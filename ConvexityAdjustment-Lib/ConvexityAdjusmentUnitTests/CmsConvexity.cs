using NUnit.Framework;
using ql = QLNet;

namespace ConvexityAdjustmentUnitTests
{
    public class CmsConvexity
    {
        [Test]
        public void HullWhiteCms()
        {
            // settings
            var calendar = new ql.TARGET();
            ql.DayCounter dc = new ql.Actual365Fixed();
            ql.Date startDate = ql.Settings.evaluationDate();
            
            // Hull-White parameters
            double k = 0.0007;
            double sigma = 0.015;

            // Curves
            double discountFlatRate = 0.015;
            double basisSpread = 0.01;
            ql.Quote quoteDiscountRate = new ql.SimpleQuote(discountFlatRate);
            ql.Quote quoteSpread = new ql.SimpleQuote(basisSpread);
            
            ql.Handle<ql.YieldTermStructure> discountCurve =  new ql.Handle<ql.YieldTermStructure>(new ql.FlatForward(startDate, quoteDiscountRate, dc));
            
            // Model and process
            ql.OrnsteinUhlenbeckProcess process = new ql.OrnsteinUhlenbeckProcess(k, sigma, 0.0, 0.0);
            ql.HullWhite model = new ql.HullWhite(discountCurve, k, sigma);
            
            
            // Swap product
            var floatingPeriod = new ql.Period(ql.Frequency.Semiannual);
            var payPeriod = new ql.Period(ql.Frequency.Annual);
            var tenorSwap = new ql.Period(5, ql.TimeUnit.Years);

            // Convexity test
            var periodToCompute = new ql.Period(1, ql.TimeUnit.Years);
            int numberOftimes = 25;
            
            List<ql.Date> datesToCompute = new List<ql.Date>(){calendar.advance(startDate, periodToCompute)};
            List<double> deltaTimes = new List<double>(){dc.yearFraction(startDate, datesToCompute[0])};
            List<double> f0ts = new List<double>(){discountCurve.link.forwardRate(deltaTimes[0], deltaTimes[0], ql.Compounding.Simple, ql.Frequency.NoFrequency).rate()};
            List<double> hjmAdjustment = new List<double>(){ConvexityAdjustment_Lib.HullWhite.HjmAdjustment(0, deltaTimes[0], k, sigma)};
            
            for (int i = 1; i < numberOftimes; i++)
            {
                datesToCompute.Add(calendar.advance(datesToCompute[i-1], periodToCompute));
                deltaTimes.Add(dc.yearFraction(startDate, datesToCompute[i]));
                f0ts.Add(discountCurve.link.forwardRate(deltaTimes[i], deltaTimes[i], ql.Compounding.Simple, ql.Frequency.NoFrequency).rate());
                hjmAdjustment.Add(ConvexityAdjustment_Lib.HullWhite.HjmAdjustment(deltaTimes[i-1], deltaTimes[i], k, sigma));
            }
            
            // MC
            int numberOfSimulations = 500000;
            ulong seed = 123545;
            int timeSteps = deltaTimes.Count;
            
            // Path Generator
            var rsg = (ql.InverseCumulativeRsg<ql.RandomSequenceGenerator<ql.MersenneTwisterUniformRng>, 
                    ql.InverseCumulativeNormal>)
                new ql.PseudoRandom().make_sequence_generator(timeSteps, seed);
            
            ql.PathGenerator<ql.IRNG> generator = new ql.PathGenerator<ql.IRNG>(process, deltaTimes.Last(), timeSteps, rsg, true, deltaTimes);
            
            // outputs
            var rateCmsMc = new Dictionary<int, double>();
            var convexityAdjustment = new List<double>();
            var std = new List<double>();
            var momentOrderTwo = new Dictionary<int, double>();
            var rateCmsForward = new List<double>();
            
            // Test
            // Swap rate at t=0
            for (int j = 0; j < numberOftimes; j++)
            {
                
                var endDate = calendar.advance(datesToCompute[j], tenorSwap,
                    ql.BusinessDayConvention.ModifiedFollowing);

                ql.Handle<ql.YieldTermStructure> staticForwardCurve =new ql.Handle<ql.YieldTermStructure>(
                    new ql.ForwardSpreadedTermStructure(discountCurve,  new ql.Handle<ql.Quote>(quoteSpread)));
                var staticIndex = new ql.Euribor6M(staticForwardCurve);
                        
                var discountEngine = new ql.DiscountingSwapEngine(discountCurve, null, startDate, startDate);
                        
                var staticSwap = new ql.MakeVanillaSwap(floatingPeriod, staticIndex, 0.0)
                    .withEffectiveDate(datesToCompute[j])
                    .withFloatingLegTenor(staticIndex.tenor())
                    .withFixedLegTenor(new ql.Period(1, ql.TimeUnit.Years))
                    .withTerminationDate(endDate)
                    .withPricingEngine(discountEngine).value();
                        
                        
                double payoffDiscount = discountCurve.link.discount(datesToCompute[j]) * staticSwap.fairRate();
                rateCmsForward.Add(payoffDiscount);
            }

            for (int i = 0; i < numberOfSimulations; i++)
            {
                ql.Sample<ql.IPath> sample = generator.next();
                var path = (ql.Path)sample.value;

                for (int j = 0; j < numberOftimes; j++)
                {
                    var ri = path[j + 1] + hjmAdjustment[j] + f0ts[j];
                    
                    var endDate = calendar.advance(datesToCompute[j], tenorSwap,
                        ql.BusinessDayConvention.ModifiedFollowing);
                    
                    // Pay date
                    var tP = calendar.advance(datesToCompute[j], payPeriod,
                        ql.BusinessDayConvention.ModifiedFollowing);
                    
                    // Dynamics curve
                    ql.Handle<ql.YieldTermStructure> hullWhiteDiscountCurve = new ql.Handle<ql.YieldTermStructure>(new ql.ShortRate1FYieldStructure<ql.HullWhite>(datesToCompute[j], dc, calendar, model));
                    var dt = dc.yearFraction(startDate, datesToCompute[j]);
                    ((ql.ShortRate1FYieldStructure<ql.HullWhite>)hullWhiteDiscountCurve.link).updateState(dt, ri);
                    
                    // Pricing Engines
                    var discountEngineHullWhite = new ql.DiscountingSwapEngine(hullWhiteDiscountCurve, null, datesToCompute[j], datesToCompute[j]);
                   
                    ql.Handle<ql.YieldTermStructure> forwardCurve =new ql.Handle<ql.YieldTermStructure>(
                       new ql.ForwardSpreadedTermStructure(hullWhiteDiscountCurve,  new ql.Handle<ql.Quote>(quoteSpread)));
                    
                    var index = new ql.Euribor6M(forwardCurve);
                   
                    var swap = new ql.MakeVanillaSwap(floatingPeriod, index, 0.0)
                        .withEffectiveDate(datesToCompute[j])
                        .withFloatingLegTenor(index.tenor())
                        .withFixedLegTenor(new ql.Period(1, ql.TimeUnit.Years))
                        .withTerminationDate(endDate)
                        .withPricingEngine(discountEngineHullWhite).value();

                    double payoffMc = discountCurve.link.discount(dt) * swap.fairRate();
                   

                    if (rateCmsMc.ContainsKey(datesToCompute[j].serialNumber()))
                    {
                        rateCmsMc[datesToCompute[j].serialNumber()] += payoffMc / numberOfSimulations;
                        momentOrderTwo[datesToCompute[j].serialNumber()] += (payoffMc * payoffMc / numberOfSimulations);
                    }
                    else
                    {
                        rateCmsMc.Add(datesToCompute[j].serialNumber(), payoffMc / numberOfSimulations);
                        momentOrderTwo.Add(datesToCompute[j].serialNumber(), payoffMc * payoffMc / numberOfSimulations);
                    }
                    
                }
            }

            for (int i = 0; i < numberOftimes; i++)
            {
                convexityAdjustment.Add(rateCmsMc[datesToCompute[i].serialNumber()] - rateCmsForward[i]);
                std.Add( (Math.Sqrt(momentOrderTwo[datesToCompute[i].serialNumber()] - Math.Pow(rateCmsMc[datesToCompute[i].serialNumber()] ,2.0)))/ Math.Sqrt(numberOfSimulations));
            }

        }
        
       
        
    }
}