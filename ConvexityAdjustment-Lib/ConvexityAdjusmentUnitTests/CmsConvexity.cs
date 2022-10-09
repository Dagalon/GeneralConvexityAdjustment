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
            double sigma = 0.01;

            // Curves
            double discountFlatRate = 0.015;
            double basisSpread = 0.01;
            ql.Quote quoteDiscountRate = new ql.SimpleQuote(discountFlatRate);
            ql.Quote quoteSpread = new ql.SimpleQuote(basisSpread);
            
            ql.Handle<ql.YieldTermStructure> discountCurve =  new ql.Handle<ql.YieldTermStructure>(new ql.FlatForward(startDate, quoteDiscountRate, dc));
            
            // Model and process
            ql.OrnsteinUhlenbeckProcess process = new ql.OrnsteinUhlenbeckProcess(k, sigma, 0.0, 0.0);
            ql.HullWhite model = new ql.HullWhite(discountCurve, k, sigma);
            
            // Dynamics curve
            ql.Handle<ql.YieldTermStructure> hullWhiteDiscountCurve = new ql.Handle<ql.YieldTermStructure>(new ql.ShortRate1FYieldStructure<ql.HullWhite>(model));
            
            // Swap product
            var floatingPeriod = new ql.Period(ql.Frequency.Semiannual);
            var payPeriod = new ql.Period(ql.Frequency.Annual);
            var tenorSwap = new ql.Period(5, ql.TimeUnit.Years);

            var discountEngine = new ql.DiscountingSwapEngine(discountCurve, null, startDate, startDate);
            
            // Convexity test
            var periodToCompute = new ql.Period(1, ql.TimeUnit.Years);
            int numberOftimes = 25;
            
            List<ql.Date> datesToCompute = new List<ql.Date>(){calendar.advance(startDate, periodToCompute)};
            List<double> deltaTimes = new List<double>(){dc.yearFraction(startDate, datesToCompute[0])};
            List<double> f0ts = new List<double>(){discountCurve.link.forwardRate(deltaTimes[0], deltaTimes[0], ql.Compounding.Continuous, ql.Frequency.NoFrequency).rate()};
            List<double> hjmAdjustment = new List<double>(){ConvexityAdjustment_Lib.HullWhite.HjmAdjustment(0, deltaTimes[0], k, sigma)};
            
            for (int i = 1; i < numberOftimes; i++)
            {
                datesToCompute.Add(calendar.advance(datesToCompute[i-1], periodToCompute));
                deltaTimes.Add(dc.yearFraction(startDate, datesToCompute[i]));
                f0ts.Add(discountCurve.link.forwardRate(deltaTimes[i], deltaTimes[i], ql.Compounding.Continuous, ql.Frequency.NoFrequency).rate());
                hjmAdjustment.Add(ConvexityAdjustment_Lib.HullWhite.HjmAdjustment(deltaTimes[i-1], deltaTimes[i], k, sigma));
            }
            
            // MC
            int numberOfSimulations = 10;
            ulong seed = 123545;
            int timeSteps = deltaTimes.Count;
            
            // Path Generator
            var rsg = (ql.InverseCumulativeRsg<ql.RandomSequenceGenerator<ql.MersenneTwisterUniformRng>, 
                    ql.InverseCumulativeNormal>)
                new ql.PseudoRandom().make_sequence_generator(timeSteps, seed);
            
            ql.PathGenerator<ql.IRNG> generator = new ql.PathGenerator<ql.IRNG>(process, deltaTimes.Last(), timeSteps, rsg, true, deltaTimes);
            
            // outputs
            var rateCmsMc = new Dictionary<int, double>();
            var rateCmsForward = new List<double>();

            // Test
            for (int i = 0; i < numberOfSimulations; i++)
            {
                ql.Sample<ql.IPath> sample = generator.next();
                var path = (ql.Path)sample.value;

                for (int j = 1; j < numberOftimes; j++)
                {
                    var ri = path[j] + hjmAdjustment[j - 1] + f0ts[j-1];
                    var endDate = calendar.advance(datesToCompute[j - 1], tenorSwap,
                        ql.BusinessDayConvention.ModifiedFollowing);
                    
                    // Pay date
                    var tP = calendar.advance(datesToCompute[j - 1], payPeriod,
                        ql.BusinessDayConvention.ModifiedFollowing);
                    
                    // Curves
                    var dt = dc.yearFraction(startDate, datesToCompute[j - 1]);
                    ((ql.ShortRate1FYieldStructure<ql.HullWhite>)hullWhiteDiscountCurve.link).updateState(dt, ri);
                    
                    // Pricing Engines
                    var discountEngineHullWhite = new ql.DiscountingSwapEngine(hullWhiteDiscountCurve, null, datesToCompute[j-1], datesToCompute[j-1]);
                   
                    ql.Handle<ql.YieldTermStructure> forwardCurve =new ql.Handle<ql.YieldTermStructure>(
                       new ql.ForwardSpreadedTermStructure(hullWhiteDiscountCurve,  new ql.Handle<ql.Quote>(quoteSpread)));
                    var index = new ql.Euribor6M(forwardCurve);
                   
                    var swap = new ql.MakeVanillaSwap(floatingPeriod, index, 0.0)
                        .withEffectiveDate(datesToCompute[j-1])
                        .withFloatingLegTenor(index.tenor())
                        .withFixedLegTenor(new ql.Period(1, ql.TimeUnit.Years))
                        .withTerminationDate(endDate)
                        .withPricingEngine(discountEngineHullWhite).value();

                    double payoffMc = discountCurve.link.discount(dt) * swap.fairRate();

                    if (rateCmsMc.ContainsKey(datesToCompute[j - 1].serialNumber()))
                    {
                        rateCmsMc[datesToCompute[j - 1].serialNumber()] += payoffMc / numberOfSimulations;
                    }
                    else
                    {
                        rateCmsMc.Add(datesToCompute[j - 1].serialNumber(), payoffMc / numberOfSimulations);
                    }
                    
                    if (i == 0)
                    {
                        swap.setPricingEngine(discountEngine);
                        double payoffDiscount = swap.fairRate();
                        rateCmsForward.Add(payoffDiscount);
                    }

                }
                
            }

        }
    }
}