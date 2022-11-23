using NUnit.Framework;
using ql = QLNet;

namespace ConvexityAdjusmentUnitTests
{
    public class CmsConvexity
    {
        [Test]
        public void HullWhiteCms()
        {
            // settings
            var calendar = new ql.TARGET();
            ql.DayCounter dc = new ql.Actual360();
            ql.Date startDate =  calendar.adjust(ql.Settings.evaluationDate());
            
            // Hull-White parameters
            double k = 0.0007;
            double sigma = 0.007;

            // Curves
            double discountFlatRate = 0.015;
            double basisSpread = 0.0;
            ql.Quote quoteDiscountRate = new ql.SimpleQuote(discountFlatRate);
            ql.Quote quoteSpread = new ql.SimpleQuote(basisSpread);
            
            ql.Handle<ql.YieldTermStructure> discountCurve =  new ql.Handle<ql.YieldTermStructure>(new ql.FlatForward(startDate, quoteDiscountRate, dc));
            
            // Model and process
            ql.OrnsteinUhlenbeckProcess process = new ql.OrnsteinUhlenbeckProcess(k, sigma, 0.0, 0.0);
            ql.HullWhite model = new ql.HullWhite(discountCurve, k, sigma);
            
            // Swap product
            var payPeriod = new ql.Period(ql.Frequency.Annual);
            var tenorSwap = new ql.Period(5, ql.TimeUnit.Years);

            // Convexity test
            var periodToCompute = new ql.Period(5, ql.TimeUnit.Years);
            int numberOftimes = 1;

            List<ql.Date> datesToCompute = new List<ql.Date>{calendar.advance(startDate, periodToCompute)};
            List<double> deltaTimes = new List<double>{dc.yearFraction(startDate, datesToCompute[0])};
            List<double> f0ts = new List<double>{discountCurve.link.forwardRate(deltaTimes[0], deltaTimes[0], ql.Compounding.Simple, ql.Frequency.NoFrequency).rate()};
            List<double> hjmAdjustment = new List<double>{ConvexityAdjustment_Lib.HullWhite.hjmAdjustment(deltaTimes[0], k, sigma)};
            var drift = new List<double>{ConvexityAdjustment_Lib.HullWhite.forwardMeasureAdjustment(deltaTimes[0], k, sigma)};
            
            for (int i = 1; i < numberOftimes; i++)
            {
                datesToCompute.Add(calendar.advance(datesToCompute[i-1], periodToCompute));
                deltaTimes.Add(dc.yearFraction(startDate, datesToCompute[i]));
                f0ts.Add(discountCurve.link.forwardRate(deltaTimes[i], deltaTimes[i], ql.Compounding.Simple, ql.Frequency.NoFrequency).rate());
                hjmAdjustment.Add(ConvexityAdjustment_Lib.HullWhite.hjmAdjustment(deltaTimes[i], k, sigma));
                drift.Add(ConvexityAdjustment_Lib.HullWhite.forwardMeasureAdjustment(deltaTimes[i], k, sigma));
            }
            
            // MC
            int numberOfSimulations = 100000;
            ulong seed = 123545;
            
            // outputs
            var rateCmsMc = new Dictionary<int, double>();
            var convexityAdjustmentMc = new List<double>();
            var convexityAdjustmentAnalytic = new List<double>();
            var std = new List<double>();
            var momentOrderTwo = new Dictionary<int, double>();
            var rateCmsForward = new List<double>();
            var annuity = new List<Tuple<double,double>>();
            var partialOisSwapTa = new List<double>();
            var partialOisSwapT0 = new List<double>();
            var partialVanillaSwap = new List<double>();
            var swapRates = new List<double>();
            var swapOisRates = new List<double>();
            
            // Test
            for (int j = 0; j < numberOftimes; j++)
            {
              
                ql.Handle<ql.YieldTermStructure> staticForwardCurve =new ql.Handle<ql.YieldTermStructure>(
                    new ql.ForwardSpreadedTermStructure(discountCurve,  new ql.Handle<ql.Quote>(quoteSpread)));
                var staticIndex = new ql.Euribor6M(staticForwardCurve);
                        
                var discountEngine = new ql.DiscountingSwapEngine(discountCurve, null, datesToCompute[j], startDate);
                        
                var staticSwap = new ql.MakeVanillaSwap(tenorSwap, staticIndex, 1.0)
                    .withEffectiveDate(datesToCompute[j])
                    .withFloatingLegTenor(staticIndex.tenor())
                    .withFixedLegTenor(new ql.Period(1, ql.TimeUnit.Years))
                    .withPricingEngine(discountEngine)
                    .withFixedLegRule(ql.DateGeneration.Rule.Forward)
                    .withFloatingLegRule(ql.DateGeneration.Rule.Forward)
                    .value();

                annuity.Add( model.GetAnnuity(startDate, startDate, f0ts[j], staticSwap.fixedSchedule(),
                    staticSwap.fixedDayCount()));
                
                var swapOisDerivativeT0 = model.GetPartialDerivativeSwapRate(startDate, startDate, f0ts[j],
                    staticSwap.floatingSchedule(), staticSwap.fixedSchedule(), staticSwap.floatingDayCount(),
                    staticSwap.fixedDayCount());
                
                var swapOisDerivativeTa = model.GetPartialDerivativeSwapRate(startDate, startDate, f0ts[j],
                    staticSwap.floatingSchedule(), staticSwap.fixedSchedule(), staticSwap.floatingDayCount(),
                    staticSwap.fixedDayCount());
                
                var swapVanillaDerivative = model.GetPartialDerivativeSwapRate(startDate, startDate, f0ts[j],
                    staticSwap.floatingSchedule(), staticForwardCurve, staticSwap.fixedSchedule(), staticSwap.floatingDayCount(),
                    staticSwap.fixedDayCount());
                
                partialOisSwapT0.Add(swapOisDerivativeT0);
                partialOisSwapTa.Add(swapOisDerivativeTa);
                partialVanillaSwap.Add(swapVanillaDerivative);
                swapRates.Add(staticSwap.fairRate());
                
                // ois swap rate
                var tb =  calendar.advance(datesToCompute[j], tenorSwap,
                    ql.BusinessDayConvention.ModifiedFollowing);
                var oisSwapRate = (1.0 - discountCurve.link.discount(tb) / discountCurve.link.discount(datesToCompute[j])) /  annuity[j].Item1; 
                swapOisRates.Add(oisSwapRate);
                
                var payDate = calendar.advance(datesToCompute[j], payPeriod,
                    ql.BusinessDayConvention.ModifiedFollowing);
                
                double payoffDiscount = discountCurve.link.discount(payDate) * staticSwap.fairRate();
                rateCmsForward.Add(payoffDiscount);
            }

            for (int j = 0; j < numberOftimes; j++)
            {
                // Path Generator
                var rsg = (ql.InverseCumulativeRsg<ql.RandomSequenceGenerator<ql.MersenneTwisterUniformRng>, 
                        ql.InverseCumulativeNormal>)
                    new ql.PseudoRandom().make_sequence_generator(1, seed);
            
                ql.PathGenerator<ql.IRNG> generator = new ql.PathGenerator<ql.IRNG>(process, deltaTimes[j], 1, rsg, false);

                
                // Pay date
                var tP = calendar.advance(datesToCompute[j], payPeriod,
                    ql.BusinessDayConvention.ModifiedFollowing);
                var deltaTimeTp = dc.yearFraction(startDate, tP);
                
                // Dynamics curve
                ql.Handle<ql.YieldTermStructure> hullWhiteDiscountCurve = new ql.Handle<ql.YieldTermStructure>(new ql.ShortRate1FYieldStructure<ql.HullWhite>(datesToCompute[j], dc, calendar, model));
                
                // Pricing Engines
                var discountEngineHullWhite = new ql.DiscountingSwapEngine(hullWhiteDiscountCurve, null, datesToCompute[j], datesToCompute[j]);
                
                ql.Handle<ql.YieldTermStructure> forwardCurve =new ql.Handle<ql.YieldTermStructure>(
                   new ql.ForwardSpreadedTermStructure(hullWhiteDiscountCurve,  new ql.Handle<ql.Quote>(quoteSpread)));
                
                var index = new ql.Euribor6M(forwardCurve);
                
                var swap = new ql.MakeVanillaSwap(tenorSwap, index, 1.0)
                    .withSettlementDays(0)
                    .withEffectiveDate(datesToCompute[j])
                    .withFloatingLegTenor(index.tenor())
                    .withFixedLegTenor(new ql.Period(1, ql.TimeUnit.Years))
                    .withPricingEngine(discountEngineHullWhite) 
                    .withFixedLegRule(ql.DateGeneration.Rule.Forward)
                    .withFloatingLegRule(ql.DateGeneration.Rule.Forward)
                    .value();
                
                var dfOis = discountCurve.link.discount(datesToCompute[j]);
                var expectedValue = 0.0;
                var dfMc = 0.0;
                var rtaMc = 0.0;
                var meanXi = 0.0;

                rateCmsMc[datesToCompute[j].serialNumber()] = 0.0;
                momentOrderTwo[datesToCompute[j].serialNumber()] = 0.0;
                
                for (int i = 0; i < numberOfSimulations; i++)
                {
                    ql.Sample<ql.IPath> sample = generator.next();
                    var path = (ql.Path)sample.value;
                    
                    var ri = path[1] + hjmAdjustment[j] + f0ts[j] - drift[j];
                    expectedValue += ri / numberOfSimulations; 
                    
                    ((ql.ShortRate1FYieldStructure<ql.HullWhite>)hullWhiteDiscountCurve.link).updateState(deltaTimes[j], ri);

                    var aux = model.discountBond(deltaTimes[j], deltaTimeTp, ri);
                    double pTaTp = hullWhiteDiscountCurve.link.discount(tP);
                    swap.recalculate();
                    var swapRate = swap.fairRate();
                  
                    double payoffMc = dfOis * swapRate * pTaTp;
                    dfMc += dfOis * pTaTp / numberOfSimulations;
                    rtaMc += swapRate / numberOfSimulations;
                    meanXi += path[1] / numberOfSimulations;;
                    
                    rateCmsMc[datesToCompute[j].serialNumber()] +=  payoffMc / numberOfSimulations;
                    momentOrderTwo[datesToCompute[j].serialNumber()] += (payoffMc * payoffMc / numberOfSimulations);
                }

                var error = dfMc - discountCurve.link.discount(tP);
                var errorRi = rtaMc * discountCurve.link.discount(datesToCompute[j]) -  swapRates[j];
                
                // Analytic convexity adjustment
                double ca = ConvexityAdjustment_Lib.HullWhite.convexityCms(discountCurve, startDate, datesToCompute[j],
                    swap.maturityDate(), tP, annuity[j].Item1, annuity[j].Item2, partialVanillaSwap[j], 
                    partialOisSwapTa[j], partialOisSwapT0[j], dc, k, sigma);
                
                convexityAdjustmentAnalytic.Add(ca);
            }

            for (int i = 0; i < numberOftimes; i++)
            {
                convexityAdjustmentMc.Add(rateCmsMc[datesToCompute[i].serialNumber()] - rateCmsForward[i]);
                std.Add( (Math.Sqrt(momentOrderTwo[datesToCompute[i].serialNumber()] - Math.Pow(rateCmsMc[datesToCompute[i].serialNumber()] ,2.0)))/ Math.Sqrt(numberOfSimulations));
            }
        }
    }
}