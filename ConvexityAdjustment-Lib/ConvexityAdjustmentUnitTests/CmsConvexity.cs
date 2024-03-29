﻿using NUnit.Framework;
using HullWhite = ConvexityAdjustment_Lib.HullWhite.HullWhite;
using ql = QLNet;

namespace ConvexityAdjustmentUnitTests
{
    public class CmsConvexity
    {
        [Test]
        public void hullWhiteCms()
        {
            // settings
            var calendar = new ql.TARGET();
            ql.DayCounter dc = new ql.Actual360();
            ql.Date startDate = calendar.adjust(ql.Settings.evaluationDate());

            // Hull-White parameters
            double k = 0.0007;
            double sigma = 0.01;

            // Curves
            double discountFlatRate = 0.015;
            double basisSpread = 0.0;
            ql.Quote quoteDiscountRate = new ql.SimpleQuote(discountFlatRate);
            ql.Quote quoteSpread = new ql.SimpleQuote(basisSpread);

            ql.Handle<ql.YieldTermStructure> discountCurve =
                new ql.Handle<ql.YieldTermStructure>(new ql.FlatForward(startDate, quoteDiscountRate, dc));

            // Model and process
            ql.OrnsteinUhlenbeckProcess process = new ql.OrnsteinUhlenbeckProcess(k, sigma, 0.0, 0.0);
            ql.HullWhite model = new ql.HullWhite(discountCurve, k, sigma);

            // Swap product
            var payPeriod = new ql.Period(ql.Frequency.Annual);
            var tenorSwap = new ql.Period(5, ql.TimeUnit.Years);

            // Convexity test
            var periodToCompute = new ql.Period(6, ql.TimeUnit.Months);
            int numberOfTimes = 60;

            List<ql.Date> datesToCompute = new List<ql.Date> {calendar.advance(startDate, periodToCompute)};
            List<double> deltaTimes = new List<double> {dc.yearFraction(startDate, datesToCompute[0])};
            
            List<double> f0ts = new List<double>
            {
                discountCurve.link
                    .forwardRate(deltaTimes[0], deltaTimes[0], ql.Compounding.Continuous, ql.Frequency.NoFrequency,true).rate()
            };

            for (int i = 1; i < numberOfTimes; i++)
            {
                // simulation date
                datesToCompute.Add(calendar.advance(datesToCompute[i - 1], periodToCompute));
                
                // year fracction form value date to sampling date
                deltaTimes.Add(dc.yearFraction(startDate, datesToCompute[i]));
                
                // initial forward curve
                f0ts.Add(discountCurve.link.forwardRate(deltaTimes[i], deltaTimes[i], ql.Compounding.Continuous,
                    ql.Frequency.NoFrequency).rate());
            }

            // MC
            int numberOfSimulations = 500000;
            ulong seed = 123545;

            // outputs
            var rateCmsMc = new Dictionary<int, double>();
            var convexityAdjustmentMc = new List<double>();
            var convexityAdjustmentAnalytic = new List<double>();
            var std = new List<double>();
            var momentOrderTwo = new Dictionary<int, double>();
            var rateCmsForward = new List<double>();
            
            // swap derivatives, rates ans ratio
            var partialOisSwap = new List<double>();
            
            var swapRates = new List<double>();
            var swapOisRates = new List<double>();
            var ratioSwaps = new List<double>();

            // Test
            for (int j = 0; j < numberOfTimes; j++)
            {
                ql.Handle<ql.YieldTermStructure> staticForwardCurve = new ql.Handle<ql.YieldTermStructure>(
                    new ql.ForwardSpreadedTermStructure(discountCurve, new ql.Handle<ql.Quote>(quoteSpread)));
                var staticIndex = new ql.Euribor6M(staticForwardCurve);

                var discountEngine = new ql.DiscountingSwapEngine(discountCurve, null, datesToCompute[j], startDate);

                var staticSwap = new ql.MakeVanillaSwap(tenorSwap, staticIndex, 1.0)
                    .withEffectiveDate(datesToCompute[j])
                    .withFloatingLegTenor(staticIndex.tenor())
                    .withFixedLegTenor(new ql.Period(1, ql.TimeUnit.Years))
                    .withPricingEngine(discountEngine)
                    .withFixedLegRule(ql.DateGeneration.Rule.Forward)
                    .withFloatingLegRule(ql.DateGeneration.Rule.Forward)
                    .withFloatingLegDayCount(discountCurve.link.dayCounter())
                    .value();
                
                var swapOisDerivative = model.GetPartialDerivativeSwapRate(startDate, startDate, f0ts[j],
                    staticSwap.floatingSchedule(), staticSwap.fixedSchedule(), discountCurve.link.dayCounter(),  staticSwap.floatingDayCount(),
                    staticSwap.fixedDayCount());


                var annuity = model.GetAnnuity(startDate, startDate, f0ts[j], staticSwap.fixedSchedule(),
                    discountCurve.link.dayCounter(), staticSwap.fixedDayCount());
         
                partialOisSwap.Add(swapOisDerivative);
                
                var swapOisFairRate = (discountCurve.link.discount(datesToCompute[j]) - discountCurve.link.discount(staticSwap.maturityDate())) / annuity.Item1;
                swapRates.Add(swapOisFairRate);
                
                // rateCmsForward.Add(staticSwap.fairRate());
                rateCmsForward.Add(swapOisFairRate);
            }

            for (int j = 0; j < numberOfTimes; j++)
            {
                // Path Generator
                var rsg = (ql.InverseCumulativeRsg<ql.RandomSequenceGenerator<ql.MersenneTwisterUniformRng>,
                        ql.InverseCumulativeNormal>)
                    new ql.PseudoRandom().make_sequence_generator(2, seed);

                // Pay date
                var tP = calendar.advance(datesToCompute[j], payPeriod, ql.BusinessDayConvention.ModifiedFollowing);
                
                // Dynamics curve
                ql.Handle<ql.YieldTermStructure> hullWhiteDiscountCurve =
                    new ql.Handle<ql.YieldTermStructure>(
                        new ql.ShortRate1FYieldStructure<ql.HullWhite>(datesToCompute[j], dc, calendar, model));

                // Pricing Engines
                var discountEngineHullWhite = new ql.DiscountingSwapEngine(hullWhiteDiscountCurve, null,
                    datesToCompute[j], datesToCompute[j]);

                ql.Handle<ql.YieldTermStructure> forwardCurve = new ql.Handle<ql.YieldTermStructure>(
                    new ql.ForwardSpreadedTermStructure(hullWhiteDiscountCurve, new ql.Handle<ql.Quote>(quoteSpread)));

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

                var X0Root = HullWhite.getX0(model, startDate, swap.floatingSchedule(),
                    discountCurve.link.dayCounter(), swap.fixedDayCount(), swapRates[j]);

                var annuityTa = model.GetAnnuity(datesToCompute[j], datesToCompute[j], X0Root + f0ts[j], swap.fixedSchedule(),
                    discountCurve.link.dayCounter(), swap.fixedDayCount());
                
                var annuityT0 = model.GetAnnuity(startDate, startDate, f0ts[j], swap.fixedSchedule(),
                    discountCurve.link.dayCounter(), swap.fixedDayCount());

                rateCmsMc[datesToCompute[j].serialNumber()] = 0.0;
                momentOrderTwo[datesToCompute[j].serialNumber()] = 0.0;

                // variables to check od the simulation
                var dfMC = 0.0;
                var meanI0t = 0.0;

                // Moments to simulate
                var varRt = HullWhite.getVarianceRt(deltaTimes[j], k, sigma);
                var varIt = HullWhite.getVarianceIt(deltaTimes[j], k, sigma);
                var covRtIt = HullWhite.getCovarianceRtIt(deltaTimes[j], deltaTimes[j], k, sigma);
                var driftRt = HullWhite.hjmAdjustment(deltaTimes[j], k, sigma);
                var driftIt = HullWhite.getDriftIt(deltaTimes[j], k, sigma);
                var rho = covRtIt / Math.Sqrt(varIt * varRt);
                var w1 = rho;
                var w2 = Math.Sqrt(1.0 - rho * rho);
                
                for (int i = 0; i < numberOfSimulations; i++)
                {
                    var sample = rsg.nextSequence();
                    var path = sample.value;
                    
                    // sampling r_t
                    var ri = driftRt + Math.Sqrt(varRt) * path[0] + f0ts[j];
                    
                    // Sampling It we have used Cholesky between r_t  and I_t
                    var noiseI01 = w1 * path[0] + w2 * path[1];
                    var i0Tp = driftIt + noiseI01 * Math.Sqrt(varIt);
                    var df = discountCurve.link.discount(deltaTimes[j]);
                    var dfTp = discountCurve.link.discount(tP);
                    var r0t = -Math.Log(df);
                    // var r0Tp = -Math.Log(dfTp);
                    var bt = Math.Exp(i0Tp + r0t);
                    
                    ((ql.ShortRate1FYieldStructure<ql.HullWhite>) hullWhiteDiscountCurve.link).updateState(
                        deltaTimes[j], ri);

                    double pTaTp = hullWhiteDiscountCurve.link.discount(tP);
                    swap.recalculate();
                    var swapRate = swap.fairRate();

                    // double annuity = Math.Abs(swap.fixedLegNPV());
                    // double mTa = pTaTp / annuity;

                    var annuityInfo = model.GetAnnuity(datesToCompute[j], datesToCompute[j], ri, 
                        swap.fixedSchedule(), discountCurve.link.dayCounter(), swap.fixedDayCount());

                    var dfEnd = hullWhiteDiscountCurve.link.discount(swap.maturityDate());
                    var swapRatePath = (1.0 - dfEnd) / annuityInfo.Item1;
                    
                    double payoffMc = pTaTp * swapRatePath / (bt * dfTp);

                    rateCmsMc[datesToCompute[j].serialNumber()] += payoffMc / numberOfSimulations;
                    momentOrderTwo[datesToCompute[j].serialNumber()] += (payoffMc * payoffMc / numberOfSimulations);
                    dfMC += (pTaTp / bt) / numberOfSimulations;
                    meanI0t += (i0Tp + r0t) / numberOfSimulations;
                }

                double ca = HullWhite.convexityCms(model,
                    startDate,
                    datesToCompute[j],
                    swap.maturityDate(),
                    tP,
                    1.0,
                    annuityT0.Item1,
                    annuityTa.Item1,
                    annuityTa.Item2,
                    annuityTa.Item3,
                    dc,
                    X0Root + f0ts[j]);
                
                convexityAdjustmentAnalytic.Add(ca);
            }

            for (int i = 0; i < numberOfTimes; i++)
            {
                convexityAdjustmentMc.Add(rateCmsMc[datesToCompute[i].serialNumber()] - rateCmsForward[i]);
                std.Add((Math.Sqrt(momentOrderTwo[datesToCompute[i].serialNumber()] -
                                   Math.Pow(rateCmsMc[datesToCompute[i].serialNumber()], 2.0))) /
                        Math.Sqrt(numberOfSimulations));
            }
        }
    }
}