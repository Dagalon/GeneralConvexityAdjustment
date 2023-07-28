using NUnit.Framework;
using HullWhite = ConvexityAdjustment_Lib.HullWhite.HullWhite;
using HullWhiteSampling = ConvexityAdjustment_Lib.HullWhite.Sampling;
using ql = QLNet;

namespace ConvexityAdjustmentUnitTests
{
    public class RfrCapPriceHullWhite
    {
        [Test]
        public void hullWhiteSofrConvexity()
        {
            // settings
            var calendar = new ql.TARGET();
            
            // Hull-White parameters
            double k = 0.0001;
            double sigma = 0.01;
            
            // cap floor information
            double[] strike = {0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.06, 0.07, 0.075, 0.08};
            int noStrikes = strike.Length;
            double maturity = 5.0;
            double tenorCaplet = 0.5;

            // Curves
            double flatRate = 0.04;
            ql.Quote quoteRate = new ql.SimpleQuote(flatRate);
            ql.DayCounter dc = new ql.Actual365Fixed();
            ql.Date startDate = ql.Settings.evaluationDate();
            ql.Handle<ql.YieldTermStructure> curve =  new ql.Handle<ql.YieldTermStructure>(new ql.FlatForward(startDate, quoteRate, dc));
            
            // Model and process
            ql.HullWhite model = new ql.HullWhite(curve, k, sigma);
            
            // MC
            int numberOfSimulations = 1000000;
            int seed = 123545;
            var numberOfMonths = 240;

            
            var forwards = HullWhiteSampling.getIntRtPathSpotMeasure(model, maturity, maturity + tenorCaplet, seed,
                numberOfSimulations);

            double[] capsPrice = new double[noStrikes];
            double[] analyticPrice = new double[noStrikes];
            var priceMomentOrerTwo = 0.0;
            var df = curve.link.discount(maturity, true);
            
            for (var i = 0; i < noStrikes; i++)
            {
                // Analytic Prices
                for (var j = 0; j < numberOfSimulations; j++)
                {
                    var pi = Math.Max(forwards[j, 0] - strike[i], 0.0);
                    capsPrice[i] += pi / numberOfSimulations;
                    priceMomentOrerTwo += pi *  pi / numberOfSimulations;
                }

                capsPrice[i] *= df;
            }

            
           
        }
        
    }
}