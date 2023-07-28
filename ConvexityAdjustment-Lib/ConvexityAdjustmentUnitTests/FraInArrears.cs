using ConvexityAdjustment_Lib.HullWhite;
using NUnit.Framework;
using ql = QLNet;

namespace ConvexityAdjustmentUnitTests
{
    public class ConvexityFraInArrearsNewVersion
    {
        [Test]
        public void hullWhiteFraInArrearsNewVersion()
        {
            
            // settings
            var calendar = new ql.TARGET();
            
            // Hull-White parameters
            double k = 0.0007;
            double sigma = 0.01;

            // Curves
            double flatRate = 0.01;
            ql.Quote quoteRate = new ql.SimpleQuote(flatRate);
            ql.DayCounter dc = new ql.Actual365Fixed();
            ql.Date startDate = ql.Settings.evaluationDate();
            ql.Handle<ql.YieldTermStructure> curve =  new ql.Handle<ql.YieldTermStructure>(new ql.FlatForward(startDate, quoteRate, dc));
            
            // Model and process
            ql.HullWhite model = new ql.HullWhite(curve, k, sigma);
            
            // MC
            int numberOfSimulations = 500000;
            ulong seed = 123545;
            var numberOfMonths = 240;

            List<double> t0s = new List<double>();
            List<double> caMc = new List<double>();
            List<double> caMalliavin = new List<double>();

            for (var j = numberOfMonths; j < numberOfMonths + 1; j++)
            {
                // Product
                var t1 = calendar.advance(startDate, j, ql.TimeUnit.Months);
                var t2 = calendar.advance(t1, 12, ql.TimeUnit.Months);
                 
                // Delta times
                double delta01 = dc.yearFraction(startDate, t1);
                double delta02 = dc.yearFraction(startDate, t2);
                
                var meanPayOff = 0.0;
                var momentOrderTwoMean = 0.0;
                
                
                // path simulation
                var libors = Sampling.getLiborForwardMeasurePaths(model, t1, t2, t1, dc, seed, numberOfSimulations);
                
                for (var i = 0; i < numberOfSimulations; i++)
                {
                    var ratio = curve.link.discount(delta02) / curve.link.discount(delta01);
                    var payOff = ratio * libors[i] * (1.0 + (delta02 - delta01) * libors[i]);
                    meanPayOff += payOff / numberOfSimulations;
                }
                
                // Statistics
                var mcFraInArrears = meanPayOff;
                var stdMc = Math.Sqrt(momentOrderTwoMean - meanPayOff * meanPayOff);
                var intervalConfidence = new double[]
                    { meanPayOff - stdMc / Math.Sqrt(numberOfSimulations), meanPayOff + stdMc / Math.Sqrt(numberOfSimulations) };
                
                // Convexity by MC and Malliavin
                var convexityMc = mcFraInArrears - curve.currentLink()
                    .forwardRate(delta01, delta02, ql.Compounding.Continuous, ql.Frequency.NoFrequency).rate();
                var convexityMalliavin =
                    HullWhite.convexityFraInArrears(curve, k, sigma, delta01, delta02, 0.0);
                
                // outputs
                t0s.Add(delta01);
                caMc.Add(convexityMc);
                caMalliavin.Add(convexityMalliavin);
            }
           
        }
        
    }
}