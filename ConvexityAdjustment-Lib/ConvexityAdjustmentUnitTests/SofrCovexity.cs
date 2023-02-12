using NUnit.Framework;
using HullWhite = ConvexityAdjustment_Lib.HullWhite.HullWhite;
using HullWhiteSampling = ConvexityAdjustment_Lib.HullWhite.Sampling;
using ql = QLNet;

namespace ConvexityAdjustmentUnitTests
{
    public class ConvexityAdjustmentSofr
    {
        [Test]
        public void hullWhiteSofrConvexity()
        {
            // settings
            var calendar = new ql.TARGET();
            
            // Hull-White parameters
            double k = 0.0001;
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
            int numberOfSimulations = 2000000;
            ulong seed = 123545;
            var numberOfMonths = 240;

            List<double> t0s = new List<double>();
            List<double> caMc = new List<double>();
            List<double> caAvgMc = new List<double>();
            List<double> caMalliavin = new List<double>();
            List<double> caAvgMalliavin = new List<double>();

            for (var j = 1; j < numberOfMonths + 1; j++)
            {
                
                
                // Product
                var t0 = calendar.advance(startDate, j, ql.TimeUnit.Months);
                var t1 = calendar.advance(startDate, j, ql.TimeUnit.Months);
                var t2 = calendar.advance(t1, 12, ql.TimeUnit.Months);
                
                // Delta times
                double delta00 = dc.yearFraction(startDate, t0);
                double delta01 = dc.yearFraction(startDate, t1);
                double delta02 = dc.yearFraction(startDate, t2);
                
                // paths
                var paths = HullWhiteSampling.getIntRtPathSpotMeasure(model, delta02, delta01, seed,
                    numberOfSimulations);

                // Path generator
                var mean = 0.0;
                var avgMean = 0.0;
                var momentOrderTwoMean = 0.0;

                for (var i = 0; i < numberOfSimulations; i++)
                {
                    mean += (paths[i,0] / numberOfSimulations);
                    avgMean += (paths[i,1] / numberOfSimulations);
                    momentOrderTwoMean += (paths[i,0] * paths[i,0]) / numberOfSimulations;
                }
                
                // Statistics
                var mcFuturePrice = mean;
                var stdMc = Math.Sqrt(momentOrderTwoMean - mean * mean);
                var intervalConfidence = new double[]
                    { mean - stdMc / Math.Sqrt(numberOfSimulations), mean + stdMc / Math.Sqrt(numberOfSimulations) };
                
                // Convexity future
                var dfTa = curve.link.discount(delta01);
                var dfTb = curve.link.discount(delta02);
                var forward = ((dfTa / dfTb) - 1.0) / (delta02 - delta01);
                
                // var forward =  curve.currentLink().
                //     forwardRate(delta01, delta02, ql.Compounding.Compounded).rate();

                var convexityMc = mcFuturePrice - forward;
                var convexityAvgMc = avgMean - forward;
                
                var convexityMalliavin = HullWhite.convexityOis(curve, k, sigma, delta01, delta02) - forward;
                var convexityAvgMalliavin = HullWhite.convexityAvgOis(curve, k, sigma, delta01, delta02) - forward;
                
                // outputs
                t0s.Add(delta00);
                caMc.Add(convexityMc);
                caAvgMc.Add(convexityAvgMc);
                caMalliavin.Add(convexityMalliavin);
                caAvgMalliavin.Add(convexityAvgMalliavin);
            }
           
        }
        
    }
    



}