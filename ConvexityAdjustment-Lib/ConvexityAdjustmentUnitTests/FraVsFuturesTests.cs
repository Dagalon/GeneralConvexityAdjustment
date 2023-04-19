using ConvexityAdjustment_Lib.HullWhite;
using NUnit.Framework;
using HullWhite = ConvexityAdjustment_Lib.HullWhite.HullWhite;
using ql = QLNet;

namespace ConvexityAdjustmentUnitTests
{
    public class ConvexityFutureFras
    {
        [Test]
        public void hullWhiteFraVsFutureConvexity()
        {
            
            // settings
            var calendar = new ql.TARGET();
            
            // Hull-White parameters
            double k = 0.0003;
            double sigma = 0.015;
            double spreadBasis = 0.001;

            // Curves
            double flatRate = 0.025;
            ql.Quote quoteRate = new ql.SimpleQuote(flatRate);
            ql.DayCounter dc = new ql.Actual365Fixed();
            ql.Date startDate = ql.Settings.evaluationDate();
            ql.Handle<ql.YieldTermStructure> curve =  new ql.Handle<ql.YieldTermStructure>(new ql.FlatForward(startDate, quoteRate, dc));
            
            // Model and process
            ql.OrnsteinUhlenbeckProcess process = new ql.OrnsteinUhlenbeckProcess(k, sigma, 0.0, 0.0);
            ql.HullWhite model = new ql.HullWhite(curve, k, sigma);
            
            // MC
            int numberOfSimulations = 2000000;
            ulong seed = 123545;
            var numberOfMonths = 240;

            List<double> t0s = new List<double>();
            List<double> caMc = new List<double>();
            List<double> caMalliavin = new List<double>();

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

                var libors = Sampling.getLiborPaths(model, t0, t1, t2, dc, seed, numberOfSimulations);
                
                // Path generator
                var mean = 0.0;
                // var meanXt0 = 0.0;
                // var varXt0 = 0.0;
                var momentOrderTwoMean = 0.0;
                
                for (var i = 0; i < numberOfSimulations; i++)
                {
                    mean += (libors[i] / numberOfSimulations);
                    momentOrderTwoMean += (libors[i] * libors[i]) / numberOfSimulations;
                }
                
                // Statistics
                var mcFuturePrice = mean;
                var stdMc = Math.Sqrt(momentOrderTwoMean - mean * mean);
                var intervalConfidence = new double[]
                    { mean - stdMc / Math.Sqrt(numberOfSimulations), mean + stdMc / Math.Sqrt(numberOfSimulations) };
                
                // Convexity
                var convexityMc = mcFuturePrice - curve.currentLink()
                    .forwardRate(delta01, delta02, ql.Compounding.Simple, ql.Frequency.NoFrequency).rate() - spreadBasis;
                var convexityMalliavin =
                    HullWhite.convexityFuture(curve, k, sigma, delta00, delta01, delta02, spreadBasis);
                
                // outputs
                t0s.Add(delta00);
                caMc.Add(convexityMc);
                caMalliavin.Add(convexityMalliavin);
            }
           
        }
        
    }
}
