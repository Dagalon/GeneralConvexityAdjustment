using NUnit.Framework;
using ql = QLNet;

namespace ConvexityAdjusmentUnitTests
{
    public class ConvexityAdjustmentHullWhite
    {
        [Test]
        public void HullWhiteFraVsFutureConvexity()
        {
            
            // settings
            var calendar = new ql.TARGET();
            
            // Hull-White parameters
            double k = 0.0003;
            double sigma = 0.0115;
            // double spreadBasis = 0.001;

            // Curves
            double flatRate = 0.01;
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
            // double deltaTime = 0.01;
            // int timeSteps = (int)(length / deltaTime);
            int timeSteps = 1;
            var numberOfMonths = 48;

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

                List<double> mandatoryPoints = new List<double> { delta00};
                double length = dc.yearFraction(startDate, t2);
                
                var rsg = (ql.InverseCumulativeRsg<ql.RandomSequenceGenerator<ql.MersenneTwisterUniformRng>, 
                        ql.InverseCumulativeNormal>)
                    new ql.PseudoRandom().make_sequence_generator(timeSteps, seed);

                ql.PathGenerator<ql.IRNG> generator = new ql.PathGenerator<ql.IRNG>(process, length, timeSteps, rsg, true, mandatoryPoints);
                double f0t = curve.currentLink().forwardRate(delta00, delta00, ql.Compounding.Simple, ql.Frequency.NoFrequency).rate();
                double expat = (1.0 - Math.Exp(-k * delta00)) / k ;
                double mt = (sigma * sigma) / (2.0 * k) * (expat - Math.Exp(-k * delta00) * expat); 
                
                // Path generator
                int i;
                var mean = 0.0;
                var meanXt0 = 0.0;
                var varXt0 = 0.0;
                var momentOrderTwoMean = 0.0;
                
                for (i = 0; i < numberOfSimulations; i++)
                {
                    ql.Sample<ql.IPath> sample = generator.next();
                    var path = (ql.Path)sample.value;
                    var rt0 = path[path.length() - 1] + mt + f0t;
                    var df01 = model.discountBond(delta00, delta01, rt0);
                    var df02 = model.discountBond(delta00, delta02, rt0);
                    var libor = ((df01 / df02) - 1.0) / (delta02 - delta01);
                    mean += (libor / numberOfSimulations);
                    meanXt0 += path[path.length() - 1] / numberOfSimulations;
                    varXt0 += (path[path.length() - 1] * path[path.length() - 1]) / numberOfSimulations;
                    momentOrderTwoMean += (libor * libor) / numberOfSimulations;
                }
                
                // Convexity future
                var varAnalyticXt = sigma * sigma * (1 - Math.Exp(-2.0 * k * delta00)) / (2.0 * k);
                var mcFuturePrice = mean;
                var stdMc = Math.Sqrt(momentOrderTwoMean - mean * mean);
                var intervalConfidence = new double[]
                    { mean - stdMc / Math.Sqrt(numberOfSimulations), mean + stdMc / Math.Sqrt(numberOfSimulations) };
                var convexityMc = mcFuturePrice - curve.currentLink()
                    .forwardRate(delta01, delta02, ql.Compounding.Simple, ql.Frequency.NoFrequency).rate();
                var convexityMalliavin =
                    ConvexityAdjustment_Lib.HullWhite.ConvexityFuture(curve, k, sigma, delta00, delta01, delta02, 0.0);
                
                // outputs
                t0s.Add(delta00);
                caMc.Add(convexityMc);
                caMalliavin.Add(convexityMalliavin);
            }
           
        }
        
    }
    



}
