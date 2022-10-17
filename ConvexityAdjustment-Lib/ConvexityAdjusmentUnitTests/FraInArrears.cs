using NUnit.Framework;
using ql = QLNet;

namespace ConvexityAdjustmentUnitTests
{
    public class ConvexityFraInArrears
    {
        [Test]
        public void HullWhiteFraInArrears()
        {
            
            // settings
            var calendar = new ql.TARGET();
            
            // Hull-White parameters
            double k = 0.0007;
            double sigma = 0.01;
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
            int timeSteps = 1;
            var numberOfMonths = 48;

            List<double> t0s = new List<double>();
            List<double> caMc = new List<double>();
            List<double> caMalliavin = new List<double>();

            for (var j = 1; j < numberOfMonths + 1; j++)
            {
                // Product
                var t1 = calendar.advance(startDate, j, ql.TimeUnit.Months);
                var t2 = calendar.advance(t1, 12, ql.TimeUnit.Months);
                 
                // Delta times
                double delta01 = dc.yearFraction(startDate, t1);
                double delta02 = dc.yearFraction(startDate, t2);

                List<double> mandatoryPoints = new List<double> { delta01};
                double length = dc.yearFraction(startDate, t2);
                
                var rsg = (ql.InverseCumulativeRsg<ql.RandomSequenceGenerator<ql.MersenneTwisterUniformRng>, 
                        ql.InverseCumulativeNormal>)
                    new ql.PseudoRandom().make_sequence_generator(timeSteps, seed);

                ql.PathGenerator<ql.IRNG> generator = new ql.PathGenerator<ql.IRNG>(process, length, timeSteps, rsg, true, mandatoryPoints);
                double f0t = curve.currentLink().forwardRate(delta01, delta01, ql.Compounding.Simple, ql.Frequency.NoFrequency).rate();
                double expat = (1.0 - Math.Exp(-k * delta01)) / k ;
                double mt = (sigma * sigma) / (2.0 * k) * (expat - Math.Exp(-k * delta01) * expat);
                
                // Adjustment forward measure
                double adjFwdMeasure =
                    -(sigma * sigma / k) * (ConvexityAdjustment_Lib.HullWhite.beta(0.0, delta01, k) - Math.Exp(-k * delta01) * delta01);
                
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
                    var rt0 = path[path.length() - 1] + mt + f0t + adjFwdMeasure;
                    var df12 = model.discountBond(delta01, delta02, rt0);
                    var libor = (1.0 / df12 - 1.0) / (delta02 - delta01);
                    mean += (libor / numberOfSimulations);
                    meanXt0 += path[path.length() - 1] / numberOfSimulations;
                    varXt0 += (path[path.length() - 1] * path[path.length() - 1]) / numberOfSimulations;
                    momentOrderTwoMean += (libor * libor) / numberOfSimulations;
                }
                
                
                // Convexity future
                var varAnalyticXt = sigma * sigma * (1 - Math.Exp(-2.0 * k * delta01)) / (2.0 * k);
                var mcFraInArrears = mean;
                var stdMc = Math.Sqrt(momentOrderTwoMean - mean * mean);
                var intervalConfidence = new double[]
                    { mean - stdMc / Math.Sqrt(numberOfSimulations), mean + stdMc / Math.Sqrt(numberOfSimulations) };
                var convexityMc = mcFraInArrears - curve.currentLink()
                    .forwardRate(delta01, delta02, ql.Compounding.Simple, ql.Frequency.NoFrequency).rate();
                var convexityMalliavin =
                    ConvexityAdjustment_Lib.HullWhite.convexityFraInArrears(curve, k, sigma, delta01, delta02, 0.0);
                
                // outputs
                t0s.Add(delta01);
                caMc.Add(convexityMc);
                caMalliavin.Add(convexityMalliavin);
            }
           
        }
        
    }
}
