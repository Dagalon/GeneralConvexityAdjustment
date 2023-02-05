using ConvexityAdjustment_Lib;
using ConvexityAdjustment_Lib.HullWhite;
using ql = QLNet;
using NUnit.Framework;

namespace ConvexityAdjusmentUnitTests
{
    public class ConvexityAdjustmentSofr
    {
        [Test]
        public void HullWhiteSofrConvexity()
        {
            
            // settings
            var calendar = new ql.TARGET();
            
            // Hull-White parameters
            double k = 0.0001;
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
            // double deltaTime = 0.01;
            // int timeSteps = (int)(length / deltaTime);
            // int timeSteps = 1;
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

                var rsg = (ql.InverseCumulativeRsg<ql.RandomSequenceGenerator<ql.MersenneTwisterUniformRng>,
                        ql.InverseCumulativeNormal>)
                    new ql.PseudoRandom().make_sequence_generator(numberOfSimulations, seed);

                double p0t1 = curve.currentLink().discount(delta01, true);
                double p0t2 = curve.currentLink().discount(delta02, true);
                double r0t = -Math.Log(p0t2 / p0t1);
                
                // Path generator
                int i;
                var mean = 0.0;
                var avgMean = 0.0;
                var meanXt0 = 0.0;
                var momentOrderTwoMean = 0.0;

                var zk = rsg.nextSequence();
                var mu = HullWhite.getExpectedIntegralRt(k, sigma, delta01, delta02);
                var std = Math.Sqrt(HullWhite.getVarianceIntegralRt(k, sigma, delta01, delta02));
                
                for (i = 0; i < numberOfSimulations; i++)
                {
                    var it0Tot1 = mu + std * zk.value[i] + r0t;
                    var libor = (Math.Exp(it0Tot1) - 1.0) / (delta02 - delta01);
                    var avgLibor = it0Tot1 / (delta02 - delta01);
                    mean += (libor / numberOfSimulations);
                    avgMean += (avgLibor / numberOfSimulations);
                    meanXt0 += it0Tot1 / numberOfSimulations;
                    momentOrderTwoMean += (libor * libor) / numberOfSimulations;
                }
                
                
                // Convexity future
                var mcFuturePrice = mean;
                var stdMc = Math.Sqrt(momentOrderTwoMean - mean * mean);
                var intervalConfidence = new double[]
                    { mean - stdMc / Math.Sqrt(numberOfSimulations), mean + stdMc / Math.Sqrt(numberOfSimulations) };
                
                var forward =  curve.currentLink()
                    .forwardRate(delta01, delta02, ql.Compounding.Continuous, ql.Frequency.NoFrequency).rate();

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