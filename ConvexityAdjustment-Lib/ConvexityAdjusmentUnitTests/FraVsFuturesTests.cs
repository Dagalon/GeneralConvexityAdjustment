using NUnit.Framework;
using ql = QLNet;
using System;

namespace ConvexityAdjustment_UnitTest
{

    public class ConvexityAdjustmentHullWhite
    {
        [Test]
        public void HullWhiteFraVsFutureConvexity()
        {
            
            // settings
            var calendar = new ql.TARGET();
            
            // Hull-White parameters
            double k = 0.003;
            double sigma = 0.05;
            // double spreadBasis = 0.001;

            // Curves
            double flatRate = 0.01;
            ql.Quote quoteRate = new ql.SimpleQuote(flatRate);
            ql.DayCounter dc = new ql.Actual365Fixed();
            ql.Date startDate = ql.Settings.evaluationDate();
            ql.Handle<ql.YieldTermStructure> curve =  new ql.Handle<ql.YieldTermStructure>(new ql.FlatForward(startDate, quoteRate, dc));

            // Model and process
            ql.HullWhiteProcess process = new ql.HullWhiteProcess(curve, k, sigma);
            ql.HullWhite model = new ql.HullWhite(curve, k, sigma);
            
            // Product
            var t0 = calendar.advance(startDate, 12, ql.TimeUnit.Months);
            var t1 = calendar.advance(t0, 12, ql.TimeUnit.Months);
            
            // MC
            int numberOfSimulations = 10;
            ulong seed = 123545;
            double deltaTime = 0.1;
            double length = dc.yearFraction(startDate, t1);
            int timeSteps = (int)(length / deltaTime);
            
            var rsg = (ql.InverseCumulativeRsg<ql.RandomSequenceGenerator<ql.MersenneTwisterUniformRng>, 
                       ql.InverseCumulativeNormal>)
                       new ql.PseudoRandom().make_sequence_generator(timeSteps, seed);

            ql.PathGenerator<ql.IRNG> generator = new ql.PathGenerator<ql.IRNG>(process, length, timeSteps, rsg, false);
            
            // Path generator
            int i;
            for (i = 0; i < numberOfSimulations; i++)
                generator.next();
            
            ql.Sample<ql.IPath> sample = generator.next();
            ql.Path value = sample.value as ql.Path;
        }
        
    }
    



}
