using ql = QLNet;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.Random;

namespace ConvexityAdjustment_Lib.HullWhite;

public static class Sampling
{

    public static double[] getLiborPathsSpotMeasure(ql.HullWhite model, 
        ql.Date t0, 
        ql.Date ta, 
        ql.Date tb, 
        ql.DayCounter dc, 
        ulong seed, 
        int numberOfPaths)
    {
        // path vector
        double[] paths = new double[numberOfPaths];
        
        // yearfracs
        double dT0Ta = dc.yearFraction(t0, ta);
        double dT0Tb = dc.yearFraction(t0, tb);
        double dTaTb = dc.yearFraction(ta, tb);

        var rsg = (ql.InverseCumulativeRsg<ql.RandomSequenceGenerator<ql.MersenneTwisterUniformRng>, 
            ql.InverseCumulativeNormal>) new ql.PseudoRandom().make_sequence_generator(numberOfPaths, seed);
        
        // mean, variance, ...
        var stdXt0 = model.sigma() * Math.Sqrt(HullWhite.beta(0.0, dT0Ta, 2.0 * model.a()));
        var driftXt0 = HullWhite.hjmAdjustment(dT0Ta, model.a(), model.sigma());
        
        var sample = rsg.nextSequence();
        var zSampleValue = sample.value;
        
        double f0T = model.termStructure().currentLink().forwardRate(dT0Ta, dT0Tb, ql.Compounding.Continuous, ql.Frequency.NoFrequency).rate();

        for (int i = 0; i < numberOfPaths; i++)
        {
            var rT0 = driftXt0 + f0T + stdXt0 * zSampleValue[i];
            var dfT0Ta = model.discountBond(dT0Ta, dT0Ta, rT0);
            var dfT0Tb = model.discountBond(dT0Ta, dT0Tb, rT0);

            paths[i] = (dfT0Ta / dfT0Tb - 1.0) / dTaTb;
        }

        return paths;
        
    }


    public static double[] getLiborPaths(ql.HullWhite model, 
        ql.Date t0, 
        ql.Date ta, 
        ql.Date tb, 
        ql.DayCounter dc, 
        ulong seed, 
        int numberOfPaths)
    {
        // path vector
        double[] paths = new double[numberOfPaths];
        
        // yearfracs
        double dT0Ta = dc.yearFraction(t0, ta);
        double dT0Tb = dc.yearFraction(t0, tb);
        double dTaTb = dc.yearFraction(ta, tb);

        double sqrtdT0Ta = Math.Sqrt(dT0Ta);
        
        var rsg = (ql.InverseCumulativeRsg<ql.RandomSequenceGenerator<ql.MersenneTwisterUniformRng>, 
                ql.InverseCumulativeNormal>) new ql.PseudoRandom().make_sequence_generator(numberOfPaths, seed);
        
        // mean, variance, ...
        var stdLibor = model.sigma() * HullWhite.beta(dT0Ta, dT0Tb, model.a());
        
        var sample = rsg.nextSequence();
        var zSampleValue = sample.value;
        
        // Initial z0
        double f0T = model.termStructure().currentLink().forwardRate(dT0Ta, dT0Tb, ql.Compounding.Continuous, ql.Frequency.NoFrequency).rate();
        var z0Delta = 1.0 + dTaTb * f0T;

        for (int i = 0; i < numberOfPaths; i++)
        {
            var ztDelta = z0Delta * Math.Exp(-0.5 * stdLibor * stdLibor * dT0Ta + stdLibor * sqrtdT0Ta * zSampleValue[i]);
            paths[i] = (ztDelta - 1.0) / dTaTb;
        }

        return paths;

    }

    public static double[] getRtPathSpotMeasure(ql.HullWhite model,
        double t,
        ulong seed,
        int numberOfPaths)
    {
        // path vector
        double[] paths = new double[numberOfPaths];
        
        // parameters
        double sigma = model.sigma();
        double k = model.a();
        double betaT = HullWhite.beta(0.0, t, k);

        // y_t and f0t
        double f0T = model.termStructure().currentLink().forwardRate(t, t, ql.Compounding.Simple, ql.Frequency.NoFrequency).rate();
        double yT = (sigma * sigma) / (2.0 * k) * (betaT - Math.Exp(-k * t) * betaT);
        
        // volatility
        double sigmaT = Math.Sqrt(HullWhite.beta(0.0, t, 2.0 * k));
        
        // path generator
        var rsg = (ql.InverseCumulativeRsg<ql.RandomSequenceGenerator<ql.MersenneTwisterUniformRng>, 
            ql.InverseCumulativeNormal>) new ql.PseudoRandom().make_sequence_generator(numberOfPaths, seed);
        
        var sample = rsg.nextSequence();
        var zSampleValue = sample.value;
        

        for (int i = 0; i < numberOfPaths; i++)
        {
            paths[i] = f0T + yT + sigmaT * zSampleValue[i];
        }

        return paths;
    }

    public static double[,] getRtBtPathSpotMeasure(ql.HullWhite model,
        double tR,
        double tBeta,
        ulong seed,
        int numberOfPaths)
    {
        
        // path vector
        double[,] paths = new double[numberOfPaths,2];
        
        // parameters
        double sigma = model.sigma();
        double k = model.a();
        
        // f0t
        double f0T = model.termStructure().currentLink()
            .forwardRate(tR, tR, ql.Compounding.Continuous, ql.Frequency.NoFrequency, true).rate();
        
        double df = model.termStructure().currentLink().discount(tR);
        double dfTbeta = model.termStructure().currentLink().discount(tBeta);
        var rTr = -Math.Log(df);
        var rTbeta = -Math.Log(dfTbeta);
        
        
        // Moments to simulate
        var varRt = HullWhite.getVarianceRt(tR, k, sigma);
        var varIt = HullWhite.getVarianceIt(tBeta, k, sigma);
        var covRtIt = HullWhite.getCovarianceRtIt(tR, tBeta, k, sigma);
        var driftRt = HullWhite.hjmAdjustment(tR, k, sigma);
        var driftIt = HullWhite.getDriftIt(tBeta, k, sigma);
        var rho = covRtIt / Math.Sqrt(varIt * varRt);
        var w1 = rho;
        var w2 = Math.Sqrt(1.0 - rho * rho);
        
        // path generator
        var rsg = (ql.InverseCumulativeRsg<ql.RandomSequenceGenerator<ql.MersenneTwisterUniformRng>, 
            ql.InverseCumulativeNormal>) new ql.PseudoRandom().make_sequence_generator(2, seed);


        for (int i = 0; i < numberOfPaths; i++)
        {
            
            var sample = rsg.nextSequence();
            var zSampleValue = sample.value;
            
            // sampling r_t
            paths[i, 0] = driftRt + Math.Sqrt(varRt) * zSampleValue[0] + f0T;
                    
            // Sampling It we have used Cholesky between r_t  and I_t
            var noiseI01 = w1 * zSampleValue[0] + w2 * zSampleValue[1];
            var i0Tp = driftIt + noiseI01 * Math.Sqrt(varIt);
           
           
            paths[i, 1] = Math.Exp(i0Tp + rTbeta);
        }

        return paths;

    }

    public static double[,] getIntRtPathSpotMeasure(ql.HullWhite model,
        double t1,
        double t0,
        int seed,
        int numberOfPaths)
    {
        
        // path vector
        double[,] paths = new double[numberOfPaths, 2];
        
        // parameters
        double sigma = model.sigma();
        double k = model.a();
        double d01 = t1 - t0;
        
        // random numbers with math.net
        var rng = new MersenneTwister(seed);
        double[] zSampleValue = new double[numberOfPaths];
        Normal.Samples(rng, zSampleValue, 0.0, 1.0);
        
        
        // path generator
        // var rsg = (ql.InverseCumulativeRsg<ql.RandomSequenceGenerator<ql.MersenneTwisterUniformRng>, 
        //     ql.InverseCumulativeNormal>) new ql.PseudoRandom().make_sequence_generator(numberOfPaths, seed);
        
        // var zSampleValue = rsg.nextSequence();
        var mu = HullWhite.getExpectedIntegralRt(k, sigma, t0, t1);
        var std = Math.Sqrt(HullWhite.getVarianceIntegralRt(k, sigma, t0, t1));
        
        double dfT0 = model.termStructure().currentLink().discount(t0, true);
        double dfT1 = model.termStructure().currentLink().discount(t1, true);
        double r0T = -Math.Log(dfT1 / dfT0);

        double meanAux = 0.0;
        double meanZs = 0.0;
        for (int i = 0; i < numberOfPaths; i++)
        {
            var it0Tot1 = mu + std * zSampleValue[i] + r0T;
            paths[i,0] = (Math.Exp(it0Tot1) - 1.0) / d01;
            paths[i,1] = it0Tot1 / d01;
            meanAux += paths[i, 1] / numberOfPaths;
            meanZs += zSampleValue[i] / numberOfPaths;

        }

        return paths;
    }

}