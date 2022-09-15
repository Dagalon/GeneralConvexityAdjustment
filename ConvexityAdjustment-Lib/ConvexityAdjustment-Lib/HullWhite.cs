using ql = QLNet;

namespace ConvexityAdjustment_Lib
{
    public static class HullWhite
    {

        public static double ConvexityFraInArrears(ql.Handle<ql.YieldTermStructure> curve, double k, double sigma,
            double t1, double t2,
            double basisSpread)
        {
            // df, year fracs, basis,....
            double df01 = curve.currentLink().discount(t1, true);
            double df02 = curve.currentLink().discount(t2, true);
            double delta12 = (t2 - t1);
            double h01 = Math.Exp(-basisSpread * t1);
            double h02 = Math.Exp(-basisSpread * t2);
            double df12 = (df02 * h02) / (df01 * h01);

            double m = sigma * sigma / (delta12 * df12 * k);
            double integral = Beta(t1, 2.0 * t1, k) - Beta(t2, t1 + t2, k);

            return m * integral;

        }

        public static double ConvexityFuture(ql.Handle<ql.YieldTermStructure> curve,  double k, double sigma, double t0, double t1, double t2,
            double basisSpread)
        {
            // df, year fracs, basis,....
            double df00 = curve.currentLink().discount(t0, true);
            double df01 = curve.currentLink().discount(t1, true);
            double df02 = curve.currentLink().discount(t2, true);
            double delta12 = (t2 - t1);
            double delta02 = (t2 - t0);
            double delta01 = (t1 - t0);
            double h01 = Math.Exp(-basisSpread * delta01);
            double h02 = Math.Exp(-basisSpread * delta02);
            
            // adjustment
            double g01 = (1.0 - Math.Exp(-k * delta01)) / k;
            double g02 = (1.0 - Math.Exp(-k * delta02)) / k;
            
            double m = Beta(0.0, delta12, k) * sigma * sigma * Math.Exp(-k * t0);
            double integral = (1.0 - Math.Exp(-k * t0)) / (k * k) - t0 * Math.Exp(-k * t2) / k;
            double alpha = (h01 * df01) / (df00 * (h02 * df02) * delta12) * ( df02 * g02 - df01 * g01);

            return alpha * m * integral;
        }

        
        public static double GetExpectedIntegralRt(double k, double sigma, double t0, double t1)
        {

            var alpha = 0.5 * sigma*sigma / k;
            var part1 = t1 - t0;
            var part2 = 0.5 * (Math.Exp(-2.0 * k * t0) - Math.Exp(-2.0 * k * t1)) / k;
            var part3 = (1.0 - Math.Exp(- k * (t1 - t0))) / k;
            var part4 = (Math.Exp(-k * (t0 + t1)) - Math.Exp(-2.0 * k * t1)) / k;

            return alpha * (part1 - part2 - part3 + part4);

        }
        
        public static double GetVarianceIntegralRt(double k, double sigma, double t0, double t1)
        {
            var alpha = sigma * sigma;
            var part1 = t1 - t0;
            var part2 = 0.5 * (1.0 - Math.Exp(-2.0 * k * (t1 - t0))) / k;
            var part3 = (1.0 - Math.Exp(-k * (t1 - t0))) / k;

            return alpha * (part1 + part2 - 2.0 * part3);
        }

        public static double Beta(double t0, double t1, double alpha)
        {
            double delta = t1 - t0;
            if (ql.Utils.close(alpha * delta, 0.0))
            {
                return Math.Exp(-alpha * t0) * delta;
            }

            return (Math.Exp(-alpha * t0) - Math.Exp(-alpha * t1)) / alpha;
        }

    }
}