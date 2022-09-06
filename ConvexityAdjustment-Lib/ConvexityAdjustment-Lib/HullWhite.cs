using ql = QLNet;

namespace ConvexityAdjustment_Lib
{
    public static class HullWhite
    {
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
            
            double m = sigma * sigma * Math.Exp(-k * t0);
            double integral = (1.0 - Math.Exp(-k * t0)) / (k * k) - t0 * Math.Exp(-k * t2) / k;
            double alpha = (h01 * df01) / (df00 * (h02 * df02) * delta12) * ( df02 * g02 - df01 * g01);

            return alpha * m * integral;


        }


    }
}