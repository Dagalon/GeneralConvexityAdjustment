using ql = QLNet;
    
namespace ConvexityAdjustment_Lib
{
    public static class HullWhite
    {
        
        #region Convexity
        public static double convexityFraInArrears(ql.Handle<ql.YieldTermStructure> curve, double k, double sigma,
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
            double integral = beta(t1, 2.0 * t1, k) - beta(t2, t1 + t2, k);

            return m * integral;

        }

        public static double convexityFuture(ql.Handle<ql.YieldTermStructure> curve,  double k, double sigma, double t0, double t1, double t2,
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
            
            double m = beta(0.0, delta12, k) * sigma * sigma * Math.Exp(-k * t0);
            double integral = (1.0 - Math.Exp(-k * t0)) / (k * k) - t0 * Math.Exp(-k * t2) / k;
            double alpha = (h01 * df01) / (df00 * (h02 * df02) * delta12) * ( df02 * g02 - df01 * g01);

            return alpha * m * integral;
        }

        public static double convexityOis(ql.Handle<ql.YieldTermStructure> curve, double k, double sigma, double t1,
            double t2)
        {
            double b1 = beta(t1, t2, k);
            double b2 = beta(t1, t2, 2.0 * k);
            double delta12 = t2 - t1;

            double df1 = curve.link.discount(t1, true);
            double df2 = curve.link.discount(t2, true);
            double m = 0.5 * Math.Pow(sigma / k, 2.0);
            double expected = -Math.Log(df2 / df1) + m * (delta12 - 2.0 * b1 + b2);
            double i1 = m * ( t1 * Math.Exp(- 2.0 * k * t1) + Math.Exp( - 2.0 * k * t2) * t1  - 2.0 * Math.Exp( - k * (t1 + t2)) * t1);
            double i2 = m * (b2 + Math.Exp(-2.0 * k * t2) * delta12 - 2.0 * beta(t1+t2, 2.0*t2, k));

            return (Math.Exp(expected + i1 + i2) - 1.0) / delta12;
        }
        
        public static double convexityAvgOis(ql.Handle<ql.YieldTermStructure> curve, double k, double sigma, double t1,
            double t2)
        {
            double b1 = beta(t1, t2, k);
            double b2 = beta(t1, t2, 2.0 * k);
            double delta12 = t2 - t1;

            double df1 = curve.link.discount(t1, true);
            double df2 = curve.link.discount(t2, true);
            double m = 0.5 * Math.Pow(sigma / k, 2.0);
            double i1 = m * ( t1 * Math.Exp(- 2.0 * k * t1) + Math.Exp( - 2.0 * k * t2) * t1  - 2.0 * Math.Exp( - k * (t1 + t2)) * t1);
            double i2 = m * (b2 + Math.Exp(-2.0 * k * t2) * delta12 - 2.0 * beta(t1+t2, 2.0*t2, k));

            var r01 = convexityOis(curve, k, sigma, t1, t2);
            var logR01 = Math.Log(1.0 + delta12 * r01) / delta12;
            
            return logR01 - (i1 + i2);
        }

        public static double convexityCms(ql.Handle<ql.YieldTermStructure> discountCurve,  ql.Date valueDate, ql.Date ta, ql.Date tb, ql.Date tp, double annuity,
            double partialAnnuity, double partialVanillaSwap, double partialOisSwap, double partialOisSwapT0 , ql.DayCounter dc, double k, double sigma)
        {
            var dta = dc.yearFraction(valueDate, ta);
            var dtp = dc.yearFraction(valueDate, tp);
            
            var dFtp = discountCurve.link.discount(tp);
            var dFta = discountCurve.link.discount(ta);
            var m = - dFtp  * (beta(dta, dtp, k) + partialAnnuity / annuity);

            var alpha = sigma * sigma * beta(0.0, dta, 2.0 * k);
            return m * alpha * partialVanillaSwap;

        }

        public static double ConvexityCmsNewApproach(ql.Handle<ql.YieldTermStructure> discountCurve,
            ql.Date valueDate,
            ql.Date ta,
            ql.Date tb,
            ql.Date tp,
            double c,
            double annuityOis,
            double partialAnnuityOis,
            ql.DayCounter dc,
            double k,
            double sigma)
        {
            
            // anuuity and partial_x annuity for T_a
            
            var dta = dc.yearFraction(valueDate, ta);
            var dtb = dc.yearFraction(valueDate, tb);
            var dtp = dc.yearFraction(valueDate, tp);
            
            var dFtp = discountCurve.link.discount(tp);
            var dFta = discountCurve.link.discount(ta);
            var dFtb = discountCurve.link.discount(tb);
            var dFTaTp = dFtp / dFta;
            var dFtaTb = dFtb / dFta; 
            
            // Hull-White's convexity adjustment parameter
            var alpha = sigma * sigma * beta(0.0, dta, 2.0 * k);
            var m = dFTaTp / annuityOis;
            // var partialM = - m  * (partialAnnuityOis / annuityOis  + beta(dta, dtp, k));
            var swapOisRate = (dFta - dFtb) / annuityOis;
            var partialSwapOis = (beta(0.0, dtb, k) * dFtb - beta(0.0, dta, k) * dFta) /annuityOis - swapOisRate * partialAnnuityOis / annuityOis;
            var integral = sigma * sigma * (beta(dta, 2.0 * dta, k) - 0.5 * beta(dta, 3.0 * dta, k)) / k;

            // return beta(dta, dtp, k) * c * (alpha * (partialSwapOis / dFta) + (swapOisRate / dFta) * integral);
            return beta(dta, dtp, k) * alpha * c * (partialSwapOis / dFta);

        }

        #endregion
        
        #region Process tools

        public static double LiborVariance(double t, double t1, double t2, double sigma, double k)
        {
            double gamma1 = GetGammaVariance(t1, t2, t2, k, k);
            double gamma2 = GetGammaVariance(t1, t1, t1, k, k);
            double gamma12 = GetGammaVariance(t1, t1, t2, k, k);

            return sigma * sigma * (gamma1 + gamma2 - 2.0 * gamma12);
        }

        #endregion

        #region Tools

        public static double GetGammaVariance(double t, double t1, double t2, double k1, double k2)
        {
            double m = 1.0 / (k1 * k2);
            double part1 = t - beta(t1 - t, t1, k1) - beta(t2 - t, t2, k2);
            double part2 = t2 > t1
                ? Math.Exp(-k2 * (t2 - t1)) * beta(t1 - t, t1, k1 + k2)
                : Math.Exp(-k1 * (t1 - t2)) * beta(t2 - t, t2, k1 + k2);

            return m * (part1 + part2);
        }

      
        public static double getExpectedIntegralRt(double k, double sigma, double t0, double t1)
        {

            var m = 0.5 * Math.Pow(sigma / k, 2.0);
            var delta = t1 - t0;
            return m * (delta - 2.0 * beta(t0, t1, k) + beta(t0, t1, 2.0 * k));

        }
        
        public static double getVarianceIntegralRt(double k, double sigma, double t0, double t1)
        {

            double m = Math.Pow(sigma / k, 2.0);
            double delta = t1 - t0;
            double i1 = m * (beta(0.0, t0, 2.0 * k) + beta(delta, t1, 2.0 * k) - beta(t1 - t0, t1 + t0, k));
            double i2 = m * (delta + beta(0.0,delta, 2.0 * k) - 2.0 * beta(0.0, delta, k));

            return i1 + i2;
        }

        public static double GetVarianceIt(double t, double k, double sigma)
        {
            double m = Math.Pow(sigma / k, 2.0);
            return m * (beta(0, t, 2.0 * k) + Math.Exp(-2.0 * k * t) * t -
                                    2.0 * Math.Exp(-k * t) * beta(0.0, t, k));
        }


        public static double beta(double t0, double t1, double alpha)
        {
            double delta = t1 - t0;
            if (ql.Utils.close(alpha * delta, 0.0))
            {
                return Math.Exp(-alpha * t0) * delta;
            }

            return (Math.Exp(-alpha * t0) - Math.Exp(-alpha * t1)) / alpha;
        }

        public static double hjmAdjustment(double t, double k, double sigma)
        {
            double m = 0.5 * sigma * sigma / k;
            return m * (beta(0, t, k) - beta(t, 2.0 * t, k));
        }

        public static double forwardMeasureAdjustment(double dtP, double k, double sigma)
        {
            return (sigma * sigma / k) * (dtP * Math.Exp(-k * dtP) - beta(dtP, 2.0 * dtP, k));
            // return (sigma * sigma / k) * ( - dtP * Math.Exp(-k * dtP) + beta(0.0, dtP, k));
        }

        public static double GetDriftIt(double t, double k, double sigma)
        {
            var m = Math.Pow(sigma / k, 2.0);
            var expKt = Math.Exp(-k * t);
            var betaT = beta(0.0, t, k);

            return m * (betaT - expKt * t - beta(0.0, t, 2.0 * k) + expKt * betaT);
        }

        public static double GetCovarianceRtIt(double t, double k, double sigma)
        {
            var m = sigma * sigma / k;
            return m * (Math.Exp(-k * t) * t - beta(t, 2.0 * t, k));
        }

        public static double GetVarianceRt(double t, double k, double sigma)
        {
            return sigma * sigma * beta(0.0, t, 2.0 * k);
        }

        #endregion

    }
}