using MathNet.Numerics.LinearAlgebra.Double;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Filtering.Helpers
{
    public static class KFBuilder
    {
        public static bool ParseContiniousTFString(string numString, string denString, out double[] num, out double[] denom)
        {
            var format = new System.Globalization.NumberFormatInfo
            {
                NumberDecimalSeparator = ".",
                NegativeSign = "-"
            };

            if (string.IsNullOrWhiteSpace(numString) || string.IsNullOrWhiteSpace(denString))
            {
                num = new double[0];
                denom = new double[0];
                return false;
            }

            string numStr = numString.Trim("[]".ToCharArray());
            string denomStr = denString.Trim("[]".ToCharArray());
            string[] vals = numStr.Split(" \t\n".ToCharArray(), System.StringSplitOptions.RemoveEmptyEntries);

            num = new double[vals.Length];
            for (int i = 0; i < vals.Length; i++)
            {
                double val;
                if (!double.TryParse(vals[i], System.Globalization.NumberStyles.Float, format, out val))
                {
                    denom = new double[0];
                    return false;
                }
                num[i] = val;
            }

            vals = denomStr.Split(" \t\n".ToCharArray(), System.StringSplitOptions.RemoveEmptyEntries);
            denom = new double[vals.Length];

            for (int i = 0; i < vals.Length; i++)
            {
                double val;
                if (!double.TryParse(vals[i], System.Globalization.NumberStyles.Float, format, out val))
                {
                    return false;
                }

                denom[i] = val;
            }
            return true;
        }

        public static SSF GenerateSensorModel(string numString, string denString, double ts)
        {
            double[] num, denom;
            if (ParseContiniousTFString(numString, denString, out num, out denom))
            {
                return SSF.C2DSS(new DenseVector(num), new DenseVector(denom), 1.0, ts);                
            }
            throw new ArgumentException("Numerator and/or Denominator are in invalid format");
        }

        public static KF BuildKf(double[] num, double[] den, DenseMatrix measCov, DenseMatrix procCov)
        {
            var format = new System.Globalization.NumberFormatInfo
            {
                NumberDecimalSeparator = ".",
                NegativeSign = "-"
            };

            var numStr = new StringBuilder();
            foreach (var par in num)
            {
                numStr.Append(" " + par.ToString(format));
            }
            var denStr = new StringBuilder();
            foreach (var par in den)
            {
                denStr.Append(" " + par.ToString(format));
            }

            var ssf = GenerateSensorModel(numStr.ToString(), denStr.ToString(), 0.01);

            var kf = new KF(ssf, measCov, procCov);
            return kf;
        }
    }
}
