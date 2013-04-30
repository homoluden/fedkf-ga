using System;
using System.Globalization;

using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.LinearAlgebra.Generic;

namespace Filtering
{
    public class SSF
    {
        #region Fields

        private double _ts;
        private double _k;
        private DenseMatrix _a;
        private DenseMatrix _b;
        private DenseMatrix _c;
        private double _d;
        private bool _tsSet;
        private bool _kSet;
        private bool _aSet;
        private bool _bSet;
        private bool _cSet;
        private DenseMatrix _x;

        #endregion // Fields
    
        #region Properties

        public double TimeSample
        {
            get { return _ts; }
            set
            {
                _ts = value;
                _tsSet = value.CompareTo(double.Epsilon) >= 0;
            }
        }

        public double Gain
        {
            get { return _k; }
            set
            {
                _k = value;
                _kSet = !(double.IsInfinity(value) || double.IsNaN(value));
            }
        }

        public DenseMatrix A
        {
            get { return _a; }
            set
            {
                _a = value;
                _aSet = value != null;
            }
        }

        public DenseMatrix B
        {
            get { return _b; }
            set
            {
                _b = value;
                _bSet = value != null;
            }
        }

        public DenseMatrix C
        {
            get { return _c; }
            set
            {
                _c = value;
                _cSet = value != null;
            }
        }

        public double D
        {
            get { return _d; }
            set
            {
                _d = value;
            }
        }

        public bool Initialized
        {
            get { return _aSet && _bSet && _cSet && _tsSet && _kSet; }
        }

        public DenseMatrix State { get { return _x; } }

        #endregion // Properties
	
        #region Ctors

        public SSF()
        {
        }

        //def initialize(a,b,c,d,ts,k)
        //    @a,@b,@c,@d,@ts,@t,@k = a,b,c,d,ts,0.0,k
        //    @x = NMatrix.float(@a.sizes[0],1) #Matrix.rows(Array.new(@a.row_size,[0.0]))
        //    @order = @a.sizes[1]#@a.row_size
        //end
        public SSF(DenseMatrix f, DenseMatrix r, DenseMatrix m, double d, double ts, double k)
        {
            _ts = ts;
            _k = k;
            _c = m;
            _b = r;
            _a = f;
            _d = d;

            _x = new DenseMatrix(_a.RowCount, 1);
        }

        #endregion // Ctors

        #region Public Methods
        //def step(u) # public double Step(double U, double xnoise, double fnoise)
        //    tx 	= @a*@x + @b*u
        //    f = @c*@x
        //    @x = tx.clone
        //        @t += @ts
        //    return (f*@k)[0,0]
        //end
        public double Step(double u)
        {
            var x = (_a * _x + _b * u);
            var y = (_c * _x)[0, 0];
            _x = x;
            return y * _k + _d * u;
        }

        public double Step(double u, out double x0)
        {
            var res = Step(u);
            x0 = _x[0, 0];
            return res;
        }

        
        //def SSF.align_num_denum(num,denum)
        //        d = denum.length - num.length
        //    if  d >= 0 then
        //        num = num.reverse
        //        d.downto(1) {
        //            num << 0.0
        //        }
        //        num = num.reverse
        //    else
        //        d.upto(-1) {
        //            denum << 0.0
        //        }
        //    end
        //    return num.length,num,denum
        //end
        public static int AlignNumDenom(ref Vector<double> numerator, ref Vector<double> denominator)
        {
            if (numerator == null || denominator == null) return 0;
            
            int dl = denominator.Count - numerator.Count;
            if (dl >= 0)
            {
                var newNum = new DenseVector(dl+numerator.Count, 0.0);
                for (int i = 0; i < numerator.Count; i++)
                {
                    newNum[dl + i] = numerator[i];
                }
                //numerator = newNum.Add(numerator);
                numerator = newNum;
            }
            else
            {
                denominator = denominator.Add(new DenseVector(-dl, 0.0));
            }
            return numerator.Count;
        }

        //def SSF.tf2ssd(num,denum,k,ts)# public static SSFilter TF2SSd(double[] Num, double[] Den, double k, double Ts)
        //    n,num,denum = SSF::align_num_denum(num,denum)
        //    a = []
        //    b = []
        //    c = []
        //    (0).upto(n-3) {
        //        |i|
        //        ta = NArray.float(n-1)
        //        ta[i+1] = 1.0
        //        a << ta.to_a
        //        b << [0.0]
        //    }
        //    ta = []
        //    (n-1).downto(1) {
        //        |i|
        //        ta << -denum[i]/denum[0]
        //    }
        //    a << ta
        //    b << [1.0/denum[0]]
        //    d = num[0]/denum[0]
		
        //    c = NMatrix.float(n-1,1)
        //    (n-1).downto(1) {
        //        |i|
        //        c[i-1,0] = d*(-denum[i]) + num[i]
        //    }
		
        //    am = NMatrix.rows(a)    #Matrix.rows(a)
        //    bm = NMatrix.rows(b)    #Matrix.rows(b)
        //    e =  NMatrix.float(n-1,n-1).unit       #Matrix.xdiag(n-1,1.0)
        //    f = e + am*ts + am**2*ts**2/2 + am**3*ts**3/3
        //    r = (e*ts + am*ts**2/2 + am**2*ts**3/3)*bm #e*ts*bm
        //    ssf = SSF.new(f,r,c,d,ts,k*denum.last/num.last)
        //    ssf
        //end
        public static SSF C2DSS(Vector<double> num, Vector<double> denom, double gain, double timeSample)
        {
            int n = AlignNumDenom(ref num, ref denom);
            var a = new DenseMatrix(n - 1);
            var b = new DenseMatrix(n - 1, 1, 0.0);
            for (int i = 0; i < n-2; i++)
            {
                a[i, i + 1] = 1.0;
            }
            b[n - 2, 0] = 1/denom[0];
            var d = num[0]/denom[0];
            for (int i = 0; i < n-1; i++)
            {
                a[n - 2, i] = -denom[n - 1 - i]/denom[0];
            }
            var c = new DenseMatrix(1, n - 1);
            for (int i = 0; i < n-1; i++)
            {
                c[0, i] = num[n - 1 - i] - denom[n - 1 - i] * d;
            }

            var e = new DenseMatrix(a.RowCount);
            for (int i = 0; i < a.RowCount; i++)
            {
                e[i,i] = 1.0;
            }

            var f = e + a.Multiply(timeSample) +
                        a.Multiply(a).Multiply(timeSample * timeSample * 0.5);// +
                        //a.Multiply(a).Multiply(a).Multiply(timeSample*timeSample*timeSample/3.0);
            var r = (e * timeSample + a.Multiply(timeSample * timeSample * 0.5)). // + a.Multiply(a).Multiply(Math.Pow(timeSample, 3)/3.0) ).
                        Multiply(b);
            return new SSF((DenseMatrix)f, (DenseMatrix)r, (DenseMatrix)c, d, timeSample, gain * denom[n - 1] / num[n - 1]);
        }

        public static SSF C2DSSStr(string numString, string denomString, double timeSample)
        {
            string[] numStrings = numString.Trim('[', ']').Split(" ,".ToCharArray(), StringSplitOptions.RemoveEmptyEntries);
            string[] denomStrings = denomString.Trim('[', ']').Split(" ,".ToCharArray(), StringSplitOptions.RemoveEmptyEntries);

            NumberFormatInfo format = new NumberFormatInfo
            {
                NumberDecimalSeparator = "."
            };

            var num = new double[numStrings.Length];
            var denom = new double[denomStrings.Length];

            int i = 0;
            foreach (var str in numStrings)
            {
                double val;
                num[i] = double.TryParse(str, System.Globalization.NumberStyles.Float, format, out val) ? val : 0.0;
                i++;
            }
            i = 0;
            foreach (var str in denomStrings)
            {
                double val;
                denom[i] = double.TryParse(str, System.Globalization.NumberStyles.Float, format, out val) ? val : 0.0;
                i++;
            }
            return C2DSS(new DenseVector(num), new DenseVector(denom), 1.0, timeSample);
        }

        #endregion // Public Methods
	
    }
}
