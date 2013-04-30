using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.LinearAlgebra.Generic;

namespace Filtering
{
    public class KF
    {
        SSF _model;
        Matrix<double> _p;
        DenseMatrix _mCov;
        DenseMatrix _pCov;
        DenseMatrix _e;
        DenseMatrix _at;
        DenseMatrix _ct;
        public KF(SSF model, DenseMatrix measCov, DenseMatrix procCov) {
            _model = model;
            _p = new DenseMatrix(_model.A.ColumnCount);
            _mCov = measCov;
            _pCov = procCov;
            _e = DenseMatrix.Identity(_model.A.ColumnCount);

            _at = new DenseMatrix(_model.A.Transpose().ToArray());
            _ct = new DenseMatrix(_model.C.Transpose().ToArray());
        }

        public Matrix<double> Step(DenseMatrix measurement, DenseMatrix input){
            var pProp = _model.A * _p * _at + _pCov;
            var k = pProp * _ct * (_model.C * pProp * _ct + _mCov).Inverse();
            var ekc = _e - k * _model.C;
            _p = ekc * pProp;

	        var xProp = _model.A * _model.State; // +_model.B * input;
            var x = xProp + k * (measurement - _model.C * xProp);
            x.CopyTo(_model.State);
            return _model.C * x;
        }
    }
}
