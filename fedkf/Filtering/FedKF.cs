using System.Linq;

using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.LinearAlgebra.Generic;


namespace Filtering
{
	using DAL.Records;

	public class FedKF
    {
        KF[] _filters;
        DenseMatrix _filteredSignals;
        DenseMatrix _dcm;
        Matrix<double> _dcmt;
        Matrix<double> _cInv;
        public FedKF(KF[] filters, DenseMatrix dcm, DenseMatrix covariances) {
            _filters = filters;
            _filteredSignals = new DenseMatrix(_filters.Length, 1);
            _dcm = dcm;
            _dcmt = dcm.Transpose();
            _cInv = covariances.Inverse();
        }

        public Matrix<double> Step(double[] measurements, double[] inputs) {
            var meas = measurements.AsParallel().Select(m => new DenseMatrix(1, 1, m)).ToList();
            var inps = inputs.AsParallel().Select(inp => new DenseMatrix(1, 1, inp)).ToList();

            int i = 0;
            foreach (var m in meas.ToArray())
            {
                _filteredSignals[i, 0] = _filters[i].Step(m, inps[i])[0, 0];
                i++;
            }


            return (_dcmt * _cInv * _dcm).Inverse() * _dcmt * _cInv * _filteredSignals;
        }
    }
}
