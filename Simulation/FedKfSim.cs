
namespace Simulation
{
    using System;
    using System.Collections.Concurrent;
    using System.Linq;
    using System.Threading.Tasks;

    using DAL;
    using DAL.Records;
    using Filtering;
    using Filtering.Helpers;
    using GA;
    using GA.Helpers;
    using MathNet.Numerics.LinearAlgebra.Double;

    public static class FedKfSim
    {
        public static ushort SensorsCount { get; set; }

        public static ushort NumOrder { get; set; }

        public static ushort DenOrder { get; set; }

        public static DenseMatrix SensorsOutputCovariances { get; set; }

        public static DenseMatrix ProcessCovariances { get; set; }

        public static DenseMatrix DirectCosinsMatrix { get; set; }

        public static DenseMatrix Noises { get; set; }

        public static DenseMatrix Signals { get; set; }

        public static DenseMatrix Targets { get; set; }

        public static int MaxSimLength { get; set; }

        public static bool PrintSimResults { get; set; }

        public static FedKF ToFedKf(this Specimen spec)
        {
            int totalLength = SensorsCount * (NumOrder + DenOrder);
            if (totalLength > spec.Genes.Length)
            {
                throw new ArgumentException("Genome of the specimen is too short;");
            }

            var filters = new KF[SensorsCount];
            var numDenLength = NumOrder + DenOrder;
            var partitioner = Partitioner.Create(0, totalLength, numDenLength);

            Parallel.ForEach(
                partitioner,
                (range, loopstate) =>
                {
                    int sensorIdx = range.Item1 / numDenLength;

                    var numParams = new double[NumOrder];
                    var denParams = new double[DenOrder];
                    for (int i = range.Item1, j = 0; i < range.Item2; i++, j++)
                    {
                        if (j < NumOrder)
                        {
                            numParams[j] = spec.Genes[i];
                        }
                        else
                        {
                            denParams[j - NumOrder] = spec.Genes[i];
                        }
                    }

                    filters[sensorIdx] = KFBuilder.BuildKf(numParams, denParams, new DenseMatrix(1, 1, SensorsOutputCovariances[sensorIdx, sensorIdx]), ProcessCovariances);
                });

            var fkf = new FedKF(filters, DirectCosinsMatrix, SensorsOutputCovariances);
            return fkf;
        }

        public static double Simulate(Specimen spec)
        {
            var fkf = spec.ToFedKf();
            var meas = new double[4];

            var n = DirectCosinsMatrix;
            //var nt = n.Transpose();
            //var c = SensorsOutputCovariances.Inverse();
            //var pseudoInverse = (nt * n).Inverse() * nt;
            
            var err = 0.0;
            int lng = Math.Min(Signals.RowCount, MaxSimLength);

            var results = new Vector3[lng];
            results[0] = new Vector3 { X = 0.0, Y = 0.0, Z = 0.0 };

            for (int i = 0; i < lng; i++)
            {
                var sigRow = Signals.Row(i);
                var noiseRow = Noises.Row(i);
                var targRow = Targets.Row(i);
                meas[0] = sigRow[0] + noiseRow[0];
                meas[1] = sigRow[1] + noiseRow[1];
                meas[2] = sigRow[2] + noiseRow[2];
                meas[3] = sigRow[3] + noiseRow[3];

                DenseMatrix inps;
                if (i > 0)
                {
                    inps = n * results[i - 1].ToMatrix();
                }
                else
                {
                    inps = new DenseMatrix(4, 1, 0.0);
                }
                var res = fkf.Step(meas, inps.ToColumnWiseArray());

                var errs = new double[] { res[0, 0] - targRow[0], res[1, 0] - targRow[1], res[2, 0] - targRow[2] };
                err += (errs[0] * errs[0]) + (errs[1] * errs[1]) + (errs[2] * errs[2]);
                results[i] = new Vector3 { X = res[0, 0], Y = res[1, 0], Z = res[2, 0] };

                if (PrintSimResults)
                {
                    Console.WriteLine(res.ToColumnWiseArray().Print());
                }

                if (double.IsNaN(err))
                {
                    return double.NaN;
                }
            }

            if (PrintSimResults)
            {
                FileParser.Write3ColonFile(@"Data\Evaluations1000.csv", results);
            }

            return 1/err*lng;
        }
    }
}
