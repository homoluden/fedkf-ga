using System;
using System.Globalization;
using System.IO;
using System.Linq;

namespace DAL
{
    using System.Collections.Concurrent;
    using System.Threading.Tasks;

    using DAL.Records;
    using FileHelpers;
    using MathNet.Numerics.LinearAlgebra.Double;

    public static class FileParser
    {
        private static readonly NumberFormatInfo _format;
        public static NumberFormatInfo NumberFormat
        {
            get
            {
                return _format;
            }
        }

        static FileParser()
        {
            _format = new NumberFormatInfo
            {
                NumberDecimalSeparator = ".",
                NegativeSign = "-"
            };
        }

        public static void Write3ColonFile(string path, Vector3[] data)
        {
            var engine = new FileHelperEngine(typeof(Vector3));
            
            Parallel.Invoke(() => engine.WriteFile(path, data));
        }

        public static double[,] Read3ColonFile(string path)
        {
            var engine = new FileHelperEngine(typeof(Vector3));
            var values = (Vector3[])engine.ReadFile(path);

            var result = new double[values.Length, 3];
            
            // Convert 10 items per thread
            var partitioner = Partitioner.Create(0, values.Length, 10);

            Parallel.ForEach(
                partitioner,
                (range, state) =>
                    {
                        for (int i = range.Item1; i < range.Item2; i++)
                        {
                            result[i, 0] = values[i].X;
                            result[i, 1] = values[i].Y;
                            result[i, 2] = values[i].Z;
                        }
                    });

            return result;
        }

        public static double[,] Read4ColonFile(string path)
        {
            var engine = new FileHelperEngine(typeof(Vector4));
            var values = (Vector4[])engine.ReadFile(path);

            var result = new double[values.Length, 4];
            
            // Convert 10 items per thread
            var partitioner = Partitioner.Create(0, values.Length, 10);

            Parallel.ForEach(
                partitioner,
                (range, state) =>
                    {
                        for (int i = range.Item1; i < range.Item2; i++)
                        {
                            result[i, 0] = values[i].W;
                            result[i, 1] = values[i].X;
                            result[i, 2] = values[i].Y;
                            result[i, 3] = values[i].Z;
                        }
                    });

            return result;
        }

        public static DenseMatrix ParseMatrixFromFile(string path)
        {
            string content;
            using (var fs = new FileStream(path, FileMode.Open, FileAccess.Read))
            {
                using (var reader = new StreamReader(fs))
                {
                    content = reader.ReadToEnd();
                }
            }

            var mtx = ParseMatlabMatrix(content);
            return mtx;
        }

        public static DenseMatrix ParseMatlabMatrix(string matrixString)
        {
            matrixString = matrixString.Trim("[]\r\n".ToCharArray());
            var rows = matrixString.Split(";".ToCharArray(), StringSplitOptions.RemoveEmptyEntries);
            var mtxStringData =
                rows.Select(row => row.Split(",\t \r\n".ToCharArray(), StringSplitOptions.RemoveEmptyEntries)).ToArray();

            var rowsCount = mtxStringData.Length;
            var colsCount = mtxStringData[0].Length;
            var mtxData = new double[rowsCount, colsCount];

            for (int i = 0; i < rowsCount; i++)
            {
                for (int j = 0; j < colsCount; j++)
                {
                    var success = double.TryParse(mtxStringData[i][j], NumberStyles.Number, _format, out mtxData[i, j]);
                    if (!success)
                    {
                        throw new FormatException("Wrong format of numbers string representation!");
                    }
                }
            }

            return new DenseMatrix(mtxData);
        }
    }
}
