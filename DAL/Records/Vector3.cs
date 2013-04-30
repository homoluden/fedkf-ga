using FileHelpers;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using MathNet.Numerics.LinearAlgebra.Double;

namespace DAL.Records
{
    [DelimitedRecord(";")]
    public class Vector3
    {
        public double X { get; set; }

        public double Y { get; set; }

        public double Z { get; set; }

        public DenseMatrix ToMatrix()
        {
            var res = new DenseMatrix(3, 1);
            res[0, 0] = this.X;
            res[1, 0] = this.Y;
            res[2, 0] = this.Z;
            return res;
        }
    }
}
