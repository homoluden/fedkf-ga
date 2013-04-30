using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using FileHelpers;

namespace DAL.Records
{
    [DelimitedRecord(";")]
    public class Vector4
    {
        public double W { get; set; }

        public double X { get; set; }

        public double Y { get; set; }
        
        public double Z { get; set; }
    }
}
