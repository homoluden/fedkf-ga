using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics.Distributions;

namespace GA.Helpers
{
    using System.Globalization;

    public static class SpecimenHelper
    {
        public static ContinuousUniform GeneRandom { get; set; }
        public static ContinuousUniform UniteRandom { get; private set; }
        public static Random RndInteger { get; private set; }

        public static double SimilarityThreshold = 1e1;

        static SpecimenHelper() {
            GeneRandom = new ContinuousUniform(Specimen.MinGeneValue, Specimen.MinGeneValue);
            UniteRandom = new ContinuousUniform(0.0, 1.0);
            RndInteger = new Random();
        }

        public static void GenerateGenes(ref Specimen spec)
        {
            if (spec.Length <= 0)
            {
                throw new ArgumentException("Specimen Genes Length must be greater that zero!");
            }

            double min = Specimen.MinGeneValue;
            double max = Specimen.MaxGeneValue;
            
            var genes = new double[spec.Length];

            for (int i = 0; i < spec.Length; i++)
            {
                genes[i] = GeneRandom.Sample();
            }

            spec.Genes = genes;
        }

        public static void SetGenes(ref Specimen spec, double[] genes) {
            if (genes == null || genes.Length == 0)
            {
                throw new ArgumentException("Genes array cannot be empty!");
            }

            spec.Genes = new double[genes.Length];
            genes.CopyTo(spec.Genes, 0);
            spec.Length = genes.Length;
        }


        public static List<Specimen> Crossover(this Specimen current, Specimen other)
        {
            if (current.Length != other.Length)
            {
                throw new ArgumentException("Specimens must have the same length of the genome!");
            }
            
            int pos = RndInteger.Next(current.Length - 3) + 2;

            var child1 = default(Specimen);
            child1.Fitness = double.NaN;
            
            var child2 = default(Specimen);
            child2.Fitness = double.NaN;
            
            var child3 = default(Specimen);
            child3.Fitness = double.NaN;

            var genes1 = new double[current.Length];
            var genes2 = new double[current.Length];
            var genes3 = new double[current.Length];

            for (int i = 0; i < current.Length; i++)
            {
                genes1[i] = (i < pos) ? current.Genes[i] : other.Genes[i];
                genes2[i] = (i >= pos) ? current.Genes[i] : other.Genes[i];
                genes3[i] = 0.5 * (current.Genes[i] + other.Genes[i]);
            }

            SetGenes(ref child1, genes1);
            SetGenes(ref child2, genes2);
            SetGenes(ref child3, genes3);

            return new List<Specimen> { child1, child2, child3 };
        }

        public static Specimen Mutate(this Specimen spec) 
        {
            if (spec.Genes == null || spec.Genes.Length == 0)
            {
                throw new ArgumentException("Genes array cannot be empty!");
            }

            var mutated = new Specimen();
            SetGenes(ref mutated, spec.Genes);

            for (int i = 0; i < spec.Genes.Length; i++)
            {
                if (UniteRandom.Sample() <= Specimen.MutationRate)
                {
                    mutated.Genes[i] = GeneRandom.Sample();
                }
            }
            
            return mutated;
        }

        public static string Print(this double[] array) {
            var builder = new StringBuilder();
            builder.Append("[");
            foreach (double val in array)
            {
                builder.Append(string.Format(" {0} ", val.ToString("0.0000e+00", CultureInfo.InvariantCulture)));
            }
            builder.Append(string.Format("]"));

            return builder.ToString();
        }
        
        public static string Print(this Specimen spec)
        {            
            var res = string.Format("{1} => Fitness: {0}", spec.Fitness, spec.Genes.Print());

            return res;
        }

        public static bool IsSimilar(this Specimen spec1, Specimen spec2) {
            var diff = 0.0;
            for (int i = 0; i < spec1.Length; i++)
            {
                diff += Math.Abs(spec1.Genes[i] - spec2.Genes[i]);
            }

            return diff <= SimilarityThreshold;
        }
    }
}
