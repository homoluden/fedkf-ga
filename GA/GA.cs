using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;

using GA.Helpers;

namespace GA
{
    using System.Diagnostics.CodeAnalysis;

    using MathNet.Numerics.Distributions;

	/// <summary>
    /// Genetic Algorithm class
    /// </summary>
    public class Ga
    {
        #region Fields

        private static readonly Random Rnd = new Random();
        
        private List<double> _fitnessTable = new List<double>();

        private List<Specimen> _currGeneration = new List<Specimen>();
        
        #endregion // Fields


        #region Properties

        public Func<Specimen, double> FitnessFunction { get; set; }

        public int PopulationSize { get; set; }

        public int GenerationsCount { get; set; }

        public int GenomeSize { get; set; }

        public double CrossoverRate { get; set; }

        public double MutationRate { get; set; }

        public double TotalFitness { get; set; }

        /// <summary>
        /// Keep previous generation's fittest specimen
        /// </summary>
        public bool Elitism { get; set; }

        #endregion


        #region Ctors

        /// <summary>
        /// Default constructor sets mutation rate to 5%, crossover to 80%, population to 100,
        /// and generations to 2000.
        /// </summary>
        public Ga()
        {
            InitialValues();
            this.MutationRate = 0.05;
            this.CrossoverRate = 0.80;
            this.PopulationSize = 100;
            this.GenerationsCount = 2000;
        }

        public Ga(double crossoverRate, double mutationRate, int populationSize, int generationSize, int genomeSize)
        {
            InitialValues();
            this.MutationRate = mutationRate;
            this.CrossoverRate = crossoverRate;
            this.PopulationSize = populationSize;
            this.GenerationsCount = generationSize;
            this.GenomeSize = genomeSize;
        }

        public Ga(int genomeSize)
        {
            InitialValues();
            this.GenomeSize = genomeSize;
        }

        #endregion // Ctors


        #region Public Methods

        public void InitialValues()
        {
            this.Elitism = false;
        }
        
        /// <summary>
        /// Method which starts the GA executing.
        /// </summary>
        public void Go(double maxGeneVal, double minGeneVal)
        {
            if (FitnessFunction == null)
            {
                throw new ArgumentNullException("Need to supply fitness function");
            }

            if (this.GenomeSize == 0)
            {
                throw new IndexOutOfRangeException("Genome size not set");
            }

            _currGeneration = new List<Specimen>();
            Specimen.MutationRate = this.MutationRate;
            Specimen.MaxGeneValue = maxGeneVal;
            Specimen.MinGeneValue = minGeneVal;
			SpecimenHelper.GeneRandom = new ContinuousUniform(minGeneVal, maxGeneVal);
            
            Initiation();

            for (int i = 0; i < this.GenerationsCount; i++)
            {
                this.Selection();                
                
                Console.WriteLine("=== Selection {0} END ===\n\n", i);                
            }
            
            // TODO: Add writing to logger, data storage and/or console
        }

        #endregion

        /// <summary>
        /// After ranking all the genomes by fitness, use a 'roulette wheel' selection
        /// method.  This allocates a large probability of selection to those with the 
        /// highest fitness.
        /// </summary>
        /// <returns>Random individual biased towards highest fitness</returns>
        private int RouletteSelection()
        {
            double randomFitness = Rnd.NextDouble() * TotalFitness;
            int idx = -1;
            int first = 0;
            int last = this.PopulationSize - 1;
            int mid = (last - first) / 2;

            while (idx == -1 && first <= last)
            {
                // FIXME: build _fitnessTable with 
                if (randomFitness < _fitnessTable[mid])
                {
                    last = mid;
                }
                    // FIXME: Consider to use _currentGeneration instead of separate fitness table
                else if (randomFitness > _fitnessTable[mid])
                {
                    first = mid;
                }
                else if (randomFitness == _fitnessTable[mid])
                {
                    return mid;
                }

                mid = (first + last) / 2;

                // lies between i and i+1
                if ((last - first) == 1)
                {
                    idx = last;
                }
            }

            return idx;
        }

        /// <summary>
        /// Create the *initial* genomes by repeated calling the supplied fitness function
        /// </summary>
        [SuppressMessage("StyleCop.CSharp.ReadabilityRules", "SA1101:PrefixLocalCallsWithThis", Justification = "Reviewed. Suppression is OK here.")]
        private void Initiation()
        {
            Console.WriteLine("\n=== Starting Initiation ===");
            Console.WriteLine("Population Size: {0};\nGenome Size: {1};\nGene Value Range: [{2}; {3}]\n\n",
                                this.PopulationSize, this.GenomeSize, Specimen.MinGeneValue, Specimen.MaxGeneValue);

            _currGeneration = new List<Specimen>();

            var newSpecies = Enumerable.Range(0, PopulationSize).AsParallel().Select(i =>
            {
                var newSpec = new Specimen
                {
                    Length = this.GenomeSize
                };
                SpecimenHelper.GenerateGenes(ref newSpec);
                var fitness = FitnessFunction(newSpec);

                // FIXME: Replace '1e5' magic number with well described constant
                newSpec.Fitness = double.IsNaN(fitness) ? 0 : (double.IsInfinity(fitness) ? 1e5 : fitness);

                Console.WriteLine("Specimen {1} has Fitness: {0}", newSpec.Fitness, i);

                return newSpec;
            }).OrderBy(s => s.Fitness).Take(this.PopulationSize);

            _currGeneration = newSpecies.ToList(); // Huge load starts here :)

            _fitnessTable = new List<double>();
            foreach (var spec in _currGeneration)
            {
                if (!_fitnessTable.Any())
                {
                    _fitnessTable.Add(spec.Fitness);
                }
                else
                {
                    _fitnessTable.Add(_fitnessTable.Last() + spec.Fitness);
                }                
            }

            TotalFitness = _currGeneration.Sum(spec => spec.Fitness);

            var best = _currGeneration.Last();

            Console.WriteLine("=== Initiation Result ===\n");
            Console.WriteLine("Best Specimen:\n{0}", best.Print());
        }

        private void Selection()
        {            
            var tempGenerationContainer = new ConcurrentBag<Specimen>();
            if (this.Elitism)
            {
                var elite = _currGeneration.Last();
                tempGenerationContainer.Add(elite);
            }
            
            for (int i = 0; i < this.PopulationSize / 2.5; i++)
            {
                int pidx1 = this.PopulationSize - i - 1;
                int pidx2 = pidx1;
                while (pidx1 == pidx2 || _currGeneration[pidx1].IsSimilar(_currGeneration[pidx2]))
                {
                    pidx2 = RouletteSelection();
                }

                var parent1 = _currGeneration[pidx1].Mutate();
                var parent2 = _currGeneration[pidx2].Mutate();

                //Console.WriteLine("Selected Species {0} and {1}", pidx1, pidx2);

                var children = Rnd.NextDouble() < this.CrossoverRate ? parent1.Crossover(parent2) : new List<Specimen> { _currGeneration[pidx1], _currGeneration[pidx2] };

                foreach (var ch in children.AsParallel())
                {
                    if (double.IsNaN(ch.Fitness))
                    {
                        var fitness = FitnessFunction(ch);
                        var newChild = new Specimen
                            {
                                Genes = ch.Genes,
                                Length = ch.Length,
                                Fitness = double.IsNaN(fitness) ? 0 : (double.IsInfinity(fitness) ? 1e5 : fitness)
                            };
                        tempGenerationContainer.Add(newChild);
                    }
                    else
                    {
                        tempGenerationContainer.Add(ch);
                    }
                }
            }

            _currGeneration = tempGenerationContainer.OrderByDescending(s => s.Fitness).Take(this.PopulationSize).Reverse().ToList();

            _fitnessTable = new List<double>();
            foreach (var spec in _currGeneration)
            {
                if (!_fitnessTable.Any())
                {
                    _fitnessTable.Add(spec.Fitness);
                }
                else
                {
                    _fitnessTable.Add(_fitnessTable.Last() + spec.Fitness);
                }
            }

            TotalFitness = _currGeneration.Sum(spec => spec.Fitness);

            Console.WriteLine("=== Selection Result ===\n\n");
            Console.WriteLine("\n--- Top 5 ---");

            var best = _currGeneration.Last();

            Console.WriteLine("Best Specimen:\n{0}\n", best.Print());

            int j = 1;
            foreach (var spec in _currGeneration.AsEnumerable().Reverse().Take(5).Skip(1))
            {
                Console.WriteLine("Specimen {0} has Fitness: {1}\n", j++, spec.Fitness); // FIXME: Printout Genes
            }

            Console.WriteLine("Average Fitness: {0}\n", TotalFitness/PopulationSize);            
        }
    }
}
