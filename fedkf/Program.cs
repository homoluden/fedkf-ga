using System;
using System.Collections.Generic;
using System.Configuration;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Threading.Tasks;

using DAL;

using GA;
using GA.Helpers;

using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.LinearAlgebra.Generic;
using MathNet.Numerics.Statistics;

using Simulation;

namespace fedkf
{
    class Program
    {
        static void Main(string[] args)
        {
            if (args.Length > 0)
            {
                var cmd = args[0].ToLower(CultureInfo.InvariantCulture);
                switch (cmd)
                {
                    case "simulate": 
                    case "simulation":
                        InitializeSimulator();

                        FedKfSim.PrintSimResults = true;

                        var spec = new Specimen();
                        SpecimenHelper.SetGenes(ref spec, ReadSimulationGenes());

                        FedKfSim.Simulate(spec);
                        break;

                    case "set":
                        var settingName = args[1];
                        var settingValue = args[2];
						var config = ConfigurationManager.OpenExeConfiguration(ConfigurationUserLevel.None);
						config.AppSettings.Settings[settingName].Value = settingValue;
						config.Save(ConfigurationSaveMode.Modified);
						ConfigurationManager.RefreshSection("appSettings");
						
						Console.WriteLine("'{0}' set to {1}", settingName, settingValue);
                        break;

                    case "print":
                        Console.WriteLine("Current settings:");
                        foreach (var name in ConfigurationManager.AppSettings.AllKeys.AsParallel())
                        {
                            var value = ConfigurationManager.AppSettings[name];
                            Console.WriteLine("'{0}' => {1}", name, value);
                        }
                        break;

					case "help":
					case "?":
					case "-h":
		                PrintHelp();
						break;
                    default:
						Console.Error.WriteLine(string.Format("\nARGUMENT ERROR\n'{0}' is unknown command!\n", cmd));
						PrintHelp();
						break;
                }
            }
            else
            {
                InitializeSimulator();

                var genCount = int.Parse(ConfigurationManager.AppSettings["GenerationsCount"]);
                var popSize = int.Parse(ConfigurationManager.AppSettings["PopulationSize"]);
                var crossOver = double.Parse(ConfigurationManager.AppSettings["CrossoverRate"], FileParser.NumberFormat);
                var mutRate = double.Parse(ConfigurationManager.AppSettings["MutationRate"], FileParser.NumberFormat);
                var maxGeneVal = double.Parse(ConfigurationManager.AppSettings["MaxGeneValue"], FileParser.NumberFormat);
                var minGeneVal = double.Parse(ConfigurationManager.AppSettings["MinGeneValue"], FileParser.NumberFormat);
                var genomeLength = int.Parse(ConfigurationManager.AppSettings["GenomeLength"]);

                SpecimenHelper.SimilarityThreshold = double.Parse(
                    ConfigurationManager.AppSettings["SimilarityThreshold"], FileParser.NumberFormat);

                var ga = new Ga(genomeLength)
                {
                    FitnessFunction = FedKfSim.Simulate,
                    Elitism = true,
                    GenerationsCount = genCount,
                    PopulationSize = popSize,
                    CrossoverRate = crossOver,
                    MutationRate = mutRate
                };

                FedKfSim.PrintSimResults = false;
                ga.Go(maxGeneVal, minGeneVal);
            }

            Console.ReadLine();
        }

		private static void PrintHelp()
		{
			Console.WriteLine("Available commands are:");
			Console.WriteLine("\n'simulate' - runs Federative Kalman Filter simulation one time with parameters from file specified by 'SimulationGenesPath' setting\n");
			Console.WriteLine("\n'set <SettingName> <SettingValue>' - sets the setting value\n");
			Console.WriteLine("\n'print' - prints all settings with values to console\n");
			Console.WriteLine("\n'help' - (also '?', '-h') prints this doc\n");
			Console.WriteLine("\n\nRun without any command line arguments to proceed Genetic Optimization of Federative Kalman Filter");
		}

        private static double[] ReadSimulationGenes()
        {
            var simGenesPath = ConfigurationManager.AppSettings["SimulationGenesPath"];

            var simGenesString = string.Empty;

            using (var file = File.OpenRead(simGenesPath))
            {
                using (var reader = new StreamReader(file))
                {
                    simGenesString = reader.ReadToEnd();
                }
            }

            return simGenesString.Split(" ,\t".ToCharArray(), StringSplitOptions.RemoveEmptyEntries).Select(strVal => double.Parse(strVal)).ToArray();
        }
        private static void InitializeSimulator()
        {
            ReadDCM();
            ReadPCov();
            ReadModelParams();
            ReadSignalsAndNoises();
            FedKfSim.MaxSimLength = int.Parse(ConfigurationManager.AppSettings["MaxSimLength"]);
        }

        private static void ReadDCM()
        {
            var dcmPath = ConfigurationManager.AppSettings["DCMFilePath"];
            
            FedKfSim.DirectCosinsMatrix = FileParser.ParseMatrixFromFile(dcmPath);
        }

        private static void ReadPCov()
        {
            var pCovPath = ConfigurationManager.AppSettings["ProcCovFilePath"];
            
            FedKfSim.ProcessCovariances = FileParser.ParseMatrixFromFile(pCovPath);
        }

        private static void ReadModelParams()
        {
            var numOrderStr = ConfigurationManager.AppSettings["NumOrder"];
            var denOrderStr = ConfigurationManager.AppSettings["DenOrder"];
            var sensorsCountStr = ConfigurationManager.AppSettings["SensorsCount"];

            ushort numOrder;
            ushort denOrder = 0; // To suppress warning about use of unassigned variable
            ushort sensCount = 0; // To suppress warning about use of unassigned variable
            var parseSuccess = ushort.TryParse(numOrderStr, NumberStyles.Number, FileParser.NumberFormat, out numOrder)
                            && ushort.TryParse(denOrderStr, NumberStyles.Number, FileParser.NumberFormat, out denOrder)
                            && ushort.TryParse(sensorsCountStr, NumberStyles.Number, FileParser.NumberFormat, out sensCount);

            if (!parseSuccess)
            {
                throw new FormatException("Model parameters in settings has invalid format");
            }

            FedKfSim.NumOrder = numOrder;
            FedKfSim.DenOrder = denOrder;
            FedKfSim.SensorsCount = sensCount;
        }

        private static void ReadSignalsAndNoises()
        {
            var noisesPath = ConfigurationManager.AppSettings["NoisesFilePath"];
            var signalsPath = ConfigurationManager.AppSettings["SignalsFilePath"];
            var targetsPath = ConfigurationManager.AppSettings["TargetsFilePath"];

            FedKfSim.Noises = new DenseMatrix(FileParser.Read4ColonFile(noisesPath));
            FedKfSim.Signals = new DenseMatrix(FileParser.Read4ColonFile(signalsPath));
            FedKfSim.Targets = new DenseMatrix(FileParser.Read3ColonFile(targetsPath));

            var measCov = new DenseMatrix(4);
            double c00 = 0, c01 = 0, c02 = 0, c03 = 0, c11 = 0, c12 = 0, c13 = 0, c22 = 0, c23 = 0, c33 = 0;

            Vector<double> v1 = new DenseVector(1);
            Vector<double> v2 = new DenseVector(1);
            Vector<double> v3 = new DenseVector(1);
            Vector<double> v4 = new DenseVector(1);
            var s1 = new DescriptiveStatistics(new double[1]);
            var s2 = new DescriptiveStatistics(new double[1]);
            var s3 = new DescriptiveStatistics(new double[1]);
            var s4 = new DescriptiveStatistics(new double[1]);

            var t00 = Task.Run(() =>
            {
                v1 = FedKfSim.Noises.Column(0);
                s1 = new DescriptiveStatistics(v1);
                c00 = s1.Variance;
            });

            var t11 = Task.Run(() =>
            {
                v2 = FedKfSim.Noises.Column(1);
                s2 = new DescriptiveStatistics(v2);
                c11 = s2.Variance;
            });

            var t22 = Task.Run(() =>
            {
                v3 = FedKfSim.Noises.Column(2);
                s3 = new DescriptiveStatistics(v3);
                c22 = s3.Variance;
            });

            var t33 = Task.Run(() =>
            {
                v4 = FedKfSim.Noises.Column(3);
                s4 = new DescriptiveStatistics(v4);
                c33 = s4.Variance;
            });
            Task.WaitAll(new[] { t00, t11, t22, t33 });

            var t01 = Task.Run(() => c01 = CalcVariance(v1, s1.Mean, v2, s2.Mean, FedKfSim.Noises.RowCount));
            var t02 = Task.Run(() => c02 = CalcVariance(v1, s1.Mean, v3, s3.Mean, FedKfSim.Noises.RowCount));
            var t03 = Task.Run(() => c03 = CalcVariance(v1, s1.Mean, v4, s4.Mean, FedKfSim.Noises.RowCount));

            var t12 = Task.Run(() => c12 = CalcVariance(v2, s2.Mean, v3, s3.Mean, FedKfSim.Noises.RowCount));
            var t13 = Task.Run(() => c13 = CalcVariance(v2, s2.Mean, v4, s4.Mean, FedKfSim.Noises.RowCount));

            var t23 = Task.Run(() => c23 = CalcVariance(v3, s3.Mean, v4, s4.Mean, FedKfSim.Noises.RowCount));

            Task.WaitAll(new[] { t01, t02, t03, t12, t13, t23 });

            measCov[0, 0] = c00; measCov[0, 1] = c01; measCov[0, 2] = c02; measCov[0, 3] = c03;
            measCov[1, 0] = c01; measCov[1, 1] = c11; measCov[1, 2] = c12; measCov[1, 3] = c13;
            measCov[2, 0] = c02; measCov[2, 1] = c12; measCov[2, 2] = c22; measCov[2, 3] = c23;
            measCov[3, 0] = c03; measCov[3, 1] = c13; measCov[3, 2] = c23; measCov[3, 3] = c33;

            FedKfSim.SensorsOutputCovariances = measCov;
        }

        private static double CalcVariance(IEnumerable<double> v1, double mean1, IEnumerable<double> v2, double mean2, int length)
        {
            var zipped = v1.Take(length).Zip(v2.Take(length), (i1, i2) => new[] { i1, i2 });

            var sum = zipped.AsParallel().Sum(z => (z[0] - mean1) * (z[1] - mean2));

            return sum / (length - 1);
        }
    }
}
