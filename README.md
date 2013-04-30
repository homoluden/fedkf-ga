FEDKF-GA
========

Federated Kalman Filter implementation with Genetic Algorithm for FKF parameters optimization.
Implemented Genetic Algorithm are based on "[A Simple C# Genetic Algorithm][1]" by Barry Lapthorn ([website][2]).

How to use
========
You can get basic usage info with "help" command:

    fedkf help
This command prints out basic documentation:

    Available commands are:
	'simulate' - runs Federative Kalman Filter simulation one time with parameters from file specified by 'SimulationGenesPath' setting
	'set <SettingName> <SettingValue>' - sets the setting value
	'print' - prints all settings with values to console
	'help' - (also '?', '-h') prints this doc
	Run without any command line arguments to proceed Genetic Optimization of Federative Kalman Filter
If you run program without any arguments, then Kalman Filter simulation process will be performed. Parameters for Kalman Filter will be taken from application config (**fedkf.exe.Config**):

    <?xml version="1.0" encoding="utf-8" ?>
	<configuration>
	  ...
	  <appSettings>
	    <add key="SensorsCount" value="4" />
	    <add key="NumOrder" value="2" />
	    <add key="DenOrder" value="4" />
	    <add key="DCMFilePath" value="Data\dcm.txt" />
	    <add key="NoisesFilePath" value="Data\Noises1000.csv" />
	    <add key="SignalsFilePath" value="Data\Signals1000.csv" />
	    <add key="TargetsFilePath" value="Data\Targets1000.csv" />
	    <add key="ProcCovFilePath" value="Data\pcov.txt" />
	    <add key="MaxSimLength" value="1000" />
	    <add key="GenerationsCount" value="1000" />
	    <add key="PopulationSize" value="300" />
	    <add key="CrossoverRate" value="0.8" />
	    <add key="MutationRate" value="0.09" />
	    <add key="MaxGeneValue" value="10.0" />
	    <add key="MinGeneValue" value="0.0" />
	    <add key="GenomeLength" value="24" />
	    <add key="SimilarityThreshold" value="0.0000001" />
 	   <add key="SimulationGenesPath" value="Data\SimGenes.txt" />
	  </appSettings>
	</configuration>
and from files specified by paths in this config.

***FIXME: provide detailed info about each parameter and file***

TODOs
=====

 - Parameters validation
 - Evolution statistics export

  [1]: http://www.codeproject.com/Articles/3172/A-Simple-C-Genetic-Algorithm
  [2]: http://www.lapthorn.net