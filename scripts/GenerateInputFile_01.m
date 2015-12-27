%% generate SCORE.Input.File.uds
clear;
clc;

id = 91008;
velo = 0.9; 
%1 - 0.0, 
%2 - 0.5, 
%5 - 0.8, 
%10 - 0.9
nuclei = 1000; %value is irrelevant for GBNucleation

%fdeformed = fopen('DefGrPool.CR90.GammaBCCMonoX.uds','r');
fdeformed = 'DefGrPool.CR90.Alpha12Gamma12OriSpread.uds';
DELIMITER = '\t';
HEADERLINES = 1; %12.uds
%HEADERLINES = 17; %MonoX.uds

%% dont change if necessary beyond this point here

shfname = ['score_phdmk.' num2str(id) '.sh'];
shid = fopen(shfname, 'w+');

fprintf( shid, '#!/usr/bin/env zsh\n');
fprintf( shid, '\n');
fprintf( shid, ['#BSUB -J SCORE' num2str(id) '\n']);
fprintf( shid, ['#BSUB -o SCORE' num2str(id) '.%%J']);
fprintf( shid, '\n');
fprintf( shid, '#BSUB -W 2:00\n');			
fprintf( shid, '#BSUB -M 50000\n');
fprintf( shid, '#BSUB -a openmpi\n');
fprintf( shid, '#BSUB -n 1\n');
%fprintf( shid, '### BSUB -m mpi-l\n');
%fprintf( shid, '### BSUB -R "select[hpcwork]"\n');
fprintf( shid, '#BSUB -x\n');
fprintf( shid, '#BSUB -P jara0076\n');
%fprintf( shid, '\n');
fprintf( shid, ['$MPIEXEC $FLAGS_MPI_BATCH score_intel14O2 SCORE.Input.' num2str(id) '.uds ' num2str(id) '\n']);
fclose(shid);

'Shellscript generated'

%% the corresponding input file
fname = ['SCORE.Input.' num2str(id) '.uds'];
fid = fopen(fname, 'w+');
fprintf( fid, '|| ParameterList        	|| *(ss)        || Parameter, Value\n');

fprintf( fid, '"CAEnsembleSize"                 "1"			// how many cellular automata should be devised\n');
fprintf( fid, '"XMAX"							"1.000"		// at which point should the integration scheme be stopped at the latest, define tmax in processing profile\n');
fprintf( fid, '"TIMEMAX"						"36000.0"	// at which point in time should every automaton stop the simulation\n');
fprintf( fid, '"NMAX"							"100000"	// after which number of integration steps should the scheme be stopped regardless the local X\n');
fprintf( fid, '"CellSize"						"1.0E+00"	// micron, edge length of one simulation voxel\n');
fprintf( fid, '"3DCAEdgeLengthInCellsRD"		"1200"		// how many cells in each automaton along each edge\n');
fprintf( fid, '"3DCAEdgeLengthInCellsND"		"1200"\n');
fprintf( fid, '"3DCAEdgeLengthInCellsTD"		"1200"\n');
fprintf( fid, '"MeanDefGrainSizeinRD"			"1.00E+03" 	// micron, average typical or median grain size utilized in the case of\n'); 
fprintf( fid, '"MeanDefGrainSizeinND"			"1.00E+01"\n');
fprintf( fid, '"MeanDefGrainSizeinTD"			"1.00E+02"\n');
fprintf( fid, '"DefStructureSynthesis"			"1"			// 1 - default is cuboid grains the GIA approach, 2 - poisson voronoi, 3 - import CPFEM data not yet correctly implemented\n');
fprintf( fid, '"MeanDefSizePoissonDiameter"     "24.8"		// micron, becomes interpreted in number of grains and seed points of a structure\n');
fprintf( fid, '\n');
fprintf( fid, ['"NucleationDensityLocalCSR"		"' num2str(nuclei) '"		// assuming each automaton hosting a complete spatial random point process how many points should be considered\n'] );
fprintf( fid, '"CSRNucleation"					"0"			// 1 - csr_enforced, 2 - csr_pickedrandomly, else no\n');
fprintf( fid, '"GBNucleation"					"2"			// 1 - yes all boundaries, 2 - yes only at high-angle grain boundaries, else no\n');
fprintf( fid, '"GBNucDensity2Number"			"2.0E-4"	// scales the unit area density as does SCALING_LAMBDA in pp3 analysis into a number of nuclei\n');
fprintf( fid, '"GBNucRhoDiff2Density"			"3.5436E-15"		// 10.0 over maximum dislocation density in the system\n');
fprintf( fid, '"GBNucMaximumScatter"			"5.0"		// how strongly are grain boundary nuclei at most deviating from the matrix orientation\n');
fprintf( fid, '"IncubationTimeModel"			"1"			// 1 site saturated as the default, 2 time dependent\n');
fprintf( fid, '"IncubationTimeScatter"			"0.0"		// sigma of a Rayleigh distribution according to which the nuclei become seeded\n');
fprintf( fid, '\n');
fprintf( fid, '"MobilityModel"					"2"			// 1 Sebald Gottstein is the default model, 2 Rollett Holm\n');
fprintf( fid, '\n');
fprintf( fid, '"RenderingOfMicrostructure"		"0"			// 0 - not at all, 2 - 2D zsection is implicit at half in z direction, 3 - 3D for all in Render Microstructure\n');
fprintf( fid, '"RenderingColorModel"			"2"			// 1 - disjoint grains with unsigned int IDs, 2 - IPFZ coloring\n');
fprintf( fid, '"RenderingBoundaries"			"0"			// 1 - identify which voxel is located at a boundaries and to which grain it belongs BE CAREFUL LARGE FILES\n');
fprintf( fid, '"OutputLogBoundaries"			"0"			// 1 - write detailed file containing all boundaries, else not\n');
fprintf( fid, '"OutputRXFrontStats"             "0" 		// 1 - yes, all CA output detailed loginformation on the interface cell list, else they dont\n');
fprintf( fid, '"OutputSingleGrainStats"         "0"			// 1 - yes, all CA output detailed information on the evolution of each grain,  else they dont\n');
fprintf( fid, '"OutputOfHPF"					"0"			// only for the automata of the MASTER the 3D and 2D HPF is computed\n');
fprintf( fid, '\n');
fprintf( fid, '\n');
fprintf( fid, '|| IntegratorAccuracy		|| *(ss)			|| Parameter, Value\n');
fprintf( fid, '"MAXFILLPERSTEP"                 "0.1"		// how much partial volume of a cell does the fastest boundary sweep in one timestep\n');
fprintf( fid, '"RediscretizationSteps"			"1000"		// in how many time step should the real time interval tstart to tend be rediscretized\n');
fprintf( fid, '\n');
fprintf( fid, '\n');
fprintf( fid, '|| PerformanceParameter		|| *(ss)			|| Parameter, Value\n');
fprintf( fid, '"InitialCacheSizeRXCells"		"0.10"		// elaborated guess about how many memory should be reallocated to store active cells, relative value is converted into absolute number of cells per local region\n');
fprintf( fid, '"ReCacheSizeRXCells"         	"1.0"		// elaborated guess what happens if in the local region more cells become necessary then the memory is expanded by a factor of this times already present memory\n');
fprintf( fid, '\n');
fprintf( fid, '|| MaterialProperties		|| *(ss)			|| Parameter, Value\n');
fprintf( fid, '"ZeroKelvinShearModulusG0"		"85.0e9"	// in Pascal\n');
fprintf( fid, '"FirstOrderdG0dT"				"0.0e9"	// internally overwritten to incorporate pronounced softening towards the transition temperature\n');
fprintf( fid, '"ZeroCelsiusBurgersVector"		"2.50e-10"\n');	
fprintf( fid, '"AlloyConstantThermalExpCoeff"	"1.0"		// in Willey model, C constants see Hatch JE, Aluminum, see also Touloukian\n');
fprintf( fid, '"FirstOrderThermalExpCoeff"		"0.0"\n');
fprintf( fid, '"SecondOrderThermalExpCoeff"     "0.0"\n');	
fprintf( fid, '"MeltingTemperature"             "800.0"		// in degrees celsius here chosen as the range of data accuracy\n');
fprintf( fid, '\n');
fprintf( fid, '"LAGBm0"                         "0.3"		// classical SebaldGottstein mobility model 3584m^4/Js preexponential factor intrinsic mobility of low-angle boundaries, default truncation at 15degrees\n');
fprintf( fid, '"LAGBHact"						"1.200"		// 1550eV activation enthalpy low-angle grain boundaries disori ij < 15degree eV\n');
fprintf( fid, '"HAGBm0"                         "3.0"		// alike for general high-angle grain boundaries\n');
fprintf( fid, '"HAGBHact"						"1.200"\n');
fprintf( fid, '"GSm0"							"3.0"		// boundaries close to 40deg111 in misorientation space\n');
fprintf( fid, '"GSHact"                         "1.200"\n');
fprintf( fid, '\n');
fprintf( fid, '"RHModelHAGBm0"					"300.0"		// the model often utilized by Rollett and Holm in which LAGB are practically immobile and a sharp sigmoidal transition of the maximum mobilities to HAGB occurs between 10 and 16 deg\n');
fprintf( fid, '"RHModelHAGBHact"				"2.45"\n');
fprintf( fid, ['"RHModelLAGBHAGBCut"             "' num2str(velo) '"		// one minus cut times exp minus trans multiplying the quotient disori over lagbhagb transition disorientation to the power of the exponent\n'] );
fprintf( fid, '"RHModelLAGBHAGBTrans"           "5.0"\n');
fprintf( fid, '"RHModelLAGBHAGBExponent"		"9.0"\n');
fprintf( fid, '\n');
fprintf( fid, '\n');
fprintf( fid, '|| RecoveryParameter		|| *(ss)	|| Parameter, Value\n');
fprintf( fid, '\n');
fprintf( fid, '\n');
fprintf( fid, '|| ZenerDragParameter		|| *(ss)	|| Parameter, Value\n');
fprintf( fid, '"ZenerConsider"					"0"			// 0 switched off, 1 constant, 2 timedependent in accord with fr evolution\n');
fprintf( fid, '"ZenerAlpha"                 	"1.5"		// 1.5gammaHAGB 0.324 J/m^2 from MeyersMurr Interfacial\n');
fprintf( fid, '"ZenerGamma"                     "0.324"		// f over r we take the first entry from EvoDragging Particles\n');
fprintf( fid, '\n');
fprintf( fid, '\n');
fprintf( fid, '|| IdealComponents			|| 	*(sffff)		|| Name, IdealcomponentBunge1, IdealcomponentBunge2, IdealcomponentBunge3, each component only counted one and to the closest all other are random, random is has ID zero\n');
fprintf( fid, '"Alpha001-110"			0.0			0.0			45.0		10.0\n');
fprintf( fid, '"Alpha112-110"			0.0			35.0		45.0		10.0\n');
fprintf( fid, '"Alpha111-110"			0.0			55.0		45.0		10.0\n');
fprintf( fid, '"Gamma111-112"			30.0		55.0		45.0		10.0\n');
fprintf( fid, '"Gamma111-110"			60.0		55.0		45.0		10.0\n');
fprintf( fid, '"Gamma111-112"			90.0		55.0		45.0		10.0\n');
%fprintf( fid, '"NucAlpha"				0.0			30.0		45.0		2.0\n');
%fprintf( fid, '"NucBeta"				55.0		55.0		45.0		2.0\n');
fprintf( fid, '\n');
fprintf( fid, '\n');
fprintf( fid, '|| RenderMicrostructure	||	*(sf)	|| Local recrystallized fraction, local value\n');
fprintf( fid, '"X"		0.25\n');
fprintf( fid, '"X"		0.50\n');
fprintf( fid, '"X"		0.75\n');
fprintf( fid, '"X"		1.00\n');
fprintf( fid, '\n');
fprintf( fid, '\n');
fprintf( fid, '|| RenderOnlyTheseRegions || *(i) || CAID indices start from 0 to n minus one\n');
fprintf( fid, '\n');
fprintf( fid, '\n');
fprintf( fid, '|| ProcessingSchedule		||	*(ff)			|| time, temperature linearized heating profiles\n');
fprintf( fid, '0.0								500.0\n');
fprintf( fid, '3600.0							500.0\n');
fprintf( fid, '36000.0							500.0\n');
fprintf( fid, '\n');
fprintf( fid, '\n');
fprintf( fid, '|| EvoDraggingParticles		||  *(ff)			|| time, fr\n');
fprintf( fid, '0.0								0.0E+0\n');
fprintf( fid, '3600.0							0.0E+0\n');
fprintf( fid, '\n');
fprintf( fid, '\n');

fprintf( fid, '|| RXGrainsPool				|| *(ffff)		|| Bunge e1, Bunge e2,  Bunge e3, tincub\n');
%fprintf( fid, '0.0			30.0		45.0        0.0\n');

% fprintf( fid, '0.0			30.0		45.0        0.0\n');
% fprintf( fid, '55.0		55.0		45.0        0.0\n');

% fprintf( fid, '55.0		55.0		45.0        0.0\n');
% fprintf( fid, '55.0		55.0		45.0        0.0\n');
% fprintf( fid, '55.0		55.0		45.0        0.0\n');
% fprintf( fid, '55.0		55.0		45.0        0.0\n');
% fprintf( fid, '55.0		55.0		45.0        0.0\n');
% fprintf( fid, '55.0		55.0		45.0        0.0\n');
% fprintf( fid, '55.0		55.0		45.0        0.0\n');
% fprintf( fid, '55.0		55.0		45.0        0.0\n');
% fprintf( fid, '55.0		55.0		45.0        0.0\n');
% fprintf( fid, '55.0		55.0		45.0        0.0\n');
% fprintf( fid, '55.0		55.0		45.0        0.0\n');
% fprintf( fid, '55.0		55.0		45.0        0.0\n');
% fprintf( fid, '55.0		55.0		45.0        0.0\n');
% fprintf( fid, '55.0		55.0		45.0        0.0\n');
% fprintf( fid, '55.0		55.0		45.0        0.0\n');
% fprintf( fid, '55.0		55.0		45.0        0.0\n');
% fprintf( fid, '55.0		55.0		45.0        0.0\n');
% fprintf( fid, '55.0		55.0		45.0        0.0\n');
fprintf( fid, '\n');
fprintf( fid, '\n');

%deformed grains
fprintf( fid, '|| DeformedGrainsPool		|| *(ffffff)        || Bunge e1, Bunge e2, Bunge e3, rho0, dsubav0, disoriav0\n');
%fprintf( fid, '%.3f\t\t%.3f\t\t%.3f\t\t%.3e\t\t%.3f\t\t%.3f\n', 60.0, 55.0, 45.0, 18.6*10^14, 0.0, 8.0 );


 newData1 = importdata( fdeformed, DELIMITER, HEADERLINES);
 vars = fieldnames(newData1);
 for i = 1:length(vars)
     assignin('base', vars{i}, newData1.(vars{i}));
 end
 nlines = length(data(:,1));
 for i=1:nlines
     fprintf( fid, '%.3f\t\t%.3f\t\t%.3f\t\t%.3e\t\t%.3f\t\t%.3f\n', data(i,1), data(i,3), data(i,5), data(i,7), data(i,9), data(i,11));
 end
 clearvars data



fprintf( fid, '\n');
fprintf( fid, '\n');
fprintf( fid, '|| HeuristicRXFrontListDefragmentation	||	*(sf)	|| Recrystallized fraction, local value\n');
fprintf( fid, '"X"		0.60\n');
fprintf( fid, '"X"		0.65\n');
fprintf( fid, '"X"		0.70\n');
fprintf( fid, '"X"		0.75\n');
fprintf( fid, '"X"		0.80\n');
fprintf( fid, '"X"		0.85\n');
fprintf( fid, '"X"		0.90\n');
fprintf( fid, '"X"		0.92\n');
fprintf( fid, '"X"		0.94\n');
fprintf( fid, '"X"		0.96\n');
fprintf( fid, '"X"		0.98\n');
fprintf( fid, '\n');
fprintf( fid, '\n');
fprintf( fid, '|| SingleGrainVolOverTime	||  *(sf)		|| Recrystallized fraction, local value\n');
fprintf( fid, '"X"		0.00\n');
fprintf( fid, '"X"		0.01\n');
fprintf( fid, '"X"		0.02\n');
fprintf( fid, '"X"		0.03\n');
fprintf( fid, '"X"		0.04\n');
fprintf( fid, '"X"		0.05\n');
fprintf( fid, '"X"		0.06\n');
fprintf( fid, '"X"		0.07\n');
fprintf( fid, '"X"		0.08\n');
fprintf( fid, '"X"		0.09\n');
fprintf( fid, '"X"		0.10\n');
fprintf( fid, '"X"		0.11\n');
fprintf( fid, '"X"		0.12\n');
fprintf( fid, '"X"		0.13\n');
fprintf( fid, '"X"		0.14\n');
fprintf( fid, '"X"		0.15\n');
fprintf( fid, '"X"		0.16\n');
fprintf( fid, '"X"		0.17\n');
fprintf( fid, '"X"		0.18\n');
fprintf( fid, '"X"		0.19\n');
fprintf( fid, '"X"		0.20\n');
fprintf( fid, '"X"		0.21\n');
fprintf( fid, '"X"		0.22\n');
fprintf( fid, '"X"		0.23\n');
fprintf( fid, '"X"		0.24\n');
fprintf( fid, '"X"		0.25\n');
fprintf( fid, '"X"		0.26\n');
fprintf( fid, '"X"		0.27\n');
fprintf( fid, '"X"		0.28\n');
fprintf( fid, '"X"		0.29\n');
fprintf( fid, '"X"		0.30\n');
fprintf( fid, '"X"		0.31\n');
fprintf( fid, '"X"		0.32\n');
fprintf( fid, '"X"		0.33\n');
fprintf( fid, '"X"		0.34\n');
fprintf( fid, '"X"		0.35\n');
fprintf( fid, '"X"		0.36\n');
fprintf( fid, '"X"		0.37\n');
fprintf( fid, '"X"		0.38\n');
fprintf( fid, '"X"		0.39\n');
fprintf( fid, '"X"		0.40\n');
fprintf( fid, '"X"		0.41\n');
fprintf( fid, '"X"		0.42\n');
fprintf( fid, '"X"		0.43\n');
fprintf( fid, '"X"		0.44\n');
fprintf( fid, '"X"		0.45\n');
fprintf( fid, '"X"		0.46\n');
fprintf( fid, '"X"		0.47\n');
fprintf( fid, '"X"		0.48\n');
fprintf( fid, '"X"		0.49\n');
fprintf( fid, '"X"		0.50\n');
fprintf( fid, '"X"		0.51\n');
fprintf( fid, '"X"		0.52\n');
fprintf( fid, '"X"		0.53\n');
fprintf( fid, '"X"		0.54\n');
fprintf( fid, '"X"		0.55\n');
fprintf( fid, '"X"		0.56\n');
fprintf( fid, '"X"		0.57\n');
fprintf( fid, '"X"		0.58\n');
fprintf( fid, '"X"		0.59\n');
fprintf( fid, '"X"		0.60\n');
fprintf( fid, '"X"		0.61\n');
fprintf( fid, '"X"		0.62\n');
fprintf( fid, '"X"		0.63\n');
fprintf( fid, '"X"		0.64\n');
fprintf( fid, '"X"		0.65\n');
fprintf( fid, '"X"		0.66\n');
fprintf( fid, '"X"		0.67\n');
fprintf( fid, '"X"		0.68\n');
fprintf( fid, '"X"		0.69\n');
fprintf( fid, '"X"		0.70\n');
fprintf( fid, '"X"		0.71\n');
fprintf( fid, '"X"		0.72\n');
fprintf( fid, '"X"		0.73\n');
fprintf( fid, '"X"		0.74\n');
fprintf( fid, '"X"		0.75\n');
fprintf( fid, '"X"		0.76\n');
fprintf( fid, '"X"		0.77\n');
fprintf( fid, '"X"		0.78\n');
fprintf( fid, '"X"		0.79\n');
fprintf( fid, '"X"		0.80\n');
fprintf( fid, '"X"		0.81\n');
fprintf( fid, '"X"		0.82\n');
fprintf( fid, '"X"		0.83\n');
fprintf( fid, '"X"		0.84\n');
fprintf( fid, '"X"		0.85\n');
fprintf( fid, '"X"		0.86\n');
fprintf( fid, '"X"		0.87\n');
fprintf( fid, '"X"		0.88\n');
fprintf( fid, '"X"		0.89\n');
fprintf( fid, '"X"		0.90\n');
fprintf( fid, '"X"		0.91\n');
fprintf( fid, '"X"		0.92\n');
fprintf( fid, '"X"		0.93\n');
fprintf( fid, '"X"		0.94\n');
fprintf( fid, '"X"		0.95\n');
fprintf( fid, '"X"		0.96\n');
fprintf( fid, '"X"		0.97\n');
fprintf( fid, '"X"		0.98\n');
fprintf( fid, '"X"		0.99\n');
fprintf( fid, '"X"		1.00\n');

fclose (fid);

'uds input file generated'