/* 
*	Copyright Markus Kühbach, 2014-2017
*	L. A. Barrales-Mora (quaternion library) and V. Mohles (I/O routines for reading the UDS file format)
*	contributed to the code.

*	SCORE is an MPI/OpenMP-parallel implementation of a cellular automaton model for the studying 
*	of microstructure evolution during the growth phase of static recrystallization in metallic alloys.
*	Its novelity is to solve an ensemble of independent simulation domains and to average their results
*	into an ensemble result. In comparison to the classical RVE-based approach this strategy enables 
*	studies with much higher statistical significance as orders of magnitude more grains can be studied
*	while these are solved at the same time in independent individual simulations which are thus executable
*	in parallel.
*	For this task, SCORE utilizes a two-layer data parallelism with a main layer of MPI-processes. 
*	Each of which solves for a queue of cellular automata domains. A second layer of OpenMP-thread 
*	parallelism accelerates the executing of each individual CA domain. The method is described in:

*	M. Kühbach, G. Gottstein, L. A. Barrales-Mora: A statistical ensemble cellular automaton 
*	microstructure model for primary recrystallization, Acta Materialia, Vol 107, 2016, p366
*	http://dx.doi.org/10.1016/j.actamat.2016.01.068

*	Further details, in particular to this implementation and the concept, are detailed in:
*	M. Kühbach: Efficient Recrystallization Microstructure Modeling by Utilizing Parallel Computation

*	The authors gratefully acknowledge the financial support from the Deutsche Forschungsgemeinschaft
*	(DFG) within the Reinhart Koselleck-Project (GO 335/44-1) and computing time grants kindly provided
*	by RWTH Aachen University and the FZ Jülich within the scope of the JARAHPC project JARA0076.


*	This file is part of SCORE.

*	SCORE is free software: you can redistribute it and/or modify
*	it under the terms of the GNU General Public License as published by
*	the Free Software Foundation, either version 3 of the License, or
*	(at your option) any later version.

*	SCORE is distributed in the hope that it will be useful,
*	but WITHOUT ANY WARRANTY; without even the implied warranty of
*	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*	GNU General Public License for more details.

*	You should have received a copy of the GNU General Public License
*	along with SCORE.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef __SCORE_KERNEL_H_INCLUDED__
//#define __SCORE_KERNEL_H_INCLUDED__

#include "SCORE_Kernel.h"
//#include <omp.h>


ensembleHdl::ensembleHdl()
{
	simid = DEFAULT_SIMID;

	ensRediscretization.ensRediscrTime = NULL;
	ensRediscretization.tensmin = 0.0;
	ensRediscretization.tensmax = 0.0;
	ensRediscretization.nslots = REDISCR_TIMESLOTS_DEFAULT;
	ensRediscretization.strategy = REDISCR_TIME_EQUIDISTANT;

	ensCAGeometry.nNucleiCSR = 1;
	ensCAGeometry.nboxedge_rd = CA_DIMENSIONS_MINIMUM;
	ensCAGeometry.nboxedge_td = CA_DIMENSIONS_MINIMUM;
	ensCAGeometry.nboxedge_nd = CA_DIMENSIONS_MINIMUM;
	ensCAGeometry.nboxarea_rdtd = SQR(CA_DIMENSIONS_MINIMUM);
	ensCAGeometry.nboxvol_rdtdnd = CUBE(CA_DIMENSIONS_MINIMUM);
	ensCAGeometry.cellsize = DEFAULT_CELLSIZE;
	ensCAGeometry.boxedge_rd = CA_DIMENSIONS_MINIMUM * DEFAULT_CELLSIZE;
	ensCAGeometry.boxedge_td = CA_DIMENSIONS_MINIMUM * DEFAULT_CELLSIZE;
	ensCAGeometry.boxedge_nd = CA_DIMENSIONS_MINIMUM * DEFAULT_CELLSIZE;
	ensCAGeometry.boxarea_rdtd = SQR(CA_DIMENSIONS_MINIMUM)*SQR(DEFAULT_CELLSIZE);
	ensCAGeometry.boxvol_rdtdnd = CUBE(CA_DIMENSIONS_MINIMUM)*CUBE(DEFAULT_CELLSIZE);

	ensPhysData.G = DEFAULT_SHEARMODULUS;
	ensPhysData.b = DEFAULT_BURGERSVECTOR;
	ensPhysData.G0 = DEFAULT_SHEARMODULUS;
	ensPhysData.bZeroCelsius = DEFAULT_BURGERSVECTOR;
	ensPhysData.dGdt = 0.0;
	ensPhysData.thermexp_C = DEFAULT_PUREALU_C;
	ensPhysData.thermexp_a = 0.0;
	ensPhysData.thermexp_b = 0.0;
	ensPhysData.Tmelt = DEFAULT_PUREALU_TMELT;

	ensPhysData.LAGBm0 = 1.0;
	ensPhysData.LAGBHact = DEFAULT_LAGB_HACT * echarge;
	ensPhysData.HAGBm0 = 1.0;
	ensPhysData.HAGBHact = DEFAULT_HAGB_HACT * echarge;
	ensPhysData.GSm0 = 1.0;
	ensPhysData.GSHact = DEFAULT_GS_HACT * echarge;

	ensPhysData.RH_HAGBm0 = 1.0;
	ensPhysData.RH_HAGBHact = DEFAULT_HAGB_HACT  * echarge;
	ensPhysData.RH_LAGBHAGBcut = 0.9;
	ensPhysData.RH_LAGBHAGBtrans = 5.0;
	ensPhysData.RH_LAGBHAGBexponent = 9.0;

	ensPhysData.defgmean_rd = DEFAULT_DEFGSIZE;
	ensPhysData.defgmean_td = DEFAULT_DEFGSIZE;
	ensPhysData.defgmean_nd = DEFAULT_DEFGSIZE;
	ensPhysData.defgmean_poisson = DEFAULT_DEFGSIZE;

	ensDragData.ZenerAlpha = 0.0;
	ensDragData.ZenerGamma = 0.0;
	ensDragData.fr = 0.0;
	ensDragData.ZenerConsider = DISPERSOIDDRAG_NO;

	ensRecoveryModel.RecoveryConsider = RECOVERY_NO;
	ensRecoveryModel.VacancyDiffGeo = 1.0;
	ensRecoveryModel.VacancyDiffD0 = 1.0;
	ensRecoveryModel.VacancyDiffHact = 1.0;
	ensRecoveryModel.SoluteDiffD0 = 1.0;
	ensRecoveryModel.SoluteLsPropFactor = 1.0;
	ensRecoveryModel.SoluteConcentration = 1.0;
	ensRecoveryModel.SoluteDiffHact = 1.0;
	ensRecoveryModel.NesAlpha3 = 1.0;
	ensRecoveryModel.NesKappa2 = 1.0;
	ensRecoveryModel.NesC3 = 1.0;
	ensRecoveryModel.NesC4 = 1.0;
	ensRecoveryModel.NesC5 = 1.0;

	ensNucleationModel.gbnucleation = GBNUCLEATION_NO;
	ensNucleationModel.csrnucleation = CSRNUCLEATION_NO;
	ensNucleationModel.clustnucleation = CLUSTERNUCLEATION_NO;
	ensNucleationModel.defaultnucdensity = DEFAULT_NUCDENSITY_CSR;

	maxfillperstep = DEFAULT_MAXFILLPERSTEP;
	ensMemGuard = 0.0;
	initialRelCellCaching = DEFAULT_RELCELLCACHING;
	transientRelCellRecaching = DEFAULT_TRANSRELCELL;
	ebsdstepsize = 0.0;
	XMAX = DEFAULT_XMAX;
	TMAX = DEFAULT_TMAX;
	NTSTEPSMAX = DEFAULT_NMAX;

	myRank = MASTER;
	nRanks = 1;
	nworldCAs = 0;
	nensCAs = 0;

	percolation = PERCOLATION_ANALYZE_NO;
	defmsmethod = CUBOID_DEFMS;
	outopt_rendermethod = RENDERING_MSNO;
	outopt_rendercolormodel = RENDERING_COLOR_GRAINID;
	outopt_localrenderboundaries = RENDERING_BOUNDARIES_NO;
	outopt_logboundaries = OUTPUT_LOGBND_NO;
	outopt_rxfront = OUTPUT_RXFRONTSTATS_NO;
	outopt_threadprofiling = OUTPUT_THREADPROFILING_NO;
	outopt_singlegrainevo = OUTPUT_SINGLEGRAIN_NO;
	outopt_generateDAMASKgeom = OUTPUT_DAMASK_GEOMETRYFILE_NO;
	outopt_artsemebsd = OUTPUT_SEMEBSD_NO;

	onthefly_defragmentation = false;
	ensembleSuccess = true;
	postProcessing = true;
	experimentInput = false;

	ensembleprngseed = DEFAULT_SEED;

	score_input = NULL;
	ensembleprng.init( DEFAULT_PRNG_SEED );


	prof_t0 = 0.0;
	prof_tstart = 0.0;
	prof_tend = 0.0;
}


ensembleHdl::~ensembleHdl()
{
	//vector of caHdl call constructor of the caHdl class which has to assure to free memory in loginfo pointer snippets
	//randomClass has own destructor
	delete [] ensRediscretization.ensRediscrTime; ensRediscretization.ensRediscrTime = NULL;
	ensMemGuard = ensMemGuard - (ensRediscretization.nslots * sizeof(double));
}


void ensembleHdl::init_ensprng( bool SetDissimilarSeedsForAll )
{
	if ( myRank == MASTER ) 
		cout << "Initializing PRNG..." << endl;

	//Initialize random number generators for the individual MPI processes
	//true - seed is (-1*myRank)-1, 
	if ( SetDissimilarSeedsForAll == true) {
		long newseed = (long) (-1 * myRank);
		newseed = newseed - 1; //avoid seed 0
		ensembleprng.init ( newseed );

#ifdef REPORTSTYLE_DEVELOPER
		cout << myRank << " my new ensembleprng seed is:" << newseed << endl;
#endif

		//set also the seed of my local math library, works because there is ONLY ONE INSTANCE OF AN ENSEMBLEHDL PUBLIC SHARING mathMethods class object PER PROCESS!
		//this->setprng( newseed );
		ensmath.setprng( newseed );
	}
	else {
		//false - all nodes get the same seed, therefore generate exactly the same chain of random numbers!
		ensembleprng.init( DEFAULT_SEED );
		//this->setprng( DEFAULT_SEED );
		ensmath.setprng( DEFAULT_SEED );
	}
}


void ensembleHdl::init_ensprng( int targetseed )
{
	if ( myRank == MASTER )
		cout << "Initializing PRNG with Diehl protocol i.e. pass explicit seed" << endl;

	//Initialize random number generators for the individual MPI processes
	//all nodes get the same seed, therefore generate exactly the same chain of random numbers!
	ensembleprngseed = -targetseed;
	ensembleprng.init( -targetseed );
	ensmath.setprng( -targetseed );

#ifdef REPORTSTYLE_DEVELOPER
		cout << myRank << " my new ensembleprng seed is:" << -targetseed << endl;
#endif
}



void ensembleHdl::init_mpidatatypes( void )
{
	if ( myRank == MASTER ) 
		cout << "Initializing MPI I/O Datatypes..." << endl;

	//MPI_IO_CAProfilingData
	int elementCounts[2] = {2,7};
	MPI_Aint displacements[2] = {0, 2*8};
	MPI_Datatype oldTypes[2] = {MPI_LONG, MPI_DOUBLE};
	MPI_Type_create_struct(2, elementCounts, displacements, oldTypes, &MPI_IO_CAProfilingInfoData_Type);

	//MPI_IO_CAPhysicsInfoData
	int elementCounts2[2] = {12, 1};
	MPI_Aint displacements2[2] = { 0, 12*4}; //was 12*8 for MPI_LONG
	MPI_Datatype oldTypes2[2] = {MPI_UNSIGNED, MPI_DOUBLE};
	MPI_Type_create_struct(2, elementCounts2, displacements2, oldTypes2, &MPI_IO_CAPhysicsInfoData_Type);

	//MPI_IO_FinalGSDInfoData
	int elementCounts3[2] = {2, 1};
	MPI_Aint displacements3[2] = { 0, 2*8};
	MPI_Datatype oldTypes3[2] = {MPI_DOUBLE, MPI_UNSIGNED };
	MPI_Type_create_struct(2, elementCounts3, displacements3, oldTypes3, &MPI_IO_FinalGSDInfoData_Type);

	//MPI_IO_FGSDComplete
	int elementCounts4[2] = {2, 2};
	MPI_Aint displacements4[2] = { 0, 2*8};
	MPI_Datatype oldTypes4[2] = {MPI_DOUBLE, MPI_INT };
	MPI_Type_create_struct(2, elementCounts4, displacements4, oldTypes4, &MPI_IO_FGSDComplete_Type);
	
	//commit types
	MPI_Type_commit(&MPI_IO_CAProfilingInfoData_Type);
	MPI_Type_commit(&MPI_IO_CAPhysicsInfoData_Type);
	MPI_Type_commit(&MPI_IO_FinalGSDInfoData_Type);
	MPI_Type_commit(&MPI_IO_FGSDComplete_Type);
}

bool ensembleHdl::init_parameter_runtime( void ) {
	bool status = true;

	dataBlockP runtimecontrol = readDataBlock("RuntimeControl", score_input);

	XMAX = geTReal("XMAX", runtimecontrol );
	if ( (XMAX < 0.0 || XMAX > 1.0) && myRank == MASTER ) { status = false; cout << "ERROR::Invalid argument in XMAX." << endl; }

	TMAX = geTReal("TIMEMAX", runtimecontrol );
	if ( (TMAX < 0.0 || TMAX > DEFAULT_TMAX) && myRank == MASTER ) { status = false; cout << "ERROR::Invalid argument in TMAX." << endl; }

	NTSTEPSMAX = geTInt("NMAX", runtimecontrol );
	if ( (NTSTEPSMAX < 1 || NTSTEPSMAX > DEFAULT_NMAX || NTSTEPSMAX >= INTEGER_RANGE_MAX) && myRank == MASTER ) { status = false; cout << "ERROR::Invalid argument for NMAX." << endl; }

	percolation = PERCOLATION_ANALYZE_NO;
	if ( geTInt("PERCOLATIONANALYSIS", runtimecontrol) == 1 ) percolation = PERCOLATION_ANALYZE_YES_NOSZDISTR;
	if ( geTInt("PERCOLATIONANALYSIS", runtimecontrol) == 2 ) percolation = PERCOLATION_ANALYZE_YES_SZDISTR;

	maxfillperstep = geTReal("MAXFILLPERSTEP", runtimecontrol );
	if ( (maxfillperstep < DEFAULT_MINFILLIN || maxfillperstep > DEFAULT_MAXFILLIN) && myRank == MASTER ) { maxfillperstep = DEFAULT_MAXFILLPERSTEP; cout << "WARNING::MAXFILLPERSTEP was chosen wrong, should be 0.01 < 0.2, was reset to = " << maxfillperstep << endl; }

	int rediscr_nslots = geTInt("RediscretizationSteps", runtimecontrol );
	if ( rediscr_nslots < REDISCR_TIMESLOTS_MIN || rediscr_nslots > REDISCR_TIMESLOTS_MAX || rediscr_nslots >= INTEGER_RANGE_MAX ) { cout << "WARNING::Rediscretization time slots reset to " << REDISCR_TIMESLOTS_DEFAULT << endl; }
	ensRediscretization.nslots = rediscr_nslots;

	initialRelCellCaching = geTReal("InitialCacheSizeRXCells", runtimecontrol );
	if ( (initialRelCellCaching < DEFAULT_MINRELCELL || initialRelCellCaching > DEFAULT_MAXRELCELL) && myRank == MASTER ) { initialRelCellCaching = DEFAULT_RELCELLCACHING; cout << "WARNING::InitialCacheSizeRXCells was chosen unwise, it is reset to = " << initialRelCellCaching << endl; }
	transientRelCellRecaching = geTReal("ReCacheSizeRXCells", runtimecontrol );
	if ( (transientRelCellRecaching < DEFAULT_MINTRANSRELCELL || transientRelCellRecaching > DEFAULT_MAXTRANSRELCELL) && myRank == MASTER ) { transientRelCellRecaching = DEFAULT_TRANSRELCELL; cout << "WARNING::ReCacheSizeRXCells was chosen unwise, it is reset to = " << transientRelCellRecaching << endl; }

	onthefly_defragmentation = false;
	uint32_t onthefly_defrag = 0;
	onthefly_defrag = geTInt( "OnTheFlyDefragmentation", runtimecontrol );
	if ( onthefly_defrag == 1 ) onthefly_defragmentation = true;

	//possibility to switch off postProcessing in order to utilize only the CA solver
	postProcessing = false;
	uint32_t kuehbach_su_analysis = 0;
	kuehbach_su_analysis = geTInt("SolitaryUnitPostProcessing", runtimecontrol );
	if ( kuehbach_su_analysis == 1 ) postProcessing = true;

	experimentInput = false;
	uint32_t support_ebsddata = 0;
	support_ebsddata = geTInt("ReadEBSDMapping", runtimecontrol );
	if ( support_ebsddata == 1 ) experimentInput = true;

	this->ebsdstepsize = geTReal("ReadEBSDStepsize", runtimecontrol );
	if ( ebsdstepsize <= DEFAULT_EBSDSTEPSIZE_MINIMUM ) { status = false; cout << "ERROR::Invalid argument in ReadEBSDStepsize." << endl; }


	delete runtimecontrol;

	return status;
}

bool ensembleHdl::init_parameter_output( void ) {
	bool status = true;

	dataBlockP outputopts = readDataBlock( "SimulationOutput", score_input );

	outopt_rendermethod = RENDERING_MSNO;
	long val = geTInt("RenderingOfMicrostructure", outputopts);
	if ( val == 2 ) outopt_rendermethod = RENDERING_MS2D;
	if ( val == 3 ) outopt_rendermethod = RENDERING_MS3D;

	if ( outopt_rendermethod != RENDERING_MSNO ) {
		renderwindow.xmin = geTReal( "RenderingSectionXMIN", outputopts );
		renderwindow.xmax = geTReal( "RenderingSectionXMAX", outputopts );
		renderwindow.ymin = geTReal( "RenderingSectionYMIN", outputopts );
		renderwindow.ymax = geTReal( "RenderingSectionYMAX", outputopts );
		renderwindow.zmin = geTReal( "RenderingSectionZMIN", outputopts );
		renderwindow.zmax = geTReal( "RenderingSectionZMAX", outputopts );
		bool modified = false; 
		if ( renderwindow.xmin < 0.0 ) { renderwindow.xmin = 0.0; modified = true; }
		if ( renderwindow.ymin < 0.0 ) { renderwindow.ymin = 0.0; modified = true; }
		if ( renderwindow.zmin < 0.0 ) { renderwindow.zmin = 0.0; modified = true; }
		if ( renderwindow.xmax > 1.0 ) { renderwindow.xmax = 1.0; modified = true; }
		if ( renderwindow.ymax > 1.0 ) { renderwindow.ymax = 1.0; modified = true; }
		if ( renderwindow.zmax > 1.0 ) { renderwindow.zmax = 1.0; modified = true; }
		if ( renderwindow.xmin >= renderwindow.xmax ) { renderwindow.xmin = 0.0; renderwindow.xmax = 0.0; modified = true; }
		if ( renderwindow.ymin >= renderwindow.ymax ) { renderwindow.ymin = 0.0; renderwindow.ymax = 0.0; modified = true; }
		if ( renderwindow.zmin >= renderwindow.zmax ) { renderwindow.zmin = 0.0; renderwindow.zmax = 0.0; modified = true; }
		if ( modified == true ) {
			cout << "WARNING::I modified the rendering window to [0.0,1.0]^3 as the input was faulty!" << endl;
		}
	}

	outopt_rendercolormodel = RENDERING_COLOR_GRAINID;
	val = geTInt( "RenderingColorModel", outputopts);
	if ( val == 2 ) outopt_rendercolormodel = RENDERING_COLOR_IPFZ;
	if ( val == 3 ) outopt_rendercolormodel = RENDERING_COLOR_RHO;

	outopt_renderfileformat = RENDERING_FILEFORMAT_RAW;
	if ( geTInt( "RenderingFileFormat", outputopts) == 2 )
		outopt_renderfileformat = RENDERING_FILEFORMAT_HDF5;

	outopt_localrenderboundaries = RENDERING_BOUNDARIES_NO;
	if ( geTInt("RenderingBoundaries", outputopts) == 1 ) 
		outopt_localrenderboundaries = RENDERING_BOUNDARIES_YES;

	outopt_logboundaries = OUTPUT_LOGBND_NO;
	if ( geTInt( "OutputLogBoundaries", outputopts) == 1 )
		outopt_logboundaries = OUTPUT_LOGBND_YES;

	outopt_rxfront = OUTPUT_RXFRONTSTATS_NO;
	if ( geTInt( "OutputRXFrontStats", outputopts) == 1 )
		outopt_rxfront = OUTPUT_RXFRONTSTATS_YES;

	outopt_threadprofiling = OUTPUT_THREADPROFILING_NO;
	if ( geTInt( "OutputThreadProfiling", outputopts ) == 1 )
		outopt_threadprofiling = OUTPUT_THREADPROFILING_YES;

	outopt_singlegrainevo = OUTPUT_SINGLEGRAIN_NO;
	val = geTInt( "OutputSingleGrainStats", outputopts);
	if ( val == 1 ) outopt_singlegrainevo = OUTPUT_SINGLEGRAIN_ASCII;
	if ( val == 2 ) outopt_singlegrainevo = OUTPUT_SINGLEGRAIN_BINARY;

	outopt_generateDAMASKgeom = OUTPUT_DAMASK_GEOMETRYFILE_NO;
	val = geTInt( "OutputDAMASKGeometryFile", outputopts );
	if ( val == 1 ) outopt_generateDAMASKgeom = OUTPUT_DAMASK_GEOMETRYFILE_YES;

	outopt_artsemebsd = OUTPUT_SEMEBSD_NO;
	val = geTInt( "OutputArtificialSEMEBSD", outputopts );
	if ( val == 1 ) outopt_artsemebsd = OUTPUT_SEMEBSD_YES;


	delete outputopts;

	return status;
}

bool ensembleHdl::init_parameter_ensemble( void ) {
	bool status = true;

	dataBlockP ensembledef = readDataBlock("EnsembleDefinition", score_input);

	nworldCAs = geTInt("CAEnsembleSize", ensembledef);
	if ( (nworldCAs < 1 || nworldCAs >= INTEGER_RANGE_MAX || nworldCAs > DEFAULT_MAX_NCAS || nworldCAs < nRanks ) && myRank == MASTER ) { status = false; cout << "ERROR::Invalid argument for the number of domains." << endl; }

	ensCAGeometry.cellsize = geTReal("CellSize", ensembledef) * MICRON2METER;
	if ( (ensCAGeometry.cellsize < DEFAULT_CELLSIZE_MIN || ensCAGeometry.cellsize > DEFAULT_CELLSIZE_MAX) && myRank == MASTER ) { status = false; cout << "ERROR::Invalid argument for the cellsize." << endl; }

	ensCAGeometry.nboxedge_rd = geTInt("3DCAEdgeLengthInCellsRD", ensembledef);
	//if ( (ensCAGeometry.nboxedge_rd < CA_DIMENSIONS_MINIMUM || ensCAGeometry.nboxedge_rd > CA_DIMENSIONS_MAXIMUM) && myRank == MASTER ) { status = false; cout << "ERROR::Invalid argument for automaton size RD." << endl; }
	if ( ensCAGeometry.nboxedge_rd < CA_DIMENSIONS_MINIMUM && myRank == MASTER ) { status = false; cout << "ERROR::Invalid argument for automaton size RD." << endl; }

	ensCAGeometry.nboxedge_td = geTInt("3DCAEdgeLengthInCellsTD", ensembledef);
	//if ( (ensCAGeometry.nboxedge_td < CA_DIMENSIONS_MINIMUM || ensCAGeometry.nboxedge_td > CA_DIMENSIONS_MAXIMUM) && myRank == MASTER ) { status = false; cout << "ERROR::Invalid argument for automaton size TD." << endl; }
	if ( ensCAGeometry.nboxedge_td < CA_DIMENSIONS_MINIMUM && myRank == MASTER ) { status = false; cout << "ERROR::Invalid argument for automaton size TD." << endl; }

	ensCAGeometry.nboxedge_nd = geTInt("3DCAEdgeLengthInCellsND", ensembledef);
	//if ( (ensCAGeometry.nboxedge_nd < CA_DIMENSIONS_MINIMUM || ensCAGeometry.nboxedge_nd > CA_DIMENSIONS_MAXIMUM) && myRank == MASTER ) { status = false; cout << "ERROR::Invalid argument for automaton size ND." << endl; }
	if ( ensCAGeometry.nboxedge_nd < CA_DIMENSIONS_MINIMUM && myRank == MASTER ) { status = false; cout << "ERROR::Invalid argument for automaton size ND." << endl; }

	ensCAGeometry.nboxarea_rdtd = ensCAGeometry.nboxedge_rd * ensCAGeometry.nboxedge_td;
	ensCAGeometry.nboxvol_rdtdnd = ensCAGeometry.nboxedge_rd * ensCAGeometry.nboxedge_td * ensCAGeometry.nboxedge_nd;
	if ( ensCAGeometry.nboxvol_rdtdnd > CUBE(CA_DIMENSIONS_MAXIMUM) && myRank == MASTER ) { status = false; cout << "ERROR::In total the resulting CA would contain more cells than what the implementation can handle currently!" << endl; }

	ensCAGeometry.boxedge_rd = ensCAGeometry.cellsize * ensCAGeometry.nboxedge_rd;
	ensCAGeometry.boxedge_td = ensCAGeometry.cellsize * ensCAGeometry.nboxedge_td;
	ensCAGeometry.boxedge_nd = ensCAGeometry.cellsize * ensCAGeometry.nboxedge_nd;
	ensCAGeometry.boxarea_rdtd = ensCAGeometry.nboxarea_rdtd * SQR(ensCAGeometry.cellsize);
	ensCAGeometry.boxvol_rdtdnd = ensCAGeometry.nboxvol_rdtdnd * CUBE(ensCAGeometry.cellsize);

	//modify render window
	renderwindow.xmi = (renderwindow.xmin * ((double) ensCAGeometry.nboxedge_rd));
	renderwindow.xmx = (renderwindow.xmax * ((double) ensCAGeometry.nboxedge_rd));
	renderwindow.ymi = (renderwindow.ymin * ((double) ensCAGeometry.nboxedge_td));
	renderwindow.ymx = (renderwindow.ymax * ((double) ensCAGeometry.nboxedge_td));
	renderwindow.zmi = (renderwindow.zmin * ((double) ensCAGeometry.nboxedge_nd));
	renderwindow.zmx = (renderwindow.zmax * ((double) ensCAGeometry.nboxedge_nd));
	//set limits to inclusive coordinates
	if ( renderwindow.xmx == ensCAGeometry.nboxedge_rd ) renderwindow.xmx--; //does not substract from 0 because SU has dimensions of at least CA_DIMENSIONS_MINIMUM > 0
	if ( renderwindow.ymx == ensCAGeometry.nboxedge_td ) renderwindow.ymx--;
	if ( renderwindow.zmx == ensCAGeometry.nboxedge_nd ) renderwindow.zmx--;
	QUICKASSERT( renderwindow.xmi >= 0 );
	QUICKASSERT( renderwindow.xmx < ensCAGeometry.nboxedge_rd );
	QUICKASSERT( renderwindow.ymi >= 0 );
	QUICKASSERT( renderwindow.ymx < ensCAGeometry.nboxedge_td );
	QUICKASSERT( renderwindow.zmi >= 0 );
	QUICKASSERT( renderwindow.zmx < ensCAGeometry.nboxedge_nd );

	ensPhysData.defgmean_rd = geTReal("MeanDefGrainSizeinRD", ensembledef) * MICRON2METER; // RD - x
	ensPhysData.defgmean_td = geTReal("MeanDefGrainSizeinTD", ensembledef) * MICRON2METER; // TD - y
	ensPhysData.defgmean_nd = geTReal("MeanDefGrainSizeinND", ensembledef) * MICRON2METER; // ND - z
	double mindiscretization = DEFAULT_MIN_GRAINDISCR * ensCAGeometry.cellsize;
	if ( (ensPhysData.defgmean_rd < mindiscretization || ensPhysData.defgmean_td < mindiscretization || ensPhysData.defgmean_nd < mindiscretization) && myRank == MASTER) { status = false; cout << "ERROR::Invalid argument for the average size of the cuboid deformed grain." << endl; } //no upper bound check necessary

	defmsmethod = CUBOID_DEFMS;
	long val = geTInt("DefStructureSynthesis", ensembledef);
	if ( val == 1 ) defmsmethod = CUBOID_DEFMS;
	if ( val == 2 ) defmsmethod = POISSONVORONOI_DEFMS;
	if ( val == 3 ) defmsmethod = CPFEMDAMASK_DEFMS;
	if ( val == 4 ) defmsmethod = SEMEBSD_2D3D_COLUMNSINGLE;
	if ( val == 5 ) defmsmethod = SEMEBSD_2D3D_COLUMNSTACK;
	if ( val == 6 ) defmsmethod = SEMEBSD_2D3D_COLUMNSHIFT;

	ensPhysData.defgmean_poisson = geTReal("MeanDefSizePoissonDiameter", ensembledef) * MICRON2METER;
	if ( (ensPhysData.defgmean_poisson < mindiscretization && myRank == MASTER ) ) { status = false; cout << "ERROR::Invalid argument for the average size of deformed grain for method poisson." << endl; }

	delete ensembledef;

	return status;
}

bool ensembleHdl::init_parameter_nucleation( void ) {
	bool status = true;

	dataBlockP nucleationmodeling = readDataBlock( "NucleationModeling", score_input );

	ensNucleationModel.defaultnucdensity = geTInt("NucleationDensityLocalCSR", nucleationmodeling );
	if ( ensNucleationModel.defaultnucdensity >= INTEGER_RANGE_MAX && myRank == MASTER ) { status = false; cout << "ERROR::Too high NucleationDensityLocalCSR!" << endl; }
	if ( ensNucleationModel.defaultnucdensity < DEFAULT_NUCDENSITY_CSR && myRank == MASTER ) { status = false; cout << "ERROR::Too low NucleationDensityLocalCSR!" << endl; }
	long nnucpractical = (ensCAGeometry.nboxvol_rdtdnd / ensNucleationModel.defaultnucdensity); //cells per nuc, assuming CSR, v=const
	if ( ( nnucpractical < MIN_DISCRETIZATION_GRAIN) && myRank == MASTER ) { cout << "WARNING::Inacceptable low discretization for this high number of nuclei." << endl; }

	//CSR nucleation
	ensNucleationModel.csrnucleation = CSRNUCLEATION_NO;
	long val = geTInt("CSRNucleation", nucleationmodeling);
	if ( val == 1 ) ensNucleationModel.csrnucleation = CSRNUCLEATION_YES_ENFORCED;
	if ( val == 2 ) ensNucleationModel.csrnucleation = CSRNUCLEATION_YES_PICKRND;
	if ( val == 3 ) ensNucleationModel.csrnucleation = CSRNUCLEATION_YES_DIEHL;
	if ( val == 4 ) ensNucleationModel.csrnucleation = CSRNUCLEATION_DIEHL_RANDOMSO3;
	if ( val == 5 ) ensNucleationModel.csrnucleation = CSRNUCLEATION_DIEHL_RANDOMDEFORI;
	if ( val == 6 ) ensNucleationModel.csrnucleation = CSRNUCLEATION_DIEHL_SCATTEREXISTENT;

	//Clustering nucleation
	ensNucleationModel.clustnucleation = CLUSTERNUCLEATION_NO;
	if ( geTInt("ClusteredNucleation", nucleationmodeling) == 1 )
		ensNucleationModel.clustnucleation = CLUSTERNUCLEATION_YES_PICKRND;

	ensNucleationModel.cluster_nclust = geTInt( "ClusteredNucleationNCluster", nucleationmodeling);
	ensNucleationModel.cluster_lambda = geTReal( "ClusteredNucleationLCluster", nucleationmodeling);
	ensNucleationModel.cluster_rvesize = geTReal( "ClusteredNucleationScaling", nucleationmodeling) * MICRON2METER;
	ensNucleationModel.cluster_a = geTInt( "ClusteredNucleationExtendA", nucleationmodeling);
	ensNucleationModel.cluster_b = geTInt( "ClusteredNucleationExtendB", nucleationmodeling);
	ensNucleationModel.cluster_c = geTInt( "ClusteredNucleationExtendC", nucleationmodeling);

	//Nucleation at pre-existent grain boundaries
	ensNucleationModel.gbnucleation = GBNUCLEATION_NO;
	//at the moment no GBnucleation in OpenMP parallel version
	if ( geTInt("GBNucleation", nucleationmodeling) == 1 )
		ensNucleationModel.gbnucleation = GBNUCLEATION_YES_PHYS_LAGBHAGB;
	if ( geTInt("GBNucleation", nucleationmodeling) == 2 )
		ensNucleationModel.gbnucleation = GBNUCLEATION_YES_PHYS_ONLYHAGB;
	if ( geTInt("GBNucleation", nucleationmodeling) == 3 ) 
		ensNucleationModel.gbnucleation = GBNUCLEATION_YES_PICKRND;

	ensNucleationModel.gbnucleation_dens2num = geTReal("GBNucDensity2Number", nucleationmodeling);

	ensNucleationModel.gbnucleation_drho2dens = geTReal("GBNucRhoDiff2Density", nucleationmodeling);
	ensNucleationModel.gbnucleation_scatter = DEG2RAD(geTReal("GBNucMaximumScatter", nucleationmodeling));
	ensNucleationModel.tincub_rayleigh_sigma = geTReal("IncubationTimeScatter", nucleationmodeling);
	if ( ensNucleationModel.gbnucleation_dens2num < 0.0 || ensNucleationModel.gbnucleation_drho2dens < 0.0 || ensNucleationModel.gbnucleation_scatter < 0.0 ) { status = false; cout << "ERROR::Invalid grain boundary nucleation properties!" << endl; }
	if ( ensNucleationModel.csrnucleation == CSRNUCLEATION_NO && ensNucleationModel.gbnucleation == GBNUCLEATION_NO && ensNucleationModel.clustnucleation == CLUSTERNUCLEATION_NO ) { status = false; cout << "ERROR::Choose a valid nucleation model!" << endl; }

	//Incubation time model
	ensNucleationModel.tincubmodel = TINCUB_SITESATURATION;
	if ( geTInt("IncubationTimeModel", nucleationmodeling) != 1 ) {
		status = false; cout << "ERROR::Only site-saturated nucleation implemented in the MPI/OpenMP model!" << endl;
	}
	/*if ( geTInt("IncubationTimeModel", nucleationmodeling) == 2 ) {
		ensNucleationModel.tincubmodel = TINCUB_TIMEDEPENDENT;
		if ( ensNucleationModel.tincub_rayleigh_sigma < MINIMUM_RAYLEIGH_SIGMA ) {
			status = false; cout << "ERROR::Invalid low incubation time scatter!" << endl;
		}
	}*/

	delete nucleationmodeling;

	return status;
}

bool ensembleHdl::init_parameter_material( void ) {
	bool status = true;

	dataBlockP matpara = readDataBlock("MaterialProperties", score_input);

	ensPhysData.G0 = geTReal("ZeroKelvinShearModulusG0", matpara);
	ensPhysData.dGdt = geTReal("FirstOrderdG0dT", matpara);
	ensPhysData.bZeroCelsius = geTReal("ZeroCelsiusBurgersVector", matpara);
	ensPhysData.thermexp_C = geTReal("AlloyConstantThermalExpCoeff", matpara); //thermal lattice expansion model
	ensPhysData.thermexp_a = geTReal("FirstOrderThermalExpCoeff", matpara);
	ensPhysData.thermexp_b = geTReal("SecondOrderThermalExpCoeff", matpara);
	ensPhysData.Tmelt = geTReal("MeltingTemperature", matpara) + TOFFSET;
	if ( ensPhysData.G0 <= 0.0 || ensPhysData.bZeroCelsius <= 0.0 || ensPhysData.Tmelt <= 0.0 ) { status = false; cout << "ERROR::Inrealistic material parameter!" << endl; }

	delete matpara;

	return status;
}

bool ensembleHdl::init_parameter_mobilities( void ) {
	bool status = true;

	dataBlockP mobilitymodeling = readDataBlock( "BoundaryMobilityModel", score_input );

	mobilitymodel = MOBILITYMODEL_SEBALDGOTTSTEIN;
	if ( geTInt("MobilityModel", mobilitymodeling) == 2 )
		mobilitymodel = MOBILITYMODEL_ROLLETTHOLM;

	//parameter of the Sebald-Gottstein inspired mobility model in which the mobility is solely dependent on the disorientation of the two crystals
	//adjointing the boundary, when the disorientation < LAGB_TO_HAGB_TRANSITION the boundary is essentially immobile, otherwise mobile as a HAGB
	//it is considered particularly mobile when the disorientation is close to a 40deg111 in misorientation space
	ensPhysData.LAGBm0 = geTReal("LAGBm0", mobilitymodeling);
	ensPhysData.LAGBHact = geTReal("LAGBHact", mobilitymodeling) * echarge; //all other codes work with 1.602e-19
	ensPhysData.HAGBm0 = geTReal("HAGBm0", mobilitymodeling);
	ensPhysData.HAGBHact = geTReal("HAGBHact", mobilitymodeling) * echarge;
	ensPhysData.GSm0 = geTReal("GSm0", mobilitymodeling);
	ensPhysData.GSHact = geTReal("GSHact", mobilitymodeling) * echarge;

	//on the other hand Rollett and Holm inspired by Humphreys for instance often consider the LAGB_TO_HAGB transition not as sharp and despite the possibility of 
	//close to 40deg111 boundaries modelled additionally more mobile than random high-angle grain boundaries (HAGB). They consider a sigmoidal function of the form
	//m = 1.0 - cut*exp(-trans * (disori/LAGB_TO_HAGB_TRANSITION)^exponent) also inspired by exp. results on tilt boundaries in Al performed by Winning and Molodov
	//MK::however often the parameter RHModelLAGBHAGBExponent is set to 4, this is physically unrealistic because then LAGB are sufficiently 
	//mobile against aforementioned experimental observations
	ensPhysData.RH_HAGBm0 = geTReal("RHModelHAGBm0", mobilitymodeling );
	ensPhysData.RH_HAGBHact = geTReal("RHModelHAGBHact", mobilitymodeling ) * echarge;
	ensPhysData.RH_LAGBHAGBcut = geTReal("RHModelLAGBHAGBCut", mobilitymodeling );
	ensPhysData.RH_LAGBHAGBtrans = geTReal("RHModelLAGBHAGBTrans", mobilitymodeling );
	ensPhysData.RH_LAGBHAGBexponent = geTReal("RHModelLAGBHAGBExponent", mobilitymodeling );

	if ( mobilitymodel == MOBILITYMODEL_SEBALDGOTTSTEIN )
		if ( ensPhysData.LAGBm0 <= 0.0 || ensPhysData.LAGBHact <= 0.0 || ensPhysData.HAGBm0 <= 0.0 || ensPhysData.HAGBHact <= 0.0 || ensPhysData.GSm0 <= 0.0 || ensPhysData.GSHact <= 0 ) { status = false; cout << "ERROR::Unrealistic mobility model parameterization!" << endl; }
	if ( mobilitymodel == MOBILITYMODEL_ROLLETTHOLM )
		if ( ensPhysData.RH_HAGBm0 <= 0.0 || ensPhysData.RH_HAGBHact <= 0.0 || ensPhysData.RH_LAGBHAGBcut <= 0.0 || ensPhysData.RH_LAGBHAGBtrans <= 0.0 || ensPhysData.RH_LAGBHAGBexponent <= 0.0 ) { status = false; cout << "ERROR::Unrealistic mobility model parameterization!" << endl; }

	delete mobilitymodeling;

	return status;
}

bool ensembleHdl::init_parameter_recovery( void ) {
	bool status = true;

	dataBlockP recoverypara = readDataBlock("RecoveryParameter", score_input);

	ensRecoveryModel.RecoveryConsider = RECOVERY_NO;
	int recoverymodel = geTInt( "RecoveryConsider", recoverypara );
	if ( recoverymodel == 1 ) ensRecoveryModel.RecoveryConsider = RECOVERY_NES_VACCOREDIFF;
	if ( recoverymodel == 2 ) ensRecoveryModel.RecoveryConsider = RECOVERY_NES_SOLUTEDRAG;
	if ( recoverymodel == 3 ) ensRecoveryModel.RecoveryConsider = RECOVERY_NES_MICHALAKPAXTON;
	ensRecoveryModel.VacancyDiffGeo = geTReal( "RecoveryVacancyDiffGeometry", recoverypara );
	ensRecoveryModel.VacancyDiffD0 = geTReal( "RecoveryVacancyDiffPreexp", recoverypara );
	ensRecoveryModel.VacancyDiffHact = geTReal( "RecoveryVacancyDiffHact", recoverypara ) * echarge;
	ensRecoveryModel.SoluteDiffD0 = geTReal( "RecoverySoluteDiffPreexp", recoverypara );
	ensRecoveryModel.SoluteDiffHact = geTReal( "RecoverySoluteDiffHact", recoverypara ) * echarge;
	ensRecoveryModel.SoluteLsPropFactor = geTReal( "RecoverySoluteLsPropFactor", recoverypara );
	ensRecoveryModel.SoluteConcentration = geTReal ( "RecoverySoluteConcentration", recoverypara );
	ensRecoveryModel.NesAlpha3 = geTReal( "RecoveryParameterAlpha3", recoverypara );
	ensRecoveryModel.NesKappa2 = geTReal( "RecoveryParameterKappa2", recoverypara );
	ensRecoveryModel.NesC3 = geTReal( "RecoveryParameterC3", recoverypara );
	ensRecoveryModel.NesC4 = geTReal( "RecoveryParameterC4", recoverypara );
	ensRecoveryModel.NesC5 = geTReal( "RecoveryParameterC5", recoverypara );
	//##MK::implement consistency check

	delete recoverypara;

	return status;
}

bool ensembleHdl::init_parameter_drag( void ) {
	bool status = true;

	dataBlockP dragpara = readDataBlock("DragParameter", score_input );

	ensDragData.ZenerConsider = DISPERSOIDDRAG_NO;
	int zenermodel = geTInt( "ZenerConsider", dragpara );
	if ( zenermodel == 1 ) ensDragData.ZenerConsider = DISPERSOIDDRAG_CONSTANT;
	if ( zenermodel == 2 ) ensDragData.ZenerConsider = DISPERSOIDDRAG_TIMEDEP;
	ensDragData.ZenerAlpha = geTReal( "ZenerAlpha", dragpara );
	ensDragData.ZenerGamma = geTReal( "ZenerGamma", dragpara );

	delete dragpara;

	dataBlockP microchemevo = readDataBlock( "EvoMicrochemistry" , score_input );

	dataLineP mline = microchemevo->first;
	long nmlines = microchemevo->lineCount;
	if ( myRank == MASTER && nmlines < 2 ) { cout << "ERROR::An insufficient number (<2) of Particle Time Dispersion data have been provided to perform the simulation!" << endl; }

	QUICKASSERT( nmlines <= DEFAULT_MAX_INPUTDATALEN );

	double ti = 0.0;
	double fri = 0.0;

	for ( long ts = 0; ts < nmlines; ts++ ) {
		ti = getReal(mline, 1);
		fri = getReal(mline, 2);

		dispersoidtime.push_back( ti );
		dispersoidfr.push_back( fri );

#ifdef REPORTSTYLE_DEVELOPER
		if ( myRank == MASTER ) cout << "Node;" << this->myRank << ";ti;fri;" << ti << ";" << fri << "in the vector instead Time[ts];Dispersion[ts];" << "\ts=" << ts << "\t" << dispersoidtime[ts] << ";" << dispersoidfr[ts] << endl;
#endif
		mline = mline->next;
	}

	delete mline;
	delete microchemevo;
	if ( dispersoidtime.size() != dispersoidfr.size() ) { status = false; cout << "ERROR::Inconsistency in reading microchemistry evolution!" << endl; }

	return status;
}

bool ensembleHdl::init_parameter_idealtexture( void ) {
	bool status = true;

	//import ideal texture components
	dataBlockP idealcomp = readDataBlock( "IdealComponents", score_input);
	dataLineP iline = idealcomp->first;
	long nilines = idealcomp->lineCount;

	QUICKASSERT( nilines <= DEFAULT_MAX_INPUTDATALEN );

	double bunge_ori[3], quaternion_ori[4], scatt;

	for ( long ii = 0; ii < nilines; ii++ ) {
		struct ideal icomp;

		bunge_ori[0] = DEG2RAD(getReal( iline, 2 ));
		bunge_ori[1] = DEG2RAD(getReal( iline, 3 ));
		bunge_ori[2] = DEG2RAD(getReal( iline, 4 ));
		scatt  = DEG2RAD(getReal( iline, 5 ));

		if ( bunge_ori[0] < BUNGE_MIN || bunge_ori[0] > BUNGE_MAX || bunge_ori[1] < BUNGE_MIN || bunge_ori[1] > BUNGE_MAX || bunge_ori[2] < BUNGE_MIN || bunge_ori[2] > BUNGE_MAX ) { status = false; cout << "ERROR::Invalid Euler angles in ideal orientation defintion!" << endl; break; } 

		ensmath.euler2quat( bunge_ori, quaternion_ori ); //interally everything with quaternions

		icomp.bunge1 = bunge_ori[0];
		icomp.bunge2 = bunge_ori[1];
		icomp.bunge3 = bunge_ori[2];
		icomp.q0 = quaternion_ori[0];
		icomp.q1 = quaternion_ori[1];
		icomp.q2 = quaternion_ori[2];
		icomp.q3 = quaternion_ori[3];
		icomp.scatter = scatt;

		standardlagen.push_back( icomp );

#ifdef REPORTSTYLE_DEVELOPER
			cout << "Idealcomponent imported::bunge123;q0123;scatt;" << standardlagen[ii].bunge1 << ";" << standardlagen[ii].bunge2 << ";" << standardlagen[ii].bunge3 << ";" << standardlagen[ii].q0 << ";" << standardlagen[ii].q1 << ";" << standardlagen[ii].q2 << ";"<< standardlagen[ii].q3 << ";" << standardlagen[ii].scatter << endl;
#endif

		iline = iline->next;
	}
	delete iline;
	delete idealcomp;

	return status;
}

bool ensembleHdl::init_parameter_processing( void ) {
	bool status = true;

	dataBlockP timetemp = readDataBlock( "AnnealingSchedule", score_input );
	dataLineP tline = timetemp->first;
	long ntlines = timetemp->lineCount;
	if ( myRank == MASTER && ntlines < 2 ) { status = false; cout << "ERROR::An insufficient number (<2) of Time-Temperature data have been provided to perform the simulation!" << endl; }

	QUICKASSERT( ntlines <= DEFAULT_MAX_INPUTDATALEN );

	double ti_old = -1.0;
	double ti_now = 0.0;
	double Ti = (25.0 + TOFFSET);

	for ( long ts = 0; ts < ntlines; ts++ ) {
		ti_now = getReal(tline, 1);
		Ti = getReal(tline, 2);
		Ti += TOFFSET;

		if ( ti_now < 0.0 || Ti < 0.0 || ti_old >= ti_now ) { status = false; cout << "ERROR::Invalid annealing schedule data (<0.0 or not consecutive in time)!" << endl; break; }
		ti_old = ti_now;

		ttime.push_back( ti_now );
		ttemperature.push_back( Ti );

#ifdef REPORTSTYLE_DEVELOPER
			cout << "Node;" << this->myRank << ";ti;Ti;" << ti_now << ";" << Ti << "in the vector instead ttime[ts];ttemperature[ts];" << "\ts=" << ts << "\t" << ttime[ts] << ";" << ttemperature[ts] << endl;
#endif

		tline = tline->next;
	}

	delete tline;
	delete timetemp;
	if ( ttime.size() != ttemperature.size() ) { status = false; cout << "ERROR::Inconsistency in readingf the annealing schedule!" << endl; }
	if ( ttime.size() < 2 ) { status = false; cout << "ERROR::Annealing schedule is empty!" << endl; }

	return status;
}

bool ensembleHdl::init_parameter_defgpool( void ) {
	bool status = true;

	dataBlockP defgpara = readDataBlock("DeformedGrainsPool", score_input);
	dataLineP dline = defgpara->first;
	long ndefgpool = defgpara->lineCount;

	if ( ndefgpool > DEFAULT_MAX_DEFORMEDGRAINS ) { 
		cout << "WARNING::The pool of deformed grains is too large, only the first " << DEFAULT_MAX_DEFORMEDGRAINS << " are read." << endl;
		ndefgpool = DEFAULT_MAX_DEFORMEDGRAINS;
	}

	worldoripool.reserve(ndefgpool);
	worlddefgpool.reserve(ndefgpool);

	uint32_t id = 0;
	double bunge_ori[3] = {0.0, 0.0, 0.0};

	for (long dg = 0; dg < ndefgpool; dg++) {
		struct defg agrain;

		bunge_ori[0] = DEG2RAD(getReal(dline, 1));
		bunge_ori[1] = DEG2RAD(getReal(dline, 2));
		bunge_ori[2] = DEG2RAD(getReal(dline, 3));

		if ( bunge_ori[0] < BUNGE_MIN || bunge_ori[0] > BUNGE_MAX || bunge_ori[1] < BUNGE_MIN || bunge_ori[1] > BUNGE_MAX || bunge_ori[2] < BUNGE_MIN || bunge_ori[2] > BUNGE_MAX ) { status = false; cout << "ERROR::Invalid Euler angles in deformed grain orientation defintion!" << endl; break; } 

		id = ens_check_disjunctness_io( bunge_ori );

		agrain.ori = id;
		agrain.rho0 = getReal(dline, 4);
		agrain.rho = agrain.rho0;

		if ( agrain.rho <= RHO_RECRYSTALLIZED_MATERIAL ) { status = false; cout << "ERROR::Unphysical low dislocation density for deformed material!" << endl; break; }

		worlddefgpool.push_back(agrain);

		//cout << "dg;ori;e1;e2;e3;q0;q1;q2;q3;rho0;rho " << dg << ";" << worlddefgpool[dg].ori << ";" << worldoripool[id].bunge1 << ";" << worldoripool[id].bunge2 << ";" << worldoripool[id].bunge3 << ";" << worldoripool[id].q0 << ";" << worldoripool[id].q1 << ";" << worldoripool[id].q2 << ";" << worldoripool[id].q3 << ";" << worlddefgpool[dg].rho0 << ";" << worlddefgpool[dg].rho << endl;

		dline = dline->next;
	} // for all grains in the pool

	delete defgpara;

	return status;
}

bool ensembleHdl::init_parameter_rxgpool( void ) {
	bool status = true;

	dataBlockP rxgpara = readDataBlock("RXGrainsPool", score_input);
	dataLineP rline = rxgpara->first;
	long nrxgpool = rxgpara->lineCount;

	if ( nrxgpool > DEFAULT_MAX_RXGRAINS ) {
		cout << "WARNING::The pool of recrystallized grains is too large, only the first " << DEFAULT_MAX_RXGRAINS << " are read." << endl;
		nrxgpool = DEFAULT_MAX_RXGRAINS;
	}
	//there are already many orientations in the oripool most likely some are not significantly disjunct from the already existing ones
	worldrxgpool.reserve(nrxgpool);

	uint32_t id = 0;
	double bunge_ori[3] = {0.0, 0.0, 0.0};

	for ( long rg = 0; rg < nrxgpool; rg++) {
		struct rxg agrain;

		bunge_ori[0] = DEG2RAD(getReal(rline, 1));
		bunge_ori[1] = DEG2RAD(getReal(rline, 2));
		bunge_ori[2] = DEG2RAD(getReal(rline, 3));

		if ( bunge_ori[0] < BUNGE_MIN || bunge_ori[0] > BUNGE_MAX || bunge_ori[1] < BUNGE_MIN || bunge_ori[1] > BUNGE_MAX || bunge_ori[2] < BUNGE_MIN || bunge_ori[2] > BUNGE_MAX ) { status = false; cout << "ERROR::Invalid Euler angles in recrystallized orientation defintion!" << endl; break; } 

		id = ens_check_disjunctness_io( bunge_ori );

		agrain.ori = id;
		agrain.tincub = getReal(rline, 4);

		if ( agrain.tincub < 0.0 ) { status = false; cout << "ERROR::Invalid incubation time!" << endl; break; }

		worldrxgpool.push_back(agrain);

		//cout << "rg;ori;e1;e2;e3;q0;q1;q2;q3;tincub;" << rg << ";" << worldrxgpool[rg].ori << ";" << worldoripool[id].bunge1 << ";" << worldoripool[id].bunge2 << ";" << worldoripool[id].bunge3 << ";" << worldoripool[id].q0 << ";" << worldoripool[id].q1 << ";" << worldoripool[id].q2 << ";" << worldoripool[id].q3 << ";" << worldrxgpool[rg].tincub << endl;

		rline = rline->next;
	} // for all grains in the pool

	delete rxgpara;

	return status;
}

bool ensembleHdl::init_parameter_additional( void ) {
	bool status = true;

	if ( outopt_rendermethod == RENDERING_MS2D || outopt_rendermethod == RENDERING_MS3D || outopt_localrenderboundaries == RENDERING_BOUNDARIES_YES ) {
		dataBlockP renderlogpoints = readDataBlock("RenderMicrostructure", score_input );

		dataLineP rline = renderlogpoints->first;
		long nrlines = renderlogpoints->lineCount;

		QUICKASSERT( nrlines < DEFAULT_MAX_INPUTDATALEN );

		double Xri = 0.0;
		for ( long r = 0; r < nrlines; r++ ) { 
			Xri = getReal( rline, 2 );
			if ( Xri < 0.0 || Xri > 1.0 ) { status = false; cout << "ERROR::Invalid recrystallized volume fraction X in RenderMicrostructure!" << endl; break; }

			UserDefLogPoint_MS_Rendering.push_back( Xri );

			rline = rline->next;
		}
		delete rline;
		delete renderlogpoints;

		std::sort ( UserDefLogPoint_MS_Rendering.begin(), UserDefLogPoint_MS_Rendering.end(), SortDoubleAscending );

		//import which simulations to render
		dataBlockP candidatelist = readDataBlock("RenderOnlyTheseSolitaryUnits", score_input );

		dataLineP cline = candidatelist->first;
		long nclines = candidatelist->lineCount;
		int id = 0;

//cout << myRank << "Render only these regions = " << nclines << ";" << cline << ";" << id << endl;
		for ( long r = 0; r < nclines; r++ ) { 
			id = getInt ( cline, 1 );
			if ( id < 0 || id >= nworldCAs ) { status = false; cout << "ERROR::Invalid simulation index listed that should become rendered!" << endl; break; }
//cout << "\t\t" << r << ";" <<  myRank << "Render only these regions = " << nclines << ";" << cline << ";" << id << endl;

			UserDefLogPoint_WhichCAtoOutput.push_back( (uint32_t) id );

			cline = cline->next;
		}
		delete cline;
		delete candidatelist;

		std::sort( UserDefLogPoint_WhichCAtoOutput.begin(), UserDefLogPoint_WhichCAtoOutput.end(), SortIntAscending );

		//import which ZSections to render
		dataBlockP zpositions = readDataBlock("RenderTheseSections", score_input );

		dataLineP zline = candidatelist->first;
		long nzlines = candidatelist->lineCount;
		double zid = 0;
		for ( long z = 0; z < nzlines; z++ ) {
			zid = getReal( zline, 2 );
			if ( zid < (0.0 - 1.0e-5) || zid > (1.0 + 1.0e-5) ) { status = false; cout << "ERROR::Invalid relative zsection coordinate provided in RenderTheseSections!" << endl; break; }

			UserDefLogPoint_MS_RenderZPos.push_back( zid );

			zline = zline->next;
		}
		delete zline;
		delete zpositions;

		std::sort( UserDefLogPoint_MS_RenderZPos.begin(), UserDefLogPoint_MS_RenderZPos.end(), SortDoubleAscending );
	}

	//##MK::remove duplicates! currently they will simply be overwritten as they get the same filename....


	//import defragmentation points
	dataBlockP defragpoints = readDataBlock( "HeuristicRXFrontListDefragmentation", score_input);
	dataLineP dline = defragpoints->first;
	long ndlines = defragpoints->lineCount;

	QUICKASSERT( ndlines < DEFAULT_MAX_INPUTDATALEN );

	double Xdi = 0.0;
	for ( long d = 0; d < ndlines; d++ ) { 
		Xdi = getReal( dline, 2 );
		if ( Xdi < 0.0 || Xdi > 1.0 ) { status = false; cout << "ERROR::Invalid recrystallized volume fraction X at which to defragment list!" << endl; break; }

		UserDefLogPoint_X_CellListDefragment.push_back( Xdi );

		dline = dline->next;
	}
	delete dline;
	delete defragpoints;

	std::sort ( UserDefLogPoint_X_CellListDefragment.begin(), UserDefLogPoint_X_CellListDefragment.end(), SortDoubleAscending );

	//import logpoints that are interpreted locally in each caHdl dependent on the local recrystallized fraction
	dataBlockP logpoints = readDataBlock( "SingleGrainVolOverTime", score_input);
	dataLineP pline = logpoints->first;
	long nplines = logpoints->lineCount;

	QUICKASSERT( nplines < DEFAULT_MAX_INPUTDATALEN );

	double Xi = 0.0;
	for ( long pp = 0; pp < nplines; pp++ ) { 
		Xi = getReal( pline, 2 );
		if ( Xi < 0.0 || Xi > 1.0 ) { status = false; cout << "ERROR::Invalid recrystallized volume fraction X at which to log locally all grains!" << endl; break; }

		UserDefLogPoint_X_Output.push_back( Xi ); 	//no reserve because single page already allocated

		pline = pline->next;
	}
	delete pline;
	delete logpoints;

	std::sort ( UserDefLogPoint_X_Output.begin(), UserDefLogPoint_X_Output.end(), SortDoubleAscending );

	return status;
}

bool ensembleHdl::init_ebsdmap_read( const char* ebsdfname ) {
	//reads three headered EBSD file of format x,y,gid,phi1,Phi,phi2, kam, optional RGB
	std::cout << "Reading SEM/EBSD data sequentially..." << std::endl;

	ifstream ebsdfile;
	string ebsdline;
	istringstream line;
	string datapiece;

	ebsdfile.open( ebsdfname );

	//map<uint32_t,uint32_t> UniqueGrainIDs;

	if ( ebsdfile.is_open() == true ) {
		//inspect header bool inheader = true; char keychar[1] = { char(35) }; //"#" character

		//kick three header lines
		getline( ebsdfile, ebsdline );
		getline( ebsdfile, ebsdline );
		getline( ebsdfile, ebsdline );

		//read data
		while ( ebsdfile.good() == true ) {
			getline( ebsdfile, ebsdline );
			if ( ebsdline.size() > 6 ) { //required when attempting to read x, y, gid, bunge1,2,3, kam
				istringstream line( ebsdline );

				struct ebsdpoint ap;
				getline(line, datapiece, '\t');	ap.x = fabs(atof( datapiece.c_str() )); //##MK::experimental x and y flipped
				getline(line, datapiece, '\t');	ap.y = fabs(atof( datapiece.c_str() )); //fabs because sometimes in Matlab the sign bit becomes set...
				getline(line, datapiece, '\t');	ap.gid = atoi( datapiece.c_str() ); //##MK::assumes indices to be positive on [0,2^32-1]
				//##MK::IDs need no relabeling SCORE will do so...

				double bunge_ori[3]; //, bunge_mod[3]; //mean ori computed from all measurement points inside the grain
				getline(line, datapiece, '\t'); bunge_ori[0] = atof( datapiece.c_str() ); //MK::*.ang file is already in rad
				getline(line, datapiece, '\t');	bunge_ori[1] = atof( datapiece.c_str() );
				getline(line, datapiece, '\t');	bunge_ori[2] = atof( datapiece.c_str() );

				//ebsdori2score( bunge_ori, modifiedbunge );
				ap.bunge1 = bunge_ori[0];
				ap.bunge2 = bunge_ori[1];
				ap.bunge3 = bunge_ori[2];

				double quaternion_mod[4];
				ensmath.euler2quat( bunge_ori, quaternion_mod );
				ap.q0 = quaternion_mod[0];
				ap.q1 = quaternion_mod[1];
				ap.q2 = quaternion_mod[2];
				ap.q3 = quaternion_mod[3];

				//getline(line, datapiece, '\t');		ap.iq = 0.0; //##MK::overwritten atof( datapiece.c_str() );
				//getline(line, datapiece, '\t');		ap.ci = 0.0; //##MK::overwritten atof( datapiece.c_str() );

				getline(line, datapiece, '\t');	ap.kam = DEG2RAD(atof( datapiece.c_str() )); //KAM needs to be in radiant

				/*if ( ebsdline.size() > 6+3 ) { //read optional RGB color code for the orientation
					int val[3] = {0,0,0};
					getline(line, datapiece, '\t');	val[RED] = atoi( datapiece.c_str() );
					getline(line, datapiece, '\t'); val[GREEN] = atoi( datapiece.c_str() );
					getline(line, datapiece, '\t'); val[BLUE] = atoi( datapiece.c_str() );
					if (val[RED] < 0 || val[RED] > RGB_MAX) { std::cout << "ERROR::Invalid color value for RED measurement point " << expData.size() << endl; return false; }
					if (val[GREEN] < 0 || val[GREEN] > RGB_MAX) { std::cout << "ERROR::Invalid color value for GREEN measurement point " << expData.size() << endl; return false; }
					if (val[BLUE] < 0 || val[BLUE] > RGB_MAX) { std::cout << "ERROR::Invalid color value for BLUE measurement point " << expData.size() << endl; return false; }
					ap.r = val[RED];
					ap.g = val[GREEN];
					ap.b = val[BLUE];
					//cout << ap.r << ";" << ap.g << ";" << ap.b << endl;
				}*/

				//collect measurement point
				expData.push_back( ap );

				//MK::identify disjoint grain ids but compute pixel averaged quantities of these later, hence scanning of shorter arrays
				auto it = expUniqueGrainIDs.find( ap.gid );
				if ( it != expUniqueGrainIDs.end() ) { //if found, the most likely case the more grains have already been added as NumberOfUniqueGrains << NumberOfPixel
					continue;
				}
				else {
					expUniqueGrainIDs.insert( pair<uint32_t,uint32_t>(ap.gid, ap.gid) );
					//expGrainIDs.push_back( ap.gid );
					cout << "SEM/EBSD grain id/phi1/Phi/phi2//KAM = " << ap.gid << ";" << ap.bunge1 << ";" << ap.bunge2 << ";" << ap.bunge3 << "\t\t" << ap.kam << "\t\t" << ap.q0 << ";" << ap.q1 << ";" << ap.q2 << ";" << ap.q3 << endl;
				}

				/*
				if ( expGrainIDs.size() > 0 ) {
					bool found = false;
					for ( uint32_t g = 0; g < expGrainIDs.size(); g++ ) {
						//if ( expGrainIDs[g] == ap.gid ) { found = true; break; }
						if ( expGrainIDs[g] != ap.gid )
							continue;

						//expGrainIDs[g] == ap.gid, already existent...
						found = true;
						break;
					}

					if ( found == false ) { //only collect unknown userGrainID
						expGrainIDs.push_back( ap.gid );
cout << "SEM/EBSD grain id/phi1/Phi/phi2//KAM = " << ap.gid << ";" << ap.bunge1 << ";" << ap.bunge2 << ";" << ap.bunge3 << "\t\t" << ap.kam << "\t\t" << ap.q0 << ";" << ap.q1 << ";" << ap.q2 << ";" << ap.q3 << endl;
					}
				}
				else { //always add the first userGrainID
					expGrainIDs.push_back( ap.gid );
cout << "SEM/EBSD grain id/phi1/Phi/phi2//KAM = " << ap.gid << ";" << ap.bunge1 << ";" << ap.bunge2 << ";" << ap.bunge3 << "\t\t" << ap.kam << "\t\t" << ap.q0 << ";" << ap.q1 << ";" << ap.q2 << ";" << ap.q3 << endl;
				}
				*/
			}
		}

		ebsdfile.close();
		cout << "Reading of " << ebsdfname << " was successful with npx/nids = " << expData.size() << ";" << expUniqueGrainIDs.size() << endl; //expGrainIDs.size() << endl;
		return true;
	}
	else { cout << "File " << ebsdfname << " either not existent or cannot be opened!" << endl; return false; }
}

bool ensembleHdl::init_parameter_defgpool_exp( void )
{
	if ( expUniqueGrainIDs.size() > DEFAULT_MAX_DEFORMEDGRAINS ) { cout << "WARNING::The pool of EBSD grains is very large." << endl; }

	//worldoripool.reserve(ndefgpool);
	//worlddefgpool.reserve(ndefgpool);
	/*std::sort( expGrainIDs.begin(), expGrainIDs.end(), SortIntAscending );*/
	//##MK::expUniqueGrainIDs is a map of uint32_t it keeps automatic track of an ascendingly order set of uint32_t keys

	//grain indexing in EBSD mapping may start with any number and be non-contiguous, so assign contiguous labels on [1,N]

	/*vector<uint32_t>* newlabels = NULL;
	newlabels = new vector<uint32_t>;
	uint32_t label = 1;
	for ( size_t g = 0; g < expGrainIDs.size(); g++ ) { //define new labels
		newlabels->push_back( label );
		label++;
	}*/
	uint32_t label = 1;
	for ( auto it = expUniqueGrainIDs.begin(); it != expUniqueGrainIDs.end(); ++it) {
		it->second = label;
		label++; //newlabel is now a hash for the new labels to assign the grains with old labels in expGrainIDs
	}
	//QUICKASSERT( label == expGrainIDs.size() );
	for ( size_t px = 0; px < expData.size(); px++ ) { //relabeling
		uint32_t old_id = expData[px].gid;
		uint32_t pos = UINT32T_MAX;
		auto it = expUniqueGrainIDs.find( old_id ); //first keeps the old ID second has the new ID
		if ( it != expUniqueGrainIDs.end() )
			expData[px].gid = it->second; //if assign the newlabel
		else {
			cerr << "Labeling inconsistence in init_parameter_defgpool" << endl;
			return false;
		}
		/*for ( size_t i = 0; i < expGrainIDs.size(); i++ ) {
			if ( expGrainIDs[i] != old_id ) 
				continue;
			//else
			pos = i; break;
		}
		QUICKASSERT( pos != UINT32T_MAX );
		expData[px].gid = newlabels->at(pos);*/
	}
	/*//now all measured points in the mapping have labels from [1,label)
	for ( size_t g = 0; g < expGrainIDs.size(); g++ ) { //relabel also in expGrainIDs
		expGrainIDs[g] = newlabels->at(g);
	}
	delete newlabels; newlabels = NULL;*/

	//uint32_t largestid = expGrainIDs[expGrainIDs.size()-1];
	uint32_t largestid = label-1;
	//QUICKASSERT( largestid == label );
	if ( largestid > DEFAULT_MAX_DEFORMEDGRAINS_EBSD ) { cout << "ERROR::The pool of EBSD grains is too large!" << endl; return false; }

	//build hashtable-type summation arrays for utilization in mean computations
	vector<uint32_t>* cnt = NULL;
	cnt = new vector<uint32_t>;
	vector<double>* xyzsum = NULL;
	xyzsum = new vector<double>;
	vector<double>* kamsum = NULL;
	kamsum = new vector<double>;
	vector<euler313>* qb = NULL;
	qb = new vector<euler313>;
	//vector<quat>* qori = NULL;
	//qori = new vector<quat>;
	//vector<vector<quat>*>* qsum = NULL;
	//qsum = new vector<vector<quat>*>;

	for ( size_t hk = 0; hk < (1+largestid); hk++ ) {
		cnt->push_back(0);
		xyzsum->push_back(0.0); xyzsum->push_back(0.0);
		kamsum->push_back(0.0);
		struct euler313 el;
		qb->push_back(el);
	}
	//qori->reserve(1+largestid);
	//struct quat aq;
	//for ( size_t hk = 0; hk < (1+largestid); hk++ )
	//	qori->push_back(aq);
	//qsum->reserve(1+largestid);
	//for ( size_t hk = 0; hk < (1+largestid); hk++ ) 
	//	qsum->push_back( NULL );
	//kamsum->reserve(1+largestid);
	//std::fill (kamsum->begin(), kamsum->end(), 0.0 );

	//MK::accumulate pixel values before averaging profits for thresholded structures as there are usually an order of magnitude fewer grains than px in an EBSD mapping
	uint32_t hashkey;
	for ( size_t px = 0; px < expData.size(); px++ ) {
		hashkey = expData[px].gid;
		cnt->at(hashkey) += 1;
		xyzsum->at(2*hashkey+0) += expData[px].x;
		xyzsum->at(2*hashkey+1) += expData[px].y;
		//xyzsum->at(3*hashkey+2) += 0.0;
		kamsum->at(hashkey) += expData[px].kam;
		//collect the quaternions of grain hashkey's pixel
		//if ( qsum->at(hashkey) == NULL )
		//	qsum->at(hashkey) = new vector<quat>;

		//double pxbunge[3] = { expData[px].bunge1, expData[px].bunge2, expData[px].bunge3 };
		//double pxquaternion[4];
		//ensmath.euler2quat( pxbunge, pxquaternion );
		qb->at(hashkey).bunge1 = expData[px].bunge1;
		qb->at(hashkey).bunge2 = expData[px].bunge2;
		qb->at(hashkey).bunge3 = expData[px].bunge3;

		//qori->at(hashkey).q0 = pxquaternion[0];
		//qori->at(hashkey).q1 = pxquaternion[1];
		//qori->at(hashkey).q2 = pxquaternion[2];
		//qori->at(hashkey).q3 = pxquaternion[3];
		//q.q0 = pxquaternion[0];
		//q.q1 = pxquaternion[1];
		//q.q2 = pxquaternion[2];
		//q.q3 = pxquaternion[3];
		//qsum->at(hashkey)->push_back( q );
	}

	/*
	//all grains consistent,i.e. cnts > 0?
	for ( size_t id = 0; id < expGrainIDs.size(); id++ ) { //imprint ascending order of ids
		if ( cnt->at(hashkey) == 0 ) {
			cout << "ERROR::Inconsistency! Grain without associated pixel!" << endl;
			delete cnt;
			delete xyzsum;
			delete kamsum;
			delete qb;
			//delete qori;
			//for ( size_t hk = 0; hk < (1+largestid); hk++ ) 
			//	if ( qsum->at(hk) != NULL ) 
			//		delete qsum->at(hk);
			//delete qsum;
			return false;
		}
	}
	*/

	//perform averaging process
	/*for ( size_t id = 0; id < expGrainIDs.size(); id++ ) { //create thresholded grains from EBSD mapping only for existent IDs*/
	for ( auto it = expUniqueGrainIDs.begin(); it != expUniqueGrainIDs.end(); it++ ) {
		/*hashkey = expGrainIDs.at(id);*/
		hashkey = it->second; //the new_id !
		struct ebsdgrain agr;
		agr.baryx = xyzsum->at(2*hashkey+0) / cnt->at(hashkey);
		agr.baryy = xyzsum->at(2*hashkey+1) / cnt->at(hashkey);
		//agr.baryz = xyzsum->at(3*hashkey+2) / cnt->at(hashkey);
		agr.kam = kamsum->at(hashkey) / cnt->at(hashkey);
		agr.gid = hashkey; //expGrains get same label as EBSD pixel have, labels on [1,N]

		//double meanquat[4] = { qori->at(hashkey).q0, qori->at(hashkey).q1, qori->at(hashkey).q2, qori->at(hashkey).q3 };
		//double meanbunge[3] = {0.0, 0.0, 0.0};
		//ensmath.quat2euler( meanquat, meanbunge );
		//quaternion_averaging( qsum->at(hashkey), meanbunge );
		agr.bunge1 = qb->at(hashkey).bunge1; //meanbunge[0];
		agr.bunge2 = qb->at(hashkey).bunge2; //meanbunge[1];
		agr.bunge3 = qb->at(hashkey).bunge3; //meanbunge[2];

		//convert Bunge Euler to quaternion only once!
		double bunge_ori[3] = { agr.bunge1, agr.bunge2, agr.bunge3 };
		double quaternion_ori[4] = {1.0, 0.0, 0.0, 0.0};
		ensmath.euler2quat( bunge_ori, quaternion_ori );

		agr.q0 = quaternion_ori[0];
		agr.q1 = quaternion_ori[1];
		agr.q2 = quaternion_ori[2];
		agr.q3 = quaternion_ori[3];

		expGrains.push_back( agr ); //organization of expGrains is linear expGrains[].gid runs from [1,N] ascendingly

cout << "Internal grain old_id/new_id phi1/Phi/phi2//KAM " << it->first << ";" << it->second << "\t\t" << agr.gid << ";" << agr.bunge1 << ";" << agr.bunge2 << ";" << agr.bunge3 << ";" << agr.kam << endl;
	}

	//as in Martin Diehl's project grains will be generated based on the expGrains we have to add to worlddefgpool and worldoripool later

	//tidy up
	delete cnt;
	delete xyzsum;
	delete kamsum;
	delete qb;
	//for ( size_t hk = 0; hk < (1+largestid); hk++ ) 
	//	if ( qsum->at(hk) != NULL ) 
	//		delete qsum->at(hk);
	//delete qsum;

	return true;
}

bool ensembleHdl::init_parameter_rxgpool_exp( void )
{
	//in Martin Diehl's project nuclei will be generated from the structure so no importing of a priori information on nucleation spectrum required
	return true;
}

void ensembleHdl::report_successful_setup( void )
{
	cout << "\t\tEnsemble size " << nworldCAs << "\n";
	cout << "\t\tSU-Geometry " << ensCAGeometry.nboxedge_rd << ";" << ensCAGeometry.nboxedge_td << ";" << ensCAGeometry.nboxedge_nd << " is RD x ND x TD CA geometry with " << ensCAGeometry.cellsize << " cellsize" << " m\n";
	cout << "\t\tGrain-Geometry " << ensPhysData.defgmean_rd << "," << ensPhysData.defgmean_td << "," << ensPhysData.defgmean_nd << " m\n";
	cout << "\t\tStandardlagen " << standardlagen.size() << "\n";
	cout << "\t\tDefgpool " << this->worlddefgpool.size() << "\n";
	cout << "\t\tRxgpool " << this->worldrxgpool.size() << "\n";
	cout << "\t\tOripool " << this->worldoripool.size() << endl;

	//cout << UserDefLogPoint_X_CellListDefragment.size() << ", logpoints for RXFront defragmentation were initialized as desired and sorted ascending." << endl;
	//cout << UserDefLogPoint_X_Output.size() << ", logpoints for X were initialized as desired and sorted in ascending order." << endl;
	//cout << XMAX << ", is the upper limit for RXFRACTION numerical integration scheme." << endl;
	//cout << TMAX << ", is the upper limit for REALTIME numerical integration scheme." << endl;
	//cout << initialRelCellCaching << ", is the initial recaching value." << endl;
	//cout << ensPhysData.G0 << "," << ensPhysData.bZeroCelsius << ", ShearModulus/Burgersvector." << endl;
	//cout << ensPhysData.LAGBm0 << "," << ensPhysData.HAGBm0 << "," << ensPhysData.GSm0 << ",  LAGB/HAGB/GS m0." << endl;
	//cout << ensPhysData.LAGBHact << "," << ensPhysData.HAGBHact << "," << ensPhysData.GSHact << ", LAGB/HAGB/GS Hact=" << endl;
}

void ensembleHdl::report_datastructures( void ) 
{
	if ( myRank == MASTER ) { 
		cout << "The size of the SCORE datastructures is as follows..." << endl;
		cout << "\t\tdefcell=" << sizeof(defcell) << endl;
		cout << "\t\tcell=" << sizeof(cell) << endl;
		cout << "\t\trxg=" << sizeof(rxg) << endl;
		cout << "\t\tdefg=" << sizeof(defg) << endl;
		cout << "\t\tori=" << sizeof(ori) << endl;
		cout << "\t\tqsymmfcc=" << sizeof(qsymm_fcc) << endl;
		cout << "\t\tcadefg=" << sizeof(cadefg) << endl;
		cout << "\t\tcarxg=" << sizeof(carxg) << endl;
		cout << "\t\tideal=" << sizeof(ideal) << endl;
		cout << "\t\tloginfophys=" << sizeof(loginfo_ca_physics) << endl;
		cout << "\t\tdefgseed=" << sizeof(defgseed) << endl;
		cout << "\t\tcellsBndFast=" << sizeof(cellsBndFast) << endl;
		cout << "\t\tbndFaceFast=" << sizeof(bndFaceFast) << endl;
		cout << "\t\tbndColumnFast=" << sizeof(bndColumnFast) << endl;
		cout << "\t\tloginfogrevo=" << sizeof(loginfo_grainevo_ca) << endl;
		cout << "\t\tloginfo_rxfront=" << sizeof(loginfo_rxfrontstats_ca) << endl;
		cout << "\t\tensHdl=" << sizeof(ensembleHdl) << endl;
		cout << "\t\tcaHdl=" << sizeof(caHdl) << endl;
	}
}


void ensembleHdl::compare_two_ang_files( const char* inputang, const char* snpang )
{
	//debug: stupidly compare disorientation between two Bunge triplets from a 1.) SCORE*.ang file against a 2.) snapshot *.ang file line1-by-line1
	ifstream scoreangfile, snpangfile;
	string scoreangline, snpangline;
	istringstream line1, line2;
	string dp1, dp2;

	scoreangfile.open( inputang );
	snpangfile.open( snpang );
	if ( scoreangfile.is_open() == true && snpangfile.is_open() == true ) {
		for ( unsigned int h1 = 0; h1 < 3; h1++ ) //kick three header lines from inputang
			getline( scoreangfile, scoreangline );
		for ( unsigned int h2 = 0; h2 < 6; h2++ ) 
			getline ( snpangfile, snpangline );

		double bunge1[3] = {0.0, 0.0, 0.0};
		double quat1[4] = {1.0, 0.0, 0.0, 0.0};

		double bunge2[3] = {0.0, 0.0, 0.0};
		double quat2[4] = {1.0, 0.0, 0.0, 0.0};

		for ( unsigned int px = 0; px < 526680; px++ ) { //##hacky
			//better reset
			bunge1[0] = 0.0; bunge1[1] = 0.0; bunge1[2] = 0.0;
			quat1[0] = 1.0; quat1[1] = 0.0; quat1[2] = 0.0; quat1[3] = 0.0;
			bunge2[0] = 0.0; bunge2[1] = 0.0; bunge2[2] = 0.0;
			quat2[0] = 1.0; quat2[1] = 0.0; quat2[2] = 0.0; quat2[3] = 0.0;

			getline( scoreangfile, scoreangline );
			if ( scoreangline.size() > 6 ) { //required when attempting to read x, y, gid, bunge1,2,3, kam
				istringstream line1( scoreangline );

				getline(line1, dp1, '\t'); getline(line1, dp1, '\t'); getline(line1, dp1, '\t');
				getline(line1, dp1, '\t');		bunge1[0] = atof( dp1.c_str() );
				getline(line1, dp1, '\t');		bunge1[1] = atof( dp1.c_str() );
				getline(line1, dp1, '\t');		bunge1[2] = atof( dp1.c_str() );
				getline(line1, dp1, '\t');
			}

			getline ( snpangfile, snpangline );
			if ( snpangline.size() > 3 ) {
				istringstream line2( snpangline );
				getline(line2, dp2, ';');		bunge2[0] = atof( dp2.c_str() );
				getline(line2, dp2, ';');		bunge2[1] = atof( dp2.c_str() );
				getline(line2, dp2, ';');		bunge2[2] = atof( dp2.c_str() );
			}

			ensmath.euler2quat( bunge1, quat1 );
			ensmath.euler2quat( bunge2, quat2 );

			//compare
			/*double mis12 = (this->disori_angle_oriquat_cubic( quat1, quat2 ) * 180.0 / _PI_);
			double mis21 = (this->disori_angle_oriquat_cubic( quat2, quat1 ) * 180.0 / _PI_);
			//if ( mis12 < 1.0e-5 ) 
			//	continue;

			//compare consistence of back-conversion
			double bunge1r[3];
			double bunge2r[3];
			ensmath.quat2euler( quat1, bunge1r );
			ensmath.quat2euler( quat2, bunge2r );

			double d1 = fabs(bunge1[0]-bunge1r[0]) + fabs(bunge1[1]-bunge1r[1]) + fabs(bunge1[2]-bunge1r[2]);
			double d2 = fabs(bunge2[0]-bunge2r[0]) + fabs(bunge2[1]-bunge2r[1]) + fabs(bunge2[2]-bunge2r[2]);
			*/

			//get passive interpretation of presumably actively defined orientation quaternion
			double quat1pas[4] = { quat1[0], -1.0*quat1[1], -1.0*quat1[2], -1.0*quat1[3] };

			double bunge1pas[3];
			ensmath.quat2euler( quat1pas, bunge1pas );

			//double dpas1 = fabs(bunge1[0]-bunge1pas[0]) + fabs(bunge1[1]-bunge1pas[1]) + fabs(bunge1[2]-bunge1pas[2]);

			//cout << px+1 << ";" << setprecision(32) << mis12 << ";" << setprecision(32) << mis21 << ";" << setprecision(32) << d1 << ";" << setprecision(32) << d2 << endl;
			cout << px+1 << ";" << setprecision(32) << bunge1[0] << ";" << setprecision(32) << bunge1[1] << ";" << setprecision(32) << bunge1[2] << ";" << setprecision(32) << bunge1pas[0] << ";" << setprecision(32) << bunge1pas[1] << ";" << setprecision(32) << bunge1pas[2] << endl;
		}
	}

	scoreangfile.close();
	snpangfile.close();
}


void ensembleHdl::ang_file_coloring( const char* inputang, unsigned int nskip )
{
	//MK::testing function to export ang file color coding with internal routines
	//reads three headered EBSD file of format x,y,gid,phi1,Phi,phi2, kam, optional RGB
	std::cout << "Reading SEM/EBSD data sequentially and colorize according to IMM scheme for test purposes..." << std::endl;

	//allocate memory to store the resulting image, ##MK::currently hardcoded dimensions as dummy function
	uint32_t imgx = 440;
	uint32_t imgy = 1197;

	unsigned char* img_rgba = NULL;
	img_rgba = new unsigned char [4*imgx*imgy];

	uint32_t c = 0;
	for ( uint32_t y = 0; y < imgy; y++ ) {
		for ( uint32_t x = 0; x < imgx; x++ ) {
			c = 4*(x + (y*imgx));
			img_rgba[c + REDCHAN ] = BLACK;
			img_rgba[c + GREENCHAN ] = BLACK;
			img_rgba[c + BLUECHAN ] = BLACK;
			img_rgba[c + ALPHACHAN ] = UCHAR_RANGE_MAX;
		}
	}

	string fname;
	fname = "SCORE.DebugANGFileIPF001.png";

	ifstream ebsdfile;
	string ebsdline;
	istringstream line;
	string datapiece;

	ebsdfile.open( inputang );

	c = 0; //4*(imgx*imgy - 1 ); //work backwards thus explicit vertical flip
	if ( ebsdfile.is_open() == true ) {
		//user format currently no header, format is phi1, PHI, phi2 from TSL
		//for ( unsigned int h = 0; h < 135; h++ ) getline( ebsdfile, ebsdline );
		//user format 6-lines header thaen bunge1,bunge2,Bunge3 in rad
		for ( unsigned int h = 0; h < nskip; h++ ) getline( ebsdfile, ebsdline );
cout << "I have now skipped " << nskip << " header lines!" << endl;

		while ( ebsdfile.good() == true ) { //read data
			getline( ebsdfile, ebsdline );
			if ( ebsdline.size() > 6 ) { //required when attempting to read x, y, gid, bunge1,2,3, kam
				istringstream line( ebsdline );

				double bunge_ori[3];
				getline(line, datapiece, '\t'); bunge_ori[0] = atof( datapiece.c_str() ); //MK::*.input already in rad
				getline(line, datapiece, '\t');	bunge_ori[1] = atof( datapiece.c_str() );
				getline(line, datapiece, '\t');	bunge_ori[2] = atof( datapiece.c_str() );

				double quaternion_mod[4];
				ensmath.euler2quat( bunge_ori, quaternion_mod );

				double ax[3] = {0.0, 0.0, 1.0}; //z ||ND for instance
				unsigned char rgb[3] = { UCHAR_RANGE_MIN, UCHAR_RANGE_MIN, UCHAR_RANGE_MIN };
				double locationinsst[2] = {FAIL_MYMATH_NUMERIC, FAIL_MYMATH_NUMERIC};

				ensmath.oriquat2sst_cubic_ipf_quat2pos( quaternion_mod, ax, locationinsst );
				ensmath.oriquat2sst_cubic_ipf_pos2rgb( locationinsst, rgb, IPF_COLOR_STRETCH_R, IPF_COLOR_STRETCH_G, IPF_COLOR_STRETCH_B );

				img_rgba[c+REDCHAN] = rgb[0];
				img_rgba[c+GREENCHAN] = rgb[1];
				img_rgba[c+BLUECHAN] = rgb[2];
				//c = c - 4;
				c = c + 4;
				//##MK::potentially reverse image build up along y direction...
			}
		}

		ebsdfile.close();
	}
	else { cout << "File " << inputang << " either not existent or cannot be opened!" << endl; }

	//encode the PNG file by calling the lodePNG library
	lodepng::encode( fname.c_str(), img_rgba , imgx, imgy );

	delete [] img_rgba; img_rgba = NULL;
}


void ensembleHdl::ipf_colormap( double ipf_stretch_r, double ipf_stretch_g, double ipf_stretch_b, unsigned int id )
{
	//generates an IPFZ colormap
	//allocate memory to store the resulting image
	uint32_t imgsz = 1000;
	uint32_t imgx = imgsz;
	uint32_t imgy = imgsz;
	double xysc = 1.0 / ((double) imgsz) * (pow(2.0, 0.5) - 1.0);
	double SQRR = SQR(pow(2.0, 0.5));

	unsigned char* img_rgba = NULL;
	img_rgba = new unsigned char [4*imgx*imgy];
	uint32_t c = 0;
	for ( uint32_t y = 0; y < imgy; y++ ) {
		uint32_t yoff = y*imgx;
		double yy = y;
		yy *= xysc;

		if ( yy > (0.5*pow(3.0, 0.5) - 0.5) ) //too high up
			continue;

		for ( uint32_t x = 0; x < imgx; x++ ) {
			double xx = x;
			xx *= xysc;

			//inside the stereographic standard triangle?
			c = 4*(x + yoff);
			img_rgba[c + REDCHAN ] = BLACK;
			img_rgba[c + GREENCHAN ] = BLACK;
			img_rgba[c + BLUECHAN ] = BLACK;
			img_rgba[c + ALPHACHAN ] = UCHAR_RANGE_MAX;

			if ( yy > xx ) 
				continue; //upper triangle is not inside, remains black

			//beyond projected great circle segment, boundary between (011) and (111) is segment of a circle radius sqrt(2) with center xc,yc = -1,0
			if ( (SQR(xx + 1.0) + SQR(yy)) > SQRR ) 
				continue;

			unsigned char rgb[3] = { UCHAR_RANGE_MIN, UCHAR_RANGE_MIN, UCHAR_RANGE_MIN };
			double locationinsst[2] = {xx, yy};

			ensmath.oriquat2sst_cubic_ipf_pos2rgb( locationinsst, rgb, ipf_stretch_r, ipf_stretch_g, ipf_stretch_b );

			img_rgba[c + REDCHAN ] = rgb[0];
			img_rgba[c + GREENCHAN ] = rgb[1];
			img_rgba[c + BLUECHAN ] = rgb[2];
		}
	}

	string fname;
	fname = "SCORE.IPFColorCoding.ID." + std::to_string(id) + ".png";

	//encode the PNG file by calling the lodePNG library
	lodepng::encode( fname.c_str(), img_rgba , imgx, imgy );

	delete [] img_rgba; img_rgba = NULL;
}


bool ensembleHdl::init_parameter( void )
{
	//read user input, master checks data integrity
	if ( myRank == MASTER ) 
		cout << "Initializing parameter/physics..." << endl;

	bool good = true;
	if ( good == true )	{
		good = init_parameter_runtime();
		cout << "Initialized runtime..." << "\n";
	}
	if ( good == true )	 {
		good = init_parameter_output();
		cout << "Initialized output..." << "\n";
	}
	if ( good == true )	{
		good = init_parameter_ensemble();
		cout << "Initialized ensemble..." << "\n";
	}
	if ( good == true )	 {
		good = init_parameter_nucleation();
		cout << "Initialized nucleation..." << "\n";
	}
	if ( good == true )	{
		good = init_parameter_material();
		cout << "Initialized material..." << "\n";
	}
	if ( good == true )	{
		good = init_parameter_mobilities();
		cout << "Initialized mobilites..." << "\n";
	}
	if ( good == true )	{
		good = init_parameter_recovery();
		cout << "Initialized recovery..." << "\n";
	}
	if ( good == true )	{
		good = init_parameter_drag();
		cout << "Initialized drag..." << "\n";
	}
	if ( good == true )	{
		good = init_parameter_idealtexture();
		cout << "Initialized ideal texture..." << "\n";
	}
	if ( good == true )	{
		good = init_parameter_processing();
		cout << "Initialized processing..." << "\n";
	}
	if ( good == true ) {
		if ( experimentInput == false ) {
			good = init_parameter_defgpool();
			cout << "Initialized defgpool..." << "\n";
		}
		else {
			good = init_parameter_defgpool_exp();
			cout << "Initialized defgpoolexp..." << "\n";
		}
	}
	if ( good == true ) {
		if ( experimentInput == false ) {
			good = init_parameter_rxgpool();
			cout << "Initialized rxgpool..." << "\n";
		}
		else {
			good = init_parameter_rxgpool_exp();
			cout << "Initialized rxgpoolexp..." << "\n";
		}
	}
	if ( good == true )	{
		good = init_parameter_additional();
		cout << "Initialized additional..." << "\n";
	}

	int mystatus = 0;
	if ( good == true ) mystatus = 1;
	int worldstatus = 0;
	MPI_Allreduce( &mystatus, &worldstatus, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
	if ( worldstatus == nRanks ) { //have all ranks been able to load the input?
		if ( myRank == MASTER ) { report_successful_setup(); }
		return true;
	}
	else {
		return false;
	}
}

void ensembleHdl::init_distributeWorkOnRanks( void )
{
	//partition ensemble of solitary units on the ranks in MPI_COMM_WORLD
	if ( myRank == MASTER ) {
		cout << "Partitioning CA ensemble on MPI processes..." << endl;
		if ( (nworldCAs % nRanks) != 0 ) { cout << "WARNING::Choose number of nuclei as a multiple number of nRanks to prevent excessive microstopping before barriers." << endl; }
	}

	//work partitioning at the moment simple a modulo round robin operation ensuring that ALL RANKS come to the same deterministic partitioning scheme
	//MK::modulo on uint % int is implementation specific and possibly non positive when the divisor int is negative or zero then even undefined, however nRanks >= 1
	uint32_t luckyRank;
	for (uint32_t ca = 0; ca < nworldCAs; ca++) {

		luckyRank = ca % nRanks; //number % x with x == 0 is undefined, but nRanks >= 1

		WhichCAonWhichRank.push_back(luckyRank);

		if ( myRank == luckyRank ) {
			myIDs.push_back ( ca );
		}
	}
	cout << "\t\t" << myRank << " got " << myIDs.size() << " CAs to process" << endl;
}

void ensembleHdl::init_distributeWorkOnThreads( void )
{
	//MK::currently each process works sequentially through a queue (myIDs) of SU domains
	//each in turn is executed with all available OpenMP threads
	//while this may have slightly more overhead than handing even more coarsely one SU to each thread
	//this choice has the benefit of enabling parallel execution of also very large CAs
	for ( uint32_t tid = 0; tid < ALL_THREADS_WORK_ONONECA; tid++) {
		vector<caHdlP> emptyBucketThreadLocalCAs;
		myCAs.push_back( emptyBucketThreadLocalCAs);
	}
	//nested vector-vector construct okay because handles only pointer to caHdl class objects, which internally spawn threads independently

	for ( uint32_t myca = 0; myca < myIDs.size(); myca++) {
		caHdlP aca = NULL;
		aca = new caHdl;
		QUICKASSERT( aca != NULL );

		ensembleHdl::myCAs[DEBUG_UTILIZE_ONLY_ONE_THREAD].push_back( aca );
		ensMemGuard = ensMemGuard + sizeof(caHdl);

		//MK::assignment after push_back because aca is not a trivial class object
		caHdlP tmpaccessthatca = this->myCAs[DEBUG_UTILIZE_ONLY_ONE_THREAD][myCAs[DEBUG_UTILIZE_ONLY_ONE_THREAD].size()-1];
		tmpaccessthatca->myensHdl = this;
		tmpaccessthatca->jobid = this->myIDs[myca]; //assure this is the global ID
		tmpaccessthatca->myensRank = this->myRank;
		tmpaccessthatca->nRanks = this->nRanks;

		//check if this results from this SU should be rendered
		tmpaccessthatca->renderingForThisCA = false;
		for ( uint32_t c = 0; c < ensembleHdl::UserDefLogPoint_WhichCAtoOutput.size(); c++ ) {
			if ( ensembleHdl::UserDefLogPoint_WhichCAtoOutput[c] == tmpaccessthatca->jobid ) {
				tmpaccessthatca->renderingForThisCA = true;
				break;
			}
		}

#ifdef REPORTSTYLE_DEVELOPER
			cout << tmpaccessthatca->myensRank << "\t\t" << myRank << "\t\t" << myca << "\t\t" << tmpaccessthatca->myensHdl << "\t\t" << tmpaccessthatca->jobid << "\t\t" << "\ttmpaccessthatca->myensRank;myRank;myca;myensHdlAddress;jobid" << endl;
#endif
	} //for all the solitary unit CA domains of the MPI process
}

void ensembleHdl::SIMULATE_myCAs( void )
{
	//MPI rank works through queue of solitary unit domains, each domain is solved with OMP_NUM_THREADS
	for ( uint32_t tid = 0; tid < ALL_THREADS_WORK_ONONECA; tid++) {
		for ( uint32_t mycathreaded = 0; mycathreaded < myCAs[tid].size(); mycathreaded++ ) {
			//analyze the automaton theca with private properties theensHdl  //ensembleHdlP theens = theca->myensHdl; have to be private temporaries to the thread
			struct profilingData prof;
			double timer[7] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
			timer[0] = MPI_Wtime();

			caHdlP theca = myCAs[tid][mycathreaded];

			theca->solve_INITIALIZATION();

			theca->solve_INITIALIZATION_THREADED_DOMAIN_DECOMPOSITION();
			timer[1] = MPI_Wtime();
			prof.MPIWTimeInitialization = timer[1] - timer[0];


			theca->solve_SYNTHESIZE_DEFORMATIONSTRUCTURE();
			timer[2] = MPI_Wtime();
			prof.MPIWTimeSpendDefMS = timer[2] - timer[1];


			if ( theca->myNucleationModel.gbnucleation != GBNUCLEATION_NO ) {
				theca->solve_DETECT_GRAINBOUNDARIES();
			}
			timer[3] = MPI_Wtime();
			prof.MPIWTimeSpendGBDetection = timer[3] - timer[2];


			//##MK::nucleation does not run OpenMP parallel
			if ( theca->myNucleationModel.csrnucleation != CSRNUCLEATION_NO || theca->myNucleationModel.gbnucleation != GBNUCLEATION_NO ) {
				theca->solve_NUCLEATIONMODELING();
			}
			timer[4] = MPI_Wtime();
			prof.MPIWTimeSpendNucleation = timer[4] - timer[3];


			theca->solve_REPLACE_CA_STRUCTURE();
			//MK::now everything is prepared for placing nuclei from myrxgpool let them grow

			theca->solve_RXGROWTH();

			timer[5] = MPI_Wtime();
			prof.MPIWTimeSpendGrowthSim = timer[5] - timer[4];

			theca->solve_FINAL_IO();

			timer[6] = MPI_Wtime();
			prof.MPIWTimeSpendFinalIO = timer[6] - timer[5];

			//commit profiling
			prof.JobID = theca->jobid; //##MK::consider change to UNSIGNED
			prof.Iterations = theca->step;
			prof.tend = theca->t;
			myCAProfiler.push_back ( prof );
			struct loginfo_ca_physics simphysics;
			theca->log_ca_physics( &simphysics );
			myCAPhysics.push_back ( simphysics );

		} //execute next SU with all threads
	}
}

void ensembleHdl::PERFORM_SOLITARYUNIT_MODELING_ANALYSIS( void ) 
{
	if ( myRank == MASTER ) 
		cout << "Performing solitary unit analysis..." << endl;

	postprocess_write_mycas_mpislow(); //sequential
	//postprocess_write_mycas_mpifast(); //gatherv

	if ( postProcessing == true ) {
		postprocess_init();
		postprocess_rediscr_kinetics();
		postprocess_rediscr_macrotexture();
		postprocess_rediscr_finalgrainsizedistribution();
	}
}


void ensembleHdl::destroy_myCAs( void )
{
	//necessary because the ensembleHdl member vector<caHdlP> can only deconstruct the vector but not the heap memory the individual element pointers are pointing to
	for ( uint32_t tid = 0; tid < myCAs.size(); tid++) {
		//cout << myRank << "\t\t" << tid << "\t\tmyRank;tid deallocating..." << endl;
		for ( uint32_t tidlocalca = 0; tidlocalca < myCAs[tid].size(); tidlocalca++) {
			//for ( uint32_t nreg = 0; myCAs[tid][tidlocalca]->regions.size(); nreg++ ) {
			//	myCAs[tid][tidlocalca]->regions[nreg]->~caregionMemHdl();
			//}
			myCAs[tid][tidlocalca]->cleanMyGrainEvolution();
			myCAs[tid][tidlocalca]->cleanBookkeeping();

			/*##MK::memory flaw
			delete myCAs[tid][tidlocalca]; //this is expected, but calls the constructor of the caHdl class that automatically destroys memory used for representing the automaton and storage of intermediate results
			myCAs[tid][tidlocalca] = NULL;
			ensMemGuard = ensMemGuard - sizeof(caHdl);*/
			//cout << myRank << "\t\t" << tid << "\t\t" << tidlocalca << "\t\tmyRank;tid;tidlocalca successful." << endl;
		}
	}
}


caregionMemHdl::caregionMemHdl()
{
	mycaHdl = NULL;
	mythreadid = MASTER;

	myGeom.nreg_rdmin = 0;
	myGeom.nreg_rdmax = 0;
	myGeom.nreg_tdmin = 0;
	myGeom.nreg_tdmax = 0;
	myGeom.nreg_ndmin = 0;
	myGeom.nreg_ndmax = 0;
	myGeom.nreg_rd = 1;
	myGeom.nreg_td = 1;
	myGeom.nreg_nd = 1;
	myGeom.nregarea_rdtd = 1 * 1;
	myGeom.nregvol_rdtdnd = 1 * 1 * 1;

	myGeom.nedge_global_rd = 1;
	myGeom.nedge_global_td = 1;
	myGeom.nedge_global_nd = 1;

	myGeom.cellsize = DEFAULT_CELLSIZE;

	thePartitioning.nreg_rdx = 1;
	thePartitioning.nreg_tdy = 1;
	thePartitioning.nreg_ndz = 1;

	regMemGuard = 0.0;
	dtCalcGrowthInside = 0.0;
	dtCalcGrowthBorder = 0.0;
	dtUpdateInside = 0.0;
	dtUpdateBorder = 0.0;

	myhalos = NULL;
	mynborhalos = NULL;
	mynborregions = NULL;

	mycellgrid = NULL;

	//3D CA VOXELIZING
	mySeedFrontInside = NULL;
	mySeedFrontBorder = NULL;
	myFullSeedInside = NULL;
	myFullSeedBorder = NULL;
	myRecyclingSeedInside = NULL;
	myRecyclingSeedBorder = NULL;

	ntotalSeedFrontInside = 0;
	ntotalSeedFrontBorder = 0;
	nextSlotNeverActiveSeedInside = 0;
	nextSlotNeverActiveSeedBorder = 0;
	ncurrentActiveSeedInside = 0;
	ncurrentActiveSeedBorder = 0;

	ntotalFullSeedInside = 0;
	ntotalFullSeedBorder = 0;
	nextSlotToFullSeedInside = 0;
	nextSlotToFullSeedBorder = 0;

	ntotalRecyclingSeedInside = 0;
	ntotalRecyclingSeedBorder = 0;
	nextToRecycleSeedInside = 0;
	nextToRecycleSeedBorder = 0;
	firstNotRecycleSeedInside = 0;
	firstNotRecycleSeedBorder = 0;


	//3D CA RX SIM
	myRXFrontInside = NULL;
	myRXFrontBorder = NULL;
	myFullRXListInside = NULL;
	myFullRXListBorder = NULL;
	myRecyclingListInside = NULL;
	myRecyclingListBorder = NULL;

	ntotalRXFrontInside = 0;
	ntotalRXFrontBorder = 0;
	nextSlotNeverActiveRXFrontInside = 0;
	nextSlotNeverActiveRXFrontBorder = 0;
	nCurrentlyActiveInside = 0;
	nCurrentlyActiveBorder = 0;

	ntotalFullRXListInside = 0;
	ntotalFullRXListBorder = 0;
	nextSlotToFullRXInside = 0;
	nextSlotToFullRXBorder = 0;

	ntotalRecyclingListInside = 0;
	ntotalRecyclingListBorder = 0;
	nextSlotThatBecomesRecycledInside = 0;
	nextSlotThatBecomesRecycledBorder = 0;
	firstNotRecycledYetInside = 0;
	firstNotRecycledYetBorder = 0;

	reg_nmydefgpool = 0;

	reg_nCurrActiveInside = 0;
	reg_nCurrActiveBorder = 0;
	reg_SvInside = 0;
	reg_SvBorder = 0;
	reg_dXstepInside = 0.0;
	reg_dXstepBorder = 0.0;

	reg_myMobilityWeightMax = DEFAULT_PMAX;
	reg_nmynuclei = 0;

	reg_defrag_cnt = 0;

	reg_myrxgp_counts = NULL;
	reg_nmyrxgp_counts = 0;
	reg_mydfgp_counts = NULL;
	reg_nmydfgp_counts = 0;
}


caregionMemHdl::~caregionMemHdl()
{
/*
	delete [] myRXFrontInside;
	myRXFrontInside = NULL;
	regMemGuard = regMemGuard - (ntotalRXFrontInside * sizeof(cellP));
	ntotalRXFrontInside = 0;
	nextSlotNeverActiveRXFrontInside = 0;
	nCurrentlyActiveInside = 0;

	delete [] myRXFrontBorder;
	myRXFrontBorder = NULL;
	regMemGuard = regMemGuard - (ntotalRXFrontBorder * sizeof(cellP));
	ntotalRXFrontBorder = 0;
	nextSlotNeverActiveRXFrontBorder = 0;
	nCurrentlyActiveBorder = 0;

	delete [] myFullRXListInside;
	myFullRXListInside = NULL;
	regMemGuard = regMemGuard - (ntotalFullRXListInside * sizeof(uint32_t));
	ntotalFullRXListInside = 0;
	nextSlotToFullRXInside = 0;

	delete [] myFullRXListBorder;
	myFullRXListBorder = NULL;
	regMemGuard = regMemGuard - (ntotalFullRXListBorder * sizeof(uint32_t));
	ntotalFullRXListBorder = 0;
	nextSlotToFullRXBorder = 0;

	delete [] myRecyclingListInside;
	myRecyclingListInside = NULL;
	regMemGuard = regMemGuard - (ntotalRecyclingListInside  * sizeof(uint32_t));
	ntotalRecyclingListInside = 0;
	nextSlotThatBecomesRecycledInside = 0;
	firstNotRecycledYetInside = 0;

	delete [] myRecyclingListBorder;
	myRecyclingListBorder = NULL;
	regMemGuard = regMemGuard - (ntotalRecyclingListBorder  * sizeof(uint32_t));
	ntotalRecyclingListBorder = 0;
	nextSlotThatBecomesRecycledBorder = 0;
	firstNotRecycledYetBorder = 0;
*/

	//reg_myrxgp_counts, reg_mydfgp_counts, profiling_growthmachine are vectors...

	//##MK::allocation at the moment problematic!
	//delete [] mynborregions;
	//mynborregions = NULL;

	//own deallocator for mycellgrid

	/*delete [] reg_myrxgp_counts;
	reg_myrxgp_counts = NULL;
	reg_nmyrxgp_counts = 0;

	delete [] reg_mydfgp_counts;
	reg_mydfgp_counts = NULL;
	reg_nmydfgp_counts = 0;*/
}


void caregionMemHdl::ompshar_init_mycellgrid ( void )
{
	//CALLED FROM WITHIN PARALLEL REGION

	//allocate memory threadlocally for the voxelized deformed structure, because
	//local variables in encapsulated function are private by default
	mycellgrid = NULL;
	mycellgrid = new uint32_t[myGeom.nregvol_rdtdnd];
	QUICKASSERT ( mycellgrid != NULL );
	regMemGuard = regMemGuard + ( myGeom.nregvol_rdtdnd * sizeof(uint32_t) );

	//first touch in parallel to assure data locality if possible (and threads are pinned!)...
	for ( uint32_t cxyz = 0; cxyz < myGeom.nregvol_rdtdnd; cxyz++ ) {
		mycellgrid[cxyz] = NO_GRAIN_ASSIGNED;
	}
}


uint32_t caregionMemHdl::get_nbortid( int regid, short dx, short dy, short dz, int npx, int npy, int npz )
{
	//CALLED FROM WITHIN PARALLEL REGION
	//determines ID field of neighbors in an implicit, topology of regions identified from 0, 1, 2, omp_get_num_threads() - 1 that are
	//arranged on a party * partz region grid with party <= partz
	//old way::
	//addressing is first in z at fixed y then in y+1 along z, y+2 along z and so forth...
	//so party = 2 and partz = 4 generates: 0, 1,2,3 with y = 0 and 4,5,6,7 with y = 1

	//addressing is first in y at fixed z then in z+1 along +y, z+2 along +y and so forth...
	//so party = 3 and partz = 4 generates: 0,1,2 with z=0, 3,4,5 at z=1, 6,7,8 at z=2, 9,10,11 at z=3

	//get coordinates of regid
	int ix = 0; //no partitioning along x
	int iy = regid % npy; //neither regid nor npy negative
	int iz = regid / npy;

	//get new coordinate
	ix += dx;
	iy += dy;
	iz += dz;

	//link memory regions periodic boundary conditions
	if ( ix < 0 )		ix += npx;
	if ( ix >= npx )	ix -= npx;
	if ( iy < 0 )		iy += npy;
	if ( iy >= npy )	iy -= npy;
	if ( iz < 0 )		iz += npz;
	if ( iz >= npz )	iz -= npz;
	//MK::ix, iy, and iz are assured >= 0, hence implicit cast to uint32_t is safe

	return ( (iz * npy) + iy );
}


void caregionMemHdl::ompshar_init_neighbortopology( uint32_t party, uint32_t partz )
{
	//CALLED FROM WITHIN PARALLEL REGION
	//determines ID field of neighbors in an implicitly linearly stored but explicitly defined 3D topology 
	//of camemoryregions identified from 0, 1, 2, omp_get_num_threads() - 1
	//these regions are arranged on a party * partz region grid with party <= partz
	//addressing is first in y at fixed z then in z+1 along +y, z+2 along +y and so forth...
	//so party = 2 and partz = 4 generates: 0,1 with z=0, 2,3 at z=1, 4,5 at z=2, 6,7 at z=3
	//in this manner regions are periodic so in that example

	//a Moore stencil kernel can issue infections into at most NUMBER_OF_NEIGHBORS adjacent memory regions
	mynborregions = NULL;
	mynborregions = new uint32_t[(MYSELF+NUMBER_OF_NEIGHBORS)*(IDS_AND_LIMITS)];
	QUICKASSERT ( mynborregions != NULL );

	//the order of this array is arranged in such a manner that the target region can be found fastest on average
	//in fact, we utilize that it is much more likely that an infection into a neighboring region directs into the
	//face than the edge than the corners...

	//myself is the first entry
	mynborregions[(0*IDS_AND_LIMITS)+THE_IDS] = this->mythreadid; //thread ids are positive or zero

	//six faces
	mynborregions[(1*IDS_AND_LIMITS)+THE_IDS] = this->get_nbortid( this->mythreadid, -1, 0, 0, 1, party, partz );
	mynborregions[(2*IDS_AND_LIMITS)+THE_IDS] = this->get_nbortid( this->mythreadid, +1, 0, 0, 1, party, partz );
	mynborregions[(3*IDS_AND_LIMITS)+THE_IDS] = this->get_nbortid( this->mythreadid, 0, -1, 0, 1, party, partz );
	mynborregions[(4*IDS_AND_LIMITS)+THE_IDS] = this->get_nbortid( this->mythreadid, 0, +1, 0, 1, party, partz );
	mynborregions[(5*IDS_AND_LIMITS)+THE_IDS] = this->get_nbortid( this->mythreadid, 0, 0, -1, 1, party, partz );
	mynborregions[(6*IDS_AND_LIMITS)+THE_IDS] = this->get_nbortid( this->mythreadid, 0, 0, +1, 1, party, partz );

	//twelve edges
	mynborregions[(7*IDS_AND_LIMITS)+THE_IDS] = this->get_nbortid( this->mythreadid, -1, 0, -1, 1, party, partz );
	mynborregions[(8*IDS_AND_LIMITS)+THE_IDS] = this->get_nbortid( this->mythreadid, +1, 0, -1, 1, party, partz );
	mynborregions[(9*IDS_AND_LIMITS)+THE_IDS] = this->get_nbortid( this->mythreadid, -1, 0, +1, 1, party, partz );
	mynborregions[(10*IDS_AND_LIMITS)+THE_IDS] = this->get_nbortid( this->mythreadid, +1, 0, +1, 1, party, partz );

	mynborregions[(11*IDS_AND_LIMITS)+THE_IDS] = this->get_nbortid( this->mythreadid, -1, -1, 0, 1, party, partz );
	mynborregions[(12*IDS_AND_LIMITS)+THE_IDS] = this->get_nbortid( this->mythreadid, -1, +1, 0, 1, party, partz );
	mynborregions[(13*IDS_AND_LIMITS)+THE_IDS] = this->get_nbortid( this->mythreadid, +1, -1, 0, 1, party, partz );
	mynborregions[(14*IDS_AND_LIMITS)+THE_IDS] = this->get_nbortid( this->mythreadid, +1, +1, 0, 1, party, partz );

	mynborregions[(15*IDS_AND_LIMITS)+THE_IDS] = this->get_nbortid( this->mythreadid, 0, -1, -1, 1, party, partz );
	mynborregions[(16*IDS_AND_LIMITS)+THE_IDS] = this->get_nbortid( this->mythreadid, 0, -1, +1, 1, party, partz );
	mynborregions[(17*IDS_AND_LIMITS)+THE_IDS] = this->get_nbortid( this->mythreadid, 0, +1, -1, 1, party, partz );
	mynborregions[(18*IDS_AND_LIMITS)+THE_IDS] = this->get_nbortid( this->mythreadid, 0, +1, +1, 1, party, partz );

	//eight corners, eight corners you have ...
	mynborregions[(19*IDS_AND_LIMITS)+THE_IDS] = this->get_nbortid( this->mythreadid, -1, -1, -1, 1, party, partz );
	mynborregions[(20*IDS_AND_LIMITS)+THE_IDS] = this->get_nbortid( this->mythreadid, +1, -1, -1, 1, party, partz );
	mynborregions[(21*IDS_AND_LIMITS)+THE_IDS] = this->get_nbortid( this->mythreadid, -1, +1, -1, 1, party, partz );
	mynborregions[(22*IDS_AND_LIMITS)+THE_IDS] = this->get_nbortid( this->mythreadid, +1, +1, -1, 1, party, partz );
	mynborregions[(23*IDS_AND_LIMITS)+THE_IDS] = this->get_nbortid( this->mythreadid, -1, -1, +1, 1, party, partz );
	mynborregions[(24*IDS_AND_LIMITS)+THE_IDS] = this->get_nbortid( this->mythreadid, +1, -1, +1, 1, party, partz );
	mynborregions[(25*IDS_AND_LIMITS)+THE_IDS] = this->get_nbortid( this->mythreadid, -1, +1, +1, 1, party, partz );
	mynborregions[(26*IDS_AND_LIMITS)+THE_IDS] = this->get_nbortid( this->mythreadid, +1, +1, +1, 1, party, partz );
	//MK::the farther apart in the array positioned are the least likely infected region
}


void caregionMemHdl::ompshar_initialize_halo( uint32_t whichhalo, uint32_t extx, uint32_t exty, uint32_t extz, uint32_t originx, uint32_t originy,uint32_t originz )
{
	//CALLED FROM WITHIN PARALLEL REGION
	myhalos[whichhalo].extendx = extx;
	myhalos[whichhalo].extendy = exty;
	myhalos[whichhalo].extendz = extz;
	myhalos[whichhalo].extendxy = extx * exty;
	myhalos[whichhalo].extendxyz = extx * exty * extz;

	myhalos[whichhalo].originx = originx;
	myhalos[whichhalo].originy = originy;
	myhalos[whichhalo].originz = originz;

	myhalos[whichhalo].inregion = mynborregions[(whichhalo*IDS_AND_LIMITS)+THE_IDS];

	//do not allocate memory if extendxyz is zero, then the halo is considered a dummy
	//it does strictly require no storage but keeps the setup of the halos for each region consistent
	uint32_t* refs = NULL;
	halocellP bucket = NULL;

	uint32_t xyz = extx * exty * extz;
	if ( xyz > 0 ) { //IF NOT A DUMMY HALO!
		refs = new uint32_t[xyz];
		QUICKASSERT( refs != NULL);
		//MK::initialize the halo with a specific marker to enable the identification of potential logic flaws, this is however not necessary!
		for (uint32_t c = 0; c < xyz; c++ ) { refs[c] = INVALID_ADDRESS; }

		bucket = new halocell[xyz];
		QUICKASSERT( bucket != NULL);

		//initialize fresh halo cells
		for (uint32_t c = 0; c < xyz; c++ ) {
			bucket[c].status = HALOCELL_NEVER_VISITED;
			bucket[c].infector = 26;
			bucket[c].ix = -1;
			bucket[c].iy = -1;
			bucket[c].iz = -1;
			bucket[c].rxFrac = 0.0;
			bucket[c].myrxgid = NO_GRAIN_ASSIGNED;
		}
	}

	myhalos[whichhalo].nextFreeHaloRefs = 0; //NO_HALO_CELL_INFECTED;
	myhalos[whichhalo].theHaloRefs = refs;
	myhalos[whichhalo].theHaloCells = bucket;

//##MK::#pragma omp critical{cout << "HALO-INITIALIZED;byWhom;whichHalo;inWhichRegion;extx/y/z/xy/xyz/;originx/y/z = " << mythreadid << ";" << whichhalo << ";" << "-->" << myhalos[whichhalo].inregion << "<--" << myhalos[whichhalo].extendx << "__" << myhalos[whichhalo].extendy << "__" << myhalos[whichhalo].extendz << "--" << myhalos[whichhalo].extendxy << "--" << myhalos[whichhalo].extendxyz << "__" << myhalos[whichhalo].originx << "--" << myhalos[whichhalo].originy << "--" << myhalos[whichhalo].originz <<  endl;}
}


void caregionMemHdl::determine_origins( uint32_t regid, short dx, short dy, short dz, int npx, int npy, int npz, uint32_t * origin )
{
	//THIS FUNCTION IDENTIFIES THE ORIGIN OF HALOREGIONS THAT ARE NOT DUMMIES I.E: whose extendxyz > 0
	//IT MUST NOT BE CALLED IF THIS IS NOT ASSURED!

	//CALLED FROM WITHIN PARALLEL REGION, determines ID field of neighbors in an implicit
	//topology of regions identified from 0, 1, 2, omp_get_num_threads() - 1 that are
	//arranged on a party * partz region grid with party <= partz
	//addressing is first in z at fixed y then in y+1 along z, y+2 along z and so forth...
	//so party = 2 and partz = 4 generates: 0, 1,2,3 with y = 0 and 4,5,6,7 with y = 1

	//get coordinates of regid
	int ix = 0;
	int iy = regid % npy;
	//MK::modulo uint % int  could become a) int --> uint so uint % uint --> int (safe! because regid and npy << INT32_MAX)
	//b) uint --> int (safe!, because regid << INT32_MAX) and then int % int --> int (safe)
	int iz = regid / npy;

	//get new coordinate, //MK::int --> short safe here because all coordinates positive and << SHORT_MAX
	ix += dx;
	iy += dy;
	iz += dz;

	//get automaton geometry
	int cax = myGeom.nedge_global_rd; //MK::implicit uint32_t 2 int not a problem because nedge_global_... > 0 and << INT32_RANGE
	int cay = myGeom.nedge_global_td;
	int caz = myGeom.nedge_global_nd;

	//determine origin
	int oox = UNKNOWN_ORIGIN;
	int ooy = UNKNOWN_ORIGIN;
	int ooz = UNKNOWN_ORIGIN;

	//x first, check periodic boundaries
	if ( dx == -1 ) {
		oox = myGeom.nreg_rdmin - 1; //still located in this region minus one takes me in next region
		if ( oox < 0 )
			oox = oox + cax;
	}
	if ( dx == 0 ) {
		oox = myGeom.nreg_rdmin; //cannot direct out of automaton
	}
	if ( dx == +1 ) {
		oox = myGeom.nreg_rdmax + 1;
		if ( oox >= cax )
			oox = oox - cax;
	}

	//y next
	if ( dy == -1 ) {
		ooy = myGeom.nreg_tdmin - 1;
		if ( ooy < 0 )
			ooy = ooy + cay;
	}
	if ( dy == 0 ) {
		ooy = myGeom.nreg_tdmin;
	}
	if ( dy == +1 ) {
		ooy = myGeom.nreg_tdmax + 1;
		if ( ooy >= cay )
			ooy = ooy - cay;
	}

	//z last
	if ( dz == -1 ) {
		ooz = myGeom.nreg_ndmin - 1;
		if ( ooz < 0 )
			ooz = ooz + caz;
	}
	if ( dz == 0 ) {
		ooz = myGeom.nreg_ndmin;
	}
	if ( dz == +1 ) {
		ooz = myGeom.nreg_ndmax + 1;
		if ( ooz >= caz )
			ooz = ooz - caz;
	}

	//##MK::DEBUG bounds check
	QUICKASSERT ( oox >= 0 && oox < UINT32T_MAX );
	QUICKASSERT ( ooy >= 0 && ooy < UINT32T_MAX );
	QUICKASSERT ( ooz >= 0 && ooz < UINT32T_MAX );

	origin[0] = oox;
	origin[1] = ooy;
	origin[2] = ooz;
}


void caregionMemHdl::ompshar_init_haloregions( uint32_t party, uint32_t partz )
{
	//CALLED FROM WITHIN PARALLEL REGION
	//serves to initialize the bookkeeping and management of the haloregions and to allocate memory for the storage of halocells
	myhalos = NULL;   //infection with Moore kernel on boundary of domain can direct at most into NUMBER_OF_NEIGHBORS disjoint regions
	myhalos = new halo[(MYSELF+NUMBER_OF_NEIGHBORS)];
	QUICKASSERT ( myhalos != NULL );

	//MAKE SURE FUNCTION IS CALLED WHEN mynborregions knows already the neighbors, get my own dimensions
	uint32_t ex = myGeom.nreg_rd;
	uint32_t ey = myGeom.nreg_td;
	uint32_t ez = myGeom.nreg_nd;

	//halo regions with extend 0 are dummies!, ##MK::one could safe such tiny memory consumption of a few KB but...
	//only dummy halo interface not required because periodic boundary conditions map back infections in border
	
	uint32_t o[3] = {INVALID_COORDINATE, INVALID_COORDINATE, INVALID_COORDINATE};
	//myself is the first entry //extend, origin
	ompshar_initialize_halo(0,      0, 0, 0,        0, 0, 0);

	//MK::the silent conversion of uin32_t to int could be dangerous but party and partz are >= 0 and << INT32_RANGE !
	//six faces
	uint32_t mythrid = this->mythreadid; //MK::comparisons of the kind int comp uint will always have the int argument casted to unsigned which can result in an undefined uint if int is negative however threadids are always >=0
	if ( mynborregions[(1*IDS_AND_LIMITS)+THE_IDS] == mythrid ) { //such if clause checks whether periodic boundaries cause that the halo is again in the same region
		ompshar_initialize_halo(1,  0,0,0,          0,0,0); //if so define a dummy
	} else {
		o[0] = INVALID_COORDINATE; o[1] = INVALID_COORDINATE; o[2] = INVALID_COORDINATE;
		determine_origins( mythrid, -1, 0, 0,		1, party, partz, o ); //if not define a real halo, i.e. with extend > 0
		ompshar_initialize_halo(1,  1,ey,ez,        o[0],o[1],o[2]);
	}

	if ( mynborregions[(2*IDS_AND_LIMITS)+THE_IDS] == mythrid ) {
		ompshar_initialize_halo(2,  0,0,0,          0,0,0);
	} else {
		o[0] = INVALID_COORDINATE; o[1] = INVALID_COORDINATE; o[2] = INVALID_COORDINATE;
		determine_origins( mythrid, +1, 0, 0,		1, party, partz, o );
		ompshar_initialize_halo(2,  1,ey,ez,        o[0],o[1],o[2]);
	}

	if ( mynborregions[(3*IDS_AND_LIMITS)+THE_IDS] == mythrid ) {
		ompshar_initialize_halo(3,  0,0,0,          0,0,0);
	} else {
		o[0] = INVALID_COORDINATE; o[1] = INVALID_COORDINATE; o[2] = INVALID_COORDINATE;
		determine_origins( mythrid, 0, -1, 0,		1, party, partz, o);
		ompshar_initialize_halo(3,  ex,1,ez,        o[0],o[1],o[2]);
	}

	if ( mynborregions[(4*IDS_AND_LIMITS)+THE_IDS] == mythrid ) {
		ompshar_initialize_halo(4,  0,0,0,          0,0,0);
	} else {
		o[0] = INVALID_COORDINATE; o[1] = INVALID_COORDINATE; o[2] = INVALID_COORDINATE;
		determine_origins( mythrid, 0, +1, 0,		1, party, partz, o);
		ompshar_initialize_halo(4,  ex,1,ez,        o[0],o[1],o[2]);
	}

	if ( mynborregions[(5*IDS_AND_LIMITS)+THE_IDS] == mythrid ) {
		ompshar_initialize_halo(5,  0,0,0,          0,0,0);
	} else {
		o[0] = INVALID_COORDINATE; o[1] = INVALID_COORDINATE; o[2] = INVALID_COORDINATE;
		determine_origins( mythrid, 0, 0, -1,		1, party, partz, o);
		ompshar_initialize_halo(5,  ex,ey,1,        o[0],o[1],o[2]);
	}

	if ( mynborregions[(6*IDS_AND_LIMITS)+THE_IDS] == mythrid ) {
		ompshar_initialize_halo(6,  0,0,0,          0,0,0);
	} else {
		o[0] = INVALID_COORDINATE; o[1] = INVALID_COORDINATE; o[2] = INVALID_COORDINATE;
		determine_origins( mythrid, 0, 0, +1,		1, party, partz, o);
		ompshar_initialize_halo(6,  ex,ey,1,        o[0],o[1],o[2]);
	}

	//twelve edges
	if ( mynborregions[(7*IDS_AND_LIMITS)+THE_IDS] == mythrid ) {
		ompshar_initialize_halo(7,  0,0,0,          0,0,0);
	} else {
		o[0] = INVALID_COORDINATE; o[1] = INVALID_COORDINATE; o[2] = INVALID_COORDINATE;
		determine_origins( mythrid, -1, 0, -1,			1, party, partz, o);
		ompshar_initialize_halo(7,  1,ey,1,         o[0],o[1],o[2]);
	}

	if ( mynborregions[(8*IDS_AND_LIMITS)+THE_IDS] == mythrid ) {
		ompshar_initialize_halo(8,  0,0,0,          0,0,0);
	} else {
		o[0] = INVALID_COORDINATE; o[1] = INVALID_COORDINATE; o[2] = INVALID_COORDINATE;
		determine_origins( mythrid, +1, 0, -1,		1, party, partz, o);
		ompshar_initialize_halo(8,  1,ey,1,         o[0],o[1],o[2]);
	}

	if ( mynborregions[(9*IDS_AND_LIMITS)+THE_IDS] == mythrid ) {
		ompshar_initialize_halo(9,  0,0,0,          0,0,0);
	} else {
		o[0] = INVALID_COORDINATE; o[1] = INVALID_COORDINATE; o[2] = INVALID_COORDINATE;
		determine_origins( mythrid, -1, 0, +1,		1, party, partz, o);
		ompshar_initialize_halo(9,  1,ey,1,         o[0],o[1],o[2]);
	}

	if ( mynborregions[(10*IDS_AND_LIMITS)+THE_IDS] == mythrid ) {
		ompshar_initialize_halo(10,  0,0,0,          0,0,0);
	} else {
		o[0] = INVALID_COORDINATE; o[1] = INVALID_COORDINATE; o[2] = INVALID_COORDINATE;
		determine_origins( mythrid, +1, 0, +1,		1, party, partz, o);
		ompshar_initialize_halo(10, 1,ey,1,         o[0],o[1],o[2]);
	}

	if ( mynborregions[(11*IDS_AND_LIMITS)+THE_IDS] == mythrid ) {
		ompshar_initialize_halo(11,  0,0,0,          0,0,0);
	} else {
		o[0] = INVALID_COORDINATE; o[1] = INVALID_COORDINATE; o[2] = INVALID_COORDINATE;
		determine_origins( mythrid, -1, -1, 0,		1, party, partz, o);
		ompshar_initialize_halo(11, 1,1,ez,         o[0],o[1],o[2]);
	}

	if ( mynborregions[(12*IDS_AND_LIMITS)+THE_IDS] == mythrid ) {
		ompshar_initialize_halo(12,  0,0,0,          0,0,0);
	} else {
		o[0] = INVALID_COORDINATE; o[1] = INVALID_COORDINATE; o[2] = INVALID_COORDINATE;
		determine_origins( mythrid, -1, +1, 0,		1, party, partz, o);
		ompshar_initialize_halo(12, 1,1,ez,         o[0],o[1],o[2]);
	}

	if ( mynborregions[(13*IDS_AND_LIMITS)+THE_IDS] == mythrid ) {
		ompshar_initialize_halo(13,  0,0,0,          0,0,0);
	} else {
		o[0] = INVALID_COORDINATE; o[1] = INVALID_COORDINATE; o[2] = INVALID_COORDINATE;
		determine_origins( mythrid, +1, -1, 0,		1, party, partz, o);
		ompshar_initialize_halo(13, 1,1,ez,         o[0],o[1],o[2]);
	}

	if ( mynborregions[(14*IDS_AND_LIMITS)+THE_IDS] == mythrid ) {
		ompshar_initialize_halo(14,  0,0,0,          0,0,0);
	} else {
		o[0] = INVALID_COORDINATE; o[1] = INVALID_COORDINATE; o[2] = INVALID_COORDINATE;
		determine_origins( mythrid, +1, +1, 0,		1, party, partz, o);
		ompshar_initialize_halo(14, 1,1,ez,         o[0],o[1],o[2]);
	}

	if ( mynborregions[(15*IDS_AND_LIMITS)+THE_IDS] == mythrid ) {
		ompshar_initialize_halo(15,  0,0,0,          0,0,0);
	} else {
		o[0] = INVALID_COORDINATE; o[1] = INVALID_COORDINATE; o[2] = INVALID_COORDINATE;
		determine_origins( mythrid, 0, -1, -1,		1, party, partz, o);
		ompshar_initialize_halo(15, ex,1,1,         o[0],o[1],o[2]);
	}

	if ( mynborregions[(16*IDS_AND_LIMITS)+THE_IDS] == mythrid ) {
		ompshar_initialize_halo(16,  0,0,0,          0,0,0);
	} else {
		o[0] = INVALID_COORDINATE; o[1] = INVALID_COORDINATE; o[2] = INVALID_COORDINATE;
		determine_origins( mythrid, 0, -1, +1,		1, party, partz, o);
		ompshar_initialize_halo(16, ex,1,1,         o[0],o[1],o[2]);
	}

	if ( mynborregions[(17*IDS_AND_LIMITS)+THE_IDS] == mythrid ) {
		ompshar_initialize_halo(17,  0,0,0,          0,0,0);
	} else {
		o[0] = INVALID_COORDINATE; o[1] = INVALID_COORDINATE; o[2] = INVALID_COORDINATE;
		determine_origins( mythrid, 0, +1, -1,		1, party, partz, o);
		ompshar_initialize_halo(17, ex,1,1,         o[0],o[1],o[2]);
	}

	if ( mynborregions[(18*IDS_AND_LIMITS)+THE_IDS] == mythrid ) {
		ompshar_initialize_halo(18,  0,0,0,          0,0,0);
	} else {
		o[0] = INVALID_COORDINATE; o[1] = INVALID_COORDINATE; o[2] = INVALID_COORDINATE;
		determine_origins( mythrid, 0, +1, +1,		1, party, partz, o);
		ompshar_initialize_halo(18, ex,1,1,         o[0],o[1],o[2]);
	}

	//eight corners, eight corners you have ...
		if ( mynborregions[(19*IDS_AND_LIMITS)+THE_IDS] == mythrid ) {
			ompshar_initialize_halo(19,  0,0,0,          0,0,0);
		} else {
			o[0] = INVALID_COORDINATE; o[1] = INVALID_COORDINATE; o[2] = INVALID_COORDINATE;
			determine_origins( mythrid, -1, -1, -1,	1, party, partz, o);
			ompshar_initialize_halo(19, 1,1,1,          o[0],o[1],o[2]);
		}

		if ( mynborregions[(20*IDS_AND_LIMITS)+THE_IDS] == mythrid ) {
			ompshar_initialize_halo(20,  0,0,0,          0,0,0);
		} else {
			o[0] = INVALID_COORDINATE; o[1] = INVALID_COORDINATE; o[2] = INVALID_COORDINATE;
			determine_origins( mythrid, +1, -1, -1,	1, party, partz, o);
			ompshar_initialize_halo(20, 1,1,1,          o[0],o[1],o[2]);
		}

		if ( mynborregions[(21*IDS_AND_LIMITS)+THE_IDS] == mythrid ) {
			ompshar_initialize_halo(21,  0,0,0,          0,0,0);
		} else {
			o[0] = INVALID_COORDINATE; o[1] = INVALID_COORDINATE; o[2] = INVALID_COORDINATE;
			determine_origins( mythrid, -1, +1, -1,	1, party, partz, o);
			ompshar_initialize_halo(21, 1,1,1,          o[0],o[1],o[2]);
		}

		if ( mynborregions[(22*IDS_AND_LIMITS)+THE_IDS] == mythrid ) {
			ompshar_initialize_halo(22,  0,0,0,          0,0,0);
		} else {
			o[0] = INVALID_COORDINATE; o[1] = INVALID_COORDINATE; o[2] = INVALID_COORDINATE;
			determine_origins( mythrid, +1, +1, -1,	1, party, partz, o);
			ompshar_initialize_halo(22, 1,1,1,          o[0],o[1],o[2]);
		}

		if ( mynborregions[(23*IDS_AND_LIMITS)+THE_IDS] == mythrid ) {
			ompshar_initialize_halo(23,  0,0,0,          0,0,0);
		} else {
			o[0] = INVALID_COORDINATE; o[1] = INVALID_COORDINATE; o[2] = INVALID_COORDINATE;
			determine_origins( mythrid, -1, -1, +1,	1, party, partz, o);
			ompshar_initialize_halo(23, 1,1,1,          o[0],o[1],o[2]);
		}

		if ( mynborregions[(24*IDS_AND_LIMITS)+THE_IDS] == mythrid ) {
			ompshar_initialize_halo(24,  0,0,0,          0,0,0);
		} else {
			o[0] = INVALID_COORDINATE; o[1] = INVALID_COORDINATE; o[2] = INVALID_COORDINATE;
			determine_origins( mythrid, +1, -1, +1,	1, party, partz, o);
			ompshar_initialize_halo(24, 1,1,1,          o[0],o[1],o[2]);
		}

		if ( mynborregions[(25*IDS_AND_LIMITS)+THE_IDS] == mythrid ) {
			ompshar_initialize_halo(25,  0,0,0,          0,0,0);
		} else {
			o[0] = INVALID_COORDINATE; o[1] = INVALID_COORDINATE; o[2] = INVALID_COORDINATE;
			determine_origins( mythrid, -1, +1, +1,	1, party, partz, o);
			ompshar_initialize_halo(25, 1,1,1,          o[0],o[1],o[2]);
		}

		if ( mynborregions[(26*IDS_AND_LIMITS)+THE_IDS] == mythrid ) {
			ompshar_initialize_halo(26,  0,0,0,          0,0,0);
		} else {
			o[0] = INVALID_COORDINATE; o[1] = INVALID_COORDINATE; o[2] = INVALID_COORDINATE;
			determine_origins( mythrid, +1, +1, +1,	1, party, partz, o ); 
			ompshar_initialize_halo(26, 1,1,1,          o[0],o[1],o[2]);
		}
}


double caregionMemHdl::reg_calc_mobilityweight( uint32_t rgpoolid, uint32_t dgpoolid )
{
	//centralize this function for easier maintenance passing theta to another function with implements the mobility models
	//Sebald Gottstein model
	double q1[4], q2[4], theta;

//#pragma omp critical
//{
//	for (int i = 0; i < this->mycaHdl->myrxgpool.size(); i++) { cout << "tid;rxg;i = " << this->mythreadid << ";" << i << "--" << this->mycaHdl->myrxgpool[i].caori << endl; }
//	for (int i = 0; i < this->mycaHdl->mydefgpool.size(); i++) { cout << "tid;defg;i = " << this->mythreadid << ";" << i << "--" << this->mycaHdl->mydefgpool[i].caori << endl; }
//}

	uint32_t rgori = this->mycaHdl->myrxgpool[rgpoolid].caori;
	uint32_t dgori = this->mycaHdl->mydefgpool[dgpoolid].caori;

	//now calculate on the fly...
	q1[0] = this->mycaHdl->myoripool[rgori].q0;
	q1[1] = this->mycaHdl->myoripool[rgori].q1;
	q1[2] = this->mycaHdl->myoripool[rgori].q2;
	q1[3] = this->mycaHdl->myoripool[rgori].q3;

	q2[0] = this->mycaHdl->myoripool[dgori].q0;
	q2[1] = this->mycaHdl->myoripool[dgori].q1;
	q2[2] = this->mycaHdl->myoripool[dgori].q2;
	q2[3] = this->mycaHdl->myoripool[dgori].q3;

	//P=weightedMob for all HAGB already defined
	double weightedMob = 0.0;

	theta = this->mycaHdl->localmath.disori_angle_oriquat_cubic( q1, q2 );

	if ( mycaHdl->mobilitymodel == MOBILITYMODEL_ROLLETTHOLM ) {
		//P value allows to calculate the mobility as often occurring in works by Rollett and Holm via P * mRHHAGB
		weightedMob = 1.0 - (mycaHdl->myPhysData.RH_LAGBHAGBcut * exp( -1.0 * mycaHdl->myPhysData.RH_LAGBHAGBtrans * pow( (theta/MAXDISORI_LAGB2HAGB), mycaHdl->myPhysData.RH_LAGBHAGBexponent ) ));

		return weightedMob;
	}

	//else MOBILITYMODEL_SEBALDGOTTSTEIN

	//assign categorical intrinsic boundary mobility to disorientation among two orientation
	//LAGB detected, overwrite default value
	if( theta <= MAXDISORI_LAGB2HAGB ) {
		weightedMob = -1.0;

		return weightedMob;
	}

	double maxDev40_111 = MAXDISORI_TO_40DEG111;
	double _sqrt3 = 1 / sqrt( 3.0 );
	double oneNinth = 1.0 / 9.0; 
	//the quaternion that describes a 40deg<111> misorientation ##MK::we have to test against all possible ones...
	double qdis[4] = {1.0, 0.0, 0.0, 0.0};
	this->mycaHdl->localmath.disori_fromoriquat_cubic_imm( q1, q2, qdis );

	//##MK::m40_111 must be a misorientation quaternion because it describes a particular boundary situation between two grains A and B in Bunge orientation gA gB
	//##MK::qdis however is also a misorientation quaternion, hence we are here assuming that the method to compute the minimum difference between two disorientation quaternions considering crystal symmetries works for orientation quaternions as for misorientation quaternions...
	double m40_111[4] = { cos( oneNinth * _PI_ ), _sqrt3 * sin( oneNinth * _PI_ ), _sqrt3 * sin( oneNinth * _PI_ ), _sqrt3 * sin( oneNinth * _PI_ ) };

	//obviously a HAGB but also eligible for GS? check proximity to 40deg111
	double dev_40_111 = this->mycaHdl->localmath.disori_angle_oriquat_cubic( qdis, m40_111 );

	if( dev_40_111 < maxDev40_111 ) {
		weightedMob = SQR( cos( 0.5 * _PI_ * dev_40_111 / maxDev40_111 ) );
	}

	return weightedMob;
}


void caregionMemHdl::ompshar_init_mySeedFrontInside( void )
{
	uint32_t ninitialSeedFront = (uint32_t) ( mycaHdl->initialRelCellCaching * (double) myGeom.nregvol_rdtdnd );
	QUICKASSERT ( ninitialSeedFront < CA_ALLOCATION_MAXIMUM );
	mySeedFrontInside = NULL;
	mySeedFrontInside = new struct defcell[ninitialSeedFront];
	QUICKASSERT( mySeedFrontInside != NULL );
	regMemGuard = regMemGuard + (ninitialSeedFront * sizeof(defcell));

	ntotalSeedFrontInside = ninitialSeedFront;
	nextSlotNeverActiveSeedInside = 0;
	ncurrentActiveSeedInside = 0;

	//initialize already associated cells
	for ( uint32_t cfr = 0; cfr < ninitialSeedFront; ++cfr ) {
		mySeedFrontInside[cfr].activity = INACTIVE;
		mySeedFrontInside[cfr].infector = 26;
		mySeedFrontInside[cfr].ix = -1;
		mySeedFrontInside[cfr].iy = -1;
		mySeedFrontInside[cfr].iz = -1;
		mySeedFrontInside[cfr].frac = NO_INFECTION;
		mySeedFrontInside[cfr].mydefgseedid = NO_GRAIN_ASSIGNED;
	}


	//allocate memory to store RecycledCells
	uint32_t ninitialRecyclingList = (uint32_t) ( mycaHdl->initialRelCellCaching * (double) myGeom.nregvol_rdtdnd * mycaHdl->maxfillperstep * FULLRECYCLING_ATT_FACTOR);
	QUICKASSERT ( ninitialRecyclingList < CA_ALLOCATION_MAXIMUM );
	myRecyclingSeedInside = NULL;
	myRecyclingSeedInside = new uint32_t[ninitialRecyclingList];
	QUICKASSERT( myRecyclingSeedInside != NULL );
	regMemGuard = regMemGuard + ( ninitialRecyclingList * sizeof(uint32_t) );

	ntotalRecyclingSeedInside = ninitialRecyclingList;
	nextToRecycleSeedInside = NOTHING_TO_RECYCLE; //not assigned as there is no cell at all yet  that can be recycled 
	firstNotRecycleSeedInside = NOTHING_TO_RECYCLE;

	//myRecyclingSeedInside needs no initialization because it is assured to have no gaps and only read [0;nextToRecycleSeedInside)

	//allocate maxfillperstep the initial memory for cell that complete their infection in a timestep
	//this list serves as a guide to reduce the amount of ifs during finding recyclable sites during the infection phase as the RXFrontList becomes for sure fragmented and thus INACTIVE CELLs will likely be encountered
	uint32_t ninitialFullSeedList = (uint32_t) ( mycaHdl->initialRelCellCaching * (double) myGeom.nregvol_rdtdnd * mycaHdl->maxfillperstep * FULLRXCACHING_ATT_FACTOR);
	QUICKASSERT ( ninitialFullSeedList < CA_ALLOCATION_MAXIMUM );
	myFullSeedInside = NULL;
	myFullSeedInside = new uint32_t[ninitialFullSeedList];
	QUICKASSERT( myFullSeedInside != NULL );
	regMemGuard = regMemGuard + ( ninitialFullSeedList * sizeof(uint32_t) );

	ntotalFullSeedInside = ninitialFullSeedList;
	nextSlotToFullSeedInside = 0; //as above not assigned yet

	//no initialization as well assured gap-free read-only [0;nextSlotToFullSeedInside)

	//cout << this->jobid << "\t\tMemGuard (Byte)= " << this->myMemGuard << endl;
}


void caregionMemHdl::ompshar_init_mySeedFrontBorder( void )
{
	//MK::it can happen that a thread infects from its border into the border of an adjacent thread memory region,
	//thus memory in the mySeedFrontBorder of thread B is required to store the cell, however that array might not have been allocated large enough
	//then thread A would be required to issue or request at least a reallocation such that thread B enlarges its mySeedFrontBorder array
	//to avoid this we preallocate in advance already for all possible cells in the outer shell and as such the list of border cells does not require a reallocation
	//this comes at the cost of more main memory utilization but potentially better cache linearity..., however, for large SUs, i.e. where this may be of relevance, 
	//the interface area to volume is low consult PhD thesis by M. K\"uhbach for finding substantiated arguments as to why this allocation is not of a significant problem!

	uint32_t nSeedFront = 2*( (myGeom.nreg_rd * myGeom.nreg_td) + (myGeom.nreg_rd * myGeom.nreg_nd) + (myGeom.nreg_td * myGeom.nreg_nd) ); //##MK::an exact calculation could save a few cells here...
	QUICKASSERT ( nSeedFront < CA_ALLOCATION_MAXIMUM );
	mySeedFrontBorder = NULL;
	mySeedFrontBorder = new struct defcell[nSeedFront];
	QUICKASSERT( mySeedFrontBorder != NULL );
	regMemGuard = regMemGuard + (nSeedFront * sizeof(defcell));

	ntotalSeedFrontBorder = nSeedFront;
	nextSlotNeverActiveSeedBorder = 0;
	ncurrentActiveSeedBorder = 0;

	//initialize already associated cells
	for ( uint32_t cfr = 0; cfr < nSeedFront; ++cfr ) {
		mySeedFrontBorder[cfr].activity = INACTIVE;
		mySeedFrontBorder[cfr].infector = 26;
		mySeedFrontBorder[cfr].ix = -1;
		mySeedFrontBorder[cfr].iy = -1;
		mySeedFrontBorder[cfr].iz = -1;
		mySeedFrontBorder[cfr].frac = NO_INFECTION;
		mySeedFrontBorder[cfr].mydefgseedid = NO_GRAIN_ASSIGNED;
	}


	//allocate memory to store RecycledCells, list does not grow!
	uint32_t nRecyclingList = (uint32_t) ( (double) nSeedFront * mycaHdl->maxfillperstep * FULLRECYCLING_ATT_FACTOR);
	QUICKASSERT ( nRecyclingList < CA_ALLOCATION_MAXIMUM );
	myRecyclingSeedBorder = NULL;
	myRecyclingSeedBorder = new uint32_t[nRecyclingList];
	QUICKASSERT( myRecyclingSeedBorder != NULL );
	regMemGuard = regMemGuard + ( nRecyclingList * sizeof(uint32_t) );

	ntotalRecyclingSeedBorder = nRecyclingList;
	nextToRecycleSeedBorder = NOTHING_TO_RECYCLE; //not assigned as there is no cell at all yet  that can be recycled 
	firstNotRecycleSeedBorder = NOTHING_TO_RECYCLE;

	//myRecyclingSeedBorder needs no initialization because it is assured to have no gaps and only read [0;nextToRecycleSeedBorder)

	//allocate maxfillperstep the initial memory for cell that complete their infection in a timestep
	//this list serves as a guidance to reduce the amount of ifs during the infection phase as the list is for sure fragmented and thus INACTIVE CELLs are  met
	uint32_t nFullSeedList = (uint32_t) ( (double) nSeedFront * mycaHdl->maxfillperstep * FULLRXCACHING_ATT_FACTOR);
	QUICKASSERT ( nFullSeedList < CA_ALLOCATION_MAXIMUM );
	myFullSeedBorder = NULL;
	myFullSeedBorder = new uint32_t[nFullSeedList];
	QUICKASSERT( myFullSeedBorder != NULL );
	regMemGuard = regMemGuard + ( nFullSeedList * sizeof(uint32_t) );

	ntotalFullSeedBorder = nFullSeedList;
	nextSlotToFullSeedBorder = 0; //as above not assigned yet

	//no initialization as well  assured gap-free read-only [0;nextSlotToFullSeedBorder)

	//cout << this->jobid << "\t\tMemGuard (Byte)= " << this->myMemGuard << endl;
}


uint32_t caregionMemHdl::omp_getNextFreeSlotInSeedFrontInside( bool how )
{
	//CALLED FROM WITHIN PARALLEL REGION, MUST NOT BE CALLED WHEN NOTHING WAS RECYCLED DURING A CALCGROWTHSTEP
	//see more detailed comments in omp_getNextFreeSlotInRXFrontInside

	bool strategy = how;
	uint32_t place = INVALID_ADDRESS;

	//if recycling is desired
	if ( strategy == CELLRECYCLE ) {
		if ( nextToRecycleSeedInside < firstNotRecycleSeedInside ) {
			place = myRecyclingSeedInside[nextToRecycleSeedInside];
			nextToRecycleSeedInside++;
			return place;
		}

		//hasnt returned yet, so recycling list should be utilized but was already exhausted so change the strategy to get memory
		strategy = CELLAPPEND;
	}

	//now strategy is for sure CELLAPPEND
	if ( nextSlotNeverActiveSeedInside < ntotalSeedFrontInside ) {
		place = nextSlotNeverActiveSeedInside;
		nextSlotNeverActiveSeedInside++;
		return place;
	}

	//obviously	, recycling list already exhausted and precached places all occupied. Hence, we need to expand mySeedFrontInside!
	size_t oldSize = ntotalSeedFrontInside;
	size_t newSize = oldSize + (mycaHdl->transientRelCellRecaching * oldSize);
	QUICKASSERT( newSize < CA_ALLOCATION_MAXIMUM );

//#ifdef REPORTSTYLE_DEVELOPER
//	cout << "\t\tALLOCATION for new_mySeedFront;oldSize;newSize\t" << ((int) oldSize) << "\t" << ((int) newSize) << endl;
//#endif

	defcellP new_mySeedFrontInside = NULL;
	new_mySeedFrontInside = new struct defcell[newSize];
	QUICKASSERT( new_mySeedFrontInside != NULL );
	regMemGuard = regMemGuard + (newSize * sizeof(defcell));

	//MK::all indices in FullSeed and Recycling list stay valid as they are referencing relative to begin of mySeedFrontInside array
	for ( uint32_t i = 0; i < oldSize; ++i ) {
		new_mySeedFrontInside[i].activity = mySeedFrontInside[i].activity;
		new_mySeedFrontInside[i].infector = mySeedFrontInside[i].infector;
		new_mySeedFrontInside[i].ix = mySeedFrontInside[i].ix;
		new_mySeedFrontInside[i].iy = mySeedFrontInside[i].iy;
		new_mySeedFrontInside[i].iz = mySeedFrontInside[i].iz;
		new_mySeedFrontInside[i].frac = mySeedFrontInside[i].frac;
		new_mySeedFrontInside[i].mydefgseedid = mySeedFrontInside[i].mydefgseedid;
	}

	delete [] mySeedFrontInside;
	regMemGuard = regMemGuard - (oldSize * sizeof(defcell));
	mySeedFrontInside = new_mySeedFrontInside;
	ntotalSeedFrontInside = newSize;

	//initialize this new memory snippet, access with nc = oldSize is to 0-th element of appended section
	for ( uint32_t nc = oldSize; nc < ntotalSeedFrontInside; ++nc) {
		mySeedFrontInside[nc].activity = INACTIVE;
		mySeedFrontInside[nc].infector = 26;
		mySeedFrontInside[nc].ix = -1;
		mySeedFrontInside[nc].iy = -1;
		mySeedFrontInside[nc].iz = -1;
		mySeedFrontInside[nc].frac = NO_INFECTION;
		mySeedFrontInside[nc].mydefgseedid = NO_GRAIN_ASSIGNED;
	}

	//now we have new space, so append
	place = nextSlotNeverActiveSeedInside;
	nextSlotNeverActiveSeedInside++;
	return place;
}


uint32_t caregionMemHdl::omp_getNextFreeSlotInSeedFrontBorder( bool how )
{
	//CALLED FROM WITHIN PARALLEL REGION, MUST NOT BE CALLED WHEN NOTHING WAS RECYCLED DURING A CALCGROWTHSTEP
	//MK::contrary to getNextFreeSlotInSeedInside we are now guaranteed to have a large enough container mySeedFrontBorder to avoid a potential reallocation
	//algorithm similar to omp_getNextFreeSlotInSeedFrontInside(), so see further comments there...
	bool strategy = how;
	uint32_t place = INVALID_ADDRESS;

	if ( strategy == CELLRECYCLE ) {
		if ( nextToRecycleSeedBorder < firstNotRecycleSeedBorder ) {
			place = myRecyclingSeedBorder[nextToRecycleSeedBorder];
			nextToRecycleSeedBorder++;
			return place;
		}

		//hasnt returned yet, so recycling list should be utilized but was already exhausted so change the strategy to get memory
		strategy = CELLAPPEND;
	}

	//now for sure strategy is CELLAPPEND
	if ( nextSlotNeverActiveSeedBorder < ntotalSeedFrontBorder ) {
		place = nextSlotNeverActiveSeedBorder;
		nextSlotNeverActiveSeedBorder++;
	}

	//return an address, maybe one that is marked as invalid
	return place;
}


void caregionMemHdl::ompshar_voxelize_growthStepInside( void )
{
	double tprof = omp_get_wtime();
	double shapefactors[27] = { FGEODIAG, FGEOEDGE, FGEODIAG, FGEOEDGE, FGEOFACE, FGEOEDGE, FGEODIAG, FGEOEDGE, FGEODIAG, 
						FGEOEDGE, FGEOFACE, FGEOEDGE, FGEOFACE,      FGEOFACE, FGEOEDGE, FGEOFACE, FGEOEDGE,	
						FGEODIAG, FGEOEDGE, FGEODIAG, FGEOEDGE, FGEOFACE, FGEOEDGE, FGEODIAG, FGEOEDGE, FGEODIAG, (1.0) }; //layers xz stacked in positive y direction...

	uint32_t IdentifiedAsActive = 0;

	//MK::reset collector array for already transformed cells
	nextSlotToFullSeedInside = 0;

	//MK::reset RecyclingCounter as int-values from previous step are now corrupted, up to now nothing has been recycled
	nextToRecycleSeedInside = 0;
	firstNotRecycleSeedInside = 0;

	//scan entire mySeedFrontInside on [0;nextSlotNeverActiveSeedInside) to 
	//	a) to find ACTIVE cells that require simulation of migration increment
	//	b) to find INACTIVE cells which we can recycle

	//MK::DEVELOPER NOTE:when nucleation should be time-dependent it is vital to have the nuclei insert here AT THE LATEST, otherwise nextSlotNeverActiveRX will become updated and as such corrupted...
	double dVrx = 0.0;
	double v = 0.0;
	for ( uint32_t c = 0; c < nextSlotNeverActiveSeedInside; c++ ) {

		//MK::ACTIVE is the most likely case in a well defragmented list
		if ( mySeedFrontInside[c].activity == ACTIVE ) {
			IdentifiedAsActive++;

			//unit time infection fill increment, scale according to infection direction
			v = mycaHdl->maxfillperstep * shapefactors[mySeedFrontInside[c].infector];
			mySeedFrontInside[c].frac = mySeedFrontInside[c].frac + (v * 1.0);
			dVrx = dVrx + (v * 1.0);

			//does this fill the cell up?
			if ( mySeedFrontInside[c].frac >= 1.0 ) { //is the cell finally full?
				this->append_to_fullseed_inside_list ( c );
			}
			continue;
		}

		//else cell seems inactive, okay remember index c in RecyclingList
		this->append_to_recycseed_inside_list ( c );

	} //analyze all cells in [0, nextSlotNeverActiveSeedInside)


	reg_dXstepInside = dVrx;
	reg_SvInside = IdentifiedAsActive;
	reg_nCurrActiveInside = IdentifiedAsActive;
	dtCalcGrowthInside = omp_get_wtime() - tprof;
}


void caregionMemHdl::ompshar_voxelize_growthStepBorder( void )
{
	double tprof = omp_get_wtime();
	double shapefactors[27] = { FGEODIAG, FGEOEDGE, FGEODIAG, FGEOEDGE, FGEOFACE, FGEOEDGE, FGEODIAG, FGEOEDGE, FGEODIAG, 
						FGEOEDGE, FGEOFACE, FGEOEDGE, FGEOFACE,      FGEOFACE, FGEOEDGE, FGEOFACE, FGEOEDGE,	
						FGEODIAG, FGEOEDGE, FGEODIAG, FGEOEDGE, FGEOFACE, FGEOEDGE, FGEODIAG, FGEOEDGE, FGEODIAG, (1.0) }; //layers xz stacked in positive y direction...

	uint32_t IdentifiedAsActive = 0;
	nextSlotToFullSeedBorder = 0;

	//MK::reset RecyclingCounter as int-values from previous step are now corrupted, up to now nothing has been recycled
	nextToRecycleSeedBorder = 0;
	firstNotRecycleSeedBorder = 0;

	//strategy as for the inside
	double dVrx = 0.0;
	double v = 0.0;
	for ( uint32_t c = 0; c < nextSlotNeverActiveSeedBorder; c++ ) {
		if ( mySeedFrontBorder[c].activity == ACTIVE ) {
			IdentifiedAsActive++;

			v = mycaHdl->maxfillperstep * shapefactors[mySeedFrontBorder[c].infector];
			mySeedFrontBorder[c].frac = mySeedFrontBorder[c].frac + (v * 1.0);
			dVrx = dVrx + (v * 1.0);

			if ( mySeedFrontBorder[c].frac >= 1.0 ) {
				this->append_to_fullseed_border_list ( c );
			}
			continue;
		} 
		//else cell == INACTIVE
		this->append_to_recycseed_border_list ( c );

	} //analyze all cells in [0, nextSlotNeverActiveSeedBorder)

	reg_dXstepBorder = dVrx;
	reg_SvBorder = IdentifiedAsActive;
	reg_nCurrActiveBorder = IdentifiedAsActive;
	dtCalcGrowthBorder = omp_get_wtime() - tprof;
}


void caregionMemHdl::ompshar_voxelize_updateFullInside( void )
{
	double tprof = omp_get_wtime();
	//CALLED FROM WITHIN PARALLEL REGION
	//myFullSeedInside in the interval [0;nextSlotToFullSeedInside) is a compact list of entries 
	//dereferencing cells from mySeedFrontInside that in this time step have completely transformed, can infect and be INACTIVED thereafter
	//but it is required their categorization if the infected cell is still in the inside of a region or now handled by the border

	//repartitioning of recrystallized volume, original factors for velocity fgeo 1 in <100>, 1/2^0.5 <110>, 1/3^0.5 <111>
	double overFac1 = 0.0682;
	double overFac2 = 0.0542;
	double overFac3 = 0.0484;
	double overshoot, carryOver1, carryOver2, carryOver3;

	uint32_t nUpdatedCells = 0;
	uint32_t xmi = myGeom.nreg_rdmin;
	uint32_t xmx = myGeom.nreg_rdmax;
	uint32_t ymi = myGeom.nreg_tdmin;
	uint32_t ymx = myGeom.nreg_tdmax;
	uint32_t zmi = myGeom.nreg_ndmin;
	uint32_t zmx = myGeom.nreg_ndmax;
	uint32_t xx = myGeom.nreg_rd;
	uint32_t xxyy = myGeom.nregarea_rdtd;

	//MK::if a cell is located inside the automaton domain and a fixed kernel, e.g. MOORE is utilized
	//in almost all cases the infected cell needs no checking of periodic boundary conditions
	//MK::but we have to distinguish in every case into which region the infection goes
	bool InfectIntoBorder; //if true placing the cell in mySeedFrontBorder, otherwise in mySeedFrontInside
	uint32_t crxfrontid, rgpid, cxyz;

	for ( uint32_t c = 0; c < nextSlotToFullSeedInside; c++ ) {
		crxfrontid = myFullSeedInside[c];

		rgpid = mySeedFrontInside[crxfrontid].mydefgseedid;

		QUICKASSERT( mySeedFrontInside[crxfrontid].frac >= 1.0 );

		overshoot = mySeedFrontInside[crxfrontid].frac - 1.0;
		carryOver1 = overshoot * overFac1; //partition normalized overshoot
		carryOver2 = overshoot * overFac2;
		carryOver3 = overshoot * overFac3;

		//##MK::update potential shape axes-aligned bounding boxes for shape tracking --> inject code here

		//infect all untransformed neighbors
		//as all these infecting cells are located at the inside of the domain they point at most into a border but not a periodic image
		//hence it is not necessary a boundary check for periodicity
		InfectIntoBorder = true;
		if ( mySeedFrontInside[crxfrontid].ix > (xmi + 1) && mySeedFrontInside[crxfrontid].ix < (xmx - 1) && mySeedFrontInside[crxfrontid].iy > (ymi + 1) && mySeedFrontInside[crxfrontid].iy < (ymx - 1) && mySeedFrontInside[crxfrontid].iz > (zmi + 1) && mySeedFrontInside[crxfrontid].iz < (zmx - 1) ) {
			InfectIntoBorder = false;
		}

		this->omp_infect_OutOfSeedFrontInside ( rgpid, crxfrontid, InfectIntoBorder, -1, -1, 0, 16, carryOver2 );
		this->omp_infect_OutOfSeedFrontInside ( rgpid, crxfrontid, InfectIntoBorder, 0, -1, 0, 15, carryOver1 );
		this->omp_infect_OutOfSeedFrontInside ( rgpid, crxfrontid, InfectIntoBorder, +1, -1, 0, 14, carryOver2 );
		this->omp_infect_OutOfSeedFrontInside ( rgpid, crxfrontid, InfectIntoBorder, -1, 0, 0, 13, carryOver1 );
		this->omp_infect_OutOfSeedFrontInside ( rgpid, crxfrontid, InfectIntoBorder, +1, 0, 0, 12, carryOver1 );
		this->omp_infect_OutOfSeedFrontInside ( rgpid, crxfrontid, InfectIntoBorder, -1, +1, 0, 11, carryOver2 );
		this->omp_infect_OutOfSeedFrontInside ( rgpid, crxfrontid, InfectIntoBorder, 0, +1, 0, 10, carryOver1 );
		this->omp_infect_OutOfSeedFrontInside ( rgpid, crxfrontid, InfectIntoBorder, +1, +1, 0, 9, carryOver2 );

		this->omp_infect_OutOfSeedFrontInside ( rgpid, crxfrontid, InfectIntoBorder, -1, +1, -1, 19, carryOver3 );
		this->omp_infect_OutOfSeedFrontInside ( rgpid, crxfrontid, InfectIntoBorder, 0, +1, -1, 18, carryOver2 );
		this->omp_infect_OutOfSeedFrontInside ( rgpid, crxfrontid, InfectIntoBorder, +1, +1, -1, 17, carryOver3 );
		this->omp_infect_OutOfSeedFrontInside ( rgpid, crxfrontid, InfectIntoBorder, -1, -1, -1, 25, carryOver3 );
		this->omp_infect_OutOfSeedFrontInside ( rgpid, crxfrontid, InfectIntoBorder, 0, -1, -1, 24, carryOver2 );
		this->omp_infect_OutOfSeedFrontInside ( rgpid, crxfrontid, InfectIntoBorder, +1, -1, -1, 23, carryOver3 );
		this->omp_infect_OutOfSeedFrontInside ( rgpid, crxfrontid, InfectIntoBorder, -1, 0, -1, 22, carryOver2 );
		this->omp_infect_OutOfSeedFrontInside ( rgpid, crxfrontid, InfectIntoBorder, 0, 0, -1, 21, carryOver1 );
		this->omp_infect_OutOfSeedFrontInside ( rgpid, crxfrontid, InfectIntoBorder, +1, 0, -1, 20, carryOver2 );

		this->omp_infect_OutOfSeedFrontInside ( rgpid, crxfrontid, InfectIntoBorder, -1, -1, +1, 8, carryOver3 );
		this->omp_infect_OutOfSeedFrontInside ( rgpid, crxfrontid, InfectIntoBorder, 0, -1, +1, 7, carryOver2 );
		this->omp_infect_OutOfSeedFrontInside ( rgpid, crxfrontid, InfectIntoBorder, +1, -1, +1, 6, carryOver3 );
		this->omp_infect_OutOfSeedFrontInside ( rgpid, crxfrontid, InfectIntoBorder, -1, 0, +1, 5, carryOver2 );
		this->omp_infect_OutOfSeedFrontInside ( rgpid, crxfrontid, InfectIntoBorder, 0, 0, +1, 4, carryOver1 );
		this->omp_infect_OutOfSeedFrontInside ( rgpid, crxfrontid, InfectIntoBorder, +1, 0, +1, 3, carryOver2 );
		this->omp_infect_OutOfSeedFrontInside ( rgpid, crxfrontid, InfectIntoBorder, -1, +1, +1, 2, carryOver3 );
		this->omp_infect_OutOfSeedFrontInside ( rgpid, crxfrontid, InfectIntoBorder, 0, +1, +1, 1, carryOver2 );
		this->omp_infect_OutOfSeedFrontInside ( rgpid, crxfrontid, InfectIntoBorder, +1, +1, +1, 0, carryOver3 );
		//MK::following code NEEDS to be executed AFTER all these infection attempts

		//deactivate cell 
		cxyz = (mySeedFrontInside[crxfrontid].ix - xmi) + ((mySeedFrontInside[crxfrontid].iy - ymi) * xx) + ((mySeedFrontInside[crxfrontid].iz - zmi) * xxyy);

		//mark representor in the cellgrid as FULLY_RECRYSTALLIZED by assigning id of the recrystallizing
		mycellgrid[cxyz] = rgpid; //##MK::this design limits the maximum number of grains possible

		//inject local texture bookkeeping here

		//cell flagged as "Freiwild" and thought of Green Mile healed from all infections...
		mySeedFrontInside[crxfrontid].activity = INACTIVE;
		mySeedFrontInside[crxfrontid].frac = 0.0;

//#ifdef REPORTSTYLE_CELLCYCLES
//			cout << "\t\tstep->updt_afterinfects;crxfrontid;activity;ACTIVEwouldbeMarkedAs;nextSlotNeverActiveSeedInside\t" << mycaHdl->step << "\t" << crxfrontid << "\t" << mySeedFrontInside[crxfrontid].activity << "\t" << ACTIVE << "\t" << nextSlotNeverActiveSeedInside << endl;
//#endif

		nUpdatedCells++;
	} //for all cells

	reg_nUpdatedCellsInside = nUpdatedCells;
	dtUpdateInside = omp_get_wtime() - tprof;
}


void caregionMemHdl::ompshar_voxelize_updateFullBorder( void )
{
	double tprof = omp_get_wtime();
	//CALLED FROM WITHIN PARALLEL REGION
	//MK::it is principal design the same as the updateFullRXInside, however now we must potentially direct infections 
	//out of our own memory region by offloading them into our halo such that the neighboring threads can take care about them

	double overFac1 = 0.0682;
	double overFac2 = 0.0542;
	double overFac3 = 0.0484;
	double overshoot, carryOver1, carryOver2, carryOver3;

	uint32_t nUpdatedCells = 0;
	uint32_t xmi = myGeom.nreg_rdmin;
	uint32_t ymi = myGeom.nreg_tdmin;
	uint32_t zmi = myGeom.nreg_ndmin;
	uint32_t xx = myGeom.nreg_rd;
	uint32_t xxyy = myGeom.nregarea_rdtd;
	uint32_t crxfrontid, rgpid, cxyz;

//##DEBUG#pragma omp critical { cout << "\t\tEntering updateFullSeedBorder;step;nextSlotToFullSeedBorder\t" << mycaHdl->step << "\t" << nextSlotToFullSeedBorder << endl; }

	for ( uint32_t c = 0; c < nextSlotToFullSeedBorder; c++ ) {
		crxfrontid = myFullSeedBorder[c];

		rgpid = mySeedFrontBorder[crxfrontid].mydefgseedid;

//##DEBUG#pragma omp critical { cout << "\t\t\t->step->Infecting BORDER with;step;c;crxfrontid;mySeedFrontBorder[crxfrontid].frac;nextSlotToFullSeedBorder\t" << mycaHdl->step << "\t" << c << "\t" << crxfrontid << "\t" << mySeedFrontBorder[crxfrontid].frac << "\t" << nextSlotToFullSeedBorder << endl; }

		QUICKASSERT( mySeedFrontBorder[crxfrontid].frac >= 1.0 );

		overshoot = mySeedFrontBorder[crxfrontid].frac - 1.0;
		carryOver1 = overshoot * overFac1; //partition normalized overshoot
		carryOver2 = overshoot * overFac2;
		carryOver3 = overshoot * overFac3;

		//##MK::update potential shape axes-aligned bounding boxes for shape tracking --> inject code here

		//three possible cases where the newly infected cells is located:
		//in the inside of this->region, in the border of this->region or in another region

		this->omp_infect_OutOfSeedFrontBorder ( rgpid, crxfrontid, -1, -1, 0, 16, carryOver2 );
		this->omp_infect_OutOfSeedFrontBorder ( rgpid, crxfrontid, 0, -1, 0, 15, carryOver1 );
		this->omp_infect_OutOfSeedFrontBorder ( rgpid, crxfrontid, +1, -1, 0, 14, carryOver2 );
		this->omp_infect_OutOfSeedFrontBorder ( rgpid, crxfrontid, -1, 0, 0, 13, carryOver1 );
		this->omp_infect_OutOfSeedFrontBorder ( rgpid, crxfrontid, +1, 0, 0, 12, carryOver1 );
		this->omp_infect_OutOfSeedFrontBorder ( rgpid, crxfrontid, -1, +1, 0, 11, carryOver2 );
		this->omp_infect_OutOfSeedFrontBorder ( rgpid, crxfrontid, 0, +1, 0, 10, carryOver1 );
		this->omp_infect_OutOfSeedFrontBorder ( rgpid, crxfrontid, +1, +1, 0, 9, carryOver2 );

		this->omp_infect_OutOfSeedFrontBorder ( rgpid, crxfrontid, -1, +1, -1, 19, carryOver3 );
		this->omp_infect_OutOfSeedFrontBorder ( rgpid, crxfrontid, 0, +1, -1, 18, carryOver2 );
		this->omp_infect_OutOfSeedFrontBorder ( rgpid, crxfrontid, +1, +1, -1, 17, carryOver3 );
		this->omp_infect_OutOfSeedFrontBorder ( rgpid, crxfrontid, -1, -1, -1, 25, carryOver3 );
		this->omp_infect_OutOfSeedFrontBorder ( rgpid, crxfrontid, 0, -1, -1, 24, carryOver2 );
		this->omp_infect_OutOfSeedFrontBorder ( rgpid, crxfrontid, +1, -1, -1, 23, carryOver3 );
		this->omp_infect_OutOfSeedFrontBorder ( rgpid, crxfrontid, -1, 0, -1, 22, carryOver2 );
		this->omp_infect_OutOfSeedFrontBorder ( rgpid, crxfrontid, 0, 0, -1, 21, carryOver1 );
		this->omp_infect_OutOfSeedFrontBorder ( rgpid, crxfrontid, +1, 0, -1, 20, carryOver2 );

		this->omp_infect_OutOfSeedFrontBorder ( rgpid, crxfrontid, -1, -1, +1, 8, carryOver3 );
		this->omp_infect_OutOfSeedFrontBorder ( rgpid, crxfrontid, 0, -1, +1, 7, carryOver2 );
		this->omp_infect_OutOfSeedFrontBorder ( rgpid, crxfrontid, +1, -1, +1, 6, carryOver3 );
		this->omp_infect_OutOfSeedFrontBorder ( rgpid, crxfrontid, -1, 0, +1, 5, carryOver2 );
		this->omp_infect_OutOfSeedFrontBorder ( rgpid, crxfrontid, 0, 0, +1, 4, carryOver1 );
		this->omp_infect_OutOfSeedFrontBorder ( rgpid, crxfrontid, +1, 0, +1, 3, carryOver2 );
		this->omp_infect_OutOfSeedFrontBorder ( rgpid, crxfrontid, -1, +1, +1, 2, carryOver3 );
		this->omp_infect_OutOfSeedFrontBorder ( rgpid, crxfrontid, 0, +1, +1, 1, carryOver2 );
		this->omp_infect_OutOfSeedFrontBorder ( rgpid, crxfrontid, +1, +1, +1, 0, carryOver3 );
		//MK::following code NEEDS to be executed AFTER all these infection attempts

		//deactivate cell, cell is definately within this->mycellgrid
		cxyz = (mySeedFrontBorder[crxfrontid].ix - xmi) + ((mySeedFrontBorder[crxfrontid].iy - ymi) * xx) + ((mySeedFrontBorder[crxfrontid].iz - zmi) * xxyy);

		//mark representor in the cellgrid as FULLY_RECRYSTALLIZED by assigning id of RXgrains
		mycellgrid[cxyz] = rgpid;

		//###WRONG::MK::local texture bookkeepingreg_myrxgp_counts[myRXFrontBorder[crxfrontid].myrxgid] += 1;reg_mydfgp_counts[myRXFrontBorder[crxfrontid].mydefgid] -= 1;

		//cell flagged as "Freiwild" and thought of Green Mile healed from all infections...
		mySeedFrontBorder[crxfrontid].activity = INACTIVE;
		mySeedFrontBorder[crxfrontid].frac = 0.0;

//#ifdef REPORTSTYLE_CELLCYCLES
//			cout << "\t\tstep->updt_afterinfects;crxfrontid;activity;ACTIVEwouldbeMarkedAs;nextSeedSlotNeverActiveBorder\t" << mycaHdl->step << "\t" << crxfrontid << "\t" << mySeedFront[crxfrontid].activity << "\t" << ACTIVE << "\t" << nextSlotNeverActiveSeedBorder << endl;
//#endif

		nUpdatedCells++;
	} //for all cells

	reg_nUpdatedCellsBorder = nUpdatedCells;
	dtUpdateBorder = omp_get_wtime() - tprof;
}


void caregionMemHdl::ompshar_synchronize_haloregions_def( void )
{
	double tprof = omp_get_wtime();

	//scan over all myhalo regions that are expected at all updated
	uint32_t hal, nbortid;
	caregionMemHdlP apartreg;
	uint32_t HowManyToSynchronize;

	uint32_t xmi = myGeom.nreg_rdmin; //inclusive
	uint32_t xmx = myGeom.nreg_rdmax;
	uint32_t ymi = myGeom.nreg_tdmin;
	uint32_t ymx = myGeom.nreg_tdmax;
	uint32_t zmi = myGeom.nreg_ndmin;
	uint32_t zmx = myGeom.nreg_ndmax;
	uint32_t xx = myGeom.nreg_rd;
	uint32_t yy = myGeom.nreg_td;
	uint32_t xxyy = xx*yy;
	uint32_t freeplace;
	uint32_t cxyz;
	uint32_t cxyz_initialstate;
	uint32_t gx, gy, gz;
	//uint32_t oox, ooy, ooz;

	//sweep through my neighbors and synchronize with their halos
	for ( uint32_t nb = 0; nb < (MYSELF + NUMBER_OF_NEIGHBORS); nb++ ) {
		nbortid = mynborregions[(nb*IDS_AND_LIMITS)+THE_IDS];
		apartreg = mycaHdl->regions[nbortid]; //MK::threads can access each others data because of shared-memory!

		hal = apartreg->mynborhalos[nb];
		if ( hal != HALO_NOT_EXISTENT ) { //is there halo information?

//#pragma omp critical
//{ cout << "synchaloinfo threadid;nb;nbortid;hal" << this->mythreadid << ";" << nb << ";" << nbortid << ";" << hal << endl; }

			//get the regionid of the neighbor that has this halo
			HowManyToSynchronize = apartreg->myhalos[hal].nextFreeHaloRefs;
			if ( HowManyToSynchronize < 1 ) { continue; }

			//not continued, so there was sth in the halo requiring sync
			//hrefs = apartreg->myhalos[hal].theHaloRefs;
			//hcells = apartreg->myhalos[hal].theHaloCells;

			//oox = apartreg->myhalos[hal].originx;
			//ooy = apartreg->myhalos[hal].originy;
			//ooz = apartreg->myhalos[hal].originz;
 
			uint32_t target;
			for ( uint32_t s = 0; s < HowManyToSynchronize; s++ ) {
				target = apartreg->myhalos[hal].theHaloRefs[s];

				//read out halo coordinate
				gx = apartreg->myhalos[hal].theHaloCells[target].ix;
				gy = apartreg->myhalos[hal].theHaloCells[target].iy;
				gz = apartreg->myhalos[hal].theHaloCells[target].iz;

				QUICKASSERT( gx >= xmi && gx <= xmx); //##MK::is this position really in the calling threads local memoryregion?
				QUICKASSERT( gy >= ymi && gy <= ymx);
				QUICKASSERT( gz >= zmi && gz <= zmx);
				//transform into local coordinate
				cxyz = (gx - xmi) + ((gy - ymi)*xx) + ((gz - zmi)*xxyy);
				cxyz_initialstate = mycellgrid[cxyz];

				//eine Anekdote zum Thema, es sind nur zwei Zeilen Code.. tausche in folgender Zeile einmal || durch && und debugge mit Totalview, viel Spass
				if ( cxyz_initialstate != NOT_ASSIGNED_YET || cxyz_initialstate == CURRENTLY_INFECTED ) {
					continue;
				}

				//not continued so halocell requires infection border
				mycellgrid[cxyz] = CURRENTLY_INFECTED;

				freeplace = this->omp_getNextFreeSlotInSeedFrontBorder( CELLRECYCLE );
				QUICKASSERT( freeplace != INVALID_ADDRESS );

				mySeedFrontBorder[freeplace].activity = ACTIVE;
				mySeedFrontBorder[freeplace].infector = apartreg->myhalos[hal].theHaloCells[target].infector;
				mySeedFrontBorder[freeplace].ix = apartreg->myhalos[hal].theHaloCells[target].ix;
				mySeedFrontBorder[freeplace].iy = apartreg->myhalos[hal].theHaloCells[target].iy;
				mySeedFrontBorder[freeplace].iz = apartreg->myhalos[hal].theHaloCells[target].iz;
				mySeedFrontBorder[freeplace].frac = apartreg->myhalos[hal].theHaloCells[target].rxFrac; //##MK::identifier naming may be misleading, however halos are such we can reutilize the halo cells
				mySeedFrontBorder[freeplace].mydefgseedid = apartreg->myhalos[hal].theHaloCells[target].myrxgid;

				//in the region where the halo was filled updateFullSeedBorder already marked the infecting cells as recrystallized
				//MK::other bookkeep is not necessary --- in fact --- it were errornous
			} //next cell in the halo
		} //halo handled
	} //check whether next neighbor has a halo to sync for me

	dtSyncHalo = omp_get_wtime() - tprof;
}


void caregionMemHdl::omp_infect_OutOfSeedFrontInside( uint32_t seedid, uint32_t seedfrontid, bool intoborder, short dx, short dy, short dz, unsigned char direction, double carryOver )
{
	//CALLED FROM WITHIN PARALLEL REGION BUT OPERATING ON caregionMemHdl class object's local cellgrid
	short x = mySeedFrontInside[seedfrontid].ix;
	short y = mySeedFrontInside[seedfrontid].iy;
	short z = mySeedFrontInside[seedfrontid].iz;

	//add relative coordinate increment to identify position of infection target in the CA, global coordinates!
	x += dx;
	y += dy;
	z += dz;
	//MK::no bounds check necessary as either the infection is directed into the bulk of the domain (intoborder == false) or its border (intoborder == true)
	//hence x,y,z are positive and thus the promotion from short to uint32_t is not a problem

	//convert ix,iy,iz in implicit 3D coordinate to access mycellgrid
	//MK::WORKS ONLY WITH SHORT IS POSITIVE!
	uint32_t cxyz = (x - myGeom.nreg_rdmin) + ((y - myGeom.nreg_tdmin) * myGeom.nreg_rd) + ((z - myGeom.nreg_ndmin) * myGeom.nregarea_rdtd);

	uint32_t cxyz_initialstate = mycellgrid[cxyz];
	if ( cxyz_initialstate != NOT_ASSIGNED_YET || cxyz_initialstate == CURRENTLY_INFECTED ) {
//##DEBUG#pragma omp critical { cout << "\t\t\tREJECT infection cxyz;x;y;z;dx;dy;dz;deformedstate;reg_nmydefgpool\t" << cxyz << "\t" << x << ";" << y << ";" << z << "\t" << dx << ";" << dy << ";" << dz << "\t" << cxyz_initialstate << "\t" << reg_nmydefgpool << endl; }
		return;
	}

	//obviously cell is free, hehe, go for it!
	mycellgrid[cxyz] = CURRENTLY_INFECTED;

	//MK::now the infection is directed either into the bulk or the border but the infection remains in threadlocal memory
	uint32_t freeplace = INVALID_ADDRESS;
	if ( intoborder == false ) { //most likely case is to infect into the bulk because domains usually cuboids with millions of cells

		freeplace = this->omp_getNextFreeSlotInSeedFrontInside ( CELLRECYCLE ); //MK::strategy is RECYCLING because OutofSeedFront infections are only called from updateFullSeed
//##DEBUG#pragma omp critical { 	cout << "\t\t\tACCEPT-INSIDE at freeplace\t\tcxyz;x;y;z;dx;dy;dz\t" << freeplace << "\t\t" << cxyz << "\t" << x << ";" << y << ";" << z << "\t" << dx << ";" << dy << ";" << dz << endl; }

		mySeedFrontInside[freeplace].activity = ACTIVE;
		mySeedFrontInside[freeplace].infector = direction;
		mySeedFrontInside[freeplace].ix = x;
		mySeedFrontInside[freeplace].iy = y;
		mySeedFrontInside[freeplace].iz = z;
		mySeedFrontInside[freeplace].frac = carryOver;
		mySeedFrontInside[freeplace].mydefgseedid = seedid;

		//MK::assessor function get_NextFreeSlotInSeedFront() ASSURES freeplace < nextSlotNeverActiveSeedFront <= ntotalSeedFront!
		return;
	}

	//not returned yet, okay infection aims for the mySeedFrontBorder
	freeplace = this->omp_getNextFreeSlotInSeedFrontBorder ( CELLRECYCLE ); //MK::strategy is RECYCLING because OutofSeedFront infections are only called from updateFullSeed
//##DEBUG#pragma omp critical { cout << "\t\t\tACCEPT-BORDER at freeplace\t\tcxyz;x;y;z;dx;dy;dz\t" << freeplace << "\t\t" << cxyz << "\t" << x << ";" << y << ";" << z << "\t" << dx << ";" << dy << ";" << dz << endl; }

	mySeedFrontBorder[freeplace].activity = ACTIVE;
	mySeedFrontBorder[freeplace].infector = direction;
	mySeedFrontBorder[freeplace].ix = x;
	mySeedFrontBorder[freeplace].iy = y;
	mySeedFrontBorder[freeplace].iz = z;
	mySeedFrontBorder[freeplace].frac = carryOver;
	mySeedFrontBorder[freeplace].mydefgseedid = seedid;
}


void caregionMemHdl::omp_infect_OutOfSeedFrontBorder( uint32_t seedid, uint32_t seedfrontid, short dx, short dy, short dz, unsigned char direction, double carryOver )
{
	//CALLED FROM WITHIN PARALLEL REGION
	//Three possible locations where the target cell is located:: inside, border, or in another region
	//additionally::we need periodic bound checks to the global domain!

	short x = mySeedFrontBorder[seedfrontid].ix;
	short y = mySeedFrontBorder[seedfrontid].iy;
	short z = mySeedFrontBorder[seedfrontid].iz;

	//add relative coordinate increment to identify position of infection target in the CA, global coordinates!
	x += dx;
	y += dy;
	z += dz;

	//apply periodic boundary conditions at solitary unit global scale
	if ( x < 0 )
		x += myGeom.nedge_global_rd;
	if ( x >= myGeom.nedge_global_rd ) //MK::the short x will for the comparison be implicitly promoted to an unsigned int, which is only what we intend when x is non-negative therefore we have to test x < 0 beforehand!
		x -= myGeom.nedge_global_rd;
	if ( y < 0 )
		y += myGeom.nedge_global_td;
	if ( y >= myGeom.nedge_global_td )
		y -= myGeom.nedge_global_td;
	if ( z < 0 )
		z += myGeom.nedge_global_nd;
	if ( z >= myGeom.nedge_global_nd )
		z -= myGeom.nedge_global_nd;

//##DEBUG#pragma omp critical { cout << x << "\t" << y << "\t" << z << "\t" << seedid << endl; }

	//MK::there are several techniques to identify in which region the infection is, I have chosen for a
	//compromise between speed and clarity instead of making even more complicated branching conditions
	//the problem to solve is: how to get the threadid of the memregion in which the global coordinated x,y,z is located?
	//the strategy: inspect first the most likely, then the less likely and last the least likely, i.e. the regions along the corners and edges of my own region

	//COMPLETELY INSIDE for a cell on the inside of a face of a region this occurs in 9 directions
	if ( x > myGeom.nreg_rdmin && x < myGeom.nreg_rdmax && y > myGeom.nreg_tdmin && y < myGeom.nreg_tdmax && z > myGeom.nreg_ndmin && z < myGeom.nreg_ndmax ) {
		//infection directs into this region, hence requires the activation of a cell from the mySeedFrontInside list

		//convert x,y,z in implicit local 3D coordinates to access mycellgrid, works only because the short's are all positive...
		uint32_t cxyz = (x - myGeom.nreg_rdmin) + ((y - myGeom.nreg_tdmin) * myGeom.nreg_rd) + ((z - myGeom.nreg_ndmin) * myGeom.nregarea_rdtd);

		uint32_t cxyz_initialstate = this->mycellgrid[cxyz];

//##DEBUG#pragma omp critical { cout << "\t\t\t" << omp_get_thread_num() << ";" << "SAMEREGION-INSIDE x,y,z,dx,dy,dz,cxyz,cxyz_initialstate" << x << ";" << y << ";" << z << ";" << cxyz << ";" << cxyz_initialstate << endl; }

		if ( cxyz_initialstate != NOT_ASSIGNED_YET || cxyz_initialstate == CURRENTLY_INFECTED ) { //get out here, cell is already infected
//##DEBUG#pragma omp critical{ cout << "\t\t\t" << omp_get_thread_num() << ";" << "SAMEREGION-INSIDE returned" << endl; }
			return;
		}

		this->mycellgrid[cxyz] = CURRENTLY_INFECTED; 
		//MK::thread-safe because non-overlapping regions plus sequential processing of border cells, for
		//##MK::these cases in an even more elaborated preprocessing the sequential OpenMP overhead could be reduced further...

		uint32_t freeplace = this->omp_getNextFreeSlotInSeedFrontInside ( CELLRECYCLE ); //MK::strategy is RECYCLING because OutofSeedFront infections are only called from updateFullSeed

		mySeedFrontInside[freeplace].activity = ACTIVE;
		mySeedFrontInside[freeplace].infector = direction;
		mySeedFrontInside[freeplace].ix = x;
		mySeedFrontInside[freeplace].iy = y;
		mySeedFrontInside[freeplace].iz = z;
		mySeedFrontInside[freeplace].frac = carryOver;
		mySeedFrontInside[freeplace].mydefgseedid = seedid;
//##DEBUG#pragma omp critical{ cout << "\t\t\t" << omp_get_thread_num() << ";" << "SAMEREGION-BORDER success" << endl; }
		return;
	} 

	//infection directed NOT INSIDE SO NOW NECESSARY TO CHECK AGAINST LIMITS OF ALL POSSIBLE NEIGHBORS, still however it could be an infection that stays within the border layer
	//THIS REGION encloses all infections that direct not out of the shell, i.e. 8 of 26 directions
	//albeit there is not a simple test, like ^sum di^2 == 2 to check for this as the cells can be also located on the edge of the region
	//but lets assume for a moment that the infection is still on this border, it cannot be inside the region as that case was already checked for
	//so check this and confirm that infection is not in another region
	uint32_t nbortid = UNKNOWN_NEIGHBOR;
	//uint32_t nborpos = UNKNOWN_POSITION;

	for ( uint32_t nb = 0; nb < (MYSELF + NUMBER_OF_NEIGHBORS); nb++ ) {
		//MK::be careful when considering to change these seemingly naive if loop into a supposedly "more performant" inversion of the test logic, because
		//the mynborregions array is order such that when the cell is on its own region's boundary the first iteration is directly a hit
		//otherwise the more likely cases (face, corners) are tested for first, hence we almost never run this loop entirely for MYSELF + NUMBER_OF_NEIGHBORS iterations...
		if ( x >= mynborregions[(nb*IDS_AND_LIMITS)+THE_XMIN] && x <= mynborregions[(nb*IDS_AND_LIMITS)+THE_XMAX] && y >= mynborregions[(nb*IDS_AND_LIMITS)+THE_YMIN] && y <= mynborregions[(nb*IDS_AND_LIMITS)+THE_YMAX] && z >= mynborregions[(nb*IDS_AND_LIMITS)+THE_ZMIN] && z <= mynborregions[(nb*IDS_AND_LIMITS)+THE_ZMAX] )	{
				nbortid = mynborregions[(nb*IDS_AND_LIMITS)+THE_IDS];
				//nborpos = nb;
				break;
		}
	}

	//now nbortid holds into which border list we have to write
	//still in my own region but on the boundary?
	if ( nbortid == mynborregions[(0*IDS_AND_LIMITS)+THE_IDS] ) {

		uint32_t cxyz = (x - myGeom.nreg_rdmin) + ((y - myGeom.nreg_tdmin) * myGeom.nreg_rd) + ((z - myGeom.nreg_ndmin) * myGeom.nregarea_rdtd);

		uint32_t cxyz_initial = mycellgrid[cxyz];

//##DEBUG#pragma omp critical { cout << "\t\t\t" << omp_get_thread_num() << ";" << "SAMEREGION-BORDER x,y,z,dx,dy,dz,cxyz,cxyz_deformedstate;nbortid;" << x << ";" << y << ";" << z << ";" << dx << ";" << dy << ";" << dz << ";" << cxyz << ";" << cxyz_deformedstate << ";" << nbortid << endl; }

		if ( cxyz_initial != NOT_ASSIGNED_YET || cxyz_initial == CURRENTLY_INFECTED ) { //get out here, cell is already infected
//##DEBUG#pragma omp critical { cout << "\t\t\t" << omp_get_thread_num() << ";" << "SAMEREGION-BORDER returned" << endl; }
			return;
		}

		mycellgrid[cxyz] = CURRENTLY_INFECTED;

		uint32_t freeplace = this->omp_getNextFreeSlotInSeedFrontBorder ( CELLRECYCLE ); //MK::strategy is RECYCLING because OutofRXFront infections are only called from updateFullRX

		mySeedFrontBorder[freeplace].activity = ACTIVE;
		mySeedFrontBorder[freeplace].infector = direction;
		mySeedFrontBorder[freeplace].ix = x;
		mySeedFrontBorder[freeplace].iy = y;
		mySeedFrontBorder[freeplace].iz = z;
		mySeedFrontBorder[freeplace].frac = carryOver;
		mySeedFrontBorder[freeplace].mydefgseedid = seedid;
//##DEBUG#pragma omp critical{ cout << "\t\t\t" << omp_get_thread_num() << ";" << "SAMEREGION-BORDER success" << endl; }
		return;
	}


	//not returned yet, so the infection directs into a halo check in which and whether the corresponding location is 
	//free nborpos must not be utilized! because multiple halos can lay in different domains!
	//find which halo of mine to place the infection
	uint32_t whichhalo = UNKNOWN_HALO;
	uint32_t hxmi, hxmx, hymi, hymx, hzmi, hzmx;
	for ( uint32_t h = 1; h < (MYSELF + NUMBER_OF_NEIGHBORS); h++ ) {
		if ( myhalos[h].inregion == nbortid ) {
			hxmi = myhalos[h].originx;
			hymi = myhalos[h].originy;
			hzmi = myhalos[h].originz;
			hxmx = hxmi + myhalos[h].extendx;
			hymx = hymi + myhalos[h].extendy;
			hzmx = hzmi + myhalos[h].extendz;
			if ( x >= hxmi && x <= hxmx && y >= hymi && y <= hymx && z >= hzmi && z <= hzmx ) { //##MK::himx are exclusive
				whichhalo = h;
				break;
			}
		}
	}
	QUICKASSERT( whichhalo != UNKNOWN_HALO );

	//is there a halo region for this neighbor?
	//QUICKASSERT ( myhalos[nborpos].inregion == mynborregions[(nborpos*IDS_AND_LIMITS)+THE_IDS] );
	//QUICKASSERT ( myhalos[nborpos].extendxyz > 0 );
	//there is a region so transform global coordinate into local halo coordinate

	short hx = x - myhalos[whichhalo].originx;
	short hy = y - myhalos[whichhalo].originy;
	short hz = z - myhalos[whichhalo].originz;
 
//#pragma omp critical
//{ cout << "issuing from my/nbtid/nbpos = " << this->mythreadid << "--" << nbortid << "--" << whichhalo << "==" << x << ";" << y << ";"<< z << "\t\t" << myhalos[whichhalo].originx << ";" << myhalos[whichhalo].originy << ";" << myhalos[whichhalo].originz << "\t\t" << myhalos[whichhalo].extendx << ";" << myhalos[whichhalo].extendy << ";" << myhalos[whichhalo].extendz << "\t\t" << hx << ";" << hy << ";" << hz << endl; }

	//##DEBUG
	QUICKASSERT ( hx >= 0 && hx < myhalos[whichhalo].extendx );
	QUICKASSERT ( hy >= 0 && hy < myhalos[whichhalo].extendy );
	QUICKASSERT ( hz >= 0 && hz < myhalos[whichhalo].extendz );

	uint32_t hid = hx + (hy * myhalos[whichhalo].extendx) + (hz * myhalos[whichhalo].extendxy);

	//this region does not know what is at the other side unless it checks in the memory shared with the thread
	//but as we are in sharedmemory one can check whether the halo was already infected
	if ( myhalos[whichhalo].theHaloCells[hid].status == HALOCELL_ALREADY_INFECTED ) {
		return; //infection unsuccessful
	}

	//HALOCELL is free mark the halocell, there is enough space and no reallocation necessary because the halo was preallocated
	myhalos[whichhalo].theHaloCells[hid].status = HALOCELL_ALREADY_INFECTED;
	myhalos[whichhalo].theHaloCells[hid].infector = direction;
	myhalos[whichhalo].theHaloCells[hid].ix = x;
	myhalos[whichhalo].theHaloCells[hid].iy = y;
	myhalos[whichhalo].theHaloCells[hid].iz = z;
	myhalos[whichhalo].theHaloCells[hid].rxFrac = carryOver;
	myhalos[whichhalo].theHaloCells[hid].myrxgid = seedid; //MK::identifier myrxgid is misleading only because we utilize the halo twice!

	//make known that in this halo an infection was done
	uint32_t freeh = myhalos[whichhalo].nextFreeHaloRefs;
	myhalos[whichhalo].nextFreeHaloRefs = myhalos[whichhalo].nextFreeHaloRefs + 1;
	myhalos[whichhalo].theHaloRefs[freeh] = hid;

//#pragma omp critical{ cout << "\t\t\t" << omp_get_thread_num() << ";" << "OTHERREGION success" << endl; }
}


void caregionMemHdl::append_to_fullseed_inside_list( uint32_t cid )
{
	if ( nextSlotToFullSeedInside < ntotalFullSeedInside ) {		//FullSeed list still large enough, append
		myFullSeedInside[nextSlotToFullSeedInside] = cid;			//indexing relative to first element in array with displacement sizeof(int)
		nextSlotToFullSeedInside++;
		return;
	}

	//FullSeed list is too small but it is necessary to store all cells that updateFullSeed should work on, so increase size of FullSeed list
	size_t oldSize = ntotalFullSeedInside;
	size_t newSize = oldSize + (mycaHdl->transientRelCellRecaching * oldSize);
	QUICKASSERT ( newSize < CA_ALLOCATION_MAXIMUM );

//#ifdef REPORTSTYLE_DEVELOPER
//	cout << "\t\tALLOCATION for fullseed list, oldSize, newSize\t" << ((int) oldSize) << "\t" << ((int) newSize) << endl;
//#endif

	uint32_t* new_myFullSeedInside = NULL;
	new_myFullSeedInside = new uint32_t[newSize];
	QUICKASSERT( new_myFullSeedInside != NULL );
	regMemGuard = regMemGuard + (newSize * sizeof(uint32_t));

	//##MK::potentially unsafe to use memcpy memcpy( new_myFullSeedInside, myFullSeedInside, oldSize * sizeof(int) );
	for ( uint32_t i = 0; i < oldSize; ++i ) {
		new_myFullSeedInside[i] = myFullSeedInside[i];
	}

	delete [] myFullSeedInside;
	regMemGuard = regMemGuard - ( oldSize * sizeof(uint32_t) );
	myFullSeedInside = new_myFullSeedInside;
	ntotalFullSeedInside = newSize;

	QUICKASSERT ( nextSlotToFullSeedInside < ntotalFullSeedInside );

	//more space again, so append
	myFullSeedInside[nextSlotToFullSeedInside] = cid;
	nextSlotToFullSeedInside++;
}


void caregionMemHdl::append_to_recycseed_inside_list( uint32_t cid )
{
	//firstNotRecycledYetPoints to position in myRecyclingList in [0, ntotalRecyclingListInside) no cell reference to myRXFrontInside was stored yet
	if ( firstNotRecycleSeedInside < ntotalRecyclingSeedInside ) {
		myRecyclingSeedInside[firstNotRecycleSeedInside] = cid;
		firstNotRecycleSeedInside++;
	}
	//else -> currently I do not recycle more than ntotalRecyclingList elements, also this list is already as large as 0.15*ntotalcells, user can optimize this structure when knowing his microstructural path function Sv(X) in more details to save additional memory therefore not further memory utilized
}


void caregionMemHdl::append_to_fullseed_border_list( uint32_t cid )
{
	if ( nextSlotToFullSeedBorder < ntotalFullSeedBorder ) {		//FullSeed list still large enough, append
		myFullSeedBorder[nextSlotToFullSeedBorder] = cid;		//indexing relative to first element in array with displacement sizeof(int)
		nextSlotToFullSeedBorder++;
		return;
	}

	//FullSeed list is too small but it is necessary to store all cells that updateFullSeed should work on, so increase size of FullRX list
	size_t oldSize = ntotalFullSeedBorder;
	size_t newSize = oldSize + (mycaHdl->transientRelCellRecaching * oldSize);
	QUICKASSERT ( newSize < CA_ALLOCATION_MAXIMUM );

//#ifdef REPORTSTYLE_DEVELOPER
	cout << "\t\tALLOCATION for fullseed border list, oldSize, newSize\t" << ((int) oldSize) << "\t" << ((int) newSize) << endl;
//#endif

	uint32_t* new_myFullSeedBorder = NULL;
	new_myFullSeedBorder = new uint32_t[newSize];
	QUICKASSERT( new_myFullSeedBorder != NULL );
	regMemGuard = regMemGuard + (newSize * sizeof(uint32_t));

	for ( uint32_t i = 0; i < oldSize; ++i ) {
		new_myFullSeedBorder[i] = myFullSeedBorder[i];
	}

	delete [] myFullSeedBorder;
	regMemGuard = regMemGuard - ( oldSize * sizeof(uint32_t) );
	myFullSeedBorder = new_myFullSeedBorder;
	ntotalFullSeedBorder = newSize;

	QUICKASSERT ( nextSlotToFullSeedBorder < ntotalFullSeedBorder );

	//more space again, so append
	myFullSeedBorder[nextSlotToFullSeedBorder] = cid;
	nextSlotToFullSeedBorder++;
}


void caregionMemHdl::append_to_recycseed_border_list( uint32_t cid )
{
	if ( firstNotRecycleSeedBorder < ntotalRecyclingSeedBorder ) {
		myRecyclingSeedBorder[firstNotRecycleSeedBorder] = cid;
		firstNotRecycleSeedBorder++;
	}
	//else -> currently I do not recycle more than ntotalRecyclingList elements, also this list is already as large as 0.15*ntotalcells, user can optimize this structure when knowing his microstructural path function Sv(X) in more details to save additional memory therefore not further memory utilized
}


void caregionMemHdl::omp_resetHaloRegions( void )
{
	//reset the halos to the initial stage such that for the solve_RXGROWTH automaton they can be reutilized
	for (uint32_t h = 0; h < (MYSELF + NUMBER_OF_NEIGHBORS); h++ ) { //MK::YOU MUST NOT DELETE THE HALOS AS THEY ARE REUTILIZED IN solve_RXGROWTH!
		for ( uint32_t i = 0; i < myhalos[h].extendxyz; i++ ) { //reset content of all non-dummy halos because a dummy halo has extendxyz == 0
			myhalos[h].theHaloCells[i].status = HALOCELL_NEVER_VISITED;
			myhalos[h].theHaloCells[i].infector = 26;
			myhalos[h].theHaloCells[i].ix = -1;
			myhalos[h].theHaloCells[i].iy = -1;
			myhalos[h].theHaloCells[i].iz = -1;
			myhalos[h].theHaloCells[i].rxFrac = NO_INFECTION;
			myhalos[h].theHaloCells[i].myrxgid = NO_GRAIN_ASSIGNED;
		}
		//resetting of theHaloRefs not necessary, the array myhalos[h].theHaloRefs and myhalos[h].theHaloCells ONE MUST NOT DELETE! but as they are 
		myhalos[h].nextFreeHaloRefs = 0;
		//ONE MUST NOT OVERWRITE THE VALUES FOR extend.. AND origin.. OR inregion AS THE PARTITIONING in solve_RXGROWTH remains the same!
	}
}


void caregionMemHdl::omp_resetInternalCounters( void )
{
	mycaHdl->nmynuclei = 0;

	reg_nCurrActiveInside = 0;
	reg_nCurrActiveBorder = 0;
	reg_SvInside = 0;
	reg_SvBorder = 0;
	reg_dXstepInside = 0.0;
	reg_dXstepBorder = 0.0;

	reg_myMobilityWeightMax = DEFAULT_PMAX;
	reg_nmynuclei = 0;
	mycaHdl->myMobilityWeightMax = DEFAULT_PMAX;


	//MK::THESE MUST NOT BE DELETED!
	/*delete [] reg_myrxgp_counts;
	reg_myrxgp_counts = NULL;
	reg_nmyrxgp_counts = 0;

	delete [] reg_mydfgp_counts;
	reg_mydfgp_counts = NULL;
	reg_nmydfgp_counts = 0;*/
}


void caregionMemHdl::omp_voxelize_memoryCleanup( void ) 
{
	if ( mySeedFrontInside != NULL ) {
		delete [] mySeedFrontInside;
		mySeedFrontInside = NULL;
		ntotalSeedFrontInside = 0;
		nextSlotNeverActiveSeedInside = 0;
		ncurrentActiveSeedInside = 0;
	}

	if ( mySeedFrontBorder != NULL ) {
		delete [] mySeedFrontBorder;
		mySeedFrontBorder = NULL;
		ntotalSeedFrontBorder = 0;
		nextSlotNeverActiveSeedBorder = 0;
		ncurrentActiveSeedBorder = 0;
	}

	if ( myFullSeedInside != NULL ) {
		delete [] myFullSeedInside;
		myFullSeedInside = NULL;
		ntotalFullSeedInside = 0;
		nextSlotToFullSeedInside = 0;
	}

	if ( myFullSeedBorder != NULL ) {
		delete [] myFullSeedBorder;
		myFullSeedBorder = NULL;
		ntotalFullSeedBorder = 0;
		nextSlotToFullSeedBorder = 0;
	}

	if ( myRecyclingSeedInside != NULL ) {
		delete [] myRecyclingSeedInside;
		myRecyclingSeedInside = NULL;
		ntotalRecyclingSeedInside = 0;
		nextToRecycleSeedInside = 0;
		firstNotRecycleSeedInside = 0;
	}

	if ( myRecyclingSeedBorder != NULL ) {
		delete [] myRecyclingSeedBorder;
		myRecyclingSeedBorder = NULL;
		ntotalRecyclingSeedBorder = 0;
		nextToRecycleSeedBorder = 0;
		firstNotRecycleSeedBorder = 0;
	}
}




void caregionMemHdl::ompshar_init_myRXFrontInside( void )
{
	uint32_t ninitialRXFront = (uint32_t) ( mycaHdl->initialRelCellCaching * (double) myGeom.nregvol_rdtdnd );
	QUICKASSERT ( ninitialRXFront < CA_ALLOCATION_MAXIMUM );
	myRXFrontInside = NULL;
	myRXFrontInside = new struct cell[ninitialRXFront];
	QUICKASSERT( myRXFrontInside != NULL );
	regMemGuard = regMemGuard + (ninitialRXFront * sizeof(cell));

	ntotalRXFrontInside = ninitialRXFront;
	nextSlotNeverActiveRXFrontInside = 0;
	nCurrentlyActiveInside = 0;

	//initialize already associated cells
	for ( uint32_t cfr = 0; cfr < ninitialRXFront; ++cfr ) {
		myRXFrontInside[cfr].activity = INACTIVE;
		myRXFrontInside[cfr].infector = 26;
		myRXFrontInside[cfr].ix = -1;
		myRXFrontInside[cfr].iy = -1;
		myRXFrontInside[cfr].iz = -1;
		myRXFrontInside[cfr].rxFrac = NO_INFECTION;
		myRXFrontInside[cfr].P = -1.0;
		myRXFrontInside[cfr].mydefgid = NO_GRAIN_ASSIGNED;
		myRXFrontInside[cfr].myrxgid = NO_GRAIN_ASSIGNED;
	}


	//allocate memory to store RecycledCells
	uint32_t ninitialRecyclingList = (uint32_t) ( mycaHdl->initialRelCellCaching * (double) myGeom.nregvol_rdtdnd * mycaHdl->maxfillperstep * FULLRECYCLING_ATT_FACTOR);
	QUICKASSERT ( ninitialRecyclingList < CA_ALLOCATION_MAXIMUM );
	myRecyclingListInside = NULL;
	myRecyclingListInside = new uint32_t[ninitialRecyclingList];
	QUICKASSERT( myRecyclingListInside != NULL );
	regMemGuard = regMemGuard + ( ninitialRecyclingList * sizeof(uint32_t) );


	ntotalRecyclingListInside = ninitialRecyclingList;
	nextSlotThatBecomesRecycledInside = NOTHING_TO_RECYCLE; //not assigned as there is no cell at all yet  that can be recycled 
	firstNotRecycledYetInside = NOTHING_TO_RECYCLE;

	//myRecyclingListInside needs no initialization because it is assured to have no gaps and only read [0;nextSlotThatBecomesRecycledInside)

	//allocate maxfillperstep the initial memory for cell that complete their infection in a timestep
	//this list serves as a guide to reduce the amount of ifs during finding recyclable sites during the infection phase as the RXFrontList becomes for sure fragmented and thus INACTIVE CELLs will likely be encountered
	uint32_t ninitialFullRXList = (uint32_t) ( mycaHdl->initialRelCellCaching * (double) myGeom.nregvol_rdtdnd * mycaHdl->maxfillperstep * FULLRXCACHING_ATT_FACTOR);
	QUICKASSERT ( ninitialFullRXList < CA_ALLOCATION_MAXIMUM );
	myFullRXListInside = NULL;
	myFullRXListInside = new uint32_t[ninitialFullRXList];
	QUICKASSERT( myFullRXListInside != NULL );
	regMemGuard = regMemGuard + ( ninitialFullRXList * sizeof(uint32_t) );

	ntotalFullRXListInside = ninitialFullRXList;
	nextSlotToFullRXInside = 0; //as above not assigned yet

	//no initialization as well assured gap-free read-only [0;nextSlotToFullRXInside)

	//cout << this->jobid << "\t\tMemGuard (Byte)= " << this->myMemGuard << endl;
}


void caregionMemHdl::ompshar_init_myRXFrontBorder( void )
{
	//MK::it can happen that a thread infects from its border into the border of an adjacent thread memory region,
	//thus memory in the myRXFrontBorder of thread B is required to store the cell, however that array might not have been allocated large enough
	//then thread A would be required to issue or request at least a reallocation such that thread B enlarges its myRXFrontBorder array
	//to avoid this we preallocate in advance already for all possible cells in the outer shell and as such the list of border cells does not require a reallocation
	//this comes at the cost of more main memory utilization but potentially better cache linearity..., however, for large SUs, i.e. where this may be of relevance, 
	//the interface area to volume is low consult PhD thesis by M. K\"uhbach for finding substantiated arguments as to why this allocation is not of a significant problem!

	uint32_t nRXFront = 2*( (myGeom.nreg_rd * myGeom.nreg_td) + (myGeom.nreg_rd * myGeom.nreg_nd) + (myGeom.nreg_td * myGeom.nreg_nd) ); //##MK::an exact calculation could save a few cells here...
	QUICKASSERT ( nRXFront < CA_ALLOCATION_MAXIMUM );
	myRXFrontBorder = NULL;
	myRXFrontBorder = new struct cell[nRXFront];
	QUICKASSERT( myRXFrontBorder != NULL );
	regMemGuard = regMemGuard + (nRXFront * sizeof(cell));

	ntotalRXFrontBorder = nRXFront;
	nextSlotNeverActiveRXFrontBorder = 0;
	nCurrentlyActiveBorder = 0;

	//initialize already associated cells
	for ( uint32_t cfr = 0; cfr < nRXFront; ++cfr ) {
		myRXFrontBorder[cfr].activity = INACTIVE;
		myRXFrontBorder[cfr].infector = 26;
		myRXFrontBorder[cfr].ix = -1;
		myRXFrontBorder[cfr].iy = -1;
		myRXFrontBorder[cfr].iz = -1;
		myRXFrontBorder[cfr].rxFrac = NO_INFECTION;
		myRXFrontBorder[cfr].P = -1.0;
		myRXFrontBorder[cfr].mydefgid = NO_GRAIN_ASSIGNED;
		myRXFrontBorder[cfr].myrxgid = NO_GRAIN_ASSIGNED;
	}


	//allocate memory to store RecycledCells, list does not grow!
	uint32_t nRecyclingList = (uint32_t) ( (double) nRXFront * mycaHdl->maxfillperstep * FULLRECYCLING_ATT_FACTOR);
	QUICKASSERT ( nRecyclingList < CA_ALLOCATION_MAXIMUM );
	myRecyclingListBorder = NULL;
	myRecyclingListBorder = new uint32_t[nRecyclingList];
	QUICKASSERT( myRecyclingListBorder != NULL );
	regMemGuard = regMemGuard + ( nRecyclingList * sizeof(uint32_t) );

	ntotalRecyclingListBorder = nRecyclingList;
	nextSlotThatBecomesRecycledBorder = NOTHING_TO_RECYCLE; //not assigned as there is no cell at all yet  that can be recycled 
	firstNotRecycledYetBorder = NOTHING_TO_RECYCLE;

	//myRecyclingListBorder needs no initialization because it is assured to have no gaps and only read [0;nextSlotThatBecomesRecycledBorder)

	//allocate maxfillperstep the initial memory for cell that complete their infection in a timestep
	//this list serves as a guidance to reduce the amount of ifs during the infection phase as the list is for sure fragmented and thus INACTIVE CELLs are  met
	uint32_t nFullRXList = (uint32_t) ( (double) nRXFront * mycaHdl->maxfillperstep * FULLRXCACHING_ATT_FACTOR);
	QUICKASSERT ( nFullRXList < CA_ALLOCATION_MAXIMUM );
	myFullRXListBorder = NULL;
	myFullRXListBorder = new uint32_t[nFullRXList];
	QUICKASSERT( myFullRXListBorder != NULL );
	regMemGuard = regMemGuard + ( nFullRXList * sizeof(uint32_t) );

	ntotalFullRXListBorder = nFullRXList;
	nextSlotToFullRXBorder = 0; //as above not assigned yet

	//no initialization as well  assured gap-free read-only [0;nextSlotToFullRXBorder)

	//cout << this->jobid << "\t\tMemGuard (Byte)= " << this->myMemGuard << endl;
}


void caregionMemHdl::ompshar_init_serialsectioning()
{
	//2d implicit array for all my ND perpendicular RDTD slices a bookkeep how much
	//specific grains intrude
	rdtd_serialsection = vector<uint32_t>( static_cast<size_t>(myGeom.nreg_nd), 0 );
}


uint32_t caregionMemHdl::omp_getNextFreeSlotBorderAnotherRegion ( uint32_t nbortid, bool how )
{
	//CALLED FROM WITHIN PARALLEL REGION BUT EXECUTION ORDER STRICTLY PREVENTS MULTIPLE EXECUTION
	//get a free place in the myRXFront List of the nbortids region!
	//utilize that this function issues an infection from the border of region from thread tid into the border of region from thread nbortid
	//MK::contrary to getNextFreeSlotInRXFrontInside it is a guaranteed always a large enough container so as well in nbortid to avoid necessity for reallocation
	//MUST NOT BE CALLED WHEN NOTHING WAS RECYCLED DURING A CALCGROWTHSTEP

	bool strategy = how;
	uint32_t place = INVALID_ADDRESS;
	caregionMemHdlP targetregion = mycaHdl->regions[nbortid]; //reflexive reference

	if ( strategy == CELLRECYCLE ) {
		if ( targetregion->nextSlotThatBecomesRecycledBorder < targetregion->firstNotRecycledYetBorder ) {
			place = targetregion->myRecyclingListBorder[targetregion->nextSlotThatBecomesRecycledBorder];
			targetregion->nextSlotThatBecomesRecycledBorder++;

			return place;
		}

		//hasnt returned yet, so recycling list should be utilized but was already exhausted so change the strategy to get memory
		strategy = CELLAPPEND;
	}

	//now strategy is/or was CELLAPPEND
	//if ( strategy == CELLAPPEND ) {
	if ( targetregion->nextSlotNeverActiveRXFrontBorder < targetregion->ntotalRXFrontBorder ) {
		place = targetregion->nextSlotNeverActiveRXFrontBorder;
		targetregion->nextSlotNeverActiveRXFrontBorder++;
	}

	//return an address, maybe one that is marked as invalid
	return place;
}


uint32_t caregionMemHdl::omp_getNextFreeSlotInRXFrontInside( bool how ) //MK::how == 0 == CELLRECYCLE ONLY CALL FOR FROM OUT OF FRONT
{
	//CALLED FROM WITHIN PARALLEL REGION, MUST NOT BE CALLED WHEN NOTHING WAS RECYCLED DURING A CALCGROWTHSTEP
	//MK::two extreme memory fetch strategies possible
	//1::always get more memory (almost no additional scans for places but increasingly costly cache hit during calcGrowthStep owing to accumulation of INACTIVE cells along the list of cells)
	//2::recycle as best as possible (reduced total amount of memory surplus less scans on average: very likely improved cache performance than strategy 1)
	//optimal would be to maintain a compact cell list (meaning no cell in the bucket is INACTIVE), however this will incurr overhead to ensure this or
	//or require to utilize a linked list structures which to iterate over --- however --- is orders of magnitude slower for millions of entries if naive or slower at least if smartly 
	//implemented as there will be pointer overhead

	//MK::the present implement selects most beneficial compromise strategy automatically
	bool strategy = how;
	uint32_t place = INVALID_ADDRESS;

	//if recycling is desired
	if ( strategy == CELLRECYCLE ) {
		if ( nextSlotThatBecomesRecycledInside < firstNotRecycledYetInside ) {
			place = myRecyclingListInside[nextSlotThatBecomesRecycledInside];
			nextSlotThatBecomesRecycledInside++;
			return place;
		}

		//hasnt returned yet, so recycling list should be utilized but was already exhausted so change the strategy to get memory
		strategy = CELLAPPEND;
	}

	//now strategy is/or was CELLAPPEND
	//if ( strategy == CELLAPPEND ) {
		//was the myRXFrontInside container precached still enough to allow the ACTIVATION OF AN INACTIVE cell?
	if ( nextSlotNeverActiveRXFrontInside < ntotalRXFrontInside ) {
		place = nextSlotNeverActiveRXFrontInside;
		nextSlotNeverActiveRXFrontInside++;
		return place;
	}

	//obviously	, recycling list already exhausted, precached places all occupied, damn it
	//no cache left, and ohh, yes mind that accesses such as myRXFrontInside[nextSlotNeverActiveRXFrontInside] for >= ntotalRXFrontInside is asking for trouble
	//so we need to get new memory!
	size_t oldSize = ntotalRXFrontInside;
	size_t newSize = oldSize + (mycaHdl->transientRelCellRecaching * oldSize);
	QUICKASSERT( newSize < CA_ALLOCATION_MAXIMUM );

//#ifdef REPORTSTYLE_DEVELOPER
//	cout << "\t\tALLOCATION for new_myRXFront;oldSize;newSize\t" << ((int) oldSize) << "\t" << ((int) newSize) << endl;
//#endif

	cellP new_myRXFrontInside = NULL;
	new_myRXFrontInside = new struct cell[newSize];
	QUICKASSERT( new_myRXFrontInside != NULL );
	regMemGuard = regMemGuard + (newSize * sizeof(cell));

	//potentially unsafe to use memcpy because of overlapping issues memcpy( new_myRXFront, myRXFront, oldSize * sizeof(cell) ); 
	//MK::all indices in FullRX and Recycling list stay valid as they are referencing relative to first array position
	for ( uint32_t i = 0; i < oldSize; ++i ) {
		new_myRXFrontInside[i].activity = myRXFrontInside[i].activity;
		new_myRXFrontInside[i].infector = myRXFrontInside[i].infector;
		new_myRXFrontInside[i].ix = myRXFrontInside[i].ix;
		new_myRXFrontInside[i].iy = myRXFrontInside[i].iy;
		new_myRXFrontInside[i].iz = myRXFrontInside[i].iz;
		new_myRXFrontInside[i].rxFrac = myRXFrontInside[i].rxFrac;
		new_myRXFrontInside[i].P = myRXFrontInside[i].P;
		new_myRXFrontInside[i].mydefgid = myRXFrontInside[i].mydefgid;
		new_myRXFrontInside[i].myrxgid = myRXFrontInside[i].myrxgid;
	}

	delete [] myRXFrontInside;
	regMemGuard = regMemGuard - (oldSize * sizeof(cell));
	myRXFrontInside = new_myRXFrontInside;
	ntotalRXFrontInside = newSize;

	//initialize this new memory snippet, access with nc = oldSize is to 0-th element of appended section
	for ( uint32_t nc = oldSize; nc < ntotalRXFrontInside; ++nc) {
		myRXFrontInside[nc].activity = INACTIVE;
		myRXFrontInside[nc].infector = 26;
		myRXFrontInside[nc].ix = -1;
		myRXFrontInside[nc].iy = -1;
		myRXFrontInside[nc].iz = -1;
		myRXFrontInside[nc].rxFrac = NO_INFECTION;
		myRXFrontInside[nc].P = -1.0;
		myRXFrontInside[nc].mydefgid = NO_GRAIN_ASSIGNED;
		myRXFrontInside[nc].myrxgid = NO_GRAIN_ASSIGNED;
	}

	//now we have new space, so append
	place = nextSlotNeverActiveRXFrontInside;
	nextSlotNeverActiveRXFrontInside++;
	return place;
}


uint32_t caregionMemHdl::omp_getNextFreeSlotInRXFrontBorder( bool how )
{
	//CALLED FROM WITHIN PARALLEL REGION, MUST NOT BE CALLED WHEN NOTHING WAS RECYCLED DURING A CALCGROWTHSTEP
	//MK::contrary to getNextFreeSlotInRXFrontInside we are now guaranteed to have a large enough container myRXFrontBorder to avoid a potential reallocation
	//algorithm similar to omp_getNextFreeSlotInRXFrontInside, see further comments there...
	bool strategy = how;
	uint32_t place = INVALID_ADDRESS;

	if ( strategy == CELLRECYCLE ) {
		if ( nextSlotThatBecomesRecycledBorder < firstNotRecycledYetBorder ) {
			place = myRecyclingListBorder[nextSlotThatBecomesRecycledBorder];
			nextSlotThatBecomesRecycledBorder++;
			return place;
		}

		//hasnt returned yet, so recycling list should be utilized but was already exhausted so change the strategy to get memory
		strategy = CELLAPPEND;
	}

	//now strategy is/or was CELLAPPEND
	if ( nextSlotNeverActiveRXFrontBorder < ntotalRXFrontBorder ) {
		place = nextSlotNeverActiveRXFrontBorder;
		nextSlotNeverActiveRXFrontBorder++;
	}

	//return an address, maybe one that is marked as invalid
	return place;
}


void caregionMemHdl::ompshar_defragmentation( void ) 
{
	//CALLED FROM WITHIN PARALLEL REGION!
	//MK::usually the function is kicked in shortly after the peak of the Sv(X) function
	//because by then new infection requests decrease strongly while at the same time ACTIVE cells progressively become INACTIVE
	//hence we can pack the working arrays myRXFrontInside and Border more compactly in contiguous memory by leaving as few INACTIVE cells as possible
	double tprof = omp_get_wtime();

	//cout << "Defragmenting..." << endl;

	//first the inside list
	uint32_t firstinactive = 0;
	while ( myRXFrontInside[firstinactive].activity == ACTIVE ) {
		firstinactive++;
	}
	//extend initial contiguous block by copying over all active cells
	uint32_t contiguous = firstinactive;
	//"popping out inactive cells",i.e. contracting the active cells
	//progressively less write cache efficient as contiguous end falls behind actual position c, however the reads in c remain cache efficient
	for ( uint32_t c = (firstinactive + 1); c < nextSlotNeverActiveRXFrontInside; c++ ) { //the +1 offset avoids selfassignment
		if ( myRXFrontInside[c].activity == ACTIVE ) {
			myRXFrontInside[contiguous].activity = myRXFrontInside[c].activity;
			myRXFrontInside[contiguous].infector = myRXFrontInside[c].infector;
			myRXFrontInside[contiguous].ix = myRXFrontInside[c].ix;
			myRXFrontInside[contiguous].iy = myRXFrontInside[c].iy;
			myRXFrontInside[contiguous].iz = myRXFrontInside[c].iz;
			myRXFrontInside[contiguous].rxFrac = myRXFrontInside[c].rxFrac;
			myRXFrontInside[contiguous].P = myRXFrontInside[c].P;
			myRXFrontInside[contiguous].mydefgid = myRXFrontInside[c].mydefgid;
			myRXFrontInside[contiguous].myrxgid = myRXFrontInside[c].myrxgid;
			contiguous++;
		}
	}
	nextSlotNeverActiveRXFrontInside = contiguous;

	//now the border
	firstinactive = 0;
	while ( myRXFrontBorder[firstinactive].activity == ACTIVE ) {
		firstinactive++;
	}	
	contiguous = firstinactive;
	for ( uint32_t c = (firstinactive + 1); c < nextSlotNeverActiveRXFrontBorder; c++ ) {
		if ( myRXFrontBorder[c].activity == ACTIVE ) {
			myRXFrontBorder[contiguous].activity = myRXFrontBorder[c].activity;
			myRXFrontBorder[contiguous].infector = myRXFrontBorder[c].infector;
			myRXFrontBorder[contiguous].ix = myRXFrontBorder[c].ix;
			myRXFrontBorder[contiguous].iy = myRXFrontBorder[c].iy;
			myRXFrontBorder[contiguous].iz = myRXFrontBorder[c].iz;
			myRXFrontBorder[contiguous].rxFrac = myRXFrontBorder[c].rxFrac;
			myRXFrontBorder[contiguous].P = myRXFrontBorder[c].P;
			myRXFrontBorder[contiguous].mydefgid = myRXFrontBorder[c].mydefgid;
			myRXFrontBorder[contiguous].myrxgid = myRXFrontBorder[c].myrxgid;
			contiguous++;
		}
	}
	nextSlotNeverActiveRXFrontBorder = contiguous;

	dtDefragment = omp_get_wtime() - tprof;
/*	//MK::may only be called after resetting of recycling and FullRX list or prior to calcGrowthStep!
	uint32_t which = 0;
	uint32_t moved = nextSlotNeverActiveRXFrontInside - 1; //right-most/last valid entry accessible in myRXFront [0,nextSlotNeverActiveRXFront)
	//therewith automatically no defragmentation of a single entry or empty list
	//the idea is the following copy state of cells from high-cnt to low-cnt index
	while ( which < moved ) {
		if ( myRXFrontInside[which].activity == INACTIVE ) {
			myRXFrontInside[which].activity = myRXFrontInside[moved].activity;
			myRXFrontInside[which].infector = myRXFrontInside[moved].infector;
			myRXFrontInside[which].ix = myRXFrontInside[moved].ix;
			myRXFrontInside[which].iy = myRXFrontInside[moved].iy;
			myRXFrontInside[which].iz = myRXFrontInside[moved].iz;
			myRXFrontInside[which].rxFrac = myRXFrontInside[moved].rxFrac;
			myRXFrontInside[which].P = myRXFrontInside[moved].P;
			myRXFrontInside[which].mydefgid = myRXFrontInside[moved].mydefgid;
			myRXFrontInside[which].myrxgid = myRXFrontInside[moved].myrxgid;

			//moved might reference inactive cell, waste operation, but as we are filling up myRXFront active from left to right most entries are ACTIVE
			moved = moved - 1;

			//##MK::however the algorithm is still provoking almost ever a cache miss for the access to the moved-th element
			//therefore work through moved' = moved - L1_CACHE_LENGTH and increment moved'++ and a toggle++ if toggle == L1_CACHE_LENGTH...
			//thereby the first access to moved' after a toggle operation is more likely to hit only a cold-miss but then advancing linearly 
			//through one cache block while the which-th element is in another cache block, clearly
			//MK::the details are expected strongly dependent on the set associativity and other cache characteristics
		}
		//points to next potential INACTIVE cell or cache-local ACTIVE cell, anyways
		which = which + 1;
	}

	//assured deterministical defragmentation process is completed
	nextSlotNeverActiveRXFrontInside = moved + 1;


	//defragment now the border!
	uint32_t whichb = 0;
	uint32_t movedb = nextSlotNeverActiveRXFrontBorder - 1;

	while ( whichb < movedb ) {
		if ( myRXFrontBorder[whichb].activity == INACTIVE ) {
			myRXFrontBorder[whichb].activity = myRXFrontBorder[movedb].activity;
			myRXFrontBorder[whichb].infector = myRXFrontBorder[movedb].infector;
			myRXFrontBorder[whichb].ix = myRXFrontBorder[movedb].ix;
			myRXFrontBorder[whichb].iy = myRXFrontBorder[movedb].iy;
			myRXFrontBorder[whichb].iz = myRXFrontBorder[movedb].iz;
			myRXFrontBorder[whichb].rxFrac = myRXFrontBorder[movedb].rxFrac;
			myRXFrontBorder[whichb].P = myRXFrontBorder[movedb].P;
			myRXFrontBorder[whichb].mydefgid = myRXFrontBorder[movedb].mydefgid;
			myRXFrontBorder[whichb].myrxgid = myRXFrontBorder[movedb].myrxgid;

			movedb = movedb - 1;
		}
		whichb = whichb + 1;
	}
	//assured deterministical defragmentation process is completed
	nextSlotNeverActiveRXFrontBorder = movedb + 1;*/
}


void caregionMemHdl::append_to_fullrx_inside_list( uint32_t cid )
{
	if ( nextSlotToFullRXInside < ntotalFullRXListInside ) {		//FullRX list still large enough, append
		myFullRXListInside[nextSlotToFullRXInside] = cid;		//indexing relative to first element in array with displacement sizeof(int)
		nextSlotToFullRXInside++;
		return;
	}

	//FullRX list is too small but it is necessary to store all cells that updateFullRX should work on, so increase size of FullRX list
	size_t oldSize = ntotalFullRXListInside;
	size_t newSize = oldSize + (mycaHdl->transientRelCellRecaching * oldSize);
	QUICKASSERT ( newSize < CA_ALLOCATION_MAXIMUM );

//#ifdef REPORTSTYLE_DEVELOPER
//	cout << "\t\tALLOCATION for fullrx list, oldSize, newSize\t" << ((int) oldSize) << "\t" << ((int) newSize) << endl;
//#endif

	uint32_t* new_myFullRXListInside = NULL;
	new_myFullRXListInside = new uint32_t[newSize];
	QUICKASSERT( new_myFullRXListInside != NULL );
	regMemGuard = regMemGuard + (newSize * sizeof(uint32_t));

	//##MK::potentially unsafe to use memcpy memcpy( new_myFullRXList, myFullRXList, oldSize * sizeof(int) );
	for ( uint32_t i = 0; i < oldSize; ++i ) {
		new_myFullRXListInside[i] = myFullRXListInside[i];
	}

	delete [] myFullRXListInside;
	regMemGuard = regMemGuard - ( oldSize * sizeof(uint32_t) );
	myFullRXListInside = new_myFullRXListInside;
	ntotalFullRXListInside = newSize;

	QUICKASSERT ( nextSlotToFullRXInside < ntotalFullRXListInside );

	//more space again, so append
	myFullRXListInside[nextSlotToFullRXInside] = cid;
	nextSlotToFullRXInside++;
}


void caregionMemHdl::append_to_recycling_inside_list( uint32_t cid )
{
	//firstNotRecycledYetPoints to position in myRecyclingList in [0, ntotalRecyclingListInside) no cell reference to myRXFrontInside was stored yet
	if ( firstNotRecycledYetInside < ntotalRecyclingListInside ) {
		myRecyclingListInside[firstNotRecycledYetInside] = cid;
		firstNotRecycledYetInside++;
	}
	//else -> currently I do not recycle more than ntotalRecyclingList elements, also this list is already as large as 0.15*ntotalcells, user can optimize this structure when knowing his microstructural path function Sv(X) in more details to save additional memory therefore not further memory utilized
}



void caregionMemHdl::ompshar_sim_myCA_calcGrowthStepInside( void )
{
	double tprof = omp_get_wtime();

	//##consider to mount them fixed in the classHdl
	double shapefactors[27] = { FGEODIAG, FGEOEDGE, FGEODIAG, FGEOEDGE, FGEOFACE, FGEOEDGE, FGEODIAG, FGEOEDGE, FGEODIAG, 
						FGEOEDGE, FGEOFACE, FGEOEDGE, FGEOFACE,      FGEOFACE, FGEOEDGE, FGEOFACE, FGEOEDGE,	
						FGEODIAG, FGEOEDGE, FGEODIAG, FGEOEDGE, FGEOFACE, FGEOEDGE, FGEODIAG, FGEOEDGE, FGEODIAG, (1.0) }; //layers xz stacked in positive y direction...

	uint32_t IdentifiedAsActive = 0;

	//MK::reset collector array for already transformed cells
	nextSlotToFullRXInside = 0;

	//MK::reset RecyclingCounter as int-values from previous step are now corrupted, up to now nothing has been recycled
	nextSlotThatBecomesRecycledInside = 0;
	firstNotRecycledYetInside = 0;

	//scan the whole myRXFront on the interval [0;nextSlotNeverActiveRXFront) to 
	//	a) to find ACTIVE cells that require simulation of boundary migration increments during annealing real time step dt
	//	b) to find INACTIVE cells that can be recycled
//##DEBUG#pragma omp critical { cout << "\t\tEntering calcGrowthStepInside;step;nextSlotNeverActiveRXFront;myca->dt;myca->cellsize\t" << mycaHdl->step << "\t" << nextSlotNeverActiveRXFrontInside << "\t" << mycaHdl->dt << "\t" << mycaHdl->_cellsize << endl; }

	//MK::DEVELOPER NOTE:when nucleation should be time-dependent it is vital to have the nuclei insert here AT THE LATEST, otherwise nextSlotNeverActiveRX will become updated and as such corrupted...
	double dVrx = 0.0;
	double v = 0.0;

	for ( uint32_t c = 0; c < nextSlotNeverActiveRXFrontInside; c++ ) {

		//MK::ACTIVE is the most likely case in a well defragmented list
		if ( myRXFrontInside[c].activity == ACTIVE ) {
			IdentifiedAsActive++;

			//read cached mobilityweight from local cell container that is assumed fixed since the infection
			//v = mp Turnbull linearized rate equation model
			double m = mycaHdl->get_currentintrinsicmobility ( myRXFrontInside[c].P );
			double p = (mycaHdl->Gbhalfsq*(mycaHdl->get_rho(myRXFrontInside[c].mydefgid) - RHO_RECRYSTALLIZED_MATERIAL)) - mycaHdl->get_zener();
			if ( p > DOUBLE_ACCURACY ) { //most likely there is driving force

				//compute migration velocity
				v = m * p;

				//scale velocity according to infection direction and stencil kernel function to assure self-similar growth
				v *= shapefactors[myRXFrontInside[c].infector];

				//eliminate unnecessary backreferences to mycaHdl
				myRXFrontInside[c].rxFrac = myRXFrontInside[c].rxFrac + (v * mycaHdl->dt * mycaHdl->_cellsize); //v*(max*cellsize/vmax)*v

//##DEBUG#pragma omp critical { 	cout << "\t\t\t\tstep->calc_ACTIVE;step;c;v;mob;[c]->P;[c]->rxFrac\t" << mycaHdl->step << "\t" << c << "\t" << v << "\t" << setprecision(14) <<  myRXFrontInside[c].P << "\t" << setprecision(14) << myRXFrontInside[c].rxFrac << endl; }

				dVrx = dVrx + (v * mycaHdl->dt * mycaHdl->_cellsize);

				if ( myRXFrontInside[c].rxFrac >= 1.0 ) { //does this fill the cell completely?
					//cout << "APPENDTOFULLRX\t\t" << c << ";" << myRXFrontInside[c].rxFrac << endl;

//##DEBUG#pragma omp critical { cout << "\t\tstep->calcACTIVE appending to recrystallized list" << c << "-->rxFrac=" << myRXFrontInside[c].rxFrac << endl; }

					this->append_to_fullrx_inside_list ( c );
				}
				continue;
			}
			// p almost zero
			//cout << "WARNING::Migration speed is almost zero\n" << endl;
			continue;
		} 
		//else {
			//else cell seems inactive, okay remember index c in RecyclingList

			this->append_to_recycling_inside_list ( c );

//##DEBUG#pragma omp critical { cout << "\t\t\t\tstep->calc_INACTIVE cell appended;c\t" << c << endl; }


		//}

	} //analyze all cells in [0, nextSlotNeverActiveRXFrontInside)


	reg_dXstepInside = dVrx;
	reg_SvInside = IdentifiedAsActive;
	reg_nCurrActiveInside = IdentifiedAsActive;

//##DEBUG#pragma omp critical{ cout << "\t\tLeaving calcGrowthStepInside;step;dXstep;Sv;Active\t" << mycaHdl->step << "\t" << reg_dXstepInside << "\t" << reg_SvInside << "\t" << reg_nCurrActiveInside << endl; }
	dtCalcGrowthInside = omp_get_wtime() - tprof;
}


void caregionMemHdl::append_to_fullrx_border_list( uint32_t cid )
{
	if ( nextSlotToFullRXBorder < ntotalFullRXListBorder ) {		//FullRX list still large enough, append
		myFullRXListBorder[nextSlotToFullRXBorder] = cid;		//indexing relative to first element in array with displacement sizeof(int)
		nextSlotToFullRXBorder++;
		return;
	}

	//FullRX list is too small but it is necessary to store all cells that updateFullRX should work on, so increase size of FullRX list
	size_t oldSize = ntotalFullRXListBorder;
	size_t newSize = oldSize + (mycaHdl->transientRelCellRecaching * oldSize);
	QUICKASSERT ( newSize < CA_ALLOCATION_MAXIMUM );

//#ifdef REPORTSTYLE_DEVELOPER
	cout << "\t\tALLOCATION for fullrx list, oldSize, newSize\t" << ((int) oldSize) << "\t" << ((int) newSize) << endl;
//#endif

	uint32_t* new_myFullRXListBorder = NULL;
	new_myFullRXListBorder = new uint32_t[newSize];
	QUICKASSERT( new_myFullRXListBorder != NULL );
	regMemGuard = regMemGuard + (newSize * sizeof(uint32_t));

	//##MK::potentially unsafe to use memcpy memcpy( new_myFullRXList, myFullRXList, oldSize * sizeof(int) );
	for ( uint32_t i = 0; i < oldSize; ++i ) {
		new_myFullRXListBorder[i] = myFullRXListBorder[i];
	}

	delete [] myFullRXListBorder;
	regMemGuard = regMemGuard - ( oldSize * sizeof(uint32_t) );
	myFullRXListBorder = new_myFullRXListBorder;
	ntotalFullRXListBorder = newSize;

	QUICKASSERT ( nextSlotToFullRXBorder < ntotalFullRXListBorder );

	//more space again, so append
	myFullRXListBorder[nextSlotToFullRXBorder] = cid;
	nextSlotToFullRXBorder++;
}


void caregionMemHdl::append_to_recycling_border_list( uint32_t cid )
{
	if ( firstNotRecycledYetBorder < ntotalRecyclingListBorder ) {
		myRecyclingListBorder[firstNotRecycledYetBorder] = cid;
		firstNotRecycledYetBorder++;
	}
	//else -> currently I do not recycle more than ntotalRecyclingList elements, also this list is already as large as 0.15*ntotalcells, user can optimize this structure when knowing his microstructural path function Sv(X) in more details to save additional memory therefore not further memory utilized
}


void caregionMemHdl::ompshar_sim_myCA_calcGrowthStepBorder( void )
{
	double tprof = omp_get_wtime();

	//##consider to mount them fixed in the classHdl
	double shapefactors[27] = { FGEODIAG, FGEOEDGE, FGEODIAG, FGEOEDGE, FGEOFACE, FGEOEDGE, FGEODIAG, FGEOEDGE, FGEODIAG, 
						FGEOEDGE, FGEOFACE, FGEOEDGE, FGEOFACE,      FGEOFACE, FGEOEDGE, FGEOFACE, FGEOEDGE,	
						FGEODIAG, FGEOEDGE, FGEODIAG, FGEOEDGE, FGEOFACE, FGEOEDGE, FGEODIAG, FGEOEDGE, FGEODIAG, (1.0) }; //layers xz stacked in positive y direction...

	uint32_t IdentifiedAsActive = 0;

	//MK::reset collector array for already transformed cells
	nextSlotToFullRXBorder = 0;

	//MK::reset RecyclingCounter as int-values from previous step are now corrupted, up to now nothing has been recycled
	nextSlotThatBecomesRecycledBorder = 0;
	firstNotRecycledYetBorder = 0;


	//scan the whole myRXFront on the interval [0;nextSlotNeverActiveRXFront) to 
	//	a) to find ACTIVE cells for which boundary migration is simulated during dt
	//	b) to find INACTIVE cells
	//all comments from caHdl::ompshar_sim_myCA_calcGrowthStepInside( void ) apply as well

	double dVrx = 0.0;
	double v = 0.0;
	for ( uint32_t c = 0; c < nextSlotNeverActiveRXFrontBorder; c++ ) {
		if ( myRXFrontBorder[c].activity == ACTIVE ) {
			IdentifiedAsActive++;

			//v = mp Turnbull linearized rate equation model
			double m = mycaHdl->get_currentintrinsicmobility ( myRXFrontBorder[c].P );
			double p = ( mycaHdl->Gbhalfsq * (mycaHdl->get_rho(myRXFrontBorder[c].mydefgid) - RHO_RECRYSTALLIZED_MATERIAL) ) - mycaHdl->get_zener();
			if ( p >  DOUBLE_ACCURACY ) {

				//compute migration velocity
				v = m * p;

				//scale velocity
				v *= shapefactors[myRXFrontBorder[c].infector];

				myRXFrontBorder[c].rxFrac = myRXFrontBorder[c].rxFrac + (v * mycaHdl->dt * mycaHdl->_cellsize);

				dVrx = dVrx + (v * mycaHdl->dt * mycaHdl->_cellsize);

				//has the cell volume recrystallized almost entirely?
				if ( myRXFrontBorder[c].rxFrac >= 1.0 ) {
					this->append_to_fullrx_border_list ( c );
				}
				continue;
			}
			//p almost zero
			//cout << "WARNING::Migration speed is almost zero\n" << endl;
			continue;
		}
		//else cell == INACTIVE
		this->append_to_recycling_border_list ( c );

	} //analyze all cells in [0, nextSlotNeverActiveRXFrontBorder)

	reg_dXstepBorder = dVrx;
	reg_SvBorder = IdentifiedAsActive;
	reg_nCurrActiveBorder = IdentifiedAsActive;

//##DEBUG#pragma omp critical{ cout << "\t\tLeaving calcGrowthStepBorder;step;dXstep;Sv;Active\t" << mycaHdl->step << "\t" << reg_dXstepBorder << "\t" << reg_SvBorder << "\t" << reg_nCurrActiveBorder << endl; }
	dtCalcGrowthBorder = omp_get_wtime() - tprof;
}


void caregionMemHdl::omp_infect_OutOfRXFrontInside ( uint32_t rxgpoolid, uint32_t rxfrontid, bool intoborder, short dx, short dy, short dz, unsigned char direction, double carryOver )
{
	//CALLED FROM WITHIN PARALLEL REGION BUT OPERATING ON caregionMemHdl class object's local cellgrid
	
	short x = myRXFrontInside[rxfrontid].ix;
	short y = myRXFrontInside[rxfrontid].iy;
	short z = myRXFrontInside[rxfrontid].iz;

	//add relative coordinate increment to identify position of infection target in the CA, global coordinates!
	x += dx;
	y += dy;
	z += dz;
	//MK::no bounds check necessary as either the infection is directed into the bulk of the domain (intoborder == false) or its border (intoborder == true)
	//hence x,y,z are positive and thus the promotion from short to uint32_t is not a problem

	//convert ix,iy,iz in implicit 3D coordinate to access mycellgrid
	//MK::WORKS ONLY WITH SHORT IS POSITIVE!
	uint32_t cxyz = (x - myGeom.nreg_rdmin) + ((y - myGeom.nreg_tdmin) * myGeom.nreg_rd) + ((z - myGeom.nreg_ndmin) * myGeom.nregarea_rdtd);

	uint32_t cxyz_deformedstate = mycellgrid[cxyz];

	if ( cxyz_deformedstate >= reg_nmydefgpool || cxyz_deformedstate == CURRENTLY_INFECTED ) {
//##DEBUG#pragma omp critical { cout << "\t\t\tREJECT infection cxyz;x;y;z;dx;dy;dz;deformedstate;reg_nmydefgpool\t" << cxyz << "\t" << x << ";" << y << ";" << z << "\t" << dx << ";" << dy << ";" << dz << "\t" << cxyz_deformedstate << "\t" << reg_nmydefgpool << endl; }
		return;
	}

	//obviously cell is free, hehe, go for it!
	mycellgrid[cxyz] = CURRENTLY_INFECTED;

	//MK::now the infection is directed either into the bulk or the border but the infection remains in threadlocal memory
	uint32_t freeplace = INVALID_ADDRESS;
	if ( intoborder == false ) { //most likely case is to infect into the bulk because domains usually cuboids with millions of cells

		freeplace = this->omp_getNextFreeSlotInRXFrontInside ( CELLRECYCLE ); //MK::strategy is RECYCLING because OutofRXFront infections are only called from updateFullRX
//##DEBUG#pragma omp critical { 	cout << "\t\t\tACCEPT-INSIDE at freeplace\t\tcxyz;x;y;z;dx;dy;dz\t" << freeplace << "\t\t" << cxyz << "\t" << x << ";" << y << ";" << z << "\t" << dx << ";" << dy << ";" << dz << endl; }

		myRXFrontInside[freeplace].activity = ACTIVE;
		myRXFrontInside[freeplace].infector = direction;
		myRXFrontInside[freeplace].ix = x;
		myRXFrontInside[freeplace].iy = y;
		myRXFrontInside[freeplace].iz = z;
		myRXFrontInside[freeplace].rxFrac = carryOver;


		//first assume most likely case - infection continues consuming the same deformed grain id, then P is the same
		myRXFrontInside[freeplace].P = myRXFrontInside[rxfrontid].P;

		//compromise growth of the grain from fully infected deformed regionn mydefgid into a different deformed grain
		if ( myRXFrontInside[rxfrontid].mydefgid != cxyz_deformedstate ) {
			myRXFrontInside[freeplace].P = reg_calc_mobilityweight ( rxgpoolid, cxyz_deformedstate );
		}

		//myMobilityWeightMax is only relevant for the global adaptive dynamic time stepping
		if ( myRXFrontInside[freeplace].P >= reg_myMobilityWeightMax ) {
			reg_myMobilityWeightMax = myRXFrontInside[freeplace].P;
		}

		//cell carries consumed orientation along
		myRXFrontInside[freeplace].myrxgid = rxgpoolid;
		myRXFrontInside[freeplace].mydefgid = cxyz_deformedstate;

		//MK::assessor function get_NextFreeSlotInRXFrontAppendOrRecycle() ASSURES freeplace < nextSlotNeverActiveRXFront <= ntotalRXFront!
		return;
	}

	//not returned yet, okay infection aims for the myRXFrontBorder
	freeplace = this->omp_getNextFreeSlotInRXFrontBorder ( CELLRECYCLE ); //MK::strategy is RECYCLING because OutofRXFront infections are only called from updateFullRX
//##DEBUG#pragma omp critical { cout << "\t\t\tACCEPT-BORDER at freeplace\t\tcxyz;x;y;z;dx;dy;dz\t" << freeplace << "\t\t" << cxyz << "\t" << x << ";" << y << ";" << z << "\t" << dx << ";" << dy << ";" << dz << endl; }

	myRXFrontBorder[freeplace].activity = ACTIVE;
	myRXFrontBorder[freeplace].infector = direction;
	myRXFrontBorder[freeplace].ix = x;
	myRXFrontBorder[freeplace].iy = y;
	myRXFrontBorder[freeplace].iz = z;
	myRXFrontBorder[freeplace].rxFrac = carryOver;

	//first assume most likely case - infection continues consuming the same deformed grain id, then P is the same
	myRXFrontBorder[freeplace].P = myRXFrontInside[rxfrontid].P;

	//compromise growth of the grain from fully infected deformed region mydefgid into a different deformed grain
	if ( myRXFrontInside[rxfrontid].mydefgid != cxyz_deformedstate ) {
		myRXFrontBorder[freeplace].P = reg_calc_mobilityweight ( rxgpoolid, cxyz_deformedstate );
	}

	//myMobilityWeightMax is only relevant for the global adaptive dynamic time stepping
	if ( myRXFrontBorder[freeplace].P >= reg_myMobilityWeightMax ) {
		reg_myMobilityWeightMax = myRXFrontBorder[freeplace].P;
	}

	//cell carries consumed orientation along
	myRXFrontBorder[freeplace].myrxgid = rxgpoolid;
	myRXFrontBorder[freeplace].mydefgid = cxyz_deformedstate;
}


void caregionMemHdl::omp_infect_OutOfRXFrontBorder ( uint32_t rxgpoolid, uint32_t rxfrontid, short dx, short dy, short dz, unsigned char direction, double carryOver )
{
	//CALLED FROM WITHIN PARALLEL REGION
	//Three possible locations where the target cell is located:: inside, border, or in another region on top of that periodic bound checks necessary!

	short x = myRXFrontBorder[rxfrontid].ix;
	short y = myRXFrontBorder[rxfrontid].iy;
	short z = myRXFrontBorder[rxfrontid].iz;

	//add relative coordinate increment to identify position of infection target in the CA, global coordinates!
	x += dx;
	y += dy;
	z += dz;

	//apply periodic boundary conditions at solitary unit global scale
	if ( x < 0 )
		x += myGeom.nedge_global_rd;
	if ( x >= myGeom.nedge_global_rd )
		x -= myGeom.nedge_global_rd;
	if ( y < 0 )
		y += myGeom.nedge_global_td;
	if ( y >= myGeom.nedge_global_td )
		y -= myGeom.nedge_global_td;
	if ( z < 0 ) {
		return; //MK::modification for in-situ experiments of Martin Diehl that should not have periodicity in the +-ND direction
		//z += myGeom.nedge_global_nd;
	}
	if ( z >= myGeom.nedge_global_nd ) {
		return; //MK::modification for in-situ experiments of Martin Diehl that should not have periodicity in the +-ND direction
		//z -= myGeom.nedge_global_nd;
	}

//##DEBUG#pragma omp critical { cout << x << "\t" << y << "\t" << z << "\t" << rxgpoolid << endl; }

	//MK::there are several techniques to identify in which region the infection is, I have chosen for a
	//compromise between speed and clarity instead of making even more complicated branching conditions
	//the problem to solve is: how to get the threadid of the memregion in which the global coordinated x,y,z is located?
	//the strategy: inspect first the most likely, then the less likely and last the least like, i.e. the corners and edges of the region

	//COMPLETELY INSIDE for a cell on the inside of a face of a region this occurs in 9 directions
	if ( x > myGeom.nreg_rdmin && x < myGeom.nreg_rdmax && y > myGeom.nreg_tdmin && y < myGeom.nreg_tdmax && z > myGeom.nreg_ndmin && z < myGeom.nreg_ndmax ) {
		//infection directs into this this region, hence requires the activation of a cell from the myRXFrontInside list

		//convert x,y,z in implicit local 3D coordinates to access mycellgrid, works only because the short's are all positive...
		uint32_t cxyz = (x - myGeom.nreg_rdmin) + ((y - myGeom.nreg_tdmin) * myGeom.nreg_rd) + ((z - myGeom.nreg_ndmin) * myGeom.nregarea_rdtd);

		uint32_t cxyz_deformedstate = this->mycellgrid[cxyz];

//##DEBUG#pragma omp critical { cout << "\t\t\t" << omp_get_thread_num() << ";" << "SAMEREGION-INSIDE x,y,z,dx,dy,dz,cxyz,cxyz_deformedstate" << x << ";" << y << ";" << z << ";" << cxyz << ";" << cxyz_deformedstate << endl; }

		if ( cxyz_deformedstate >= reg_nmydefgpool || cxyz_deformedstate == CURRENTLY_INFECTED ) { //get out here, cell is already infected
//##DEBUG#pragma omp critical{ cout << "\t\t\t" << omp_get_thread_num() << ";" << "SAMEREGION-INSIDE returned" << endl; }
			return;
		}

		this->mycellgrid[cxyz] = CURRENTLY_INFECTED; 
		//MK::thread-safe because non-overlapping regions plus sequential processing of border cells, for
		//##MK::these cases in an even more elaborated preprocessing the sequential OpenMP overhead could be reduced further...

		uint32_t freeplace = this->omp_getNextFreeSlotInRXFrontInside ( CELLRECYCLE ); //MK::strategy is RECYCLING because OutofRXFront infections are only called from updateFullRX

		myRXFrontInside[freeplace].activity = ACTIVE;
		myRXFrontInside[freeplace].infector = direction;
		myRXFrontInside[freeplace].ix = x;
		myRXFrontInside[freeplace].iy = y;
		myRXFrontInside[freeplace].iz = z;
		myRXFrontInside[freeplace].rxFrac = carryOver;

		myRXFrontInside[freeplace].P = myRXFrontBorder[rxfrontid].P;
		if ( myRXFrontBorder[rxfrontid].mydefgid != cxyz_deformedstate ) { //recalculate mobility if encountering new orientation
			myRXFrontInside[freeplace].P = reg_calc_mobilityweight ( rxgpoolid, cxyz_deformedstate );
		}
		if ( myRXFrontInside[freeplace].P >= reg_myMobilityWeightMax ) { //dynamic timestepping
			reg_myMobilityWeightMax = myRXFrontInside[freeplace].P;
		}

		myRXFrontInside[freeplace].myrxgid = rxgpoolid;
		myRXFrontInside[freeplace].mydefgid = cxyz_deformedstate;
//##DEBUG#pragma omp critical{ cout << "\t\t\t" << omp_get_thread_num() << ";" << "SAMEREGION-BORDER success" << endl; }
		return;
	} 

	//infection directed NOT INSIDE SO NOW NECESSARY TO CHECK AGAINST LIMITS OF ALL POSSIBLE NEIGHBORS, still however it could be an infection that stays within the border layer
	//THIS REGION encloses all infections that direct not out of the shell, i.e. 8 of 26 directions
	//albeit there is not a simple test, like ^sum di^2 == 2 to check for this as the cells can be also located on the edge of the region
	//but lets assume for a moment that the infection is still on this border, it cannot be inside the region as that case was already checked for
	//so check this and confirm that infection is not in another region
	uint32_t nbortid = UNKNOWN_NEIGHBOR;
	//uint32_t nborpos = UNKNOWN_POSITION;
	//bool logicalfail = false;

	for ( uint32_t nb = 0; nb < (MYSELF + NUMBER_OF_NEIGHBORS); nb++ ) {
		//MK::carefully reconsider to change these ifs be the seemingly "more performant" inversion of the test logic, because
		//the mynborregions array is order such that when the cell is on its own region's boundary the first iteration is directly a hit
		//otherwise the more likely cases (face, corners) are tested for first, hence we almost never run this loop entirely for MYSELF + NUMBER_OF_NEIGHBORS iterations...
		if ( x >= mynborregions[(nb*IDS_AND_LIMITS)+THE_XMIN] && x <= mynborregions[(nb*IDS_AND_LIMITS)+THE_XMAX] && y >= mynborregions[(nb*IDS_AND_LIMITS)+THE_YMIN] && y <= mynborregions[(nb*IDS_AND_LIMITS)+THE_YMAX] && z >= mynborregions[(nb*IDS_AND_LIMITS)+THE_ZMIN] && z <= mynborregions[(nb*IDS_AND_LIMITS)+THE_ZMAX] )	{
				nbortid = mynborregions[(nb*IDS_AND_LIMITS)+THE_IDS];
				//nborpos = nb;
				break;
		}
	}

	//now nbortid holds into which border list we have to write
	//still in my own region but on the boundary?
	if ( nbortid == mynborregions[(0*IDS_AND_LIMITS)+THE_IDS] ) {

		uint32_t cxyz = (x - myGeom.nreg_rdmin) + ((y - myGeom.nreg_tdmin) * myGeom.nreg_rd) + ((z - myGeom.nreg_ndmin) * myGeom.nregarea_rdtd);

		uint32_t cxyz_deformedstate = mycellgrid[cxyz];

//##DEBUG#pragma omp critical { cout << "\t\t\t" << omp_get_thread_num() << ";" << "SAMEREGION-BORDER x,y,z,dx,dy,dz,cxyz,cxyz_deformedstate;nbortid;" << x << ";" << y << ";" << z << ";" << dx << ";" << dy << ";" << dz << ";" << cxyz << ";" << cxyz_deformedstate << ";" << nbortid << endl; }

		if ( cxyz_deformedstate >= reg_nmydefgpool || cxyz_deformedstate == CURRENTLY_INFECTED ) { //get out here, cell is already infected
//##DEBUG#pragma omp critical { cout << "\t\t\t" << omp_get_thread_num() << ";" << "SAMEREGION-BORDER returned" << endl; }
			return;
		}

		mycellgrid[cxyz] = CURRENTLY_INFECTED;

		uint32_t freeplace = this->omp_getNextFreeSlotInRXFrontBorder ( CELLRECYCLE ); //MK::strategy is RECYCLING because OutofRXFront infections are only called from updateFullRX

		myRXFrontBorder[freeplace].activity = ACTIVE;
		myRXFrontBorder[freeplace].infector = direction;
		myRXFrontBorder[freeplace].ix = x;
		myRXFrontBorder[freeplace].iy = y;
		myRXFrontBorder[freeplace].iz = z;
		myRXFrontBorder[freeplace].rxFrac = carryOver;

		//first assume most likely case - infection continues consuming the same deformed grain id, then P is the same
		myRXFrontBorder[freeplace].P = myRXFrontBorder[rxfrontid].P;

		//compromise growth of the grain from fully infected deformed regionn mydefgid into a different deformed grain
		if ( myRXFrontBorder[rxfrontid].mydefgid != cxyz_deformedstate ) {
			myRXFrontBorder[freeplace].P = reg_calc_mobilityweight ( rxgpoolid, cxyz_deformedstate );
		}
		if ( myRXFrontBorder[freeplace].P >= reg_myMobilityWeightMax ) {
			reg_myMobilityWeightMax = myRXFrontBorder[freeplace].P;
		}

		//cell carries consumed orientation along
		myRXFrontBorder[freeplace].myrxgid = rxgpoolid;
		myRXFrontBorder[freeplace].mydefgid = cxyz_deformedstate;
//##DEBUG#pragma omp critical{ cout << "\t\t\t" << omp_get_thread_num() << ";" << "SAMEREGION-BORDER success" << endl; }
		return;
	}


	//not returned yet, so the infection directs into a halo check in which and whether the corresponding location is 
	//free nborpos must not be utilized because multiple halos can lay in different domains!
	//find which halo of mine to place the infection
	uint32_t whichhalo = UNKNOWN_HALO;
	uint32_t hxmi, hxmx, hymi, hymx, hzmi, hzmx;
	for ( uint32_t h = 1; h < (MYSELF + NUMBER_OF_NEIGHBORS); h++ ) {
		if ( myhalos[h].inregion == nbortid ) {
			hxmi = myhalos[h].originx;
			hymi = myhalos[h].originy;
			hzmi = myhalos[h].originz;
			hxmx = hxmi + myhalos[h].extendx;
			hymx = hymi + myhalos[h].extendy;
			hzmx = hzmi + myhalos[h].extendz;
			if ( x >= hxmi && x <= hxmx && y >= hymi && y <= hymx && z >= hzmi && z <= hzmx ) {
				whichhalo = h;
				break;
			}
		}
	}
	QUICKASSERT( whichhalo != UNKNOWN_HALO );

	//is there a halo region for this neighbor?
	//QUICKASSERT ( myhalos[nborpos].inregion == mynborregions[(nborpos*IDS_AND_LIMITS)+THE_IDS] );
	//QUICKASSERT ( myhalos[nborpos].extendxyz > 0 );
	//there is a region so transform global coordinate into local halo coordinate

	short hx = x - myhalos[whichhalo].originx;
	short hy = y - myhalos[whichhalo].originy;
	short hz = z - myhalos[whichhalo].originz;
 
//#pragma omp critical
//{ cout << "issuing from my/nbtid/nbpos = " << this->mythreadid << "--" << nbortid << "--" << whichhalo << "==" << x << ";" << y << ";"<< z << "\t\t" << myhalos[whichhalo].originx << ";" << myhalos[whichhalo].originy << ";" << myhalos[whichhalo].originz << "\t\t" << myhalos[whichhalo].extendx << ";" << myhalos[whichhalo].extendy << ";" << myhalos[whichhalo].extendz << "\t\t" << hx << ";" << hy << ";" << hz << endl; }

	//##DEBUG
	QUICKASSERT ( hx >= 0 && hx < myhalos[whichhalo].extendx );
	QUICKASSERT ( hy >= 0 && hy < myhalos[whichhalo].extendy );
	QUICKASSERT ( hz >= 0 && hz < myhalos[whichhalo].extendz );

	uint32_t hid = hx + (hy * myhalos[whichhalo].extendx) + (hz * myhalos[whichhalo].extendxy);

	//this region does not know what is at the other side unless it checks in the memory shared with the thread
	//but as we are in sharedmemory one can check whether the halo was already infected
	if ( myhalos[whichhalo].theHaloCells[hid].status == HALOCELL_ALREADY_INFECTED ) {
		return; //infection unsuccessful
	}

	//HALOCELL is free mark the halocell, there is enough space and no reallocation necessary because the halo was preallocated
	myhalos[whichhalo].theHaloCells[hid].status = HALOCELL_ALREADY_INFECTED;
	myhalos[whichhalo].theHaloCells[hid].infector = direction;
	myhalos[whichhalo].theHaloCells[hid].ix = x;
	myhalos[whichhalo].theHaloCells[hid].iy = y;
	myhalos[whichhalo].theHaloCells[hid].iz = z;
	myhalos[whichhalo].theHaloCells[hid].rxFrac = carryOver;
	//during synchronization P is recalculated and then the neighbor finds potentially an even higher P
	myhalos[whichhalo].theHaloCells[hid].myrxgid = rxgpoolid;

	//make known that in this halo an infection was done
	uint32_t freeh = myhalos[whichhalo].nextFreeHaloRefs;
	myhalos[whichhalo].nextFreeHaloRefs = myhalos[whichhalo].nextFreeHaloRefs + 1;
	myhalos[whichhalo].theHaloRefs[freeh] = hid;

//#pragma omp critical{ cout << "\t\t\t" << omp_get_thread_num() << ";" << "OTHERREGION success" << endl; }
}


void caregionMemHdl::ompshar_sim_myCA_updateFullRXInside( void )
{
	double tprof = omp_get_wtime();
	//CALLED FROM WITHIN PARALLEL REGION
	//myFullRXListInside in the interval [0;nextSlotToFullRXInside) is a compact list of entries 
	//dereferencing cells from myRXFrontInside that in this time step have completely recrystallized, can infect and be INACTIVED thereafter
	//but it is required their categorization if the infected cell is still in the inside of a region or now handled by the border

	//repartitioning of recrystallized volume, original factors for velocity fgeo 1 in <100>, 1/2^0.5 <110>, 1/3^0.5 <111>
	double overFac1 = 0.0682;
	double overFac2 = 0.0542;
	double overFac3 = 0.0484;
	double overshoot, carryOver1, carryOver2, carryOver3;

	uint32_t nUpdatedCells = 0;
	uint32_t xmi = myGeom.nreg_rdmin;
	uint32_t xmx = myGeom.nreg_rdmax;
	uint32_t ymi = myGeom.nreg_tdmin;
	uint32_t ymx = myGeom.nreg_tdmax;
	uint32_t zmi = myGeom.nreg_ndmin;
	uint32_t zmx = myGeom.nreg_ndmax;
	uint32_t xxyy = myGeom.nregarea_rdtd;
	uint32_t xx = myGeom.nreg_rd;

	//MK::if a cell is located inside the automaton domain and a fixed kernel, e.g. MOORE is utilized
	//in almost all cases the infected cell needs no checking of periodic boundary conditions
	//MK::but we have to distinguish in every case into which region the infection goes
	bool InfectIntoBorder; //if true placing the cell in myRXFrontBorder, otherwise in myRXFrontInside
	uint32_t crxfrontid, rgpid, cxyz;

	for ( uint32_t c = 0; c < nextSlotToFullRXInside; c++ ) {
		crxfrontid = myFullRXListInside[c];

		rgpid = myRXFrontInside[crxfrontid].myrxgid;

		QUICKASSERT( myRXFrontInside[crxfrontid].rxFrac >= 1.0 );

		overshoot = myRXFrontInside[crxfrontid].rxFrac - 1.0;
		carryOver1 = overshoot * overFac1; //partition normalized overshoot
		carryOver2 = overshoot * overFac2;
		carryOver3 = overshoot * overFac3;

		//##MK::update potential shape axes-aligned bounding boxes for shape tracking --> inject code here

		//infect all untransformed neighbors
		//as all these infecting cells are located at the inside it is not necessary a boundary check for periodicity
		//the omp_infect_OutOfRXFrontInside is guided by infectBorder
		InfectIntoBorder = true;
		if ( myRXFrontInside[crxfrontid].ix > (xmi + 1) && myRXFrontInside[crxfrontid].ix < (xmx - 1) && myRXFrontInside[crxfrontid].iy > (ymi + 1) && myRXFrontInside[crxfrontid].iy < (ymx - 1) && myRXFrontInside[crxfrontid].iz > (zmi + 1) && myRXFrontInside[crxfrontid].iz < (zmx - 1) ) {
			InfectIntoBorder = false;
		}

		this->omp_infect_OutOfRXFrontInside ( rgpid, crxfrontid, InfectIntoBorder, -1, -1, 0, 16, carryOver2 );
		this->omp_infect_OutOfRXFrontInside ( rgpid, crxfrontid, InfectIntoBorder, 0, -1, 0, 15, carryOver1 );
		this->omp_infect_OutOfRXFrontInside ( rgpid, crxfrontid, InfectIntoBorder, +1, -1, 0, 14, carryOver2 );
		this->omp_infect_OutOfRXFrontInside ( rgpid, crxfrontid, InfectIntoBorder, -1, 0, 0, 13, carryOver1 );
		this->omp_infect_OutOfRXFrontInside ( rgpid, crxfrontid, InfectIntoBorder, +1, 0, 0, 12, carryOver1 );
		this->omp_infect_OutOfRXFrontInside ( rgpid, crxfrontid, InfectIntoBorder, -1, +1, 0, 11, carryOver2 );
		this->omp_infect_OutOfRXFrontInside ( rgpid, crxfrontid, InfectIntoBorder, 0, +1, 0, 10, carryOver1 );
		this->omp_infect_OutOfRXFrontInside ( rgpid, crxfrontid, InfectIntoBorder, +1, +1, 0, 9, carryOver2 );

		this->omp_infect_OutOfRXFrontInside ( rgpid, crxfrontid, InfectIntoBorder, -1, +1, -1, 19, carryOver3 );
		this->omp_infect_OutOfRXFrontInside ( rgpid, crxfrontid, InfectIntoBorder, 0, +1, -1, 18, carryOver2 );
		this->omp_infect_OutOfRXFrontInside ( rgpid, crxfrontid, InfectIntoBorder, +1, +1, -1, 17, carryOver3 );
		this->omp_infect_OutOfRXFrontInside ( rgpid, crxfrontid, InfectIntoBorder, -1, -1, -1, 25, carryOver3 );
		this->omp_infect_OutOfRXFrontInside ( rgpid, crxfrontid, InfectIntoBorder, 0, -1, -1, 24, carryOver2 );
		this->omp_infect_OutOfRXFrontInside ( rgpid, crxfrontid, InfectIntoBorder, +1, -1, -1, 23, carryOver3 );
		this->omp_infect_OutOfRXFrontInside ( rgpid, crxfrontid, InfectIntoBorder, -1, 0, -1, 22, carryOver2 );
		this->omp_infect_OutOfRXFrontInside ( rgpid, crxfrontid, InfectIntoBorder, 0, 0, -1, 21, carryOver1 );
		this->omp_infect_OutOfRXFrontInside ( rgpid, crxfrontid, InfectIntoBorder, +1, 0, -1, 20, carryOver2 );

		this->omp_infect_OutOfRXFrontInside ( rgpid, crxfrontid, InfectIntoBorder, -1, -1, +1, 8, carryOver3 );
		this->omp_infect_OutOfRXFrontInside ( rgpid, crxfrontid, InfectIntoBorder, 0, -1, +1, 7, carryOver2 );
		this->omp_infect_OutOfRXFrontInside ( rgpid, crxfrontid, InfectIntoBorder, +1, -1, +1, 6, carryOver3 );
		this->omp_infect_OutOfRXFrontInside ( rgpid, crxfrontid, InfectIntoBorder, -1, 0, +1, 5, carryOver2 );
		this->omp_infect_OutOfRXFrontInside ( rgpid, crxfrontid, InfectIntoBorder, 0, 0, +1, 4, carryOver1 );
		this->omp_infect_OutOfRXFrontInside ( rgpid, crxfrontid, InfectIntoBorder, +1, 0, +1, 3, carryOver2 );
		this->omp_infect_OutOfRXFrontInside ( rgpid, crxfrontid, InfectIntoBorder, -1, +1, +1, 2, carryOver3 );
		this->omp_infect_OutOfRXFrontInside ( rgpid, crxfrontid, InfectIntoBorder, 0, +1, +1, 1, carryOver2 );
		this->omp_infect_OutOfRXFrontInside ( rgpid, crxfrontid, InfectIntoBorder, +1, +1, +1, 0, carryOver3 );
		//MK::following code NEEDS to be executed AFTER all these infection attempts

		//deactivate cell 
		cxyz = (myRXFrontInside[crxfrontid].ix - xmi) + ((myRXFrontInside[crxfrontid].iy - ymi) * xx) + ((myRXFrontInside[crxfrontid].iz - zmi) * xxyy);

		//mark representor in the cellgrid as FULLY_RECRYSTALLIZED by assigning id of the recrystallizing
		mycellgrid[cxyz] = (reg_nmydefgpool + rgpid);

		//local texture bookkeeping
		reg_myrxgp_counts[myRXFrontInside[crxfrontid].myrxgid] += 1;
		reg_mydfgp_counts[myRXFrontInside[crxfrontid].mydefgid] -= 1;

		//cell flagged as "Freiwild" and thought of Green Mile healed from all infections...
		myRXFrontInside[crxfrontid].activity = INACTIVE;
		myRXFrontInside[crxfrontid].rxFrac = 0.0;

//#ifdef REPORTSTYLE_CELLCYCLES
//			cout << "\t\tstep->updt_afterinfects;crxfrontid;activity;ACTIVEwouldbeMarkedAs;nextRXSlotNeverActive\t" << mycaHdl->step << "\t" << crxfrontid << "\t" << myRXFront[crxfrontid].activity << "\t" << ACTIVE << "\t" << nextSlotNeverActiveRXFront << endl;
//#endif

		nUpdatedCells++;
	} //for all cells

	reg_nUpdatedCellsInside = nUpdatedCells;
	dtUpdateInside = omp_get_wtime() - tprof;
}


void caregionMemHdl::ompshar_clear_halo_bookkeeping( void )
{
	//CALLED FROM WITHIN PARALLEL REGION
	//it is necessary before a node can collect new infections in each integration step to clear the infection buffer
	for (uint32_t h = 0; h < (MYSELF+NUMBER_OF_NEIGHBORS); ++h ) {
		//simply set collector to zero, overwrites temporary buffers
		myhalos[h].nextFreeHaloRefs = 0;
	}
}


void caregionMemHdl::ompshar_sim_myCA_updateFullRXBorder( void )
{
	double tprof = omp_get_wtime();
	//CALLED FROM WITHIN PARALLEL REGION
	//MK::it is principal design the same as the updateFullRXInside, however now we must potentially direct infections 
	//out of our own memory region by offloading them into our halos

	double overFac1 = 0.0682;
	double overFac2 = 0.0542;
	double overFac3 = 0.0484;
	double overshoot, carryOver1, carryOver2, carryOver3;

	uint32_t nUpdatedCells = 0;
	uint32_t xmi = myGeom.nreg_rdmin;
	uint32_t ymi = myGeom.nreg_tdmin;
	uint32_t zmi = myGeom.nreg_ndmin;
	uint32_t xx = myGeom.nreg_rd;
	uint32_t xxyy = myGeom.nregarea_rdtd;
	uint32_t crxfrontid, rgpid, cxyz;

//##DEBUG#pragma omp critical { cout << "\t\tEntering updateFullRXBorder;step;nextSlotToFullRXBorder\t" << mycaHdl->step << "\t" << nextSlotToFullRXBorder << endl; }

	for ( uint32_t c = 0; c < nextSlotToFullRXBorder; c++ ) {
		crxfrontid = myFullRXListBorder[c];

		rgpid = myRXFrontBorder[crxfrontid].myrxgid; //runs from [0;myrxgpool.size) !

//##DEBUG#pragma omp critical { cout << "\t\t\t->step->Infecting BORDER with;step;c;crxfrontid;myRXFrontBorder[crxfrontid].rxFrac;nextSlotToFullRXBorder\t" << mycaHdl->step << "\t" << c << "\t" << crxfrontid << "\t" << myRXFrontBorder[crxfrontid].rxFrac << "\t" << nextSlotToFullRXBorder << endl; }

		QUICKASSERT( myRXFrontBorder[crxfrontid].rxFrac >= 1.0 );

		overshoot = myRXFrontBorder[crxfrontid].rxFrac - 1.0;
		carryOver1 = overshoot * overFac1; //partition normalized overshoot
		carryOver2 = overshoot * overFac2;
		carryOver3 = overshoot * overFac3;

		//##MK::update potential shape axes-aligned bounding boxes for shape tracking --> inject code here

		//three possible cases where the newly infected cells is located:
		//in the inside of this->region, in the border of this->region or in another region

		this->omp_infect_OutOfRXFrontBorder ( rgpid, crxfrontid, -1, -1, 0, 16, carryOver2 );
		this->omp_infect_OutOfRXFrontBorder ( rgpid, crxfrontid, 0, -1, 0, 15, carryOver1 );
		this->omp_infect_OutOfRXFrontBorder ( rgpid, crxfrontid, +1, -1, 0, 14, carryOver2 );
		this->omp_infect_OutOfRXFrontBorder ( rgpid, crxfrontid, -1, 0, 0, 13, carryOver1 );
		this->omp_infect_OutOfRXFrontBorder ( rgpid, crxfrontid, +1, 0, 0, 12, carryOver1 );
		this->omp_infect_OutOfRXFrontBorder ( rgpid, crxfrontid, -1, +1, 0, 11, carryOver2 );
		this->omp_infect_OutOfRXFrontBorder ( rgpid, crxfrontid, 0, +1, 0, 10, carryOver1 );
		this->omp_infect_OutOfRXFrontBorder ( rgpid, crxfrontid, +1, +1, 0, 9, carryOver2 );

		this->omp_infect_OutOfRXFrontBorder ( rgpid, crxfrontid, -1, +1, -1, 19, carryOver3 );
		this->omp_infect_OutOfRXFrontBorder ( rgpid, crxfrontid, 0, +1, -1, 18, carryOver2 );
		this->omp_infect_OutOfRXFrontBorder ( rgpid, crxfrontid, +1, +1, -1, 17, carryOver3 );
		this->omp_infect_OutOfRXFrontBorder ( rgpid, crxfrontid, -1, -1, -1, 25, carryOver3 );
		this->omp_infect_OutOfRXFrontBorder ( rgpid, crxfrontid, 0, -1, -1, 24, carryOver2 );
		this->omp_infect_OutOfRXFrontBorder ( rgpid, crxfrontid, +1, -1, -1, 23, carryOver3 );
		this->omp_infect_OutOfRXFrontBorder ( rgpid, crxfrontid, -1, 0, -1, 22, carryOver2 );
		this->omp_infect_OutOfRXFrontBorder ( rgpid, crxfrontid, 0, 0, -1, 21, carryOver1 );
		this->omp_infect_OutOfRXFrontBorder ( rgpid, crxfrontid, +1, 0, -1, 20, carryOver2 );

		this->omp_infect_OutOfRXFrontBorder ( rgpid, crxfrontid, -1, -1, +1, 8, carryOver3 );
		this->omp_infect_OutOfRXFrontBorder ( rgpid, crxfrontid, 0, -1, +1, 7, carryOver2 );
		this->omp_infect_OutOfRXFrontBorder ( rgpid, crxfrontid, +1, -1, +1, 6, carryOver3 );
		this->omp_infect_OutOfRXFrontBorder ( rgpid, crxfrontid, -1, 0, +1, 5, carryOver2 );
		this->omp_infect_OutOfRXFrontBorder ( rgpid, crxfrontid, 0, 0, +1, 4, carryOver1 );
		this->omp_infect_OutOfRXFrontBorder ( rgpid, crxfrontid, +1, 0, +1, 3, carryOver2 );
		this->omp_infect_OutOfRXFrontBorder ( rgpid, crxfrontid, -1, +1, +1, 2, carryOver3 );
		this->omp_infect_OutOfRXFrontBorder ( rgpid, crxfrontid, 0, +1, +1, 1, carryOver2 );
		this->omp_infect_OutOfRXFrontBorder ( rgpid, crxfrontid, +1, +1, +1, 0, carryOver3 );
		//MK::following code NEEDS to be executed AFTER all these infection attempts

		//deactivate cell, cell is definately within this->mycellgrid
		cxyz = (myRXFrontBorder[crxfrontid].ix - xmi) + ((myRXFrontBorder[crxfrontid].iy - ymi) * xx) + ((myRXFrontBorder[crxfrontid].iz - zmi) * xxyy);

		//mark representor in the cellgrid as FULLY_RECRYSTALLIZED by assigning id of RXgrains
		mycellgrid[cxyz] = (reg_nmydefgpool + rgpid);

		//local texture bookkeeping
		reg_myrxgp_counts[myRXFrontBorder[crxfrontid].myrxgid] += 1;
		reg_mydfgp_counts[myRXFrontBorder[crxfrontid].mydefgid] -= 1;

		//cell flagged as "Freiwild" and thought of Green Mile healed from all infections...
		myRXFrontBorder[crxfrontid].activity = INACTIVE;
		myRXFrontBorder[crxfrontid].rxFrac = 0.0;

//#ifdef REPORTSTYLE_CELLCYCLES
//			cout << "\t\tstep->updt_afterinfects;crxfrontid;activity;ACTIVEwouldbeMarkedAs;nextRXSlotNeverActive\t" << mycaHdl->step << "\t" << crxfrontid << "\t" << myRXFront[crxfrontid].activity << "\t" << ACTIVE << "\t" << nextSlotNeverActiveRXFront << endl;
//#endif

		nUpdatedCells++;
	} //for all cells

	reg_nUpdatedCellsBorder = nUpdatedCells;

//##DEBUG#pragma omp critical { cout << "\t\tLeaving updateFullRXBorder nUpdatedCells;\t" << reg_nUpdatedCellsBorder << endl; }
	dtUpdateBorder = omp_get_wtime() - tprof;
}


void caregionMemHdl::ompshar_synchronize_haloregions_rx( void )
{
	double tprof = omp_get_wtime();

	//scan over all myhalo regions that are expected at all updated
	uint32_t hal, nbortid;
	caregionMemHdlP apartreg;
	uint32_t HowManyToSynchronize;

	uint32_t xmi = myGeom.nreg_rdmin;
	uint32_t xmx = myGeom.nreg_rdmax;
	uint32_t ymi = myGeom.nreg_tdmin;
	uint32_t ymx = myGeom.nreg_tdmax;
	uint32_t zmi = myGeom.nreg_ndmin;
	uint32_t zmx = myGeom.nreg_ndmax;
	uint32_t xx = myGeom.nreg_rd;
	uint32_t yy = myGeom.nreg_td;
	uint32_t xxyy = xx*yy;
	uint32_t freeplace;
	uint32_t cxyz;
	uint32_t cxyz_defstate;
	uint32_t gx, gy, gz;
	//uint32_t oox, ooy, ooz;

	//sweep through my neighbors and synchronize with their halos
	for ( uint32_t nb = 0; nb < (MYSELF + NUMBER_OF_NEIGHBORS); nb++ ) {
		nbortid = mynborregions[(nb*IDS_AND_LIMITS)+THE_IDS];
		apartreg = mycaHdl->regions[nbortid]; //MK::threads can access each others data because of shared-memory!

		hal = apartreg->mynborhalos[nb];
		if ( hal != HALO_NOT_EXISTENT ) { //is there halo information?

//#pragma omp critical
//{ cout << "synchaloinfo threadid;nb;nbortid;hal" << this->mythreadid << ";" << nb << ";" << nbortid << ";" << hal << endl; }

			//get the regionid of the neighbor that has this halo
			HowManyToSynchronize = apartreg->myhalos[hal].nextFreeHaloRefs;
			if ( HowManyToSynchronize < 1 ) { continue; }

			//not continued, so there was sth in the halo requiring sync
			//hrefs = apartreg->myhalos[hal].theHaloRefs;
			//hcells = apartreg->myhalos[hal].theHaloCells;

			//oox = apartreg->myhalos[hal].originx;
			//ooy = apartreg->myhalos[hal].originy;
			//ooz = apartreg->myhalos[hal].originz;
 
			uint32_t target;
			for ( uint32_t s = 0; s < HowManyToSynchronize; s++ ) {
				target = apartreg->myhalos[hal].theHaloRefs[s];

				//read out halo coordinate
				gx = apartreg->myhalos[hal].theHaloCells[target].ix;
				gy = apartreg->myhalos[hal].theHaloCells[target].iy;
				gz = apartreg->myhalos[hal].theHaloCells[target].iz;

				
				QUICKASSERT( gx >= xmi && gx <= xmx);//##DEBUG
				QUICKASSERT( gy >= ymi && gy <= ymx);
				QUICKASSERT( gz >= zmi && gz <= zmx);
				//transform in local coordinate
				cxyz = (gx - xmi) + ((gy - ymi)*xx) + ((gz - zmi)*xxyy);
				cxyz_defstate = mycellgrid[cxyz];

				//eine Anekdote zum Thema, es sind nur zwei Zeilen Code.. tausche in folgender Zeile einmal || durch && und debugge mit Totalview, viel Spass
				if ( cxyz_defstate >= reg_nmydefgpool || cxyz_defstate == CURRENTLY_INFECTED ) {
					continue;
				}

				//not continued so halocell requires infection border
				mycellgrid[cxyz] = CURRENTLY_INFECTED;

				freeplace = this->omp_getNextFreeSlotInRXFrontBorder( CELLRECYCLE );
				QUICKASSERT( freeplace != INVALID_ADDRESS );

				myRXFrontBorder[freeplace].activity = ACTIVE;
				myRXFrontBorder[freeplace].infector = apartreg->myhalos[hal].theHaloCells[target].infector;
				myRXFrontBorder[freeplace].ix = apartreg->myhalos[hal].theHaloCells[target].ix;
				myRXFrontBorder[freeplace].iy = apartreg->myhalos[hal].theHaloCells[target].iy;
				myRXFrontBorder[freeplace].iz = apartreg->myhalos[hal].theHaloCells[target].iz;
				myRXFrontBorder[freeplace].rxFrac = apartreg->myhalos[hal].theHaloCells[target].rxFrac;
				myRXFrontBorder[freeplace].P = reg_calc_mobilityweight( apartreg->myhalos[hal].theHaloCells[target].myrxgid, cxyz_defstate);
				if ( myRXFrontBorder[freeplace].P >= reg_myMobilityWeightMax ) {
					reg_myMobilityWeightMax = myRXFrontBorder[freeplace].P;
				}
				myRXFrontBorder[freeplace].myrxgid = apartreg->myhalos[hal].theHaloCells[target].myrxgid;
				myRXFrontBorder[freeplace].mydefgid = cxyz_defstate;

				//in the region where the halo was filled updateFullRXBorder already marked the infecting cells as recrystallized
				//other bookkeep is not necessary --- in fact --- it were errornous
			} //next cell in the halo
		} //halo handled
	} //check whether next neighbor has a halo to sync for me

	dtSyncHalo = omp_get_wtime() - tprof;
}



void caregionMemHdl::ompshar_sim_myCA_clarify_status( double threshold )
{
	//called from within parallel region!
	//MK::SCORE model is implemented such that cell carries a) index of deformed grain, b) flag index CURRENTLY_INFECTED showing infected, c) index of rxgrain
	//thus, once a cell is infected only the myRXFrontLists know the initial and future assignment, therefore we require a function
	//which collects for all such cells the current grain assignment
	//most naively one could for each pixel to print when hitting a CURRENTLY_INFECTED scan all lists in all memregions to find that particular cell carrying the information
	//however this is O(N^2) complexity, better, compile and cache the assignments once and thereafter operate over all these cells to
	//modify within the image all those pixels at which an infection happened, resulting in O(N) complexity

	//first of all delete all old region entries, because the structure evolves, hence the assignments change
	myCellStatii.clear();

	//now interpret the assignments for each region; master analyzes regions sequentially could also be parallelized ##MK
	uint32_t ndefg = mydefgpool.size();
	uint32_t nbx = myCAGeometry.nboxedge_rd; //avoid the cache collision when calling this in the for loop to select the points
	uint32_t nby = myCAGeometry.nboxedge_td;
	uint32_t nbxy = nbx * nby;
	struct cellstatus ac;

	//analyze cells inside
	uint32_t ncandidates = nextSlotNeverActiveRXFrontInside;
	cellP front = myRXFrontInside;
	for ( uint32_t c = 0; c < ncandidates; c++ ) {
		if ( front[c].activity == ACTIVE ) { //most likely case for densely populated list
			ac.ixyz = front[c].ix + (front[c].iy * nbx) + (front[c].iz * nbxy); //ix, ... are short but positive! so short to uint32_t no problem
			if ( front[c].rxFrac < threshold ) {
				ac.gid = front[c].mydefgid;
				myCellStatii.push_back( ac );
				continue;
			}
			//else, rxed
			ac.gid = ndefg + front[c].myrxgid;
			myCellStatii.push_back( ac );
		}
	}

	//analyze cells border
	ncandidates = nextSlotNeverActiveRXFrontBorder;
	front = myRXFrontBorder;
	for ( uint32_t c = 0; c < ncandidates; c++ ) {
		if ( front[c].activity == ACTIVE ) { //most likely case for densely populated list
			ac.ixyz = front[c].ix + (front[c].iy * nbx) + (front[c].iz * nbxy); //ix, ... are short but positive! so short to uint32_t no problem
			if ( front[c].rxFrac < threshold ) {
				ac.gid = front[c].mydefgid;
				myCellStatii.push_back( ac );
				continue;
			}
			//else, rxed
			ac.gid = ndefg + front[c].myrxgid;
			myCellStatii.push_back( ac );
		}
	}
}


void caregionMemHdl::ompshar_sectioning_analysis_rx()
{
	//monitor how much area fraction is recrystallized in zpos used to control output in Martin/Markus MMM2018 model
	//corrected for excluding all rx grains whose sectioned area  is smaller than threshold

	uint32_t rx = myGeom.nreg_rd; //local coordinates [0, ri) so ri exclusive
	uint32_t ry = myGeom.nreg_td;
	uint32_t rz = myGeom.nreg_nd;
	uint32_t rxy = rx*ry;

	uint32_t imgxy = myCAGeometry.nboxedge_rd * myCAGeometry.nboxedge_td;
	uint32_t ndg = mydefgpool.size();
	uint32_t nrxg = myrxgpool.size();

	//clean buffer content and reset to zero counts
	rdtd_serialsection.clear();
	for( uint32_t zz = 0; zz < rz; zz++ )
		for ( uint32_t rxg = 0; rxg < nrxg; ++rxg )
			rdtd_serialsection.push_back( 0 );

	//add contribution from BULK CELLS
	uint32_t c, zroff, yzroff;
	for ( uint32_t zz = 0; zz < rz; zz++ ) {
		//for every of my z section we account now the sectioned size of the RX grains as well and filter for threshold
		//local cell coordinates of automaton
		zroff = zz*rxy;
		for ( uint32_t y = 0; y < ry; ++y ) {
			yzroff = y*rx + zroff;
			for ( uint32_t x = 0; x < rx; ++x ) {
				c = mycellgrid[x+yzroff]; //read cell index, recrystallized?
				if ( c <= CA_GRAINIDS_MAXIMUM ) { //so not CURRENTLY_INFECTED
					if ( c >= ndg ) {
						uint32_t rxgid = c - ndg;
						rdtd_serialsection[zz*nrxg+rxgid]++;
					}
				}
			} //scan +x line
		} //scan +y line
	} //scan +z lines


	//add contribution from CURRENTLY_INFECTED cells
	uint32_t zmi = myGeom.nreg_ndmin;
	uint32_t zmx = myGeom.nreg_ndmax;
	for( auto it = myCellStatii.begin(); it != myCellStatii.end(); ++it ) {
		//in which z layer is the infected pixel?
		uint32_t gridvalue = it->gid;

		//offset from global z to threadlocal z
		uint32_t zz = (it->ixyz / imgxy) - zmi;

		if ( gridvalue >= ndg ) {
			uint32_t rxgid = gridvalue - ndg;
			rdtd_serialsection[zz*nrxg+rxgid]++;
		}
	} //done checking all CURRENTLY_INFECTED candidates of the region
}


void caregionMemHdl::omp_log_synthfrontstats( void )
{
	//CALLED FROM WITHIN PARALLEL REGION
	struct omp_log_synthfstats memca;

	memca.localX = mycaHdl->X;
	memca.tCalcGrowthInside = this->dtCalcGrowthInside;
	memca.tCalcGrowthBorder = this->dtCalcGrowthBorder;
	memca.tUpdateInside = this->dtUpdateInside;
	memca.tUpdateBorder = this->dtUpdateBorder;
	memca.tSyncHalos = this->dtSyncHalo;
	memca.tSeqOverhead = this->dtSeqOverhead;

	memca.localstep = mycaHdl->step;
	memca.regionSvInside = this->reg_SvInside;
	memca.regionSvBorder = this->reg_SvBorder;
	memca.ntotalSeedInside = this->ntotalSeedFrontInside;
	memca.ntotalSeedBorder = this->ntotalSeedFrontBorder;
	memca.firstNeverActiveSeedInside = this->nextSlotNeverActiveSeedInside;
	memca.firstNeverActiveSeedBorder = this->nextSlotNeverActiveSeedBorder;
	//memca.nowActiveSeedInside = this->ncurrentActiveSeedInside;
	//memca.nowActiveSeedBorder = this->ncurrentActiveSeedBorder;

	//push_back copyconstruction of this stack struct
	profiling_synthmachine.push_back ( memca );
}


void caregionMemHdl::omp_log_rxfrontstats( void )
{
	//CALLED FROM WITHIN PARALLEL REGION
	struct omp_log_rxfstats memca;

	memca.localtime = mycaHdl->t;
	memca.localX = mycaHdl->X;
	memca.localmemory = regMemGuard;
	memca.localPmax = mycaHdl->myMobilityWeightMax;
	memca.localstep = mycaHdl->step;
	memca.regionSvInside = this->reg_SvInside;
	memca.regionSvBorder = this->reg_SvBorder;

	//query the time that elapsed as in this region some measurement was queried
	//MK::is at the moment utilized to get to know the time spend in per main loop cycle
	memca.tCalcGrowthInside = this->dtCalcGrowthInside;
	memca.tCalcGrowthBorder = this->dtCalcGrowthBorder;
	memca.tUpdateInside = this->dtUpdateInside;
	memca.tUpdateBorder = this->dtUpdateBorder;
	memca.tSyncHalos = this->dtSyncHalo;
	memca.tDefragment = this->dtDefragment;

	memca.ntotalRXFrontInside = this->ntotalRXFrontInside;
	memca.ntotalRXFrontBorder = this->ntotalRXFrontBorder;
	memca.nextSlotNeverActiveRXFrontInside = this->nextSlotNeverActiveRXFrontInside;
	memca.nextSlotNeverActiveRXFrontBorder = this->nextSlotNeverActiveRXFrontBorder;
	memca.nCurrentlyActiveInside = this->nCurrentlyActiveInside;
	memca.nCurrentlyActiveBorder = this->nCurrentlyActiveBorder;

	memca.ntotalFullRXListInside = this->ntotalFullRXListInside;
	memca.ntotalFullRXListBorder = this->ntotalFullRXListBorder;
	memca.nextSlotToFullRXInside = this->nextSlotToFullRXInside;

	memca.ntotalRecyclingListInside = this->ntotalRecyclingListInside;
	memca.ntotalRecyclingListBorder = this->ntotalRecyclingListBorder;
	memca.nextSlotThatBecomesRecycledInside = this->nextSlotThatBecomesRecycledInside;
	memca.nextSlotThatBecomesRecycledBorder = this->nextSlotThatBecomesRecycledBorder;
	memca.firstNotRecycledYetInside = this->firstNotRecycledYetInside;
	memca.firstNotRecycledYetBorder = this->firstNotRecycledYetBorder;

	//push_back copyconstruction of this stack struct
	profiling_growthmachine.push_back ( memca );
	//reset old values
	this->dtDefragment = 0.0;
}


void caregionMemHdl::omp_cleanupMemoryUsedForManagement( void )
{
	delete [] reg_myrxgp_counts;
	reg_myrxgp_counts = NULL;
	reg_nmyrxgp_counts = 0;

	delete [] reg_mydfgp_counts;
	reg_mydfgp_counts = NULL;
	reg_nmydfgp_counts = 0;
}


void caregionMemHdl::omp_cleanupMemoryUsedForGrowthMachine( void )
{
	//MK::when the solve_RXGROWTH has been executed the mygrainevolution contains the intermediate stage of all grains but 
	//the memory for the cellgrid and the dynamic memory allocated for managing cells in no longer necessary so delete it

	delete [] mycellgrid;
	mycellgrid = NULL;
	regMemGuard = regMemGuard - ( myGeom.nregvol_rdtdnd * sizeof(uint32_t));

	delete [] myFullRXListInside;
	myFullRXListInside = NULL;
	regMemGuard = regMemGuard - (ntotalFullRXListInside * sizeof(uint32_t));

	delete [] myFullRXListBorder;
	myFullRXListBorder = NULL;
	regMemGuard = regMemGuard - (ntotalFullRXListBorder * sizeof(uint32_t));

	delete [] myRecyclingListInside;
	myRecyclingListInside = NULL;
	regMemGuard = regMemGuard - (ntotalRecyclingListInside * sizeof(uint32_t));

	delete [] myRecyclingListBorder;
	myRecyclingListBorder = NULL;
	regMemGuard = regMemGuard - (ntotalRecyclingListBorder * sizeof(uint32_t));

	delete [] myRXFrontInside;
	myRXFrontInside = NULL;
	regMemGuard = regMemGuard - (ntotalRXFrontInside * sizeof(cell));

	delete [] myRXFrontBorder;
	myRXFrontBorder = NULL;
	regMemGuard = regMemGuard - (ntotalRXFrontBorder * sizeof(cell));


	delete [] mynborregions;
	mynborregions = NULL;
}


void caregionMemHdl::omp_cleanupMemoryUsedForHaloRegion( void )
{
	//first free memory in the halos
	for (uint32_t h = 0; h < (NUMBER_OF_NEIGHBORS + MYSELF); h++ ) {
		delete [] myhalos[h].theHaloRefs;
		myhalos[h].theHaloRefs = NULL;
		delete [] myhalos[h].theHaloCells;
		myhalos[h].theHaloCells = NULL;
	}

	//now the halos
	delete [] myhalos;
	myhalos = NULL;

	delete [] mynborhalos;
	mynborhalos = NULL;
}


caHdl::caHdl()
{
	//myensHdl nothing to free is only a reference to an address
	myensHdl = NULL;	//upon construction no manager has been assigned to that object yet!, no jobid = worldca-id as well
	jobid = 0;

	myCAGeometry.nNucleiCSR = 1;
	myCAGeometry.nboxedge_rd = CA_DIMENSIONS_MINIMUM;
	myCAGeometry.nboxedge_td = CA_DIMENSIONS_MINIMUM;
	myCAGeometry.nboxedge_nd = CA_DIMENSIONS_MINIMUM;
	myCAGeometry.nboxarea_rdtd = SQR(CA_DIMENSIONS_MINIMUM);
	myCAGeometry.nboxvol_rdtdnd = CUBE(CA_DIMENSIONS_MINIMUM);

	myCAGeometry.cellsize = DEFAULT_CELLSIZE;
	myCAGeometry.boxedge_rd = CA_DIMENSIONS_MINIMUM * DEFAULT_CELLSIZE;
	myCAGeometry.boxedge_td = CA_DIMENSIONS_MINIMUM * DEFAULT_CELLSIZE;
	myCAGeometry.boxedge_nd = CA_DIMENSIONS_MINIMUM * DEFAULT_CELLSIZE;
	myCAGeometry.boxarea_rdtd = SQR(CA_DIMENSIONS_MINIMUM)*SQR(DEFAULT_CELLSIZE);
	myCAGeometry.boxvol_rdtdnd = CUBE(CA_DIMENSIONS_MINIMUM)*CUBE(DEFAULT_CELLSIZE);

	myPhysData.G = DEFAULT_SHEARMODULUS;
	myPhysData.b = DEFAULT_BURGERSVECTOR;
	myPhysData.G0 = DEFAULT_SHEARMODULUS;
	myPhysData.dGdt = 0.0;
	myPhysData.bZeroCelsius = DEFAULT_BURGERSVECTOR;
	myPhysData.thermexp_C = DEFAULT_PUREALU_C;
	myPhysData.thermexp_a = 0.0;
	myPhysData.thermexp_b = 0.0;
	myPhysData.Tmelt = DEFAULT_PUREALU_TMELT;

	myPhysData.LAGBm0 = 1.0;
	myPhysData.LAGBHact = DEFAULT_LAGB_HACT * echarge;
	myPhysData.HAGBm0 = 1.0;
	myPhysData.HAGBHact = DEFAULT_HAGB_HACT * echarge;
	myPhysData.GSm0 = 1.0;
	myPhysData.GSHact = DEFAULT_GS_HACT * echarge;

	myPhysData.RH_HAGBm0 = 1.0;
	myPhysData.RH_HAGBHact = DEFAULT_HAGB_HACT * echarge;
	myPhysData.RH_LAGBHAGBcut = 0.9;
	myPhysData.RH_LAGBHAGBtrans = 5.0;
	myPhysData.RH_LAGBHAGBexponent = 9.0;


	myPhysData.defgmean_rd = DEFAULT_DEFGSIZE;
	myPhysData.defgmean_td = DEFAULT_DEFGSIZE;
	myPhysData.defgmean_nd = DEFAULT_DEFGSIZE;
	myPhysData.defgmean_poisson = DEFAULT_DEFGSIZE;

	myRecoveryModel.RecoveryConsider = RECOVERY_NO;
	myRecoveryModel.VacancyDiffGeo = 1.0;
	myRecoveryModel.VacancyDiffD0 = 1.0;
	myRecoveryModel.VacancyDiffHact = 1.0;
	myRecoveryModel.SoluteDiffD0 = 1.0;
	myRecoveryModel.SoluteDiffHact = 1.0;
	myRecoveryModel.SoluteLsPropFactor = 1.0;
	myRecoveryModel.SoluteConcentration = 1.0;
	myRecoveryModel.NesAlpha3 = 1.0;
	myRecoveryModel.NesKappa2 = 1.0;
	myRecoveryModel.NesC3 = 1.0;
	myRecoveryModel.NesC4 = 1.0;
	myRecoveryModel.NesC5 = 1.0;

	myDragData.ZenerConsider = DISPERSOIDDRAG_NO;
	myDragData.ZenerAlpha = 0.0;
	myDragData.ZenerGamma = 0.0;
	myDragData.fr = 0.0;

	myCADefMS.defmstype = CUBOID_DEFMS;
	myCADefMS.ngrx = 1;
	myCADefMS.ngry = 1;
	myCADefMS.ngrz = 1;
	myCADefMS.ngrxy = 1*1;
	myCADefMS.ngrxyz = 1*1*1;
	myCADefMS.u_xrd = 0.0;
	myCADefMS.v_ytd = 0.0;
	myCADefMS.w_znd = 0.0;

	myNucleationModel.gbnucleation = GBNUCLEATION_NO;
	myNucleationModel.csrnucleation = CSRNUCLEATION_NO;
	myNucleationModel.clustnucleation = CLUSTERNUCLEATION_NO;

	tmpdefgseeds = NULL;
	ndefgseeds = 0;

	mycellgrid = NULL;


	tprocessingend = 0.0;



	SvDeformed = 0;
	DefMicrotextureClasses = 0;
	NucleiMicrotextureClasses = 0;
	DefMicrotexture = NULL;
	NucleiMicrotexture = NULL;
	myMemGuard = 0.0;
	StoredEnergy = 0.0;
	CurrentG = 0.0;
	CurrentBurgersVector = DEFAULT_BURGERSVECTOR;
	CurrentTemperature = (25.0 + TOFFSET);
	t = 0.0;
	X = 0.0;
	Xcells = 0.0;
	dXstep = 0.0;
	dt = DEFAULT_DELTATIME;
	tsimend = 0.0;
	stepsimend = 0;
	step = 0;
	XMAX = DEFAULT_XMAX;
	TMAX = DEFAULT_TMAX;
	NTSTEPSMAX = DEFAULT_NMAX;
	myMobilityWeightMax = DEFAULT_PMAX;	//start assuming always HAGB, system will detect whether there are faster boundaries that need finer time discretization
	myrhomax = RHOMAX_WELLANNEALED;
	myrhomax0 = RHO_DEFORMED_MATERIAL_MAX;
	myrhomin0 = RHO_RECRYSTALLIZED_MATERIAL;
	mypzmin = 0.0;


	Gbhalfsq = 0.5 * myPhysData.G0 * SQR(myPhysData.bZeroCelsius);
	//_kT = (1.0 / (kboltzman * CurrentTemperature));
	mLAGB = myPhysData.LAGBm0 * exp( -1.0 *  DEFAULT_LAGB_HACT / (kboltzman * CurrentTemperature) );
	mHAGB = myPhysData.HAGBm0 * exp( -1.0 *  DEFAULT_HAGB_HACT / (kboltzman * CurrentTemperature) );
	mGS = myPhysData.GSm0 * exp( -1.0 *  DEFAULT_GS_HACT / (kboltzman * CurrentTemperature) );
	mRHHAGB = myPhysData.RH_HAGBm0 * exp( -1.0 * myPhysData.RH_HAGBHact / (kboltzman * CurrentTemperature) );

	nmydefgpool = 0;
	nmyrxgpool = 0;
	nmynuclei = 0;

	loginfo_rxfrontstats_cnt = 0;
	loginfo_mpm_cnt = 0;
	loginfo_grainevo_cnt = 0;
	loginfo_rendering_cnt = 0;
	loginfo_rendering_sectionbased_cnt = vector<uint32_t>( 3, 0 );
	loginfo_defrag_cnt = 0;
	loginfo_percolation_cnt = 0;

	myensRank = MASTER;
	nRanks = 1;

	outopt_localrenderhow = RENDERING_MSNO;
	outopt_localrendercolor = RENDERING_COLOR_GRAINID;
	outopt_localrenderfileformat = RENDERING_FILEFORMAT_RAW;
	outopt_localrenderboundaries = RENDERING_BOUNDARIES_NO;
	outopt_logboundaries = OUTPUT_LOGBND_NO;
	outopt_localrxfront = OUTPUT_RXFRONTSTATS_NO;
	outopt_localthreadprof = OUTPUT_THREADPROFILING_NO;
	outopt_localsinglegrainevo = OUTPUT_SINGLEGRAIN_NO;
	outopt_localgenDAMASKgeom = OUTPUT_DAMASK_GEOMETRYFILE_NO;
	outopt_localartsemebsd = OUTPUT_SEMEBSD_NO;

	mySuccess = true;
	onthefly_defragmentation = false;
	renderingForThisCA = false;

	localprng.init( DEFAULT_PRNG_SEED );
}


caHdl::~caHdl()
{
	/* ##MK::check for deallocation of dynmaic memory to store status of the grains
	//because caHdl does not know about dynamic allocation of memory from within the struct
	for ( uint32_t s = 0; s < mygrainevolution.size(); s++) {
		delete [] mygrainevolution[s].localdatabucket;
		mygrainevolution[s].localdatabucket = NULL;
		myMemGuard = myMemGuard - ( (double) mygrainevolution[s].nlocaldata * sizeof(uint32_t) );
	}
	//vector vector does allocate self*/
}


//now all caHdl helper functions
void caHdl::update_system_fixedtime( double time )
{
	//##MK::implement here a presolver that evolves the microstructure
	//to the earliest point in time when nucleation occurs to safe integration steps and as such reduce the simulation time
}


void caHdl::update_temperature( void )
{
	uint32_t i = 0;
	//scan vector t and perform linear interpolation to get Ti(ti)
	while (ttime[i] < t && i < ttime.size() ) {
		i++;
	}

	if (i == 0) {
		CurrentTemperature = ttemperature[0];
		return;
	}
	if (i >= ttime.size() ) {
		CurrentTemperature = ttemperature[ttime.size() - 1];
		return;
	}

	CurrentTemperature = ttemperature[i-1] + ( (ttemperature[i] - ttemperature[i-1]) / (ttime[i] - ttime[i-1]) ) * (t - ttime[i-1]);
}


void caHdl::update_atomisticproperties( void )
{
	CurrentG = get_shearmodulus( CurrentTemperature );
	myPhysData.G = CurrentG;

	CurrentBurgersVector = get_burgersvector( CurrentTemperature );
	myPhysData.b = CurrentBurgersVector;

	//##add here further any temperature dependent values which need update, e.g. recovery parameters

	//cache Gbhalfsq
	Gbhalfsq = 0.5 * myPhysData.G * SQR( myPhysData.b );
}


void caHdl::update_intrinsicmobilities( void )
{
	//intrinsic grain boundary mobilities
	mLAGB = myPhysData.LAGBm0 * exp( -1.0 *  myPhysData.LAGBHact / (kboltzman * CurrentTemperature) );
	mHAGB = myPhysData.HAGBm0 * exp( -1.0 *  myPhysData.HAGBHact / (kboltzman * CurrentTemperature) );
	mGS = myPhysData.GSm0 * exp( -1.0 *  myPhysData.GSHact / (kboltzman * CurrentTemperature) );

	mRHHAGB = myPhysData.RH_HAGBm0 * exp( -1.0 * myPhysData.RH_HAGBHact / (kboltzman * CurrentTemperature) );

	QUICKASSERT ( mLAGB > 0.0 );
	QUICKASSERT ( mHAGB >= mLAGB );
	QUICKASSERT ( mGS >= mHAGB );
	QUICKASSERT ( mRHHAGB > 0.0 );
}


void caHdl::update_microchemistry( void )
{
	//update the Zener Drag
	if ( myDragData.ZenerConsider == DISPERSOIDDRAG_NO ) {
		myDragData.fr = 0.0;
		this->mypzmin = 0.0;
		return;
	}

	if ( myDragData.ZenerConsider == DISPERSOIDDRAG_CONSTANT ) {
		//MK::convention take the first entry from EvoDraggingParticles
		myDragData.fr = zenerdispersion[0];
		this->mypzmin = myDragData.ZenerAlpha * myDragData.ZenerGamma * myDragData.fr;
		return;
	}


	//DISPERSOIDDRAG_TIMEDEP
	uint32_t i = 0; //scan vector t and get current linear fr
	while ( zenertime[i] < t && i < zenertime.size() ) {
		i++;
	}

	if ( i == 0 ) {
		myDragData.fr = 0.0;
	}

	//linear interpolate
	if (i < zenertime.size() ) {
		double dfrdt = ( (zenerdispersion[i] - zenerdispersion[i-1]) / (zenertime[i] - zenertime[i-1]) );
		myDragData.fr = zenerdispersion[i-1] + dfrdt * ( this->t - zenertime[i-1] );
	}


	this->mypzmin = myDragData.ZenerAlpha * myDragData.ZenerGamma * myDragData.fr;
	//##MK::-->inject code for update of solute drag or any other parts of a microchemistry model here
}


void caHdl::update_deformedsubstructure( void )
{
	//##MK::inject the recovery model here which updates rho
	//evolution of the statistical dislocation density affecting mydefgpool[i].rho
	unsigned char recoverymodel = myRecoveryModel.RecoveryConsider;

	if ( recoverymodel == RECOVERY_NO ) 
		return;

	if ( recoverymodel == RECOVERY_NES_VACCOREDIFF ) {
		nes_networkgrowthmodel_vacancycorediff();
		return;
	}


	if ( recoverymodel == RECOVERY_NES_SOLUTEDRAG ) {
		nes_networkgrowthmodel_solutedrag();
		return;
	}

	if ( recoverymodel == RECOVERY_NES_MICHALAKPAXTON ) {
		nes_michalakpaxton();
		return;
	}

	//evolution of the average subgrain size ...

	//evolution of the average misorientation in the deformed grain affecting mydefgpool[i].dgav
}



double caHdl::get_rho ( uint32_t mydefgpoolid )
{
	//MK::interpret mydefgpool rho
	return mydefgpool[mydefgpoolid].rho;
}


double caHdl::get_zener( void )
{
	if ( myDragData.ZenerConsider == DISPERSOIDDRAG_NO ) 
		return 0.0;

	//DISPERSOIDDRAG_CONSTANT || DISPERSOIDDRAG_TIMEDEP
	//MK::at the moment constant drag in each grain
	return ( 1.0 * myDragData.ZenerAlpha * myDragData.ZenerGamma * myDragData.fr );
}


double caHdl::get_currentintrinsicmobility( double Pvalue )
{
	if ( mobilitymodel == MOBILITYMODEL_ROLLETTHOLM )
		return ( Pvalue * mRHHAGB );

	//else MOBILITYMODEL_SEBALDGOTTSTEIN
	if (Pvalue < 0.0) {
		return mLAGB;
	}
	return ( ( Pvalue * mGS ) + ( (1.0 - Pvalue) * mHAGB ) );
	//most likely case Pvalue >= 0.0 but then two comparisons...
}


uint32_t caHdl::characterize_cellstatus( uint32_t r, uint32_t lx, uint32_t ly, uint32_t lz )
{
	//finds interface cell at local coordinates lx,ly,lz and global coordinates
	//lx,ly, lz are [0,regionlimit)
	short x = ((short) this->regions[r]->myGeom.nreg_rdmin) + ((short) lx); //first inclusive location in memory region r in x direction
	short y = ((short) this->regions[r]->myGeom.nreg_tdmin) + ((short) ly);
	short z = ((short) this->regions[r]->myGeom.nreg_ndmin) + ((short) lz);
	//MK::conversion from uint32_t li to short possible as extend of automaton well below a short's positive range...

	//scan the two container's for active interface cells, inside and border
	cellP bucket = this->regions[r]->myRXFrontInside;
	uint32_t ncand = this->regions[r]->nextSlotNeverActiveRXFrontInside;

	for ( uint32_t c = 0; c < ncand; c++ ) {
		if ( bucket[c].activity == ACTIVE ) {
			if ( bucket[c].ix != x ) continue;
			if ( bucket[c].iy != y ) continue;
			if ( bucket[c].iz != z ) continue;
			//not continued so cell was found identify status

			//##MK::assure that the following definition is always consistent with the one in caHdl::binarize_partiallyrxed_microstructure()
			//CELL_IS_INFECTED, i.e. recrystallized volume
			return (nmydefgpool + bucket[c].myrxgid);
		}
	}

	bucket = this->regions[r]->myRXFrontBorder;
	ncand = this->regions[r]->nextSlotNeverActiveRXFrontBorder;
	for ( uint32_t c = 0; c < ncand; c++ ) {
		if ( bucket[c].activity == ACTIVE ) {
			if ( bucket[c].ix != x ) continue;
			if ( bucket[c].iy != y ) continue;
			if ( bucket[c].iz != z ) continue;
			//not continued so cell was found identify status

			//##MK::assure that the following definition is always consistent with the one in caHdl::binarize_partiallyrxed_microstructure()
			//CELL_IS_INFECTED, i.e. recrystallized volume
			return (nmydefgpool + bucket[c].myrxgid);
		}
	}

	//not returned yet, so cell was not found!
	return NO_GRAIN_ASSIGNED;
}


double caHdl::get_dtmax_instantslope_cells( void ) //double when )
{
	double dtmax = INFINITE;

	if ( mobilitymodel == MOBILITYMODEL_ROLLETTHOLM ) {
		double mRHHAGBmax = myPhysData.RH_HAGBm0 * exp( - 1.0 * myPhysData.RH_HAGBHact / (kboltzman * CurrentTemperature) );

		double vRHmax = mRHHAGBmax * ((Gbhalfsq * myrhomax) - mypzmin);

		dtmax = ((maxfillperstep * myCAGeometry.cellsize) / vRHmax);

		return dtmax;
	}


	//MOBILITYMODEL_SEBALDGOTTSTEIN
	//maximum intrinsic mobility
	//current not absolute somewhen, as we can utilize explicit time-synchronization
	double mgs = myPhysData.GSm0 * exp( -1.0 * myPhysData.GSHact / (kboltzman * CurrentTemperature) );
	double mhagb = myPhysData.HAGBm0 * exp ( -1.0 * myPhysData.HAGBHact / (kboltzman * CurrentTemperature) );

	//but there might not at all be a nucleus-defg achieving this high mobility to completely or there are only
	//only LAGB boundaries in the system
	double m_max = ( myMobilityWeightMax * mgs ) + ( (1.0 - myMobilityWeightMax) * mhagb );

	//only LAGB in the system
	if ( myMobilityWeightMax < (0.0 - DEFAULT_SMALL_NUMBER) ) {
		m_max = myPhysData.LAGBm0 * exp( -1.0 * myPhysData.LAGBHact / (kboltzman * CurrentTemperature) );
	}

	//maximum migration velocity
	double v_max = m_max * ((Gbhalfsq * myrhomax) - mypzmin);

	//transform only a fraction of a cell in a time step
	dtmax = ((maxfillperstep * myCAGeometry.cellsize) / v_max);

	//#ifdef DETAILED_PROMPTScout << this->jobid << ";jobid;" << m_max << ";m_max;" << v_max << ";vmax;" << tmin << ";tmin" << endl,#endif
	return dtmax;
}


double caHdl::get_dtmax_instantslope_temp( void )
{
	uint32_t i = 0; //scan vector t and get current linear heating rate
	double dtmax = INITIAL_DELTAT;

	while ( ttime[i] < t && i < ttime.size() ) {
		i++;
	}

	if ( i == 0 ) {
		return dtmax;
	}

	if (i < ttime.size() ) {
		double dTdt = ( fabs(ttemperature[i] - ttemperature[i-1]) / (ttime[i] - ttime[i-1]) ); 

		if (dTdt <= SMALL_HEATRATE) dTdt = SMALL_HEATRATE; //avoid division by zero error

		dtmax = (SMALL_HEAT / dTdt);

		//but not more than the distance of the time step intervals!
		double dtintvl = ttime[i] - ttime[i-1];
		if ( dtmax >= dtintvl )
			dtmax = dtintvl;

		return dtmax;
	}

	//if not already returned, then processing schedule ends, temperature stays constant with last value, nevertheless so dT/dt = 0
	dtmax = INFINITE;
	return dtmax;
}


double caHdl::get_dtmax_instantslope_rho( void )
{
	unsigned char recoverymodel = myRecoveryModel.RecoveryConsider;

	if ( recoverymodel == RECOVERY_NO )
		return INFINITE;

	double drdtmax, dtmax;
	if ( recoverymodel == RECOVERY_NES_VACCOREDIFF ) {
		double Bx = _PI_ * myRecoveryModel.NesC3 * myRecoveryModel.NesC5 * myRecoveryModel.VacancyDiffD0 * exp ( -1.0 * myRecoveryModel.VacancyDiffHact / ( kboltzman * CurrentTemperature ));
		Bx = Bx / (myRecoveryModel.NesKappa2 * ( 1.0 / pow(myrhomax, 0.5) ) );

		double exponent = SQR(CurrentBurgersVector) / (kboltzman * CurrentTemperature);

		drdtmax = Bx * exp ( myRecoveryModel.NesAlpha3 * myRecoveryModel.NesKappa2 * CurrentG * CurrentBurgersVector * exponent );

		if ( drdtmax < (1.0 / pow(SMALL_RECOVERYRATE, 0.5)) ) {
			drdtmax = (1.0 / pow(SMALL_RECOVERYRATE, 0.5));
		}

		dtmax = (1.0 / pow(SMALL_RHO, 0.5) ) / drdtmax;
		return dtmax;
	}

	if ( recoverymodel == RECOVERY_NES_SOLUTEDRAG ) {
		double Bsa =  myRecoveryModel.NesC3 * CurrentBurgersVector * myRecoveryModel.SoluteDiffD0 * exp ( -1.0 * myRecoveryModel.SoluteDiffHact / ( kboltzman * CurrentTemperature ));

		double exponent = SQR(CurrentBurgersVector) / (kboltzman * CurrentTemperature);
		double ls = myRecoveryModel.SoluteLsPropFactor * pow( myRecoveryModel.SoluteConcentration, (-2.0 / 3.0) );
		drdtmax = Bsa * exp ( myRecoveryModel.NesAlpha3 * ls * CurrentG * CurrentBurgersVector * exponent / (1.0 / SQR(this->myrhomax)) );

		if ( drdtmax < (1.0 / pow(SMALL_RECOVERYRATE, 0.5)) ) {
			drdtmax = (1.0 / pow(SMALL_RECOVERYRATE, 0.5));
		}

		dtmax = (1.0 / pow(SMALL_RHO, 0.5) ) / drdtmax;
		return dtmax;
	}


	//unrecognized recovery model
	return INFINITE;
}


double caHdl::get_dtmax_instantslope_zener( void )
{
	if ( myDragData.ZenerConsider != DISPERSOIDDRAG_TIMEDEP ) { //when there is no ZenerDrag simulated or the drag is not time-dependent there is no necessity to adjust the numerical time integration of the CA
		return INFINITE;
	}

	uint32_t i = 0;
	double dtmax = INITIAL_DELTAT;

	while ( zenertime[i] < t && i < zenertime.size() ) {
		i++;
	}

	if ( i == 0 ) {
		return dtmax;
	}

	if (i < zenertime.size() ) {
		double dpzdt = myDragData.ZenerAlpha * myDragData.ZenerGamma * ( fabs(zenerdispersion[i] - zenerdispersion[i-1]) / (zenertime[i] - zenertime[i-1]) );

		if ( fabs(dpzdt) <= SMALL_DRAGGINGRATE) dpzdt = SMALL_DRAGGINGRATE; //avoid division by zero error

		dtmax = (SMALL_ZENERFORCE / dpzdt);
		return dtmax;
	}

	//if not already returned, then processing schedule has ended
	dtmax = INFINITE;
	return dtmax;
}


double caHdl::get_dtmax_instantslope_nucleation( void )
{
	return INFINITE;

	/*if ( this->myNucleationModel.tincubmodel != TINCUB_TIMEDEPENDENT )
		return INFINITE;

	//nucleation is time dependent
	double dtmax = INFINITE;
	double dNdt = (double) myrxgpool.size() * pdf_rayleigh( t ); //##MK::all planned nuclei, MK::not t+dt because function is being called after this->t was incremented!
	if ( dNdt <= SMALL_NUCLEATIONRATE ) dNdt = SMALL_NUCLEATIONRATE; //avoid division by zero error
	dtmax = (SMALL_NUMBEROFNUCLEI / dNdt);
	return dtmax;*/
}


double caHdl::get_dtmax_minimum( void )
{
	double min_dtmax = INFINITE;
	double dtmax_model = INFINITE;

	//analyze all physical mechanisms that need proper time discretization during the integration scheme
	dtmax_model = get_dtmax_instantslope_temp();
	if ( dtmax_model < min_dtmax ) 
		min_dtmax = dtmax_model;

	dtmax_model = get_dtmax_instantslope_rho();
	if ( dtmax_model < min_dtmax ) 
		min_dtmax = dtmax_model;

	dtmax_model = get_dtmax_instantslope_zener();
	if ( dtmax_model < min_dtmax ) 
		min_dtmax = dtmax_model;

	double dtmax_migr_discr = get_dtmax_instantslope_cells();
	if ( dtmax_migr_discr < min_dtmax )
		min_dtmax = dtmax_migr_discr;

	//##MK::modification
	double dtmax_nucrate_discr = get_dtmax_instantslope_nucleation();

	//initially min_dtmax was scaled such that not more than one nucleus
	//difference per integration step was seeded, however, this voted for orders of 
	//magnitude shorter time steps for timedependent nucleation, however, if the migration distance of a cell is limited, one could also 
	//without significant effect on the simulation results seed all nuclei at once, as they cannot migrate more than a partial step of a cell
	if ( dtmax_nucrate_discr >=  dtmax_migr_discr ) {
		if ( dtmax_migr_discr < min_dtmax )
			min_dtmax = dtmax_migr_discr;
	}
	//not returned indicated nucleation rate is very high and requires finer

	return min_dtmax;
}





void caHdl::init_cahdlprng( void )
{
	long localseed = (long) -1*( this->jobid  + 1); //+1 because jobid can be 0, MK::do not set seeds like myensHdl->myRank because then all SUs from the rank would reinitialize the same PRNG stream!
	if ( myensHdl->ensembleprngseed == DEFAULT_SEED ) {
		localseed = localseed - 10;
		cout << "This jobid " << this->jobid << " has localprng seed " << localseed << endl;
	}
	else {
		localseed = myensHdl->ensembleprngseed;
		cout << "This jobid " << this->jobid << " has localprng seed " << localseed << endl;
	}

	this->localprng.init( localseed );
	//set as well the seed for the local math library
	//thus each caDomain initializes different random number streams per solitary unit
	this->localmath.r.init( localseed );
}


void caHdl::init_parameter( void )
{
	//beneficial, because I/O only for each process, not each automaton...
	ensembleHdlP myens = this->myensHdl; //MK::allows class object caHdl to refer to his executing MPI rank!

	XMAX = myens->XMAX;
	TMAX = myens->TMAX;
	NTSTEPSMAX = myens->NTSTEPSMAX;

	myCAGeometry.cellsize = myens->ensCAGeometry.cellsize; //MK::within my thesis all SUs of the same geometry
	this->_cellsize = (1.0 / myCAGeometry.cellsize);
	myCAGeometry.nNucleiCSR = myens->ensCAGeometry.nNucleiCSR;
	myCAGeometry.nboxedge_rd = myens->ensCAGeometry.nboxedge_rd;
	myCAGeometry.nboxedge_td = myens->ensCAGeometry.nboxedge_td;
	myCAGeometry.nboxedge_nd = myens->ensCAGeometry.nboxedge_nd;
	myCAGeometry.nboxarea_rdtd = myens->ensCAGeometry.nboxarea_rdtd;
	myCAGeometry.nboxvol_rdtdnd = myens->ensCAGeometry.nboxvol_rdtdnd;
	myCAGeometry.boxedge_rd = myens->ensCAGeometry.boxedge_rd;
	myCAGeometry.boxedge_td = myens->ensCAGeometry.boxedge_td;
	myCAGeometry.boxedge_nd = myens->ensCAGeometry.boxedge_nd;
	myCAGeometry.boxarea_rdtd = myens->ensCAGeometry.boxarea_rdtd;
	myCAGeometry.boxvol_rdtdnd = myens->ensCAGeometry.boxvol_rdtdnd;

	myPhysData.G = myens->ensPhysData.G;
	myPhysData.b = myens->ensPhysData.b;
	myPhysData.G0 = myens->ensPhysData.G0;
	myPhysData.dGdt = myens->ensPhysData.dGdt;
	myPhysData.bZeroCelsius = myens->ensPhysData.bZeroCelsius;
	myPhysData.thermexp_C = myens->ensPhysData.thermexp_C;
	myPhysData.thermexp_a = myens->ensPhysData.thermexp_a;
	myPhysData.thermexp_b = myens->ensPhysData.thermexp_b;
	myPhysData.Tmelt = myens->ensPhysData.Tmelt;

	myPhysData.LAGBm0 =	myens->ensPhysData.LAGBm0;
	myPhysData.LAGBHact = myens->ensPhysData.LAGBHact;
	myPhysData.HAGBm0 = myens->ensPhysData.HAGBm0;
	myPhysData.HAGBHact = myens->ensPhysData.HAGBHact;
	myPhysData.GSm0 = myens->ensPhysData.GSm0;
	myPhysData.GSHact = myens->ensPhysData.GSHact;

	myPhysData.RH_HAGBm0 = myens->ensPhysData.RH_HAGBm0;
	myPhysData.RH_HAGBHact = myens->ensPhysData.RH_HAGBHact;
	myPhysData.RH_LAGBHAGBcut = myens->ensPhysData.RH_LAGBHAGBcut;
	myPhysData.RH_LAGBHAGBtrans = myens->ensPhysData.RH_LAGBHAGBtrans;
	myPhysData.RH_LAGBHAGBexponent = myens->ensPhysData.RH_LAGBHAGBexponent;

	myPhysData.defgmean_rd = myens->ensPhysData.defgmean_rd;
	myPhysData.defgmean_td = myens->ensPhysData.defgmean_td;
	myPhysData.defgmean_nd = myens->ensPhysData.defgmean_nd;
	myPhysData.defgmean_poisson = myens->ensPhysData.defgmean_poisson;

	myNucleationModel.gbnucleation = myens->ensNucleationModel.gbnucleation;
	myNucleationModel.csrnucleation = myens->ensNucleationModel.csrnucleation;
	myNucleationModel.clustnucleation = myens->ensNucleationModel.clustnucleation;
	myNucleationModel.tincubmodel = myens->ensNucleationModel.tincubmodel;
	myNucleationModel.defaultnucdensity = myens->ensNucleationModel.defaultnucdensity;

	myNucleationModel.cluster_nclust = myens->ensNucleationModel.cluster_nclust;
	myNucleationModel.cluster_lambda = myens->ensNucleationModel.cluster_lambda;
	myNucleationModel.cluster_rvesize = myens->ensNucleationModel.cluster_rvesize;
	myNucleationModel.cluster_a = myens->ensNucleationModel.cluster_a;
	myNucleationModel.cluster_b = myens->ensNucleationModel.cluster_b;
	myNucleationModel.cluster_c = myens->ensNucleationModel.cluster_c;
	myNucleationModel.gbnucleation_dens2num = myens->ensNucleationModel.gbnucleation_dens2num;
	myNucleationModel.gbnucleation_drho2dens = myens->ensNucleationModel.gbnucleation_drho2dens;
	myNucleationModel.gbnucleation_scatter = myens->ensNucleationModel.gbnucleation_scatter;
	myNucleationModel.tincub_rayleigh_sigma = myens->ensNucleationModel.tincub_rayleigh_sigma;

	percolation = myens->percolation;
	mobilitymodel = myens->mobilitymodel;
	nucleationmodel_status = NUCSITES_STILL_SOME_FREE;
	outopt_localrenderhow = myens->outopt_rendermethod;
	outopt_localrendercolor = myens->outopt_rendercolormodel;
	outopt_localrenderfileformat = myens->outopt_renderfileformat;
	outopt_localrenderboundaries = myens->outopt_localrenderboundaries;
	outopt_logboundaries = myens->outopt_logboundaries;
	outopt_localthreadprof = myens->outopt_threadprofiling;
	outopt_localrxfront = myens->outopt_rxfront;
	outopt_localsinglegrainevo = myens->outopt_singlegrainevo;
	outopt_localgenDAMASKgeom = myens->outopt_generateDAMASKgeom;
	outopt_localartsemebsd = myens->outopt_artsemebsd;

	maxfillperstep = myens->maxfillperstep;
	initialRelCellCaching = myens->initialRelCellCaching;
	transientRelCellRecaching = myens->transientRelCellRecaching;
	ebsdstepsize = myens->ebsdstepsize;

	myRecoveryModel.RecoveryConsider = myens->ensRecoveryModel.RecoveryConsider;
	myRecoveryModel.VacancyDiffGeo = myens->ensRecoveryModel.VacancyDiffGeo;
	myRecoveryModel.VacancyDiffD0 = myens->ensRecoveryModel.VacancyDiffD0;
	myRecoveryModel.VacancyDiffHact = myens->ensRecoveryModel.VacancyDiffHact;
	myRecoveryModel.SoluteDiffD0 = myens->ensRecoveryModel.SoluteDiffD0;
	myRecoveryModel.SoluteDiffHact = myens->ensRecoveryModel.SoluteDiffHact;
	myRecoveryModel.SoluteLsPropFactor = myens->ensRecoveryModel.SoluteLsPropFactor;
	myRecoveryModel.SoluteConcentration = myens->ensRecoveryModel.SoluteConcentration;
	myRecoveryModel.NesAlpha3 = myens->ensRecoveryModel.NesAlpha3;
	myRecoveryModel.NesKappa2 = myens->ensRecoveryModel.NesKappa2;
	myRecoveryModel.NesC3 = myens->ensRecoveryModel.NesC3;
	myRecoveryModel.NesC4 = myens->ensRecoveryModel.NesC4;
	myRecoveryModel.NesC5 = myens->ensRecoveryModel.NesC5;

	myDragData.ZenerConsider = myens->ensDragData.ZenerConsider;
	myDragData.ZenerAlpha = myens->ensDragData.ZenerAlpha;
	myDragData.ZenerGamma = myens->ensDragData.ZenerGamma;

	myCADefMS.defmstype = myens->defmsmethod;

	onthefly_defragmentation = myens->onthefly_defragmentation;

	//MK::ideal orientations are only relevant at the ensembleHdl level at which I/O proceeds but reference possible with myens

	for ( uint32_t r = 0; r < myens->UserDefLogPoint_MS_Rendering.size(); r++ ) {
		rendering_atthisX.push_back ( myens->UserDefLogPoint_MS_Rendering[r] );
	}
	for ( uint32_t z = 0; z < myens->UserDefLogPoint_MS_RenderZPos.size(); z++ ) {
		rendering_atthisZPos.push_back( myens->UserDefLogPoint_MS_RenderZPos[z] );
	}
	for ( uint32_t d = 0; d < myens->UserDefLogPoint_X_CellListDefragment.size(); d++) {
		defragRXFront_atthisX.push_back ( myens->UserDefLogPoint_X_CellListDefragment[d] );
	}
	for ( uint32_t s = 0; s < myens->UserDefLogPoint_X_Output.size(); s++ ) {
		output_atthisX.push_back ( myens->UserDefLogPoint_X_Output[s] );
	}

	myrenderwindow.xmin = myens->renderwindow.xmin;
	myrenderwindow.xmax = myens->renderwindow.xmax;
	myrenderwindow.ymin = myens->renderwindow.ymin;
	myrenderwindow.ymax = myens->renderwindow.ymax;
	myrenderwindow.zmin = myens->renderwindow.zmin;
	myrenderwindow.zmax = myens->renderwindow.zmax;

	myrenderwindow.xmi = myens->renderwindow.xmi;
	myrenderwindow.xmx = myens->renderwindow.xmx;
	myrenderwindow.ymi = myens->renderwindow.ymi;
	myrenderwindow.ymx = myens->renderwindow.ymx;
	myrenderwindow.zmi = myens->renderwindow.zmi;
	myrenderwindow.zmx = myens->renderwindow.zmx;
}


void caHdl::init_processing( void )
{
	ensembleHdlP myens = this->myensHdl;
	tprocessingend = 0.0;

	for ( uint32_t t = 0; t < myens->ttime.size(); t++) {
		this->ttime.push_back ( myens->ttime[t] );
		this->ttemperature.push_back ( myens->ttemperature[t] );

		if ( ttime[t] >= tprocessingend ) {
			tprocessingend = ttime[t];
		}
	}
}


void caHdl::init_zenerdrag( void ) 
{
	ensembleHdlP myens = this->myensHdl;
	for ( uint32_t t = 0; t < myens->dispersoidtime.size(); t++) {
		this->zenertime.push_back ( myens->dispersoidtime[t] );
		this->zenerdispersion.push_back ( myens->dispersoidfr[t] );
	}
}


void caHdl::solve_INITIALIZATION( void )
{
	cout << "\t\tSolitary unit initializing..." << endl;
	init_cahdlprng();
	init_parameter();
	init_processing();
	init_zenerdrag();
}
























void caHdl::ompshar_initialize_region_myrxgpools( int regionid )
{
	//CALLED FROM WITHIN PARALLEL REGION
	//only call once for each region! from disjoint threadids regionid, 
	//such that memory for bookkeeping array is local
	uint32_t nrxg = myrxgpool.size();
	caregionMemHdlP theregion = regions[regionid];

	theregion->reg_myrxgp_counts = NULL;
	theregion->reg_myrxgp_counts = new uint32_t[nrxg];
	QUICKASSERT ( theregion->reg_myrxgp_counts != NULL );
	theregion->reg_nmyrxgp_counts = nrxg;

	for ( uint32_t rg = 0; rg < nrxg; rg++ ) {
		theregion->reg_myrxgp_counts[rg] = 0;
	}
}


void caHdl::ompshar_initialize_region_mydfgpools( int regionid )
{
	//CALLED FROM WITHIN PARALLEL REGION
	//only call once per region from thread regionid, such that memory for bookkeeping array is local
	uint32_t ndg = mydefgpool.size();
	caregionMemHdlP theregion = regions[regionid];

	theregion->reg_mydfgp_counts = NULL;
	theregion->reg_mydfgp_counts = new uint32_t[ndg];
	QUICKASSERT ( theregion->reg_mydfgp_counts != NULL );
	theregion->reg_nmydfgp_counts = ndg;

	for ( uint32_t dg = 0; dg < ndg; dg++ ) {
		theregion->reg_mydfgp_counts[dg] = 0;
	}
}


void caHdl::synchronize_mydefpool_counts( void )
{
	//CALLED FROM WITHIN PARALLEL REGION
	//##MK::clear current content of the mydefgpool.cellcount
	uint32_t ndgp = this->mydefgpool.size();
	for ( uint32_t dg = 0; dg < ndgp; ++dg ) {
		mydefgpool[dg].cellcount = 0;
	}


	//now collect cell counts from all regions
	for ( uint32_t nreg = 0; nreg < regions.size(); nreg++ ) { //not omp_get_num_threads() as it evaluates to 1 outside a parallel region!
		caregionMemHdlP theregion = regions[nreg];
		uint32_t nreg_dgp = theregion->reg_nmydfgp_counts;

		QUICKASSERT ( nreg_dgp == mydefgpool.size() ); //##DEBUG

		for ( uint32_t d = 0; d < nreg_dgp; ++d ) {
			mydefgpool[d].cellcount += theregion->reg_mydfgp_counts[d];
		}
	}
}


void caHdl::synchronize_myrxgpool_counts( void )
{
	//CALLED FROM WITHIN PARALLEL REGION
	//##MK::clear current content of the myrxgpool.cellcount
	uint32_t nrxgp = this->myrxgpool.size();
	for ( uint32_t rg = 0; rg < nrxgp; ++rg ) {
		myrxgpool[rg].cellcount = 0;
	}


	//now collect cell counts from all regions
	for ( uint32_t nreg = 0; nreg < regions.size(); nreg++ ) {
		caregionMemHdlP theregion = regions[nreg];
		uint32_t nreg_rxgp = theregion->reg_nmyrxgp_counts; //reg_myrxgp_counts.size();

		QUICKASSERT ( nreg_rxgp == myrxgpool.size() ); //##DEBUG

		for ( uint32_t r = 0; r < nreg_rxgp; ++r ) {
			myrxgpool[r].cellcount += theregion->reg_myrxgp_counts[r]; //reg_myrxgp_counts[r]; //arrays are synchronized...
		}
	}
}


void caHdl::determine_initial_nucleitexture_myrxgpoolids( void )
{
	//cache deformation structure
	uint32_t* nucmicrotex = NULL;
	uint32_t nnucmicrotex = myensHdl->standardlagen.size() + 1;
	nucmicrotex = new uint32_t[nnucmicrotex];
	QUICKASSERT ( nucmicrotex != NULL );

	for (uint32_t i = 0; i < nnucmicrotex; i++) {
		nucmicrotex[i] = 0;
	}

	//scan all RX nuclei
	uint32_t ocategory;
	uint32_t nrxgp = myrxgpool.size();
	for ( uint32_t nuc = 0; nuc < nrxgp; nuc++) {
		ocategory = myoripool[myrxgpool[nuc].caori].closestideal;
		nucmicrotex[ocategory] = nucmicrotex[ocategory] + 1;
	}

//##DEBUGfor ( uint32_t oc = 0; oc < nnucmicrotex; oc++) {cout << "NucleationTexture component " << oc << "\t\t" << nucmicrotex[oc] << endl;}

	NucleiMicrotextureClasses = nnucmicrotex;
	NucleiMicrotexture = nucmicrotex;
}


void caHdl::solve_INITIALIZATION_THREADED_DOMAIN_DECOMPOSITION( void )
{
	//MK::this function exemplifies how to initialize in parallel the CA region construct and an associated deformed microstructure representation
	//MK::first threads need initalization and allocation and first-touch of NUMA-local memory, then critical functions can be executed after which
	//start parallel use all OMP_NUM_THREADS, if work should only be done by one thread
	//MK::the CONVENTION is this is the MASTER THREAD to keep compatibility with MPI libraries that require master threads for FUNNELED hybrid operation
	#pragma omp parallel	//by default now all variables are shared among the threads in the team
	{
		#pragma omp master
		{
			//cout << "\t\tInitializing OMP domain decomposition..." << endl;
			for ( int t = 0; t < omp_get_num_threads(); t++ ) {
				regions.push_back ( NULL );
				cellStatii.push_back( NULL );
			}
		}
		#pragma omp barrier //necessary, team relies on this ca decomposition

		ompshar_init_regions( omp_get_thread_num() );

		ompshar_ca2regions();
		#pragma omp barrier
		//barrier with implicit flush necessary as all threads need to have the limits established and stored before other threads can read them out safely //#pragma omp flush
		//also the halos have to be defined before they can be identified by neighboring threads

		ompshar_init_regionlimits( omp_get_thread_num() ); //executed in parallel, thread ids are disjoint and hence initialize limits in disjoint regions!

		ompshar_link_haloregions( omp_get_thread_num() ); //serves to make each region now where it can find the pointers to the haloregions in the adjacent regions

		#pragma omp master
		{
			cout << "myRank " << this->myensHdl->myRank << "; JobID " << this->jobid << "; myRegion " << omp_get_thread_num() << " voxelization was successful." << endl;
			caregionMemHdlP re = regions[omp_get_thread_num()];
			cout << "\t\tLimits of myRegion " << omp_get_thread_num() << " are xmimx/ymimx/zmimx = " << re->myGeom.nreg_rdmin << ";" << re->myGeom.nreg_rdmax << ";" << re->myGeom.nreg_tdmin << ";" << re->myGeom.nreg_tdmax << ";" << re->myGeom.nreg_ndmin << ";" << re->myGeom.nreg_ndmax << endl;
			cout << "\t\tLimits of myRegion " << omp_get_thread_num() << " are x/y/z/xy/xyz/cellsize/xglo/yglo/zglo = " << re->myGeom.nreg_rd << ";" << re->myGeom.nreg_td << ";" << re->myGeom.nreg_nd << ";" << re->myGeom.nregarea_rdtd << ";" << re->myGeom.nregvol_rdtdnd << ";" << re->myGeom.cellsize << ";" << re->myGeom.nedge_global_rd << ";" << re->myGeom.nedge_global_td << ";" << re->myGeom.nedge_global_nd << endl;
			/*cout << "\t\tNeighbors of myRegion " << omp_get_thread_num() << " are = ";
			for ( uint32_t r = 0; r < (NUMBER_OF_NEIGHBORS+MYSELF); r++ ) {
				cout << re->mynborregions[(r*IDS_AND_LIMITS)+THE_IDS] << "__" << re->mynborregions[(r*IDS_AND_LIMITS)+THE_XMIN] << "__" << re->mynborregions[(r*IDS_AND_LIMITS)+THE_XMAX] << "__" << re->mynborregions[(r*IDS_AND_LIMITS)+THE_YMIN] << "__" << re->mynborregions[(r*IDS_AND_LIMITS)+THE_YMAX] << "__" << re->mynborregions[(r*IDS_AND_LIMITS)+THE_ZMIN] << "__" << re->mynborregions[(r*IDS_AND_LIMITS)+THE_ZMAX] << ";";
			} cout << endl;*/
		}
	} //implicit barrier
}


void caHdl::solve_SYNTHESIZED_DEFORMEDSTRUCTURE_CUBOIDBLOCKS( void )
{
	#pragma omp parallel
	{
		#pragma omp master
		{
			determine_polycrystalline_geometry();
			pick_deformedgrains();
		}
		#pragma omp barrier //necessary because pragma omp construct has no implied barrier

		ompshar_voxelize_defms();

		regions[omp_get_thread_num()]->reg_nmydefgpool = mydefgpool.size();

		//##MK::add further functionality here

	} //implicit barrier
}


void caHdl::solve_SYNTHESIZED_DEFORMEDSTRUCTURE_FROMEBSDDATA( void )
{
	/*ebsd2polycrystal();

	//mycell memory was already allocated so now we can utilize an OpenMP parallel CA to voxelize the synthetic structure from the seeds with unit speed
	grow_deformedgrains_voxelize_omp();
	*/

	ebsd2polycrystal2(); //##MK::only single-threaded

	for ( uint32_t mr = 0; mr < regions.size(); mr++ ) {
		regions[mr]->reg_nmydefgpool = mydefgpool.size();
	}

	//cout << "Thread " << i << "-->" << regions[i]->reg_nmydefgpool << endl;
	//##MK::add further functionality here
}


void caHdl::solve_SYNTHESIZE_DEFORMATIONSTRUCTURE( void )
{
	//init log file
	if ( renderingForThisCA == true && outopt_localrenderhow == RENDERING_MS3D ) {
		if ( outopt_localrenderfileformat == RENDERING_FILEFORMAT_HDF5 && (outopt_localrendercolor == RENDERING_COLOR_GRAINID || outopt_localrendercolor == RENDERING_COLOR_IPFZ || outopt_localrendercolor == RENDERING_COLOR_RHO) ) {
			string h5fname = "SCORE." + std::to_string( this->myensHdl->simid ) + "." + "VoxelData.h5";
			bool h5status = true;
			//omp_set_lock(&h5lock);
			h5status = hdf5_write_coordinates( h5fname, "/coordinates" );
			//omp_unset_lock(&h5lock);
		}
	}

	if ( myCADefMS.defmstype == CUBOID_DEFMS ) {
		solve_SYNTHESIZED_DEFORMEDSTRUCTURE_CUBOIDBLOCKS();
	}
	else if ( 	myCADefMS.defmstype == SEMEBSD_2D3D_COLUMNSINGLE ||
				myCADefMS.defmstype == SEMEBSD_2D3D_COLUMNSTACK ||
				myCADefMS.defmstype == SEMEBSD_2D3D_COLUMNSHIFT  ) {
		solve_SYNTHESIZED_DEFORMEDSTRUCTURE_FROMEBSDDATA();
	}
	else {
		cout << "ERROR::Unable to detect the deformation microstructure synthesis method!" << endl;
	}

	nmydefgpool = mydefgpool.size();

//##DEBUGcout << "mydefgpool.size();nmydefgpool = " << mydefgpool.size() << "\t\t" << nmydefgpool << endl;
}


void caHdl::solve_DETECT_GRAINBOUNDARIES( void ) 
{
	double timer = MPI_Wtime();
	bool stillgood = true;

	#pragma omp parallel
	{
		double myt0 = omp_get_wtime();

		//each thread has distinct memory region he works on
		int tid = omp_get_thread_num();
		int nthreads = omp_get_num_threads();

		regions[tid]->myjunctions.mycaHdl = this;
		regions[tid]->myjunctions.mythreadid = tid; //thread takes care of statically scheduled domain part...
		regions[tid]->myjunctions.mymemreg = this->regions[tid];
		regions[tid]->myjunctions.juncMemGuard = 0.0;

		regions[tid]->myjunctions.ompshar_init_junctions();
		regions[tid]->myjunctions.ompshar_detect_junctions();
		#pragma omp barrier

		#pragma omp master
		{
			//check if all buckets have the same length
			stillgood = true;
			unsigned int mastersize = regions[MASTER]->myjunctions.fjuncSkeletonSize;
			for ( int t = (MASTER+1); t < nthreads; t++ ) {
				if ( regions[t]->myjunctions.fjuncSkeletonSize != mastersize ) {
					stillgood = false;
					break;
				}
			} //by entry into parallel region the variable is shared by default, no race however as we are in a omp master construct!
		}
		#pragma omp barrier //necessary! as there is no implicit one after a omp master construct...

		/*if ( this->outopt_logboundaries == OUTPUT_LOGBND_YES ) { //##MK::output all boundaries
			write_junction_log_gb( tid, false );
		}*/

		if ( stillgood == true ) { //synchronize bookkeeping information over all junctions
			//##MK::get overview of junctions
			regions[omp_get_thread_num()]->myjunctions.ompshar_report_junctions();

			regions[omp_get_thread_num()]->myjunctions.ompshar_consolidate_junctions();
			#pragma omp barrier //necessary because we must not delete memory which other thread may still need to synchronize...!

			regions[omp_get_thread_num()]->myjunctions.ompshar_calc_junction_properties();
			regions[omp_get_thread_num()]->myjunctions.ompshar_cleanup_foreign_junctions();

			//consecutive labeling of consolidated boundaries
			#pragma omp master
			{
				uint32_t consecutiveGBID = 0;
				uint32_t nfj = this->regions[MASTER]->myjunctions.fjuncSkeletonSize;
				uint32_t nthr = nthreads;
				for ( uint32_t pmax = 0; pmax < nfj; pmax++ ) { //fjuncSize the same on all threads after consolidation
					uint32_t thr = workPartitioning(pmax, nthr); //NEEDS TO BE THE SAME STATIC WORK PARTITIONING AS IN CONSOLIDATE JUNCTIONS
					uint32_t nf = this->regions[thr]->myjunctions.fjuncSkeleton[pmax].len;
					for ( uint32_t f = 0; f < nf; f++ ) {
						this->regions[thr]->myjunctions.fjuncSkeleton[pmax].thefaces[f].gbID = consecutiveGBID;
						consecutiveGBID++;
					}
				}
			} //progress prompting can continue in parallel
			#pragma omp critical
			{
				std::cout << "\t\tThread ID " << omp_get_thread_num() << " junction detection " << (omp_get_wtime() - myt0) << std::endl;
			}
		}
	} //implicit barrier at end of parallel region

	if ( stillgood == false ) { std::cout << "ERROR::Junction bucket inconsistency!" << std::endl; return; }

	if ( renderingForThisCA == true ) {
		if ( outopt_localrenderboundaries == RENDERING_BOUNDARIES_YES ) {
			write_junction_skeleton( OUTPUT_JUNCTIONS_GBTHREADCOLORING );
			write_junction_skeleton( OUTPUT_JUNCTIONS_GBFACES );
			write_junction_skeleton( OUTPUT_JUNCTIONS_TRIJUNCTIONS );
			write_junction_skeleton( OUTPUT_JUNCTIONS_HIGHERORDER );
		}
	}

	if ( this->outopt_logboundaries == OUTPUT_LOGBND_YES ) { //##MK::output all boundaries
		write_junction_log_gb( MASTER, true );
		write_junction_log_tj();
		write_junction_log_hj();
	}

	cout << "Grain boundary detection in " << ( (double) MPI_Wtime() - timer ) << " seconds!" << endl;
}


//////////////NUCLEATION MODELING at TRACKED DEFORMED GRAIN BOUNDARIES UTILIZING AUTOMATON WITH SEEDIDS

void caHdl::solve_GRAINBOUNDARYNUCLEATION( void )
{
	//##MK in M. Kuehbach's Phd thesis the specific implementation of grain boundary nucleated reactions which were exercised within 
	//C. Haase's HighMn SFB761 Twip steel project aided to model site-saturated grain boundary nucleated RX processes
	//so it sufficies to populate the myrxgpool vector with picking physically the desired per se disjoint boundary cells,
	//at which nuclei should be generated, then solve_NUCLEATIONMODELING must not be executed 
	//instead solve_RXGROWTH places the nuclei without prior knowledge that the mycellgrid assignment came from a DEFMS automaton

	if( this->myNucleationModel.gbnucleation == GBNUCLEATION_YES_PICKRND ) {
		solve_nucmodeling_gbnuc_pickrandomly();
		return;
	}

	if ( this->myNucleationModel.gbnucleation == GBNUCLEATION_YES_PHYS_LAGBHAGB || this->myNucleationModel.gbnucleation == GBNUCLEATION_YES_PHYS_ONLYHAGB ) {
		//solve_nucmodeling_gbnuc_physicalmodelFast();
		cout << "ERROR::For grain boundary nucleation at the currently only option 3 is implemented!" << endl;
		return;
	}

	//##MK::we could empty the boundaries at that point, otherwise the remain until destruction of the this caHdl class object
}




uint32_t caHdl::getDeformedGrainID( uint32_t xglo, uint32_t yglo, uint32_t zglo )
{
	//find in which memory region the global position x,y,z is and read tmpdefgseed ID
	//get region limits
	uint32_t xmi, xmx, ymi, ymx, zmi, zmx, xx, xxyy; //inclusive CA global bounds in the memory region
	for ( uint32_t mr = 0; mr < this->regions.size(); mr++ ) {
		xmi = this->regions[mr]->myGeom.nreg_rdmin;
		xmx = this->regions[mr]->myGeom.nreg_rdmax;
		ymi = this->regions[mr]->myGeom.nreg_tdmin;
		ymx = this->regions[mr]->myGeom.nreg_tdmax;
		zmi = this->regions[mr]->myGeom.nreg_ndmin;
		zmx = this->regions[mr]->myGeom.nreg_ndmax;
		xx = this->regions[mr]->myGeom.nreg_rd;
		xxyy = this->regions[mr]->myGeom.nregarea_rdtd;
		if ( xglo < xmi ) continue;
		if ( xglo > xmx ) continue;
		if ( yglo < ymi ) continue;
		if ( yglo > ymx ) continue;
		if ( zglo < zmi ) continue;
		if ( zglo > zmx ) continue;
		//not continued we found the domain point can only be in one memoryregions as these are nonoverlapping but SU domain filling
		return regions[mr]->mycellgrid[xglo-xmi+(yglo-ymi)*xx+(zglo-zmi)*xxyy];
	}

	return 0;
}


void caHdl::solve_NUCLEATIONMODELING( void ) 
{
	if ( myNucleationModel.gbnucleation != GBNUCLEATION_NO ) {
		solve_GRAINBOUNDARYNUCLEATION();
		return;
	}

	if ( myNucleationModel.csrnucleation == CSRNUCLEATION_YES_ENFORCED ) {
		solve_nucmodeling_csr_enforced();
		return;
	}

	if ( myNucleationModel.csrnucleation == CSRNUCLEATION_YES_PICKRND ) {
		solve_nucmodeling_csr_pickrandomly();
		return;
	}

	//if ( myNucleationModel.csrnucleation == CSRNUCLEATION_YES_DIEHL
	if ( 	myNucleationModel.csrnucleation == CSRNUCLEATION_DIEHL_RANDOMSO3 ||
			myNucleationModel.csrnucleation == CSRNUCLEATION_DIEHL_RANDOMDEFORI	||
			myNucleationModel.csrnucleation == CSRNUCLEATION_DIEHL_SCATTEREXISTENT  ) {
		solve_nucmodeling_csrdiehl();
	}

	if ( myNucleationModel.clustnucleation == CLUSTERNUCLEATION_YES_PICKRND ) {
		//##MK::verify for use with OpenMP
		solve_nucmodeling_ellipsoidalcluster_pickrandomly();
		return;
	}
}


void caHdl::solve_REPLACE_CA_STRUCTURE( void ) 
{
	//now we know all potential nuclei, so no new orientations are added further because only existing ones grow
	nmyrxgpool = myrxgpool.size();

	//MK::allocate memory to store the deformation texture, only once and sequential
	DefMicrotextureClasses = myensHdl->standardlagen.size() + 1;
	DefMicrotexture = NULL;
	DefMicrotexture = new uint32_t[DefMicrotextureClasses];
	QUICKASSERT ( DefMicrotexture != NULL );
	for ( uint32_t d = 0; d < DefMicrotextureClasses; d++ ) { this->DefMicrotexture[d] = 0; }

	int threadid;
	uint32_t c, nxyz, seedid, mydgid;
	double reg_storedenergy, scaler;

	#pragma omp parallel private(threadid, c, nxyz, seedid, mydgid, reg_storedenergy, scaler)
	{
		threadid = omp_get_thread_num(); //MK::hereby it is dictacted that each thread works on disjoint data!
		nxyz = regions[threadid]->myGeom.nregvol_rdtdnd;

		//threads initialize local copy of the myrxgpool and the mydefgpoolcell counters
		//such that during updating the cells in the growth simulation, the bookkeeping is race free and thread-locality maximized
		ompshar_initialize_region_myrxgpools(threadid);
		ompshar_initialize_region_mydfgpools(threadid);

		//cache locally the deformation structure
		uint32_t* reg_defmicrotexture = NULL;
		uint32_t ndefmicrotexture = myensHdl->standardlagen.size() + 1; //required the same for all threads!
		reg_defmicrotexture = new uint32_t[ndefmicrotexture];
		for (uint32_t i = 0; i < ndefmicrotexture; i++) { reg_defmicrotexture[i] = 0; }

		//MK::still the mycellgrid is populated with seedids, additionally for the growth of the RX nuclei the boundaries among the deformed grains 
		//are no longer of interest, so replace mycellgrid from indices referencing entries in tmpdefgseeds to indices referencing entries in mydefgpool
		//first replace the CA structure, then bookkeep
		for ( c = 0; c < nxyz; c++ ) {
			seedid = regions[threadid]->mycellgrid[c];
			mydgid = tmpdefgseeds[seedid].mydefgpoolid;
			regions[threadid]->mycellgrid[c] = mydgid;
//##DEBUGcout << "c;seedid;mydgid;grid\t\t" << c << ";" << seedid << ";" << mydgid << ";" << regions[threadid]->mycellgrid[c] << endl;
		}

		//now  the cells know to which defgIDs they were assigned, bookkeep texture over regions
		reg_storedenergy = 0.0;

		uint32_t cadefgid;
		uint32_t ocategory;
		scaler = (1.0 / SCALING_DEFAULT_DISDENS);
		for ( c = 0; c < nxyz; c++ ) {
			cadefgid = regions[threadid]->mycellgrid[c];

			regions[threadid]->reg_mydfgp_counts[cadefgid] += 1;

			reg_storedenergy = reg_storedenergy + (mydefgpool[cadefgid].rho0 * scaler);
		}

		//beneficially, as compared to MPI only version now instead of nxyz cache missed to get the ocategory, much fewer accesses to the myoripool
		uint32_t ndgid = regions[threadid]->reg_nmydfgp_counts;
		for ( cadefgid = 0; cadefgid < ndgid; cadefgid++ ) {
			ocategory = myoripool[mydefgpool[cadefgid].caori].closestideal; //because random is 0
			reg_defmicrotexture[ocategory] = reg_defmicrotexture[ocategory] + regions[threadid]->reg_mydfgp_counts[cadefgid];
		}

		//out of the perspective of each thread reg_defmicrotexture was locally updated for all ocategories, but now
		//threads need to reduce information back to the caHdl...
		#pragma omp critical
		{
			StoredEnergy = StoredEnergy + reg_storedenergy;
			for ( uint32_t d = 0; d < DefMicrotextureClasses; d++ ) {
				DefMicrotexture[d] = DefMicrotexture[d] + reg_defmicrotexture[d];
			}
		}

		//thread local deformation texture no longer of interest
		delete [] reg_defmicrotexture;
		reg_defmicrotexture = NULL;

	} //implicit barrier and flush

	//counter pool bookkeeping, must not be executed in parallel otherwise race conditions!
	synchronize_mydefpool_counts();
	synchronize_myrxgpool_counts();

	//nuclei bookkeeping
	determine_initial_nucleitexture_myrxgpoolids();

	//clean-up tmdefgseeds as it is not longer needed
	delete [] tmpdefgseeds;
	tmpdefgseeds = NULL;
	myMemGuard = myMemGuard - ( ndefgseeds * sizeof(struct defgseed) );

	if (myensHdl->myRank == MASTER ) { cout << "myRank " << this->myensHdl->myRank << "; JobID " << this->jobid << " CA structure replaced successfully" << endl; }
}


void caHdl::cleanMyGrainEvolution( void )
{
	//because caHdl does not know about dynamic allocation of memory from within the struct
	for ( uint32_t s = 0; s < mygrainevolution.size(); s++) {
		delete [] mygrainevolution[s].localdatabucket;
		mygrainevolution[s].localdatabucket = NULL;
		myMemGuard = myMemGuard - ( (double) mygrainevolution[s].nlocaldata * sizeof(uint32_t) );
	}
	//vector vector does allocate self
	
	for ( uint32_t s = 0; s < myrxprofile.size(); s++ ) {
		delete [] myrxprofile[s].localdatabucket;
		myrxprofile[s].localdatabucket = NULL;
	}
}


void caHdl::cleanBookkeeping( void )
{
	delete [] DefMicrotexture;
	DefMicrotexture = NULL;

	delete [] NucleiMicrotexture;
	NucleiMicrotexture = NULL;

//#pragma omp parallel 
//{
	for ( uint32_t r = 0; r < regions.size(); r++ ) {
		//if ( omp_get_thread_num() == r ) {
			if ( regions.at(r) != NULL ) {
				delete regions[r];
				regions[r] = NULL;
			}
		//}
			if ( cellStatii.at(r) != NULL ) {
				delete cellStatii[r];
				cellStatii[r] = NULL;
			}
	}
//}

/*	regionsuint32_t nregions = regions.size();
	for ( uint32_t t = 0; t < nregions; t++ ) {
		delete [] regions[t];
	}
*/
}


void caHdl::solve_RXGROWTH( void )
{
	log_initialization();

	unsigned char status = OMP_OPERATIVE;
	uint32_t integrationstep = 0;
	double myomptimer = 0.0;
	t = 0.0;
	X = 0.0;
	Xcells = 0.0;
	dXstep = 0.0;
	Sv = 0;
	step = 0;

	//MK::you desire initial I/O prior to simulation? --> inject the code then here

	update_mydgoptimize_rho();
	update_mydgoptimize_cnt();
	update_myrhomax();
	log_initrho_scalebar();

	update_temperature();			//not OpenMP parallelized because writes to caHdl and of low op intensity
	update_atomisticproperties();
	//update_microchemistry();
	update_intrinsicmobilities();
	//update_deformedsubstructure();

#ifdef REPORTSTYLE_USER
	if ( this->myensRank == MASTER ) 
		cout << this->jobid << "\tstep;t;X;dt;CurrTemperature;Sv\t\t" << this->step << "\t\t" << this->t << "\t\t" << this->X  << "\t\t" << this->dt << "\t\t" << this->CurrentTemperature << "\t\t" << this->Sv << endl;
#endif

	//get initial delta time for the numerical integration scheme
	dt = get_dtmax_minimum();

	//mechanism to enable the writing of HDF5 files (though sequentially) from within a threaded region
	omp_init_lock(&(this->h5lock)); //only the calling thread who sets a lock is allowed to unlock it again, so outside a parallel region locker is the master thread

	#pragma omp parallel private( integrationstep, status, myomptimer )
	{
		#pragma omp master
		{
			ompcrit_clarify_status_infectedcells( RXFRACTION_THRESHOLD ); //caching of partially infected cells
		}
		//barrier necessary 
		#pragma omp barrier

		//initial output
		if ( outopt_localrenderhow == RENDERING_MS2D ) {
			if ( outopt_localrenderfileformat == RENDERING_FILEFORMAT_RAW ) {
				if ( outopt_localrendercolor == RENDERING_COLOR_IPFZ ) {
					ompcrit_write_zsection_looping_ipfz( 0.0, COLORIZE_RX_LEAVE_DEFORMED_BLACK );
					ompcrit_write_zsection_looping_ipfz( 0.0, COLORIZE_DEF_LEAVE_RX_BLACK );
					ompcrit_write_zsection_looping_ipfz( 0.0, COLORIZE_RX_RHOGREYSCALE_DEFORMED );
				}
				else
					cout << "ERROR::2D output in RAW format supports currently only IPFZ coloring!" << endl;
			}
			if ( outopt_localrenderfileformat == RENDERING_FILEFORMAT_HDF5 ) 
				cout << "ERROR::2D output of complete 2D structure with HDF5 currently not supported!" << endl;
		}

		//##MK:for study with Martin Diehl, MMM2018 Proceedings Paper
		//switched off initial ang file
		/*if ( outopt_localartsemebsd == OUTPUT_SEMEBSD_YES ) {
			//##MK::requires ompcrit_clarify_status_infectedcells to have been executed
			ompcrit_write_semebsd_looping( DEFAULT_SEMEBSD_MODE );
		}*/
		#pragma omp barrier

		integrationstep = 0;
		status = OMP_OPERATIVE;

		QUICKASSERT ( regions.size() == omp_get_num_threads() );

		regions[omp_get_thread_num()]->ompshar_init_myRXFrontInside();
		regions[omp_get_thread_num()]->ompshar_init_myRXFrontBorder();

		//MK::a barrier is necessary!, no placing of nuclei before RX front cell management structures have been initialized
		#pragma omp barrier

		//add nuclei into the simulation, at the moment ##MK::MPI/OMP only site-saturated --> inject time-dependent nucleation here...
		vector<rxgcommit>* thrsafesync = NULL;	thrsafesync = new vector<rxgcommit>;
		ompshar_sim_myCA_sitesaturatedNucleation( thrsafesync );

		regions[omp_get_thread_num()]->ompshar_init_serialsectioning();

		#pragma omp barrier

		#pragma omp critical
		{
			for ( uint32_t s = 0; s < thrsafesync->size(); s++ ) {
				myrxgpool[thrsafesync->at(s).nrxgid].nucsite = thrsafesync->at(s).nucsite;
				myrxgpool[thrsafesync->at(s).nrxgid].tincub = thrsafesync->at(s).tincub;
			}
		}
		delete thrsafesync; thrsafesync = NULL;

		//necessary, not possible to start before the threads seeded all nuclei surplus barrier implicitly flushes
		#pragma omp barrier

		#pragma omp master
		{
			for (uint32_t r = 0; r < regions.size(); r++) {
				if ( regions[r]->reg_myMobilityWeightMax >= myMobilityWeightMax ) { myMobilityWeightMax = regions[r]->reg_myMobilityWeightMax; }
				nmynuclei = nmynuclei + regions[r]->reg_nmynuclei;
			}
			//cout << "Starting RX Growth with myMobilityWeightMax = " << myMobilityWeightMax << " and " << nmynuclei << " nuclei." << endl;
		}
		//barrier necessary, as otherwise the threads would start iterating even though myMobilityWeightMax is not yet consistent because pragma omp master construct has no implicit barrier
		#pragma omp barrier

		#pragma omp critical
		{
			cout << "Thread " << omp_get_thread_num() << " enters main loop" << endl;
		}

		//MAIN LOOP needed transformation as do/while not allowed in parallel region do {//write the correct value into the loop
		for ( integrationstep = 0; integrationstep < DEFAULT_NMAX; integrationstep++ ) {

			status = OMP_OPERATIVE;
			//each thread reads globally synchronized variable and checks whether or not there is still something to do, because we are not allowed to break out uncontrolled of a parallel for construct
			if ( (Xcells >= myCAGeometry.nboxvol_rdtdnd) || ( X > XMAX ) || (integrationstep > NTSTEPSMAX) || (t > tprocessingend) || (t > TMAX) ) {
				status = OMP_COMPLETED;
			}

			if ( status == OMP_OPERATIVE ) {
				//to assure the myRXFront list to be as compact filled with ACTIVE cells between [0, nextRX...) as possible
				if ( onthefly_defragmentation == true ) {
					if ( (Xcells >= defragRXFront_atthisX[regions[omp_get_thread_num()]->reg_defrag_cnt]) && (regions[omp_get_thread_num()]->reg_defrag_cnt < defragRXFront_atthisX.size()) ) {
						regions[omp_get_thread_num()]->ompshar_defragmentation();
						regions[omp_get_thread_num()]->reg_defrag_cnt++;
					}
				}

				//MK::no time dependent nucleation model for OpenMP at the moment--> if desired inject it here!

				//fill ACTIVE cells while keeping track of INACTIVE cells in myRecyclingList, totally local and independent of others
				regions[omp_get_thread_num()]->ompshar_sim_myCA_calcGrowthStepInside();
				regions[omp_get_thread_num()]->ompshar_sim_myCA_calcGrowthStepBorder();

				//let the fully recrystallized ACTIVE cells now infect their neighbors
				//MK::inside is independent within the cellregion of one thread requiring no synchronization, as inside cells cannot in Moore neighborhood protrude beyond local mycellgrids at most only into local border
				//in this function cells on myRXFront become activated again, but during updating in one thread all interesting positions NOT ADDRESSES! are already known in the FullRX Lists
				//so it cannot happen that reactivated cells infect out of cycle...
				regions[omp_get_thread_num()]->ompshar_sim_myCA_updateFullRXInside();

				//possible infections outside a region are stored in the halos or rejected directly by the halos
				//thus, the halos in turn keep track of whether volume in a neighboring region is still infectable
				regions[omp_get_thread_num()]->ompshar_clear_halo_bookkeeping();
				regions[omp_get_thread_num()]->ompshar_sim_myCA_updateFullRXBorder();
				#pragma omp barrier
				//MK::necessary! because all threads should finish populating the halo regions before we can exchange halo pieces of information
				//otherwise we require overhead to ensure a conflictless overlap of reading out halos for pairs of already finished threads and those not,
				//very likely the overhead is in the order of the savings...
				//instead, after having enforced the barrier the threads can read out the halos in parallel and synchronize
				//they are guided by neighboring regions' haloref list which stores which positions were newly infected by the neighboring thread in each halo
				regions[omp_get_thread_num()]->ompshar_synchronize_haloregions_rx();

				#pragma omp critical
				{
					cout << "Thread " << omp_get_thread_num() << " about to perform serial sectioning" << endl;
				}

				if ( renderingForThisCA == true ) {

					regions[omp_get_thread_num()]->ompshar_sim_myCA_clarify_status( RXFRACTION_THRESHOLD );
					#pragma omp critical
					{
						cout << "Thread " << omp_get_thread_num() << " clarified status" << endl;
					}
					regions[omp_get_thread_num()]->ompshar_sectioning_analysis_rx();
					#pragma omp critical
					{
						cout << "Thread " << omp_get_thread_num() << " performed threadlocal serial sectioning" << endl;
					}
				}

				#pragma omp barrier

				if ( outopt_localthreadprof == OUTPUT_THREADPROFILING_YES ) {
					regions[omp_get_thread_num()]->omp_log_rxfrontstats();
					#pragma omp critical
					{
						cout << "Thread " << omp_get_thread_num() << " profiled thread" << endl;
					}
				}
			} //growth step completed...
			myomptimer = omp_get_wtime();

			//...now all synchronization work on the master as writing out global variables in memory shared by all threads
			#pragma omp master
			{
				if ( status == OMP_OPERATIVE ) {
					dXstep = 0.0;
					Sv = 0;
					double mymobmax = DEFAULT_PMAX;
					for ( uint32_t tid = 0; tid < regions.size(); tid++ ) { //accumulate local results
						if ( regions[tid]->reg_myMobilityWeightMax >= mymobmax ) {
							mymobmax = regions[tid]->reg_myMobilityWeightMax;
						}
						dXstep = dXstep + ( regions[tid]->reg_dXstepInside + regions[tid]->reg_dXstepBorder );
						Sv = Sv + (regions[tid]->reg_SvInside + regions[tid]->reg_SvBorder);
						Xcells = Xcells + ( (double) (regions[tid]->reg_nUpdatedCellsInside + regions[tid]->reg_nUpdatedCellsBorder) ); 
					}
					myMobilityWeightMax = mymobmax;
					X = (Xcells / ((double) myCAGeometry.nboxvol_rdtdnd));

					synchronize_myrxgpool_counts();
					synchronize_mydefpool_counts();

					if ( myRecoveryModel.RecoveryConsider != RECOVERY_NO ) {
						update_mydgoptimize_rho(); //conservative checking of maximum rho left, ##MK::once implementing recovery model utilize additionally
					}
					update_mydgoptimize_cnt(); //dont call prior to having syncrhonized cnts!
					update_myrhomax();

					cout << "Thread MASTER about to enter log_OUTPUT_Diehl()" << endl;

					//log_OUTPUT();
					log_OUTPUT_Diehl();

					log_PERCOLATION();

					//--> up to here it is possible to account always a new for the myMobMax in order to stretch the time towards the end of the simulation in a 
					//scenario where only slowly migrating boundaries remain and the estimate of the maximum Pvalue that was EVER AT SOME POINT obtained is too conservative
					//prepare for next timestep, same integration order than START1D/3DCA and COReV3

					t += dt; 

					update_temperature();
					update_atomisticproperties();
					update_microchemistry();
					update_intrinsicmobilities();
					update_deformedsubstructure(); //can reduce rho thus after this recovery call the actual maximum rho in the system can be less than myrhomax, however assuming myrhomax is larger, thus more driving force available we are conservative and assume that the time increment has to be slightly shorter as in reality it could...

					//updateFullRX was the last one to determined Pmax in this timestep, now consider potentially higher mobility to refine time stepping
					dt = get_dtmax_minimum();

					step++;

//#ifdef REPORTSTYLE_USER
					if ( this->myensRank == MASTER ) 
						cout << "\tstep;t;X;dt;CurrTemperature;Sv\t\t" << this->step << "\t\t" << this->t << "\t\t" << this->X  << "\t\t" << this->dt << "\t\t" << this->CurrentTemperature << "\t\t" << this->Sv << endl;
//#endif
				}
			}
			#pragma omp barrier //##MOVE this barrier behind the loop all threads identify the state on shared variables all values of which are the same for each thread

			regions[omp_get_thread_num()]->profiling_growthmachine[regions[omp_get_thread_num()]->profiling_growthmachine.size()-1].tSeqOverhead = omp_get_wtime() - myomptimer;
		}
		//OMP_COMPLETED

		#pragma omp master 
		{
			tsimend = t;
			stepsimend = step;

			//final output
			if ( this->renderingForThisCA == true ) {
				if ( outopt_localrenderfileformat == RENDERING_FILEFORMAT_RAW ) { //final structure as raw
					cout << "ERROR::Output of final structure currently supported only with HDF5, chose instead 1.0 as an additional output log point in the input file!" << endl;
				}
				if ( outopt_localrenderfileformat == RENDERING_FILEFORMAT_HDF5 ) { //final structure and finalize XDMF file
				
					if ( outopt_localrendercolor == RENDERING_COLOR_GRAINID ) 
						ompcrit_write_voxeldata_h5( RENDERING_COLOR_GRAINID, "/GROWTH3DGrainID", true, "dummy", this->t, this->X, this->step );
					else if ( outopt_localrendercolor == RENDERING_COLOR_IPFZ )
						ompcrit_write_voxeldata_h5( RENDERING_COLOR_IPFZ, "/GROWTH3DGrainID", true, "dummy", this->t, this->X, this->step );
					else if ( outopt_localrendercolor == RENDERING_COLOR_RHO )
						ompcrit_write_voxeldata_h5( RENDERING_COLOR_RHO, "/GROWTH3DGrainID", true, "dummy", this->t, this->X, this->step );
					else
						cout << "ERROR::HDF5 output desired but utilizing unspecified color model!" << endl;
				}
			}

			//##DEBUG::ompcrit_write_voxeldata_coloring_grainids( "GROWTH.3DGrainID" );
			if ( this->outopt_localgenDAMASKgeom == OUTPUT_DAMASK_GEOMETRYFILE_YES ) 
				write_damask_geom_h5( "NetworkFinalState" );
		}
		#pragma omp barrier

		//all threads clean their data in parallel
		regions[omp_get_thread_num()]->omp_cleanupMemoryUsedForManagement();
		regions[omp_get_thread_num()]->omp_cleanupMemoryUsedForGrowthMachine();
		regions[omp_get_thread_num()]->omp_cleanupMemoryUsedForHaloRegion();

	} //implicit barrier and flush parallel execution done, go out when completely recrystallized, recrystallized to use defined value, when too many integration steps, when processing out
	omp_destroy_lock(&(this->h5lock));

	//output complete growth path of all grains in plain MPI I/O binary files
	if ( outopt_localsinglegrainevo != OUTPUT_SINGLEGRAIN_NO && outopt_localsinglegrainevo == OUTPUT_SINGLEGRAIN_BINARY )
		write_grainevolution_binary();

	cout << myensHdl->myRank << " myRank,jobid = " << this->jobid << " local CA simulation was successfully executed, tsimend = " << tsimend << ", stepsimend = " << stepsimend << endl;
}

void caHdl::solve_FINAL_IO( void )
{
	//add all I/O functionalities after the simulation executed
	if ( outopt_localrxfront == OUTPUT_RXFRONTSTATS_YES )
		write_rxfrontstats();

	if ( outopt_localthreadprof == OUTPUT_THREADPROFILING_YES )
		write_ThreadProfilingSummaryGrowthMachine();

	if ( outopt_localsinglegrainevo != OUTPUT_SINGLEGRAIN_NO ) {
		if ( outopt_localsinglegrainevo == OUTPUT_SINGLEGRAIN_ASCII ) 
			write_grainevolution_ascii();
		if ( outopt_localsinglegrainevo == OUTPUT_SINGLEGRAIN_BINARY )
			write_grainevolution_binary();
	}

	if ( percolation != PERCOLATION_ANALYZE_NO ) {
		write_PercolationProfilingSummary();
	}


	write_rxareaprofile_ascii();

	//##MK::currently not supported with MPI/OpenMP
	/*if ( renderingForThisCA == true ) {
		if ( outopt_localrenderhow == RENDERING_MS2D ) {
			if ( outopt_localrendercolor == RENDERING_COLOR_IPFZ )
				ompcrit_write_zsection_coloring_ipfz( DEFAULT_ZSECTIONING_ZPOS );
		}

		if ( outopt_localrenderhow == RENDERING_MS3D ) {
			if ( outopt_localrendercolor == RENDERING_COLOR_GRAINID ) {
				ompcrit_write_voxeldata_coloring_grainids( "FINAL.3DGrainID" );
				//ompcrit_write_voxeldata_coloring_regions();
			}

			if ( outopt_localrendercolor == RENDERING_COLOR_IPFZ )
				//##MK::currently deactivated with OpenMP
				//ompcrit_write_voxeldata_coloring_ipfz();
			}

			loginfo_rendering_cnt++;
		}
	}*/
}

#endif
