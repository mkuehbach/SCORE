//SCORE automaton developed by M K{\"u}hbach, 2014/2015, for questions and details contact markus.kuehbach@rwth-aachen.de


#include "SCORE_Kernel.h"
#include "lodepng.h"


inline bool SortDoubleAscending(const double &dbl1 , const double &dbl2) {
	return dbl1 < dbl2;
}

inline bool SortLoadFactorJobTicketAscending(const jobticket &jt1 , const jobticket &jt2) {
	return jt1.loadfactor < jt2.loadfactor;
}

inline bool SortGSDAscending( const agrain &ag1, const agrain &ag2 ) {
	return ag1.vol < ag2.vol;
}

inline bool SortFGSDCompleteAsc( const MPI_IO_FGSDComplete &fg1, const MPI_IO_FGSDComplete &fg2 ) {
	return fg1.finalvol < fg2.finalvol;
}

inline bool SortIntAscending( const uint32_t &it1, const uint32_t &it2 )
{
	return it1 < it2;
}



ensembleHdl::ensembleHdl()
{
	simid = DEFAULT_SIMID;

	ensRediscrTime = NULL;
	ensRediscrWindow.tensmin = 0.0;
	ensRediscrWindow.tensmax = 0.0;
	ensRediscrWindow.nslots = REDISCR_TIMESLOTS_DEFAULT;
	ensRediscrWindow.strategy = REDISCR_TIME_EQUIDISTANT;

	ensCAGeometry.nNucleiCSR = 1;
	ensCAGeometry.nboxedge_rd = CA_DIMENSIONS_MINIMUM;
	ensCAGeometry.nboxedge_nd = CA_DIMENSIONS_MINIMUM;
	ensCAGeometry.nboxedge_td = CA_DIMENSIONS_MINIMUM;
	ensCAGeometry.nboxarea_rdnd = SQR(CA_DIMENSIONS_MINIMUM);
	ensCAGeometry.nboxvol_rdndtd = CUBE(CA_DIMENSIONS_MINIMUM);
	ensCAGeometry.cellsize = DEFAULT_CELLSIZE;
	ensCAGeometry.boxedge_rd = CA_DIMENSIONS_MINIMUM * DEFAULT_CELLSIZE;
	ensCAGeometry.boxedge_nd = CA_DIMENSIONS_MINIMUM * DEFAULT_CELLSIZE;
	ensCAGeometry.boxedge_td = CA_DIMENSIONS_MINIMUM * DEFAULT_CELLSIZE;
	ensCAGeometry.boxarea_rdnd = SQR(CA_DIMENSIONS_MINIMUM)*SQR(DEFAULT_CELLSIZE);
	ensCAGeometry.boxvol_rdndtd = CUBE(CA_DIMENSIONS_MINIMUM)*CUBE(DEFAULT_CELLSIZE);

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
	ensPhysData.defgmean_nd = DEFAULT_DEFGSIZE;
	ensPhysData.defgmean_td = DEFAULT_DEFGSIZE;
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

	ensNucleationModel.gbnucleation = NO_GBNUCLEATION;
	ensNucleationModel.csrnucleation = NO_CSRNUCLEATION;
	ensNucleationModel.clustnucleation = NO_CLUSTNUCLEATION;
	ensNucleationModel.defaultnucdensity = DEFAULT_NUCDENSITY_CSR;

	maxfillperstep = DEFAULT_MAXFILLPERSTEP;
	ensMemGuard = 0.0;
	initialRelCellCaching = DEFAULT_RELCELLCACHING;
	transientRelCellRecaching = DEFAULT_TRANSRELCELL;
	XMAX = DEFAULT_XMAX;
	TMAX = DEFAULT_TMAX;
	NTSTEPSMAX = DEFAULT_NMAX;
	//Pmax = DEFAULT_PMAX;

	myRank = MASTER;
	nRanks = 1;
	nworldCAs = 0;
	nensCAs = 0;

	defmsmethod = CUBOID_DEFMS;
	outopt_rendermethod = RENDERING_MSNO;
	outopt_rendercolormodel = RENDERING_COLOR_GRAINID;
	outopt_localrenderboundaries = RENDERING_BOUNDARIES_NO;
	outopt_logboundaries = OUTPUT_LOGBND_NO;
	outopt_rxfront = OUTPUT_RXFRONTSTATS_NO;
	outopt_singlegrainevo = OUTPUT_SINGLEGRAIN_NO;

	ensembleSuccess = true;
	outopt_hpf3d2donmaster = false;

	ensembleprng.init( DEFAULT_PRNG_SEED );

	prof_t0 = 0.0;
	prof_tstartpp = 0.0;
	prof_tend = 0.0;
}


ensembleHdl::~ensembleHdl()
{
	//vector of caHdl call constructor of the caHdl class which has to assure to free memory in loginfo pointer snippets
	//randomClass has own destructor
	delete [] ensRediscrTime;
	ensRediscrTime = NULL;
	ensMemGuard = ensMemGuard - (ensRediscrWindow.nslots * sizeof(double));
}


uint32_t ensembleHdl::get_closest_standardlage( double  * quat )
{
	double closest_disori = MAX_FCC;
	uint32_t closest_ideal = CATEGORIZED_AS_RANDOM;
	double q10 = quat[0];
	double q11 = quat[1];
	double q12 = quat[2];
	double q13 = quat[3];


	for (uint32_t cand = 0; cand < standardlagen.size(); cand++) {
		double disori = misorientationCubicQxQ ( q10, q11, q12, q13,  standardlagen[cand].q0, standardlagen[cand].q1, standardlagen[cand].q2, standardlagen[cand].q3 );

		//scatter within scattering range and closer as to all other candidates?
		if ( disori <= standardlagen[cand].scatter && disori < closest_disori ) {
			closest_disori = disori;
			closest_ideal = cand;
		}
		//MK::dont break here for all possible standardlagen candidates to find the closest match even when scatter range overlaps
	}

	//MK::within the code orientations categorized as RANDOM are at the zero-th position

	if ( closest_ideal == CATEGORIZED_AS_RANDOM )
		return 0;

	//else
		return closest_ideal + 1;
}


uint32_t ensembleHdl::check_disjunctness( double * bunge )
{
	QUICKASSERT ( worldoripool.size() < MAXIMUM_DISJOINT_ORIS );

	uint32_t closestid = UNKNOWN_ORIENTATION;
	double disori = 2 * _PI_;
	double closestdisori = RESOLUTION_SO3GRID;
	double qbunge[4], qcand[4];

	euler2quaternion( bunge, qbunge );

	//check disorientation to all other components in worldoripool
	for (uint32_t cand = 0; cand < worldoripool.size(); cand++ ) {
		qcand[0] = worldoripool[cand].q0;
		qcand[1] = worldoripool[cand].q1;
		qcand[2] = worldoripool[cand].q2;
		qcand[3] = worldoripool[cand].q3;
		disori = misorientationCubicQxQ( qbunge[0], qbunge[1], qbunge[2], qbunge[3], qcand[0], qcand[1], qcand[2], qcand[3]);

		if (disori <= closestdisori) {
			closestdisori = disori;
			closestid = cand;
		} //MK::dont break here! find really the closest in the set!
	}

	if ( closestid == UNKNOWN_ORIENTATION ) {
		struct ori anori;
		anori.bunge1 = bunge[0];
		anori.bunge2 = bunge[1];
		anori.bunge3 = bunge[2];

		anori.q0 = qbunge[0];
		anori.q1 = qbunge[1];
		anori.q2 = qbunge[2];
		anori.q3 = qbunge[3];

		anori.closestideal = get_closest_standardlage( qbunge ); //RANDOM or one of our components

		anori.RGBA[RED] = UCHAR_RANGE_MIN;
		anori.RGBA[GREEN] = UCHAR_RANGE_MIN;
		anori.RGBA[BLUE] = UCHAR_RANGE_MIN;
		anori.RGBA[ALPHA] = UCHAR_RANGE_MAX;

		worldoripool.push_back( anori );

		//#ifdef DETAILED_PROMPTS
#ifdef REPORTSTYLE_DEVELOPER
		if (myRank == MASTER) {
			cout << "ADD\t\t" << (worldoripool.size() - 1) << "\t\t" << anori.bunge1 << "\t\t" << anori.bunge2 << "\t\t" << anori.bunge3;
			cout << "\t\t" << anori.q0 << "\t\t" << anori.q1 << "\t\t" << anori.q2 << "\t\t" << anori.q3 << endl;
			if ( anori.closestideal == 0 ) { cout << "Node;" << this->myRank << ", I have categorized close to RANDOM." << endl; }
			else { cout << "Node;" << this->myRank << ", I have categorized close to listentry " << (anori.closestideal + 1) << endl; }
		}
		//#endif DETAILED_PROMPTS
#endif

		return (worldoripool.size() - 1);
	}

	//else present closest already existent orientation id
#ifdef REPORTSTYLE_DEVELOPER
	if (myRank == MASTER) { cout << "GET\t\t" << closestid << "\t\t" << closestdisori << endl; }
#endif

	return closestid;
}


void ensembleHdl::init_ensprng( bool SetDissimilarSeedsForAll )
{
	//MK::all nodes running in parallel would see the same seed, thus generate exactly the same chain of random numbers!
	if ( SetDissimilarSeedsForAll == true) {
		long newseed = (long) (-1 * myRank);
		newseed = newseed - 1;
		ensembleprng.init ( newseed );

#ifdef REPORTSTYLE_DEVELOPER
		cout << myRank << " my new ensembleprng seed is:" << newseed << endl;
#endif

		//adjust mathPRNG, MK::works because there is ONLY ONE INSTANCE OF AN ENSEMBLEHDL PUBLIC SHARING mathMethods functionalities PER PROCESS!
		this->setprng( newseed );


		return;
	}
	else {
		ensembleprng.init ( DEFAULT_SEED );
	}
}


void ensembleHdl::init_mpidatatypes( void )
{
	//MPI_IO_CAProfilingData
	int elementCounts[2] = {2,5};
	MPI_Aint displacements[2] = {0, 2*8};
	MPI_Datatype oldTypes[2] = {MPI_LONG, MPI_DOUBLE};
	MPI_Type_create_struct(2, elementCounts, displacements, oldTypes, &MPI_IO_CAProfilingInfoData_Type);

	int elementCounts2[2] = {9, 1};
	MPI_Aint displacements2[2] = { 0, 9*8};
	MPI_Datatype oldTypes2[2] = {MPI_LONG, MPI_DOUBLE};
	MPI_Type_create_struct(2, elementCounts2, displacements2, oldTypes2, &MPI_IO_CAPhysicsInfoData_Type);


	int elementCounts3[2] = {2, 1};
	MPI_Aint displacements3[2] = { 0, 2*8};
	MPI_Datatype oldTypes3[2] = {MPI_DOUBLE, MPI_UNSIGNED };
	MPI_Type_create_struct(2, elementCounts3, displacements3, oldTypes3, &MPI_IO_FinalGSDInfoData_Type);

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


void ensembleHdl::init_parameter( void )
{
	dataBlockP parameterlist = readDataBlock("ParameterList", score_input);

		//only MASTER prompts
		nworldCAs = geTInt("CAEnsembleSize", parameterlist);
		if ( (nworldCAs < 1 || nworldCAs >= INTEGER_RANGE_MAX || nworldCAs > DEFAULT_MAX_NCAS) && myRank == MASTER ) { cout << "Invalid argument for the number of domains." << endl; }

		XMAX = geTReal("XMAX", parameterlist);
		if ( (XMAX < 0.0 || XMAX > 1.0) && myRank == MASTER ) { cout << "Invalid argument in XMAX." << endl; }

		TMAX = geTReal("TIMEMAX", parameterlist);
		if ( (TMAX < 0.0 || TMAX > DEFAULT_TMAX) && myRank == MASTER ) { cout << "Invalid argument in TMAX." << endl; }

		NTSTEPSMAX = geTInt("NMAX", parameterlist);
		if ( (NTSTEPSMAX < 1 || NTSTEPSMAX > DEFAULT_NMAX || NTSTEPSMAX >= INTEGER_RANGE_MAX) && myRank == MASTER ) { cout << "Invalid argument for NMAX." << endl; }

		//###currently implicit cubic domains
		ensCAGeometry.cellsize = geTReal("CellSize", parameterlist) * MICRON2METER;
		if ( (ensCAGeometry.cellsize < DEFAULT_CELLSIZE_MIN || ensCAGeometry.cellsize > DEFAULT_CELLSIZE_MAX) && myRank == MASTER ) { cout << "Invalid argument for the cellsize." << endl; }

		ensCAGeometry.nboxedge_rd = geTInt("3DCAEdgeLengthInCellsRD", parameterlist);
		if ( (ensCAGeometry.nboxedge_rd < CA_DIMENSIONS_MINIMUM || ensCAGeometry.nboxedge_rd > CA_DIMENSIONS_MAXIMUM) && myRank == MASTER ) { cout << "Invalid argument for automaton size RD." << endl; } 

		ensCAGeometry.nboxedge_nd = geTInt("3DCAEdgeLengthInCellsND", parameterlist);
		if ( (ensCAGeometry.nboxedge_nd < CA_DIMENSIONS_MINIMUM || ensCAGeometry.nboxedge_nd > CA_DIMENSIONS_MAXIMUM) && myRank == MASTER ) { cout << "Invalid argument for automaton size ND." << endl; } 

		ensCAGeometry.nboxedge_td = geTInt("3DCAEdgeLengthInCellsTD", parameterlist);
		if ( (ensCAGeometry.nboxedge_td < CA_DIMENSIONS_MINIMUM || ensCAGeometry.nboxedge_td > CA_DIMENSIONS_MAXIMUM) && myRank == MASTER ) { cout << "Invalid argument for automaton size TD." << endl; } 

		ensCAGeometry.nboxarea_rdnd = ensCAGeometry.nboxedge_rd * ensCAGeometry.nboxedge_nd;
		ensCAGeometry.nboxvol_rdndtd = ensCAGeometry.nboxedge_rd * ensCAGeometry.nboxedge_nd * ensCAGeometry.nboxedge_td;

		ensCAGeometry.boxedge_rd = ensCAGeometry.cellsize * ensCAGeometry.nboxedge_rd;
		ensCAGeometry.boxedge_nd = ensCAGeometry.cellsize * ensCAGeometry.nboxedge_nd;
		ensCAGeometry.boxedge_td = ensCAGeometry.cellsize * ensCAGeometry.nboxedge_td;
		ensCAGeometry.boxarea_rdnd = ensCAGeometry.nboxarea_rdnd * SQR(ensCAGeometry.cellsize);
		ensCAGeometry.boxvol_rdndtd = ensCAGeometry.nboxvol_rdndtd * CUBE(ensCAGeometry.cellsize); //##MK::accuracy issues

		ensPhysData.defgmean_rd = geTReal("MeanDefGrainSizeinRD", parameterlist) * MICRON2METER; // RD - x
		ensPhysData.defgmean_nd = geTReal("MeanDefGrainSizeinND", parameterlist) * MICRON2METER; // ND - y
		ensPhysData.defgmean_td = geTReal("MeanDefGrainSizeinTD", parameterlist) * MICRON2METER; // TD - z
		double mindiscretization = DEFAULT_MIN_GRAINDISCR * ensCAGeometry.cellsize;
		if ( (ensPhysData.defgmean_rd < mindiscretization || ensPhysData.defgmean_nd < mindiscretization || ensPhysData.defgmean_td < mindiscretization) && myRank == MASTER) { cout << "Invalid argument for the average size of the cuboid deformed grain." << endl; } //no upper bound check necessary

		defmsmethod = CUBOID_DEFMS;
		if ( geTInt("DefStructureSynthesis", parameterlist) == 1 ) defmsmethod = CUBOID_DEFMS;
		if ( geTInt("DefStructureSynthesis", parameterlist) == 2 ) defmsmethod = POISSONVORONOI_DEFMS;
		if ( geTInt("DefStructureSynthesis", parameterlist) == 3 ) defmsmethod = CPFEMDAMASK_DEFMS;

		ensPhysData.defgmean_poisson = geTReal("MeanDefSizePoissonDiameter", parameterlist) * MICRON2METER;
		if ( (ensPhysData.defgmean_poisson < mindiscretization && myRank == MASTER ) ) { cout << "Invalid argument for the average size of deformed grain for method poisson." << endl; }


		ensNucleationModel.defaultnucdensity = geTInt("NucleationDensityLocalCSR", parameterlist);
		if ( ensNucleationModel.defaultnucdensity >= INTEGER_RANGE_MAX && myRank == MASTER ) { cout << "Invalid NucleationDensityLocalCSR!" << endl; }
		long nnucpractical = (ensCAGeometry.nboxvol_rdndtd / ensNucleationModel.defaultnucdensity); //cells per nuc, assuming CSR, v=const
		if ( (ensNucleationModel.defaultnucdensity < DEFAULT_NUCDENSITY_CSR && myRank == MASTER) || (nnucpractical < MIN_DISCRETIZATION_GRAIN && myRank == MASTER) ) { cout << "Inacceptable low discretization for this high number of nuclei." << endl; }

		ensNucleationModel.csrnucleation = NO_CSRNUCLEATION;
		if ( geTInt("CSRNucleation", parameterlist) == 1 )
			ensNucleationModel.csrnucleation = CSR_ENFORCED;
		if ( geTInt("CSRNucleation", parameterlist) == 2 )
			ensNucleationModel.csrnucleation = CSR_PICKRANDOMLY;

		ensNucleationModel.clustnucleation = NO_CLUSTNUCLEATION;
		if ( geTInt("ClusteredNucleation", parameterlist) == 1 )
			ensNucleationModel.clustnucleation = CLUSTNUCLEATION_PICKRANDOMLY;

		ensNucleationModel.gbnucleation = NO_GBNUCLEATION;
		if ( geTInt("GBNucleation", parameterlist) == 1 )
			ensNucleationModel.gbnucleation = GB_ALLBOUNDARIES;
		if ( geTInt("GBNucleation", parameterlist) == 2 )
			ensNucleationModel.gbnucleation = GB_ONLYATHAGB;

		if ( ensNucleationModel.csrnucleation == NO_CSRNUCLEATION && ensNucleationModel.gbnucleation == NO_GBNUCLEATION && ensNucleationModel.clustnucleation == NO_CLUSTNUCLEATION ) {
			cout << "ERROR::Choose a valid nucleation model!" << endl;
			ensNucleationModel.csrnucleation = CSR_ENFORCED;
		}

		//clustering modeling
		ensNucleationModel.cluster_nclust = geTInt( "ClusteredNucleationNCluster", parameterlist);
		ensNucleationModel.cluster_lambda = geTReal( "ClusteredNucleationLCluster", parameterlist);
		ensNucleationModel.cluster_rvesize = geTReal( "ClusteredNucleationScaling", parameterlist) * MICRON2METER;
		ensNucleationModel.cluster_a = geTInt( "ClusteredNucleationExtendA", parameterlist);
		ensNucleationModel.cluster_b = geTInt( "ClusteredNucleationExtendB", parameterlist);
		ensNucleationModel.cluster_c = geTInt( "ClusteredNucleationExtendC", parameterlist);

		//nucleation modeling
		ensNucleationModel.gbnucleation_dens2num = geTReal("GBNucDensity2Number", parameterlist);
		ensNucleationModel.gbnucleation_drho2dens = geTReal("GBNucRhoDiff2Density", parameterlist);
		ensNucleationModel.gbnucleation_scatter = geTReal("GBNucMaximumScatter", parameterlist) / 180.0 * _PI_;
		ensNucleationModel.tincub_rayleigh_sigma = geTReal("IncubationTimeScatter", parameterlist);
		ensNucleationModel.tincubmodel = TINCUB_SITESATURATION;
		if ( geTInt("IncubationTimeModel", parameterlist) == 2 ) {
			ensNucleationModel.tincubmodel = TINCUB_TIMEDEPENDENT;
			if ( ensNucleationModel.tincub_rayleigh_sigma < MINIMUM_RAYLEIGH_SIGMA ) {
				cout << "ERROR::Invalid low incubation time scatter!" << endl;
			}
		}

		mobilitymodel = MOBILITYMODEL_SEBALDGOTTSTEIN;
		if ( geTInt("MobilityModel", parameterlist) == 2 )
			mobilitymodel = MOBILITYMODEL_ROLLETTHOLM;

		//additional output
		outopt_rendermethod = RENDERING_MSNO;
		if ( geTInt("RenderingOfMicrostructure", parameterlist) == 2 )
			outopt_rendermethod = RENDERING_MS2D;
		if ( geTInt("RenderingOfMicrostructure", parameterlist) == 3 )
			outopt_rendermethod = RENDERING_MS3D;

		outopt_rendercolormodel = RENDERING_COLOR_GRAINID;
		if ( geTInt( "RenderingColorModel", parameterlist) == 2 ) 
			outopt_rendercolormodel = RENDERING_COLOR_IPFZ;

		outopt_localrenderboundaries = RENDERING_BOUNDARIES_NO;
		if ( geTInt("RenderingBoundaries", parameterlist) == 1 ) 
			outopt_localrenderboundaries = RENDERING_BOUNDARIES_YES;

		outopt_logboundaries = OUTPUT_LOGBND_NO;
		if ( geTInt( "OutputLogBoundaries", parameterlist) == 1 )
			outopt_logboundaries = OUTPUT_LOGBND_YES;

		outopt_rxfront = OUTPUT_RXFRONTSTATS_NO;
		if ( geTInt( "OutputRXFrontStats", parameterlist) == 1 )
			outopt_rxfront = OUTPUT_RXFRONTSTATS_YES;

		outopt_singlegrainevo = OUTPUT_SINGLEGRAIN_NO;
		if ( geTInt( "OutputSingleGrainStats", parameterlist) == 1 ) 
			outopt_singlegrainevo = OUTPUT_SINGLEGRAIN_YES;

		outopt_hpf3d2donmaster = false;
		if ( geTInt("OutputOfHPF", parameterlist) == 1 ) { outopt_hpf3d2donmaster = true; }


	delete parameterlist;


	dataBlockP  integratoraccuracy = readDataBlock("IntegratorAccuracy", score_input);

		maxfillperstep = geTReal("MAXFILLPERSTEP", integratoraccuracy);
		if ( (maxfillperstep < DEFAULT_MINFILLIN || maxfillperstep > DEFAULT_MAXFILLIN) && myRank == MASTER ) { maxfillperstep = DEFAULT_MAXFILLPERSTEP; cout << "MAXFILLPERSTEP was chosen wrong, should be 0.01 < 0.2, was reset to = " << maxfillperstep << endl; }
		//define further values relevant for the numerical integration scheme

		int rediscr_nslots = geTInt("RediscretizationSteps", integratoraccuracy);
		if ( rediscr_nslots < REDISCR_TIMESLOTS_MIN || rediscr_nslots > REDISCR_TIMESLOTS_MAX || rediscr_nslots >= INTEGER_RANGE_MAX ) { cout << "WARNING::Rediscretization time slots reset to " << REDISCR_TIMESLOTS_DEFAULT << endl; }
		ensRediscrWindow.nslots = rediscr_nslots;

	delete integratoraccuracy;


	dataBlockP performancepara = readDataBlock("PerformanceParameter", score_input);

		initialRelCellCaching = geTReal("InitialCacheSizeRXCells", performancepara);
		if ( (initialRelCellCaching < DEFAULT_MINRELCELL || initialRelCellCaching > DEFAULT_MAXRELCELL) && myRank == MASTER ) { initialRelCellCaching = DEFAULT_RELCELLCACHING; cout << "InitialCacheSizeRXCells was chosen unwise, it is reset to = " << initialRelCellCaching << endl; }
		transientRelCellRecaching = geTReal("ReCacheSizeRXCells", performancepara);
		if ( (transientRelCellRecaching < DEFAULT_MINTRANSRELCELL || transientRelCellRecaching > DEFAULT_MAXTRANSRELCELL) && myRank == MASTER ) { transientRelCellRecaching = DEFAULT_TRANSRELCELL; cout << "ReCacheSizeRXCells was chosen unwise, it is reset to = " << transientRelCellRecaching << endl; }

	delete performancepara;


	dataBlockP matpara = readDataBlock("MaterialProperties", score_input);
		//no boundary checks for these values!
		ensPhysData.G0 = geTReal("ZeroKelvinShearModulusG0", matpara);
		ensPhysData.dGdt = geTReal("FirstOrderdG0dT", matpara);
		ensPhysData.bZeroCelsius = geTReal("ZeroCelsiusBurgersVector", matpara);
		ensPhysData.thermexp_C = geTReal("AlloyConstantThermalExpCoeff", matpara); //thermal lattice expansion model
		ensPhysData.thermexp_a = geTReal("FirstOrderThermalExpCoeff", matpara);
		ensPhysData.thermexp_b = geTReal("SecondOrderThermalExpCoeff", matpara);
		ensPhysData.Tmelt = geTReal("MeltingTemperature", matpara) + TOFFSET;

		//parameter of the Sebald-Gottstein inspired mobility model in which the mobility is solely dependent on the disorientation of the two crystals
		//adjointing the boundary, when the disorientation < LAGB_TO_HAGB_TRANSITION the boundary is essentially immobile, otherwise mobile as a HAGB
		//it is considered particularly mobile when the disorientation is close to a 40deg111 in misorientation space
		ensPhysData.LAGBm0 = geTReal("LAGBm0", matpara);
		ensPhysData.LAGBHact = geTReal("LAGBHact", matpara) * echarge; //all other codes work with 1.602e-19
		ensPhysData.HAGBm0 = geTReal("HAGBm0", matpara);
		ensPhysData.HAGBHact = geTReal("HAGBHact", matpara) * echarge;
		ensPhysData.GSm0 = geTReal("GSm0", matpara);
		ensPhysData.GSHact = geTReal("GSHact", matpara) * echarge;


		//on the other hand Rollett and Holm for instance often consider the LAGB_TO_HAGB transition not as sharp and despite the possibility of 
		//close to 40deg111 boundaries modelled additionally more mobile than random high-angle grain boundaries (HAGB) they consider a 
		//sigmoidal function of the form m = 1.0 - cut*exp(-trans * (disori/LAGB_TO_HAGB_TRANSITION)^exponent) also inspired by exp. results on tilt boundaries in Al performed by Winning and Molodov
		ensPhysData.RH_HAGBm0 = geTReal("RHModelHAGBm0", matpara );
		ensPhysData.RH_HAGBHact = geTReal("RHModelHAGBHact", matpara ) * echarge;
		ensPhysData.RH_LAGBHAGBcut = geTReal("RHModelLAGBHAGBCut", matpara );
		ensPhysData.RH_LAGBHAGBtrans = geTReal("RHModelLAGBHAGBTrans", matpara );
		ensPhysData.RH_LAGBHAGBexponent = geTReal("RHModelLAGBHAGBExponent", matpara );


		if ( ensPhysData.Tmelt < 0.0  &&  myRank == MASTER ) { cout << "Tmelting is invalid." << endl; }

	delete matpara;

	//caching variables Gbhalfsq are becoming updated in the main loops

	//RECOVERY
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
	delete recoverypara;


	ensDragData.ZenerConsider = DISPERSOIDDRAG_NO;

	dataBlockP zenerpara = readDataBlock("ZenerDragParameter", score_input );
		int zenermodel = geTInt( "ZenerConsider", zenerpara );
		if ( zenermodel == 1 ) ensDragData.ZenerConsider = DISPERSOIDDRAG_CONSTANT;
		if ( zenermodel == 2 ) ensDragData.ZenerConsider = DISPERSOIDDRAG_TIMEDEP;
		ensDragData.ZenerAlpha = geTReal( "ZenerAlpha", zenerpara );
		ensDragData.ZenerGamma = geTReal( "ZenerGamma", zenerpara );

	delete zenerpara;

	
	init_zenerdrag();


	init_idealorientations();


	if ( outopt_rendermethod == RENDERING_MS2D || outopt_rendermethod == RENDERING_MS3D || outopt_localrenderboundaries == RENDERING_BOUNDARIES_YES ) {
		dataBlockP renderlogpoints = readDataBlock("RenderMicrostructure", score_input );

			dataLineP rline = renderlogpoints->first;
			long nrlines = renderlogpoints->lineCount;

			QUICKASSERT( nrlines < DEFAULT_MAX_INPUTDATALEN );

			double Xri = 0.0;
			for ( long r = 0; r < nrlines; r++ ) { 
				Xri = getReal( rline, 2 );

				UserDefLogPoint_MS_Rendering.push_back( Xri );

				rline = rline->next;
			}
			delete rline;

		delete renderlogpoints;

		std::sort ( UserDefLogPoint_MS_Rendering.begin(), UserDefLogPoint_MS_Rendering.end(), SortDoubleAscending );

		//import which simulations to render
		dataBlockP candidatelist = readDataBlock("RenderOnlyTheseRegions", score_input );

			dataLineP cline = candidatelist->first;
			long nclines = candidatelist->lineCount;

			int id = 0;

//cout << myRank << "Render only these regions = " << nclines << ";" << cline << ";" << id << endl;
			for ( long r = 0; r < nclines; r++ ) { 
				id = getInt ( cline, 1 );

//cout << "\t\t" << r << ";" <<  myRank << "Render only these regions = " << nclines << ";" << cline << ";" << id << endl;

				if ( id < 0 || id >= nworldCAs ) { cout << "Invalid simulation index listed that should become rendered!" << endl; return; }

				UserDefLogPoint_WhichCAtoOutput.push_back( (uint32_t) id );

				cline = cline->next;
			}
			delete cline;

		delete candidatelist;

		std::sort( UserDefLogPoint_WhichCAtoOutput.begin(), UserDefLogPoint_WhichCAtoOutput.end(), SortIntAscending );
	}

	init_processing();

	//import defragmentation points
	dataBlockP defragpoints = readDataBlock( "HeuristicRXFrontListDefragmentation", score_input);
		dataLineP dline = defragpoints->first;
		long ndlines = defragpoints->lineCount;

		QUICKASSERT( ndlines < DEFAULT_MAX_INPUTDATALEN );

		double Xdi = 0.0;
		for ( long d = 0; d < ndlines; d++ ) { 
			Xdi = getReal( dline, 2 );

			UserDefLogPoint_X_CellListDefragment.push_back( Xdi );

			dline = dline->next;
		}
		delete dline;
	delete defragpoints;


	//better sort if user has forgotten to do so
	std::sort ( UserDefLogPoint_X_CellListDefragment.begin(), UserDefLogPoint_X_CellListDefragment.end(), SortDoubleAscending );


	//import logpoints that are interpreted locally in each caHdl dependent on the local recrystallized fraction
	dataBlockP logpoints = readDataBlock( "SingleGrainVolOverTime", score_input);
		dataLineP pline = logpoints->first;
		long nplines = logpoints->lineCount;

		QUICKASSERT( nplines < DEFAULT_MAX_INPUTDATALEN );

		double Xi = 0.0;
		for ( long pp = 0; pp < nplines; pp++ ) { 
			Xi = getReal( pline, 2 );

			UserDefLogPoint_X_Output.push_back( Xi ); 	//no reserve because single page already allocated

			pline = pline->next;
		}
		delete pline;
	delete logpoints;

	//better sort the list in this vector if the user has forgotten to do so, what if two elements are numerically indistinguishable, then yes well what about the universe would be a coffee cup...
	std::sort ( UserDefLogPoint_X_Output.begin(), UserDefLogPoint_X_Output.end(), SortDoubleAscending );


	//output
	if (myRank == MASTER) {
		cout << "I have been made familiar with " << standardlagen.size() << " ideal texture components that are used for texture volume analyses."  << endl;

		cout << UserDefLogPoint_X_CellListDefragment.size() << ", logpoints for RXFront defragmentation were initialized as desired and sorted ascending." << endl;
		cout << UserDefLogPoint_X_Output.size() << ", logpoints for X were initialized as desired and sorted in ascending order." << endl;
		cout << XMAX << ", is the upper limit for RXFRACTION numerical integration scheme." << endl;
		cout << TMAX << ", is the upper limit for REALTIME numerical integration scheme." << endl;
		cout << initialRelCellCaching << ", is the initial recaching value." << endl;
		cout << nworldCAs << ", is the total number of automata in the MPI_COMM_WORLD." << endl;
		cout << ensCAGeometry.cellsize << ", is the cellsize." << endl;
		cout << maxfillperstep << ", is the MAXFILLPERSTEP." << endl;
		cout << ensCAGeometry.nboxedge_rd << ";" << ensCAGeometry.nboxedge_nd << ";" << ensCAGeometry.nboxedge_td << ", is RD x ND x TD CA geometry." << endl;
		cout << ensPhysData.defgmean_rd << "," << ensPhysData.defgmean_nd << "," << ensPhysData.defgmean_td << ", defgmedian_rdndtd." << endl;
		cout << ensPhysData.G0 << "," << ensPhysData.bZeroCelsius << ", ShearModulus/Burgersvector." << endl;
		cout << ensPhysData.LAGBm0 << "," << ensPhysData.HAGBm0 << "," << ensPhysData.GSm0 << ",  LAGB/HAGB/GS m0." << endl;
		cout << ensPhysData.LAGBHact << "," << ensPhysData.HAGBHact << "," << ensPhysData.GSHact << ", LAGB/HAGB/GS Hact=" << endl;

		cout << "Initial parameter have been loaded successfully." << endl;
	}
}


void ensembleHdl::init_idealorientations( void )
{
	//import ideal texture components
	dataBlockP idealcomp = readDataBlock( "IdealComponents", score_input);
		dataLineP iline = idealcomp->first;
		long nilines = idealcomp->lineCount;

		QUICKASSERT( nilines <= DEFAULT_MAX_INPUTDATALEN );

		double eul[3], qbunge[4], scatt;

		for ( long ii = 0; ii < nilines; ii++ ) {
			struct ideal icomp;

			eul[0] = getReal( iline, 2 ) / 180.0 * _PI_;
			eul[1] = getReal( iline, 3 ) / 180.0 * _PI_;
			eul[2] = getReal( iline, 4 ) / 180.0 * _PI_;
			scatt  = getReal( iline, 5 ) / 180.0 * _PI_;

			euler2quaternion ( eul, qbunge );

			icomp.bunge1 = eul[0];
			icomp.bunge2 = eul[1];
			icomp.bunge3 = eul[2];
			icomp.q0 = qbunge[0];
			icomp.q1 = qbunge[1];
			icomp.q2 = qbunge[2];
			icomp.q3 = qbunge[3];
			icomp.scatter = scatt;

			standardlagen.push_back( icomp );

			#ifdef REPORTSTYLE_DEVELOPER
				cout << "Idealcomponent imported::bunge123;q0123;scatt;" << standardlagen[ii].bunge1 << ";" << standardlagen[ii].bunge2 << ";" << standardlagen[ii].bunge3 << ";" << standardlagen[ii].q0 << ";" << standardlagen[ii].q1 << ";" << standardlagen[ii].q2 << ";"<< standardlagen[ii].q3 << ";" << standardlagen[ii].scatter << endl;
			#endif

			iline = iline->next;
		}
		delete iline;
	delete idealcomp;
}


void ensembleHdl::init_processing( void )
{
	dataBlockP timetemp = readDataBlock( "ProcessingSchedule", score_input );
		dataLineP tline = timetemp->first;
		long ntlines = timetemp->lineCount;
		if ( myRank == MASTER && ntlines < 2 ) { cout << "ERROR::An insufficient number (<2) of Time-Temperature data have been provided to perform the simulation!" << endl; }

		QUICKASSERT( ntlines <= DEFAULT_MAX_INPUTDATALEN );

		double ti = 0.0;
		double Ti = (25.0 + TOFFSET);

		for ( long ts = 0; ts < ntlines; ts++ ) {
			ti = getReal(tline, 1);
			Ti = getReal(tline, 2);
			Ti += TOFFSET;

			ttime.push_back( ti );
			ttemperature.push_back( Ti );

			#ifdef REPORTSTYLE_DEVELOPER
				cout << "Node;" << this->myRank << ";ti;Ti;" << ti << ";" << Ti << "in the vector instead ttime[ts];ttemperature[ts];" << "\ts=" << ts << "\t" << ttime[ts] << ";" << ttemperature[ts] << endl;
			#endif

			tline = tline->next;
		}

		delete tline;
	delete timetemp;
}


void ensembleHdl::init_zenerdrag( void )
{
	dataBlockP particleevo = readDataBlock( "EvoDraggingParticles" , score_input );
	dataLineP zline = particleevo->first;
		long nzlines = particleevo->lineCount;
		if ( myRank == MASTER && nzlines < 2 ) { cout << "ERROR::An insufficient number (<2) of Particle Time Dispersion data have been provided to perform the simulation!" << endl; }

		QUICKASSERT( nzlines <= DEFAULT_MAX_INPUTDATALEN );

		double ti = 0.0;
		double fri = 0.0;

		for ( long ts = 0; ts < nzlines; ts++ ) {
			ti = getReal(zline, 1);
			fri = getReal(zline, 2);

			dispersoidtime.push_back( ti );
			dispersoidfr.push_back( fri );

			#ifdef REPORTSTYLE_DEVELOPER
				cout << "Node;" << this->myRank << ";ti;fri;" << ti << ";" << fri << "in the vector instead Time[ts];Dispersion[ts];" << "\ts=" << ts << "\t" << dispersoidtime[ts] << ";" << dispersoidfr[ts] << endl;
			#endif

			zline = zline->next;
		}

		delete zline;
	delete particleevo;
}


void ensembleHdl::init_defgpool( void )
{
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
		double bunge[3] = {0.0, 0.0, 0.0};
		double q[4] = {1.0, 0.0, 0.0, 0.0};

		for (long dg = 0; dg < ndefgpool; dg++) {
			struct defg agrain;

			bunge[0] = getReal(dline, 1) / 180.0 * _PI_;
			bunge[1] = getReal(dline, 2) / 180.0 * _PI_;
			bunge[2] = getReal(dline, 3) / 180.0 * _PI_;

			id = check_disjunctness( bunge );

			agrain.ori = id;
			agrain.rho0 = getReal(dline, 4);
			agrain.rho = getReal(dline, 4);
			agrain.dav0 = getReal(dline, 5); //##MK::currently not ingrain ori spread or subgrain size distribution
			agrain.dav = getReal(dline, 5);
			agrain.avdg0 = getReal(dline, 6) / 180.0 * _PI_;
			agrain.avdg = getReal(dline, 6) / 180.0 * _PI_;


			worlddefgpool.push_back(agrain);

			//cout << "dg;ori;e1;e2;e3;q0;q1;q2;q3;rho0;rho;avdg;" << dg << ";" << worlddefgpool[dg].ori << ";" << worldoripool[id].bunge1 << ";" << worldoripool[id].bunge2 << ";" << worldoripool[id].bunge3 << ";" << worldoripool[id].q0 << ";" << worldoripool[id].q1 << ";" << worldoripool[id].q2 << ";" << worldoripool[id].q3 << ";" << worlddefgpool[dg].rho0 << ";" << worlddefgpool[dg].rho << ";" << worlddefgpool[dg].avdg << endl;

			dline = dline->next;
		} // for all grains in the pool

	delete defgpara;

	if (myRank == MASTER) { cout << "All deformed grains have been successfully imported, highest dislocation density detected  is." << endl; }
}


void ensembleHdl::init_rxgpool( void )
{
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
		double bunge[3] = {0.0, 0.0, 0.0};
		double q[4] = {1.0, 0.0, 0.0, 0.0};

		for ( long rg = 0; rg < nrxgpool; rg++) {
			struct rxg agrain;

			bunge[0] = getReal(rline, 1) / 180.0 * _PI_;
			bunge[1] = getReal(rline, 2) / 180.0 * _PI_;
			bunge[2] = getReal(rline, 3) / 180.0 * _PI_;

			id = check_disjunctness( bunge );

			agrain.ori = id;
			agrain.tincub = getReal(rline, 4);

			worldrxgpool.push_back(agrain);

			//cout << "rg;ori;e1;e2;e3;q0;q1;q2;q3;tincub;" << rg << ";" << worldrxgpool[rg].ori << ";" << worldoripool[id].bunge1 << ";" << worldoripool[id].bunge2 << ";" << worldoripool[id].bunge3 << ";" << worldoripool[id].q0 << ";" << worldoripool[id].q1 << ";" << worldoripool[id].q2 << ";" << worldoripool[id].q3 << ";" << worldrxgpool[rg].tincub << endl;

			rline = rline->next;
		} // for all grains in the pool

	delete rxgpara;

	if (myRank == MASTER) { cout << nrxgpool << " recrystallizing grain orientations have been successfully imported." << endl; }
}


inline double ensembleHdl::calculateWorkload( uint32_t ncells, double averagefillperstep )
{
	//strategy is fastest boundary swept maxfillperstep - so two CAs with same ncells but different max1 > max2 -> max1 is faster and has less loop iteration overhead
	//workload of an automaton foremost dependent on the amount the amount of cells O(n^3) and the path heterogeneity/velocity heterogeneity

	//###there is an approximate measure for the time-dependent averagefillperstep
	double workload = (double) ncells * (1.0 / averagefillperstep);

	//if (workload >= worldmaxworkload) { worldmaxworkload = workload; }

	return workload;
}


inline double ensembleHdl::estimateMemory( uint32_t ncells, double frontcelloverhead, uint32_t avnuclei, uint32_t logsteps )
{
	//strategy for memory is:
	double estMemory = 0.0;
	//static memory: deformed cells
	estMemory += (ncells * DEFAULT_MEMPERCELL_DEF);
	//dynamic memory: front cells
	estMemory += (ncells * frontcelloverhead * DEFAULT_MEMPERCELL_RX);
	//dynamic memory: for intermediate stages log
	estMemory += (avnuclei * logsteps * DEFAULT_MEMPERGRAIN_LOG);
	//arbitrary additional overhead, e.g. texture and loginfo, local variable cache copies
	estMemory += (1.0 * DEFAULT_MEMPERCA_OVERHEAD);

	return estMemory;
}


void ensembleHdl::init_distributeWorkOnRanks( void )
{
	//##MK::MIND THAT IT IS NECESSARY TO HAVE A MPI_Barrier ( MPI_COMM_WORLD) if load should be distributed iteratively

	//calculate expected load for each automaton
	/*double ActiveCellOverhead = initialRelCellCaching;

	vector<jobticket> expectedloadCAs;	//##maybe it is not necessary to have all nodes know about all jobs
	expectedloadCAs.reserve( nworldCAs );

	//no Barrier necessary when distributeWorkOnRanksOnRanks separates work without iteratively packing thus all process come to the same partitioning plan
	for (int ca = 0; ca < nworldCAs; ca++) {
		struct jobticket jt;

		jt.jobid = ca;
		jt.estimatedMemory = estimateMemory (ensCAGeometry.nboxvol_rdndtd, ActiveCellOverhead, defaultnucleationdensity, UserDefLogPoint_X_Output.size() );
		jt.loadfactor = calculateWorkload (ensCAGeometry.nboxvol_rdndtd, maxfillperstep);

		expectedloadCAs.push_back ( jt );
	}

	//sort after load, all come to this result
	std::sort ( expectedloadCAs.begin(), expectedloadCAs.end(), SortLoadFactorJobTicketAscending );*/


	//partition work according to load
	if ( myRank == MASTER ) {
		if ( (nworldCAs % nRanks) != 0 ) { cout << "WARNING::Choose number of nuclei as a multiple number of nRanks to prevent excessive microstopping before barriers." << endl; }
	}

	//work partitioning at the moment simple a modulo round robin operation ensuring that ALL RANKS come to the same deterministic partitioning scheme
	//MK::modulo on uint % int is implementation specific and possibly non positive when int is negative or zero then even undefined, however nRanks > 0
	uint32_t luckyRank;

	for (uint32_t ca = 0; ca < nworldCAs; ca++) {

		luckyRank = ca % nRanks; //number % x with x = 0 is undefined, but nRanks >= 1

		WhichCAonWhichRank.push_back(luckyRank);

		if ( myRank == luckyRank ) {
			myIDs.push_back ( ca );

			//cout << myRank << "\t\t" << ca << "\t\tmyRank;ca" << endl;
		}
	} //for all nuclei


	#ifdef REPORTSTYLE_DEVELOPER
		cout << "myRank " << myRank << " has been assigned " << myIDs.size() << " work pieces." << endl;
	#endif

	if (myRank == MASTER) { cout << "The work has been partitioned onto the processes." << endl; }
}


//##just an example of otherwise fixed external constant ENVIRONMENT variable OMP_NUM_THREADS
#define OMP_NUM_THREADS					1
#define DEBUG_UTILIZE_ONLY_ONE_THREAD	0


void ensembleHdl::init_distributeWorkOnStartedThreads( void )
{
	uint32_t luckyThread = 0; //omp is global
	int omp_get_num_threads = OMP_NUM_THREADS; //MK >0

	//sequential partitioning of the work, then the thread team is initialized
	for ( uint32_t myca = 0; myca < myIDs.size(); myca++ ) {
		//##here starting with simple partitioning on local threads by a modulo operation
		//##though, better approach would be to utilize estimateWorkload to devide further in a guided manner

		luckyThread = myca % omp_get_num_threads;

		WhichMyCAonWhichThread.push_back(luckyThread);
//cout << myRank << "\t\t" << luckyThread << "\t\t" << myIDs[myca] << "\t\tmyRank/luckyThread/myIDs[myca]" << endl;
	}

	//#pragma omp parallel ... spawn a group of threads here ...
	//only first thread does this in MPI way
	//MK::IT WILL BE CRUCIAL FOR CCNUMA TYPE ARCHITECTURES HERE TO HAVE THE POINTER TO THE CAHDLs and these class objects in the local memory of socket and core-bound threads!

	for ( uint32_t tid = 0; tid < omp_get_num_threads; tid++) {
		vector<caHdlP> emptyBucketThreadLocalCAs;
		myCAs.push_back( emptyBucketThreadLocalCAs);
	}
	//distribute the pages or rely on interconnect bus (ICB) delivering... vector-vector construct but only for <= myIDs.size() caHdl pointer, so should be acceptable...

	//##MK::pragma omp parallel for 
	int omp_get_thread_num = DEBUG_UTILIZE_ONLY_ONE_THREAD;
	for ( uint32_t myca = 0; myca < myIDs.size(); myca++) {
		if ( WhichMyCAonWhichThread[myca] == omp_get_thread_num ) {
			//shall I work on this automaton - initialize and firsttouch in thread local memory!
			caHdlP aca = NULL;
			aca = new caHdl;
			QUICKASSERT( aca != NULL );

			ensembleHdl::myCAs[omp_get_thread_num].push_back( aca );
			ensMemGuard = ensMemGuard + sizeof(caHdl);

			//MK::assignment after push_back because aca is not a trivial object
			caHdlP tmpaccessthatca = myCAs[omp_get_thread_num][myCAs[omp_get_thread_num].size()-1];
			tmpaccessthatca->myensHdl = this;
			tmpaccessthatca->ccnumatid = WhichMyCAonWhichThread[myca];
			tmpaccessthatca->jobid = this->myIDs[myca]; //ensembleHdl::myIDs[myca];
			tmpaccessthatca->myensRank = this->myRank; //ensembleHdl::myRank;
			tmpaccessthatca->nRanks = this->nRanks; //ensembleHdl::nRanks;

//cout << "myRank;cajobid;myIDs[myca];camyensRank;canRanks;" << myRank << ";" << tmpaccessthatca->jobid << ";" << myIDs[myca] << ";" << tmpaccessthatca->myensHdl << ";" << tmpaccessthatca->nRanks << endl;

			//rendering in aca?
			tmpaccessthatca->renderingForThisCA = false;
			uint32_t ncandmax = ensembleHdl::UserDefLogPoint_WhichCAtoOutput.size();
//cout << "RenderingBoundaries ncandmax=" << ncandmax << endl;
			for ( uint32_t c = 0; c < ncandmax; c++ ) {
//cout << "RenderingBoundariesYES/NO\t\t" << myRank << "\t\t" << c << "\t\t" << ensembleHdl::UserDefLogPoint_WhichCAtoOutput[c] << "\t\t" << tmpaccessthatca->jobid << "Check..." << endl;
				if ( ensembleHdl::UserDefLogPoint_WhichCAtoOutput[c] == tmpaccessthatca->jobid ) {
//cout << "RenderingBoundariesYES/NO\t\t" << myRank << "\t\t" << c << "\t\t" << ensembleHdl::UserDefLogPoint_WhichCAtoOutput[c] << "\t\t" << tmpaccessthatca->jobid << "YES!" << endl;
					tmpaccessthatca->renderingForThisCA = true;
					break;
				}
			}

#ifdef REPORTSTYLE_DEVELOPER
			cout << tmpaccessthatca->myensRank << "\t\t" << myRank << "\t\t" << myca << "\t\t" << tmpaccessthatca->myensHdl << "\t\t" << tmpaccessthatca->ccnumatid << "\t\t" << tmpaccessthatca->jobid << "\t\t" << "\ttmpaccessthatca->myensRank;myRank;myca;myensHdlAddress;ccnumatid;jobid" << endl;
#endif

		}
	} //for all CAs in the process nail onto bounded! threads
}


void ensembleHdl::SIMULATE_myCAs( void )
{
	//now as distributeWorkOnStartedThreads has created ccNUMA-aware thread local ensembleHdl caHdl objects
	//##OPENMP Task parallelism on a vector of pointer to caHdl class objects each of which hosts one cellular automaton
	//##ON CCNUMA ARCHITECTURES IT IS CRUCIAL TO HAVE THE INDIVIDUAL THREADS INITIALIZE AND FIRSTTOUCH THE AUTOMATONS ACCORDING TO THEIR THREADID
	//##SO THAT ONLY THREAD LOCAL MEMORY IS DOMINANTLY USED AND THE LOAD ON THE INTERCONNECTION BUS STAYS MINIMAL
	//##current strategy -> each thread one automaton at a time other myCAs of the thread sequential
	//fine-scaled parallelism is thread localized
	uint32_t omp_get_num_threads = OMP_NUM_THREADS;

	//outer loop that tasks out sequential at the moment only the attributed automata
	uint32_t omp_get_thread_num = DEBUG_UTILIZE_ONLY_ONE_THREAD;

	//#pragma omp parallel out
	for ( uint32_t tid = 0; tid < omp_get_num_threads; tid++) {
		for ( uint32_t mycathreaded = 0; mycathreaded < myCAs[tid].size(); mycathreaded++ ) {
			//analyze the automaton theca with private properties theensHdl  //ensembleHdlP theens = theca->myensHdl; have to be private temporaries to the thread
			struct profilingData prof;
			double timer[5] = { 0.0, 0.0, 0.0, 0.0, 0.0 };
			timer[0] = MPI_Wtime();


			caHdlP theca = myCAs[tid][mycathreaded];

			theca->solve_INITIALIZATION();


			theca->solve_SYNTHETIZE_DEFORMATIONSTRUCTURE();
			timer[1] = MPI_Wtime();
			prof.MPIWTimeSpendDefMS = timer[1] - timer[0];

			if ( theca->myNucleationModel.gbnucleation != NO_GBNUCLEATION ) {
				theca->solve_DETECT_GRAINBOUNDARIES(); 
			}
			timer[2] = MPI_Wtime();
			prof.MPIWTimeSpendGBDetection = timer[2] - timer[1];


			//nucleation modeling in the deformation structure
			if ( theca->myNucleationModel.gbnucleation != NO_GBNUCLEATION ) {
				theca->solve_GRAINBOUNDARYNUCLEATION();
			}

			if ( theca->myNucleationModel.csrnucleation != NO_CSRNUCLEATION || theca->myNucleationModel.clustnucleation != NO_CLUSTNUCLEATION ) {
				if ( theca->myNucleationModel.csrnucleation == CSR_ENFORCED || theca->myNucleationModel.csrnucleation == CSR_PICKRANDOMLY || theca->myNucleationModel.clustnucleation == CLUSTNUCLEATION_PICKRANDOMLY ) {
					theca->solve_NUCLEATIONMODELING();
				}
			}
			timer[3] = MPI_Wtime();
			prof.MPIWTimeSpendNucleation = timer[3] - timer[2];

			theca->solve_REPLACE_CA_STRUCTURE();

			//MK::now everything is prepared for placing nuclei from myrxgpool and to let them grow
			theca->solve_RXGROWTH();


			timer[4] = MPI_Wtime();
			prof.MPIWTimeSpendGrowthSim = timer[4] - timer[3];
			prof.JobID = theca->jobid;
			prof.Iterations = theca->step;
			prof.tend = theca->t;
			myCAProfiler.push_back ( prof );

			struct loginfo_ca_physics simphysics;
			theca->log_ca_physics( &simphysics );
			myCAPhysics.push_back ( simphysics );


			if ( this->outopt_rxfront == OUTPUT_RXFRONTSTATS_YES ) {
				theca->write_rxfrontstats();
			}

			if ( this->outopt_singlegrainevo == OUTPUT_SINGLEGRAIN_YES ) {
				theca->write_grainevolution();
			}


			if ( (outopt_hpf3d2donmaster == true) && (myRank == MASTER) ) { 
				theca->write_hpf3d(); 
				theca->write_hpf2d( 0 );
			}

#ifdef FRAGMENTATION_RADAR
			if ( theca->myrxfrontfragmentation.size() > 0 ) {
				theca->write_fragmentation_radar();
			}
#endif

		} //execute all automata on a particular thread
	} //execute ##(in parallel) all of the threadteam
}


void ensembleHdl::postprocess_initrediscrtimes( void )
{
	//acceptable time frame investigated
	double dtworld = ensRediscrWindow.tensmax - ensRediscrWindow.tensmin;


	//construct rediscretized time-annealing scheme
	ensRediscrTime = NULL;
	ensRediscrTime = new double[ensRediscrWindow.nslots];
	QUICKASSERT( ensRediscrTime != NULL );
	ensMemGuard = ensMemGuard + (ensRediscrWindow.nslots * sizeof(double));


	if ( ensRediscrWindow.strategy == REDISCR_TIME_EQUIDISTANT ) {
		double ttmin = ensRediscrWindow.tensmin;
		double ttmax = ensRediscrWindow.tensmax;
		uint32_t nmax = ensRediscrWindow.nslots;

		double dt = (ttmax - ttmin ) / ((double) (nmax - 1));
		QUICKASSERT ( dt > REDISCR_DTMIN );
		
		for ( uint32_t n = 0; n < nmax; ++n ) { //points to end of the interval
			ensRediscrTime[n] = ttmin + ((double) (n + 1) * dt);

			//cout << "->RediscrTime;myRank;n;value\t\t" << myRank << "\t" << n << "\t" << ensRediscrTime[n] << endl;
		}
	}
	//else {
	//	REDISCR_TIME_LOGDISTANT;
	//}
}


#define SENDMESPARA		98
#define SENDRECPARA		99
#define SENDMESPHYS		88
#define SENDRECPHYS		89
#define SENDMESDEFTEX	78
#define SENDRECDEFTEX	79
#define SENDMESRXTEX	68
#define SENDRECRXTEX	69

void ensembleHdl::postprocess_write_mycas( void )
{
	//collect from the workers the profiling information
	MPI_IO_CAProfilingInfoData* bucket = NULL;
	MPI_IO_CAPhysicsInfoData* physbucket = NULL;
	long* deftexture = NULL; //is long to guarantee MPI_LONG running not in range issues if a count is larger than INTEGER_RANGE_MAX
	long* nuctexture = NULL;

	//allow the MASTER to output information when the message of r have arrived
	double prof_tio;

	uint32_t thisrank_nids = this->myIDs.size();
	//##MK::quick and dirty hack nidealcomponents the same for deformation and nuclei texture components
	uint32_t nidealcomponents = DEFAULT_NSTANDARDLAGEN;

	stringstream log_profiling_fname;
	ofstream log_profiling_file;

	//first of all master outputs
	if ( this->myRank == MASTER ) {
		log_profiling_fname << "SCORE." << simid << ".ProfilingLog.csv";
		log_profiling_file.open ( log_profiling_fname.str().c_str() );
		log_profiling_file << "MyRank;JobID;tend;Iterations;MPIWTimeDefMS;MPIWTimeGBDetection;MPIWTimeNucleation;MPIWTimeGrowth;MPIWtimeSinceMasterExitedMPIInit;CAnx;CAny;CAnz;CAnxyz;CASvDeformed;CAStoredEnergy/in1E+14;CAnDefSeeds;CAnRXG;CAnIdeal;IdealComponentsDeformation;IdealComponentsNucleiIdentified\n";

		for ( uint32_t myca = 0; myca < thisrank_nids; myca++ ) {
			prof_tio = MPI_Wtime();

			log_profiling_file << myRank << ";" << myCAProfiler[myca].JobID << ";" << myCAProfiler[myca].tend << ";" << myCAProfiler[myca].Iterations << ";" << setprecision(6) << myCAProfiler[myca].MPIWTimeSpendDefMS << ";" << setprecision(6) << myCAProfiler[myca].MPIWTimeSpendGBDetection << ";" << setprecision(6) << myCAProfiler[myca].MPIWTimeSpendNucleation << ";" << setprecision(6) << myCAProfiler[myca].MPIWTimeSpendGrowthSim;
			log_profiling_file << ";" << setprecision(6) << (prof_tio - prof_t0) << ";" << myCAPhysics[myca].nx << ";" << myCAPhysics[myca].ny << ";" << myCAPhysics[myca].nz << ";" << myCAPhysics[myca].nxyz << ";" << myCAPhysics[myca].nboundarycells << ";" << setprecision(8) << myCAPhysics[myca].storedenergy << ";" << myCAPhysics[myca].ndgrseeds << ";" << myCAPhysics[myca].nrxgrseeds << ";" << myCAPhysics[myca].ndefmicrotexture;

			for ( uint32_t oi = 0; oi < myCAPhysics[myca].ndefmicrotexture; oi++ ) //deformed
				log_profiling_file << ";" << myCAPhysics[myca].defmicrotexture[oi];
			for ( uint32_t oi = myCAPhysics[myca].ndefmicrotexture; oi < nidealcomponents; oi++ )
				log_profiling_file << ";";

			for ( uint32_t ni = 0; ni < myCAPhysics[myca].nnucleimicrotexture; ni++ ) //nuclei
				log_profiling_file << ";" << myCAPhysics[myca].nucleimicrotexture[ni];
			for ( uint32_t ni = myCAPhysics[myca].nnucleimicrotexture; ni < nidealcomponents; ni++ )
				log_profiling_file << ";";

			log_profiling_file << "\n";
		}
	}

	//##MK
	MPI_Barrier( MPI_COMM_WORLD );

	for ( int r = 1; r < this->nRanks; r++ ) {

//if ( myRank == MASTER ) cout << "\t\tpostprocessing results ..." << r << endl;

		thisrank_nids = this->myIDs.size();

		MPI_Bcast( &thisrank_nids, 1, MPI_INT, r, MPI_COMM_WORLD);

//cout << myRank << "\t\t\t" << thisrank_nids << endl;

		if ( r == myRank ) {
			bucket = new MPI_IO_CAProfilingInfoData[thisrank_nids];
			QUICKASSERT ( bucket != NULL );
			physbucket = new MPI_IO_CAPhysicsInfoData[thisrank_nids];
			QUICKASSERT ( physbucket != NULL );
			deftexture = new long[thisrank_nids*nidealcomponents];
			QUICKASSERT ( deftexture != NULL );
			nuctexture = new long[thisrank_nids*nidealcomponents];
			for ( long k = 0; k < (thisrank_nids * nidealcomponents); ++k ) { 
				deftexture[k] = 0; 
				nuctexture[k] = 0;
			}

			//r fills in his data
			for ( uint32_t i = 0; i < thisrank_nids; i++ ) {
				//parallel performance
				bucket[i].JobID = this->myCAProfiler[i].JobID;
				bucket[i].Iterations = this->myCAProfiler[i].Iterations;
				bucket[i].tend = this->myCAProfiler[i].tend;
				bucket[i].MPIWTimeSpendDefMS = this->myCAProfiler[i].MPIWTimeSpendDefMS;
				bucket[i].MPIWTimeSpendGBDetection = this->myCAProfiler[i].MPIWTimeSpendGBDetection;
				bucket[i].MPIWTimeSpendNucleation = this->myCAProfiler[i].MPIWTimeSpendNucleation;
				bucket[i].MPIWTimeSpendGrowthSim = this->myCAProfiler[i].MPIWTimeSpendGrowthSim;

				//physical information
				physbucket[i].nx = this->myCAPhysics[i].nx;
				physbucket[i].ny = this->myCAPhysics[i].ny;
				physbucket[i].nz = this->myCAPhysics[i].nz;
				physbucket[i].nxyz = this->myCAPhysics[i].nxyz;
				physbucket[i].nboundarycells = this->myCAPhysics[i].nboundarycells;
				physbucket[i].ndgrseeds = this->myCAPhysics[i].ndgrseeds;
				physbucket[i].nrxgrseeds = this->myCAPhysics[i].nrxgrseeds;
				physbucket[i].ndefmicrotexture = this->myCAPhysics[i].ndefmicrotexture;
				physbucket[i].nnucleimicrotexture = this->myCAPhysics[i].nnucleimicrotexture;
				physbucket[i].storedenergy = this->myCAPhysics[i].storedenergy;

				//standardlagen deformed grains
				QUICKASSERT( myCAPhysics[i].ndefmicrotexture < nidealcomponents );
				for ( uint32_t j = 0; j < this->myCAPhysics[i].ndefmicrotexture; j++ ) {
					deftexture[(i*nidealcomponents)+j] = this->myCAPhysics[i].defmicrotexture[j];
				}

				//standardlagen rx nuclei
				QUICKASSERT ( myCAPhysics[i].nnucleimicrotexture < nidealcomponents );
				for ( uint32_t nj = 0; nj < this->myCAPhysics[i].nnucleimicrotexture; nj++ ) {
					nuctexture[(i*nidealcomponents)+nj] = this->myCAPhysics[i].nucleimicrotexture[nj];
				}
			}

//cout << r << "\t\t" << myRank << "ifrmyrank-beforesend" << endl;

			MPI_Send( bucket, thisrank_nids, MPI_IO_CAProfilingInfoData_Type, MASTER, SENDMESPARA, MPI_COMM_WORLD);

			MPI_Send( physbucket, thisrank_nids, MPI_IO_CAPhysicsInfoData_Type, MASTER, SENDMESPHYS, MPI_COMM_WORLD );

			MPI_Send( deftexture, thisrank_nids*nidealcomponents, MPI_LONG, MASTER, SENDMESDEFTEX, MPI_COMM_WORLD );

			MPI_Send ( nuctexture, thisrank_nids*nidealcomponents, MPI_LONG, MASTER, SENDMESRXTEX, MPI_COMM_WORLD );

//cout << r << "\t\t" << myRank << "ifrmyrank-aftersend" << endl;

			delete [] bucket;
			bucket = NULL;
			delete [] physbucket;
			physbucket = NULL;
			delete [] deftexture;
			deftexture = NULL;
			delete [] nuctexture;
			nuctexture = NULL;
		}

		if ( myRank == MASTER ) {
			bucket = new MPI_IO_CAProfilingInfoData[thisrank_nids];
			QUICKASSERT ( bucket != NULL );
			physbucket = new MPI_IO_CAPhysicsInfoData[thisrank_nids];
			QUICKASSERT ( physbucket != NULL );
			deftexture = new long[thisrank_nids*nidealcomponents];
			QUICKASSERT ( deftexture != NULL );
			nuctexture = new long[thisrank_nids*nidealcomponents];
			QUICKASSERT ( nuctexture != NULL );
			for ( long k = 0; k < (thisrank_nids * nidealcomponents); ++k ) { 
				deftexture[k] = 0; 
				nuctexture[k] = 0;
			}

//cout << r << "\t\t" << "MASTERifrmyrank-beforerecv" << endl;

			MPI_Recv( bucket, thisrank_nids, MPI_IO_CAProfilingInfoData_Type, r, SENDMESPARA, MPI_COMM_WORLD, MPI_STATUS_IGNORE );

			MPI_Recv( physbucket, thisrank_nids, MPI_IO_CAPhysicsInfoData_Type, r, SENDMESPHYS, MPI_COMM_WORLD, MPI_STATUS_IGNORE );

			MPI_Recv( deftexture, thisrank_nids*nidealcomponents, MPI_LONG, r, SENDMESDEFTEX, MPI_COMM_WORLD, MPI_STATUS_IGNORE );

			MPI_Recv( nuctexture, thisrank_nids*nidealcomponents, MPI_LONG, r, SENDMESRXTEX, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
//cout << r << "\t\t" << "MASTERfrmyrank-afterrecv" << endl;



			for ( uint32_t rca = 0; rca < thisrank_nids; rca++ ) {
				prof_tio = MPI_Wtime();

				log_profiling_file << r << ";" << bucket[rca].JobID << ";" << bucket[rca].tend << ";" << bucket[rca].Iterations << ";" << setprecision(6) << bucket[rca].MPIWTimeSpendDefMS << ";" << setprecision(6) << bucket[rca].MPIWTimeSpendGBDetection << ";" << setprecision(6) << bucket[rca].MPIWTimeSpendNucleation << ";" << setprecision(6) << bucket[rca].MPIWTimeSpendGrowthSim;
				log_profiling_file << ";" << setprecision(6) << (prof_tio - prof_t0) << ";" << physbucket[rca].nx << ";" << physbucket[rca].ny << ";" << physbucket[rca].nz << ";" << physbucket[rca].nxyz << ";" << physbucket[rca].nboundarycells << ";" << physbucket[rca].storedenergy << ";" << physbucket[rca].ndgrseeds << ";" << physbucket[rca].nrxgrseeds << ";" << physbucket[rca].ndefmicrotexture;
				//deftexture
				for ( uint32_t oi = 0; oi < physbucket[rca].ndefmicrotexture; oi++ ) 
					log_profiling_file << ";" << deftexture[(rca*nidealcomponents)+oi];
				for ( uint32_t oi = physbucket[rca].ndefmicrotexture; oi < nidealcomponents; oi++ )
					log_profiling_file << ";";

				//rxtexture
				for ( uint32_t ni = 0; ni < physbucket[rca].nnucleimicrotexture; ni++ ) 
					log_profiling_file << ";" << nuctexture[(rca*nidealcomponents)+ni];
				for ( uint32_t ni = physbucket[rca].nnucleimicrotexture; ni < nidealcomponents; ni++ ) 
					log_profiling_file << ";";

				log_profiling_file << "\n";
			}

			delete [] bucket;
			bucket = NULL;
			delete [] physbucket;
			physbucket = NULL;
			delete [] deftexture;
			deftexture = NULL;
			delete [] nuctexture;
			nuctexture = NULL;
		}

		MPI_Barrier( MPI_COMM_WORLD );
	} //handle next node


	if ( this->myRank == MASTER ) {
		log_profiling_file.flush();
		log_profiling_file.close();
		cout << "The MASTER has output profiling information regarding the simulation." << endl;
	}

	//cout << myRank << " exiting writing post information!" << endl;
}


void ensembleHdl::postprocess_write_mycas_mpifast( void )
{
	//collect from the workers the profiling information, worker can have different number of cellular automata domains!
	double prof_tio = 0.0;
	int thisrank_nids = this->myIDs.size();
	long nidealcomponents = DEFAULT_NSTANDARDLAGEN; //MK::required the same for all ranks!

	//first the workers except the master prepare their data that the master pipes to the file
	MPI_IO_CAProfilingInfoData* bucket = NULL;
	bucket = new MPI_IO_CAProfilingInfoData[thisrank_nids];
	QUICKASSERT ( bucket != NULL );
	MPI_IO_CAPhysicsInfoData* physbucket = NULL;
	physbucket = new MPI_IO_CAPhysicsInfoData[thisrank_nids];
	QUICKASSERT ( physbucket != NULL );
	long* texbucket = NULL; //convention has a size of nidealcomponents*2*howmany CA per process, first all deformed grain components (left part), then all recrystallized (right part)
	texbucket = new long[thisrank_nids*2*nidealcomponents]; 
	QUICKASSERT ( texbucket != NULL );
	for ( long k = 0; k < (thisrank_nids * 2 * nidealcomponents); ++k ) { texbucket[k] = 0; }

	//rank fills in his own data, all ranks in parallel
	if ( this->myRank != MASTER ) {
		for ( uint32_t i = 0; i < thisrank_nids; i++ ) {
			//parallel performance
			bucket[i].JobID = this->myCAProfiler[i].JobID;
			bucket[i].Iterations = this->myCAProfiler[i].Iterations;
			bucket[i].tend = this->myCAProfiler[i].tend;
			bucket[i].MPIWTimeSpendDefMS = this->myCAProfiler[i].MPIWTimeSpendDefMS;
			bucket[i].MPIWTimeSpendGBDetection = this->myCAProfiler[i].MPIWTimeSpendGBDetection;
			bucket[i].MPIWTimeSpendNucleation = this->myCAProfiler[i].MPIWTimeSpendNucleation;
			bucket[i].MPIWTimeSpendGrowthSim = this->myCAProfiler[i].MPIWTimeSpendGrowthSim;

			//physical information
			physbucket[i].nx = this->myCAPhysics[i].nx;
			physbucket[i].ny = this->myCAPhysics[i].ny;
			physbucket[i].nz = this->myCAPhysics[i].nz;
			physbucket[i].nxyz = this->myCAPhysics[i].nxyz;
			physbucket[i].nboundarycells = this->myCAPhysics[i].nboundarycells;
			physbucket[i].ndgrseeds = this->myCAPhysics[i].ndgrseeds;
			physbucket[i].nrxgrseeds = this->myCAPhysics[i].nrxgrseeds;
			physbucket[i].ndefmicrotexture = this->myCAPhysics[i].ndefmicrotexture;
			physbucket[i].nnucleimicrotexture = this->myCAPhysics[i].nnucleimicrotexture;
			physbucket[i].storedenergy = this->myCAPhysics[i].storedenergy;

			//standardlagen deformed grains
			QUICKASSERT( myCAPhysics[i].ndefmicrotexture < nidealcomponents );
			for ( uint32_t dj = 0; dj < this->myCAPhysics[i].ndefmicrotexture; dj++ ) {
				texbucket[(i*2*nidealcomponents)+dj] = this->myCAPhysics[i].defmicrotexture[dj];
			}

			//standardlagen rx nuclei
			QUICKASSERT ( myCAPhysics[i].nnucleimicrotexture < nidealcomponents );
			for ( uint32_t rj = 0; rj < this->myCAPhysics[i].nnucleimicrotexture; rj++ ) {
				texbucket[(i*2*nidealcomponents)+nidealcomponents+rj] = this->myCAPhysics[i].nucleimicrotexture[rj];
			}
		}
	}

	//meanwhile the master can already start I/O to write out his own data...
	stringstream log_profiling_fname;
	ofstream log_profiling_file;

	if ( this->myRank == MASTER ) {
		log_profiling_fname << "SCORE." << simid << ".ProfilingLog.csv";
		log_profiling_file.open ( log_profiling_fname.str().c_str() );
		log_profiling_file << "MyRank;JobID;tend;Iterations;MPIWTimeDefMS;MPIWTimeGBDetection;MPIWTimeNucleation;MPIWTimeGrowth;MPIWtimeSinceMasterExitedMPIInit;CAnx;CAny;CAnz;CAnxyz;CASvDeformed;CAStoredEnergy/in1E+14;CAnDefSeeds;CAnRXG;CAnIdeal;IdealComponentsDeformation;IdealComponentsNucleiIdentified\n";

		for ( uint32_t myca = 0; myca < thisrank_nids; myca++ ) {
			prof_tio = MPI_Wtime();

			log_profiling_file << myRank << ";" << myCAProfiler[myca].JobID << ";" << myCAProfiler[myca].tend << ";" << myCAProfiler[myca].Iterations << ";" << setprecision(6) << myCAProfiler[myca].MPIWTimeSpendDefMS << ";" << setprecision(6) << myCAProfiler[myca].MPIWTimeSpendGBDetection << ";" << setprecision(6) << myCAProfiler[myca].MPIWTimeSpendNucleation << ";" << setprecision(6) << myCAProfiler[myca].MPIWTimeSpendGrowthSim;
			log_profiling_file << ";" << setprecision(6) << (prof_tio - prof_t0) << ";" << myCAPhysics[myca].nx << ";" << myCAPhysics[myca].ny << ";" << myCAPhysics[myca].nz << ";" << myCAPhysics[myca].nxyz << ";" << myCAPhysics[myca].nboundarycells << ";" << setprecision(8) << myCAPhysics[myca].storedenergy << ";" << myCAPhysics[myca].ndgrseeds << ";" << myCAPhysics[myca].nrxgrseeds << ";" << myCAPhysics[myca].ndefmicrotexture;

			for ( uint32_t oi = 0; oi < myCAPhysics[myca].ndefmicrotexture; oi++ ) //deformed
				log_profiling_file << ";" << myCAPhysics[myca].defmicrotexture[oi];
			for ( uint32_t oi = myCAPhysics[myca].ndefmicrotexture; oi < nidealcomponents; oi++ )
				log_profiling_file << ";";

			for ( uint32_t ni = 0; ni < myCAPhysics[myca].nnucleimicrotexture; ni++ ) //nuclei
				log_profiling_file << ";" << myCAPhysics[myca].nucleimicrotexture[ni];
			for ( uint32_t ni = myCAPhysics[myca].nnucleimicrotexture; ni < nidealcomponents; ni++ )
				log_profiling_file << ";";

			log_profiling_file << "\n";
		}
	}

	//communication management, ##MK::be careful fails for more than 2^32-1 CAs in the ensemble!
	int* nids_per_rank = NULL;
	nids_per_rank = new int[this->nRanks];
	QUICKASSERT ( nids_per_rank != NULL );

	int* nids_per_rank_cumul = NULL;
	nids_per_rank_cumul = new int[this->nRanks];
	QUICKASSERT ( nids_per_rank_cumul != NULL );

	//##MK::limits further the maximum admissible size of domains to 2^32-1 / (2*DEFAULT_NSTANDARDLAGEN), which is <=71.5 mio solitary units !!! however very likely sufficient...
	int* tex_per_rank = NULL;
	tex_per_rank = new int[this->nRanks];
	QUICKASSERT( tex_per_rank != NULL);

	int* tex_per_rank_cumul = NULL;
	tex_per_rank_cumul = new int[this->nRanks];
	QUICKASSERT ( tex_per_rank_cumul != NULL );

	for ( uint32_t rr = 0; rr < this->nRanks; ++rr ) {
		nids_per_rank[rr] = 0;
		nids_per_rank_cumul[rr] = 0;
		tex_per_rank[rr] = 0;
		tex_per_rank_cumul[rr] = 0;
	}

	//collector container only significant for the master
	MPI_IO_CAProfilingInfoData* bucket_from_all = NULL;
	MPI_IO_CAPhysicsInfoData* physbucket_from_all = NULL;
	long* texbucket_from_all = NULL; //convention has a size of nidealcomponents*2*howmany CA per process, first all deformed grain components, then all recrystallized

	if ( this->myRank == MASTER ) {
		bucket_from_all = new MPI_IO_CAProfilingInfoData[this->nworldCAs];		QUICKASSERT ( bucket_from_all != NULL);
		physbucket_from_all = new MPI_IO_CAPhysicsInfoData[this->nworldCAs];	QUICKASSERT ( physbucket_from_all != NULL );
		texbucket_from_all = new long[this->nworldCAs * (2*nidealcomponents)];	QUICKASSERT ( texbucket_from_all != NULL );
	}

	//necessary because ranks require waiting for data preparation before I/O
	MPI_Barrier( MPI_COMM_WORLD );

	MPI_Allgather( &thisrank_nids, 1, MPI_INT, nids_per_rank, 1, MPI_INT, MPI_COMM_WORLD );


	int csumnids = 0; //calculate displacement where to place in the master the data from all workers
	int csumtex = 0;
	for ( uint32_t r = 0; r < this->nRanks; r++ ) {
		//nids_per_rank already known...
		tex_per_rank[r] = nids_per_rank[r] * 2 * nidealcomponents;

		nids_per_rank_cumul[r] = csumnids;
		tex_per_rank_cumul[r] = csumtex;

		csumnids = csumnids + nids_per_rank[r];
		csumtex = csumtex + tex_per_rank[r];
	}

//cout << "myRank=1,2,3,4--" << this->myRank << "--" << nids_per_rank[0] << ";" << nids_per_rank[1] << endl;

	MPI_Barrier( MPI_COMM_WORLD );

	//now gatherv, order of elements in recvbuffers is not important because JobID is stored in the profiling information
	MPI_Gatherv( bucket, thisrank_nids, MPI_IO_CAProfilingInfoData_Type, bucket_from_all, nids_per_rank, nids_per_rank_cumul, MPI_IO_CAProfilingInfoData_Type, MASTER, MPI_COMM_WORLD );

	MPI_Barrier ( MPI_COMM_WORLD ); //##MK::debug, safety

	MPI_Gatherv( physbucket, thisrank_nids, MPI_IO_CAPhysicsInfoData_Type, physbucket_from_all, nids_per_rank, nids_per_rank_cumul, MPI_IO_CAPhysicsInfoData_Type, MASTER, MPI_COMM_WORLD );

	MPI_Barrier ( MPI_COMM_WORLD ); //##MK::debug, safety

	MPI_Gatherv( texbucket, (thisrank_nids * 2 * nidealcomponents), MPI_LONG, texbucket_from_all, tex_per_rank, tex_per_rank_cumul, MPI_LONG, MASTER, MPI_COMM_WORLD );

	MPI_Barrier ( MPI_COMM_WORLD );

	//master writes all data into ASCII file, even faster (possibly) --> MPI I/O
	if ( this->myRank == MASTER ) {
		for ( int r = 1; r < this->nRanks; r++ ) {
			for ( uint32_t rca = 0; rca < nids_per_rank[r]; rca++ ) {
				prof_tio = MPI_Wtime();

				log_profiling_file << r << ";" << bucket_from_all[nids_per_rank_cumul[r]+rca].JobID << ";" << bucket_from_all[nids_per_rank_cumul[r]+rca].tend << ";" << bucket_from_all[nids_per_rank_cumul[r]+rca].Iterations << ";" << setprecision(6) << bucket_from_all[nids_per_rank_cumul[r]+rca].MPIWTimeSpendDefMS << ";" << setprecision(6) << bucket_from_all[nids_per_rank_cumul[r]+rca].MPIWTimeSpendGBDetection << ";" << setprecision(6) << bucket_from_all[nids_per_rank_cumul[r]+rca].MPIWTimeSpendNucleation << ";" << setprecision(6) << bucket_from_all[nids_per_rank_cumul[r]+rca].MPIWTimeSpendGrowthSim;
				log_profiling_file << ";" << setprecision(6) << (prof_tio - prof_t0) << ";" << physbucket_from_all[nids_per_rank_cumul[r]+rca].nx << ";" << physbucket_from_all[nids_per_rank_cumul[r]+rca].ny << ";" << physbucket_from_all[nids_per_rank_cumul[r]+rca].nz << ";" << physbucket_from_all[nids_per_rank_cumul[r]+rca].nxyz << ";" << physbucket_from_all[nids_per_rank_cumul[r]+rca].nboundarycells << ";" << physbucket_from_all[nids_per_rank_cumul[r]+rca].storedenergy << ";" << physbucket_from_all[nids_per_rank_cumul[r]+rca].ndgrseeds << ";" << physbucket_from_all[nids_per_rank_cumul[r]+rca].nrxgrseeds << ";" << physbucket_from_all[nids_per_rank_cumul[r]+rca].ndefmicrotexture;

				//deftexture
				for ( uint32_t oi = 0; oi < physbucket_from_all[nids_per_rank_cumul[r]+rca].ndefmicrotexture; oi++ ) 
					log_profiling_file << ";" << texbucket_from_all[(nids_per_rank_cumul[r]*2*nidealcomponents)+(rca*2*nidealcomponents)+oi]; //left part are data for DEFORMED
				for ( uint32_t oi = physbucket_from_all[nids_per_rank_cumul[r]+rca].ndefmicrotexture; oi < nidealcomponents; oi++ ) //add rest of fixed number of columns
					log_profiling_file << ";";

				//rxtexture
				for ( uint32_t ni = 0; ni < physbucket_from_all[nids_per_rank_cumul[r]+rca].nnucleimicrotexture; ni++ ) 
					log_profiling_file << ";" << texbucket_from_all[(nids_per_rank_cumul[r]*2*nidealcomponents)+(rca*2*nidealcomponents)+(nidealcomponents+ni)]; //right part are data for RX
				for ( uint32_t ni = physbucket_from_all[nids_per_rank_cumul[r]+rca].nnucleimicrotexture; ni < nidealcomponents; ni++ ) //add rest of fixed number of columns
					log_profiling_file << ";";

				log_profiling_file << "\n";
			}
		} //handle next node
	}

	//memory bookkeeping
	delete [] bucket;					bucket = NULL;
	delete [] physbucket;				physbucket = NULL;
	delete [] texbucket;				texbucket = NULL;

	delete [] nids_per_rank;			nids_per_rank = NULL;
	delete [] nids_per_rank_cumul;		nids_per_rank_cumul = NULL;
	delete [] tex_per_rank;				tex_per_rank = NULL;
	delete [] tex_per_rank_cumul;		tex_per_rank_cumul = NULL;

	delete [] bucket_from_all;			bucket_from_all = NULL;
	delete [] physbucket_from_all;		physbucket_from_all = NULL;
	delete [] texbucket_from_all;		texbucket_from_all = NULL;


	if ( this->myRank == MASTER ) {
		log_profiling_file.flush();
		log_profiling_file.close();
		cout << "The MASTER has output profiling information regarding the simulation." << endl;
	}

	//cout << myRank << " exiting writing post information!" << endl;
}


void ensembleHdl::postprocess_init( void )
{
	//ensembleHdl first determines maximum local tsimend
	double ensmin_tsimend = INFINITE;
	double ensmax_tsimend = 0.0;
	for ( uint32_t tid = 0; tid < myCAs.size(); tid++ ) {
		for ( uint32_t ca = 0; ca < myCAs[tid].size(); ca++) {
			caHdlP theca = myCAs[tid][ca];

			if ( theca->tsimend <= ensmin_tsimend ) {
				ensmin_tsimend = theca->tsimend;
			}

			if ( theca->tsimend >= ensmax_tsimend ) {
				ensmax_tsimend = theca->tsimend;
			}
		}
	}

	//provides this to the worldpool
	double worldmin_tsimend = 0.0;
	double worldmax_tsimend = 0.0;
	MPI_Allreduce( &ensmin_tsimend, &worldmin_tsimend, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );
	MPI_Allreduce( &ensmax_tsimend, &worldmax_tsimend, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );

	if (myRank == MASTER) { cout << "The overall worldmax_tsimend is between " << worldmin_tsimend << " and " << worldmax_tsimend << endl; }

	//from now on, all ensembleHdl have agreed on the same time rediscretization scheme [ensmin;ensmax]
	ensRediscrWindow.tensmin = 0.0;
	ensRediscrWindow.tensmax = worldmax_tsimend; //##MK::nslots was = 100
	ensRediscrWindow.strategy = REDISCR_TIME_EQUIDISTANT;

	postprocess_initrediscrtimes();
}


//##MK::the postprocessing functions
void ensembleHdl::postprocess_rediscr_kinetics( void )
{
	//thread based interpolation possible //##pragma but pending
	//##MK::strategy cache-friendly when, as here, each grain is analyzed along all times
	uint32_t nmax = ensRediscrWindow.nslots;

	double* rediscr_myCA_allgrains_vol = NULL;
	rediscr_myCA_allgrains_vol = new double[nmax];
	QUICKASSERT ( rediscr_myCA_allgrains_vol != NULL );
	ensMemGuard = ensMemGuard + (nmax * sizeof(double));

	//ensemble collector
	double* rediscr_ensemble_vol = NULL;
	rediscr_ensemble_vol = new double[nmax];
	QUICKASSERT ( rediscr_ensemble_vol != NULL );
	ensMemGuard = ensMemGuard + (nmax * sizeof(double));
	for ( uint32_t n = 0; n < nmax; ++n ) { rediscr_ensemble_vol[n] = 0.0; }

	//ensemble size can be larger than int range so either long or double
	double ensemble_myCAs_nboxvol = 0.0;

//handle data at ensemble level by interpolation analyses on all mycas for all rediscr time steps
	//MK::the idea of the rediscretization scheme is the following: all simulations last up to tend at that point in time,
	//the microstructure transformed completely or partially, in any case the simulation box volume - sum of the volume of cells still in the state of deformed account for the RX volume fraction of the sample
	//no guard zone necessary periodic boundary conditions are applied or their effect on ms evolution to analyze is the scope of the analysis
	for ( int tid = 0; tid < myCAs.size(); tid++ ) {
		for ( int ca = 0; ca < myCAs[tid].size(); ca++) {
			caHdlP theca = myCAs[tid][ca];

			double theca_nboxvol = theca->myCAGeometry.nboxvol_rdndtd;

			//MK::how much volume sample for the myCAs at the ensemble level?
			ensemble_myCAs_nboxvol = ensemble_myCAs_nboxvol + theca_nboxvol;

			//reset temporary collector for the theca to 0.0
			for ( uint32_t n = 0; n < nmax; ++n ) { 
				rediscr_myCA_allgrains_vol[n] = 0.0; 
			}

			//how much volume still left deformed?
			uint32_t theca_ndefg = theca->mydefgpool.size();

			for ( uint32_t dg = 0; dg < theca_ndefg; ++dg ) {
				for ( uint32_t n = 0; n < nmax; ++n ) {
					//interpolate volume of the dg grain in the mydefgpool of theca at that rediscretized time slot
					rediscr_myCA_allgrains_vol[n] = rediscr_myCA_allgrains_vol[n] + theca->get_interpCellCount( dg, ensRediscrTime[n] );
				}
			}

			//pipe this information to the local ensembleLevel
			//utilizing that always Sum(Vdef) + Sum(Vrxg) == Vnboxvol_rdndtd holds valid
			for ( uint32_t n = 0; n < nmax; ++n ) {
				rediscr_ensemble_vol[n] = rediscr_ensemble_vol[n] + ( theca_nboxvol - rediscr_myCA_allgrains_vol[n] );
			}

			//rediscr_ensemble is constantly added to but rediscr_myCA is switched back to 0.0 and recycled

		} //all CAs on one thread
	} //for all thread -> myCAs

	//cleanup temporaries
	delete [] rediscr_myCA_allgrains_vol;
	rediscr_myCA_allgrains_vol = NULL;
	ensMemGuard = ensMemGuard - (nmax * sizeof(double));

	//safer to introduce an MPI_Barrier if the MPI_Reduces are called significantly delayed among the processes owing to load imbalances
	MPI_Barrier( MPI_COMM_WORLD );

//handle data at world level output only relevant for MASTER who does I/O
	double world_allCAs_nboxvol = 0.0;
	MPI_Reduce ( &ensemble_myCAs_nboxvol, &world_allCAs_nboxvol, 1, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD );

	double* rediscr_world_vol = NULL;
	rediscr_world_vol = new double[nmax]; //nmax is the same in each process
	QUICKASSERT( rediscr_world_vol != NULL );
	ensMemGuard = ensMemGuard + (nmax * sizeof(double));
	for ( uint32_t n = 0; n < nmax; ++n ) { rediscr_world_vol[n] = 0.0; }

	MPI_Reduce( rediscr_ensemble_vol, rediscr_world_vol, nmax, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD );


//MASTER I/O
	if (myRank == MASTER) {
		stringstream log_kinetics_fname;
		ofstream log_kinetics_file;

		log_kinetics_fname << "SCORE." << simid << ".Rediscretized.Kinetics.csv";
		//cout << "File " << log_kinetics_fname.str().c_str() << " is opened now" << endl;
		log_kinetics_file.open ( log_kinetics_fname.str().c_str() );

		//header
		log_kinetics_file << "Step;RediscrTime/s;SumRediscrVolume--" << (world_allCAs_nboxvol) << "--;X;lnt;ln(-ln(1-X))\n";

		double _vref, rxfrac;
		QUICKASSERT ( world_allCAs_nboxvol > 0.0 );
		_vref = ((1.0) / world_allCAs_nboxvol);

		for ( uint32_t n = 0; n < nmax; ++n ) {
			rxfrac = rediscr_world_vol[n] * _vref;

			log_kinetics_file << n << ";" << ensRediscrTime[n] << ";" << rediscr_world_vol[n] << ";" << rxfrac;

			if ( ensRediscrTime[n] > SMALLEST_TIME_LOGARITMIZE )
				log_kinetics_file << ";" << log(ensRediscrTime[n]);
			else
				log_kinetics_file << ";";


			if ( (1.0 - rxfrac) > SMALLEST_TIME_LOGARITMIZE )
				log_kinetics_file << ";" << log(-1.0 * log( 1.0 - rxfrac ));
			else
				log_kinetics_file << ";";

			log_kinetics_file << "\n";
		}

		log_kinetics_file.flush();
		log_kinetics_file.close();
	}

	delete [] rediscr_ensemble_vol;
	rediscr_ensemble_vol = NULL;
	ensMemGuard = ensMemGuard - (nmax * sizeof(double));

	delete [] rediscr_world_vol;
	rediscr_world_vol = NULL;
	ensMemGuard = ensMemGuard - (nmax * sizeof(double));
}


void ensembleHdl::postprocess_rediscr_macrotexture( void )
{
	//thread based interpolation possible //##pragma but pending
	//##MK::strategy cache-friendly if each grain is analyzed along all times
	//macrotexture deformed and recrystallized grains
	uint32_t nmax = ensRediscrWindow.nslots;
	uint32_t nstandardlagen = standardlagen.size() + 1; //for RANDOM_ORIENTATION

	double* rediscr_myCA_allgrains_oridepvol = NULL;
	rediscr_myCA_allgrains_oridepvol = new double[nstandardlagen*nmax];
	QUICKASSERT ( rediscr_myCA_allgrains_oridepvol != NULL );
	ensMemGuard = ensMemGuard + (nstandardlagen * nmax * sizeof(double));

	//ensemble collector
	double* rediscr_ensemble_oridepvol = NULL;
	rediscr_ensemble_oridepvol = new double[nstandardlagen*nmax];
	QUICKASSERT ( rediscr_ensemble_oridepvol != NULL );
	ensMemGuard = ensMemGuard + (nstandardlagen * nmax * sizeof(double));
	for ( long n = 0; n < (nstandardlagen*nmax); ++n ) { rediscr_ensemble_oridepvol[n] = 0.0; }

	//ensemble size can be larger than int range so double
	double ensemble_myCAs_nboxvol = 0.0;

//handle data at ensemble level by interpolation analyses on all mycas for all rediscr time steps
	for ( uint32_t tid = 0; tid < myCAs.size(); tid++ ) {
		for ( uint32_t ca = 0; ca < myCAs[tid].size(); ca++) {
			caHdlP theca = myCAs[tid][ca];

			double theca_nboxvol = theca->myCAGeometry.nboxvol_rdndtd;
			//####MK::MIND FINAL RECRYSTALLIZED FRACTION
			ensemble_myCAs_nboxvol = ensemble_myCAs_nboxvol + theca_nboxvol;

			uint32_t theca_ndefg = theca->mydefgpool.size();
			uint32_t theca_nrxg = theca->myrxgpool.size();

			uint32_t nideal;
			uint32_t idx;
			//reset temporary collector for the theca to 0.0
			for ( long n = 0; n < (nstandardlagen*nmax); ++n ) { 
				rediscr_myCA_allgrains_oridepvol[n] = 0.0;
			}

			//scan all grains for theca
			for ( uint32_t dg = 0; dg < theca_ndefg; ++dg ) {
				//categorize according to standardlagen, catch RANDOM as the first element
				nideal = theca->myoripool[theca->mydefgpool[dg].caori].closestideal;

				//to which ideal orientation is dg belonging?
				for ( uint32_t n = 0; n < nmax; ++n ) {
					//interpolate volume of the dg grain in the mydefgpool of theca at that rediscretized time slot
					rediscr_myCA_allgrains_oridepvol[(nstandardlagen*n)+nideal] += theca->get_interpCellCount( dg, ensRediscrTime[n] );
				}
			}

			for ( uint32_t rxg = 0; rxg < theca_nrxg; ++rxg ) {
				//categorize according to standardlagen, catch random as the first element
				nideal = theca->myoripool[theca->myrxgpool[rxg].caori].closestideal;

				idx = theca_ndefg + rxg;

				for ( uint32_t n = 0; n < nmax; ++n ) {
					rediscr_myCA_allgrains_oridepvol[(nstandardlagen*n)+nideal] += theca->get_interpCellCount( idx, ensRediscrTime[n] );
				}
			}

			//pipe this information to the local ensembleLevel
			//contrary to kinetics here all texture information counts
			for ( long n = 0; n < (nstandardlagen*nmax); ++n ) {
				rediscr_ensemble_oridepvol[n] += rediscr_myCA_allgrains_oridepvol[n];
			}

			//rediscr_ensemble is reutilized, rediscr_myCA is switched back to 0.0 and recycled

		} //all CAs on one thread
	} //for all thread -> myCAs

	//cleanup temporaries
	delete [] rediscr_myCA_allgrains_oridepvol;
	rediscr_myCA_allgrains_oridepvol = NULL;
	ensMemGuard = ensMemGuard - (nstandardlagen * nmax * sizeof(double));

	//safer to introduce an MPI_Barrier if the MPI_Reduces are called significantly delayed among the processes owing to load imbalances
	MPI_Barrier( MPI_COMM_WORLD );

//handle data at world level output only relevant for MASTER who does I/O
	double world_allCAs_nboxvol = 0.0;
	MPI_Reduce ( &ensemble_myCAs_nboxvol, &world_allCAs_nboxvol, 1, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD );

	double* rediscr_world_oridepvol = NULL;
	rediscr_world_oridepvol = new double[nstandardlagen*nmax]; //##MK::nmax is the same in each process, so is the number of ideal components
	QUICKASSERT( rediscr_world_oridepvol != NULL );
	ensMemGuard = ensMemGuard + (nstandardlagen * nmax * sizeof(double));
	for ( long n = 0; n < (nstandardlagen*nmax); ++n ) { rediscr_world_oridepvol[n] = 0.0; }

	MPI_Reduce( rediscr_ensemble_oridepvol, rediscr_world_oridepvol, nstandardlagen*nmax , MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD );

//MASTER I/O
	if (myRank == MASTER) {
		stringstream log_texture_fname;
		ofstream log_texture_file;

		log_texture_fname << "SCORE." << simid << ".Rediscretized.Macrotexture.csv";
		//cout << "File " << log_texture_fname.str().c_str() << " is opened now" << endl;
		log_texture_file.open ( log_texture_fname.str().c_str() );

		//header
		log_texture_file << "Step;RediscrTime/s;RANDOM-phi1-PHI-phi2";
		
		for ( uint32_t id = 0; id < standardlagen.size(); ++id) { 
			log_texture_file << ";" << id << "-" << (standardlagen[id].bunge1 / _PI_ * 180.0) << "-" << (standardlagen[id].bunge2  / _PI_ * 180.0) << "-" << (standardlagen[id].bunge3  / _PI_ * 180.0);
		}
		log_texture_file << ";SumRediscrVolume=" << (world_allCAs_nboxvol) << endl;

		//pipe information to file
		for ( uint32_t n = 0; n < nmax; ++n ) {
			log_texture_file << n << ";" << ensRediscrTime[n];
			for ( uint32_t id = 0; id < nstandardlagen; ++id ) { log_texture_file << ";" << rediscr_world_oridepvol[(nstandardlagen*n)+id]; }
			log_texture_file << ";" << endl;
		}

		log_texture_file.flush();
		log_texture_file.close();
	}

	delete [] rediscr_ensemble_oridepvol;
	rediscr_ensemble_oridepvol = NULL;
	ensMemGuard = ensMemGuard - (nstandardlagen * nmax * sizeof(double));
	delete [] rediscr_world_oridepvol;
	rediscr_world_oridepvol = NULL;
	ensMemGuard = ensMemGuard + (nstandardlagen * nmax * sizeof(double));
}


#define SENDINGGSD		56

void ensembleHdl::postprocess_rediscr_finalgrainsizedistribution( void )
{
	if (myRank == MASTER ) { cout << "Sequential FinalGrainSizeDistribution collection." << endl; }
	uint32_t nmax = ensRediscrWindow.nslots;
	double now = ensRediscrWindow.tensmax;

	//each rank collects volume from all his nuclei at ensRediscrWindow.tensmax in all his automata,
	//MK::if tsimend locally is < tensmax the final volume after completion of the local simulation is taken 
	//MK::thus, assuming negligible volume change due to grain-growth
	unsigned long ensemble_myCAs_nnuclei = 0;
	if ( this->ensNucleationModel.tincubmodel == TINCUB_SITESATURATION ) {
		for ( uint32_t tid = 0; tid < myCAs.size(); tid++ ) {
			for ( uint32_t ca = 0; ca < myCAs[tid].size(); ca++) {
				caHdlP theca = myCAs[tid][ca];
				ensemble_myCAs_nnuclei += theca->myrxgpool.size();
			}
		}
	}
	else if ( this->ensNucleationModel.tincubmodel == TINCUB_TIMEDEPENDENT ) {
		for ( uint32_t tid = 0; tid < myCAs.size(); tid++ ) {
			for ( uint32_t ca = 0; ca < myCAs[tid].size(); ca++) {
				caHdlP theca = myCAs[tid][ca];
				uint32_t nrxg = theca->myrxgpool.size();
				for ( uint32_t nuc = 0; nuc < nrxg; nuc++ ) {
					if ( theca->myrxgpool[nuc].nucsite == NUCLEUS_ALREADY_PLACED ) //MK::ALREADY_CONSUMED excludes in the output grains that were planned but came never into existence and thus would alter the size distribution
						ensemble_myCAs_nnuclei++;
				}
			}
		}
	}
	else { cout << "ERR::Unknown nucleation time model!" << endl; return; }

	//handle myCAs at the ensemble level
	MPI_IO_FinalGSDInfoData* rediscr_ens_allnuclei = NULL;
	rediscr_ens_allnuclei = new MPI_IO_FinalGSDInfoData[ensemble_myCAs_nnuclei];
	QUICKASSERT ( rediscr_ens_allnuclei != NULL );
	ensMemGuard = ensMemGuard + (ensemble_myCAs_nnuclei * sizeof(MPI_IO_FinalGSDInfoData));

	//identify the orientation class and the volume of each nucleus from the process
	long nid = 0; //MK::many more nuclei possible than INTEGER_RANGE_MAX

	if ( this->ensNucleationModel.tincubmodel == TINCUB_SITESATURATION ) {
		for ( uint32_t tid = 0; tid < myCAs.size(); tid++ ) {
			for ( uint32_t ca = 0; ca < myCAs[tid].size(); ca++) {
				caHdlP theca = myCAs[tid][ca];

				uint32_t ndefg = theca->mydefgpool.size();
				uint32_t nrxg = theca->myrxgpool.size();
				uint32_t ndisjoint_grains = ndefg + nrxg;

				uint32_t nideal;
				for ( uint32_t rxg = ndefg; rxg < ndisjoint_grains; rxg++) {
					nideal = theca->myoripool[theca->myrxgpool[rxg-ndefg].caori].closestideal;
					rediscr_ens_allnuclei[nid].finalvol = theca->get_interpCellCount( rxg, now );
					rediscr_ens_allnuclei[nid].tincub = theca->myrxgpool[rxg-ndefg].tincub;
					rediscr_ens_allnuclei[nid].ideal = nideal;
					nid++;
				} //because each grains was nucleated
			}
		}
	}
	else if ( this->ensNucleationModel.tincubmodel == TINCUB_TIMEDEPENDENT ) {
		for ( uint32_t tid = 0; tid < myCAs.size(); tid++ ) {
			for ( uint32_t ca = 0; ca < myCAs[tid].size(); ca++) {
				caHdlP theca = myCAs[tid][ca];

				uint32_t ndefg = theca->mydefgpool.size();
				uint32_t nrxg = theca->myrxgpool.size();
				uint32_t ndisjoint_grains = ndefg + nrxg;

				uint32_t nideal;
				for ( uint32_t rxg = ndefg; rxg < ndisjoint_grains; rxg++) {
					if ( theca->myrxgpool[rxg-ndefg].nucsite == NUCLEUS_ALREADY_PLACED ) {
						nideal = theca->myoripool[theca->myrxgpool[rxg-ndefg].caori].closestideal;
						rediscr_ens_allnuclei[nid].finalvol = theca->get_interpCellCount( rxg, now );
						rediscr_ens_allnuclei[nid].tincub = theca->myrxgpool[rxg-ndefg].tincub;
						rediscr_ens_allnuclei[nid].ideal = nideal;
						nid++;
					} //because probably not all nuclei were placed and thus many zero sized grains would distort the grain size distribution!
				}
			}
		}
	}
	else { cout << "ERR::Unknown nucleation time model!" << endl; return; }


	long* nnuclei_planned_world = NULL;
	long* nnuclei_seeded_world = NULL;
	if (myRank == MASTER) {
		nnuclei_planned_world = new long[nRanks];
		QUICKASSERT ( nnuclei_planned_world != NULL );
		ensMemGuard = ensMemGuard + (nRanks * sizeof(long));
		nnuclei_seeded_world = new long[nRanks];
		QUICKASSERT ( nnuclei_seeded_world != NULL );
		ensMemGuard = ensMemGuard + (nRanks * sizeof(long));
	}

	MPI_Gather( &ensemble_myCAs_nnuclei, 1, MPI_LONG, nnuclei_planned_world, 1, MPI_LONG, MASTER, MPI_COMM_WORLD );
	MPI_Gather( &nid, 1, MPI_LONG, nnuclei_seeded_world, 1, MPI_LONG, MASTER, MPI_COMM_WORLD );

	//MASTER outputs sequentially all grains
	long nnuclei_allranks_world = 0;

	//as it is declared as default it utilized the default allocator class wrapping about new
	//so as desired finalgsd lives on the heap and becomes destructed by the vector destructor when exiting the function
	vector<agrain> finalgsd;
	double finalgsdsum = 0.0;
	double finalgsdmean = 0.0;
	double finalgsdmedian = 0.0;
	double finalgsdstddev = 0.0;
	double finalgsdmin = -1.0;
	double finalgsdmax = INFINITE;

	if ( this->myRank == MASTER ) {
		//how many nuclei the master has to aggregate?
		for ( int r = 0; r < this->nRanks; ++r ) {
			nnuclei_allranks_world += nnuclei_seeded_world[r];
		}

cout << "Master outputs his own results nid=" << nid << endl;

		//collect own results first
		for ( long n = 0; n < nid; n++) { 
			finalgsdsum = finalgsdsum + rediscr_ens_allnuclei[n].finalvol;
			struct agrain gr;
			gr.vol = rediscr_ens_allnuclei[n].finalvol;
			gr.tincub = rediscr_ens_allnuclei[n].tincub;
			gr.ideal = rediscr_ens_allnuclei[n].ideal;
			gr.rank = MASTER;
			finalgsd.push_back( gr );
		}
	}

	MPI_Barrier( MPI_COMM_WORLD );

	long nrecv = 0;
	long nrecvseeded = 0;
	MPI_IO_FinalGSDInfoData* recvbuf = NULL;

	for ( int r = 1; r < nRanks; r++ ) {
		if (myRank == r ) { 
			MPI_Send( rediscr_ens_allnuclei, ensemble_myCAs_nnuclei, MPI_IO_FinalGSDInfoData_Type, MASTER, SENDINGGSD, MPI_COMM_WORLD );
		}
		//MASTER always collects
		if (myRank == MASTER) { 
			nrecv = nnuclei_planned_world[r];
			nrecvseeded = nnuclei_seeded_world[r];

			recvbuf = new MPI_IO_FinalGSDInfoData[nrecv];

			MPI_Recv( recvbuf, nrecv, MPI_IO_FinalGSDInfoData_Type, r, SENDINGGSD, MPI_COMM_WORLD, MPI_STATUS_IGNORE );

			for ( long n = 0; n < nrecvseeded; n++) { //incorporate results from rank r
				finalgsdsum = finalgsdsum + recvbuf[n].finalvol;
				struct agrain gr;
				gr.vol = recvbuf[n].finalvol;
				gr.tincub = recvbuf[n].tincub;
				gr.ideal = recvbuf[n].ideal;
				gr.rank = r;
				finalgsd.push_back( gr );
			}

			delete [] recvbuf;		recvbuf = NULL;
		}

		MPI_Barrier(MPI_COMM_WORLD);
	} //collect grains over all ranks, think about further improvement with MPI_Gatherv


	//I/O
	stringstream log_finalgsd_fname;
	ofstream log_finalgsd_file;

	if ( myRank == MASTER ) {
		//sort ascending in order to calculate the median and to construct probability plots directly
		std::sort( finalgsd.begin(), finalgsd.end(), SortGSDAscending );

		//determine average grain size
		finalgsdmean = (finalgsdsum / (double) nnuclei_allranks_world);

		finalgsdmedian = -1.0;
		//determine median grain size
		if ( finalgsd.size() > 3) {
			size_t nnn = finalgsd.size() / 2;
			//std::nth_element( finalgsd.begin(), finalgsd.begin()+nnn, finalgsd.end() ); not necessary because vector has already been strictly sorted
			if (finalgsd.size() % 2 == 1) { //odd? -> take the center
				finalgsdmedian = finalgsd[nnn].vol;
			} 
			else { //even --> take the linear average in the middle, //std::nth_element( vol2median.begin(), vol2median.begin()+nnn-1, vol2median.end() ); same story...
				finalgsdmedian = 0.5*(finalgsd[nnn].vol + finalgsd[nnn-1].vol);
			}
		}

		QUICKASSERT( nnuclei_allranks_world == finalgsd.size() );
		QUICKASSERT( finalgsdmean > DOUBLE_ACCURACY );
		//Master has collected all grains, so descriptive statistics in C/C++ is faster than later in Origin and co...
		//standard deviation, min and max
		size_t nfgsd = finalgsd.size();
		if ( nfgsd > 0 ) {
			double sum = 0.0;
			double v;
			for ( uint32_t g = 0; g < nfgsd; g++ ) {
				v = finalgsd[g].vol;
				sum = sum + SQR(v - finalgsdmean);
				if ( v <= finalgsdmin ) finalgsdmin = v;
				if ( v >= finalgsdmax ) finalgsdmax = v;
			}

			finalgsdstddev = pow( (sum / ((double) nfgsd) ), 0.5 );
		}

		//open output file and create header
		log_finalgsd_fname << "SCORE." << simid << ".Rediscretized.FinalGSD.csv";
		log_finalgsd_file.open ( log_finalgsd_fname.str().c_str() );
		//write file header
		log_finalgsd_file << "Time=" << now << "s,NGrains=" << nnuclei_allranks_world << ",Vref=" << finalgsdsum << ",Mean=" << finalgsdmean << ",Median=" << finalgsdmedian << ",StdDev=" << finalgsdstddev << ",Min=" << finalgsdmin << ",Max=" << finalgsdmax << "\n";
		log_finalgsd_file << "Rank;Volume/Cells;IncubationTime/s;Ideal(RANDOMis0);ln(Vi/Vav);CumSumi/Vref*100%\n";

		double cumsum = 0.0;
		double volnorm = 0.0;
		double probability = 0.0;
		//output probability plot raw data and final grain size distribution
		for ( unsigned long n = 0; n < nnuclei_allranks_world; ++n ) { 
			cumsum = cumsum + finalgsd[n].vol;
			volnorm = finalgsd[n].vol / finalgsdmean;
			probability = cumsum / finalgsdsum * 100.0;

			log_finalgsd_file << finalgsd[n].rank << ";" << finalgsd[n].vol << ";" << finalgsd[n].tincub << ";" << finalgsd[n].ideal;

			if ( volnorm > SMALLEST_TIME_LOGARITMIZE )	log_finalgsd_file << ";" << (log(volnorm)) << ";" << probability << "\n";
			else										log_finalgsd_file << ";;\n";
		}

		log_finalgsd_file.flush();
		log_finalgsd_file.close();


		delete [] nnuclei_planned_world;
		nnuclei_planned_world = NULL;
		ensMemGuard = ensMemGuard - (nRanks * sizeof(long)); 
		delete [] nnuclei_seeded_world;
		nnuclei_seeded_world = NULL;
		ensMemGuard = ensMemGuard - (nRanks * sizeof(long)); 
	}

	delete [] rediscr_ens_allnuclei;
	rediscr_ens_allnuclei = NULL;
	ensMemGuard = ensMemGuard - (ensemble_myCAs_nnuclei * sizeof(MPI_IO_FinalGSDInfoData));
}


bool ensembleHdl::postprocess_rediscr_finalgrainsizedistribution_mpifast( void )
{
	if (myRank == MASTER ) { cout << "Fast FinalGrainSizeDistribution collection." << endl; }
	//IT IS NECESSARY A LARGE ENOUGH BUFFER IN THE MASTER TO COLLECT ALL GRAINS AT ONCE!
	uint32_t nmax = ensRediscrWindow.nslots;
	double now = ensRediscrWindow.tensmax;

	//each rank collects volume from all his nuclei at ensRediscrWindow.tensmax in all his automata,
	//MK::if tsimend locally is < tensmax the final volume after completion of the local simulation is taken 
	//MK::thus, assuming negligible volume change due to grain-growth
	int ensemble_myCAs_nnuclei = 0;
	if ( this->ensNucleationModel.tincubmodel == TINCUB_SITESATURATION ) {
		for ( uint32_t tid = 0; tid < myCAs.size(); tid++ ) {
			for ( uint32_t ca = 0; ca < myCAs[tid].size(); ca++) {
				caHdlP theca = myCAs[tid][ca];
				ensemble_myCAs_nnuclei += theca->myrxgpool.size();
			}
		}
	}
	else if ( this->ensNucleationModel.tincubmodel == TINCUB_TIMEDEPENDENT ) {
		for ( uint32_t tid = 0; tid < myCAs.size(); tid++ ) {
			for ( uint32_t ca = 0; ca < myCAs[tid].size(); ca++) {
				caHdlP theca = myCAs[tid][ca];
				uint32_t nrxg = theca->myrxgpool.size();
				for ( uint32_t nuc = 0; nuc < nrxg; nuc++ ) {
					if ( theca->myrxgpool[nuc].nucsite == NUCLEUS_ALREADY_PLACED ) //MK::ALREADY_CONSUMED excludes in the output grains that were planned but came never into existence and thus would alter the size distribution
						ensemble_myCAs_nnuclei++;
				}
			}
		}
	}
	else { cout << "ERR::Unknown nucleation time model!" << endl; return false; }

	//communication management, ##MK::be careful limits maximum number of grains to 2^31 - 1, 2.sth billion nuclei!
	int* nnuc_per_rank = NULL;
	nnuc_per_rank = new int[this->nRanks];
	QUICKASSERT ( nnuc_per_rank != NULL );
	int* nnuc_per_rank_cumul = NULL;
	nnuc_per_rank_cumul = new int[this->nRanks];
	QUICKASSERT ( nnuc_per_rank_cumul != NULL );
	for ( int r = 0; r < this->nRanks; r++ ) {
		nnuc_per_rank[r] = 0;
		nnuc_per_rank_cumul[r] = 0;
	}

	MPI_Barrier( MPI_COMM_WORLD );

	MPI_Allgather( &ensemble_myCAs_nnuclei, 1, MPI_INT, nnuc_per_rank, 1, MPI_INT, MPI_COMM_WORLD );

	//MPI_Barrier( MPI_COMM_WORLD );
//cout << "myRank=1,2,3,4--" << this->myRank << "--" << nnuc_per_rank[0] << ";" << nnuc_per_rank[1] << endl;


	//maximum number of nuclei should not be too large that MPI buffers are exceeded..., ##MK::currently working only for at most 2^31 - 1 nuclei, >2 billion nuclei
	int ensemble_all_nnuclei = 0;
	for ( int r = 0; r < this->nRanks; r++ ) { ensemble_all_nnuclei += nnuc_per_rank[r]; }
	if ( ensemble_all_nnuclei > MAXIMUM_NUMBER_OF_NUCLEI_MPI_BUFFERED ) { return false; }

	//not returned, then collect from all ranks directly the MPI_IO_FGSDComplete objects into one vector that the master sorts and outputs

	//analyze cumulative distribution of nnuc_per_rank
	int csum = 0;
	for ( int r = 0; r < this->nRanks; r++ ) {
		nnuc_per_rank_cumul[r] = csum;
		csum = csum + nnuc_per_rank[r];
	}

	MPI_IO_FGSDComplete* grbuf_rank = NULL;

	//if ( this->myRank != MASTER ) {
		grbuf_rank = new MPI_IO_FGSDComplete[ensemble_myCAs_nnuclei];
		QUICKASSERT ( grbuf_rank != NULL );
		ensMemGuard = ensMemGuard + (ensemble_myCAs_nnuclei * sizeof(MPI_IO_FGSDComplete));

		int grb = 0;
		if ( this->ensNucleationModel.tincubmodel == TINCUB_SITESATURATION ) {
			for ( uint32_t tid = 0; tid < myCAs.size(); tid++ ) {
				for ( uint32_t ca = 0; ca < myCAs[tid].size(); ca++) {
					caHdlP theca = myCAs[tid][ca];

					uint32_t ndefg = theca->mydefgpool.size();
					uint32_t nrxg = theca->myrxgpool.size();
					uint32_t ndisjoint_grains = ndefg + nrxg;

					uint32_t nideal;
					int myr = this->myRank;
					for ( uint32_t rxg = ndefg; rxg < ndisjoint_grains; rxg++) {
						nideal = theca->myoripool[theca->myrxgpool[rxg-ndefg].caori].closestideal;
						grbuf_rank[grb].finalvol = theca->get_interpCellCount( rxg, now );
						grbuf_rank[grb].tincub = theca->myrxgpool[rxg-ndefg].tincub;
						grbuf_rank[grb].ideal = nideal;
						grbuf_rank[grb].rank = myr;
						grb++;
					} //because each grains was nucleated
				}
			}
		}
		else if ( this->ensNucleationModel.tincubmodel == TINCUB_TIMEDEPENDENT ) {
			for ( uint32_t tid = 0; tid < myCAs.size(); tid++ ) {
				for ( uint32_t ca = 0; ca < myCAs[tid].size(); ca++) {
					caHdlP theca = myCAs[tid][ca];

					uint32_t ndefg = theca->mydefgpool.size();
					uint32_t nrxg = theca->myrxgpool.size();
					uint32_t ndisjoint_grains = ndefg + nrxg;

					uint32_t nideal;
					int myr = this->myRank;
					for ( uint32_t rxg = ndefg; rxg < ndisjoint_grains; rxg++) {
						if ( theca->myrxgpool[rxg-ndefg].nucsite == NUCLEUS_ALREADY_PLACED ) {
							nideal = theca->myoripool[theca->myrxgpool[rxg-ndefg].caori].closestideal;
							grbuf_rank[grb].finalvol = theca->get_interpCellCount( rxg, now );
							grbuf_rank[grb].tincub = theca->myrxgpool[rxg-ndefg].tincub;
							grbuf_rank[grb].ideal = nideal;
							grbuf_rank[grb].rank = myr;
							grb++;
						} //because probably not all nuclei were placed and thus many zero sized grains would distort the grain size distribution!
					}
				}
			}
		}
		else { cout << "ERR::Unknown nucleation time model!" << endl; return false; }
	//}
	QUICKASSERT ( grb == ensemble_myCAs_nnuclei );


	//meanwhile the master initializes a large vector to keep all the grains and their properties
	vector<MPI_IO_FGSDComplete> finalgsd;

	if ( this->myRank == MASTER ) {
		finalgsd.reserve( ensemble_all_nnuclei ); //cout << "Master allocates space for a large vector to store the grains" << endl;

		MPI_IO_FGSDComplete g;
		g.finalvol = 0.0;		g.tincub = 0.0;		g.ideal = 0;	g.rank = 0;

		for ( unsigned int k = 0; k < ensemble_all_nnuclei; ++k ) { 
			finalgsd.push_back( g ); 
		}
	} //now finalgsd[] is safe to use to drop the grain by a MPI_Gatherv call

//cout << "myRank/ensemble_myCAs_nnuclei/ensemble_all_nnuclei/finalgsd.size/finalgsd.capacity = " << this->myRank << ";" << ensemble_myCAs_nnuclei << ";" << ensemble_all_nnuclei << ";" << finalgsd.size() << ";" << finalgsd.capacity() << endl;

	//as it is declared as default it utilized the default allocator class wrapping about new
	//so as desired finalgsd lives on the heap and becomes destructed by the vector destructor when exiting the function

	MPI_Barrier ( MPI_COMM_WORLD );

	//mind that MPI_Gatherv evaluates to MPI_INT, ##MK so a single process cannot send more than 2^32 - 1 grains safely... however, no simulation with the SCORE model
	//will have as many nuclei on one process, because 4 billions grains to discretize requires at least a factor of 1000 more cells, which no process however has memory for to discretize the cellular automaton grids for
	MPI_Gatherv( grbuf_rank, nnuc_per_rank[this->myRank], MPI_IO_FGSDComplete_Type, &(finalgsd[0]), nnuc_per_rank, nnuc_per_rank_cumul, MPI_IO_FGSDComplete_Type, MASTER, MPI_COMM_WORLD );

	//all grains from all processes ready to output
	MPI_Barrier( MPI_COMM_WORLD );


	//I/O
	stringstream log_finalgsd_fname;
	ofstream log_finalgsd_file;

	double finalgsdsum = 0.0;
	double finalgsdmean = 0.0;
	double finalgsdmedian = 0.0;
	double finalgsdstddev = 0.0;
	double finalgsdmin = -1.0;
	double finalgsdmax = INFINITE;

	if ( myRank == MASTER ) {
		//sort ascending in order to calculate the median and to construct probability plots directly
		std::sort( finalgsd.begin(), finalgsd.end(), SortFGSDCompleteAsc );

		size_t nfgsd = finalgsd.size();

		for ( unsigned int i = 0; i < nfgsd; i++ ) {
			finalgsdsum += finalgsd[i].finalvol;
		}

		if ( nfgsd > 3 ) { //descriptive stats
			finalgsdmean = (finalgsdsum / (double) nfgsd );

			size_t nnn = nfgsd / 2;
			//std::nth_element( finalgsd.begin(), finalgsd.begin()+nnn, finalgsd.end() ); not necessary because vector has already been strictly sorted
			if (nfgsd % 2 == 1) { //odd? -> take the center
				finalgsdmedian = finalgsd[nnn].finalvol;
			} 
			else { //even --> take the linear average in the middle, //std::nth_element( vol2median.begin(), vol2median.begin()+nnn-1, vol2median.end() ); same story...
				finalgsdmedian = 0.5*(finalgsd[nnn].finalvol + finalgsd[nnn-1].finalvol);
			}

			QUICKASSERT( finalgsdmean > DOUBLE_ACCURACY );
			double sum = 0.0;
			for ( unsigned long n = 0; n < nfgsd; n++ ) {
				sum = sum + SQR(finalgsd[n].finalvol - finalgsdmean);
				if ( finalgsd[n].finalvol <= finalgsdmin ) finalgsdmin = finalgsd[n].finalvol;
				if ( finalgsd[n].finalvol >= finalgsdmax ) finalgsdmax = finalgsd[n].finalvol;
			}

			finalgsdstddev = pow( (sum / ((double) nfgsd) ), 0.5 );
		}

		//open output file and create header
		log_finalgsd_fname << "SCORE." << simid << ".Rediscretized.FinalGSD.csv";
		log_finalgsd_file.open ( log_finalgsd_fname.str().c_str() );
		//write file header
		log_finalgsd_file << "Time=" << now << "s,NGrains=" << nfgsd << ",Vref=" << finalgsdsum << ",Mean=" << setprecision(8) << finalgsdmean << ",Median=" << setprecision(8) << finalgsdmedian << ",StdDev=" << setprecision(8) << finalgsdstddev << ",Min=" << finalgsdmin << ",Max=" << finalgsdmax << "\n";
		log_finalgsd_file << "Rank;Volume/Cells;IncubationTime/s;Ideal(RANDOMis0);ln(Vi/Vav);CumSumi/Vref*100%\n";

		double cumsum = 0.0;
		double volnorm = 0.0;
		double probability = 0.0;
		//output probability plot raw data and final grain size distribution
		for ( unsigned long n = 0; n < nfgsd; ++n ) { 
			cumsum = cumsum + finalgsd[n].finalvol;
			volnorm = finalgsd[n].finalvol / finalgsdmean;
			probability = cumsum / finalgsdsum * 100.0;

			log_finalgsd_file << finalgsd[n].rank << ";" << finalgsd[n].finalvol << ";" << finalgsd[n].tincub << ";" << finalgsd[n].ideal;

			if ( volnorm > SMALLEST_TIME_LOGARITMIZE )	log_finalgsd_file << ";" << (log(volnorm)) << ";" << probability << "\n";
			else										log_finalgsd_file << ";;\n";
		}

		log_finalgsd_file.flush();
		log_finalgsd_file.close();
	}


	delete [] nnuc_per_rank;		nnuc_per_rank = NULL;
	delete [] nnuc_per_rank_cumul;	nnuc_per_rank_cumul = NULL;
	delete [] grbuf_rank;			grbuf_rank = NULL;
	ensMemGuard = ensMemGuard - (ensemble_myCAs_nnuclei * sizeof(MPI_IO_FGSDComplete));


	return true;
	//vector:we can only safely use operator[]() (or at()) to modify elements that are actually present, which means that they count towards size()
	//www.gotw.ca gotw 074.htm
}



void ensembleHdl::destroy_myCAs( void )
{
	//#pragma omp parallel end all threads if not already happen
	
	//necessary because the ensembleHdl member vector<caHdlP> can only deconstruct the vector but not the heap memory the individual element pointers are pointing to
	//vector<vector<caHdlP>> refers to process local memory in the same shared memory but with different distance, the OS knows which snippets of heap process ID has spawn by issuing threads
	for ( uint32_t tid = 0; tid < myCAs.size(); tid++) {

		//cout << myRank << "\t\t" << tid << "\t\tmyRank;tid deallocating..." << endl;

		for ( uint32_t tidlocalca = 0; tidlocalca < myCAs[tid].size(); tidlocalca++) {
			delete myCAs[tid][tidlocalca]; //this is expected, but calls the constructor of the caHdl class that automatically destroys memory used for representing the automaton and storage of intermediate results
			myCAs[tid][tidlocalca] = NULL;
			ensMemGuard = ensMemGuard - sizeof(caHdl);
			//cout << myRank << "\t\t" << tid << "\t\t" << tidlocalca << "\t\tmyRank;tid;tidlocalca successful." << endl;
		}
	}
}




caHdl::caHdl()
{
	//myensHdl nothing to free is only a reference to an address
	myensHdl = NULL;	//upon construction no manager has been assigned to that object yet!, no jobid = worldca-id as well
	ccnumatid = 0;
	jobid = 0;

	myCAGeometry.nNucleiCSR = 1;
	myCAGeometry.nboxedge_rd = CA_DIMENSIONS_MINIMUM;
	myCAGeometry.nboxedge_nd = CA_DIMENSIONS_MINIMUM;
	myCAGeometry.nboxedge_td = CA_DIMENSIONS_MINIMUM;
	myCAGeometry.nboxarea_rdnd = SQR(CA_DIMENSIONS_MINIMUM);
	myCAGeometry.nboxvol_rdndtd = CUBE(CA_DIMENSIONS_MINIMUM);

	myCAGeometry.cellsize = DEFAULT_CELLSIZE;
	myCAGeometry.boxedge_rd = CA_DIMENSIONS_MINIMUM * DEFAULT_CELLSIZE;
	myCAGeometry.boxedge_nd = CA_DIMENSIONS_MINIMUM * DEFAULT_CELLSIZE;
	myCAGeometry.boxedge_td = CA_DIMENSIONS_MINIMUM * DEFAULT_CELLSIZE;
	myCAGeometry.boxarea_rdnd = SQR(CA_DIMENSIONS_MINIMUM)*SQR(DEFAULT_CELLSIZE);
	myCAGeometry.boxvol_rdndtd = CUBE(CA_DIMENSIONS_MINIMUM)*CUBE(DEFAULT_CELLSIZE);

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
	myPhysData.defgmean_nd = DEFAULT_DEFGSIZE;
	myPhysData.defgmean_td = DEFAULT_DEFGSIZE;
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
	myCADefMS.v_ynd = 0.0;
	myCADefMS.w_ztd = 0.0;

	myNucleationModel.gbnucleation = NO_GBNUCLEATION;
	myNucleationModel.csrnucleation = NO_CSRNUCLEATION;
	myNucleationModel.clustnucleation = NO_CLUSTNUCLEATION;

	tmpdefgseeds = NULL;
	ndefgseeds = 0;

	myFullSeedList = NULL;
	myRecycList = NULL;
	mySeedFront = NULL;
	ntotalSeedFront = 0;
	nextSlotNeverActiveSeedFront = 0;
	nCurrentlyActiveSeedFront = 0;
	ntotalFullSeedList = 0;
	nextSlotToFullSeed = 0;
	ntotalRecycList = 0;
	nextSlotThatBecomesRecyc = 0;
	firstNotRecycYet = 0;


	//mydefggrid = NULL;
	mycellgrid = NULL;
	myFullRXList = NULL;
	myRecyclingList = NULL;
	myRXFront = NULL;
	myrxgpoolboundarycontact = NULL;

	ntotalRXFront = 0;
	nextSlotNeverActiveRXFront = 0;
	nCurrentlyActive = 0;
	ntotalFullRXList = 0;
	nextSlotToFullRX = 0;
	ntotalRecyclingList = 0;
	nextSlotThatBecomesRecycled = 0;
	firstNotRecycledYet = 0;

	tprocessingend = 0.0;

	//##MK::detect grain boundaries fast
	boundaryBucketFast = NULL;
	boundaryBucketSizeFast = 0;
	gbCountFast = 0;

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
	mypzmin = 0.0;


	Gbhalfsq = 0.5 * myPhysData.G0 * SQR(myPhysData.bZeroCelsius);
	//_kT = (1.0 / (kboltzman * CurrentTemperature));
	mLAGB = myPhysData.LAGBm0 * exp( -1.0 *  DEFAULT_LAGB_HACT / (kboltzman * CurrentTemperature) );
	mHAGB = myPhysData.HAGBm0 * exp( -1.0 *  DEFAULT_HAGB_HACT / (kboltzman * CurrentTemperature) );
	mGS = myPhysData.GSm0 * exp( -1.0 *  DEFAULT_GS_HACT / (kboltzman * CurrentTemperature) );
	mRHHAGB = myPhysData.RH_HAGBm0 * exp( -1.0 * myPhysData.RH_HAGBHact / (kboltzman * CurrentTemperature) );

	nmyoripool = 0;
	nmydefgpool = 0;
	nmyrxgpool = 0;
	nmynuclei = 0;

	loginfo_rxfrontstats_cnt = 0;
	loginfo_grainevo_cnt = 0;
	loginfo_rendering_cnt = 0;
	loginfo_defrag_cnt = 0;

	myhpfwindow.dx = DEFAULT_HPF_WINDOWSIZE;
	myhpfwindow.dy = DEFAULT_HPF_WINDOWSIZE;
	myhpfwindow.dz = DEFAULT_HPF_WINDOWSIZE;
	myhpfwindow.nx = (2 * myhpfwindow.dx) + 1;
	myhpfwindow.ny = (2 * myhpfwindow.dy) + 1;
	myhpfwindow.nz = (2 * myhpfwindow.dz) + 1;
	myhpfwindow.nxy = myhpfwindow.nx * myhpfwindow.ny;
	myhpfwindow.nxyz = myhpfwindow.nxy * myhpfwindow.nz;
	myhpfvalues = NULL;

	myensRank = MASTER;
	nRanks = 1;

	outopt_localrenderhow = RENDERING_MSNO;
	outopt_localrendercolor = RENDERING_COLOR_GRAINID;
	outopt_localrenderboundaries = RENDERING_BOUNDARIES_NO;
	outopt_logboundaries = OUTPUT_LOGBND_NO;
	outopt_localrxfront = OUTPUT_RXFRONTSTATS_NO;
	outopt_localsinglegrainevo = OUTPUT_SINGLEGRAIN_NO;

	mySuccess = true;
	renderingForThisCA = false;
	outopt_localhpf3d2d = false;

	localprng.init ( DEFAULT_PRNG_SEED );
}


caHdl::~caHdl()
{
	//because caHdl does not know about dynamic allocation of memory from within the struct
	for ( uint32_t s = 0; s < mygrainevolution.size(); s++) {
		delete [] mygrainevolution[s].localdatabucket;
		mygrainevolution[s].localdatabucket = NULL;
		myMemGuard = myMemGuard - ( (double) mygrainevolution[s].nlocaldata * sizeof(uint32_t) );
	}
	//vector vector does allocate self

	//fragmentation radar
	for ( uint32_t s = 0; s < this->myrxfrontfragmentation.size(); s++ ) {
		delete [] myrxfrontfragmentation[s].rgba;
		myrxfrontfragmentation[s].rgba = NULL;
		myMemGuard = myMemGuard - ( (double) myrxfrontfragmentation[s].nblocks * 4 * sizeof(unsigned char) );
	}

	//tmpdefgseeds already deleted after successful execution of pickedgrains
	//##DEBUG::delete [] mydefggrid;
	//##DEBUG::mydefggrid = NULL;
	//##DEBUG::myMemGuard = myMemGuard - (myCADefMS.ngrxyz * sizeof(uint32_t));
	/*delete [] mycellgrid;
	mycellgrid = NULL;
	myMemGuard = myMemGuard - (myCAGeometry.nboxvol_rdndtd * sizeof(uint32_t));
	delete [] myFullRXList;
	myFullRXList = NULL;
	myMemGuard = myMemGuard - (ntotalFullRXList * sizeof(uint32_t));
	delete [] myRecyclingList;
	myRecyclingList = NULL;
	myMemGuard = myMemGuard - (ntotalRecyclingList * sizeof(uint32_t));
	delete [] myRXFront;
	myRXFront = NULL;
	myMemGuard = myMemGuard - (ntotalRXFront * sizeof(cell));*/
	delete [] myrxgpoolboundarycontact;
	myrxgpoolboundarycontact = NULL;
	myMemGuard = myMemGuard - (myrxgpool.size() * sizeof(bool));

	delete [] myhpfvalues;
	myhpfvalues = NULL;

	delete [] DefMicrotexture;
	DefMicrotexture = NULL;

	delete [] NucleiMicrotexture;
	NucleiMicrotexture = NULL;
}


void caHdl::cleanupMemoryUsedForGrowthMachine( void )
{
	//MK::when the solve_RXGROWTH has been executed the mygrainevolution contains the intermediate stage of all grains but 
	//the memory for the cellgrid and the dynamic memory allocated for managing cells in no longer necessary so delete it

	delete [] mycellgrid;
	mycellgrid = NULL;
	myMemGuard = myMemGuard - ( myCAGeometry.nboxvol_rdndtd * sizeof(uint32_t));
	delete [] myFullRXList;
	myFullRXList = NULL;
	myMemGuard = myMemGuard - (ntotalFullRXList * sizeof(uint32_t));
	delete [] myRecyclingList;
	myRecyclingList = NULL;
	myMemGuard = myMemGuard - (ntotalRecyclingList * sizeof(uint32_t));
	delete [] myRXFront;
	myRXFront = NULL;
	myMemGuard = myMemGuard - (ntotalRXFront * sizeof(cell));

	//reduces memory consumption to generate the myMaternUniverse
//cout << "Capacity/Size of my point process is = " << mypointprocess.capacity() << ";" << mypointprocess.size() << endl;
	vector<point>().swap( mypointprocess );
//cout << "Capacity/Size of my point process is = " << mypointprocess.capacity() << ";" << mypointprocess.size() << endl;

}


uint32_t caHdl::ca_get_closest_standardlage( double  * quat )
{
	double closest_disori = MAX_FCC;
	uint32_t closest_ideal = CATEGORIZED_AS_RANDOM;

	double q10 = quat[0];
	double q11 = quat[1];
	double q12 = quat[2];
	double q13 = quat[3];

	ensembleHdlP ens = this->myensHdl;

	for ( uint32_t cand = 0; cand < ens->standardlagen.size(); cand++) {
		double disori = misorientationCubicQxQ ( q10, q11, q12, q13,  ens->standardlagen[cand].q0, ens->standardlagen[cand].q1, ens->standardlagen[cand].q2, ens->standardlagen[cand].q3 );

		//scatter within the range and closer as to all other candidates?
		if ( disori <= ens->standardlagen[cand].scatter && disori <= closest_disori ) {
			closest_disori = disori;
			closest_ideal = cand;
		}
		//for all possible standardlagen candidates find always the closest match even when scatter range of two components overlaps such that macrotexture properly adds to 1.0
	}

	//MK::within the code orientations categorized as RANDOM are at the zero-th position
	if ( closest_ideal == CATEGORIZED_AS_RANDOM )
		return 0;
	//else
		return closest_ideal + 1;
}


uint32_t caHdl::ca_check_disjunctness( double * bunge )
{
	QUICKASSERT ( myoripool.size() < MAXIMUM_DISJOINT_ORIS );

	uint32_t closestid = UNKNOWN_ORIENTATION;
	double disori = 2 * _PI_;
	double closestdisori = RESOLUTION_SO3GRID;
	double qbunge[4], qcand[4];

	euler2quaternion( bunge, qbunge );

	//check disorientation to all other components in myoripool
	for (uint32_t cand = 0; cand < this->myoripool.size(); cand++ ) {
		qcand[0] = myoripool[cand].q0;
		qcand[1] = myoripool[cand].q1;
		qcand[2] = myoripool[cand].q2;
		qcand[3] = myoripool[cand].q3;
		disori = misorientationCubicQxQ( qbunge[0], qbunge[1], qbunge[2], qbunge[3], qcand[0], qcand[1], qcand[2], qcand[3]);

		//closer than the SO3GRID resolution to any other orientation already known?
		if (disori <= closestdisori) {
			closestdisori = disori;
			closestid = cand;
		}
		//however still testing against all other to get the closest match
	}

	//most likely case
	if ( closestid == UNKNOWN_ORIENTATION ) {
		struct ori anori;
		anori.bunge1 = bunge[0];
		anori.bunge2 = bunge[1];
		anori.bunge3 = bunge[2];

		anori.q0 = qbunge[0];
		anori.q1 = qbunge[1];
		anori.q2 = qbunge[2];
		anori.q3 = qbunge[3];

		anori.closestideal = ca_get_closest_standardlage( qbunge ); //RANDOM or one of our components

		anori.RGBA[RED] = UCHAR_RANGE_MIN;
		anori.RGBA[GREEN] = UCHAR_RANGE_MIN;
		anori.RGBA[BLUE] = UCHAR_RANGE_MIN;
		anori.RGBA[ALPHA] = UCHAR_RANGE_MAX;


		myoripool.push_back( anori );

//cout << "ADD\t\t" << (myoripool.size() - 1) << "\t\t" << anori.bunge1 << "\t\t" << anori.bunge2 << "\t\t" << anori.bunge3;
//cout << "\t\t" << anori.q0 << "\t\t" << anori.q1 << "\t\t" << anori.q2 << "\t\t" << anori.q3 << endl;
//#ifdef DETAILED_PROMPTS
//if ( anori.closestideal == RANDOM_ORIENTATION ) { cout << "Node;" << this->myRank << ", I have categorized close to RANDOM." << endl; }
//else { cout << "Node;" << this->myRank << ", I have categorized close to listentry " << (anori.closestideal + 1) << endl; }
//#endif DETAILED_PROMPTS

		//make known the new orientation
		return (myoripool.size() - 1);
	}

	//else present closest already existent orientation id
	//cout << "GET\t\t" << closestid << "\t\t" << closestdisori << endl;

	//one of the goodfellas...
	return closestid;
}


void caHdl::nes_networkgrowthmodel_vacancycorediff( void )
{
	//MK::template code for a recovery according to Nes vacancy diffusion controlled themrally-activated glide
	//MK::Euler forward difference scheme with dr/dt ~ r(t+dt) - r(t) / dt
	//model parameter at this->t
	double Bx;
	double B = _PI_ * myRecoveryModel.NesC3 * myRecoveryModel.NesC5 * myRecoveryModel.VacancyDiffD0 * exp ( -1.0 * myRecoveryModel.VacancyDiffHact / ( kboltzman * CurrentTemperature ));
	B /= myRecoveryModel.NesKappa2;
	double exponent = SQR(CurrentBurgersVector) / (kboltzman * CurrentTemperature);
	double expterm = exp( myRecoveryModel.NesAlpha3 * myRecoveryModel.NesKappa2 * CurrentG * CurrentBurgersVector * exponent );

	uint32_t ndg = mydefgpool.size();
	double r_tdt, r_t;
	double currdt = this->dt;
	double currrhomax = RHOMAX_WELLANNEALED;

	for ( uint32_t dg = 0; dg < ndg; dg++ ) {
		r_t = 1.0 / pow( mydefgpool[dg].rho, 0.5 ); //>1.0

		Bx = (B / r_t) * expterm;

		r_tdt = (Bx * currdt) + r_t;

		mydefgpool[dg].rho = 1.0 / SQR(r_tdt);

		if ( mydefgpool[dg].rho > currrhomax ) 
			currrhomax = mydefgpool[dg].rho;
	}

	this->myrhomax = currrhomax;
}


void caHdl::nes_networkgrowthmodel_solutedrag( void )
{
	//#MK::not implemented yet
}


void caHdl::nes_michalakpaxton( void ) 
{
//set of parameter determined as in Michalak Paxton 1961, TransAIME, 221, 1961 p853
//and Nes, Recovery revisited
	double C5 = 2.0/3.0;
	double kb = 8.61733e-5;
	double k3a3 = 0.1;
	double Usd = 2.90;
	//double atau2 = 6.92878e-13;
	//double btau5 = 4.72858e-14;
	double f1 = 0.25;
	double f2 = 0.75;

	//double temp = CurrentG * CurrentBurgersVector;
	//double tau2 = pow( (atau2 * exp( ((C5 * Usd) - (k3a3 * temp * SQR(CurrentBurgersVector))) / (kb*CurrentTemperature) )), -1.0);
	//double tau5 = pow( (btau5 * exp( (C5 * Usd)/(kb*CurrentTemperature) )), -1.0);
	double tau2 = 0.3605; //min
	double tau5 = 5.2824; //min

	uint32_t ndg = mydefgpool.size();
	double currrhomax = RHOMAX_WELLANNEALED;
	double currtime_minutes = this->t / 60.0;

	double factor = f1 * pow( (1 + ( currtime_minutes / tau2 )), -1.0 ) + f2 * pow( (1 + ( currtime_minutes / tau5 )), -0.5 );

	//"recovery modeling" in each grain
	for ( uint32_t dg = 0; dg < ndg; dg++ ) {
		mydefgpool[dg].rho = mydefgpool[dg].rho0 * SQR(factor);

		if ( mydefgpool[dg].rho > currrhomax ) 
			currrhomax = mydefgpool[dg].rho;
	}

	this->myrhomax = currrhomax;

cout << "Michalak Paxton Recovery;tau2;tau5;factor;currrhomax\t\t" << this->t << ";" << this->X << ";" << tau2 << ";" << tau5 << ";" << factor << ";" << currrhomax << endl;
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

//cout << this->time[0] << ";" << this->temperature[0] << endl;
//cout << this->time[1] << ";" << this->temperature[1] << endl;
//cout << "i=" << i << "\ttime[i]=" << time[i] << "\tt=" << t << endl;
	//scan vector t and perform linear interpolation to get Ti(ti)

	while (ttime[i] < t && i < ttime.size() ) {
		i++;
	}

//cout << "afterwhile i is=" << i << endl;

	if (i == 0) {
//cout << "temperature[0]=" << temperature[0] << endl;
		CurrentTemperature = ttemperature[0];
		return;
	}
	if (i >= ttime.size() ) {
//cout << "temperature[ntimetempstamps-1]=" << temperature[ntimetempstamps-1] << endl;
		CurrentTemperature = ttemperature[ttime.size() - 1];
		return;
	}

	CurrentTemperature = ttemperature[i-1] + ( (ttemperature[i] - ttemperature[i-1]) / (ttime[i] - ttime[i-1]) ) * (t - ttime[i-1]);
//cout << "currT=" << currT << endl;
}


void caHdl::update_atomisticproperties( void )
{
	//temperature dependent shear modulus if desired
	//see Nadal, M-H and Le Poac, P in J. of Appl. Phys. 93 2472 2003 for further details what to do close to the melting point
	//MK::for Alu, also applied to TWIP steel
	//CurrentG = myPhysData.G0 - (myPhysData.dGdt * ( CurrentTemperature / myPhysData.Tmelt ));

	//overwritten for iron in accordance with Ledstetter, Reed, 1973
	//MK::for bcc-iron
	CurrentG = myPhysData.G0 * (1.0 - ((9e-10 * CUBE(CurrentTemperature) - 7e-7*SQR(CurrentTemperature) + 0.0003*CurrentTemperature - 0.00028)));

	//temperature dependent Burgers vector as a result of lattice elongation/contraction
	double T = CurrentTemperature - TOFFSET;
	//MK::for Alu, also applied to TWIP steel
	//CurrentBurgersVector = myPhysData.bZeroCelsius * ( 1.0 + myPhysData.thermexp_C * ((myPhysData.thermexp_a * T) + ( myPhysData.thermexp_b * SQR(T) )) * 1e-6 );

	//MK::for bcc iron up to 800deg celsius
	CurrentBurgersVector = myPhysData.bZeroCelsius;
	//##add here further any temperature dependent values which need update, e.g. recovery parameters

	//cache Gbhalfsq
	Gbhalfsq = 0.5 * CurrentG * SQR( CurrentBurgersVector );

//cout << "time/G/b = " << this->t << ";" << this->X << ";" << CurrentG << ";" << CurrentBurgersVector << endl;
}


void caHdl::update_intrinsicmobilities( void )
{
	//intrinsic grain boundary mobilities
	mLAGB = myPhysData.LAGBm0 * exp( -1.0 *  myPhysData.LAGBHact / (kboltzman * CurrentTemperature) );
	mHAGB = myPhysData.HAGBm0 * exp( -1.0 *  myPhysData.HAGBHact / (kboltzman * CurrentTemperature) );
	mGS = myPhysData.GSm0 * exp( -1.0 *  myPhysData.GSHact / (kboltzman * CurrentTemperature) );

	mRHHAGB = myPhysData.RH_HAGBm0 * exp( -1.0 * myPhysData.RH_HAGBHact / (kboltzman * CurrentTemperature) );

	QUICKASSERT ( mLAGB > 0.0 );
	QUICKASSERT ( mHAGB > mLAGB );
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

	while ( zenertime[i] < this->t && i < zenertime.size() ) {
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
	//##MK::-->inject code for update of solute drag or any other microchemistry model here
}


void caHdl::update_deformedsubstructure( void )
{
	//##MK::inject the recovery model here which updates rho, dgav, dav
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

	//evolution of the average subgrain size affecting mydefgpool[i].dav


	//evolution of the average misorientation in the deformed grain affecting mydefgpool[i].dgav
}


inline double caHdl::get_currentintrinsicmobility( double Pvalue )
{
	if ( mobilitymodel == MOBILITYMODEL_ROLLETTHOLM )
		return ( Pvalue * mRHHAGB );


	//MOBILITYMODEL_SEBALDGOTTSTEIN
	if (Pvalue < 0.0) {
		return mLAGB;
	}
	return ( ( Pvalue * mGS ) + ( (1.0 - Pvalue) * mHAGB ) );
	//most likely case Pvalue >= 0.0 but then two comparisons...
}


double caHdl::calc_mobilityweight ( uint32_t rgpoolid , uint32_t dgpoolid )
{
	//Sebald Gottstein model
	double q1[4], q2[4], qdis[4];

	uint32_t rgori = myrxgpool[rgpoolid].caori;
	uint32_t dgori = mydefgpool[dgpoolid].caori;

	//now calculate on the fly...
	q1[0] = myoripool[rgori].q0;
	q1[1] = myoripool[rgori].q1;
	q1[2] = myoripool[rgori].q2;
	q1[3] = myoripool[rgori].q3;

	q2[0] = myoripool[dgori].q0;
	q2[1] = myoripool[dgori].q1;
	q2[2] = myoripool[dgori].q2;
	q2[3] = myoripool[dgori].q3;

	//P=weightedMob for all HAGB already defined
	double weightedMob = 0.0;

	misorientationQuaternionCubic( q1, q2, qdis );
	double theta = qdis[0];
	if( theta > 1.0 ) {
		theta = (double) (int) theta;
	}
	theta = 2*acos(theta);


	if ( this->mobilitymodel == MOBILITYMODEL_ROLLETTHOLM ) {
		//P value allows to calculate the mobility as often occurring in works by Rollett and Holm via P * mRHHAGB
		weightedMob = 1.0 - (myPhysData.RH_LAGBHAGBcut * exp( -1.0 * myPhysData.RH_LAGBHAGBtrans * pow( (theta/MAXDISORI_LAGB2HAGB), myPhysData.RH_LAGBHAGBexponent ) ));

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
	//the quaternion that describes a 40°<111> misorientation
	double m40_111[4] = { cos( oneNinth * _PI_ ), _sqrt3 * sin( oneNinth * _PI_ ), _sqrt3 * sin( oneNinth * _PI_ ), _sqrt3 * sin( oneNinth * _PI_ ) };

	//obviously a HAGB but also eligible for GS? check proximity to 40deg111
	double dev_40_111 = misorientationCubicQxQ( qdis[0], qdis[1], qdis[2], qdis[3], m40_111[0], m40_111[1], m40_111[2], m40_111[3] );

	if( dev_40_111 < maxDev40_111 ) {
		weightedMob = SQR( cos( 0.5 * _PI_ * dev_40_111 / maxDev40_111 ) );
	}

	return weightedMob;
}



double caHdl::calc_mobilityweight_sebaldgottstein ( uint32_t rgpoolid , uint32_t dgpoolid )
{
	//Sebald Gottstein model
	double q1[4], q2[4], qdis[4];

	double weightedMob = 0.0;
	double maxDev40_111 = MAXDISORI_TO_40DEG111;
	double _sqrt3 = 1 / sqrt( 3.0 );
	double oneNinth = 1.0 / 9.0; 
	//the quaternion that describes a 40°<111> misorientation
	double m40_111[4] = { cos( oneNinth * _PI_ ), _sqrt3 * sin( oneNinth * _PI_ ), _sqrt3 * sin( oneNinth * _PI_ ), _sqrt3 * sin( oneNinth * _PI_ ) };

	uint32_t rgori = myrxgpool[rgpoolid].caori;
	uint32_t dgori = mydefgpool[dgpoolid].caori;

	//now calculate on the fly...
	q1[0] = myoripool[rgori].q0;
	q1[1] = myoripool[rgori].q1;
	q1[2] = myoripool[rgori].q2;
	q1[3] = myoripool[rgori].q3;

	q2[0] = myoripool[dgori].q0;
	q2[1] = myoripool[dgori].q1;
	q2[2] = myoripool[dgori].q2;
	q2[3] = myoripool[dgori].q3;

	//P=weightedMob for all HAGB already defined
	weightedMob = 0.0;

	misorientationQuaternionCubic( q1, q2, qdis );
	double theta = qdis[0];
	if( theta > 1.0 ) {
		theta = (double) (int) theta;
	}
	theta = 2*acos(theta);

	//assign categorical intrinsic boundary mobility to disorientation among two orientation
	//LAGB detected, overwrite default value
	if( theta <= MAXDISORI_LAGB2HAGB ) {
		weightedMob = -1.0;

		return weightedMob;
	}

	//obviously a HAGB but also eligible for GS? check proximity to 40deg111
	double dev_40_111 = misorientationCubicQxQ( qdis[0], qdis[1], qdis[2], qdis[3], m40_111[0], m40_111[1], m40_111[2], m40_111[3] );

	if( dev_40_111 < maxDev40_111 ) {
		weightedMob = SQR( cos( 0.5 * _PI_ * dev_40_111 / maxDev40_111 ) );
	}

	return weightedMob;
}


double caHdl::calc_mobilityweight_rollettholm ( uint32_t rgp, uint32_t dgp )
{
	double q1[4], q2[4], qdis[4];

	uint32_t rgori = myrxgpool[rgp].caori;
	uint32_t dgori = mydefgpool[dgp].caori;

	//now calculate on the fly...
	q1[0] = myoripool[rgori].q0;
	q1[1] = myoripool[rgori].q1;
	q1[2] = myoripool[rgori].q2;
	q1[3] = myoripool[rgori].q3;

	q2[0] = myoripool[dgori].q0;
	q2[1] = myoripool[dgori].q1;
	q2[2] = myoripool[dgori].q2;
	q2[3] = myoripool[dgori].q3;

	//disorientation angle
	misorientationQuaternionCubic( q1, q2, qdis );
	double theta = qdis[0];
	if( theta > 1.0 ) {
		theta = (double) (int) theta;
	}
	theta = 2*acos(theta);

	//P value allows to calculate the mobility as often occurring in works by Rollett and Holm via P * mRHHAGB
	double P = 1.0 - (myPhysData.RH_LAGBHAGBcut * exp( -1.0 * myPhysData.RH_LAGBHAGBtrans * pow( (theta/MAXDISORI_LAGB2HAGB), myPhysData.RH_LAGBHAGBexponent ) ));

	return P;
}


inline double caHdl::get_zener( void )
{
	if ( myDragData.ZenerConsider == DISPERSOIDDRAG_NO ) 
		return 0.0;

	//DISPERSOIDDRAG_CONSTANT || DISPERSOIDDRAG_TIMEDEP
	//##MK::at the moment constant drag in each grain
	double pz =  1.0 * myDragData.ZenerAlpha * myDragData.ZenerGamma * myDragData.fr;

	return pz;
}


inline double caHdl::get_rho ( uint32_t mydefgpoolid )
{
	//MK::interpret mydefgpool rho
	return mydefgpool[mydefgpoolid].rho;
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

	//but there might not at all be any nucleus-defg achieving the this high mobility to completely

	//when there are only LAGB boundaries in the system
	double m_max = ( myMobilityWeightMax * mgs ) + ( (1.0 - myMobilityWeightMax) * mhagb );

	//only LAGB in the system
	if ( myMobilityWeightMax < (0.0 - DEFAULT_SMALL_NUMBER) ) {
		m_max = myPhysData.LAGBm0 * exp( -1.0 * myPhysData.LAGBHact / (kboltzman * CurrentTemperature) );
	}

	//maximum migration velocity
	double v_max = m_max * ((Gbhalfsq * myrhomax) - mypzmin);

	//transform only a fraction of a cell in a time step
	dtmax = ((maxfillperstep * myCAGeometry.cellsize) / v_max);

	#ifdef DETAILED_PROMPTS
		cout << this->jobid << ";jobid;" << m_max << ";m_max;" << v_max << ";vmax;" << tmin << ";tmin" << endl;
	#endif

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

	//###disjunctness of time-temperature profile has to be assured by user
	if (i < ttime.size() ) {
		double dTdt = ( fabs(ttemperature[i] - ttemperature[i-1]) / (ttime[i] - ttime[i-1]) ); 

		if (dTdt <= SMALL_HEATRATE) dTdt = SMALL_HEATRATE; //avoid division by zero error

		dtmax = (SMALL_HEAT / dTdt);
		return dtmax;
	}

	//if not already returned, then processing schedule ends, temperature stays constant with last value, nevertheless so dT/dt = 0
	dtmax = INFINITE;
	return dtmax;
}


inline double caHdl::get_dtmax_instantslope_rho( void )
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


inline double caHdl::get_dtmax_instantslope_zener( void )
{
	if ( myDragData.ZenerConsider != DISPERSOIDDRAG_TIMEDEP ) { //when there is no ZenerDrag simulated or the drag is not time-dependent there is no necessity to adjust the numerical time integration of the CA
		return (INFINITE);
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

		if (dpzdt <= SMALL_DRAGGINGRATE) dpzdt = SMALL_DRAGGINGRATE; //avoid division by zero error

		dtmax = (SMALL_ZENERFORCE / dpzdt);
		return dtmax;
	}

	//if not already returned, then processing schedule has ended
	dtmax = INFINITE;

	return dtmax;
}


inline double caHdl::get_dtmax_instantslope_nucleation( void )
{
	
	if ( this->myNucleationModel.tincubmodel != TINCUB_TIMEDEPENDENT )
		return INFINITE;

	//nucleation is time dependent
	double dtmax = INFINITE;

	double dNdt = (double) myrxgpool.size() * pdf_rayleigh( t ); //##MK::all planned nuclei, MK::not t+dt because function is being called after this->t was incremented!

	if ( dNdt <= SMALL_NUCLEATIONRATE ) dNdt = SMALL_NUCLEATIONRATE; //avoid division by zero error

	dtmax = (SMALL_NUMBEROFNUCLEI / dNdt);

	return dtmax;
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

	dtmax_model = get_dtmax_instantslope_nucleation();
	if ( dtmax_model < min_dtmax )
		min_dtmax = dtmax_model;

	dtmax_model = get_dtmax_instantslope_cells();
	if ( dtmax_model < min_dtmax )
		min_dtmax = dtmax_model;


	return min_dtmax;
}



void caHdl::init_cahdlprng( void )
{
	//MK::jobid can be 0
	long localseed = (long) -1*( this->jobid  + 1); //##MK::myensHdl->myRank would render as many same jobs locally as there are jobs per process
	localseed = localseed - 10;

	this->localprng.init( localseed );


#ifdef DETAILED_PROMPTS
	cout << "myRank " << this->myensHdl->myRank << "; ThreadID " << this->ccnumatid << "; Localseed " << localseed << endl;
#endif
}


void caHdl::init_parameter( void )
{
	//beneficial, because I/O only for each process, not each automaton...
	ensembleHdlP myens = this->myensHdl;

	XMAX = myens->XMAX;
	TMAX = myens->TMAX;
	NTSTEPSMAX = myens->NTSTEPSMAX;

	//##MK::CURRENTLY ALL CAs are alike
	myCAGeometry.cellsize = myens->ensCAGeometry.cellsize;
	this->_cellsize = (1.0 / myCAGeometry.cellsize);
	myCAGeometry.nNucleiCSR = myens->ensCAGeometry.nNucleiCSR;
	myCAGeometry.nboxedge_rd = myens->ensCAGeometry.nboxedge_rd;
	myCAGeometry.nboxedge_nd = myens->ensCAGeometry.nboxedge_nd;
	myCAGeometry.nboxedge_td = myens->ensCAGeometry.nboxedge_td;
	myCAGeometry.nboxarea_rdnd = myens->ensCAGeometry.nboxarea_rdnd;
	myCAGeometry.nboxvol_rdndtd = myens->ensCAGeometry.nboxvol_rdndtd;
	myCAGeometry.boxedge_rd = myens->ensCAGeometry.boxedge_rd;
	myCAGeometry.boxedge_nd = myens->ensCAGeometry.boxedge_nd;
	myCAGeometry.boxedge_td = myens->ensCAGeometry.boxedge_td;
	myCAGeometry.boxarea_rdnd = myens->ensCAGeometry.boxarea_rdnd;
	myCAGeometry.boxvol_rdndtd = myens->ensCAGeometry.boxvol_rdndtd;

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
	myPhysData.defgmean_nd = myens->ensPhysData.defgmean_nd;
	myPhysData.defgmean_td = myens->ensPhysData.defgmean_td;
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

	mobilitymodel = myens->mobilitymodel;
	nucleationmodel_status = NUCSITES_STILL_SOME_FREE;
	outopt_localrenderhow = myens->outopt_rendermethod;
	outopt_localrendercolor = myens->outopt_rendercolormodel;
	outopt_localrenderboundaries = myens->outopt_localrenderboundaries;
	outopt_logboundaries = myens->outopt_logboundaries;
	outopt_localrxfront = myens->outopt_rxfront;
	outopt_localsinglegrainevo = myens->outopt_singlegrainevo;
	outopt_localhpf3d2d = myens->outopt_hpf3d2donmaster;

	maxfillperstep = myens->maxfillperstep;
	initialRelCellCaching = myens->initialRelCellCaching;
	transientRelCellRecaching = myens->transientRelCellRecaching;

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

	//ideal orientations only to myens

	uint32_t nrenderingmax = myens->UserDefLogPoint_MS_Rendering.size();
	for ( uint32_t r = 0; r < nrenderingmax; r++ ) {
		rendering_atthisX.push_back ( myens->UserDefLogPoint_MS_Rendering[r] );
	}

	uint32_t ndefragmax = myens->UserDefLogPoint_X_CellListDefragment.size();
	for ( uint32_t d = 0; d < ndefragmax; d++) {
		defragRXFront_atthisX.push_back ( myens->UserDefLogPoint_X_CellListDefragment[d] );
	}

	uint32_t nlogmax = myens->UserDefLogPoint_X_Output.size();
	for ( uint32_t s = 0; s < nlogmax; s++ ) {
		output_atthisX.push_back ( myens->UserDefLogPoint_X_Output[s] );
	}

	//characterization of hostgrain switching probability function
	myhpfwindow.dx = DEFAULT_HPF_WINDOWSIZE;
	myhpfwindow.dy = DEFAULT_HPF_WINDOWSIZE;
	myhpfwindow.dz = DEFAULT_HPF_WINDOWSIZE;
	myhpfwindow.nx = (2 * myhpfwindow.dx) + 1;
	myhpfwindow.ny = (2 * myhpfwindow.dy) + 1;
	myhpfwindow.nz = (2 * myhpfwindow.dz) + 1;
	myhpfwindow.nxy = myhpfwindow.nx * myhpfwindow.ny;
	myhpfwindow.nxyz = myhpfwindow.nxy * myhpfwindow.nz;
	myhpfvalues = NULL;

	/*
	#ifdef REPORTSTYLE_DEVELOPER
		cout << jobid << "\t\t" << myensRank << "\t\t" << this->ccnumatid << "jobid;myensRank;ccnumathreadid; parameter have been loaded copied into caHdl" << endl;
	#endif
	*/
}


void caHdl::init_processing( void )
{
	ensembleHdlP myens = this->myensHdl;
	uint32_t ntmax = myens->ttime.size();
	uint32_t nTmax = myens->ttemperature.size();
	QUICKASSERT ( ntmax == nTmax );

	tprocessingend = 0.0;

	for ( uint32_t t = 0; t < ntmax; t++) {
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
	uint32_t ntmax = myens->dispersoidtime.size();
	uint32_t nZdmax = myens->dispersoidfr.size();
	QUICKASSERT ( ntmax == nZdmax );

	for ( uint32_t t = 0; t < ntmax; t++) {
		this->zenertime.push_back ( myens->dispersoidtime[t] );
		this->zenerdispersion.push_back ( myens->dispersoidfr[t] );
	}
}


void caHdl::solve_INITIALIZATION( void )
{
	init_cahdlprng();
	init_parameter();
	init_processing();
	init_zenerdrag();
}


defgseedP caHdl::init_discreteWindowFromMasterCSR( uint32_t * howmanyseeds )
{
	//MK::naively confining a small population of points in a confined window assuming a density of lambda yields different characteristics than sampling at random from
	//a large realization of the point process with intensity lambda, therefore: sample in tmpdefgseeds a subpopulation each point of which does not overlap with another point
	//cut out an arbitrary window the center point is not necessarily attached to a point of the process
	defgseedP seedlist = NULL;

	double AvGrainDiameterCells = (myPhysData.defgmean_poisson / myCAGeometry.cellsize);
	double tmp = myCAGeometry.nboxvol_rdndtd / ( 4.0/3.0 * _PI_ * pow( (AvGrainDiameterCells / 2.0), 3.0) );
	QUICKASSERT ( 0.5*tmp < UINT32T_MAX );
	uint32_t NumberOfPlannedSeeds = tmp;

//cout << AvGrainDiameterCells << ";" << tmp << ";" << NumberOfPlannedSeeds << endl;

	double MaximumSeeds = (1.0/CUBE(10.0)) * (double) myCAGeometry.nboxvol_rdndtd;
	QUICKASSERT ( NumberOfPlannedSeeds < MaximumSeeds ); //heuristic approach to be sure disjoint places can be found and discretization sufficient

	//MK::use scaling of Poisson point processes, which when realized in a unit cube the Voronoi cell has a volume of 4/3PI(defgmean_poisson/2)^3 then 1e6 points occupy
	double volmastercsr = 4.0/3.0 * _PI_ * pow ( (myPhysData.defgmean_poisson / 2.0), 3.0 );
	volmastercsr = volmastercsr * MASTER_CSR_NPOINTS;
	double edgelengthmastercsr = pow( volmastercsr, (1.0/3.0) );

	//allocate container, possibly needs reallocation
	QUICKASSERT ( 2*NumberOfPlannedSeeds < CA_ALLOCATION_MAXIMUM ); //##MK::memory reallocation
	uint32_t nseedlist = 2*NumberOfPlannedSeeds;
	seedlist = new struct defgseed[nseedlist];
	QUICKASSERT ( seedlist != NULL );
	for ( uint32_t s = 0; s < nseedlist; ++s ) {
		seedlist[s].ensdefgpoolid = NOT_ASSIGNED_YET;
		seedlist[s].mydefgpoolid = I_DONT_KNOW_YET;
		seedlist[s].location = NOT_ASSIGNED_YET;
	}

	//populate a CSR point pattern on a unit cube
	pointP pp3 = NULL;
	pp3 = new point[MASTER_CSR_NPOINTS];
	QUICKASSERT ( pp3 != NULL );

	for ( uint32_t p = 0; p < MASTER_CSR_NPOINTS; ++p ) {
		pp3[p].x = localprng.leEcuyer(); //[0,1]^3
		pp3[p].y = localprng.leEcuyer();
		pp3[p].z = localprng.leEcuyer();
	}

	//find owin enclosed in unit cube
	double owin[3] = {0.0, 0.0, 0.0};
	owin[0] = (myCAGeometry.boxedge_rd / edgelengthmastercsr); if ( owin[0] > 1.0 ) { cout << "ERROR::MASTER point process is too small!" << endl; return NULL; }
	owin[1] = (myCAGeometry.boxedge_nd / edgelengthmastercsr); if ( owin[1] > 1.0 ) { cout << "ERROR::MASTER point process is too small!" << endl; return NULL; }
	owin[2] = (myCAGeometry.boxedge_td / edgelengthmastercsr); if ( owin[2] > 1.0 ) { cout << "ERROR::MASTER point process is too small!" << endl; return NULL; }

	double owinoo[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	bool inside = false;

	while (inside == false ) {
		owinoo[XMI] = localprng.leEcuyer();
		owinoo[XMX] = owinoo[XMI] + owin[0];
		owinoo[YMI] = localprng.leEcuyer();
		owinoo[YMX] = owinoo[YMI] + owin[1];
		owinoo[ZMI] = localprng.leEcuyer();
		owinoo[ZMX] = owinoo[ZMI] + owin[2];

		if ( (owinoo[XMX] < 1.0) && (owinoo[YMX] < 1.0) && (owinoo[ZMX] < 1.0) ) //all coordinates positive owinoo cannot protrude in negative direction
			inside = true;
	}

//cout << "owin[012] = " << owin[0] << ";" << owin[1] << ";" << owin[2] << "\t\t" << owinoo[XMI] << ";" << owinoo[XMX] << ";" << owinoo[YMI] << ";" << owinoo[YMX] << ";" << owinoo[ZMI] << ";" << owinoo[ZMX] << endl;

	//extract all points from this window and sample via (i - owinoo[iMI])/owin[i] * ni into an automaton coordinate
	double dx, dy, dz;
	uint32_t ix, iy, iz;
	uint32_t hole;
	uint32_t nbx = myCAGeometry.nboxedge_rd; //avoid the cache collision when calling this in the for loop to select the points
	uint32_t nby = myCAGeometry.nboxedge_nd;
	uint32_t nbz = myCAGeometry.nboxedge_td;
	bool valid = false;
	uint32_t nvalidp = 0;

	for (unsigned int p = 0; p < MASTER_CSR_NPOINTS; ++p ) {
		if ( pp3[p].x < owinoo[XMI] ) continue;
		if ( pp3[p].x > owinoo[XMX] ) continue;
		if ( pp3[p].y < owinoo[YMI] ) continue;
		if ( pp3[p].y > owinoo[YMX] ) continue;
		if ( pp3[p].z < owinoo[ZMI] ) continue;
		if ( pp3[p].z > owinoo[ZMX] ) continue;

		//p is inside owinoo so, so dx, dy, dz > 0.0 and < 1.0, discretize location
		dx = (pp3[p].x - owinoo[XMI]) / (owinoo[XMX] - owinoo[XMI]);
		dy = (pp3[p].y - owinoo[YMI]) / (owinoo[YMX] - owinoo[YMI]);
		dz = (pp3[p].z - owinoo[ZMI]) / (owinoo[ZMX] - owinoo[ZMI]);
		ix = dx * nbx;
		iy = dy * nby;
		iz = dz * nbz;

//##DEBUGcout << ix << "\t\t" << iy << "\t\t" << iz << endl;

		hole = ix + (iy * nbx) + (iz * nbx * nby);

		//disprove the assumption that the hole contains no seed, overlap of two seeds within <= gridresu not allowed
		valid = true;
		for ( uint32_t s = 0; s < (nvalidp+1); s++ ) {
			if ( hole == seedlist[s].location ) { valid = false; break; }
		}

		//still true, okay hole is empty, plant a seed
		if ( valid == true ) {

			seedlist[nvalidp].location = hole;

//##DEBUGcout << "seedlist adding ... " << nvalidp << ";" << hole << endl;

			nvalidp++;

			//remalloc of seedlist?
			if ( nvalidp >= nseedlist ) {
				defgseedP tmp = NULL;
				tmp = new struct defgseed[2*nseedlist]; //that such be fine
				QUICKASSERT( tmp != NULL );
				for ( uint32_t used = 0; used < nseedlist; used++ ) {
					tmp[used].ensdefgpoolid = seedlist[used].ensdefgpoolid;
					tmp[used].mydefgpoolid = seedlist[used].mydefgpoolid;
					tmp[used].location = seedlist[used].location;
				}
				for ( uint32_t unused = nseedlist; unused < 2*nseedlist; unused++ ) {
					tmp[unused].ensdefgpoolid = NOT_ASSIGNED_YET;
					tmp[unused].mydefgpoolid = I_DONT_KNOW_YET;
					tmp[unused].location = NOT_ASSIGNED_YET;
				}
				delete [] seedlist;
				nseedlist = 2 * nseedlist;
				seedlist = tmp;
			} //remalloced
		}
		//next point
	}


	delete [] pp3;
	pp3 = NULL;

	//trim the seedlist to free unutilized memory
	defgseedP trimmed = NULL;
	trimmed = new struct defgseed[nvalidp];
	QUICKASSERT( trimmed != NULL );
	myMemGuard = myMemGuard + (nvalidp * sizeof(struct defgseed));

	for ( uint32_t used = 0; used < nvalidp; used++ ) {
		trimmed[used].ensdefgpoolid = seedlist[used].ensdefgpoolid;
		trimmed[used].mydefgpoolid = seedlist[used].mydefgpoolid;
		trimmed[used].location = seedlist[used].location;
	}
	delete [] seedlist;


	*howmanyseeds = nvalidp;

	return trimmed;
}

////////////SYNTHETIZATION WITH AUTOMATON
void caHdl::seed_deformedgrains_poisson( void )
{
	uint32_t np[1] = {0};
	tmpdefgseeds = NULL;
	tmpdefgseeds = init_discreteWindowFromMasterCSR( np );
	ndefgseeds = np[0];


//cout << "SEEDINGDEFORMEDGRAINSPOISSON " << ndefgseeds << endl;

	//MK::random sampling of deformed grains to generate a local deformation structure from worlddefgpool resulting in an uncorrelated MODF
	ensembleHdlP ens = this->myensHdl;
	uint32_t nensdgpool = ens->worlddefgpool.size();

//cout << "WorldensemblePoolSize=" << nensdgpool << endl;

	//start assigning the places in for loop but whiling to fine the place
	for ( uint32_t sd = 0; sd < ndefgseeds; sd++) {

		uint32_t luckygr = localprng.leEcuyer() * nensdgpool; //lucky grain is an ID in ensembleHdl defgpool array thus the orientation can be obtained by dereferencing the "ens->defgpool[luckyGrain].ori"-th entry of ens->oripool...

		//ASSIGNMENT IS WITH WORLDDEFGPOOL-ID, mind this is in contrast to mydefggrid where IDs refer to this->mydefgpool !
		tmpdefgseeds[sd].ensdefgpoolid = luckygr;

		//append deformed grain only to mydefgpool if we have not picked this grain already, therefore scan already known ones
		bool DoIKnowLuckyGrain = false;
		for ( uint32_t tg = 0; tg < sd; tg++) {
			if ( luckygr == tmpdefgseeds[tg].ensdefgpoolid ) { 
				DoIKnowLuckyGrain = true;
				break; //because obviously I know the grain
			}
		}

//cout << sd << "\tdoiknow=" << DoIKnowLuckyGrain << "\tpoolsize=" << mydefgpool.size() << "\tluckygr=" << luckygr << endl;

		if ( DoIKnowLuckyGrain == false ) { //most likely this condition is met when ensHdl->defgpool is large
			struct cadefg dgr;

			//get orientation copied from ensHdl
			uint32_t ensoriid = ens->worlddefgpool[luckygr].ori;
			double ensbunge[3];
			ensbunge[0] = ens->worldoripool[ensoriid].bunge1;
			ensbunge[1] = ens->worldoripool[ensoriid].bunge2;
			ensbunge[2] = ens->worldoripool[ensoriid].bunge3;

			//orientation recategorization with known Bunge values in local caHdl to utilize that not all possible orientations are used in the local CA grid thus a more efficient hashtable can be generated however at additional memory costs...
			dgr.worlddefgid = luckygr;
			dgr.caori = ca_check_disjunctness( ensbunge );
			dgr.cellcount = 0;
			dgr.rho0 = ens->worlddefgpool[luckygr].rho0;

			if (dgr.rho0 >= this->myrhomax) { //identify highest dislocation density in the system for adaptive integration scheme
				myrhomax = dgr.rho0;
			}
			dgr.rho = ens->worlddefgpool[luckygr].rho;
			dgr.dav0 = ens->worlddefgpool[luckygr].dav0;
			dgr.dav = ens->worlddefgpool[luckygr].dav;
			dgr.avdg0 = ens->worlddefgpool[luckygr].avdg0;
			dgr.avdg = ens->worlddefgpool[luckygr].avdg;

			mydefgpool.push_back ( dgr );

			//finally assign local copy of this grain to the pool
			tmpdefgseeds[sd].mydefgpoolid = mydefgpool.size() - 1;

//cout << "F-sd;luckyGrain;DoIKnowLuckyGrain;mydefgpoolid;" << sd << ";" << luckygr << ";" << DoIKnowLuckyGrain << ";" << mydefgpool[(mydefgpool.size() - 1)].worlddefgid << endl;

			continue;
		}

		//else {//obviously, DoIKnowLuckyGrain == true -> I picked the grain already during the loop
		//well then find it in mylist utilizing that worlddefgid are >= 0 and unique
		for (int testg = 0; testg < mydefgpool.size(); testg++) {
			if ( luckygr == mydefgpool[testg].worlddefgid ) {
				tmpdefgseeds[sd].mydefgpoolid = testg;
				break;
//cout << "T-sd;luckyGrain;DoIKnowLuckyGrain;testg;" << sd << ";" << luckygr << ";" << DoIKnowLuckyGrain << ";" << testg << endl;
			}
		}

	} //for all candidates


//	for (int sd = 0; sd < NumberOfPlannedSeeds; sd++) //##DEBUG small output
//		cout << "sd;loc;enspool;defpool\t" << sd << "\t" << tmpdefgseeds[sd].location << "\t" << tmpdefgseeds[sd].ensdefgpoolid << "\t" <<  tmpdefgseeds[sd].mydefgpoolid << endl;
}

/*
inline void caHdl::id2xyz_cpfem( int id, short * coor )
{
	short z = id / CPFEM_NGRAINSXY;
	int rem = id - (z * CPFEM_NGRAINSXY);
	short y = rem / CPFEM_NGRAINSX;
	short x = rem - (y * CPFEM_NGRAINSX);

	coor[0] = x;
	coor[1] = y;
	coor[2] = z;
}


inline int caHdl::xyz2id_cpfem( short * coor )
{
	int idxyz = coor[0] + (coor[1] * CPFEM_NGRAINSX) + (coor[2] * CPFEM_NGRAINSXY);
	return idxyz;
}
*/


inline uint32_t caHdl::get_damask_idperiodic( int oid, int dnx, int dny, int dnz )
{
	//interprets from a linear array mapping implicitly a 3D aggregate of CPFEM Damask grains ordered as x chain of grains stacked in y and these xy layers stacked in z
	//oid is the position of the grain in the linear array that makes the basis leftmost, bottommost, frontmost grain
	int z = oid / CPFEM_NGRAINSXY;
	int rem = oid - (z * CPFEM_NGRAINSXY);
	int y = rem / CPFEM_NGRAINSX;
	int x = rem - (y * CPFEM_NGRAINSX);

//cout << "Firstgrain at = " << x << ";" << y << ";" << z << endl;

	//identify the target grain
	x += dnx;
	y += dny;
	z += dnz;

	//##DEBUG
	QUICKASSERT ( x >= 0 && x < SHORT_MAX && y >= 0 && y < SHORT_MAX && z >= 0 && z < SHORT_MAX );

	//apply periodic boundary conditions
	//##MK::probably modulo is not correct, verified
	if ( x < 0 )				x += CPFEM_NGRAINSX;
	if ( x >= CPFEM_NGRAINSX )	x = x % CPFEM_NGRAINSX;
	if ( y < 0 )				y += CPFEM_NGRAINSY;
	if ( y >= CPFEM_NGRAINSY )	y = y % CPFEM_NGRAINSY;
	if ( z < 0 )				z += CPFEM_NGRAINSZ;
	if ( z >= CPFEM_NGRAINSZ )	z = z % CPFEM_NGRAINSZ;

	uint32_t idxyz = x + (y * CPFEM_NGRAINSX) + (z * CPFEM_NGRAINSXY);

//cout << "\t\t\t\tidperiodic\t" << dnx << ";" << dny << ";" << dnz << "\t\t" << x << ";" << y << ";" << z << "\t\t" << idxyz << endl;

	return idxyz;
}



void caHdl::grow_deformedgrains_initialization( void )
{
	//##MK::currently no defragmentation
}


void caHdl::grow_deformedgrains_init_mySeedFront( void )
{
	//allocate initial memory to store ACTIVE cells
	uint32_t ninitialSeedFront = (uint32_t) (initialRelCellCaching * (double) myCAGeometry.nboxvol_rdndtd);
	ninitialSeedFront = ninitialSeedFront + MINLENGTH_ACTIVELIST;
	QUICKASSERT ( ninitialSeedFront < CA_ALLOCATION_MAXIMUM );
	mySeedFront = NULL;
	mySeedFront = new struct defcell[ninitialSeedFront];
	QUICKASSERT( mySeedFront != NULL );
	myMemGuard = myMemGuard + (ninitialSeedFront * sizeof(defcell));

	ntotalSeedFront = ninitialSeedFront;
	nextSlotNeverActiveSeedFront = 0;
	nCurrentlyActiveSeedFront = 0;

	//initialize already associated cells
	for ( uint32_t sfr = 0; sfr < ninitialSeedFront; ++sfr) {
		mySeedFront[sfr].activity = INACTIVE;
		mySeedFront[sfr].infector = 26;
		mySeedFront[sfr].ix = -1;
		mySeedFront[sfr].iy = -1;
		mySeedFront[sfr].iz = -1;
		mySeedFront[sfr].frac = NO_INFECTION;
		mySeedFront[sfr].mydefgseedid = NO_GRAIN_ASSIGNED;
	}


	//allocate memory to store RecycCells
	uint32_t ninitialRecycList = (uint32_t) (initialRelCellCaching * (double) myCAGeometry.nboxvol_rdndtd * maxfillperstep * FULLRECYCLING_ATT_FACTOR);
	ninitialRecycList = ninitialRecycList + MINLENGTH_RECYCLIST;
	QUICKASSERT ( ninitialRecycList < CA_ALLOCATION_MAXIMUM );
	myRecycList = NULL;
	myRecycList = new uint32_t[ninitialRecycList];
	QUICKASSERT( myRecycList != NULL );
	myMemGuard = myMemGuard + ( ninitialRecycList * sizeof(uint32_t) );

	ntotalRecycList = ninitialRecycList;
	nextSlotThatBecomesRecyc = NOTHING_TO_RECYCLE; //not assigned as there is no cell at all yet  that can be recycled 
	firstNotRecycYet = NOTHING_TO_RECYCLE;

	//myRecycList needs no initialization because it is assured to have no gaps and only read in between [0;nextSlotThatBecomesRecycled)

	//allocate maxfillperstep the initial memory for cell that complete their infection in a timestep
	//this list serves as a guidance to reduce the amount of ifs during the infection phase as the list is for sure fragmented and thus INACTIVE CELLs are  met
	//active cells * fraction that become filled * heuristic factor
	uint32_t ninitialFullSeedList = (uint32_t) (initialRelCellCaching * (double) myCAGeometry.nboxvol_rdndtd * maxfillperstep * FULLRXCACHING_ATT_FACTOR);
	ninitialFullSeedList = ninitialFullSeedList + MINLENGTH_FULLSEEDLIST;
	QUICKASSERT ( ninitialFullSeedList < CA_ALLOCATION_MAXIMUM );
	myFullSeedList = NULL;
	myFullSeedList = new uint32_t[ninitialFullSeedList];
	QUICKASSERT( myFullSeedList != NULL );
	myMemGuard = myMemGuard + ( ninitialFullSeedList * sizeof(uint32_t) );

	ntotalFullSeedList = ninitialFullSeedList;
	nextSlotToFullSeed = 0; //as above not assigned yet

	//no initialization as well read-only assured gap-free [0;nextSlotToFullSeed)

	//cout << this->jobid << "\t\tMemGuard (Byte)= " << this->myMemGuard << endl;
}


void caHdl::grow_deformedgrains_append_to_fullseed_list( uint32_t seedfrontid )
{
	if ( nextSlotToFullSeed < ntotalFullSeedList ) { //FullSeed list still large enough, append
		myFullSeedList[nextSlotToFullSeed] = seedfrontid;		//indexing relative to first element in array with displacement sizeof(int)
		nextSlotToFullSeed++;
		return;
	}

	//FullSeed list is too small but it is necessary to store all cells that updateFullSeed should work on
	size_t oldSize = ntotalFullSeedList;
	size_t newSize = oldSize + (this->transientRelCellRecaching * oldSize);
	QUICKASSERT( newSize < CA_ALLOCATION_MAXIMUM );

//#ifdef REPORTSTYLE_DEVELOPER
//	cout << "\t\tALLOCATION for fullrx list, oldSize, newSize\t" << ((int) oldSize) << "\t" << ((int) newSize) << endl;
//#endif

	uint32_t* new_myFullSeedList = NULL;
	new_myFullSeedList = new uint32_t[newSize];
	QUICKASSERT( new_myFullSeedList != NULL );
	myMemGuard = myMemGuard + ( newSize * sizeof(uint32_t) );

	for ( size_t i = 0; i < oldSize; ++i) {
		new_myFullSeedList[i] = myFullSeedList[i];
	}

	delete [] myFullSeedList;
	myMemGuard = myMemGuard - ( oldSize * sizeof(uint32_t) );
	myFullSeedList = new_myFullSeedList;
	ntotalFullSeedList = newSize;

	QUICKASSERT ( nextSlotToFullSeed < ntotalFullSeedList );

	//more space again, so append
	myFullSeedList[nextSlotToFullSeed] = seedfrontid;
	nextSlotToFullSeed++;
}


void caHdl::grow_deformedgrains_append_to_recyc_list( uint32_t seedfrontid )
{
	if ( firstNotRecycYet < ntotalRecycList ) {
		myRecycList[firstNotRecycYet] = seedfrontid;
		firstNotRecycYet++;
	}
	//else -> currently I do not recycle more than ntotalRecycList elements, also this list is already as large as 0.15*ntotalcells, user can optimize this structure when knowing his microstructural path function Sv(X) in more details to save additional memory therefore not further memory utilized
}


uint32_t caHdl::grow_deformedgrains_NextFreeSlotInSeedFrontAppendOrRecycle( bool how )
{
	//##MK::two extreme strategies possible
	//1::always get more memory (almost no additional scans for places but increasing number of cache hits)
	//2::recycle as best as possible (reduced total amount of memory and less scans on likely improved cache performance than strategy 1)
	//optimal condition invokes compact (meaning no defcell in the bucket is INACTIVE) that is only possible with linked list structures are utilized however those are inefficient to traverse along and to parallelize
	
	//##MK::algorithm is able to select most beneficial strategy automatically
	bool strategy = how;
	uint32_t place;

	//if recycling is desired, return an index to a position in the CellList that is currently INACTIVE
	if ( strategy == CELLRECYCLE ) {
		//is there some place that can be recycled?
		if ( nextSlotThatBecomesRecyc < firstNotRecycYet ) { //##DEBUG//&& (firstNotRecycYet < ntotalRecycList)
			place = myRecycList[nextSlotThatBecomesRecyc];
			nextSlotThatBecomesRecyc++;

			//#ifdef REPORTSTYLE_CELLCYCLES
			//	cout << "\t\tRECYCLE\t" << place << endl;
			//#endif

			return place;
		}

		//hasnt returned yet, so recycling list should be utilized but was already exhausted so change the strategy...
		strategy = CELLAPPEND;
	}

	//now strategy is/or was CELLAPPEND, then take the next cell that was already initialize but never used ACTIVEly
	//if ( strategy == CELLAPPEND ) {
		//mySeedFront cached large enough?
		if ( nextSlotNeverActiveSeedFront < ntotalSeedFront ) {
			place = nextSlotNeverActiveSeedFront;
			nextSlotNeverActiveSeedFront++; // = nextSlotNeverActiveSeedFront + 1;

			//#ifdef REPORTSTYLE_CELLCYCLES
			//	cout << "\t\tAPPEND\t" << place << endl;
			//#endif

			return place;
		}

		//obviously no preinitialized cell cache left and mind that accesses such as mySeedFront[nextSlotNeverActiveSeedFront] are faulty when >= ntotalSeedFront...
		size_t oldSize = ntotalSeedFront;
		size_t newSize = oldSize + (transientRelCellRecaching * oldSize);
		QUICKASSERT ( newSize < CA_ALLOCATION_MAXIMUM );

//#ifdef REPORTSTYLE_DEVELOPER
//	cout << "\t\tALLOCATION for new_mySeedFront;oldSize;newSize\t" << ((int) oldSize) << "\t" << ((int) newSize) << endl;
//#endif
		defcellP new_mySeedFront = new struct defcell[newSize];
		myMemGuard = myMemGuard + (newSize * sizeof(cell));

		for ( size_t i = 0; i < oldSize; ++i) {
			new_mySeedFront[i].activity = mySeedFront[i].activity;
			new_mySeedFront[i].infector = mySeedFront[i].infector;
			new_mySeedFront[i].ix = mySeedFront[i].ix;
			new_mySeedFront[i].iy = mySeedFront[i].iy;
			new_mySeedFront[i].iz = mySeedFront[i].iz;
			new_mySeedFront[i].frac = mySeedFront[i].frac;
			new_mySeedFront[i].mydefgseedid = mySeedFront[i].mydefgseedid;
		}

		delete [] mySeedFront;
		myMemGuard = myMemGuard - (oldSize * sizeof(defcell));
		mySeedFront = new_mySeedFront;
		ntotalSeedFront = newSize;

		//initialize the trailing new memory snippet, access with nc = oldSize is to the 0-th element of the appended section and not yet initialized section
		for ( size_t nc = oldSize; nc < ntotalSeedFront; ++nc ) {
			mySeedFront[nc].activity = INACTIVE;
			mySeedFront[nc].infector = 26;
			mySeedFront[nc].ix = -1;
			mySeedFront[nc].iy = -1;
			mySeedFront[nc].iz = -1;
			mySeedFront[nc].frac = NO_INFECTION;
			mySeedFront[nc].mydefgseedid = NO_GRAIN_ASSIGNED;
		}

		//now we have new space, so append
		place = nextSlotNeverActiveSeedFront;
		nextSlotNeverActiveSeedFront++; // = nextSlotNeverActiveSeedFront + 1;

		//#ifdef REPORTSTYLE_CELLCYCLES
		//	cout << "\t\tNEWREGION\t" << place << endl;
		//#endif

		return place;
	//}
}


void caHdl::grow_deformedgrains_infect_JuvenileNucleation( uint32_t seedid, uint32_t WhichCellToInfect, double frac0 )
{
	//identify the grid location x,y,z of the infectorcell, 
	//MK:: works because positive range of short [CA_DIMENSIONS_MINIMUM, CA_DIMENSIONS_MAXIMUM]
	uint32_t z = WhichCellToInfect / myCAGeometry.nboxarea_rdnd;
	uint32_t rem = WhichCellToInfect - (myCAGeometry.nboxarea_rdnd * z);
	uint32_t y = rem / myCAGeometry.nboxedge_rd;
	uint32_t x = rem - (myCAGeometry.nboxedge_rd * y);

	//MK::if desired check domain contact here, but way faster to do it later when structure is synthetized
	
	//get out here if cell is already infected, or particle but in case of site-saturated synthetization no cell can be infected already, no boundary tracking here
	if ( (mycellgrid[WhichCellToInfect] != NOT_ASSIGNED_YET) || (mycellgrid[WhichCellToInfect] == CURRENTLY_INFECTED) ) { //particles are discretized after packing the grains
		return;
	}

	//obviously cell is free, hehe, go for it!
	uint32_t freeplace = grow_deformedgrains_NextFreeSlotInSeedFrontAppendOrRecycle( CELLAPPEND );

	//write down all cell information
	mySeedFront[freeplace].activity = ACTIVE;
	mySeedFront[freeplace].infector = 26;
	mySeedFront[freeplace].ix = (short) x;
	mySeedFront[freeplace].iy = (short) y;
	mySeedFront[freeplace].iz = (short) z;
	mySeedFront[freeplace].frac = frac0;
	mySeedFront[freeplace].mydefgseedid = seedid; //is replaced after completed boundary construction - boundaryTracking requires disjoint grains but grains from the ensemble pool might appear multiple times in CA - //tmpdefgseeds[seedid].mydefgpoolid;

	//P obsolete

	//mark cell as infected to prevent multiple reinfections, now only the front cell carries deformation state of the voxel
	mycellgrid[WhichCellToInfect] = CURRENTLY_INFECTED;

#ifdef REPORTSTYLE_DEFCELLCYCLES
	cout << "MS-Synthesis::" << this->jobid << "::JUVENILE;seedid;WhichCellToInfect;x;y;z;frac0;freeplace\t\t" << seedid << ";" << WhichCellToInfect << ";" << x << ";" << y << ";" << z << ";" << frac0 << ";" << freeplace << endl;
#endif
}


void caHdl::grow_deformedgrains_placeAllSeeds( void )
{
	for ( uint32_t sd = 0; sd < ndefgseeds; ++sd) {
		uint32_t plantWhere = tmpdefgseeds[sd].location;

		grow_deformedgrains_infect_JuvenileNucleation( sd, plantWhere, ALMOST_FULLY_INFECTED );

//cout << "Placing " << sd << " a seed at " << plantWhere << endl;
	}

	//if ( myensHdl->myRank == MASTER ) {
	//	cout << this->jobid << " site-saturated planting of seeds happened successfully ndefgseeds" << ndefgseeds << endl;
	//}
}


void caHdl::grow_deformedgrains_defragmentation( void )
{
	//##MK::currently no defragmentation
}


void caHdl::grow_deformedgrains_calcGrowthStep( void )
{
	//##consider to mount them fixed in the classHdl
	double shapefactors[27] = { FGEODIAG, FGEOEDGE, FGEODIAG, FGEOEDGE, FGEOFACE, FGEOEDGE, FGEODIAG, FGEOEDGE, FGEODIAG, 
						FGEOEDGE, FGEOFACE, FGEOEDGE, FGEOFACE,      FGEOFACE, FGEOEDGE, FGEOFACE, FGEOEDGE,	
						FGEODIAG, FGEOEDGE, FGEODIAG, FGEOEDGE, FGEOFACE, FGEOEDGE, FGEODIAG, FGEOEDGE, FGEODIAG, (1.0) }; //layers xz stacked in positive y direction...

	double dVrx = 0.0;
	uint32_t IdentifiedAsActive = 0;

	//reset collector array for already transformed cells this array is always dense and only interpreted between [0, nextSlotToFullSeed)
	nextSlotToFullSeed = 0;

	//##MK::reset RecycCounter as int-values from previous step are now corrupted, up to now nothing has been recyc
	nextSlotThatBecomesRecyc = 0;
	firstNotRecycYet = 0;

	//scan the whole mySeedFront on the interval [0;nextSlotNeverActiveSeedFront) to 
	//	a) to find ACTIVE cells for which boundary migration is simulated during dt
	//	b) to find slots that in the last integration step have become INACTIVE and thus can be reutilized by collecting them in the RecycList in an effort to keep the mySeedFront list compact ie cache local
//#ifdef REPORTSTYLE_CELLCYCLES
//	cout << "\n\n\t\tstep->calc_beforeloop;nextSlotNeverActiveSeedFront\t" << step << "\t" << nextSlotNeverActiveSeedFront << endl;
//#endif

	for ( uint32_t c = 0; c < nextSlotNeverActiveSeedFront; c++) {

		if ( mySeedFront[c].activity == ACTIVE ) { //list is dense and full so ACTIVE the most likely case to raise branch predictor efficiency
			IdentifiedAsActive++;

			//unit time infection fill increment, scale according to infection direction
			double v = maxfillperstep * shapefactors[mySeedFront[c].infector];

			mySeedFront[c].frac = mySeedFront[c].frac + (v * 1.0);

			dVrx = dVrx + (v * 1.0);

			//does this fill up the cell completely?
			if ( mySeedFront[c].frac >= 1.0 ) {
				//cout << "APPENDTOFULLRX\t\t" << c << ";" << myRXFront[c].rxFrac << endl;
				grow_deformedgrains_append_to_fullseed_list ( c );
			}

//#ifdef REPORTSTYLE_CELLCYCLES
//			cout << "\t\tstep->calc_beforecontinue is active;rxFrac;v;nactive\t" << step << "\t" << c << "\t" << myRXFront[c].rxFrac << "\t" << v << "\t" << nCurrentlyActive << endl;
//			//QUICKASSERT ( v > 0.0 );
//#endif
			//##DEBUG else construct more expensive but for testing purposes continue in the condition should have the same effect
			continue;
		} 
		//else {
			//else cell seems inactive, okay remember index c in RecycList

		grow_deformedgrains_append_to_recyc_list ( c );
		//}

	} //for each cell


	this->dXstep = dVrx;
	this->Sv = IdentifiedAsActive;
	this->nCurrentlyActive = IdentifiedAsActive;

#ifdef REPORTSTYLE_CELLCYCLES
	cout << "\t\tstep->calc_end;step;dXstep;Sv;\t" << step << "\t" << dXstep << "\t" << Sv << endl;
#endif
}


void caHdl::grow_deformedgrains_infect_OutofSeedFront_BoundsCheck( uint32_t defgid, uint32_t seedfrontid, short dx, short dy, short dz, unsigned char Direction, double CarryOver )
{
	//read the location of the infector
	short x = mySeedFront[seedfrontid].ix;
	short y = mySeedFront[seedfrontid].iy;
	short z = mySeedFront[seedfrontid].iz;

	//add relative coordinate increment to identify position of infection target in the CA, no overflow as -1 <= dx, dy, dz <= 1
	x += dx;
	y += dy;
	z += dz;
	
	//##MK::modelling approach DEFAULT periodic boundary conditions in each automaton, switch contact flag if applicable
	//no periodic boundary condition during defms synthesis, but of course out of the box infections are still forbidden
	//MK::if periodic boundary conditions desired link in periodic boundary conditions more statements like...
	//x += myCAGeometry.nboxedge_rd;
	//myrxgpoolboundarycontact[rxgpoolid] = true;
	if ( x < 0 )
		return;
	if ( x >= myCAGeometry.nboxedge_rd )
		return;
	if ( y < 0 )
		return;
	if ( y >= myCAGeometry.nboxedge_nd )
		return;
	if ( z < 0 )
		return;
	if ( z >= myCAGeometry.nboxedge_td )
		return;

	//convert ix,iy,iz in implicit 3D coordinate to access mycellgrid
	uint32_t cxyz = x + (y * myCAGeometry.nboxedge_rd) + (z * myCAGeometry.nboxarea_rdnd);

	uint32_t seedtarget = mycellgrid[cxyz];
	uint32_t seedinfector = defgid;

	if ( seedtarget != NOT_ASSIGNED_YET ) { ////MK::particle become discretized later ##DEBUG|| seedtarget == CURRENTLY_INFECTED|| seedtarget == CELL_IS_A_PARTICLE ) { 
		//MK::the cell is either completely transformed already or still infected, anyhow a grain has already been assigned

		/*//MK::tracking boundaries on fly can be done, but is not as accurate therefore here performed afterwards 
		//#ifdef REPORTSTYLE_CELLCYCLES
		//cout << "\t\t\t\tREJECT infection cxyz;x;y;z;dx;dy;dz\t" << cxyz << "\t" << x << ";" << y << ";" << z << "\t" << dx << ";" << dy << ";" << dz << endl;
		//#endif*/

		//anyhow infection attempt is rejected
		return;
	}

	//obviously cell is free NOT_ASSIGNED_YET, hehe, go for it!
	uint32_t freeplace = grow_deformedgrains_NextFreeSlotInSeedFrontAppendOrRecycle ( CELLRECYCLE ); //##MK::strategy is RECYCLING because OutofSeedFront infection only called from updateFullSeed

	//newly infected cell at freeplace inherits properties from infector from mySeedFront location seedfrontid
	mySeedFront[freeplace].activity = ACTIVE;
	mySeedFront[freeplace].infector = Direction;
	mySeedFront[freeplace].ix = x;
	mySeedFront[freeplace].iy = y;
	mySeedFront[freeplace].iz = z;
	mySeedFront[freeplace].frac = CarryOver;
	//cell carries consumed orientation along
	mySeedFront[freeplace].mydefgseedid = defgid;


	//mark cell as infected to prevent multiple reinfections
	mycellgrid[cxyz] = CURRENTLY_INFECTED;

	//assessor function grow_deformedgrains_NextFreeSlotInSeedFrontAppendOrRecycle() ASSURES freeplace < nextSlotNeverActiveSeedFront <= ntotalSeedFront!
#ifdef REPORTSTYLE_CELLCYCLES
	cout << this->jobid << "\t\tout-of-front infection;defgid;seedfrontid;cxyz;x;y;z;freeplace" << defgid << ";" << seedfrontid << ";" << cxyz << ";" << x << ";" << y << ";" << z << ";" << freeplace << endl;
#endif
}


void caHdl::grow_deformedgrains_infect_OutofSeedFront_NoBoundsCheck( uint32_t defgid, uint32_t seedfrontid, short dx, short dy, short dz, unsigned char Direction, double CarryOver )
{
	//read the location of the infector
	short x = mySeedFront[seedfrontid].ix;
	short y = mySeedFront[seedfrontid].iy;
	short z = mySeedFront[seedfrontid].iz;

	//add relative coordinate increment to identify position of infection target in the CA, no overflow as -1 <= dx, dy, dz <= 1
	x += dx;
	y += dy;
	z += dz;

	//MK::MUSTNT BE CALLED WHEN IT CANNOT BE ASSURED THAT x,y,z stay in the automaton domain

	//convert ix,iy,iz in implicit 3D coordinate to access mycellgrid
	uint32_t cxyz = x + (y * myCAGeometry.nboxedge_rd) + (z * myCAGeometry.nboxarea_rdnd);

	uint32_t seedtarget = mycellgrid[cxyz];
	uint32_t seedinfector = defgid;

	if ( seedtarget != NOT_ASSIGNED_YET ) { ////MK::particle become discretized later ##DEBUG|| seedtarget == CURRENTLY_INFECTED|| seedtarget == CELL_IS_A_PARTICLE ) { 
		//MK::the cell is either completely transformed already or still infected, anyhow a grain has already been assigned

		/*//MK::tracking boundaries on fly can be done, but is not as accurate therefore here performed afterwards 
		//#ifdef REPORTSTYLE_CELLCYCLES
		//cout << "\t\t\t\tREJECT infection cxyz;x;y;z;dx;dy;dz\t" << cxyz << "\t" << x << ";" << y << ";" << z << "\t" << dx << ";" << dy << ";" << dz << endl;
		//#endif*/

		//anyhow infection attempt is rejected
		return;
	}

	//obviously cell is free NOT_ASSIGNED_YET, hehe, go for it!
	uint32_t freeplace = grow_deformedgrains_NextFreeSlotInSeedFrontAppendOrRecycle ( CELLRECYCLE ); //##MK::strategy is RECYCLING because OutofSeedFront infection only called from updateFullSeed

	//newly infected cell at freeplace inherits properties from infector from mySeedFront location seedfrontid
	mySeedFront[freeplace].activity = ACTIVE;
	mySeedFront[freeplace].infector = Direction;
	mySeedFront[freeplace].ix = x;
	mySeedFront[freeplace].iy = y;
	mySeedFront[freeplace].iz = z;
	mySeedFront[freeplace].frac = CarryOver;
	//cell carries consumed orientation along
	mySeedFront[freeplace].mydefgseedid = defgid;


	//mark cell as infected to prevent multiple reinfections
	mycellgrid[cxyz] = CURRENTLY_INFECTED;

	//assessor function grow_deformedgrains_NextFreeSlotInSeedFrontAppendOrRecycle() ASSURES freeplace < nextSlotNeverActiveSeedFront <= ntotalSeedFront!
#ifdef REPORTSTYLE_CELLCYCLES
	cout << this->jobid << "\t\tout-of-front infection;defgid;seedfrontid;cxyz;x;y;z;freeplace" << defgid << ";" << seedfrontid << ";" << cxyz << ";" << x << ";" << y << ";" << z << ";" << freeplace << endl;
#endif
}



void caHdl::grow_deformedgrains_updateFullCells( void )
{
	//##MK::we can rely on the fact that updateFullSeed is called after calcGrowthStep such that
	//myFullSeedList in the interval [0;nextSlotToFullSeed) is a compact list of entries dereferencing cells from mySeedFront that are completely filled, infect others and then be switched off for recycling

	//currently //compromise //originalfactors //seemed optimized factors but are not //original factors for velocity fgeo 1 in <100>, 1/2^0.5 <110>, 1/3^0.5 <111>
	const double overFac1 = 0.0682; //0.04946756; //0.03738274; //0.0523448297609848;  //###MK20121211, 26NN, 8 face, 12 edge and 6 diag neighbors each one getting a portion of the overshoot, all other overshoot is discarded
	const double overFac2 = 0.0542; //0.0387; //0.03949696; //0.04059546; //0.0370133840840478;
	const double overFac3 = 0.0484; //0.0307; //0.02865390; //0.03606975; //0.0302213015531897;
	double overshoot, carryOver1, carryOver2, carryOver3;

	uint32_t nUpdatedCells = 0;
	uint32_t mycageo_nboxedge_rd = myCAGeometry.nboxedge_rd;

	uint32_t lowerlimit = 2;
	uint32_t upperlimit = mycageo_nboxedge_rd - 2;

	uint32_t mycageo_nboxarea_rdnd = myCAGeometry.nboxarea_rdnd;
	unsigned char location;

	for ( uint32_t c = 0; c < nextSlotToFullSeed; c++) {
		uint32_t cseedfrontid = myFullSeedList[c];

		uint32_t dgid = mySeedFront[cseedfrontid].mydefgseedid;

		//####BOOKKEEP DIRECTLY WITH mydefgseedids!

		//#ifdef REPORTSTYLE_CELLCYCLES
		//	cout << "\t\t\t->step->updt_priorinfects;c;crxfrontid;myRXfront[crxfrontid].rxFrac;nextSlotToFullRX" << step << ";" << c << ";" << crxfrontid << ";" << myRXFront[crxfrontid].rxFrac << ";" << nextSlotToFullRX << endl;
		//#endif
		QUICKASSERT( mySeedFront[cseedfrontid].frac >= 1.0 ); //##DEBUG

		overshoot = mySeedFront[cseedfrontid].frac - 1.0;

		//partition normalized overshoot
		carryOver1 = overshoot * overFac1;
		carryOver2 = overshoot * overFac2;
		carryOver3 = overshoot * overFac3;

		//MK::update shape tracking if desired here...

		//infect all untransformed neighbors
		//MK::ASSURE THAT COMPILER OPTIMIZATIONS DO NOT PARALLELIZE THIS PART SO THAT IN A WORST CASE ALL INFECT FUNCTIONS TRY TO RECYCLE THE SAME CELL STRUCT IN THE MYFRONTSLIST!
		/*
		//start with von Neumann neighborhood that has "six faces, six faces is all that we have ..."
		grow_deformedgrains_infect_OutofSeedFront ( dgid, cseedfrontid, -1, 0, 0, 13, carryOver1 );
		grow_deformedgrains_infect_OutofSeedFront ( dgid, cseedfrontid, 0, -1, 0, 15, carryOver1 );
		grow_deformedgrains_infect_OutofSeedFront ( dgid, cseedfrontid, 0, 0, -1, 21, carryOver1 );
		grow_deformedgrains_infect_OutofSeedFront ( dgid, cseedfrontid, +1, 0, 0, 12, carryOver1 );
		grow_deformedgrains_infect_OutofSeedFront ( dgid, cseedfrontid, 0, +1, 0, 10, carryOver1 );
		grow_deformedgrains_infect_OutofSeedFront ( dgid, cseedfrontid, 0, 0, +1, 4, carryOver1 );

		//additional Moore neighborhood
		//"twelve edges, twelve edges, between you and me..."
		grow_deformedgrains_infect_OutofSeedFront ( dgid, cseedfrontid, +1, +1, 0, 9, carryOver2 ); //xy
		grow_deformedgrains_infect_OutofSeedFront ( dgid, cseedfrontid, -1, +1, 0, 11, carryOver2 );
		grow_deformedgrains_infect_OutofSeedFront ( dgid, cseedfrontid, +1, -1, 0, 14, carryOver2 );
		grow_deformedgrains_infect_OutofSeedFront ( dgid, cseedfrontid, -1, -1, 0, 16, carryOver2 );

		grow_deformedgrains_infect_OutofSeedFront ( dgid, cseedfrontid, 0, +1, +1, 1, carryOver2 ); //yz
		grow_deformedgrains_infect_OutofSeedFront ( dgid, cseedfrontid, 0, -1, +1, 7, carryOver2 );
		grow_deformedgrains_infect_OutofSeedFront ( dgid, cseedfrontid, 0, +1, -1, 18, carryOver2 );
		grow_deformedgrains_infect_OutofSeedFront ( dgid, cseedfrontid, 0, -1, -1, 24, carryOver2 );

		grow_deformedgrains_infect_OutofSeedFront ( dgid, cseedfrontid, +1, 0, +1, 3, carryOver2 ); //xz
		grow_deformedgrains_infect_OutofSeedFront ( dgid, cseedfrontid, -1, 0, +1, 5, carryOver2 );
		grow_deformedgrains_infect_OutofSeedFront ( dgid, cseedfrontid, +1, 0, -1, 20, carryOver2 );
		grow_deformedgrains_infect_OutofSeedFront ( dgid, cseedfrontid, -1, 0, -1, 22, carryOver2 );

		//"eight diagonals, ten is all that's been done..."
		grow_deformedgrains_infect_OutofSeedFront ( dgid, cseedfrontid, +1, +1, +1, 0, carryOver3 ); //xyz
		grow_deformedgrains_infect_OutofSeedFront ( dgid, cseedfrontid, -1, +1, +1, 2, carryOver3 );
		grow_deformedgrains_infect_OutofSeedFront ( dgid, cseedfrontid, +1, -1, +1, 6, carryOver3 );
		grow_deformedgrains_infect_OutofSeedFront ( dgid, cseedfrontid, -1, -1, +1, 8, carryOver3 );
		grow_deformedgrains_infect_OutofSeedFront ( dgid, cseedfrontid, +1, +1, -1, 17, carryOver3 );
		grow_deformedgrains_infect_OutofSeedFront ( dgid, cseedfrontid, -1, +1, -1, 19, carryOver3 );
		grow_deformedgrains_infect_OutofSeedFront ( dgid, cseedfrontid, +1, -1, -1, 23, carryOver3 );
		grow_deformedgrains_infect_OutofSeedFront ( dgid, cseedfrontid, -1, -1, -1, 25, carryOver3 );
		//##MK::following code NEEDS to be executed AFTER all these infection attempts
		*/


		//categorize location
		location = INFECTOR_DEEP_IN_THE_CUBE; //MK::assuming cubic domain!
		if ( mySeedFront[cseedfrontid].ix < lowerlimit || mySeedFront[cseedfrontid].ix > upperlimit ) location = INFECTOR_CLOSE_TO_THE_DOMAINWALL;
		if ( mySeedFront[cseedfrontid].iy < lowerlimit || mySeedFront[cseedfrontid].iy > upperlimit ) location = INFECTOR_CLOSE_TO_THE_DOMAINWALL;
		if ( mySeedFront[cseedfrontid].iz < lowerlimit || mySeedFront[cseedfrontid].iz > upperlimit ) location = INFECTOR_CLOSE_TO_THE_DOMAINWALL;

		if ( location == INFECTOR_DEEP_IN_THE_CUBE ) { //most likely case

			//MK::restructured access pattern compared to original version for improved data locality
			//MK::any caches in the infection order of a CA algorithm affects the assignment of the voxel owing to the first-come-first serve principle of these algorithms and as such require care
			grow_deformedgrains_infect_OutofSeedFront_NoBoundsCheck ( dgid, cseedfrontid, -1, -1, -1, 25, carryOver3 );
			grow_deformedgrains_infect_OutofSeedFront_NoBoundsCheck ( dgid, cseedfrontid, 0, -1, -1, 24, carryOver2 );
			grow_deformedgrains_infect_OutofSeedFront_NoBoundsCheck ( dgid, cseedfrontid, +1, -1, -1, 23, carryOver3 );
			grow_deformedgrains_infect_OutofSeedFront_NoBoundsCheck ( dgid, cseedfrontid, -1, 0, -1, 22, carryOver2 );
			grow_deformedgrains_infect_OutofSeedFront_NoBoundsCheck ( dgid, cseedfrontid, 0, 0, -1, 21, carryOver1 );
			grow_deformedgrains_infect_OutofSeedFront_NoBoundsCheck ( dgid, cseedfrontid, +1, 0, -1, 20, carryOver2 );
			grow_deformedgrains_infect_OutofSeedFront_NoBoundsCheck ( dgid, cseedfrontid, -1, +1, -1, 19, carryOver3 );
			grow_deformedgrains_infect_OutofSeedFront_NoBoundsCheck ( dgid, cseedfrontid, 0, +1, -1, 18, carryOver2 );
			grow_deformedgrains_infect_OutofSeedFront_NoBoundsCheck ( dgid, cseedfrontid, +1, +1, -1, 17, carryOver3 );

			grow_deformedgrains_infect_OutofSeedFront_NoBoundsCheck ( dgid, cseedfrontid, -1, -1, 0, 16, carryOver2 );
			grow_deformedgrains_infect_OutofSeedFront_NoBoundsCheck ( dgid, cseedfrontid, 0, -1, 0, 15, carryOver1 );
			grow_deformedgrains_infect_OutofSeedFront_NoBoundsCheck ( dgid, cseedfrontid, +1, -1, 0, 14, carryOver2 );
			grow_deformedgrains_infect_OutofSeedFront_NoBoundsCheck ( dgid, cseedfrontid, -1, 0, 0, 13, carryOver1 );
			grow_deformedgrains_infect_OutofSeedFront_NoBoundsCheck ( dgid, cseedfrontid, +1, 0, 0, 12, carryOver1 );
			grow_deformedgrains_infect_OutofSeedFront_NoBoundsCheck ( dgid, cseedfrontid, -1, +1, 0, 11, carryOver2 );
			grow_deformedgrains_infect_OutofSeedFront_NoBoundsCheck ( dgid, cseedfrontid, 0, +1, 0, 10, carryOver1 );
			grow_deformedgrains_infect_OutofSeedFront_NoBoundsCheck ( dgid, cseedfrontid, +1, +1, 0, 9, carryOver2 );

			grow_deformedgrains_infect_OutofSeedFront_NoBoundsCheck ( dgid, cseedfrontid, -1, -1, +1, 8, carryOver3 );
			grow_deformedgrains_infect_OutofSeedFront_NoBoundsCheck ( dgid, cseedfrontid, 0, -1, +1, 7, carryOver2 );
			grow_deformedgrains_infect_OutofSeedFront_NoBoundsCheck ( dgid, cseedfrontid, +1, -1, +1, 6, carryOver3 );
			grow_deformedgrains_infect_OutofSeedFront_NoBoundsCheck ( dgid, cseedfrontid, -1, 0, +1, 5, carryOver2 );
			grow_deformedgrains_infect_OutofSeedFront_NoBoundsCheck ( dgid, cseedfrontid, 0, 0, +1, 4, carryOver1 );
			grow_deformedgrains_infect_OutofSeedFront_NoBoundsCheck ( dgid, cseedfrontid, +1, 0, +1, 3, carryOver2 );
			grow_deformedgrains_infect_OutofSeedFront_NoBoundsCheck ( dgid, cseedfrontid, -1, +1, +1, 2, carryOver3 );
			grow_deformedgrains_infect_OutofSeedFront_NoBoundsCheck ( dgid, cseedfrontid, 0, +1, +1, 1, carryOver2 );
			grow_deformedgrains_infect_OutofSeedFront_NoBoundsCheck ( dgid, cseedfrontid, +1, +1, +1, 0, carryOver3 );

			//deactivate cell 
			uint32_t cxyz = mySeedFront[cseedfrontid].ix + (mySeedFront[cseedfrontid].iy * mycageo_nboxedge_rd) + (mySeedFront[cseedfrontid].iz * mycageo_nboxarea_rdnd);

			//mark representor in the cellgrid as FULLY_SYNTHETIZED by assigning id of seed deformed grain
			mycellgrid[cxyz] = dgid;

			//texture bookkeeping not necessary as it completely synthetized structure is scanned later...

			//cell flagged as "Freiwild" and thought of Green Mile healed from all infections...
			mySeedFront[cseedfrontid].activity = INACTIVE;
			mySeedFront[cseedfrontid].frac = 0.0;

	//#ifdef REPORTSTYLE_CELLCYCLES
	//		cout << "\t\tstep->updt_afterinfects;cseedfrontid;activity;ACTIVEwouldbeMarkedAs;nextSeedSlotNeverActive\t" << step << "\t" << cseedfrontid << "\t" << mySeedFront[cseedfrontid].activity << "\t" << ACTIVE << "\t" << nextSlotNeverActiveSeedFront << endl;
	//#endif

			nUpdatedCells++;

			continue;
		}

		//INFECTOR_CLOSE_TO_THE_DOMAINWALL
			//MK::restructured access pattern compared to original version for improved data locality
			//MK::any caches in the infection order of a CA algorithm affects the assignment of the voxel owing to the first-come-first serve principle of these algorithms and as such require care
			grow_deformedgrains_infect_OutofSeedFront_BoundsCheck ( dgid, cseedfrontid, -1, -1, -1, 25, carryOver3 );
			grow_deformedgrains_infect_OutofSeedFront_BoundsCheck ( dgid, cseedfrontid, 0, -1, -1, 24, carryOver2 );
			grow_deformedgrains_infect_OutofSeedFront_BoundsCheck ( dgid, cseedfrontid, +1, -1, -1, 23, carryOver3 );
			grow_deformedgrains_infect_OutofSeedFront_BoundsCheck ( dgid, cseedfrontid, -1, 0, -1, 22, carryOver2 );
			grow_deformedgrains_infect_OutofSeedFront_BoundsCheck ( dgid, cseedfrontid, 0, 0, -1, 21, carryOver1 );
			grow_deformedgrains_infect_OutofSeedFront_BoundsCheck ( dgid, cseedfrontid, +1, 0, -1, 20, carryOver2 );
			grow_deformedgrains_infect_OutofSeedFront_BoundsCheck ( dgid, cseedfrontid, -1, +1, -1, 19, carryOver3 );
			grow_deformedgrains_infect_OutofSeedFront_BoundsCheck ( dgid, cseedfrontid, 0, +1, -1, 18, carryOver2 );
			grow_deformedgrains_infect_OutofSeedFront_BoundsCheck ( dgid, cseedfrontid, +1, +1, -1, 17, carryOver3 );

			grow_deformedgrains_infect_OutofSeedFront_BoundsCheck ( dgid, cseedfrontid, -1, -1, 0, 16, carryOver2 );
			grow_deformedgrains_infect_OutofSeedFront_BoundsCheck ( dgid, cseedfrontid, 0, -1, 0, 15, carryOver1 );
			grow_deformedgrains_infect_OutofSeedFront_BoundsCheck ( dgid, cseedfrontid, +1, -1, 0, 14, carryOver2 );
			grow_deformedgrains_infect_OutofSeedFront_BoundsCheck ( dgid, cseedfrontid, -1, 0, 0, 13, carryOver1 );
			grow_deformedgrains_infect_OutofSeedFront_BoundsCheck ( dgid, cseedfrontid, +1, 0, 0, 12, carryOver1 );
			grow_deformedgrains_infect_OutofSeedFront_BoundsCheck ( dgid, cseedfrontid, -1, +1, 0, 11, carryOver2 );
			grow_deformedgrains_infect_OutofSeedFront_BoundsCheck ( dgid, cseedfrontid, 0, +1, 0, 10, carryOver1 );
			grow_deformedgrains_infect_OutofSeedFront_BoundsCheck ( dgid, cseedfrontid, +1, +1, 0, 9, carryOver2 );

			grow_deformedgrains_infect_OutofSeedFront_BoundsCheck ( dgid, cseedfrontid, -1, -1, +1, 8, carryOver3 );
			grow_deformedgrains_infect_OutofSeedFront_BoundsCheck ( dgid, cseedfrontid, 0, -1, +1, 7, carryOver2 );
			grow_deformedgrains_infect_OutofSeedFront_BoundsCheck ( dgid, cseedfrontid, +1, -1, +1, 6, carryOver3 );
			grow_deformedgrains_infect_OutofSeedFront_BoundsCheck ( dgid, cseedfrontid, -1, 0, +1, 5, carryOver2 );
			grow_deformedgrains_infect_OutofSeedFront_BoundsCheck ( dgid, cseedfrontid, 0, 0, +1, 4, carryOver1 );
			grow_deformedgrains_infect_OutofSeedFront_BoundsCheck ( dgid, cseedfrontid, +1, 0, +1, 3, carryOver2 );
			grow_deformedgrains_infect_OutofSeedFront_BoundsCheck ( dgid, cseedfrontid, -1, +1, +1, 2, carryOver3 );
			grow_deformedgrains_infect_OutofSeedFront_BoundsCheck ( dgid, cseedfrontid, 0, +1, +1, 1, carryOver2 );
			grow_deformedgrains_infect_OutofSeedFront_BoundsCheck ( dgid, cseedfrontid, +1, +1, +1, 0, carryOver3 );

			//deactivate cell 
			uint32_t cxyz = mySeedFront[cseedfrontid].ix + (mySeedFront[cseedfrontid].iy * mycageo_nboxedge_rd) + (mySeedFront[cseedfrontid].iz * mycageo_nboxarea_rdnd);

			//mark representor in the cellgrid as FULLY_SYNTHETIZED by assigning id of seed deformed grain
			mycellgrid[cxyz] = dgid;

			//texture bookkeeping not necessary as it completely synthetized structure is scanned later...

			//cell flagged as "Freiwild" and thought of Green Mile healed from all infections...
			mySeedFront[cseedfrontid].activity = INACTIVE;
			mySeedFront[cseedfrontid].frac = 0.0;

	//#ifdef REPORTSTYLE_CELLCYCLES
	//		cout << "\t\tstep->updt_afterinfects;cseedfrontid;activity;ACTIVEwouldbeMarkedAs;nextSeedSlotNeverActive\t" << step << "\t" << cseedfrontid << "\t" << mySeedFront[cseedfrontid].activity << "\t" << ACTIVE << "\t" << nextSlotNeverActiveSeedFront << endl;
	//#endif

			nUpdatedCells++;

	} //for all cells

	Xcells += (double) nUpdatedCells;
	X = (Xcells / ((double) myCAGeometry.nboxvol_rdndtd));
}


void caHdl::grow_deformedgrains_cleanCellManagementDefMS( void )
{
	//MK::mycellgrid MUST NOT BE DELETED it is required for the subsequent growth simulation!

	delete [] myFullSeedList;
	myFullSeedList = NULL;
	myMemGuard = myMemGuard - ( ntotalFullSeedList * sizeof(uint32_t) );
	ntotalFullSeedList = 0;

	delete [] myRecycList;
	myRecycList = NULL;
	myMemGuard = myMemGuard - ( ntotalRecycList * sizeof(uint32_t) );
	ntotalRecycList = 0;

	delete [] mySeedFront;
	myRXFront = NULL;
	myMemGuard = myMemGuard - ( ntotalSeedFront * sizeof(defcell) );
	ntotalSeedFront = 0;
}


void caHdl::grow_deformedgrains_voxelize_poisson( void )
{
	//generate list of deformed grains and where the seeds as barycenters should be placed disjointly in mycellgrid
	seed_deformedgrains_poisson(); //either poisson (seeds located in the CA managed by tmpdefgseeds, dyn. GB tracking) or damask (voxelizing, subsequent GB tracking)

	grow_deformedgrains_initialization();
	grow_deformedgrains_init_mySeedFront();

	t = 0.0;
	X = 0.0;
	Xcells = 0.0;
	dXstep = 0.0;
	Sv = 0;
	step = 0;

	//MK::I/O prior to simulation desired? --> inject code here
	//no physical properties and mobilities necessary because constant velocity distance tessellation of space
	dt = 1.0;

	//add nuclei into the simulation, per se always site saturated
	grow_deformedgrains_placeAllSeeds();

	do {
		dXstep = 0.0;

		grow_deformedgrains_defragmentation();

		//fill ACTIVE cells while collecting INACTIVE cells in myRecycList or Matthias Loeck's equivalent garbage collection...
		grow_deformedgrains_calcGrowthStep();

		//let ACTIVE infect new cells
		grow_deformedgrains_updateFullCells();

#ifdef REPORTSTYLE_DETAILS
	cout << "DeformationSynthesis-Rank;" << myensHdl->myRank << "-Jobid;" << this->jobid << "-Step;X\t\t" << step << "\t" << X << endl;
#endif

		//prepare for next timestep, same integration order than START1D/3DCA and COReV3
		t += dt; 

		//###nothing to update physically here
		step++;
	}
	while ( (Xcells < myCAGeometry.nboxvol_rdndtd) && (step < DEFAULT_NTSTEPSMAX_DEFSYNTH) );

	//MK::optional final output? inject code here

	cout << myensHdl->myRank << " myRank,jobid = " << this->jobid << " local CA deformation structure POISSONVORONOI was synthetized with " << ndefgseeds  << endl;

	//##MK::output deformation structure
	//write_voxeldata_coloring_grainids();


	grow_deformedgrains_cleanCellManagementDefMS();
}


void caHdl::grow_deformedgrains_voxelize_damask( void )
{
	//MK::initialize the grain which should be considered
	//populates tmpdefgseeds as a implicit grid of deformed grains that can subsequently be voxelized directly without using an automaton
	//basic idea sample a complete block of grains applying periodic boundary conditions starting with a random grain from the input list
	//the input list is required to encode the data as xlines of length CPFEM_NGRAINSX stacking in CPFEM_NGRAINSY in y direction forming NGRAINSZ stack in positive z direction
	//pack 3D aggregates starting from the left-and-front-and-bottommost edge of the cuboid grain t displaced about the CA origin thus embedding the CA volume/box completely in the cuboid aggregate
	double tx = localprng.leEcuyer() * -1.0 * CPFEM_GRAINSIZEX;
	double ty = localprng.leEcuyer() * -1.0 * CPFEM_GRAINSIZEY;
	double tz = localprng.leEcuyer() * -1.0 * CPFEM_GRAINSIZEZ;

	ensembleHdlP ens = this->myensHdl;
	uint32_t nensdgpool = ens->worlddefgpool.size();

	int firstgrain = localprng.leEcuyer() * nensdgpool;

	//how many grains are necessary?
	int ngrx = ((myCAGeometry.boxedge_rd - tx) / CPFEM_GRAINSIZEX);
	ngrx++; //will be zero if defg+u larger than box, so always add one grain
	int ngry = ((myCAGeometry.boxedge_nd - ty) / CPFEM_GRAINSIZEY);
	ngry++;
	int ngrz = ((myCAGeometry.boxedge_td - tz) / CPFEM_GRAINSIZEZ);
	ngrz++;

	QUICKASSERT( ngrx < SHORT_MAX && ngry < SHORT_MAX && ngrz < SHORT_MAX );

cout << "DAMASK=tx;ty;tz;ngrx;ngry;ngrz\t\t" << tx << ";" << ty << ";" << tz << ";" << ngrx << ";" << ngry << ";" << ngrz << "\tfirstgrain=" << firstgrain << endl;

	int ngrxy = ngrx * ngry;
	int ngrxyz = ngrxy * ngrz;
	ndefgseeds = ngrxyz;
	QUICKASSERT ( ndefgseeds < CA_GRAINIDS_MAXIMUM );

	//allocate memory for these seeds
	tmpdefgseeds = NULL;
	tmpdefgseeds = new struct defgseed[ndefgseeds];
	QUICKASSERT ( tmpdefgseeds != NULL );
	myMemGuard = myMemGuard + (ndefgseeds * sizeof(struct defgseed));

	for ( uint32_t s = 0; s < ndefgseeds; ++s) {
		tmpdefgseeds[s].ensdefgpoolid = NOT_ASSIGNED_YET;
		tmpdefgseeds[s].mydefgpoolid = I_DONT_KNOW_YET;
		tmpdefgseeds[s].location = NOT_ASSIGNED_YET;
	}

	int sd = 0;
	uint32_t luckygr;
	for ( int zz = 0; zz < ngrz; zz++ ) {
		for ( int yy = 0; yy < ngry; yy++ ) {
			for ( int xx = 0; xx < ngrx; xx++ ) {

				luckygr = get_damask_idperiodic( firstgrain,  xx, yy, zz );

				//as the 3D aggregate of CPFEM grains is larger than and enclosing completely the CA simulation domain
				//starting distance tessellation from the starting point location does not work as the barycenter of some grains
				//might not lay in the domain at all

				//ASSIGNMENT IS WITH WORLDDEFGPOOL-ID, mind this is in contrast to mydefggrid where IDs refer to this->mydefgpool !
				tmpdefgseeds[sd].ensdefgpoolid = luckygr;
				tmpdefgseeds[sd].location = SUBSEQ_VOXELIZATION;

				//append deformed grain only to mydefgpool if we have not picked this grain already, therefore scan already known ones
				bool DoIKnowLuckyGrain = false;
				for (int tg = 0; tg < sd; tg++) {
					if ( luckygr == tmpdefgseeds[tg].ensdefgpoolid ) {
						DoIKnowLuckyGrain = true; 
						break;
					}
				}

//cout << "sd/zz/yy/xx/luckygr/known/ = " << sd << "--" << zz << ";" << yy << ";" << xx << "--" << luckygr << ";" << DoIKnowLuckyGrain << endl;

				if ( DoIKnowLuckyGrain == false ) { //most likely this condition is met when ensHdl->defgpool is large
					struct cadefg dgr;

					//get orientation copied from ensHdl
					uint32_t ensoriid = ens->worlddefgpool[luckygr].ori;
					double ensbunge[3];
					ensbunge[0] = ens->worldoripool[ensoriid].bunge1;
					ensbunge[1] = ens->worldoripool[ensoriid].bunge2;
					ensbunge[2] = ens->worldoripool[ensoriid].bunge3;

					//orientation recategorization with known Bunge values in local caHdl to utilize that not all possible orientations are used in the local CA grid thus a more efficient hashtable can be generated however at additional memory costs...
					dgr.worlddefgid = luckygr;
					dgr.caori = ca_check_disjunctness( ensbunge );
					dgr.cellcount = 0;
					dgr.rho0 = ens->worlddefgpool[luckygr].rho0;

					if (dgr.rho0 >= this->myrhomax) { //identify highest dislocation density in the system for adaptive integration scheme
						myrhomax = dgr.rho0;
					}
					dgr.rho = ens->worlddefgpool[luckygr].rho;
					dgr.avdg0 = ens->worlddefgpool[luckygr].avdg0;
					dgr.avdg = ens->worlddefgpool[luckygr].avdg;
					dgr.dav0 = ens->worlddefgpool[luckygr].dav0;
					dgr.dav = ens->worlddefgpool[luckygr].dav;

					mydefgpool.push_back ( dgr );

					//finally assign local copy of this grain to the pool
					tmpdefgseeds[sd].mydefgpoolid = mydefgpool.size() - 1;
					sd++;
//cout << "F-sd;luckyGrain;DoIKnowLuckyGrain;mydefgpoolid;" << sd << ";" << luckygr << ";" << DoIKnowLuckyGrain << ";" << mydefgpool[(mydefgpool.size() - 1)].worlddefgid << endl;
					continue;
				}

				//obviously, DoIKnowLuckyGrain == true -> I picked the grain already some loops before
				for ( uint32_t testg = 0; testg < mydefgpool.size(); testg++) {
					if ( luckygr == mydefgpool[testg].worlddefgid ) {
						tmpdefgseeds[sd].mydefgpoolid = testg;
						sd++;
						break;
//cout << "T-sd;luckyGrain;DoIKnowLuckyGrain;testg;" << sd << ";" << luckygr << ";" << DoIKnowLuckyGrain << ";" << testg << endl;
					}
				}
			} //next seed please
		}
	}

	QUICKASSERT ( sd == ngrxyz );

cout << "DAMASK Tmpdefgseeds list generated." << endl;

	//voxelize, packing strategy: align positive x stack x lines positive in y and stack xy slices in positive z direction
	uint32_t mycanx = myCAGeometry.nboxedge_rd;
	uint32_t mycany = myCAGeometry.nboxedge_nd;
	uint32_t mycanz = myCAGeometry.nboxedge_td;
	double cs = myCAGeometry.cellsize;
	double cs2 = 0.5 * cs;

	uint32_t cxyz = 0;
	uint32_t gz, gz_ngrxy, gy, gy_ngrx, gx, gxyz_sd;

	for ( uint32_t iz = 0; iz < mycanz; ++iz) {
		gz = ( ((iz * cs) + cs2 - tz) / CPFEM_GRAINSIZEZ );
		gz_ngrxy = gz * ngrxy;

		for ( uint32_t iy = 0; iy < mycany; ++iy) {
			gy = ( ((iy * cs) + cs2 - ty) / CPFEM_GRAINSIZEY );
			gy_ngrx = gy * ngrx;

			for ( uint32_t ix = 0; ix < mycanx; ++ix) {
				gx = ( ((ix * cs) + cs2 - tx) / CPFEM_GRAINSIZEX );
				gxyz_sd = gx + gy_ngrx + gz_ngrxy;

				//###DEBUG ONLY
				QUICKASSERT ( gxyz_sd < ndefgseeds); 

				mycellgrid[cxyz] = gxyz_sd; //MK::the grain boundary detection algorithm relies on the fact that the deformed grains are disjoint even though the same in properties;

				cxyz++;
			}
		} //stack xlines in y
	} //stack xy slices in z

	if (myensHdl->myRank == MASTER) {
		cout << "myRank " << this->myensHdl->myRank << "; JobID " << this->jobid << " voxelization damask was successful." << endl;
	}
}


//SYNTHETIZATION WITH POINT PROCESS + DISTANCE TRANSFORMATION
void caHdl::solve_SYNTHETIZED_DEFORMEDSTRUCTURE_AUTOMATON( void )
{
	//allocate memory for the voxelized deformed structure = the automaton
	mycellgrid = NULL;
	mycellgrid = new uint32_t[myCAGeometry.nboxvol_rdndtd];
	QUICKASSERT ( mycellgrid != NULL );
	myMemGuard = myMemGuard + (myCAGeometry.nboxvol_rdndtd * sizeof(uint32_t));
	//now memory state is arbitrary but set explicitly in the following functions to a defined vacuum state against which the CA can test

	//MK::the assignment in mycellgrid refers to tmpdefgseeds not of mydefgpool IDs!
	//MK::the synthetization includes definition of seed points + tessellation of the structure by distance transformation with a cellular automaton
	//there is a functionality with which to track explicitly the boundaries among impinging boundaries, though this was shown to be less accurate than to reidentify the boundaries afterwards

	if ( myCADefMS.defmstype == POISSONVORONOI_DEFMS ) {
		//strictly speaking only here we have to initialize mycellgrid[c] to NOT_ASSIGNED_YET
		uint32_t nxyz = myCAGeometry.nboxvol_rdndtd;
		for ( uint32_t c = 0; c < nxyz; ++c ) {
			mycellgrid[c] = NOT_ASSIGNED_YET;
		}

		grow_deformedgrains_voxelize_poisson();

		//###place particles
		return;
	}

	if ( myCADefMS.defmstype == CPFEMDAMASK_DEFMS ) {
		grow_deformedgrains_voxelize_damask();

		//###place particles
		return;
	}
}


//SYNTHETIZATION CLASSICALLY PACKING GRAINS DIRECTLY WITH MYCELLGRID IDS
void caHdl::pick_deformedgrains( void )
{
	//##MK::everything is possible, GIA-aggregate packing, MODF optimization or simply random picking
	//##MK::here we start with sampling by picking randomly from the worlddefgpool grains resulting in an uncorrelated MODF
	ensembleHdlP ens = this->myensHdl;
	uint32_t nensdgpool = ens->worlddefgpool.size();
	uint32_t luckyGrain;

//cout << "ndefgseeds=" << ndefgseeds << endl;

	for ( uint32_t dg = 0; dg < ndefgseeds; dg++ ) {
		luckyGrain = localprng.leEcuyer() * nensdgpool; //lucky grain is an ID in ensembleHdl defgpool array thus the orientation can be obtained by dereferencing the "ens->defgpool[luckyGrain].ori"-th entry of ens->oripool...

		//ASSIGNMENT IS WITH WORLDDEFGPOOL-ID, mind this is in contrast to mydefggrid where IDs refer to this->mydefgpool !
		tmpdefgseeds[dg].ensdefgpoolid = luckyGrain;

//cout << "1-dg;luckyGrain\t\t" << dg << ";" << luckyGrain << ";" << endl;

		//append only to mydefgpool if we have not picked this grain already, therefore scan known ones
		bool DoIKnowLuckyGrain = false;
		for ( uint32_t tg = 0; tg < dg; tg++) {
			if ( luckyGrain == tmpdefgseeds[tg].ensdefgpoolid ) {
				DoIKnowLuckyGrain = true;
				break;
			}
		}

//cout << "2-dg;luckyGrain;DoIKnowLuckyGrain\t\t" << dg << ";" << luckyGrain << ";" << DoIKnowLuckyGrain << endl;

		if ( DoIKnowLuckyGrain == false ) { //most likely this condition is met when ensHdl->defgpool is large
			struct cadefg dgr;

			//get orientation copied from ensHdl
			uint32_t ensoriid = ens->worlddefgpool[luckyGrain].ori;
			double ensbunge[3];
			ensbunge[0] = ens->worldoripool[ensoriid].bunge1;
			ensbunge[1] = ens->worldoripool[ensoriid].bunge2;
			ensbunge[2] = ens->worldoripool[ensoriid].bunge3;

			//orientation recategorization with known Bunge values in local caHdl to utilize that not all possible orientations are used in the local CA grid thus a more efficient hashtable can be generated however at additional memory costs...
			dgr.worlddefgid = luckyGrain;
			dgr.caori = ca_check_disjunctness( ensbunge );
			dgr.cellcount = 0;
			dgr.rho0 = ens->worlddefgpool[luckyGrain].rho0;

			if (dgr.rho0 >= this->myrhomax) { //identify highest dislocation density in the system for adaptive integration scheme
				myrhomax = dgr.rho0;
			}
			dgr.rho = ens->worlddefgpool[luckyGrain].rho;
			dgr.dav0 = ens->worlddefgpool[luckyGrain].dav0;
			dgr.dav = ens->worlddefgpool[luckyGrain].dav;
			dgr.avdg0 = ens->worlddefgpool[luckyGrain].avdg0;
			dgr.avdg = ens->worlddefgpool[luckyGrain].avdg;

			mydefgpool.push_back ( dgr );

			//as I dont know this grain assign the copy of the grain local to caHdl to the pool
			tmpdefgseeds[dg].mydefgpoolid = (mydefgpool.size() - 1);

//cout << "3-dg;luckyGrain;UNKNOWN" << dg << ";" << luckyGrain << ";" << DoIKnowLuckyGrain << ";" << tmpdefgseeds[dg].ensdefgpoolid << "__" << tmpdefgseeds[dg].mydefgpoolid << endl;
			continue;
		}


		//obviously, DoIKnowLuckyGrain == true -> I picked the grain already during the loop
		//well then find it in mylist utilizing that worlddefgid are >= 0 and unique
		for ( uint32_t testg = 0; testg < mydefgpool.size(); ++testg) {
			if ( luckyGrain == mydefgpool[testg].worlddefgid ) {
				tmpdefgseeds[dg].mydefgpoolid = testg;
//cout << "4-dg;luckyGrain;KNOWN" << dg << ";" << luckyGrain << ";" << DoIKnowLuckyGrain << ";" << "testgrain=" << testg << "\t\t" << tmpdefgseeds[dg].ensdefgpoolid << "__" << tmpdefgseeds[dg].mydefgpoolid << endl;
				break;
			}
		}
	}

	//##calculate local MODF if desired here
#ifdef REPORTSTYLE_DEVELOPER
	cout << this->jobid << "\t\t" << this->myensHdl->myRank << "\t\t" << myMemGuard << "\t\t" << mydefgpool.size() << "\t\t" << "\t\tjobid;myensHdlmyRank;myMemGuard (Byte);mydefgpool.size() successfully packed GIA grains." << endl;
#endif
}


void caHdl::voxelize_defms( void )
{
	//MK::a memory marking intensive function + number crunching
	//in which deformed grain the center of the ix,iy,iz -th comes to rest
	//packing strategy: align positive x stack x lines positive in y and stack xy slices in positive z direction
	uint32_t mycanz = myCAGeometry.nboxedge_td;
	uint32_t mycany = myCAGeometry.nboxedge_nd;
	uint32_t mycanx = myCAGeometry.nboxedge_rd;
	uint32_t ngrx = myCADefMS.ngrx;
	uint32_t ngry = myCADefMS.ngry;
	uint32_t ngrz = myCADefMS.ngrz;
	uint32_t ngrxy = myCADefMS.ngrxy;
	uint32_t ngrxyz = myCADefMS.ngrxyz;

	double cs = myCAGeometry.cellsize;
	double cs2 = 0.5 * cs;
	double ux = myCADefMS.u_xrd;	//already absolute displacement of both grids!
	double uy = myCADefMS.v_ynd;
	double uz = myCADefMS.w_ztd;

	double dgrd = myPhysData.defgmean_rd;
	double dgnd = myPhysData.defgmean_nd;
	double dgtd = myPhysData.defgmean_td;


	//now everything is between the stack and two linear arrays...
	uint32_t cxyz = 0;
	uint32_t gz, gz_ngrxy, gy, gy_ngrx, gx, gxyz;

//cout << "ngrxyz\t\t" << ngrx << ";" << ngry << ";" << ngrz << ";" << ngrxy << ";" << ngrxyz << endl;

	for (uint32_t iz = 0; iz < mycanz; ++iz) {
		gz = (((iz * cs) + cs2 - uz) / dgtd);
		gz_ngrxy = gz * ngrxy;

		for (uint32_t iy = 0; iy < mycany; ++iy) {
			gy = (((iy * cs) + cs2 - uy) / dgnd);
			gy_ngrx = gy * ngrx;

			for (uint32_t ix = 0; ix < mycanx; ++ix) {
				gx = (((ix * cs) + cs2 - ux) / dgrd);
				gxyz = gx + gy_ngrx + gz_ngrxy;

				QUICKASSERT ( gxyz < ndefgseeds);
				//MKlocation of the tmpdefgseeds is implicit in sampling the cuboid aggregate

				mycellgrid[cxyz] = gxyz;

//cout << "iz;iy;ix;gz;gy;gx;gz_ngrxy;gy_ngrx;gxyz;tmpdefgsseeds;cxyzmycellgridcxyz\t\t" << iz << ";" << iy << ";" << ix << ";" << gz << ";" << gy << ";" << gx << ";" << gz_ngrxy << ";" << gy_ngrx << ";" << gxyz << ";" << tmpdefgseeds[gxyz].mydefgpoolid << ";" << cxyz << ";" << mycellgrid[cxyz] << endl;

				cxyz++;
			}
		} //stack xlines
	} //stack xy slices

	QUICKASSERT ( cxyz == (mycanx * mycany * mycanz) );

	if (myensHdl->myRank == MASTER) {
		cout << "myRank " << this->myensHdl->myRank << "; JobID " << this->jobid << " voxelization was successful." << endl;
	}
}


void caHdl::determine_initial_deformationtexture_mydefgpoolids( void )
{
	//traverse implicit 3D array in the cachelines
	uint32_t mycanxyz = myCAGeometry.nboxvol_rdndtd;

	StoredEnergy = 0.0;

	//cache deformation structure
	uint32_t* defmicrotexture = NULL;
	uint32_t ndefmicrotexture = myensHdl->standardlagen.size() + 1;
	defmicrotexture = new uint32_t[ndefmicrotexture];
	QUICKASSERT ( defmicrotexture != NULL );
	for (uint32_t i = 0; i < ndefmicrotexture; i++) { 
		defmicrotexture[i] = 0; 
	}

	//scan entire automaton
	uint32_t cadefgid;
	uint32_t ocategory;
	double scaler = (1.0 / SCALING_DEFAULT_DISDENS);
	for ( uint32_t cxyz = 0; cxyz < mycanxyz; ++cxyz ) {
		cadefgid = mycellgrid[cxyz];

		mydefgpool[cadefgid].cellcount += 1;

		StoredEnergy = StoredEnergy + (mydefgpool[cadefgid].rho0 * scaler);

//cout << cxyz << ";" << StoredEnergy << ";" << myoripool[mydefgpool[cadefgid].caori].bunge1 << ";" <<  myoripool[mydefgpool[cadefgid].caori].bunge2 << ";" <<  myoripool[mydefgpool[cadefgid].caori].bunge3 << ";" << mydefgpool[cadefgid].rho0 << endl;

		ocategory = myoripool[mydefgpool[cadefgid].caori].closestideal; //because random is 0
		defmicrotexture[ocategory] = defmicrotexture[ocategory] + 1;

//cout << "cxyz;cadefgid;cellcnt\t\t" << cxyz << "\t\t" << cadefgid << "\t\t" << mydefgpool[cadefgid].cellcount << endl;
	}

	//carry over this characterization of the deformation microstructure
	DefMicrotextureClasses = ndefmicrotexture;
	DefMicrotexture = defmicrotexture;

	if (myensHdl->myRank == MASTER ) { 
		cout << "myRank " << this->myensHdl->myRank << "; JobID " << this->jobid << " initial texture determined." << endl;
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

	NucleiMicrotextureClasses = nnucmicrotex;
	NucleiMicrotexture = nucmicrotex;
}


void caHdl::solve_SYNTHETIZED_DEFORMEDSTRUCTURE_CUBOIDBLOCKS( void )
{
	//determine relative displacement with localprng
	double u = -1.0 * localprng.leEcuyer(); //relative to the grain size
	double v = -1.0 * localprng.leEcuyer();
	double w = -1.0 * localprng.leEcuyer();

	//therewith, determine required size of the 3D grain grid
	uint32_t nx = ( (myCAGeometry.boxedge_rd - (u * myPhysData.defgmean_rd)) / myPhysData.defgmean_rd ); //will be zero if defg+u larger than box, so always add one grain
	nx++; 
	uint32_t ny = ( (myCAGeometry.boxedge_nd - (v * myPhysData.defgmean_nd)) / myPhysData.defgmean_nd );
	ny++;
	uint32_t nz = ( (myCAGeometry.boxedge_td - (w * myPhysData.defgmean_td)) / myPhysData.defgmean_td );
	nz++;

//cout << "MyRank;JobID;uvw;nxyz\t\t" << myensHdl->myRank << ";" << this->jobid << "\t\t" << u << ";" << v << ";" << w << "\t\t" << nx << ";" << ny << ";" << nz << endl;

	myCADefMS.ngrx = nx;
	myCADefMS.ngry = ny;
	myCADefMS.ngrz = nz;
	myCADefMS.ngrxy = nx * ny;
	myCADefMS.ngrxyz = nx * ny * nz;

	ndefgseeds = myCADefMS.ngrxyz;
	QUICKASSERT ( ndefgseeds < CA_ALLOCATION_MAXIMUM );

	myCADefMS.u_xrd = u * myPhysData.defgmean_rd;
	myCADefMS.v_ynd = v * myPhysData.defgmean_nd;
	myCADefMS.w_ztd = w * myPhysData.defgmean_td;

	//MK::allocate memory for the topology of the deformed cuboid-shaped grains, as grains from the input can be used multiple times but
	//the boundary detection requires grains to be of disjoint IDs it is absolutely necessary to keep a list of the seeds to which any grain if all the same can be linked
	tmpdefgseeds = NULL;
	tmpdefgseeds = new struct defgseed[ndefgseeds];
	QUICKASSERT ( tmpdefgseeds != NULL );
	myMemGuard = myMemGuard + (ndefgseeds * sizeof(struct defgseed));
	for ( uint32_t sd = 0; sd < ndefgseeds; ++sd ) {
		tmpdefgseeds[sd].ensdefgpoolid = NOT_ASSIGNED_YET;
		tmpdefgseeds[sd].mydefgpoolid = I_DONT_KNOW_YET;
		tmpdefgseeds[sd].location = NOT_ASSIGNED_YET;
	}


	pick_deformedgrains();


	//allocate memory for the voxelized deformed structure
	mycellgrid = NULL;
	mycellgrid = new uint32_t[myCAGeometry.nboxvol_rdndtd];
	QUICKASSERT ( mycellgrid != NULL );
	myMemGuard = myMemGuard + ( myCAGeometry.nboxvol_rdndtd * sizeof(uint32_t) );


	voxelize_defms();

	//##MK::add, if desired the explicit placement of second phase constituent particles

	//##MK::for Indranil Basu, assign some of the deformed grains here to form contiguous regions of macroscopic shear bands
}


void caHdl::solve_SYNTHETIZE_DEFORMATIONSTRUCTURE( void )
{
	if ( myCADefMS.defmstype == CUBOID_DEFMS ) {
		solve_SYNTHETIZED_DEFORMEDSTRUCTURE_CUBOIDBLOCKS();
	}
	else if ( myCADefMS.defmstype == POISSONVORONOI_DEFMS || myCADefMS.defmstype == CPFEMDAMASK_DEFMS ) {
		solve_SYNTHETIZED_DEFORMEDSTRUCTURE_AUTOMATON();
	}
	else {
		cout << "ERROR::Unable to detect the deformation microstructure synthesis method!" << endl; 
	}

	nmydefgpool = mydefgpool.size();
}


void caHdl::solve_DETECT_GRAINBOUNDARIES( void ) 
{
	initializeBoundariesFast();

	detectGrainBoundaryCellsFast();

	trimGBMemoryAndCalculateDisoriFast();

	if ( this->outopt_localrenderboundaries == RENDERING_BOUNDARIES_YES && renderingForThisCA == true ) {
		visualizeGBTopologyFast();
	}
}


//////////////NUCLEATION MODELING at TRACKED DEFORMED GRAIN BOUNDARIES UTILIZING AUTOMATON WITH SEEDIDS
void caHdl::solve_GRAINBOUNDARYNUCLEATION( void )
{
	//##MK in M. Kuehbach's Phd thesis the specific implementation of grain boundary nucleated reactions that
	//was exercised with C. Haase's HighMn SFB761 Twip steel project aided to model site-saturated grain boundary nucleated RX processes
	//so it sufficies to populate the myrxgpool vector with picking physically the desired per se disjoint boundary cells,
	//at which nuclei should be generated, then solve_NUCLEATIONMODELING must not be executed 
	//instead solve_RXGROWTH places the nuclei without prior knowledge that the mycellgrid assignment came from a DEFMS automaton

	solve_nucmodeling_gbnuc_physicalmodelFast();

	emptyBoundariesFast();
}


void caHdl::solve_nucmodeling_csrenforced( void )
{
	//MK::place exactly myNucleationModel.defaultnucdensity in each automaton working on implicit 3D coordinates, THIS IS THE OLD SCORE NUCLEATION MODEL
	uint32_t nuc = 0;
	uint32_t maxattempts = 0;
	uint32_t ncellstotal = myCAGeometry.nboxvol_rdndtd;
	ensembleHdlP ens = this->myensHdl;
	uint32_t nworldrxgpool = ens->worldrxgpool.size();

	do {
		//pick a side at random
		uint32_t luckySite = localprng.leEcuyer() * ncellstotal;

		bool AlreadyOccupied = false;
		for ( uint32_t n = 0; n < nuc; n++ ) {
			if ( luckySite == myrxgpool[n].nucsite ) {
				AlreadyOccupied == true;
				break;
			}
		}

		if ( AlreadyOccupied == false ) { //obviously juvenile place
			//###currently, pick nucleus at random from worldrxgpool 
			uint32_t luckyNucleus = localprng.leEcuyer() * nworldrxgpool;

			//contrary to the deformed grains in mydefgpool, each nucleus requires entry in myrxgpool as then incubation time can be different and furthermore also the site and single-nucleus resolved kinetics should be stored
			struct carxg crxgr;

			//get orientation
			uint32_t ensoriid = ens->worldrxgpool[luckyNucleus].ori;

			double ensbunge[3];
			ensbunge[0] = ens->worldoripool[ensoriid].bunge1;
			ensbunge[1] = ens->worldoripool[ensoriid].bunge2;
			ensbunge[2] = ens->worldoripool[ensoriid].bunge3;

			//orientation recategorization utilizes that not all possible orientations are used in the local CA grid thus a more efficient hashtable can be generated however at additional memory costs...
			crxgr.caori = ca_check_disjunctness( ensbunge );
			crxgr.cellcount = 0;
			crxgr.nucsite = luckySite; //(200^3 center = 4020100); //(100^3 center=505050); //mluckySite;
			crxgr.startingsite = luckySite;
			crxgr.tincub = ens->worldrxgpool[ensoriid].tincub;

			myrxgpool.push_back ( crxgr );

#ifdef REPORTSTYLE_DEVELOPER
			cout << "jobid;nuc;att;nucsite;luckyNucleusWorldRXGPool\t\t" << this->jobid << ";" << nuc << ";" << maxattempts << ";" << myrxgpool[myrxgpool.size()-1].startingsite << ";" << luckyNucleus << endl;
#endif
			nuc++;
			continue;
		}

		//damn, collision - AlreadyOccupied == true 
		maxattempts++;
	} 
	while (nuc < myNucleationModel.defaultnucdensity && maxattempts < MAXATTEMPTS_NUCSITE );

	//characterize the probabilities in the CA that nuclei are able to grow in their parent deformed grain in the surrounding x,y,z minding periodic boundary conditions

	//characterize_nucsite_environment(); //##MK::no HPF for C.Haase

#ifdef DETAILED_PROMPTS
	if ( myensHdl->myRank == MASTER ) { cout << this->jobid << "\t\tNucleation completed." << endl; }
#endif
}


void caHdl::solve_nucmodeling_csrpickrandomly( void )
{
	//in a large CSR point pattern with a unit density what is the size of the observation window in generalized coordinates?
	//with myNucleationModel.defaultnucdensity the limit density for infinite draws of owins
	//the difference compared to enforce is that when an obserrvation volume is small defaultnucdensity is only the expectation value
	//but never always the total amount of points one may find in the box
	double scale = pow( (((double) myNucleationModel.defaultnucdensity) / ((double) MASTER_CSR_NPOINTS)), (1.0/3.0) );

	double owin[3];
	owin[0] = scale;	if ( owin[0] > 1.0 ) { cout << "ERROR::MASTER point process is too small!" << endl; return; }
	owin[1] = scale;	if ( owin[1] > 1.0 ) { cout << "ERROR::MASTER point process is too small!" << endl; return; }
	owin[2] = scale;	if ( owin[2] > 1.0 ) { cout << "ERROR::MASTER point process is too small!" << endl; return; }


	//populate a CSR point pattern on a unit cube from which to take a sample
	pointP pp3 = NULL;
	pp3 = new point[MASTER_CSR_NPOINTS];
	QUICKASSERT ( pp3 != NULL );

	for ( uint32_t p = 0; p < MASTER_CSR_NPOINTS; ++p ) {
		pp3[p].x = localprng.leEcuyer(); //[0,1]^3
		pp3[p].y = localprng.leEcuyer();
		pp3[p].z = localprng.leEcuyer();
	}

	//sample an owin from pp3 on the unit cube
	double owinoo[6];
	bool inside = false;

	while (inside == false ) {
		owinoo[XMI] = localprng.leEcuyer();
		owinoo[XMX] = owinoo[XMI] + owin[0];
		owinoo[YMI] = localprng.leEcuyer();
		owinoo[YMX] = owinoo[YMI] + owin[1];
		owinoo[ZMI] = localprng.leEcuyer();
		owinoo[ZMX] = owinoo[ZMI] + owin[2];

		if ( owinoo[XMX] < 1.0 && owinoo[YMX] < 1.0 && owinoo[ZMX] < 1.0) //all coordinates positive owinoo cannot protrude in negative direction
			inside = true;
	}

	//extract all points from this window and sample via (i - owinoo[iMI])/owin[i] * ni into an automaton coordinate
	double dx, dy, dz;
	uint32_t ix, iy, iz, luckysite;
	uint32_t nbx = myCAGeometry.nboxedge_rd; //avoid the cache collision when calling this in the for loop to select the points
	uint32_t nby = myCAGeometry.nboxedge_nd;
	uint32_t nbz = myCAGeometry.nboxedge_td;
	bool valid = false;

	for ( uint32_t p = 0; p < MASTER_CSR_NPOINTS; ++p ) {
		if ( pp3[p].x < owinoo[XMI] ) continue;
		if ( pp3[p].x > owinoo[XMX] ) continue;
		if ( pp3[p].y < owinoo[YMI] ) continue;
		if ( pp3[p].y > owinoo[YMX] ) continue;
		if ( pp3[p].z < owinoo[ZMI] ) continue;
		if ( pp3[p].z > owinoo[ZMX] ) continue;

		//p is inside owinoo so, discretize location, dx, dy, dz > 0 < 1
		dx = (pp3[p].x - owinoo[XMI]) / (owinoo[XMX] - owinoo[XMI]);
		dy = (pp3[p].y - owinoo[YMI]) / (owinoo[YMX] - owinoo[YMI]);
		dz = (pp3[p].z - owinoo[ZMI]) / (owinoo[ZMX] - owinoo[ZMI]);
		ix = dx * nbx;
		iy = dy * nby;
		iz = dz * nbz;
		luckysite = ix + (iy * nbx) + (iz * nbx * nby);

//cout << "csrpickrnd;luckysite = " << luckysite << "\t\t" << ix << ";" << iy << ";" << iz << endl;

		//disprove the assumption that the hole contains no seed, overlap of two seeds within <= gridresu not allowed
		valid = true;
		for ( uint32_t s = 0; s < myrxgpool.size(); s++ ) {
			if ( luckysite == myrxgpool[s].nucsite ) { valid = false; break; }
		}

		//still true, okay hole is empty, plant a seed
		if ( valid == true ) {

			//assign orientation and a nucleus
			uint32_t luckyNucleus = localprng.leEcuyer() * myensHdl->worldrxgpool.size();

			//contrary to the deformed grains in mydefgpool, each nucleus requires entry in myrxgpool as then incubation time can be different and furthermore also the site and single-nucleus resolved kinetics should be stored
			struct carxg crxgr;

			//get orientation
			uint32_t ensoriid = myensHdl->worldrxgpool[luckyNucleus].ori;

			double ensbunge[3];
			ensbunge[0] = myensHdl->worldoripool[ensoriid].bunge1;
			ensbunge[1] = myensHdl->worldoripool[ensoriid].bunge2;
			ensbunge[2] = myensHdl->worldoripool[ensoriid].bunge3;

			//orientation recategorization utilizes that not all possible orientations are used in the local CA grid thus a more efficient hashtable can be generated however at additional memory costs...
			crxgr.caori = ca_check_disjunctness( ensbunge );
			crxgr.cellcount = 0;
			crxgr.nucsite = luckysite;
			crxgr.startingsite = luckysite;
			crxgr.tincub = myensHdl->worldrxgpool[ensoriid].tincub;

			myrxgpool.push_back ( crxgr );
		} //next point
	}

	delete [] pp3;
	pp3 = NULL;
}


void caHdl::init_myMaternUniverse( void )
{
	//a very large Matern cluster process on a cubic RVE of size rvesize in micron is squeezed on a unit cube in \mathbb{R}^3
	//point process is on cubic RVE with myNucleationModel.cluster_rvesize^3 micron each cluster radius is cluster_extenda * cellsize
	double scalingradius = myNucleationModel.cluster_a * myCAGeometry.cellsize / myNucleationModel.cluster_rvesize;
	QUICKASSERT ( scalingradius < DEFAULT_MASTERMATERN_MAXRADIUS );

	mypointprocess.reserve ( (long) (myNucleationModel.cluster_nclust * myNucleationModel.cluster_lambda) );

cout << this->jobid << " -- GenerateMaternClusterUniverse NClust/MeanDensityPerClust/ScalingRVESize/RadiusClust/ExtendABC = " << myNucleationModel.cluster_nclust << ";" << myNucleationModel.cluster_lambda << ";" << myNucleationModel.cluster_rvesize << ";" << scalingradius << ";" << myNucleationModel.cluster_a << ";" << myNucleationModel.cluster_b << ";" << myNucleationModel.cluster_c << endl;

	double nx, ny, nz, len;
	//##MK::distortion currently not implemented
	//##MK::assumption is that number of clusters is large (at least 1000!) to sample random 3D poisson point process

	double cx, cy, cz;
	long npc;
	for (unsigned long c = 0; c < myNucleationModel.cluster_nclust; c++) {
		//construction of Matern: define first a cluster center that itself does not count as a point
		cx = localprng.leEcuyer();
		cy = localprng.leEcuyer();
		cz = localprng.leEcuyer();

		//how many points to locate in one cluster, evaluate poisson distribution, naive cumulant method
		long k = 0;
		double lambda = myNucleationModel.cluster_lambda;
		double unifrnd = localprng.leEcuyer();
		double P = exp(-lambda);
		double sum = P;

		npc = 0; //handles case sum >= unifrnd
		if ( sum < unifrnd ) {
			for ( k = 1; k < DEFAULT_KMAX_POISSRND; k++) {
				P *= lambda / (double) k;
				sum += P;
				if ( sum >= unifrnd ) break;
			}
		}
		npc = k;


		//populate the process
		for (unsigned long p = 0; p < npc; p++) {
			//random rotor instead of predictor/branch/reject
			nx = 2*localprng.leEcuyer() - 1.0; //(-1.0,1.0) a random 3D vector pointing into space...
			ny = 2*localprng.leEcuyer() - 1.0;
			nz = 2*localprng.leEcuyer() - 1.0;
			len = ( 1.0 / pow( ( SQR(nx)+SQR(ny)+SQR(nz) ) , 0.5) ) * (localprng.leEcuyer() * scalingradius); //normalized component value * scale(1.0 / ... ) normalizes define scaling of random rotor inside a sphere of maximum radius mmatr
			//scale direction vector
			nx *= len;
			ny *= len;
			nz *= len;
			//mount scaled vector pointing in nxyz direction vector into the cluster center
			nx += cx;
			ny += cy;
			nz += cz;

			//apply PBC in the universe
			if ( nx < 0.0 ) nx += 1.0;
			if ( nx > 1.0 ) nx -= 1.0;
			if ( ny < 0.0 ) ny += 1.0;
			if ( ny > 1.0 ) ny -= 1.0;
			if ( nz < 0.0 ) nz += 1.0;
			if ( nz > 1.0 ) nz -= 1.0;

			struct point ap;
			ap.x = nx;
			ap.y = ny;
			ap.z = nz;

//cout << "\t\t" << c << "\t\t" << npc << "\t\t" << p << "\t\t" << lambda << ";" << unifrnd << ";" << P << ";" << sum << ";" << cx << ";" << cy << ";" << cz << "--" << len << ";" << nx << ";" << ny << ";" << nz << endl;

			mypointprocess.push_back(ap);
		} //populate the cluster
	} //for all clusters
}


void caHdl::pickrandomly_myMaternUniverse( void )
{
	//the solitary unit domain has a specific size (edge length in microns) that is a fraction of myNucleationModel.clust_rvesize
	double owinx = myCAGeometry.boxedge_rd / myNucleationModel.cluster_rvesize;
	double owiny = myCAGeometry.boxedge_nd / myNucleationModel.cluster_rvesize;
	double owinz = myCAGeometry.boxedge_td / myNucleationModel.cluster_rvesize;

	double owinoo[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	bool found = false;

	//MK::to center the domain on a particular point view PPANALYSIS source code...
	//cut myself a window from mypointprocess
	found = false;
	while ( found == false ) {
		owinoo[XMI] = localprng.leEcuyer();
		owinoo[XMX] = owinoo[XMI] + owinx;
		owinoo[YMI] = localprng.leEcuyer();
		owinoo[YMX] = owinoo[YMI] + owiny;
		owinoo[ZMI] = localprng.leEcuyer();
		owinoo[ZMX] = owinoo[ZMI] + owinz;

		//inside?
		if (owinoo[XMI] < 0.0) 	continue;
		if (owinoo[XMX] > 1.0)	continue;
		if (owinoo[YMI] < 0.0) 	continue;
		if (owinoo[YMX] > 1.0)	continue;
		if (owinoo[ZMI] < 0.0) 	continue;
		if (owinoo[ZMX] > 1.0)	continue;

		found = true;
	}

//cout << "Owinxyz = " << owinx << ";" << owiny << ";" << owinz << endl;
cout << "The window = " << owinoo[XMI] << ";" << owinoo[XMX] << ";" << owinoo[YMI] << ";" << owinoo[YMX] << ";" << owinoo[ZMI] << ";" << owinoo[ZMX] << endl;

	QUICKASSERT ( found == true );

	//now test which points from mypointprocess are in this window and plan them for the simulation
	double afftransx = (1.0 / (owinoo[XMX] - owinoo[XMI])) * myCAGeometry.boxedge_rd;
	double afftransy = (1.0 / (owinoo[YMX] - owinoo[YMI])) * myCAGeometry.boxedge_nd;
	double afftransz = (1.0 / (owinoo[ZMX] - owinoo[ZMI])) * myCAGeometry.boxedge_td;

	struct point sp;
	uint32_t nrdnd = myCAGeometry.nboxarea_rdnd;
	uint32_t nrd = myCAGeometry.nboxedge_rd;
	uint32_t nnd = myCAGeometry.nboxedge_nd;
	uint32_t ntd = myCAGeometry.nboxedge_td;
	uint32_t ix, iy, iz;
	double box_rd = myCAGeometry.boxedge_rd;
	double box_nd = myCAGeometry.boxedge_nd;
	double box_td = myCAGeometry.boxedge_td;
	double cellsize = myCAGeometry.cellsize;

	uint32_t nworldrxgpool = this->myensHdl->worldrxgpool.size();

	long n = mypointprocess.size();
	bool occupied = false;

	for (unsigned long p = 0; p < n; p++) { //all nuclei
		if (mypointprocess[p].x < owinoo[XMI]) continue;
		if (mypointprocess[p].x > owinoo[XMX]) continue;
		if (mypointprocess[p].y < owinoo[YMI]) continue;
		if (mypointprocess[p].y > owinoo[YMX]) continue;
		if (mypointprocess[p].z < owinoo[ZMI]) continue;
		if (mypointprocess[p].z > owinoo[ZMX]) continue;

		//location inside the CA domain MK:: mypointprocess[p].i >= owinoo[iMX]
		sp.x = (mypointprocess[p].x - owinoo[XMI]) * afftransx;
		sp.y = (mypointprocess[p].y - owinoo[YMI]) * afftransy;
		sp.z = (mypointprocess[p].z - owinoo[ZMI]) * afftransz;

		//map on discrete CA grid
		ix = sp.x / cellsize;
		iy = sp.y / cellsize;
		iz = sp.z / cellsize;

		uint32_t luckySite = ix + (iy * nrd) + (iz * nrdnd);

#ifdef REPORTSTYLE_DEVELOPER
	cout << "\t\tJobID;p;luckySite;ix;iy;iz\t\t" << this->jobid << "\t\t" << p << "\t\t" << luckySite << ";" << ix << ";" << iy << ";" << iz << endl; // << "--" << sp.x << ";" << sp.y << ";" << sp.z <<  endl;
#endif

		//site still unoccupied?
		occupied = false;
		uint32_t nucSite = INVALID_NUCSITE;

		//disprove expectation that the site is still free...
		for ( uint32_t i = 0; i < myrxgpool.size(); i++ ) {
			if ( luckySite == myrxgpool[i].nucsite ) {
				occupied = true;
				break;
			}
		}

		//placement only if site still unoccupied
		if ( occupied == false ) {
			//###currently, pick nucleus at random from worldrxgpool 
			uint32_t luckyNucleus = localprng.leEcuyer() * nworldrxgpool;

			//contrary to the deformed grains in mydefgpool, each nucleus requires entry in myrxgpool as then incubation time can be different and furthermore also the site and single-nucleus resolved kinetics should be stored
			struct carxg crxgr;

			//get orientation
			uint32_t ensoriid = myensHdl->worldrxgpool[luckyNucleus].ori;

			double ensbunge[3];
			ensbunge[0] = myensHdl->worldoripool[ensoriid].bunge1;
			ensbunge[1] = myensHdl->worldoripool[ensoriid].bunge2;
			ensbunge[2] = myensHdl->worldoripool[ensoriid].bunge3;

			//orientation recategorization utilizes that not all possible orientations are used in the local CA grid thus a more efficient hashtable can be generated however at additional memory costs...
			crxgr.caori = ca_check_disjunctness( ensbunge );
			crxgr.cellcount = 0;
			crxgr.nucsite = luckySite;
			crxgr.startingsite = luckySite;
			crxgr.tincub = myensHdl->worldrxgpool[ensoriid].tincub;

			myrxgpool.push_back ( crxgr );

#ifdef REPORTSTYLE_DEVELOPER
	cout << "\t\tJobID;p;spx;spy;spz;ix;iy;iz;luckySite;nucsite\t\t" << this->jobid << "\t\t" << p << "\t\t" << sp.x << ";" << sp.y << ";" << sp.z << "--" << ix << ";" << iy << ";" << iz << "--" << myrxgpool[myrxgpool.size()-1].nucsite << ";" << myrxgpool[myrxgpool.size()-1].startingsite << endl;
#endif

		} //a point was found and planned for nucleation

	} //test next point

	cout << "Clustering process performed with so many entries = " << myrxgpool.size() << endl;
}


void caHdl::solve_nucmodeling_ellipsoidalcluster_pickrandomly( void )
{
	
cout << "Pick cluster randomly" << endl;

	//generate a large Matern cluster process on the unit cube [0,1]^3
	init_myMaternUniverse();

	//pick a random window from this process and map all points in the automaton domain
	pickrandomly_myMaternUniverse();

//#ifdef REPORTSTYLE_DEVELOPERcout << "\t\tcluster;nuc;jobid;att;nucsite;luckyNucleusWorldRXGPool\t\t" << c << ";" << k << ";" << this->jobid << ";" << ntrials << ";" << myrxgpool[myrxgpool.size()-1].startingsite << ";" << luckyNucleus << endl;#endif

	//characterize the probabilities in the CA that nuclei are able to grow in their parent deformed grain in the surrounding x,y,z minding periodic boundary conditions

	//characterize_nucsite_environment(); //##MK::no HPF for C.Haase

#ifdef DETAILED_PROMPTS
	if ( myensHdl->myRank == MASTER ) { cout << this->jobid << "\t\tMatern type ellipsoid nucleation completed." << endl; }
#endif
}


void caHdl::solve_NUCLEATIONMODELING( void ) 
{
	if ( myNucleationModel.csrnucleation == CSR_ENFORCED ) {
		solve_nucmodeling_csrenforced();
	}

	if ( myNucleationModel.csrnucleation == CSR_PICKRANDOMLY ) {
		solve_nucmodeling_csrpickrandomly();
	}

	if ( myNucleationModel.clustnucleation == CLUSTNUCLEATION_PICKRANDOMLY ) {
		solve_nucmodeling_ellipsoidalcluster_pickrandomly();
	}
}


void caHdl::solve_REPLACE_CA_STRUCTURE( void ) 
{
	//now we know all potential nuclei, so no new orientations are added further
	nmyrxgpool = myrxgpool.size();

	//MK::still the mycellgrid is populated with seedids
	//and for the growth of the RX nuclei the boundaries among the deformed grains are no longer very interesting
	//therefore replace mycellgrid from indices referencing entries in tmpdefgseeds to indices referencing entries in mydefgpool


	uint32_t nxyz = myCAGeometry.nboxvol_rdndtd;
	uint32_t seedid, mydgid;


	for ( uint32_t c = 0; c < nxyz; c++) {
		seedid = mycellgrid[c];
		mydgid = tmpdefgseeds[seedid].mydefgpoolid;
		mycellgrid[c] = mydgid;
	}

	//texture bookkeeping
	determine_initial_deformationtexture_mydefgpoolids();

	//nuclei bookkeeping
	determine_initial_nucleitexture_myrxgpoolids();


	//clean-up tmdefgseeds as it is not longer needed
	delete [] tmpdefgseeds;
	tmpdefgseeds = NULL;
	myMemGuard = myMemGuard - ( ndefgseeds * sizeof(struct defgseed) );
}


void caHdl::init_myRXFront( void )
{
	//allocate initial memory to store ACTIVE cells
	uint32_t ninitialRXFront = (initialRelCellCaching * (double) myCAGeometry.nboxvol_rdndtd);
	QUICKASSERT ( ninitialRXFront < CA_ALLOCATION_MAXIMUM );
	myRXFront = NULL;
	myRXFront = new struct cell[ninitialRXFront];
	QUICKASSERT( myRXFront != NULL );
	myMemGuard = myMemGuard + (ninitialRXFront * sizeof(cell));

	ntotalRXFront = ninitialRXFront;
	nextSlotNeverActiveRXFront = 0;
	nCurrentlyActive = 0;

	//initialize already associated cells
	for ( uint32_t cfr = 0; cfr < ninitialRXFront; ++cfr ) {
		myRXFront[cfr].activity = INACTIVE;
		myRXFront[cfr].infector = 26;
		myRXFront[cfr].ix = -1;
		myRXFront[cfr].iy = -1;
		myRXFront[cfr].iz = -1;
		myRXFront[cfr].rxFrac = NO_INFECTION;
		myRXFront[cfr].P = -1.0;
		myRXFront[cfr].mydefgid = NO_GRAIN_ASSIGNED;
		myRXFront[cfr].myrxgid = NO_GRAIN_ASSIGNED;
	}


	//allocate memory to store RecycledCells
	uint32_t ninitialRecyclingList = (uint32_t) (initialRelCellCaching * (double) myCAGeometry.nboxvol_rdndtd * maxfillperstep * FULLRECYCLING_ATT_FACTOR);
	QUICKASSERT ( ninitialRecyclingList < CA_ALLOCATION_MAXIMUM );
	myRecyclingList = NULL;
	myRecyclingList = new uint32_t[ninitialRecyclingList];
	QUICKASSERT( myRecyclingList != NULL );
	myMemGuard = myMemGuard + ( ninitialRecyclingList * sizeof(uint32_t) );


	ntotalRecyclingList = ninitialRecyclingList;
	nextSlotThatBecomesRecycled = NOTHING_TO_RECYCLE; //not assigned as there is no cell at all yet  that can be recycled 
	firstNotRecycledYet = NOTHING_TO_RECYCLE;

	//myRecyclingList needs no initialization because it is assured to have no gaps and only read [0;nextSlotThatBecomesRecycled)


	//allocate maxfillperstep the initial memory for cell that complete their infection in a timestep
	//this list serves as a guidance to reduce the amount of ifs during the infection phase as the list is for sure fragmented and thus INACTIVE CELLs are  met
	uint32_t ninitialFullRXList = (uint32_t) (initialRelCellCaching * (double) myCAGeometry.nboxvol_rdndtd * maxfillperstep * FULLRXCACHING_ATT_FACTOR);
	QUICKASSERT ( ninitialFullRXList < CA_ALLOCATION_MAXIMUM );
	myFullRXList = NULL;
	myFullRXList = new uint32_t[ninitialFullRXList];
	QUICKASSERT( myFullRXList != NULL );
	myMemGuard = myMemGuard + ( ninitialFullRXList * sizeof(uint32_t) );

	ntotalFullRXList = ninitialFullRXList;
	nextSlotToFullRX = 0; //as above not assigned yet

	//no initialization as well  assured gap-free read-only [0;nextSlotToFullRX)


	//extra container to check boundary contact
	myrxgpoolboundarycontact = NULL;
	myrxgpoolboundarycontact = new bool[myrxgpool.size()];
	QUICKASSERT( myrxgpoolboundarycontact != NULL );
	myMemGuard = myMemGuard + (myrxgpool.size() * sizeof(bool));
	for ( uint32_t rg = 0; rg < myrxgpool.size(); ++rg ) { 
		myrxgpoolboundarycontact[rg] = false; 
	}

	//cout << this->jobid << "\t\tMemGuard (Byte)= " << this->myMemGuard << endl;
}


uint32_t caHdl::get_NextFreeSlotInRXFrontAppendOrRecycle( bool how ) //IS ONLY CALLED FOR RECYCLING FROM OutOfRXFront
{
	//##MK::two extreme strategies possible
	//1::always get more memory (almost no additional scans for places but increasing cache hit during calcGrowthStep owing to accumulation of INACTIVE cells)
	//2::recycle as best as possible (reduced total amount of memory and less scans on likely improved cache performance than strategy 1)
	//optimal condition invokes compact (meaning no cell in the bucket is INACTIVE) that is only to obtain by array management inducing overhead 
	//or utilizing a linked list structures that to iterate over however is orders of magnitude slower for millions of entries and not practical to become parallelized
	
	//##MK::algorithm is able to select most beneficial strategy automatically
	bool strategy = how;
	uint32_t place;

	//if recycling is desired
	if ( strategy == CELLRECYCLE ) {
		if ( nextSlotThatBecomesRecycled < firstNotRecycledYet ) { //##DEBUG&& (firstNotRecycledYet < ntotalRecyclingList) && 
			place = myRecyclingList[nextSlotThatBecomesRecycled];
			nextSlotThatBecomesRecycled++;

#ifdef REPORTSTYLE_CELLCYCLES
	cout << "\t\tRECYCLE\t" << place << endl;
#endif

			return place;
		}

		//hasnt returned yet, so recycling list should be utilized but was already exhausted so change the strategy
		strategy = CELLAPPEND;
	}

	//now strategy is/or was CELLAPPEND
	//if ( strategy == CELLAPPEND ) {
		//myRXFront cached large enough?
		if ( nextSlotNeverActiveRXFront < ntotalRXFront ) {
			place = nextSlotNeverActiveRXFront;
			nextSlotNeverActiveRXFront++; // ##DEBUG = nextSlotNeverActiveRXFront + 1;

#ifdef REPORTSTYLE_CELLCYCLES
	cout << "\t\tAPPEND\t" << place << endl;
#endif

			return place;
		}

		//obviously no cache left and mind that accesses such as myRXFront[nextSlotNeverActiveRXFront] are faulty when >= ntotalRXFront...
		size_t oldSize = ntotalRXFront;
		size_t newSize = oldSize + (transientRelCellRecaching * oldSize); //##DEBUGmyCAGeometry.nboxvol_rdndtd);
		QUICKASSERT( newSize < CA_ALLOCATION_MAXIMUM );

#ifdef REPORTSTYLE_DEVELOPER
	cout << "\t\tALLOCATION for new_myRXFront;oldSize;newSize\t" << ((int) oldSize) << "\t" << ((int) newSize) << endl;
#endif
		cellP new_myRXFront = new struct cell[newSize];
		myMemGuard = myMemGuard + (newSize * sizeof(cell));

		//potentially unsafe to use memcpy because of overlapping issues memcpy( new_myRXFront, myRXFront, oldSize * sizeof(cell) ); 
		//MK::all indices in FullRX and Recycling list stay valid as they are referencing relative to first array position
		for ( uint32_t i = 0; i < oldSize; ++i ) {
			new_myRXFront[i].activity = myRXFront[i].activity;
			new_myRXFront[i].infector = myRXFront[i].infector;
			new_myRXFront[i].ix = myRXFront[i].ix;
			new_myRXFront[i].iy = myRXFront[i].iy;
			new_myRXFront[i].iz = myRXFront[i].iz;
			new_myRXFront[i].rxFrac = myRXFront[i].rxFrac;
			new_myRXFront[i].P = myRXFront[i].P;
			new_myRXFront[i].mydefgid = myRXFront[i].mydefgid;
			new_myRXFront[i].myrxgid = myRXFront[i].myrxgid;
		}

		delete [] myRXFront;
		myMemGuard = myMemGuard - (oldSize * sizeof(cell));
		myRXFront = new_myRXFront;
		ntotalRXFront = newSize;

		//initialize this new memory snippet, access with nc = oldSize is to 0-th element of appended section
		for ( uint32_t nc = oldSize; nc < ntotalRXFront; ++nc) {
			myRXFront[nc].activity = INACTIVE;
			myRXFront[nc].infector = 26;
			myRXFront[nc].ix = -1;
			myRXFront[nc].iy = -1;
			myRXFront[nc].iz = -1;
			myRXFront[nc].rxFrac = NO_INFECTION;
			myRXFront[nc].P = -1.0;
			myRXFront[nc].mydefgid = NO_GRAIN_ASSIGNED;
			myRXFront[nc].myrxgid = NO_GRAIN_ASSIGNED;
		}

		//now we have new space, so append
		place = nextSlotNeverActiveRXFront;
		nextSlotNeverActiveRXFront++; // = nextSlotNeverActiveRXFront + 1;

		#ifdef REPORTSTYLE_CELLCYCLES
			cout << "\t\tNEWREGION\t" << place << endl;
		#endif

		return place;
	//}
}


void caHdl::sim_myCA_infect_JuvenileNucleation( uint32_t rxgpoolid, uint32_t WhichCelltoInfect, double rxfrac0 )
{
	//identify the grid location x,y,z of the infectorcell
	uint32_t z = WhichCelltoInfect / myCAGeometry.nboxarea_rdnd;
	uint32_t rem = WhichCelltoInfect - (myCAGeometry.nboxarea_rdnd * z);
	uint32_t y = rem / myCAGeometry.nboxedge_rd;
	uint32_t x = rem - (myCAGeometry.nboxedge_rd * y);

	//switch or not whether rxgpoolid already has domain wall contact
	if ( x <= 0 || x >= (myCAGeometry.nboxedge_rd - 1) ) { myrxgpoolboundarycontact[rxgpoolid] = true; }
	if ( y <= 0 || y >= (myCAGeometry.nboxedge_nd - 1) ) { myrxgpoolboundarycontact[rxgpoolid] = true; }
	if ( z <= 0 || z >= (myCAGeometry.nboxedge_td - 1) ) { myrxgpoolboundarycontact[rxgpoolid] = true; }

	//get out here, cell is already infected, but in case of site-saturated nucleation no cell can be infected already
	if ( (mycellgrid[WhichCelltoInfect] >= nmydefgpool) || (mycellgrid[WhichCelltoInfect] == CURRENTLY_INFECTED) ) { //|| (mycellgrid[WhichCelltoInfect] == CELL_IS_A_PARTICLE) ) { 
		myrxgpool[rxgpoolid].nucsite = SITE_ALREADY_CONSUMED;
		//.startingsite is must not be changed
		myrxgpool[rxgpoolid].tincub = -1.0;
		return;
	}

	//obviously cell is free, hehe, go for it!
	uint32_t freeplace = get_NextFreeSlotInRXFrontAppendOrRecycle( CELLAPPEND );

	//request initial state of the cell and then mark cell as infected to prevent multiple reinfections, do right here to localize memory access
	uint32_t mydgid = mycellgrid[WhichCelltoInfect];
	mycellgrid[WhichCelltoInfect] = CURRENTLY_INFECTED;
	myRXFront[freeplace].activity = ACTIVE;
	myRXFront[freeplace].infector = 26;
	myRXFront[freeplace].ix = (short) x;
	myRXFront[freeplace].iy = (short) y;
	myRXFront[freeplace].iz = (short) z;
	myRXFront[freeplace].rxFrac = rxfrac0;
	//cell carries consumed orientation along
	myRXFront[freeplace].mydefgid = mydgid;
	myRXFront[freeplace].myrxgid = rxgpoolid;

	myRXFront[freeplace].P = calc_mobilityweight( rxgpoolid, mydgid );

	//myMobilityWeightMax, which factorizes the global adaptive dynamic time stepping
	if ( myRXFront[freeplace].P >= myMobilityWeightMax ) {
		myMobilityWeightMax = myRXFront[freeplace].P;
	}

	//one more nucleus growing in this automaton
	nmynuclei++;

	myrxgpool[rxgpoolid].nucsite = NUCLEUS_ALREADY_PLACED;
	//.startingsite must not be changed
	myrxgpool[rxgpoolid].tincub = this->t;

	//MK::get_NextFreeSlotInRXFront() ASSURES freeplace <= nextSlotNeverActiveRXFront < ntotalRXFront!

#ifdef REPORTSTYLE_CELLCYCLES
	cout << this->jobid << "\t\tjuvenile infection;rxgpoolid;WhichCelltoInfect;x;y;z;rxfrac0;freeplace" << rxgpoolid << ";" << WhichCelltoInfect << ";" << x << ";" << y << ";" << z << ";" << rxfrac0 << ";" << freeplace << endl;
#endif
}


void caHdl::sim_myCA_sitesaturatedNucleation( void )
{
	//this function expects the mycellgrid positions to be marked with mydefgpool ids not seedids!
	uint32_t nucleateWhere;

	for ( uint32_t nuc = 0; nuc < myrxgpool.size(); nuc++) {
		nucleateWhere = myrxgpool[nuc].nucsite;

		sim_myCA_infect_JuvenileNucleation( nuc, nucleateWhere, ALMOST_FULLY_INFECTED );
	}

	if ( myensHdl->myRank == MASTER ) { 
		cout << "myRank " << myensHdl->myRank << "; JobID " << this->jobid << " nucleated site-saturated attempting to place " << myrxgpool.size() << " nuclei." << endl;
	}
}



inline double caHdl::pdf_rayleigh( double realtime )
{
	//return the nucleation rate at time realtime
	return ( realtime / SQR(myNucleationModel.tincub_rayleigh_sigma) * exp( -0.5 * SQR(realtime)/SQR(myNucleationModel.tincub_rayleigh_sigma) ) );
}


inline double caHdl::cdf_rayleigh( double realtime )
{
	//return the expected number of nuclei at that time
	return ( 1.0 - exp ( -0.5 * SQR(realtime)/SQR(myNucleationModel.tincub_rayleigh_sigma) ) );
}



void caHdl::sim_myCA_timedependentNucleation( void )
{
	if ( nucleationmodel_status == NUCSITES_ALL_EXHAUSTED ) {

#ifdef REPORTSTYLE_TDEPNUCLEATION
	cout << "All sites exhausted!" << endl;
#endif

		return;
	}


	//##MK::to begin we investigate a formal rate model according to a Rayleigh distribution formal rate theory CSR nucleation from any of the nuclei
	//then it is at least necessary to get the current time and dt and how many nuclei are expected
	uint32_t nrxg = myrxgpool.size();
	double dnrxg = nrxg; //##MK::for tolerance estimate

	double nExpectedNumberOfNuclei = dnrxg * cdf_rayleigh( this->t ); //could be 25 * 0.9999 so 24.999 which however is cast down indicating that never all nuclei were placed
	double nCurrentNumberOfNuclei = nmynuclei;

	//how many have to be become placed?
	double nDifference = nExpectedNumberOfNuclei - nCurrentNumberOfNuclei;
	uint32_t nNucleiToPlaceNow = nDifference;

	//MK::principal problem myrxgpool might have many entries most of which in an uncorrelated and fragmented way however have already been placed
	//now how to find efficiently the still available sites without preferences for a particular position?

#ifdef REPORTSTYLE_TDEPNUCLEATION
	cout << "TimeDependentNucleationModelling-deltaN=" << nExpectedNumberOfNuclei << "\t\t" << nCurrentNumberOfNuclei << "\t\t" << nDifference << endl;
#endif

	if ( nNucleiToPlaceNow > 0 ) {
		//create structure over which to identify shorter in caches the still available places
		uint32_t nextFreeSiteReference = 0;
		uint32_t* nucsite_avail = NULL;
		nucsite_avail = new uint32_t[nrxg];
		if ( nucsite_avail == NULL ) { cout << "Memory allocation in time dependent nucleation!" << endl; return; }
		for ( uint32_t n = 0; n < nrxg; n++ ) {
			if ( myrxgpool[n].nucsite != NUCLEUS_ALREADY_PLACED && myrxgpool[n].nucsite != SITE_ALREADY_CONSUMED ) {
				nucsite_avail[nextFreeSiteReference] = n;
				nextFreeSiteReference++;
			}
		}
		//as of now at least nrxg ifs when nucleation is not yet complete but now searching for free places on much shorter cache-linear array...

		//havent found any place at all then in all future iterations we dont have to run through this function
		if ( nextFreeSiteReference == 0 ) {
			nucleationmodel_status = NUCSITES_ALL_EXHAUSTED;
		}

#ifdef REPORTSTYLE_TDEPNUCLEATION
	cout << "TimeDependentNucleationModelling-freesites=" << nextFreeSiteReference << "\t\t" << nNucleiToPlaceNow << endl;
#endif

		if ( nNucleiToPlaceNow > nextFreeSiteReference ) //do not try to seek more places then available
			nNucleiToPlaceNow = nextFreeSiteReference;

		//MK::uncorrelated picking of nuclei from still available sites by picking from nucsite_avail references to myrxgpool
		uint32_t luckysite;
		uint32_t repeat;
		bool found;
		for ( uint32_t nuc = 0; nuc < nNucleiToPlaceNow; nuc++ ) {

			found = false;
			repeat = 0;
			while ( found == false && repeat < (2*nextFreeSiteReference) ) { //when many free sites many possible trials but at the same time usually also more chances to find free place
				luckysite = localprng.leEcuyer() * nextFreeSiteReference; //[0,1]
				if ( nucsite_avail[luckysite] != SITE_ALREADY_CONSUMED ) { //most likely case
					found = true;
					repeat = 0;
					break;
				}

				repeat++;
			}

			if ( found == true ) {
				//place myrxgpool[nucsite_avail[luckysite]]

#ifdef REPORTSTYLE_TDEPNUCLEATION
	cout << "\t\t--placing " << nuc << "\t\t" << luckysite << "\t\t" << repeat << endl;
#endif
				sim_myCA_infect_JuvenileNucleation( nucsite_avail[luckysite], myrxgpool[nucsite_avail[luckysite]].nucsite, ALMOST_FULLY_INFECTED );
				nucsite_avail[luckysite] = SITE_ALREADY_CONSUMED;
			}
		}

		delete [] nucsite_avail;
		nucsite_avail = NULL;
	}
}


void caHdl::append_to_fullrx_list( uint32_t rxfrontid )
{
	if ( nextSlotToFullRX < ntotalFullRXList ) {		//FullRX list still large enough, append
		myFullRXList[nextSlotToFullRX] = rxfrontid;		//indexing relative to first element in array with displacement sizeof(int)
		nextSlotToFullRX++;
		return;
	}

	//FullRX list is too small but it is necessary to store all cells that updateFullRX should work on, so increase size of FullRX list
	size_t oldSize = ntotalFullRXList;
	size_t newSize = oldSize + (this->transientRelCellRecaching * oldSize);
	QUICKASSERT ( newSize < CA_ALLOCATION_MAXIMUM );

#ifdef REPORTSTYLE_DEVELOPER
	cout << "\t\tALLOCATION for fullrx list, oldSize, newSize\t" << ((int) oldSize) << "\t" << ((int) newSize) << endl;
#endif

	uint32_t* new_myFullRXList = NULL;
	new_myFullRXList = new uint32_t[newSize];
	QUICKASSERT( new_myFullRXList != NULL );
	myMemGuard = myMemGuard + (newSize * sizeof(uint32_t));

	//##MK::potentially unsafe to use memcpy memcpy( new_myFullRXList, myFullRXList, oldSize * sizeof(int) );
	for ( uint32_t i = 0; i < oldSize; ++i ) {
		new_myFullRXList[i] = myFullRXList[i];
	}

	delete [] myFullRXList;
	myMemGuard = myMemGuard - ( oldSize * sizeof(uint32_t) );
	myFullRXList = new_myFullRXList;
	ntotalFullRXList = newSize;

	QUICKASSERT ( nextSlotToFullRX < ntotalFullRXList );

	//more space again, so append
	myFullRXList[nextSlotToFullRX] = rxfrontid;
	nextSlotToFullRX++;
}


void caHdl::append_to_recycling_list( uint32_t rxfrontid )
{
	if ( firstNotRecycledYet < ntotalRecyclingList ) {
		myRecyclingList[firstNotRecycledYet] = rxfrontid;
		firstNotRecycledYet++;
	}
	//else -> currently I do not recycle more than ntotalRecyclingList elements, also this list is already as large as 0.15*ntotalcells, user can optimize this structure when knowing his microstructural path function Sv(X) in more details to save additional memory therefore not further memory utilized
}


void caHdl::sim_myCA_infect_OutofRXFront_BoundsCheck ( uint32_t rxgpoolid, uint32_t rxfrontid, short dx, short dy, short dz, unsigned char direction, double rxfrac0 )
{
	//read the location of the infector
	short x = myRXFront[rxfrontid].ix;
	short y = myRXFront[rxfrontid].iy;
	short z = myRXFront[rxfrontid].iz;
	//uint32_t dgid = myRXFront[rxfrontid].mydefgid;
	//uint32_t rgid = myRXFront[rxfrontid].myrxgid;

	//add relative coordinate increment to identify position of infection target in the CA
	x += dx;
	y += dy;
	z += dz;

	//MK::modelling approach DEFAULT periodic boundary conditions in each automaton, switch contact flag if applicable
	if ( x < 0 ) {
		x += myCAGeometry.nboxedge_rd;
		myrxgpoolboundarycontact[rxgpoolid] = true;
	}
	if ( x >= myCAGeometry.nboxedge_rd ) {
		x -= myCAGeometry.nboxedge_rd;
		myrxgpoolboundarycontact[rxgpoolid] = true;
	}
	if ( y < 0 ) {
		y += myCAGeometry.nboxedge_nd;
		myrxgpoolboundarycontact[rxgpoolid] = true;
	}
	if ( y >= myCAGeometry.nboxedge_nd ) {
		y -= myCAGeometry.nboxedge_nd;
		myrxgpoolboundarycontact[rxgpoolid] = true;
	}
	if ( z < 0 ) {
		z += myCAGeometry.nboxedge_td;
		myrxgpoolboundarycontact[rxgpoolid] = true;
	}
	if ( z >= myCAGeometry.nboxedge_td ) {
		z -= myCAGeometry.nboxedge_td;
		myrxgpoolboundarycontact[rxgpoolid] = true;
	}

	//convert ix,iy,iz in implicit 3D coordinate to access mycellgrid
	uint32_t cxyz = x + (y * myCAGeometry.nboxedge_rd) + (z * myCAGeometry.nboxarea_rdnd);

	uint32_t cxyz_deformedstate = mycellgrid[cxyz];

	//MK::reconsider initialization value in the int range because currently mydefgpool.size() is by definition > 0
	//##DEBUG was >= mydefgpool.size()
	if ( cxyz_deformedstate >= nmydefgpool || cxyz_deformedstate == CURRENTLY_INFECTED ) { //|| cxyz_deformedstate  == CELL_IS_A_PARTICLE ) { 
		//get out here, cell is already infected

#ifdef REPORTSTYLE_CELLCYCLES
	cout << "\t\t\t\tREJECT infection cxyz;x;y;z;dx;dy;dz\t" << cxyz << "\t" << x << ";" << y << ";" << z << "\t" << dx << ";" << dy << ";" << dz << endl;
#endif

		return;
	}

	//obviously cell is free, hehe, go for it!
	//mark cell as infected in order to avoid multiple reinfections and further cache thrashing
	mycellgrid[cxyz] = CURRENTLY_INFECTED;


	uint32_t freeplace = get_NextFreeSlotInRXFrontAppendOrRecycle ( CELLRECYCLE ); //MK::strategy is RECYCLING because OutofRXFront infections are only called from updateFullRX

	myRXFront[freeplace].activity = ACTIVE;
	myRXFront[freeplace].infector = direction;
	myRXFront[freeplace].ix = x;
	myRXFront[freeplace].iy = y;
	myRXFront[freeplace].iz = z;
	myRXFront[freeplace].rxFrac = rxfrac0;


	//first assume most likely case - infection continues consuming the same deformed grain id, then P is the same
	myRXFront[freeplace].P = myRXFront[rxfrontid].P;

	//compromise growth of the grain from fully infected deformed regionn mydefgid into a different deformed grain
	if ( myRXFront[rxfrontid].mydefgid != cxyz_deformedstate ) {
		myRXFront[freeplace].P = calc_mobilityweight ( rxgpoolid, cxyz_deformedstate );
	}

	//myMobilityWeightMax is only relevant for the global adaptive dynamic time stepping
	if ( myRXFront[freeplace].P >= myMobilityWeightMax ) {
		myMobilityWeightMax = myRXFront[freeplace].P;
	}

	//cell carries consumed orientation along
	myRXFront[freeplace].myrxgid = rxgpoolid;
	myRXFront[freeplace].mydefgid = cxyz_deformedstate;


	//MK::assessor function get_NextFreeSlotInRXFrontAppendOrRecycle() ASSURES freeplace < nextSlotNeverActiveRXFront <= ntotalRXFront!
#ifdef REPORTSTYLE_CELLCYCLES
	cout << this->jobid << "\t\tout-of-front infection;rxgpoolid;rxfrontid;cxyz;x;y;z;freeplace" << rxgpoolid << ";" << rxfrontid << ";" << cxyz << ";" << x << ";" << y << ";" << z << ";" << freeplace << endl;
#endif
}


void caHdl::sim_myCA_infect_OutofRXFront_NoBoundsCheck ( uint32_t rxgpoolid, uint32_t rxfrontid, short dx, short dy, short dz, unsigned char direction, double rxfrac0 )
{
	//read the location of the infector
	short x = myRXFront[rxfrontid].ix;
	short y = myRXFront[rxfrontid].iy;
	short z = myRXFront[rxfrontid].iz;

	//add relative coordinate increment to identify position of infection target in the CA
	x += dx;
	y += dy;
	z += dz;

	//MK::THIS FUNCTION MUST NOT BE CALLED WHEN IT CANNOT BE ASSURED THAT i + di is located inside the automaton domain!

	//convert ix,iy,iz in implicit 3D coordinate to access mycellgrid
	uint32_t cxyz = x + (y * myCAGeometry.nboxedge_rd) + (z * myCAGeometry.nboxarea_rdnd);

	uint32_t cxyz_deformedstate = mycellgrid[cxyz];

	//MK::reconsider initialization value in the int range because currently mydefgpool.size() is by definition > 0
	if ( cxyz_deformedstate >= nmydefgpool || cxyz_deformedstate == CURRENTLY_INFECTED ) { //|| cxyz_deformedstate  == CELL_IS_A_PARTICLE ) { 
		//get out here, cell is already infected

#ifdef REPORTSTYLE_CELLCYCLES
	cout << "\t\t\t\tREJECT infection cxyz;x;y;z;dx;dy;dz\t" << cxyz << "\t" << x << ";" << y << ";" << z << "\t" << dx << ";" << dy << ";" << dz << endl;
#endif

		return;
	}

	//obviously cell is free, hehe, go for it!
	//mark cell as infected in order to avoid multiple reinfections and further cache thrashing
	mycellgrid[cxyz] = CURRENTLY_INFECTED;

	uint32_t freeplace = get_NextFreeSlotInRXFrontAppendOrRecycle ( CELLRECYCLE ); //MK::strategy is RECYCLING because OutofRXFront infections are only called from updateFullRX

	myRXFront[freeplace].activity = ACTIVE;
	myRXFront[freeplace].infector = direction;
	myRXFront[freeplace].ix = x;
	myRXFront[freeplace].iy = y;
	myRXFront[freeplace].iz = z;
	myRXFront[freeplace].rxFrac = rxfrac0;


	//first assume most likely case - infection continues consuming the same deformed grain id, then P is the same
	myRXFront[freeplace].P = myRXFront[rxfrontid].P;

	//compromise growth of the grain from fully infected deformed regionn mydefgid into a different deformed grain
	if ( myRXFront[rxfrontid].mydefgid != cxyz_deformedstate ) {
		myRXFront[freeplace].P = calc_mobilityweight ( rxgpoolid, cxyz_deformedstate );
	}

	//myMobilityWeightMax is only relevant for the global adaptive dynamic time stepping
	if ( myRXFront[freeplace].P >= myMobilityWeightMax ) {
		myMobilityWeightMax = myRXFront[freeplace].P;
	}

	//cell carries consumed orientation along
	myRXFront[freeplace].myrxgid = rxgpoolid;
	myRXFront[freeplace].mydefgid = cxyz_deformedstate;


	//MK::assessor function get_NextFreeSlotInRXFrontAppendOrRecycle() ASSURES freeplace < nextSlotNeverActiveRXFront <= ntotalRXFront!
#ifdef REPORTSTYLE_CELLCYCLES
	cout << this->jobid << "\t\tout-of-front infection;rxgpoolid;rxfrontid;cxyz;x;y;z;freeplace" << rxgpoolid << ";" << rxfrontid << ";" << cxyz << ";" << x << ";" << y << ";" << z << ";" << freeplace << endl;
#endif
}





void caHdl::sim_myCA_calcGrowthStep( void )
{
	//##consider to mount them fixed in the classHdl
	double shapefactors[27] = { FGEODIAG, FGEOEDGE, FGEODIAG, FGEOEDGE, FGEOFACE, FGEOEDGE, FGEODIAG, FGEOEDGE, FGEODIAG, 
						FGEOEDGE, FGEOFACE, FGEOEDGE, FGEOFACE,      FGEOFACE, FGEOEDGE, FGEOFACE, FGEOEDGE,	
						FGEODIAG, FGEOEDGE, FGEODIAG, FGEOEDGE, FGEOFACE, FGEOEDGE, FGEODIAG, FGEOEDGE, FGEODIAG, (1.0) }; //layers xz stacked in positive y direction...

	uint32_t IdentifiedAsActive = 0;

	//MK::reset collector array for already transformed cells
	nextSlotToFullRX = 0;

	//MK::reset RecyclingCounter as int-values from previous step are now corrupted, up to now nothing has been recycled
	nextSlotThatBecomesRecycled = 0;
	firstNotRecycledYet = 0;


	//scan the whole myRXFront on the interval [0;nextSlotNeverActiveRXFront) to 
	//	a) to find ACTIVE cells for which boundary migration is simulated during dt
	//	b) to find slots that have become INACTIVE and thus can be referenced into the RecyclingList
#ifdef REPORTSTYLE_CELLCYCLES
	cout << "\n\n\t\tstep->calc_beforeloop;nextSlotNeverActiveRXFront\t" << step << "\t" << nextSlotNeverActiveRXFront << endl;
#endif

	//MK::when nucleation should be time-dependent it is vital to insert
	//the nuclei here at the latest, otherwise nextSlotNeverActiveRX would become updated and as such corrupted...
	double dVrx = 0.0;
	double v = 0.0;

#ifdef FRAGMENTATION_RADAR
	uint32_t HowManyBlocks = EXPECTEDSHARE * (double) myCAGeometry.nboxvol_rdndtd / BLOCKLENGTH;
	uint32_t FirstUnconsideredRXFrontCell = HowManyBlocks * BLOCKLENGTH; //operating on [0, First...)
	uint32_t CurrentBlockID = 0;
	uint32_t HowManyCellsInTotalInTheBlock = 0;
	uint32_t HowManyCellsActiveInTheBlock = 0;
	unsigned char* RGBA = NULL;
	RGBA = new unsigned char[4*HowManyBlocks];
	QUICKASSERT ( RGBA != NULL );
	for ( uint32_t b = 0; b < HowManyBlocks; b++ ) {
		RGBA[4*b+0] = 255;		RGBA[4*b+1] = 255;		RGBA[4*b+2] = 255;		RGBA[4*b+3] = 255; //white default value
	}

	struct loginfo_fragmentation_radar ping;
	ping.localtstep = this->step;
	ping.blocklength = BLOCKLENGTH;
	ping.nblocks = HowManyBlocks;
	ping.rgba = RGBA;
	this->myrxfrontfragmentation.push_back ( ping );
#endif


	for ( uint32_t c = 0; c < nextSlotNeverActiveRXFront; c++ ) {

#ifdef FRAGMENTATION_RADAR
		CurrentBlockID = c / BLOCKLENGTH;
		HowManyCellsInTotalInTheBlock++;
#endif

		//MK::ACTIVE is the most likely case in a well defragmented list, as then most cells are active
		if ( myRXFront[c].activity == ACTIVE ) {
			IdentifiedAsActive++;

#ifdef FRAGMENTATION_RADAR
			HowManyCellsActiveInTheBlock++;
#endif

			//read cached mobilityweight from local cell container that is assumed fixed since the infection
			//v = mp Turnbull linearized rate equation model
			v = get_currentintrinsicmobility ( myRXFront[c].P );

#ifdef REPORTSTYLE_CELLCYCLES
			cout << "\t\tstep->calc_beforecontinue is active;v;nactive\t" << step << "\t" << c << "\t" << v << "\t" << nCurrentlyActive << "\tP=" << (myRXFront[c].P) << endl;
#endif

			//compute migration velocity
			v = v * ( (Gbhalfsq * get_rho ( myRXFront[c].mydefgid )) - get_zener () );

			//scale velocity according to infection direction and kernel function to assure self-similar growth
			v *= shapefactors[myRXFront[c].infector];


			myRXFront[c].rxFrac = myRXFront[c].rxFrac + (v * dt * _cellsize); //v*(max*cellsize/vmax)*v

			dVrx = dVrx + (v * dt * _cellsize);

			//is the cell finally full?
			if ( myRXFront[c].rxFrac >= 1.0 ) {
				//cout << "APPENDTOFULLRX\t\t" << c << ";" << myRXFront[c].rxFrac << endl;
				append_to_fullrx_list ( c );
			}

#ifdef REPORTSTYLE_CELLCYCLES
			cout << "\t\tstep->calc_beforecontinue is active;rxFrac;v;nactive\t" << step << "\t" << c << "\t" << myRXFront[c].rxFrac << "\t" << v << "\t" << nCurrentlyActive << endl;
			//QUICKASSERT ( v > 0.0 );
#endif
			//MK:: else construct to expensive

#ifdef FRAGMENTATION_RADAR
			if ( HowManyCellsInTotalInTheBlock >= BLOCKLENGTH ) { //write-out and reset
				RGBA[4*CurrentBlockID+0] = 255.0 * (1.0 - pow( ( (double) HowManyCellsActiveInTheBlock / (double) HowManyCellsInTotalInTheBlock ), REDSTRETCH ));
				RGBA[4*CurrentBlockID+1] = 0;
				RGBA[4*CurrentBlockID+2] = 0;
//cout << 255.0 * (1.0 - pow( ( (double) HowManyCellsActiveInTheBlock / (double) HowManyCellsInTotalInTheBlock ), REDSTRETCH )) << "-" << HowManyCellsActiveInTheBlock << "-" << HowManyCellsInTotalInTheBlock << ";";
				HowManyCellsActiveInTheBlock = 0;
				HowManyCellsInTotalInTheBlock = 0;
			}
#endif


			continue;
		} 
		//else {
			//else cell seems inactive, okay remember index c in RecyclingList

#ifdef FRAGMENTATION_RADAR
			if ( HowManyCellsInTotalInTheBlock >= BLOCKLENGTH ) { //write-out and reset
				RGBA[4*CurrentBlockID+0] = 255.0 * (1.0 - pow( ( (double) HowManyCellsActiveInTheBlock / (double) HowManyCellsInTotalInTheBlock ), REDSTRETCH ));
				RGBA[4*CurrentBlockID+1] = 0;
				RGBA[4*CurrentBlockID+2] = 0;
//cout << 255.0 * (1.0 - pow( ( (double) HowManyCellsActiveInTheBlock / (double) HowManyCellsInTotalInTheBlock ), REDSTRETCH )) << "-" << HowManyCellsActiveInTheBlock << "-" << HowManyCellsInTotalInTheBlock << ";";
				HowManyCellsActiveInTheBlock = 0;
				HowManyCellsInTotalInTheBlock = 0;
			}
#endif


			append_to_recycling_list ( c );
		//}

	} //analyze all cells in [0, nextSlotNeverActiveRXFront)

#ifdef FRAGMENTATION_RADAR
	//enforce export of last block
//cout << 255.0 * (1.0 - pow( ( (double) HowManyCellsActiveInTheBlock / (double) HowManyCellsInTotalInTheBlock ), REDSTRETCH )) << ";" << endl;
	RGBA[4*CurrentBlockID+0] = 255.0 * (1.0 - pow( ( (double) HowManyCellsActiveInTheBlock / (double) HowManyCellsInTotalInTheBlock ), REDSTRETCH ));
	RGBA[4*CurrentBlockID+1] = 0;
	RGBA[4*CurrentBlockID+2] = 0;
cout << "Pinging fragmentation radar data step/blocklength/nblocks\t\t" << ping.localtstep << ";" << ping.blocklength << ";" << ping.nblocks << "\t\tnextSlotNeverActiveRXFront;InactiveFractionSoBranchmispredict\t\t" << nextSlotNeverActiveRXFront << ";" << ( (double) (nextSlotNeverActiveRXFront - IdentifiedAsActive) / (double) nextSlotNeverActiveRXFront ) << endl;
#endif


	this->dXstep = dVrx;
	this->Sv = IdentifiedAsActive;
	this->nCurrentlyActive = IdentifiedAsActive;

#ifdef REPORTSTYLE_CELLCYCLES
	cout << "\t\tstep->calc_end;step;dXstep;Sv;\t" << step << "\t" << dXstep << "\t" << Sv << endl;
#endif
}


void caHdl::sim_myCA_updateFullRX( void )
{
	//MK::CALL ONLY AFTER calcGrowthStep
	//myFullRXList in the interval [0;nextSlotToFullRX) is a compact list of entries dereferencing cells from myRXFront that are completely recrystallized and can be switched off

	//currently //compromise //originalfactors //seemed optimized factors but are not //original factors for velocity fgeo 1 in <100>, 1/2^0.5 <110>, 1/3^0.5 <111>
	const double overFac1 = 0.0682; //0.04946756; //0.03738274; //0.0523448297609848;  //###MK20121211, 26NN, 8 face, 12 edge and 6 diag neighbors each one getting a portion of the overshoot, all other overshoot is discarded
	const double overFac2 = 0.0542; //0.0387; //0.03949696; //0.04059546; //0.0370133840840478;
	const double overFac3 = 0.0484; //0.0307; //0.02865390; //0.03606975; //0.0302213015531897;
	double overshoot, carryOver1, carryOver2, carryOver3;

	uint32_t nUpdatedCells = 0;
	uint32_t mycageo_nboxedge_rd = myCAGeometry.nboxedge_rd;

	//MK::most importantly now: if the cell is located inside the automaton domain and a fixed kernel, ie. VONNEUMANN OR MOORE is utilized
	//in almost all cases the infected cell needs not mapping of periodic boundary conditions, this can be exploited by first categorizing
	//the infector cell and then to branch either in the NoBoundsCheck avoid 6 branches per infection per direction! and the BoundsChecking version for the
	//few cases of cells located in the outer shell of the automaton in a layer of two cells thickness for VONNEUMANN AND MOORE kernels!
	uint32_t lowerlimit = 2;
	uint32_t upperlimit = mycageo_nboxedge_rd - 2;


	uint32_t mycageo_nboxarea_rdnd = myCAGeometry.nboxarea_rdnd;
	unsigned char location;

	for ( uint32_t c = 0; c < nextSlotToFullRX; c++ ) {
		uint32_t crxfrontid = myFullRXList[c];

		uint32_t rgpid = myRXFront[crxfrontid].myrxgid; //runs from [0;myrxgpool.size) !

		#ifdef REPORTSTYLE_CELLCYCLES
			cout << "\t\t\t->step->updt_priorinfects;c;crxfrontid;myRXfront[crxfrontid].rxFrac;nextSlotToFullRX" << step << ";" << c << ";" << crxfrontid << ";" << myRXFront[crxfrontid].rxFrac << ";" << nextSlotToFullRX << endl;
		#endif
		QUICKASSERT( myRXFront[crxfrontid].rxFrac >= 1.0 ); //##DEBUG

		overshoot = myRXFront[crxfrontid].rxFrac - 1.0;
		//partition normalized overshoot
		carryOver1 = overshoot * overFac1;
		carryOver2 = overshoot * overFac2;
		carryOver3 = overshoot * overFac3;

		//##MK::update potential shape axes-aligned bounding boxes for shape tracking --> inject code here

		//infect all untransformed neighbors
		//##MK::ASSURE THAT COMPILER OPTIMIZATIONS DO NOT PARALLELIZE THIS PART SO THAT IN A WORST CASE ALL INFECT FUNCTIONS TRY TO RECYCLE THE SAME CELL STRUCT IN THE MYFRONTSLIST!

		location = INFECTOR_DEEP_IN_THE_CUBE; //MK::assuming cubic domain!
		if ( myRXFront[crxfrontid].ix < lowerlimit || myRXFront[crxfrontid].ix > upperlimit ) location = INFECTOR_CLOSE_TO_THE_DOMAINWALL;
		if ( myRXFront[crxfrontid].iy < lowerlimit || myRXFront[crxfrontid].iy > upperlimit ) location = INFECTOR_CLOSE_TO_THE_DOMAINWALL;
		if ( myRXFront[crxfrontid].iz < lowerlimit || myRXFront[crxfrontid].iz > upperlimit ) location = INFECTOR_CLOSE_TO_THE_DOMAINWALL;

		if ( location == INFECTOR_DEEP_IN_THE_CUBE ) { //most likely case

			//MK::changed infection order according to LM,KL to help increasing cache locality on mycellgrid
			sim_myCA_infect_OutofRXFront_NoBoundsCheck ( rgpid, crxfrontid, -1, -1, 0, 16, carryOver2 );
			sim_myCA_infect_OutofRXFront_NoBoundsCheck ( rgpid, crxfrontid, 0, -1, 0, 15, carryOver1 );
			sim_myCA_infect_OutofRXFront_NoBoundsCheck ( rgpid, crxfrontid, +1, -1, 0, 14, carryOver2 );
			sim_myCA_infect_OutofRXFront_NoBoundsCheck ( rgpid, crxfrontid, -1, 0, 0, 13, carryOver1 );
			sim_myCA_infect_OutofRXFront_NoBoundsCheck ( rgpid, crxfrontid, +1, 0, 0, 12, carryOver1 );
			sim_myCA_infect_OutofRXFront_NoBoundsCheck ( rgpid, crxfrontid, -1, +1, 0, 11, carryOver2 );
			sim_myCA_infect_OutofRXFront_NoBoundsCheck ( rgpid, crxfrontid, 0, +1, 0, 10, carryOver1 );
			sim_myCA_infect_OutofRXFront_NoBoundsCheck ( rgpid, crxfrontid, +1, +1, 0, 9, carryOver2 );

			sim_myCA_infect_OutofRXFront_NoBoundsCheck ( rgpid, crxfrontid, -1, +1, -1, 19, carryOver3 );
			sim_myCA_infect_OutofRXFront_NoBoundsCheck ( rgpid, crxfrontid, 0, +1, -1, 18, carryOver2 );
			sim_myCA_infect_OutofRXFront_NoBoundsCheck ( rgpid, crxfrontid, +1, +1, -1, 17, carryOver3 );
			sim_myCA_infect_OutofRXFront_NoBoundsCheck ( rgpid, crxfrontid, -1, -1, -1, 25, carryOver3 );
			sim_myCA_infect_OutofRXFront_NoBoundsCheck ( rgpid, crxfrontid, 0, -1, -1, 24, carryOver2 );
			sim_myCA_infect_OutofRXFront_NoBoundsCheck ( rgpid, crxfrontid, +1, -1, -1, 23, carryOver3 );
			sim_myCA_infect_OutofRXFront_NoBoundsCheck ( rgpid, crxfrontid, -1, 0, -1, 22, carryOver2 );
			sim_myCA_infect_OutofRXFront_NoBoundsCheck ( rgpid, crxfrontid, 0, 0, -1, 21, carryOver1 );
			sim_myCA_infect_OutofRXFront_NoBoundsCheck ( rgpid, crxfrontid, +1, 0, -1, 20, carryOver2 );

			sim_myCA_infect_OutofRXFront_NoBoundsCheck ( rgpid, crxfrontid, -1, -1, +1, 8, carryOver3 );
			sim_myCA_infect_OutofRXFront_NoBoundsCheck ( rgpid, crxfrontid, 0, -1, +1, 7, carryOver2 );
			sim_myCA_infect_OutofRXFront_NoBoundsCheck ( rgpid, crxfrontid, +1, -1, +1, 6, carryOver3 );
			sim_myCA_infect_OutofRXFront_NoBoundsCheck ( rgpid, crxfrontid, -1, 0, +1, 5, carryOver2 );
			sim_myCA_infect_OutofRXFront_NoBoundsCheck ( rgpid, crxfrontid, 0, 0, +1, 4, carryOver1 );
			sim_myCA_infect_OutofRXFront_NoBoundsCheck ( rgpid, crxfrontid, +1, 0, +1, 3, carryOver2 );
			sim_myCA_infect_OutofRXFront_NoBoundsCheck ( rgpid, crxfrontid, -1, +1, +1, 2, carryOver3 );
			sim_myCA_infect_OutofRXFront_NoBoundsCheck ( rgpid, crxfrontid, 0, +1, +1, 1, carryOver2 );
			sim_myCA_infect_OutofRXFront_NoBoundsCheck ( rgpid, crxfrontid, +1, +1, +1, 0, carryOver3 );
			//MK::following code NEEDS to be executed AFTER all these infection attempts

			//deactivate cell 
			uint32_t cxyz = myRXFront[crxfrontid].ix + (myRXFront[crxfrontid].iy * mycageo_nboxedge_rd) + (myRXFront[crxfrontid].iz * mycageo_nboxarea_rdnd);

			//mark representor in the cellgrid as FULLY_RECRYSTALLIZED by assigning id of RXgrains
			mycellgrid[cxyz] = (nmydefgpool + rgpid); //##mydefgpool.size() + rgpid);

			//texture bookkeeping
			myrxgpool[myRXFront[crxfrontid].myrxgid].cellcount += 1;
			mydefgpool[myRXFront[crxfrontid].mydefgid].cellcount -= 1;

			//cell flagged as "Freiwild" and thought of Green Mile healed from all infections...
			myRXFront[crxfrontid].activity = INACTIVE;
			myRXFront[crxfrontid].rxFrac = 0.0;

#ifdef REPORTSTYLE_CELLCYCLES
			cout << "\t\tstep->updt_afterinfects;crxfrontid;activity;ACTIVEwouldbeMarkedAs;nextRXSlotNeverActive\t" << step << "\t" << crxfrontid << "\t" << myRXFront[crxfrontid].activity << "\t" << ACTIVE << "\t" << nextSlotNeverActiveRXFront << endl;
#endif

			nUpdatedCells++;

			continue;
		}

		//INFECTOR_CLOSE_TO_THE_DOMAINWALL

			//MK::changed infection order according to LM,KL to help increasing cache locality on mycellgrid
			sim_myCA_infect_OutofRXFront_BoundsCheck ( rgpid, crxfrontid, -1, -1, 0, 16, carryOver2 );
			sim_myCA_infect_OutofRXFront_BoundsCheck ( rgpid, crxfrontid, 0, -1, 0, 15, carryOver1 );
			sim_myCA_infect_OutofRXFront_BoundsCheck ( rgpid, crxfrontid, +1, -1, 0, 14, carryOver2 );
			sim_myCA_infect_OutofRXFront_BoundsCheck ( rgpid, crxfrontid, -1, 0, 0, 13, carryOver1 );
			sim_myCA_infect_OutofRXFront_BoundsCheck ( rgpid, crxfrontid, +1, 0, 0, 12, carryOver1 );
			sim_myCA_infect_OutofRXFront_BoundsCheck ( rgpid, crxfrontid, -1, +1, 0, 11, carryOver2 );
			sim_myCA_infect_OutofRXFront_BoundsCheck ( rgpid, crxfrontid, 0, +1, 0, 10, carryOver1 );
			sim_myCA_infect_OutofRXFront_BoundsCheck ( rgpid, crxfrontid, +1, +1, 0, 9, carryOver2 );

			sim_myCA_infect_OutofRXFront_BoundsCheck ( rgpid, crxfrontid, -1, +1, -1, 19, carryOver3 );
			sim_myCA_infect_OutofRXFront_BoundsCheck ( rgpid, crxfrontid, 0, +1, -1, 18, carryOver2 );
			sim_myCA_infect_OutofRXFront_BoundsCheck ( rgpid, crxfrontid, +1, +1, -1, 17, carryOver3 );
			sim_myCA_infect_OutofRXFront_BoundsCheck ( rgpid, crxfrontid, -1, -1, -1, 25, carryOver3 );
			sim_myCA_infect_OutofRXFront_BoundsCheck ( rgpid, crxfrontid, 0, -1, -1, 24, carryOver2 );
			sim_myCA_infect_OutofRXFront_BoundsCheck ( rgpid, crxfrontid, +1, -1, -1, 23, carryOver3 );
			sim_myCA_infect_OutofRXFront_BoundsCheck ( rgpid, crxfrontid, -1, 0, -1, 22, carryOver2 );
			sim_myCA_infect_OutofRXFront_BoundsCheck ( rgpid, crxfrontid, 0, 0, -1, 21, carryOver1 );
			sim_myCA_infect_OutofRXFront_BoundsCheck ( rgpid, crxfrontid, +1, 0, -1, 20, carryOver2 );

			sim_myCA_infect_OutofRXFront_BoundsCheck ( rgpid, crxfrontid, -1, -1, +1, 8, carryOver3 );
			sim_myCA_infect_OutofRXFront_BoundsCheck ( rgpid, crxfrontid, 0, -1, +1, 7, carryOver2 );
			sim_myCA_infect_OutofRXFront_BoundsCheck ( rgpid, crxfrontid, +1, -1, +1, 6, carryOver3 );
			sim_myCA_infect_OutofRXFront_BoundsCheck ( rgpid, crxfrontid, -1, 0, +1, 5, carryOver2 );
			sim_myCA_infect_OutofRXFront_BoundsCheck ( rgpid, crxfrontid, 0, 0, +1, 4, carryOver1 );
			sim_myCA_infect_OutofRXFront_BoundsCheck ( rgpid, crxfrontid, +1, 0, +1, 3, carryOver2 );
			sim_myCA_infect_OutofRXFront_BoundsCheck ( rgpid, crxfrontid, -1, +1, +1, 2, carryOver3 );
			sim_myCA_infect_OutofRXFront_BoundsCheck ( rgpid, crxfrontid, 0, +1, +1, 1, carryOver2 );
			sim_myCA_infect_OutofRXFront_BoundsCheck ( rgpid, crxfrontid, +1, +1, +1, 0, carryOver3 );
			//MK::following code NEEDS to be executed AFTER all these infection attempts

			//deactivate cell 
			uint32_t cxyz = myRXFront[crxfrontid].ix + (myRXFront[crxfrontid].iy * mycageo_nboxedge_rd) + (myRXFront[crxfrontid].iz * mycageo_nboxarea_rdnd);

			//mark representor in the cellgrid as FULLY_RECRYSTALLIZED by assigning id of RXgrains
			mycellgrid[cxyz] = (nmydefgpool + rgpid); //##mydefgpool.size() + rgpid);

			//texture bookkeeping
			myrxgpool[myRXFront[crxfrontid].myrxgid].cellcount += 1;
			mydefgpool[myRXFront[crxfrontid].mydefgid].cellcount -= 1;

			//cell flagged as "Freiwild" and thought of Green Mile healed from all infections...
			myRXFront[crxfrontid].activity = INACTIVE;
			myRXFront[crxfrontid].rxFrac = 0.0;

#ifdef REPORTSTYLE_CELLCYCLES
			cout << "\t\tstep->updt_afterinfects;crxfrontid;activity;ACTIVEwouldbeMarkedAs;nextRXSlotNeverActive\t" << step << "\t" << crxfrontid << "\t" << myRXFront[crxfrontid].activity << "\t" << ACTIVE << "\t" << nextSlotNeverActiveRXFront << endl;
#endif

			nUpdatedCells++;

	} //for all cells


	Xcells += (double) nUpdatedCells;
	X = (Xcells / ((double) myCAGeometry.nboxvol_rdndtd));
}



void caHdl::log_colorize_myoripool_ipfz( void )
{
	uint32_t nori = this->myoripool.size();

	for ( uint32_t oi = 0; oi < myoripool.size(); oi++) {
		unsigned char rgb[3] = { UCHAR_RANGE_MIN, UCHAR_RANGE_MIN, UCHAR_RANGE_MIN };
		double locationinsst[2] = {0.0, 0.0};

		double p1 = myoripool[oi].bunge1;
		double p2 = myoripool[oi].bunge2;
		double p3 = myoripool[oi].bunge3;

		bunge2ipfz( p1, p2, p3, rgb, locationinsst );

		myoripool[oi].RGBA[RED] = rgb[0];
		myoripool[oi].RGBA[GREEN] = rgb[1];
		myoripool[oi].RGBA[BLUE] = rgb[2];
		myoripool[oi].RGBA[ALPHA] = UCHAR_RANGE_MAX;
//cout << "Bunge2IPFZ=" << oi << "\t" << (int) rgb[RED] << ";" << (int) rgb[GREEN] << (int) rgb[BLUE] << endl;
	}
}


void caHdl::log_initialization( void )
{
	//user defined recrystallized fraction of desired logpoints normalized, so multiply with boxdimensions
	double boxvolume = myCAGeometry.nboxvol_rdndtd;

	for ( uint32_t d = 0; d < this->defragRXFront_atthisX.size(); ++d ) {
		defragRXFront_atthisX[d] = defragRXFront_atthisX[d] * boxvolume;
	}

	for ( uint32_t s = 0; s < this->output_atthisX.size(); ++s ) {
		output_atthisX[s] = output_atthisX[s] * boxvolume;
	}

	for ( uint32_t r = 0; r < this->rendering_atthisX.size(); ++r ) {
		rendering_atthisX[r] = rendering_atthisX[r] * boxvolume;
	}


	//initialize the color coding
	if ( outopt_localrenderhow == RENDERING_MS2D || outopt_localrenderhow == RENDERING_MS3D ) {
		if ( this->outopt_localrendercolor == RENDERING_COLOR_IPFZ )
			log_colorize_myoripool_ipfz();
	}

	//implement intelligent precaching of logcontainers rxfrontstats and grainevo
}


void caHdl::log_rxfrontstats( void )
{
	//collect current information
	struct loginfo_rxfrontstats_ca memca;

	memca.localtime = t;
	memca.localX = X;
	memca.localmemory = myMemGuard;
	memca.localPmax = myMobilityWeightMax;
	memca.localstep = step;
	memca.localSv = Sv;

	memca.ntotalRXFront = this->ntotalRXFront;
	memca.nextSlotNeverActiveRXFront = this->nextSlotNeverActiveRXFront;
	memca.nCurrentlyActive = this->nCurrentlyActive;

	memca.ntotalFullRXList = this->ntotalFullRXList;
	memca.nextSlotToFullRX = this->nextSlotToFullRX;

	memca.ntotalRecyclingList = this->ntotalRecyclingList;
	memca.nextSlotThatBecomesRecycled = this->nextSlotThatBecomesRecycled;
	memca.firstNotRecycledYet = this->firstNotRecycledYet;


	//push_back copyconstruction of this stack struct
	myrxfrontstatus.push_back( memca );

	loginfo_rxfrontstats_cnt++;
}


void caHdl::log_grainevo( void )
{
	//collect current microstructure state
	//##MK::further optimization potential by collecting only changes of the structure!
	uint32_t ndisjoint_grains = mydefgpool.size() + myrxgpool.size();

	//temporary stack object invoking vector copy constructor of uninitialized grevo
	struct loginfo_grainevo_ca grevo;
	mygrainevolution.push_back ( grevo );

	uint32_t i = mygrainevolution.size() - 1;

	mygrainevolution[i].localtime = t;
	mygrainevolution[i].localX = X;
	mygrainevolution[i].locallogcnt = loginfo_grainevo_cnt;
	mygrainevolution[i].localtstep = step;
	mygrainevolution[i].localSv = Sv;
	mygrainevolution[i].nlocaldata = ndisjoint_grains;

	uint32_t* cellcount_bucket = NULL;
	cellcount_bucket = new uint32_t[ndisjoint_grains];
	QUICKASSERT ( cellcount_bucket != NULL );
	myMemGuard = myMemGuard + ( ndisjoint_grains * sizeof(uint32_t) );

	//log discrete volume of deformed grains
	uint32_t ndefg = mydefgpool.size();
	for ( uint32_t dg = 0; dg < ndefg; dg++ ) {
		cellcount_bucket[dg] = mydefgpool[dg].cellcount;
	}

	//log discrete volume of rx grains
	for ( uint32_t rg = ndefg; rg < ndisjoint_grains; rg++ ) {
		cellcount_bucket[rg] = myrxgpool[rg-ndefg].cellcount;
	}

	//here pointer to heap segment is carried over
	mygrainevolution[i].localdatabucket = cellcount_bucket;

#ifdef REPORTSTYLE_DEVELOPER
	cout << this->jobid << " logged grain evolution at step = " << this->step << " for the cnt = " << loginfo_grainevo_cnt << "-th time at desired X = " << UserDefLogPoint_X_Output[loginfo_grainevo_cnt] << endl;
#endif
	loginfo_grainevo_cnt++;
}

void caHdl::write_zsection_coloring_grainids( void )
{
	//##MK
	cout << "WARNING::Returning CURRENTLY there is not option to colorcode grains other than IPFZ" << endl;
}


void caHdl::write_zsection_coloring_ipfz( void ) 
{
	//render a xy section sliced at z = DEFAULT_ZSECTIONING
	uint32_t imgx = this->myCAGeometry.nboxedge_rd;
	uint32_t imgy = this->myCAGeometry.nboxedge_nd;
	uint32_t imgxy = imgx * imgy;
	uint32_t imgz = (uint32_t) ( DEFAULT_ZSECTIONING_ZPOS * ((double) myCAGeometry.nboxedge_td));

	//allocate memory to store the resulting image
	unsigned char* img_rgba = new unsigned char [imgx*imgy*4];

	ostringstream fname;
	fname << "SCORE." << this->myensHdl->simid << ".JobID." << this->jobid << ".XYIter." << this->step << ".X.";
	uint32_t XX = (X * 1000.0);
	if ( XX < 10 )					fname << "000" << XX;
	if ( XX >= 10 && XX < 100 )		fname << "00" << XX;
	if ( XX >= 100 && XX < 1000 )	fname << "0" << XX;
	if ( XX >= 1000 )				fname << XX;
	fname << ".IPFZ.png";

cout << "Rendering a section__" << fname.str().c_str() << "__ at " << imgx << ";" << imgy << ";" << imgz << endl;
	//render microstructure
	uint32_t cellgridoff;
	uint32_t imggridpoint;
	uint32_t gid;
	uint32_t nmydgpool = this->nmydefgpool;
	int oid;

//fill image with white black for deformed structure
	for ( uint32_t y = 0; y < imgy; y++ ) {
		for ( uint32_t x = 0; x < imgx; x++ ) {
			imggridpoint = 4*(x + (y*imgx));
			img_rgba[imggridpoint + REDCHAN ] = BLACK;
			img_rgba[imggridpoint + GREENCHAN ] = BLACK;
			img_rgba[imggridpoint + BLUECHAN ] = BLACK;
			img_rgba[imggridpoint + ALPHACHAN ] = UCHAR_RANGE_MAX;
		}
	}


	for ( uint32_t y = 0; y < imgy; y++ ) {
		cellgridoff = (y * imgx) + (imgz * imgxy);

		for ( uint32_t x = 0; x < imgx; x++ ) {
			gid = mycellgrid[x+cellgridoff];

//cout << "x/y/gid/nmydefgpool\t\t" << x << "\t" << y << "\t" << gid << "\t" << this->nmydefgpool << endl;

			if ( gid >= nmydgpool && gid != CURRENTLY_INFECTED ) { //|| gid != CELL_IS_A_PARTICLE ) { 
				oid = myrxgpool[gid].caori;

				imggridpoint = 4*(x + (y*imgx));

//cout << y << "\t\t" << x << "\t\t" << oid << "\t\t" << myoripool[oid].RGBA[REDCHAN] << ";"<< myoripool[oid].RGBA[GREENCHAN] << ";"<< myoripool[oid].RGBA[BLUECHAN] << ";"<< myoripool[oid].RGBA[ALPHACHAN] << endl;

				img_rgba[imggridpoint + REDCHAN ] = myoripool[oid].RGBA[REDCHAN];
				img_rgba[imggridpoint + GREENCHAN ] = myoripool[oid].RGBA[GREENCHAN];
				img_rgba[imggridpoint + BLUECHAN ] = myoripool[oid].RGBA[BLUECHAN];
				img_rgba[imggridpoint + ALPHACHAN ] = myoripool[oid].RGBA[ALPHACHAN];
			}
		}
	}



	//encode the PNG file by calling the lodePNG library
	lodepng::encode( fname.str().c_str(), img_rgba , imgx, imgy );


	delete [] img_rgba;
}


//##DEBUG::function currently improperly embedded in network of MPI ensembleHdl...
void caHdl::write_voxeldata_coloring_grainids( void )
{
	//TARGET FORMAT IS INT FOR PARAVIEW OR AVIZO ONLY RECRYSTALLIZED GRAINS ARE VISUALIZED CAN BE KEYED OUT BY ZERO
	//CELL_IS_PARTICLE	0
	//CELL_IS_INFECTED	1
	//DEFORMED 			1+mydefgid
	//RECRYSTALLIZED	1+mydefgid.size()+myrxgid
	//<MPI_INT>uint

	MPI_File msFileHdl;
	MPI_Status msFileStatus;

	// create C-consistent file name for MPI I/O
	stringstream msFileName;
	msFileName << "SCORE." << myensHdl->simid << ".JobID." <<  this->jobid << ".Iter." << this->step << ".X.";
	uint32_t XX = (X * 1000.0);
	if ( XX < 10 )					msFileName << "000" << XX;
	if ( XX >= 10 && XX < 100 )		msFileName << "00" << XX;
	if ( XX >= 100 && XX < 1000 )	msFileName << "0" << XX;
	if ( XX >= 1000 )				msFileName << XX;
	msFileName << ".MS3DGrainID.raw";

	int msFileNameLength = msFileName.str().size();
	char* CmsFileName = new char[msFileNameLength+1];
	strcpy(CmsFileName, msFileName.str().c_str());

	//working all nodes open the file in create and write-only mode
	MPI_File_open(MPI_COMM_SELF, CmsFileName, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &msFileHdl);

	//MK:: with 64-bit architectures, offsets calculations can be performed as longs 2^64 bit thus on petabyte files...
	//in mpi.h MPI_Offset is defined as an __int64 which is long long, thus we can jump much more than 2^32 directly when unsigned int would be utilized
	//__int64 totalOffset = 0;
	int totalOffset = 0;
	MPI_File_seek( msFileHdl, totalOffset, MPI_SEEK_SET );

	//single automaton smaller than uint32_t range on shared memory in one process!
	uint32_t cxyz = 0;
	int mpi_io_cell_size = sizeof(int); //4
	uint32_t gridvalue = RAW_PARTICLE;
	uint32_t nx = myCAGeometry.nboxedge_rd;
	uint32_t ny = myCAGeometry.nboxedge_nd;
	uint32_t nz = myCAGeometry.nboxedge_td;
	uint32_t nxy = myCAGeometry.nboxarea_rdnd;
	uint32_t nxyz = myCAGeometry.nboxvol_rdndtd;


	//MPI_IO_Cell *cellLine = new MPI_IO_Cell[this->xPer]; //all lines in one node have the same size!, cell.OriIndex is now a grain index from grains.all!
	uint32_t* rawdata = NULL;
	rawdata = new uint32_t[nxyz];
	QUICKASSERT ( rawdata != NULL );

	//assign values cache aware of data layout in memory
	for (uint32_t z = 0; z < nz; ++z) {
		uint32_t offsetxy = z * nxy;
		for (uint32_t y = 0; y < ny; ++y) {
			uint32_t offsety = y * nx;
			for (uint32_t x = 0; x < nx; ++x) {
				cxyz = x + offsety + offsetxy;	//could be made omitting nx further additions

				gridvalue = mycellgrid[cxyz];

				//cout << cxyz << "\t\t" << gridvalue << endl;

				//mutually exclusive
				if ( gridvalue == CELL_IS_A_PARTICLE ) { 
					rawdata[cxyz] = 0; //RAW_PARTICLE;
					continue;
				}
				if ( gridvalue == CURRENTLY_INFECTED ) { 
					rawdata[cxyz] = 1; //RAW_INFECTED;
					continue;
				}
				//obviously deformed=untouched or recrystallized already
				rawdata[cxyz] = 2 + gridvalue;
				//rawdata[cxyz] += gridvalue; //that is > filevalue = ####defg or rx;
			}
		}
	}

	//push in file at once, runtime environment will optimize, further potential by nowait multiple write from the thread, but tricky...
	//MPI_File_write_at(msFileHdl, totalOffset, rawdata, nxyz, MPI_INT, &msFileStatus);
	MPI_File_write(msFileHdl, rawdata, nxyz, MPI_INT, &msFileStatus);

	delete [] rawdata;
	delete [] CmsFileName;

//cout << myensHdl->myRank << " successful MPI I/O in timestep = " << this->step << endl;

	MPI_File_close(&msFileHdl); //no Barrier as MPI_COMM_SELF
}


void caHdl::write_voxeldata_coloring_ipfz( void )
{
	//TARGET FORMAT IS INT FOR PARAVIEW OR AVIZO ONLY RECRYSTALLIZED GRAINS (rxFrac < 1.0) ARE colored white
	//CELL_IS_PARTICLE	0
	//CELL_IS_INFECTED	1
	//DEFORMED 			1+mydefgid
	//RECRYSTALLIZED	1+mydefgid.size()+myrxgid
	//<MPI_INT>uint

	MPI_File msFileHdl;
	MPI_Status msFileStatus;

	// create C-consistent file name for MPI I/O
	stringstream msFileName;
	msFileName << "SCORE." << myensHdl->simid << ".JobID." <<  this->jobid << ".Iter." << this->step << ".X.";
	uint32_t XX = (X * 1000.0);
	if ( XX < 10 )					msFileName << "000" << XX;
	if ( XX >= 10 && XX < 100 )		msFileName << "00" << XX;
	if ( XX >= 100 && XX < 1000 )	msFileName << "0" << XX;
	if ( XX >= 1000 )				msFileName << XX;
	msFileName << ".MS3DIPFZ.raw";
	
	int msFileNameLength = msFileName.str().size();
	char* CmsFileName = new char[msFileNameLength+1];
	strcpy(CmsFileName, msFileName.str().c_str());

	//working all nodes open the file in create and write-only mode
	MPI_File_open(MPI_COMM_SELF, CmsFileName, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &msFileHdl);

	//MK:: with 64-bit architectures, offsets calculations can be performed as longs 2^64 bit thus on petabyte files...
	//in mpi.h MPI_Offset is defined as an __int64 which is long long, thus we can jump much more than 2^32 directly when unsigned int would be utilized
	int totalOffset = 0;
	MPI_File_seek( msFileHdl, totalOffset, MPI_SEEK_SET );

	//single automaton smaller than uint32_t range on shared memory in one process!
	uint32_t cxyz = 0;
	uint32_t gridvalue = RAW_PARTICLE;
	uint32_t nx = myCAGeometry.nboxedge_rd;
	uint32_t ny = myCAGeometry.nboxedge_nd;
	uint32_t nz = myCAGeometry.nboxedge_td;
	uint32_t nxy = myCAGeometry.nboxarea_rdnd;
	uint32_t nxyz = myCAGeometry.nboxvol_rdndtd;


	unsigned char* rawdata = NULL;
	rawdata = new unsigned char[RGB * nxyz];
	QUICKASSERT ( rawdata != NULL );

	uint32_t rxgid;
	uint32_t oi;
	uint32_t ndg = this->nmydefgpool;

	//assign values cache aware of data layout in memory
	for (uint32_t z = 0; z < nz; ++z) {
		uint32_t offsetxy = z * nxy;
		for (uint32_t y = 0; y < ny; ++y) {
			uint32_t offsety = y * nx;
			for (uint32_t x = 0; x < nx; ++x) {
				cxyz = x + offsety + offsetxy;	//could be made omitting nx further additions

				gridvalue = mycellgrid[cxyz];

				//mutually exclusive
				/*
				if ( gridvalue == CELL_IS_A_PARTICLE ) { 
					rawdata[(cxyz*RGB)+RED] = UCHAR_RANGE_MIN; //RAW_PARTICLE they are black
					rawdata[(cxyz*RGB)+GREEN] = UCHAR_RANGE_MIN; //RAW_PARTICLE they are black
					rawdata[(cxyz*RGB)+BLUE] = UCHAR_RANGE_MIN; //RAW_PARTICLE they are black
					continue;
				}
				*/
				if ( gridvalue < mydefgpool.size() || gridvalue == CURRENTLY_INFECTED ) { 
					//rawdata[(cxyz*RGB)+RED] = UCHAR_RANGE_MAX; //RAW_INFECTED they are white
					//rawdata[(cxyz*RGB)+GREEN] = UCHAR_RANGE_MAX; //RAW_INFECTED they are white
					//rawdata[(cxyz*RGB)+BLUE] = UCHAR_RANGE_MAX; //RAW_INFECTED they are white

					rawdata[(cxyz*RGB)+RED] = UCHAR_RANGE_MIN; //RAW_INFECTED they are black
					rawdata[(cxyz*RGB)+GREEN] = UCHAR_RANGE_MIN; //RAW_INFECTED they are black
					rawdata[(cxyz*RGB)+BLUE] = UCHAR_RANGE_MIN; //RAW_INFECTED they are black

					continue;
				}
				//obviously >= mydefgpool.size() so completely recrystallized already
				//if ( gridvalue >= mydefgpool.size() ) {
					rxgid = gridvalue - ndg;
					oi = myrxgpool[rxgid].caori;

					rawdata[(cxyz*RGB)+RED] = myoripool[oi].RGBA[RED];
					rawdata[(cxyz*RGB)+GREEN] = myoripool[oi].RGBA[GREEN];
					rawdata[(cxyz*RGB)+BLUE] = myoripool[oi].RGBA[BLUE];
				//}
			} //along x
		} //xlines stacked in y
	} //xy layers

	//push in file at once, runtime environment will optimize, further potential by nowait multiple write from the thread, but tricky...
	//MPI_File_write_at(msFileHdl, totalOffset, rawdata, nxyz, MPI_INT, &msFileStatus);
	MPI_File_write(msFileHdl, rawdata, (RGB*nxyz), MPI_CHAR, &msFileStatus);

	delete [] rawdata;
	delete [] CmsFileName;

cout << myensHdl->myRank << " successful MPI I/O in timestep = " << this->step << endl;

	MPI_File_close(&msFileHdl); //no Barrier as MPI_COMM_SELF
}


#ifdef FRAGMENTATION_RADAR
	void caHdl::write_fragmentation_radar( void )
	{
		//render a xy section sliced at z = DEFAULT_ZSECTIONING
		uint32_t imgx = this->myrxfrontfragmentation[0].nblocks;	//status along x
		uint32_t imgy = this->myrxfrontfragmentation.size();		//time along y
		uint32_t imggridpoint;

		//allocate memory to store the resulting image
		unsigned char* img_rgba = new unsigned char [imgx*imgy*4];

		ostringstream fname;
		fname << "SCORE." << this->myensHdl->simid << ".JobID." << this->jobid << ".FragmentationRadar.png";

cout << "Rendering the fragmentation radar__" << fname.str().c_str() << "__ at " << imgx << ";" << imgy << endl;
		unsigned char* thevalues = NULL;
		for ( uint32_t y = 0; y < imgy; y++ ) {
			thevalues = myrxfrontfragmentation[y].rgba;

			for ( uint32_t x = 0; x < imgx; x++ ) {

				imggridpoint = (4 * y * imgx) + 4*x;

//cout << y << "\t\t" << x << "\t\t" << oid << "\t\t" << myoripool[oid].RGBA[REDCHAN] << ";"<< myoripool[oid].RGBA[GREENCHAN] << ";"<< myoripool[oid].RGBA[BLUECHAN] << ";"<< myoripool[oid].RGBA[ALPHACHAN] << endl;
				img_rgba[imggridpoint + REDCHAN ] = thevalues[4*x+REDCHAN];
				img_rgba[imggridpoint + GREENCHAN ] = thevalues[4*x+GREENCHAN];
				img_rgba[imggridpoint + BLUECHAN ] = thevalues[4*x+BLUECHAN];
				img_rgba[imggridpoint + ALPHACHAN ] = thevalues[4*x+ALPHACHAN];
			}
		}

		//encode the PNG file by calling the lodePNG library
		lodepng::encode( fname.str().c_str(), img_rgba , imgx, imgy );

		delete [] img_rgba;


		//construct scale bar
		//render a xy section sliced at z = DEFAULT_ZSECTIONING
		imgx = 100;
		imgy = 1000;
		imggridpoint;
		unsigned char* scale_rgba = new unsigned char [imgx*imgy*4];
		ostringstream fnamescale;
		fnamescale << "SCORE." << this->myensHdl->simid << ".JobID." << this->jobid << ".FragmentationRadarScaleBar.png";
		double fragmentation;
		unsigned char redvalue;
		for ( uint32_t y = 0; y < imgy; y++ ) {
			fragmentation = y / 10 * 0.01;
			redvalue = 255.0 * (1.0 - pow( fragmentation, REDSTRETCH ));

			for ( uint32_t x = 0; x < imgx; x++ ) {
				imggridpoint = (4 * y * imgx) + 4*x;
				scale_rgba[imggridpoint + REDCHAN ] = redvalue;
				scale_rgba[imggridpoint + GREENCHAN ] = 0;
				scale_rgba[imggridpoint + BLUECHAN ] = 0;
				scale_rgba[imggridpoint + ALPHACHAN ] = 255;
			}
		}
		//encode the scale bar as a PNG file by calling the lodePNG library
		lodepng::encode( fnamescale.str().c_str(), scale_rgba , imgx, imgy );
		delete [] scale_rgba;
	}
#endif


void caHdl::defragment_myRXFront( void )
{
	//MK::usually the function is kicked in shortly after the peak of the Sv(X) function
	//because the number of newly infected starts to decrease strongly while ACTIVE cells progressively become INACTIVE
	//always critical and not thread safe at the level of the caHdl
	//may only be call after resetting of recycling and FullRX list or prior to calcGrowthStep

	long which = 0;
	long moved = nextSlotNeverActiveRXFront - 1; //right-most/last valid entry accessible in myRXFront [0,nextSlotNeverActiveRXFront)
	//therewith automatically no defragmentation of a single entry or empty list

	//the idea is the following copy state of cells from high-cnt to low-cnt index
	while ( which < moved ) {

		if ( myRXFront[which].activity == INACTIVE ) {
			//slightly cache inefficient because which << moved but moved is not expected to be referencing start of memory page or cacheline
			myRXFront[which].activity = myRXFront[moved].activity;
			myRXFront[which].infector = myRXFront[moved].infector;
			myRXFront[which].ix = myRXFront[moved].ix;
			myRXFront[which].iy = myRXFront[moved].iy;
			myRXFront[which].iz = myRXFront[moved].iz;
			myRXFront[which].rxFrac = myRXFront[moved].rxFrac;
			myRXFront[which].P = myRXFront[moved].P;
			myRXFront[which].mydefgid = myRXFront[moved].mydefgid;
			myRXFront[which].myrxgid = myRXFront[moved].myrxgid;

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
	nextSlotNeverActiveRXFront = moved + 1;


	loginfo_defrag_cnt++;

	//cout << jobid << " I have defragmented myRXFront list - which;moved " << which << ";" << moved << endl;
	//from now on shorter overall if == ACTIVE comparisons in the main list
}


void caHdl::defragmentation( void ) 
{
	if ( Xcells >= defragRXFront_atthisX[loginfo_defrag_cnt] && loginfo_defrag_cnt < defragRXFront_atthisX.size() ) {
		defragment_myRXFront();
	}
}


void caHdl::log_OUTPUT( void )
{
	//based on current values of t,CurrentTemperature,X,Sv,step and so forth
	if ( outopt_localrxfront == OUTPUT_RXFRONTSTATS_YES ) {
		log_rxfrontstats();
	}


	//2D and 3D sections
	if ( renderingForThisCA == true ) {
		if ( Xcells >= rendering_atthisX[loginfo_rendering_cnt] && loginfo_rendering_cnt < rendering_atthisX.size() ) {

			if ( outopt_localrenderhow == RENDERING_MS2D ) {
				if ( outopt_localrendercolor == RENDERING_COLOR_GRAINID )
					write_zsection_coloring_grainids();

				if ( outopt_localrendercolor == RENDERING_COLOR_IPFZ )
					write_zsection_coloring_ipfz();
			}

			if ( outopt_localrenderhow == RENDERING_MS3D ) {
				if ( outopt_localrendercolor == RENDERING_COLOR_GRAINID )
					write_voxeldata_coloring_grainids(); 

				if ( outopt_localrendercolor == RENDERING_COLOR_IPFZ )
					write_voxeldata_coloring_ipfz();
			}

			loginfo_rendering_cnt++;
		}
	}


	//only if total RX fraction >= desired value
	if ( Xcells >= output_atthisX[loginfo_grainevo_cnt] && loginfo_grainevo_cnt < output_atthisX.size() ) {
		log_grainevo();
	}
}


void caHdl::solve_RXGROWTH( void )
{
	log_initialization();
	init_myRXFront();

	//count recrystallized volume in cells
	double boxvolume = myCAGeometry.nboxvol_rdndtd;

	t = 0.0;
	X = 0.0;
	Xcells = 0.0;
	dXstep = 0.0;
	Sv = 0;
	step = 0;

	//##MK::initial I/O prior to simulation? --> inject code here

	update_temperature();
	update_atomisticproperties();
	update_microchemistry();
	update_intrinsicmobilities();
	update_deformedsubstructure();

#ifdef REPORTSTYLE_USER
	cout << this->jobid << "\tstep;t;X;dXstep;dt;CurrTemperature;Sv;maxfillperstep\t\t" << this->step << "\t\t" << this->t << "\t\t" << this->X  << "\t\t" << this->dXstep << "\t\t" << this->dt << "\t\t" << this->CurrentTemperature << "\t\t" << this->Sv << "\t\t" << maxfillperstep << endl;
#endif

	//add nuclei into the simulation
	unsigned char incubationtimemodel = this->myNucleationModel.tincubmodel;

	if ( incubationtimemodel == TINCUB_TIMEDEPENDENT ) {
		update_system_fixedtime( 0.0 ); //##MK::currently only a dummy -->inject code here if desired
	}
	else { //SITE_SATURATION
		sim_myCA_sitesaturatedNucleation();
	}

	//get initial delta time for the numerical integration scheme
	dt = this->get_dtmax_minimum();

	do {
		dXstep = 0.0;

		//to assure the myRXFront list to be as compact filled with ACTIVE cells between [0, nextRX...) as possible
		defragmentation();

		//##MK::inject a optional time dependent nucleation model here
		if ( incubationtimemodel == TINCUB_TIMEDEPENDENT ) {
			sim_myCA_timedependentNucleation();
		}

		//fill ACTIVE cells while keep track of INACTIVE cells in myRecyclingList
		sim_myCA_calcGrowthStep();

		//let ACTIVE cells now infect new cells
		sim_myCA_updateFullRX();

		log_OUTPUT();

		//--> up to here it is possible to account always a new for the myMobMax in order to stretch the time towards the end of the simulation in a 
		//scenario where only slowly migrating boundaries remain and the estimate of the maximum Pvalue that was EVER AT SOME POINT obtained is too conservative
		//prepare for next timestep, same integration order than START1D/3DCA and COReV3
		t += dt; 

		update_temperature();
		update_atomisticproperties();
		update_microchemistry();
		update_intrinsicmobilities();
		update_deformedsubstructure();

		//updateFullRX was the last one to determined Pmax in this timestep, now consider potentially higher mobility to refine time stepping
		dt = this->get_dtmax_minimum();

		step++;

#ifdef REPORTSTYLE_USER
		cout << this->jobid << "\tstep;t;X;dXstep;dt;CurrTemperature;Sv;maxfillperstep\t\t" << this->step << "\t\t" << this->t << "\t\t" << this->X  << "\t\t" << this->dXstep << "\t\t" << this->dt << "\t\t" << this->CurrentTemperature << "\t\t" << this->Sv << "\t\t" << maxfillperstep << endl;
#endif

	} //go out when completely recrystallized, recrystallized to use defined value, when too many integration steps, when processing out
	while ( (Xcells < myCAGeometry.nboxvol_rdndtd) && (X < XMAX) && (step < this->NTSTEPSMAX) && (t < tprocessingend) && (t < TMAX) );
	//no MPI_Barriers whatsoever necessary...

	this->tsimend = t;
	this->stepsimend = step;

	//##MK::optional final output

	//from now on the memory used to voxelize the RX microstructure and manage cells is no longer necessary
	this->cleanupMemoryUsedForGrowthMachine();


	cout << myensHdl->myRank << " myRank,jobid = " << this->jobid << " local CA simulation was successfully executed, tsimend = " << tsimend << ", stepsimend = " << stepsimend << endl;
}



void caHdl::log_ca_physics( loginfo_ca_physicsP  container )
{
	container->nx = this->myCAGeometry.nboxedge_rd;
	container->ny = this->myCAGeometry.nboxedge_nd;
	container->nz = this->myCAGeometry.nboxedge_td;
	container->nxyz = this->myCAGeometry.nboxvol_rdndtd;

	container->nboundarycells = this->SvDeformed;
	container->ndgrseeds = this->ndefgseeds;
	container->nrxgrseeds = this->nmyrxgpool;

	container->ndefmicrotexture = this->DefMicrotextureClasses;
	container->nnucleimicrotexture = this->NucleiMicrotextureClasses;
	container->defmicrotexture = this->DefMicrotexture;
	container->nucleimicrotexture = this->NucleiMicrotexture;
	container->storedenergy = this->StoredEnergy;
}


double caHdl::get_interpCellCount( uint32_t localid, double when )
{
	//assumes correct localid
	//##MK::mygrainevolution has to contain information taken at t = 0.0 when X=0.0 and t = tsimend when X = XMAX
	double interpCellCount = 0.0;

	if ( when >= this->tsimend ) { //MK::when RX was already complete local I the take final value stored in mydefgpool

		if ( localid < mydefgpool.size() ) { 
			interpCellCount = mydefgpool[localid].cellcount;
			return interpCellCount;
		}
		//localid >= mydefgpool.size() -> it is an id of a nucleus
		interpCellCount = myrxgpool[(localid - mydefgpool.size())].cellcount;
		return interpCellCount;
	}

	//scan from mygrainevolution
	uint32_t i = 0;

	uint32_t nmygrainevo = mygrainevolution.size();
	while ( mygrainevolution[i].localtime < when && i < nmygrainevo ) { //framing the interval ([i-1].localtime, [i].localtime)
		i++;
	}

	//MK::implicitly assuming the microstructure unchanged
	//MK::lift load from branch predictor by planning for the most likely case
	if ( i > 0 && i < nmygrainevo ) {
		//###disjunctness of rediscretization scheme has been assured by ensRediscrWindow
		//only possible case left if (i < nmygrainevo ) {
		//MK::here almost for sure cache-thrashing occuring as [i-1] and [i] are pointers to disjoint heap regions and localid on top of this is usually > L1 size
		double ti1 = mygrainevolution[i-1].localtime;
		double cnti1 = mygrainevolution[i-1].localdatabucket[localid];

		double ti = mygrainevolution[i].localtime;
		double cnti = mygrainevolution[i].localdatabucket[localid];

		interpCellCount = cnti1 + ((when - ti1) * ((cnti - cnti1) / (ti - ti1)));

		//cout << "----->when,i,cnti,cnti1,ti,ti1,deltaccnt;intpcnt\t" << when << "\t" << i << "\t" << cnti << "\t" << cnti1 << "\t" << ti << "\t" << ti1 << "\t" << interpCellCount << endl;
		return interpCellCount;
	}

	//test the seldom scenarios
	if ( i == 0 ) {
		interpCellCount = mygrainevolution[i].localdatabucket[localid];
		return interpCellCount;
	}

	//MK::implicitly assuming coarsening via grain growth to be negligible!
	if ( i >= nmygrainevo ) {
		interpCellCount = mygrainevolution[nmygrainevo-1].localdatabucket[localid];
		return interpCellCount; //necessary because mygrainevolution[i=nmygrainevo] would seg...
	}
}



void caHdl::write_rxfrontstats( void )
{
	//I/O interaction
	stringstream log_rxfrontstats_fname;
	ofstream log_rxfrontstats_file;
	
	log_rxfrontstats_fname << "SCORE." << myensHdl->simid << ".MyRank." << myensHdl->myRank << ".JobID." << this->jobid << ".RXFrontStats.csv";
#ifdef REPORTSTYLE_DEVELOPER
	cout << "File " << log_rxfrontstats_fname.str().c_str() << " is opened now" << endl;
#endif
	log_rxfrontstats_file.open ( log_rxfrontstats_fname.str().c_str() );

	//header
	log_rxfrontstats_file << "Step;Time/s;X;Sv;myMemory/Byte;currPmax;ntotalRXFront;nextSlotNeverActiveRXFront;nCurrentlyActive;ntotalFullRXList;nextSlotToFullRX;ntotalRecyclingList;nextSlotThatBecomesRecycled;firstNotRecycledYet" << endl;

	for ( uint32_t s = 0; s < myrxfrontstatus.size(); ++s ) {
		log_rxfrontstats_file << myrxfrontstatus[s].localstep << ";" << myrxfrontstatus[s].localtime << ";" << myrxfrontstatus[s].localX << ";" << myrxfrontstatus[s].localSv << ";" << myrxfrontstatus[s].localmemory << ";" << myrxfrontstatus[s].localPmax << ";" << myrxfrontstatus[s].ntotalRXFront << ";" << myrxfrontstatus[s].nextSlotNeverActiveRXFront << ";" << myrxfrontstatus[s].nCurrentlyActive << ";" << myrxfrontstatus[s].ntotalFullRXList << ";" << myrxfrontstatus[s].nextSlotToFullRX << ";" << myrxfrontstatus[s].ntotalRecyclingList << ";" << myrxfrontstatus[s].nextSlotThatBecomesRecycled << ";" << myrxfrontstatus[s].firstNotRecycledYet << endl;
	}

	log_rxfrontstats_file.flush();
	log_rxfrontstats_file.close();
}


void caHdl::write_grainevolution( void )
{
	stringstream log_grainevo_fname;
	ofstream log_grainevo_file;
	
	log_grainevo_fname << "SCORE." << myensHdl->simid << ".MyRank." << myensHdl->myRank << ".JobID." << this->jobid << ".SingleGrainData.csv";
#ifdef REPORTSTYLE_DEVELOPER
	cout << "File " << log_grainevo_fname.str().c_str() << " is opened now" << endl;
#endif
	log_grainevo_file.open ( log_grainevo_fname.str().c_str() );

	//header
	log_grainevo_file << "GrainID/Step;;;;;;;;";
	for ( uint32_t s = 0; s < loginfo_grainevo_cnt; s++ ) { log_grainevo_file << ";" << output_atthisX[s]; }
	log_grainevo_file << endl;
	log_grainevo_file << "Step;;;;;;;;";
	for ( uint32_t s = 0; s < loginfo_grainevo_cnt; s++ ) { log_grainevo_file << ";" << mygrainevolution[s].localtstep; }
	log_grainevo_file << endl;
	log_grainevo_file << "Time;;;;;;;;";
	for ( uint32_t s = 0; s < loginfo_grainevo_cnt; s++ ) { log_grainevo_file << ";" << mygrainevolution[s].localtime; }
	log_grainevo_file << endl;
	log_grainevo_file << "X;;;;;;;;";
	for ( uint32_t s = 0; s < loginfo_grainevo_cnt; s++ ) { log_grainevo_file << ";" << mygrainevolution[s].localX; }
	log_grainevo_file << endl;
	log_grainevo_file << "Sv;;;;phi1;PHI;phi2;IdealOri;";
	for ( uint32_t s = 0; s < loginfo_grainevo_cnt; s++ ) { log_grainevo_file << ";" << mygrainevolution[s].localSv; }
	log_grainevo_file << endl;

	//each grain one row, each step a column, this matrix transpose is practical because most often the number of logsteps << number of grains
	//however this is indeed a performance hit
	uint32_t ndisjoint_defgrains = mydefgpool.size();
	uint32_t ndisjoint_grains = ndisjoint_defgrains + myrxgpool.size();
	for ( uint32_t dgr = 0; dgr < ndisjoint_defgrains; dgr++ ) {
		uint32_t oriid = mydefgpool[dgr].caori;

		log_grainevo_file << dgr << ";;;;;" << myoripool[oriid].bunge1 << ";" << myoripool[oriid].bunge2 << ";" << myoripool[oriid].bunge3 << ";" << myoripool[oriid].closestideal; 

		//yeah, it hits caches badly I know, as matrix transposes do ...
		for ( uint32_t s = 0; s < loginfo_grainevo_cnt; s++) {
			log_grainevo_file << ";" << this->mygrainevolution[s].localdatabucket[dgr];
		}
		log_grainevo_file << endl;
	}

	log_grainevo_file << ";ix;iy;iz;boundarycontact;phi1;PHI;phi2;IdealOri";
	for ( uint32_t s = 0; s < loginfo_grainevo_cnt; s++) { log_grainevo_file << ";"; }
	log_grainevo_file << endl;

	for ( uint32_t rgr = ndisjoint_defgrains; rgr < ndisjoint_grains; rgr++ ) {
		uint32_t i = rgr - ndisjoint_defgrains;

		uint32_t z = myrxgpool[i].startingsite / myCAGeometry.nboxarea_rdnd;
		uint32_t rem = myrxgpool[i].startingsite - (myCAGeometry.nboxarea_rdnd * z);
		uint32_t y = rem / myCAGeometry.nboxedge_rd;
		uint32_t x = rem - (myCAGeometry.nboxedge_rd * y);

		uint32_t oriid = myrxgpool[i].caori;

		log_grainevo_file << i << ";" << x << ";" << y << ";" << z << ";" << myrxgpoolboundarycontact[i] << ";" << myoripool[oriid].bunge1 << ";" << myoripool[oriid].bunge2 << ";" << myoripool[oriid].bunge3 << ";" << myoripool[oriid].closestideal; //optionally add << "-RXG";
		
		//yeah, it hits caches badly I know, as matrix transposes do ...
		for ( uint32_t s = 0; s < loginfo_grainevo_cnt; s++) {
			log_grainevo_file << ";" << mygrainevolution[s].localdatabucket[rgr];
		}
		log_grainevo_file << endl;
	}

	log_grainevo_file.flush();
	log_grainevo_file.close();
}


void caHdl::write_hpf3d( void )
{
	//HPF in plain text
	stringstream log_hpfvalues3d_fname;
	ofstream log_hpfvalues3d_file;
	
	log_hpfvalues3d_fname << "SCORE." << myensHdl->simid << ".MyRank." << myensHdl->myRank << ".JobID." << this->jobid << ".HPF3D.csv";
#ifdef REPORTSTYLE_DEVELOPER
	cout << "File " << log_hpfvalues3d_fname.str().c_str() << " is opened now" << endl;
#endif
	log_hpfvalues3d_file.open ( log_hpfvalues3d_fname.str().c_str() );

	//header
	log_hpfvalues3d_file << "wx;wy;wz;counts" << endl;

	for (int wz = - myhpfwindow.dz; wz <= myhpfwindow.dz; wz++) {
		int hz = wz + myhpfwindow.dz;

		for (int wy = - myhpfwindow.dy; wy <= myhpfwindow.dy; wy++) {
			int hy = wy + myhpfwindow.dy;

			for (int wx = - myhpfwindow.dx; wx <= myhpfwindow.dx; wx++) {
				int hx = wx + myhpfwindow.dx;
				int hxyz = hx + (hy * myhpfwindow.nx) + (hz * myhpfwindow.nxy);

				log_hpfvalues3d_file << wx << ";" << wy << ";" << wz << ";" << myhpfvalues[hxyz] << endl;
			}
		}
	}


	log_hpfvalues3d_file.flush();
	log_hpfvalues3d_file.close();
}


void caHdl::write_hpf2d( int wz )
{
	//HPF in plain text
	stringstream log_hpfvalues2d_fname;
	ofstream log_hpfvalues2d_file;
	
	log_hpfvalues2d_fname << "SCORE." << myensHdl->simid << ".MyRank." << myensHdl->myRank << ".JobID." << this->jobid << ".HPF2DzPos." << wz << ".csv";
#ifdef REPORTSTYLE_DEVELOPER
	cout << "File " << log_hpfvalues2d_fname.str().c_str() << " is opened now" << endl;
#endif
	log_hpfvalues2d_file.open ( log_hpfvalues2d_fname.str().c_str() );

	//header
	log_hpfvalues2d_file << "wx;wy;counts" << endl;


	int hpfx = myhpfwindow.dx;
	int hpfy = myhpfwindow.dy;
	int hpfz = myhpfwindow.dz;
	int hpfxy = hpfx * hpfy;

	int hx, hy, hz, hyzoff, hxyz;

	hz = wz + hpfz;

	for (int wy = - hpfy; wy <= hpfy; wy++) {
		hy = wy + hpfy;

		hyzoff = (hy * hpfx) + (hz * hpfxy);
		for (int wx = - hpfx; wx <= hpfx; wx++) {
			hx = wx + hpfx;
			hxyz = hx + hyzoff;

			log_hpfvalues2d_file << wx << ";" << wy << ";" << myhpfvalues[hxyz] << endl;
		}
	}


	log_hpfvalues2d_file.flush();
	log_hpfvalues2d_file.close();
}


#define NUCSCALING_HAASE	(1.0e6)

long caHdl::gbfacenucnumbermodel_haase( double facearea, double rhodifference, double disori )
{
	//no nucleation at low-angle grain boundaries desired? then no nuclei!
	if ( myNucleationModel.gbnucleation == GB_ONLYATHAGB && disori < LAGB_TO_HAGB_TRANS )
		return NO_NUCLEI;

	//too low a difference? no nuclei here!
	if ( rhodifference <= MINIMUM_DRHO )
		return NO_NUCLEI;


	//implements the physics of the nucleation at grain boundary faces
	long gbnuc = NO_NUCLEI;

	//number density is scaled according to absolute dislocation density difference * empirical efficiency

	double lamb = rhodifference * myNucleationModel.gbnucleation_drho2dens; //rhodifference/rhodiffmax*SCALING

	double nnuc = lamb * myNucleationModel.gbnucleation_dens2num * facearea;


//##DEBUGcout << "face/rhodiff/lamb/nnuc = " << facearea << "\t\t" << rhodifference << "\t\t" << lamb << "\t\t" << nnuc << endl;

	//cast down
	gbnuc = (long) nnuc; //truncates trailing decimals

	//therefore, are nuclei or lost or not?
	nnuc = nnuc - (double) gbnuc;

	//add another one
	if ( nnuc >= localprng.leEcuyer() ) 
		gbnuc++;

	return gbnuc;
}



long caHdl::gbfacenucnumbermodel_phdmk( double facearea, double rhodifference, double disori )
{
	//no nucleation at low-angle grain boundaries desired? then no nuclei!
	if ( myNucleationModel.gbnucleation == GB_ONLYATHAGB && disori < LAGB_TO_HAGB_TRANS )
		return NO_NUCLEI;

	//too low a difference? no nuclei here!
	if ( rhodifference <= MINIMUM_DRHO )
		return NO_NUCLEI;


	//implements the physics of the nucleation at grain boundary faces
	long gbnuc = NO_NUCLEI;

	//number density is scaled according to Poisson-distribution the probability of seeding more nuclei the higher the stronger and the large drho and the large the boundary area
	//larger and higher the dislocation density difference

	double lamb = rhodifference * myNucleationModel.gbnucleation_drho2dens; //rhodifference/rhodiffmax*SCALING

//cout << "gbfacenuc\t\t" << dg1 << ";" << dg2 << "\t\t" << disori << "\t\t" << lamb << endl;

	double pcum[POISSON_CUMSUM_TABLE];
	double pkcumsum = 0.0;
	for ( long k = 0; k <= POISSON_CUMSUM_CUTOFF; k++ ) {
		pcum[k] = poissondistribution( lamb, k );

//cout << "\t\t---" << k << "\t\t" << pcum[k] << "\t\t" << pow(lamb, k) << "\t\t" << fac(k) << "\t\t" << pkcumsum << endl;

		pkcumsum = pkcumsum + pcum[k];
	}

	if ( pkcumsum < DOUBLE_ACCURACY ) { cout << "ERROR::pkcumsum is unexpectedly close to zero " << pkcumsum << endl; return NO_NUCLEI; }

	//form cumulated distribution
	double pk = 0.0;
	for ( long k = 0; k <= POISSON_CUMSUM_CUTOFF; k++) {
		pk = pk + pcum[k];
		pcum[k] = pk / pkcumsum;

//cout << "\t\t" << k << "\t\t" << pcum[k] << endl;
	}

	//pick random number
	double luckynumber = localprng.leEcuyer();

	long j = 0; //which number to take?
	while ( (luckynumber > pcum[j]) && (j < POISSON_CUMSUM_CUTOFF) ) {
		j++;
	}

//cout << luckynumber << "\t\t" << j << "\t\t" << pcum[j] << "\t\t" << SCALE_NUCDENSITY << "\t\t" << bndarea_micron << endl;

	double nnuc = j;
	nnuc = nnuc * myNucleationModel.gbnucleation_dens2num * facearea;

	//cast down
	gbnuc = (long) nnuc; //truncates trailing decimals

	//therefore, are nuclei or lost or not?
	nnuc = nnuc - (double) gbnuc;

	//add another one
	if ( nnuc >= localprng.leEcuyer() ) 
		gbnuc++;

	return gbnuc;
}


void caHdl::solve_nucmodeling_gbnuc_physicalmodelFast( void )
{
	//pack the nuclei locally into the myrxgpool vector analysis contributions one boundary after another
	long N = 0;
	uint32_t activeBoundaryCells = 0;
	uint32_t nworldrxgpool = this->myensHdl->worldrxgpool.size();


	for ( uint32_t b = 0; b < boundaryBucketSizeFast; b++ ) {
		//no boundary here
		if ( boundaryBucketFast[b].len == 0) continue;

		for ( uint32_t f = 0; f < boundaryBucketFast[b].len; f++ ) {
			bndFaceFast* aface = &(boundaryBucketFast[b].thefaces[f]);

			if ( aface->nvoxeltotal <= EMPTY ) continue;

			//what is the deformed grain with the higher dislocation density?
			uint32_t seedmax = aface->gposUp;
			uint32_t seedmin = aface->gposDown;
			//lets assume and reject that posUp has higher dislocation density

			double rhoUp = mydefgpool[tmpdefgseeds[seedmax].mydefgpoolid].rho0;
			double rhoDown = mydefgpool[tmpdefgseeds[seedmin].mydefgpoolid].rho0;

			if( rhoDown > rhoUp ) {
				seedmax = aface->gposDown;
				seedmin = aface->gposUp;
			} //okay but now seedmax has the highest density

			double drho = fabs(rhoUp - rhoDown);

			//the reservoir is the aface->voxelizedface
			activeBoundaryCells = aface->nvoxeltotal; //not nextfreevoxelslot because aface-> was trimmed by trimGB...

			N = gbfacenucnumbermodel_haase( activeBoundaryCells, drho, aface->disori );
			//##MK::PhdMK
			//N = gbfacenucnumbermodel_phdmk( activeBoundaryCells, drho, aface->disori );

			//##MK::QUICKHACKFIXEDNUMBERN = 1;

			//##DEBUG
			if ( N >= activeBoundaryCells ) { 
				cout << "WARNING::" << myensHdl->myRank << ";" << this->jobid << "\tplanned for N = " << N << " nuclei but only " << activeBoundaryCells << " available sites in the boundary so I am limiting N to (activeBoundaryCells-1)\t\tdrho;disori;" << drho << ";" << aface->disori << endl;
				N = activeBoundaryCells - 1;
			}

			//begin placement by picking randomly from the reservoir the sites but do not repeat not more than REPEAT_MAX_DRAWING
			uint32_t repeat = 0;
			uint32_t randomCell, nucsite;

//cout << "gbnucphysical\t\tb;f;N == " << b << ";" << f << ";" << N << endl;

			//MK::Orientation of the nucleus is the same as that of the grain with lower dislocation density, no scatter
			uint32_t oid = mydefgpool[tmpdefgseeds[seedmin].mydefgpoolid].caori;
			double gbnuc_dgav = this->myNucleationModel.gbnucleation_scatter;
			double seedmin_dgav0 = mydefgpool[tmpdefgseeds[seedmin].mydefgpoolid].avdg0;
			double seedminOri[3];
			seedminOri[0] = myoripool[oid].bunge1;
			seedminOri[1] = myoripool[oid].bunge2;
			seedminOri[2] = myoripool[oid].bunge3;


			while( N > 0 ) {
				repeat ++;

				if ( repeat >= REPEAT_MAX_DRAWING ) {
					repeat = 0;
					N--; 
					continue;
				}

				randomCell = (localprng.leEcuyer() * activeBoundaryCells); //select one cell from the grain boundaries

				//##DEBUG
				QUICKASSERT ( randomCell < activeBoundaryCells );

				//when the place is already evaluated or not located in a grain with higher energy chose a new site
				if ( aface->voxelizedface[randomCell].location == EMPTY || aface->voxelizedface[randomCell].seedid != seedmax ) { continue; }

				nucsite = aface->voxelizedface[randomCell].location;
				aface->voxelizedface[randomCell].location = EMPTY; //not location can be empty because even uint32_t location 0 which is voxel (ix,iy,iz = 0,0,0) has been assigned MARKER_.. + 0

				//inheriting the orientation from the seedmin


				//now with an ##MK uncorrelated scatter about the orientation with the lower dislocation density - SIBM
				double nucleusOri[3];
				//##MK::Haase randomly from the list
				uint32_t luckyNucleus = localprng.leEcuyer() * nworldrxgpool;

				uint32_t ensoriid = myensHdl->worldrxgpool[luckyNucleus].ori;

				nucleusOri[0] = myensHdl->worldoripool[ensoriid].bunge1;
				nucleusOri[1] = myensHdl->worldoripool[ensoriid].bunge2;
				nucleusOri[2] = myensHdl->worldoripool[ensoriid].bunge3;


				//##MK::giving random disorientation distribution 
				//newOrientationFromReference( seedminOri, gbnuc_dgav, nucleusOri );

				//##MK::giving close to Rayleigh distribution with sigma = seedmin_dgav0 distributed orientation spectrum
				specificallyDisorientednewOriFromReference( seedminOri, seedmin_dgav0, nucleusOri );

				uint32_t noid = ca_check_disjunctness( nucleusOri );

unsigned char rrgb[3]; double pposs[2]; bunge2ipfz( nucleusOri[0], nucleusOri[1], nucleusOri[2], rrgb, pposs );
cout << "gbnucphysical\t\tbunge123;scatter;newOri;DisoriBunge123ToNewOri\t\t" << seedminOri[0] << ";" << seedminOri[1] << ";" << seedminOri[2] << "\t\t" << myNucleationModel.gbnucleation_scatter << "\t\t" << nucleusOri[0] << ";" << nucleusOri[1] << ";" << nucleusOri[2] << "\t\t;" << this->misorientationCubic( seedminOri[0], seedminOri[1], seedminOri[2], nucleusOri[0], nucleusOri[1], nucleusOri[2] ) << "\t\t;" << pposs[0] << ";" << pposs[1] << endl;

				//##MK::or implement other models or more complex models or twinning models if there are better experimental data justifying there usage

				//recrystallized nuclei are always unique in the local automaton, other automata are not interested in local nuclei
				struct carxg newNucleus;

				newNucleus.caori = noid;
				newNucleus.cellcount = 0;
				newNucleus.nucsite = nucsite;
				newNucleus.startingsite = nucsite;

				//MK::incubation time model is site-saturated at the moment
				newNucleus.tincub = 0.0;

				myrxgpool.push_back ( newNucleus );

#ifdef REPORTSTYLE_DEVELOPER
				cout << "myRank;jobid;nucid;oid;ix;iy;iz;nucsite;cellcnt;tincub\t\t" << this->myensHdl->myRank << ";" << this->jobid << ";" << (myrxgpool.size() - 1) << ";" << oid << ";" << cb->ix << ";" << cb->iy << ";" << cb->iz << ";" << nucsite << ";" << myrxgpool[myrxgpool.size()-1].cellcount << ";" << myrxgpool[myrxgpool.size()-1].tincub << endl;
#endif

				N--;
				repeat = 0;
			}
		} //next face
	} //next face collection
}


void caHdl::characterize_nucsite_environment( void )
{
	//get memory for that correlation function
	myhpfvalues = NULL;
	myhpfvalues = new int[myhpfwindow.nxyz];
	QUICKASSERT ( myhpfvalues != NULL );

	//initialize as 0
	for ( int hpfxyz = 0; hpfxyz < myhpfwindow.nxyz; hpfxyz++) { 
		myhpfvalues[hpfxyz] = 0; 
	}

	//MK::HPF estimation window can not be larger than domain itself!
	uint32_t nbx = myCAGeometry.nboxedge_rd;
	uint32_t nby = myCAGeometry.nboxedge_nd;
	uint32_t nbz = myCAGeometry.nboxedge_td;
	uint32_t nbxy = nbx * nby;

	int hpfx = myhpfwindow.dx;
	int hpfy = myhpfwindow.dy;
	int hpfz = myhpfwindow.dz;
	int hpfxy = hpfx * hpfy;
	int hpfxyz;

	QUICKASSERT ( myhpfwindow.nx < nbx ); //the window is smaller as the whole domain
	QUICKASSERT ( myhpfwindow.ny < nby );
	QUICKASSERT ( myhpfwindow.nz < nbz );

	//analyze all places for all nuclei in the observation window 
	//MK::if window passes domain boundaries, apply periodic boundary conditions
cout << "Characterizing the HPF function ..." << endl;

	uint32_t nrxg = myrxgpool.size();

	for ( uint32_t nuc = 0; nuc < nrxg; nuc++) {
		uint32_t nucleateWhere = myrxgpool[nuc].nucsite;

		uint32_t dgid = mycellgrid[nucleateWhere];

		//implicit to explicit coordinate transformation, int is possible because the automaton is at least SQR(CA_DIMENSIONS_MINIMUM)
		int iz = nucleateWhere / myCAGeometry.nboxarea_rdnd;
		uint32_t rem = nucleateWhere - (nbxy * iz); //uint32_t necessary because iz can be 0 and then rem = nucleateWhere which can exceed the range of int!
		int iy = rem / nbx;
		int ix = rem - (nbx * iy);

		//analyze window for that nucleus
		int hz, hy, hx, hpfyzoff;
		int x, y, z; //can become negative
		uint32_t px, py, pz; //MK::always positive
		uint32_t yzoff;
		uint32_t xyz;

		for (int wz = - hpfz; wz <= hpfz; wz++) {
			z = iz + wz;
			hz = wz + hpfz;
			//some position maybe also outside the automaton

			if ( z < 0 )
				z += nbz;
			if ( z >= nbz )
				z -= nbz;

			//now definately inside so between [CA_DIMENSIONS_MINIMUM,CA_DIMENSIONS_MAXIMUM]
			//now z is positive
			pz = (uint32_t) z;

			for (int wy = - hpfy; wy <= hpfy; wy++) {
				y = iy + wy;
				hy = wy + hpfy;

				if ( y < 0 )
					y += nby;
				if ( y >= nby )
					y -= nby;

				//now y is positive
				py = (uint32_t) y;

				yzoff = (py * nbx) + (pz * nbxy);
				hpfyzoff = (hy * hpfx) + (hz * hpfxy);

				for (int wx = - hpfx; wx <= hpfx; wx++) {
					x = ix + wx;
					hx = wx + hpfx;

					//apply periodic boundary conditions
					if ( x < 0 )
						x += nbx;
					if ( x >= nbx )
						x -= nbx;

					//now x is positive
					px = (uint32_t) x;
					xyz = px + yzoff;

					//scan mycellgrid if at (x,y,z) the id is exactly dgid if so HPF function is 1 otherwise 0
					hpfxyz = hx + hpfyzoff;

					if ( mycellgrid[xyz] == dgid ) {
						myhpfvalues[hpfxyz] += 1;
					}
				} //along positive x parallel RD
			} //for xlines stacked y
		} //for xy slabs stacked in z

	} //next nucleus

cout << "Having characterized the HPF function." << endl;
}


//##MK::boundaries fast
void caHdl::initializeBoundariesFast( void )
{
	//an array pointing to a collection of faces per disjoint seed
	boundaryBucketSizeFast = ndefgseeds;

	boundaryBucketFast = NULL;
	boundaryBucketFast = new bndColumnFast[boundaryBucketSizeFast];
	QUICKASSERT ( boundaryBucketFast != NULL );
	myMemGuard = myMemGuard + boundaryBucketSizeFast * sizeof(bndColumnFast);

	for( uint32_t i = 0; i < boundaryBucketSizeFast; i++ ) {
		boundaryBucketFast[i].len = 0;
		boundaryBucketFast[i].maxLen = 0;
		boundaryBucketFast[i].thefaces = NULL;
	}
}


void caHdl::emptyBoundariesFast( void )
{
	double bmemfreed = 0.0;
	double szofcbnd = sizeof(cellsBndFast);
	double szofface = sizeof(bndFaceFast);

	for ( uint32_t b = 0; b < boundaryBucketSizeFast; b++ ) {

		if ( boundaryBucketFast[b].len == 0 ) continue;

		for ( uint32_t f = 0; f < boundaryBucketFast[b].maxLen; f++ ) {
			delete [] boundaryBucketFast[b].thefaces[f].voxelizedface;
			boundaryBucketFast[b].thefaces[f].voxelizedface = NULL;
			bmemfreed = bmemfreed + (boundaryBucketFast[b].thefaces[f].nextfreevoxelslot * szofcbnd);
		}

		delete [] boundaryBucketFast[b].thefaces;
		boundaryBucketFast[b].thefaces = NULL;
		bmemfreed = bmemfreed + (szofface * boundaryBucketFast[b].maxLen);
		boundaryBucketFast[b].len = 0;
		boundaryBucketFast[b].maxLen = 0;
	}

	delete [] boundaryBucketFast;
	boundaryBucketFast = NULL;
	boundaryBucketSizeFast = 0;

	myMemGuard = myMemGuard - bmemfreed - (this->boundaryBucketSizeFast * sizeof( bndColumnFast));
}


void caHdl::addVoxelAtBoundary( uint32_t seed0, uint32_t seedtt, uint32_t loc )
{
	//positions upsite and dwsite were identified part of two seedgrains with disjoint and positive ids seedup and seeddown
	
	//a unique hash is necessary that identifies to which boundary with which two grains the voxel is adjacent
	long pUp = seed0;
	long pDown = seedtt;

	long posmax = MAX( pUp, pDown );

	//is this the first voxelpair linked into a face identified by grain MAX(pUp, pDown)?
	if( boundaryBucketFast[posmax].maxLen == 0 ) //get memory
	{
		boundaryBucketFast[posmax].len = 0;
		boundaryBucketFast[posmax].maxLen = STDLEN;
		boundaryBucketFast[posmax].thefaces = new bndFaceFast[STDLEN];
		myMemGuard = myMemGuard + ( STDLEN * sizeof(bndFaceFast) );

		//initialize for safety
		bndFaceFastP fc = boundaryBucketFast[posmax].thefaces;
		for ( uint32_t af = 0; af < boundaryBucketFast[posmax].maxLen; af++ ) {
			fc[af].voxelizedface = NULL;
			fc[af].nvoxeltotal = EMPTY;
			fc[af].nextfreevoxelslot = THEFIRSTONE;
		}
	}

	//link in directly this voxelpair at the appropriate boundary
	uint_fast64_t hashtag = MYHASH(seed0,seedtt);

	//populate the faces at posmax

	uint32_t i;
	//if bin has just been created i is headcontrolled set to 0 but for loop not evaluated
	//find an already existing grain boundary face, utilize that MK::ids are always positive and disjoint!
	for ( i = 0; i < boundaryBucketFast[posmax].len; i++ )
	{
		if( boundaryBucketFast[posmax].thefaces[i].id == hashtag ) //if the current oriIndex is the one in the oribucket, it is increased
		{
			//enough memory?
			bndFaceFastP aface = &(boundaryBucketFast[posmax].thefaces[i]);
			//still space in the container most likely case
			if ( aface->nextfreevoxelslot < aface->nvoxeltotal ) {
				aface->voxelizedface[aface->nextfreevoxelslot].location = MARKER_TO_IDENTIFY_EXPELLED_SITES + loc;
				aface->voxelizedface[aface->nextfreevoxelslot].seedid = seed0;
				aface->nextfreevoxelslot++;
				return;
			}

			//container too small reallocate and copy
			size_t oldSize = aface->nvoxeltotal;
			size_t newSize = aface->nvoxeltotal + DEFAULT_VOXEL2FACE_REFERENCES;
			QUICKASSERT ( newSize < CA_ALLOCATION_MAXIMUM );

			cellsBndFastP newblock = NULL;
			newblock = new cellsBndFast[newSize];
			QUICKASSERT( newblock != NULL );
			myMemGuard = myMemGuard + ( newSize * sizeof(cellsBndFast) );

			//copy the old voxelpairs first
			cellsBndFastP oldblock = aface->voxelizedface;

			for ( uint32_t vxp = 0; vxp < oldSize; vxp++ ) {
				newblock[vxp].location = oldblock[vxp].location;
				newblock[vxp].seedid = oldblock[vxp].seedid;
			}

			//fill new site into newblock
			newblock[oldSize].location = MARKER_TO_IDENTIFY_EXPELLED_SITES + loc;
			newblock[oldSize].seedid = seed0;

			//free the old block
			delete [] boundaryBucketFast[posmax].thefaces[i].voxelizedface;
			myMemGuard = myMemGuard - (oldSize * sizeof(cellsBndFast) );

			//link the new block
			boundaryBucketFast[posmax].thefaces[i].voxelizedface = newblock;
			boundaryBucketFast[posmax].thefaces[i].nvoxeltotal = newSize;
			boundaryBucketFast[posmax].thefaces[i].nextfreevoxelslot++;

			return;
		}
	}

	//mhh we are initializing a completely new face
	if( i < boundaryBucketFast[posmax].maxLen )
	{

		boundaryBucketFast[posmax].thefaces[i].gposUp = MAX(seed0, seedtt);
		//QUICKASSERT ( boundaryBucketFast[posmax].thefaces[i].gposUp == pUp ); //##MK::because risky, initially promoted int to long downgraded again to unit_fast32t...

		boundaryBucketFast[posmax].thefaces[i].gposDown = MIN(seed0, seedtt); 
		//QUICKASSERT ( boundaryBucketFast[posmax].thefaces[i].gposDown == pDown );

		boundaryBucketFast[posmax].thefaces[i].id = hashtag;
		boundaryBucketFast[posmax].thefaces[i].disori = GBDISORI_NOT_REQUIRED_YET;

		boundaryBucketFast[posmax].thefaces[i].voxelizedface = NULL;
		boundaryBucketFast[posmax].thefaces[i].voxelizedface = new cellsBndFast[DEFAULT_VOXEL2FACE_REFERENCES];
		QUICKASSERT( boundaryBucketFast[posmax].thefaces[i].voxelizedface != NULL );
		myMemGuard = myMemGuard + (DEFAULT_VOXEL2FACE_REFERENCES * sizeof(cellsBndFast));

		boundaryBucketFast[posmax].thefaces[i].nvoxeltotal = DEFAULT_VOXEL2FACE_REFERENCES;
		boundaryBucketFast[posmax].thefaces[i].nextfreevoxelslot = THEFIRSTONE;

			//add the voxelpair
			boundaryBucketFast[posmax].thefaces[i].voxelizedface[0].location = MARKER_TO_IDENTIFY_EXPELLED_SITES + loc;
			boundaryBucketFast[posmax].thefaces[i].voxelizedface[0].seedid = seed0;
			boundaryBucketFast[posmax].thefaces[i].nextfreevoxelslot++;

		boundaryBucketFast[posmax].len++;
		gbCountFast++;

		return;
	}

	//##DEBUG
	QUICKASSERT ( boundaryBucketFast[posmax].len == boundaryBucketFast[posmax].maxLen );

	//add further faces into that collection
	size_t oSize = boundaryBucketFast[posmax].maxLen;
	size_t nSize = oSize + STDLEN;

	bndFaceFastP newfaces = NULL;
	newfaces = new bndFaceFast[nSize];
	QUICKASSERT( newfaces != NULL );
	myMemGuard = myMemGuard + ( nSize * sizeof(bndFaceFast) );

	//copy the old face data first
	bndFaceFastP ofaces = boundaryBucketFast[posmax].thefaces;
	for ( uint32_t oldf = 0; oldf < oSize; oldf++ ) {
		newfaces[oldf].gposUp = ofaces[oldf].gposUp;
		newfaces[oldf].gposDown = ofaces[oldf].gposDown;
		newfaces[oldf].id = ofaces[oldf].id;
		newfaces[oldf].disori = ofaces[oldf].disori;
		//carry over the pointer to memory to avoid leaks
		newfaces[oldf].voxelizedface = ofaces[oldf].voxelizedface;
		newfaces[oldf].nvoxeltotal = ofaces[oldf].nvoxeltotal;
		newfaces[oldf].nextfreevoxelslot = ofaces[oldf].nextfreevoxelslot;
	}

	//initialize the rest for safety
	for ( uint32_t newf = oSize; newf < nSize; newf++ ) {
		newfaces[newf].voxelizedface = NULL;
		newfaces[newf].nvoxeltotal = EMPTY;
		newfaces[newf].nextfreevoxelslot = THEFIRSTONE;
	}


	//add data on the i = maxLen-th new face to the boundaryBucket[posmax]
	newfaces[i].gposUp = MAX(seed0, seedtt);
	//QUICKASSERT ( newfaces[i].gposUp == pUp ); //##MK::because risky, initially promoted int to long downgraded again to unit_fast32t...

	newfaces[i].gposDown = MIN(seed0, seedtt); 
	//QUICKASSERT ( newfaces[i].gposDown == pDown );

	newfaces[i].id = hashtag;
	newfaces[i].disori = GBDISORI_NOT_REQUIRED_YET;

	newfaces[i].voxelizedface = NULL;
	newfaces[i].voxelizedface = new cellsBndFast[DEFAULT_VOXEL2FACE_REFERENCES];
	QUICKASSERT( newfaces[i].voxelizedface != NULL );
	myMemGuard = myMemGuard + (DEFAULT_VOXEL2FACE_REFERENCES * sizeof(cellsBndFast));

	newfaces[i].nvoxeltotal = DEFAULT_VOXEL2FACE_REFERENCES;
	newfaces[i].nextfreevoxelslot = THEFIRSTONE;


		//add the voxelpair in this new face
		newfaces[i].voxelizedface[0].location = MARKER_TO_IDENTIFY_EXPELLED_SITES + loc;
		newfaces[i].voxelizedface[0].seedid = seed0;
		newfaces[i].nextfreevoxelslot++;


	//handshake hand over
	delete [] boundaryBucketFast[posmax].thefaces; //the memory that was referred to in thefaces[..].voxelizedface was copied over
	myMemGuard = myMemGuard - (oSize * sizeof(bndFaceFast));

	boundaryBucketFast[posmax].thefaces = newfaces;
	boundaryBucketFast[posmax].maxLen = nSize;
	boundaryBucketFast[posmax].len++;
	gbCountFast++;
}


double caHdl::calculateBoundaryDisoriFast( uint_fast32_t seedup, uint_fast32_t seeddown )
{
	double q1[4], q2[4], qdis[4];
	uint32_t upOri = mydefgpool[tmpdefgseeds[seedup].mydefgpoolid].caori;
	uint32_t downOri = mydefgpool[tmpdefgseeds[seeddown].mydefgpoolid].caori;


	//now calculate on the fly...
	q1[0] = myoripool[upOri].q0;
	q1[1] = myoripool[upOri].q1;
	q1[2] = myoripool[upOri].q2;
	q1[3] = myoripool[upOri].q3;

	q2[0] = myoripool[downOri].q0;
	q2[1] = myoripool[downOri].q1;
	q2[2] = myoripool[downOri].q2;
	q2[3] = myoripool[downOri].q3;

	misorientationQuaternionCubic( q1, q2, qdis );
	double theta = qdis[0];
	if( theta > 1.0 ) {
		theta = (double) (int) theta;
	}
	theta = 2*acos(theta);

	return theta;
}


void caHdl::trimGBMemoryAndCalculateDisoriFast( void )
{
	uint32_t totalVoxelAtBoundaries = 0;

	for ( uint32_t pos = 0; pos < boundaryBucketSizeFast; pos++ ) {
		//no grain linked?
		if ( boundaryBucketFast[pos].len == 0 ) continue; //precached maybe even detrimental

		for ( uint32_t f = 0; f < boundaryBucketFast[pos].len; f++) {
			boundaryBucketFast[pos].thefaces[f].disori = calculateBoundaryDisoriFast( boundaryBucketFast[pos].thefaces[f].gposUp, boundaryBucketFast[pos].thefaces[f].gposDown );

			//trim boundary bucket to release memory
			cellsBndFast* toolongarray = boundaryBucketFast[pos].thefaces[f].voxelizedface;

			uint32_t toolong = boundaryBucketFast[pos].thefaces[f].nvoxeltotal;
			uint32_t exactSize = boundaryBucketFast[pos].thefaces[f].nextfreevoxelslot;
			QUICKASSERT ( exactSize < CA_ALLOCATION_MAXIMUM );

			cellsBndFast* reservoir = NULL;
			reservoir = new cellsBndFast[exactSize];
			QUICKASSERT( reservoir != NULL );
			myMemGuard = myMemGuard + (exactSize * sizeof(cellsBndFast));

			//carry over
			for ( uint32_t vxp = 0; vxp < exactSize; vxp++) {
				reservoir[vxp].location = toolongarray[vxp].location;
				reservoir[vxp].seedid = toolongarray[vxp].seedid;
			}

			//free the obsolete and too long snippet
			delete [] boundaryBucketFast[pos].thefaces[f].voxelizedface;
			myMemGuard = myMemGuard - (toolong * sizeof(cellsBndFast));
			totalVoxelAtBoundaries = totalVoxelAtBoundaries + exactSize;

			//link in the new exactly trimmed snippet such nucleation can pick from linear positions in array rather random pointer hinting in more cases to memory farther away
			boundaryBucketFast[pos].thefaces[f].voxelizedface = reservoir;
			boundaryBucketFast[pos].thefaces[f].nvoxeltotal = exactSize;
			boundaryBucketFast[pos].thefaces[f].nextfreevoxelslot = exactSize;
		}
	}

	this->SvDeformed = totalVoxelAtBoundaries;

	if ( this->outopt_logboundaries == OUTPUT_LOGBND_YES ) { //##MK::output all boundaries
		stringstream localBoundariesfname;
		ofstream localBoundaries;
		localBoundariesfname << "SCORe." << myensHdl->simid << ".JobID." << this->jobid << ".LogBoundaries.csv";
		localBoundaries.open( localBoundariesfname.str().c_str() );
		localBoundaries << "BoundaryBucketID;BoundaryFaceID;BoundaryCellCnt;Disori;gPosUpID;gPosDownID;rho0gPosUp;rho0gPosDown" << endl;

		for( uint32_t j = 0; j < boundaryBucketSizeFast; j++ ) {
			//no grain linked?
			if ( boundaryBucketFast[j].len == 0 ) continue;

			uint32_t gup, gdw;
			for ( uint32_t f = 0; f < boundaryBucketFast[j].len; f++) {
				gup = boundaryBucketFast[j].thefaces[f].gposUp;
				gdw = boundaryBucketFast[j].thefaces[f].gposDown;

				localBoundaries << j << ";" << f << ";" << boundaryBucketFast[j].thefaces[f].nvoxeltotal << ";" << boundaryBucketFast[j].thefaces[f].disori << ";";
				localBoundaries << gup << ";" << gdw << ";" << mydefgpool[tmpdefgseeds[gup].mydefgpoolid].rho0 << ";" << mydefgpool[tmpdefgseeds[gdw].mydefgpoolid].rho0 << endl;
			}
		}

		localBoundaries.flush();
		localBoundaries.close();
	}

cout << "GrainBoundary Disorientations calculated successfully = bucketSizeFast;gbCountFast;" << boundaryBucketSizeFast << ";" << gbCountFast << endl;
}


uint32_t caHdl::whichSeedAtLocationFast( uint32_t seedid0, int ix , int iy, int iz, short dx, short dy, short dz )
{
	//prepare to become pos or neg, who knows...
	int tx = ix + dx;
	int ty = iy + dy;
	int tz = iz + dz;

	//construct no boundaries over periodic domain boundaries
	if ( tx < 0 )							return seedid0;
	if ( tx >= myCAGeometry.nboxedge_rd )	return seedid0;
	if ( ty < 0 )							return seedid0;
	if ( ty >= myCAGeometry.nboxedge_nd )	return seedid0;
	if ( tz < 0 )							return seedid0;
	if ( tz >= myCAGeometry.nboxedge_td )	return seedid0;

	//so tx, ty, tz are always positive and confined in a smaller than CUBE(CA_DIMENSIONS_MAXIMUM)^3 container
	//convert ix,iy,iz in implicit 3D coordinate to access mycellgrid
	uint32_t txyz = tx + (ty * myCAGeometry.nboxedge_rd) + (tz * myCAGeometry.nboxarea_rdnd);

	uint32_t seedidt = mycellgrid[txyz]; //all indices are positive

	return seedidt;
}


void caHdl::detectGrainBoundaryCellsFast( void )
{
	//MK::collect pixel close to disjoint seedids
	//scan through the whole automaton for each cell in either direction test
	uint32_t nbx = myCAGeometry.nboxedge_rd;
	uint32_t nby = myCAGeometry.nboxedge_nd;
	uint32_t nbz = myCAGeometry.nboxedge_td;
	uint32_t nbxy = nbx * nby;
	uint32_t nbxyz = nbx * nby * nbz;

	uint32_t seed0, seedtest;
	for ( uint32_t c = 0; c < nbxyz; c++) {

		//referencing entries of tmpdefgseeds are positive and < UINT32MAX
		seed0 = mycellgrid[c]; 

		//coordinates in the automaton are always positive and ##MK wont work for 1x2 pixelx10000 cells automata!
		uint32_t sz = c / nbxy;
		uint32_t rem = c - (nbxy * sz);
		uint32_t sy = rem / nbx;
		uint32_t sx = rem - (nbx * sy);
		//which SeedAtLocationFast works with int might flow over BUT nbx, nby, nbz is carefully limited

		//improve cache locality
		seedtest = whichSeedAtLocationFast( seed0, sx, sy, sz, -1, -1, -1 );
		if ( seedtest != seed0 ) { addVoxelAtBoundary( seed0, seedtest, c ); continue; }
		seedtest = whichSeedAtLocationFast( seed0, sx, sy, sz,  0, -1, -1 );
		if ( seedtest != seed0 ) { addVoxelAtBoundary( seed0, seedtest, c ); continue; }
		seedtest = whichSeedAtLocationFast( seed0, sx, sy, sz, +1, -1, -1 );
		if ( seedtest != seed0 ) { addVoxelAtBoundary( seed0, seedtest, c ); continue; }
		seedtest = whichSeedAtLocationFast( seed0, sx, sy, sz, -1,  0, -1 );
		if ( seedtest != seed0 ) { addVoxelAtBoundary( seed0, seedtest, c ); continue; }
		seedtest = whichSeedAtLocationFast( seed0, sx, sy, sz,  0,  0, -1 );
		if ( seedtest != seed0 ) { addVoxelAtBoundary( seed0, seedtest, c ); continue; }
		seedtest = whichSeedAtLocationFast( seed0, sx, sy, sz, +1,  0, -1 );
		if ( seedtest != seed0 ) { addVoxelAtBoundary( seed0, seedtest, c ); continue; }
		seedtest = whichSeedAtLocationFast( seed0, sx, sy, sz, -1, +1, -1 );
		if ( seedtest != seed0 ) { addVoxelAtBoundary( seed0, seedtest, c ); continue; }
		seedtest = whichSeedAtLocationFast( seed0, sx, sy, sz,  0, +1, -1 );
		if ( seedtest != seed0 ) { addVoxelAtBoundary( seed0, seedtest, c ); continue; }
		seedtest = whichSeedAtLocationFast( seed0, sx, sy, sz, +1, +1, -1 );
		if ( seedtest != seed0 ) { addVoxelAtBoundary( seed0, seedtest, c ); continue; }

		seedtest = whichSeedAtLocationFast( seed0, sx, sy, sz, -1, -1, 0 );
		if ( seedtest != seed0 ) { addVoxelAtBoundary( seed0, seedtest, c ); continue; }
		seedtest = whichSeedAtLocationFast( seed0, sx, sy, sz,  0, -1, 0 );
		if ( seedtest != seed0 ) { addVoxelAtBoundary( seed0, seedtest, c ); continue; }
		seedtest = whichSeedAtLocationFast( seed0, sx, sy, sz, +1, -1, 0 );
		if ( seedtest != seed0 ) { addVoxelAtBoundary( seed0, seedtest, c ); continue; }
		seedtest = whichSeedAtLocationFast( seed0, sx, sy, sz, -1, 0, 0 );
		if ( seedtest != seed0 ) { addVoxelAtBoundary( seed0, seedtest, c ); continue; }
		seedtest = whichSeedAtLocationFast( seed0, sx, sy, sz, +1, 0, 0 );
		if ( seedtest != seed0 ) { addVoxelAtBoundary( seed0, seedtest, c ); continue; }
		seedtest = whichSeedAtLocationFast( seed0, sx, sy, sz, -1, +1, 0 );
		if ( seedtest != seed0 ) { addVoxelAtBoundary( seed0, seedtest, c ); continue; }
		seedtest = whichSeedAtLocationFast( seed0, sx, sy, sz,  0, +1, 0 );
		if ( seedtest != seed0 ) { addVoxelAtBoundary( seed0, seedtest, c ); continue; }
		seedtest = whichSeedAtLocationFast( seed0, sx, sy, sz, +1, +1, 0 );
		if ( seedtest != seed0 ) { addVoxelAtBoundary( seed0, seedtest, c ); continue; }

		seedtest = whichSeedAtLocationFast( seed0, sx, sy, sz, -1, -1, +1 );
		if ( seedtest != seed0 ) { addVoxelAtBoundary( seed0, seedtest, c ); continue; }
		seedtest = whichSeedAtLocationFast( seed0, sx, sy, sz,  0, -1, +1 );
		if ( seedtest != seed0 ) { addVoxelAtBoundary( seed0, seedtest, c ); continue; }
		seedtest = whichSeedAtLocationFast( seed0, sx, sy, sz, +1, -1, +1 );
		if ( seedtest != seed0 ) { addVoxelAtBoundary( seed0, seedtest, c ); continue; }
		seedtest = whichSeedAtLocationFast( seed0, sx, sy, sz, -1,  0, +1 );
		if ( seedtest != seed0 ) { addVoxelAtBoundary( seed0, seedtest, c ); continue; }
		seedtest = whichSeedAtLocationFast( seed0, sx, sy, sz,  0,  0, +1 );
		if ( seedtest != seed0 ) { addVoxelAtBoundary( seed0, seedtest, c ); continue; }
		seedtest = whichSeedAtLocationFast( seed0, sx, sy, sz, +1, 0, +1 );
		if ( seedtest != seed0 ) { addVoxelAtBoundary( seed0, seedtest, c ); continue; }
		seedtest = whichSeedAtLocationFast( seed0, sx, sy, sz, -1, +1, +1 );
		if ( seedtest != seed0 ) { addVoxelAtBoundary( seed0, seedtest, c ); continue; }
		seedtest = whichSeedAtLocationFast( seed0, sx, sy, sz,  0, +1, +1 );
		if ( seedtest != seed0 ) { addVoxelAtBoundary( seed0, seedtest, c ); continue; }
		seedtest = whichSeedAtLocationFast( seed0, sx, sy, sz, +1, +1, +1 );
		if ( seedtest != seed0 ) { addVoxelAtBoundary( seed0, seedtest, c ); continue; }
	}

	//ready for grain boundary nucleation
	cout << "Cells located at grain boundaries were identified with the method fast 2." << endl;
}


void caHdl::visualizeGBTopologyFast( void )
{
	uint32_t nx = myCAGeometry.nboxedge_rd;
	uint32_t ny = myCAGeometry.nboxedge_nd;
	uint32_t nz = myCAGeometry.nboxedge_td;
	uint32_t nxy = myCAGeometry.nboxarea_rdnd;
	uint32_t nxyz = myCAGeometry.nboxvol_rdndtd;

	uint32_t* seedID2raw = new uint32_t[nxyz];
	uint32_t* faceID2raw = new uint32_t[nxyz];

	//initialize rendering values of the containers
	for (uint32_t vx = 0; vx < nxyz; vx++) 
		seedID2raw[vx] = PLOTBOUNDARIES_DEFORMED;

	for (uint32_t vx = 0; vx < nxyz; vx++) 
		faceID2raw[vx] = PLOTBOUNDARIES_DEFORMED;

	//scan all grain boundary faces and mark cells as located at the boundary
	int faceid = 0;
	for (uint32_t b = 0; b < boundaryBucketSizeFast; b++) {
		if ( boundaryBucketFast[b].len == 0 ) continue;

		for ( uint32_t f = 0; f < boundaryBucketFast[b].len; f++ ) {
			cellsBndFast* vxps = boundaryBucketFast[b].thefaces[f].voxelizedface;
			uint32_t nvxp = boundaryBucketFast[b].thefaces[f].nextfreevoxelslot;
			uint32_t gup = boundaryBucketFast[b].thefaces[f].gposUp;
			uint32_t gdw = boundaryBucketFast[b].thefaces[f].gposDown;

			QUICKASSERT( vxps != NULL );

			uint32_t cxyz, seed;
			for ( uint32_t vx = 0; vx < nvxp; vx++) {

				cxyz = vxps[vx].location - MARKER_TO_IDENTIFY_EXPELLED_SITES;
				seed = vxps[vx].seedid;

				//##DEBUG
				QUICKASSERT ( seed == mycellgrid[cxyz] );
				//##DEBUG

				seedID2raw[cxyz] = PLOTBOUNDARIES_DEFORMED + seed;

				faceID2raw[cxyz] = faceid;
			}

			faceid++;
		} //all faces of a face collection
	} //all boundaries


	MPI_File msFileHdlSeedID2Voxel;
	MPI_Status msFileStatusSeedID2Voxel;
	MPI_File msFileHdlBoundaryID2Voxel;
	MPI_Status msFileStatusBoundaryID2Voxel;

	// create C-consistent file name for MPI I/O
	stringstream msFileNameSeed2Voxel;
	msFileNameSeed2Voxel << "SCORE." << myensHdl->simid << ".JobID." <<  this->jobid << ".SeedCellGBAssgn.raw";
	int msFileNameLengthSeed2Voxel = msFileNameSeed2Voxel.str().size();
	char* CmsFileNameSeed2Voxel = new char[msFileNameLengthSeed2Voxel+1];
	strcpy(CmsFileNameSeed2Voxel, msFileNameSeed2Voxel.str().c_str());

	stringstream msFileNameBnd2Voxel;
	msFileNameBnd2Voxel << "SCORE." << myensHdl->simid << ".JobID." <<  this->jobid << ".FaceIDCellAssgn.raw";
	int msFileNameLengthBnd2Voxel = msFileNameBnd2Voxel.str().size();
	char* CmsFileNameBnd2Voxel = new char[msFileNameLengthBnd2Voxel+1];
	strcpy(CmsFileNameBnd2Voxel, msFileNameBnd2Voxel.str().c_str());


	//working all nodes open the file in create and write-only mode
	MPI_File_open(MPI_COMM_SELF, CmsFileNameSeed2Voxel, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &msFileHdlSeedID2Voxel);
	MPI_File_open(MPI_COMM_SELF, CmsFileNameBnd2Voxel, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &msFileHdlBoundaryID2Voxel);

	//MPI I/O
	int totalOffset = 0;
	MPI_File_seek( msFileHdlSeedID2Voxel, totalOffset, MPI_SEEK_SET );
	MPI_File_write(msFileHdlSeedID2Voxel, seedID2raw, nxyz, MPI_INT, &msFileStatusSeedID2Voxel);

	totalOffset = 0;
	MPI_File_seek( msFileHdlBoundaryID2Voxel, totalOffset, MPI_SEEK_SET );
	MPI_File_write(msFileHdlBoundaryID2Voxel, faceID2raw, nxyz, MPI_INT, &msFileStatusBoundaryID2Voxel);

	MPI_File_close(&msFileHdlSeedID2Voxel); //no Barrier because MPI self
	MPI_File_close(&msFileHdlBoundaryID2Voxel);

	delete [] seedID2raw;
	delete [] faceID2raw;
	delete [] CmsFileNameSeed2Voxel;
	delete [] CmsFileNameBnd2Voxel;
cout << "Successful MPI I/O of boundary cell assignment." << endl;
}
