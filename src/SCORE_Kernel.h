//SCORE automaton developed by M K{\"u}hbach, 2014/2015, for questions and details contact markus.kuehbach@rwth-aachen.de
 

#ifndef __SCORE_KERNEL_H_INCLUDED__
#define __SCORE_KERNEL_H_INCLUDED__


#include "SCORE_Definitions.h"
#include "SCORE_Io.h"
#include "SCORE_Math.h"

#include <vector>
#include <algorithm>
#include <iomanip>


//fragmentation radar
//#define FRAGMENTATION_RADAR
#define BLOCKLENGTH			1000		//over how many cells are the fragmentation data averaged?
#define REDSTRETCH			5.0			//scaling the R value 255*(1.0-x^REDSTRETCH)
#define EXPECTEDSHARE		0.2



using namespace std;


class caHdl;


//the dynamic heap object that stores a cell that through its infection
struct defcell
{
	bool activity;
	unsigned char infector;
	short ix;
	short iy;
	short iz;
	double frac;
	//P and myrxgid is not necessary as mobility in the microstructure generator is always vmax and defgseeds grows into VACUUM
	uint32_t mydefgseedid;	//the grain that becomes synthetized
};
typedef struct defcell * defcellP;


//the dynamic heap object that stores the active cells
struct cell
{
	bool activity;
	//unsigned char rhoFac_RFU;
	unsigned char infector;
	short ix;
	short iy;
	short iz;
	double rxFrac;
	double P;
	//##MK::in order to safe memory it could be worthwhile to discretize P and rxFrac as int
	uint32_t mydefgid;
	uint32_t myrxgid;
};
typedef struct cell * cellP;


//the recrystallized grain
struct rxg						//rxg - recrystallizing grains
{
	uint32_t ori;				//an ori-index from the oripool
	double tincub;				//an incubation time in seconds after which the nucleus starts growing
};
typedef struct rxg * rxgP;


//a deformed grain that can be consumed
struct defg						//defg - deformed grains
{
	uint32_t ori;				//an ori index from the oripool
	double rho0;				//initial dislocation density
	double rho;					//instantaneous dislocation density
	double dav0;				//initial average subgrain size
	double dav;					//instantaneous average subgrain size
	double avdg0;				//initial average disorientation angle when the disorientation among the subgrains in this grain is described by a Rayleigh-type distribution
	double avdg;				//instantaneous average disorientation angle
};
typedef struct defg * defgP;


//a crystallographic orientation
struct ori
{
	double bunge1;				//Bunge (3,1,3) convention Euler angles
	double bunge2;
	double bunge3;

	double q0;					//unit quaternion representing this orientation
	double q1;
	double q2;
	double q3;

	uint32_t closestideal;		//the discrete closest orientation to which the grain is associated default is RANDOM_ORIENTATION

	unsigned char RGBA[4];		//a color/alpha channel information in accord with the RGB color model
};
typedef struct ori * oriP;


//a crystalsymmetry operator in quaternion notation
struct qsymm_fcc
{
	double q0;					//unit quaternion representing this orientation
	double q1;
	double q2;
	double q3;
};
typedef struct qsymm_fcc * qsymm_fccP;


struct cadefg
{
	uint32_t worlddefgid;				//where in my ensHdl I can find this grain
	uint32_t caori;						//ori index in caHdl local myoripool
	uint32_t cellcount;					//number of cells comprising this grain
	double rho0;
	double rho;
	double dav0;
	double dav;
	double avdg0;
	double avdg;
};
typedef struct cadefg * cadefgP;


struct carxg
{
	uint32_t caori;						//ori index in caHdl local myoripool
	uint32_t cellcount;					//number of cells comprising this orientation
	uint32_t nucsite;					//serves as a switch flag, initially the coordinate is placed, once placed or identified as an unoccupyable place marked with flag ,implicit 3D coordinate ix + iy*box + iz * boxxy
	uint32_t startingsite;				//implicit 3D coordinate ix + iy*box + iz * boxxy persitent can be utilized to get the original location of the nucleus
	double tincub;
};
typedef struct carxg * carxgP;


//an ideal orientation and up to which scatter it is considered as such
struct ideal
{
	double bunge1;				//Bunge (3,1,3) convention Euler angles
	double bunge2;
	double bunge3;

	double q0;					//unit quaternion representing the orientation
	double q1;
	double q2;
	double q3;

	double scatter;				//scatter range in degrees utilized
};
typedef struct ideal * idealP;


//##MK::add further if desired
struct physData
{
	double G0;					//0K Shear modulus
	double dGdt;				//linear shear modulus temperature-dependence model
	double bZeroCelsius;		//Zero Celsius Burgers vector
	double thermexp_C;			//linear thermal expansion coefficient
	double thermexp_a;			//first order value
	double thermexp_b;			//second order
	double Tmelt;				//melting temperature

	double LAGBm0;				//low-angle grain boundary preexponential factor intrinsic mobility
	double LAGBHact;			//activation enthalpy
	double HAGBm0;				//general high-angle grain boundary
	double HAGBHact;
	double GSm0;				//close to 40deg111 in disorientation space misoriented grains that are eligible for faster boundaries ie for Al
	double GSHact;

	double RH_HAGBm0;			//Rollett Holm model
	double RH_HAGBHact;
	double RH_LAGBHAGBcut;
	double RH_LAGBHAGBtrans;
	double RH_LAGBHAGBexponent;

	double defgmean_rd;
	double defgmean_nd;
	double defgmean_td;
	double defgmean_poisson;
};
typedef struct physData * physDataP;


struct dragData
{
	double ZenerAlpha;			//3/2 or close to it
	double ZenerGamma;			//interfacial energy
	double fr;					//dispersion degree

	unsigned char ZenerConsider;

	//add everything for solute drag
};
typedef struct dragData * dragDataP;


struct recoveryData
{
	double VacancyDiffGeo;
	double VacancyDiffD0;
	double VacancyDiffHact;
	double SoluteDiffD0;
	double SoluteDiffHact;
	double SoluteLsPropFactor;
	double SoluteConcentration;
	double NesAlpha3;
	double NesKappa2;
	double NesC3;
	double NesC4;
	double NesC5;

	unsigned char RecoveryConsider;

	//##MK::add further recovery models here
};
typedef struct recoveryData * recoveryDataP;


struct userCAGeom
{
	uint32_t nNucleiCSR;				//so many nuclei in case of prescribed fixed number of nuclei
	uint32_t nboxedge_rd;			//so many cells in each direction
	uint32_t nboxedge_nd;
	uint32_t nboxedge_td;
	uint32_t nboxarea_rdnd;
	uint32_t nboxvol_rdndtd;

	double cellsize;			//in (m)
	double boxedge_rd;
	double boxedge_nd;
	double boxedge_td;
	double boxarea_rdnd;
	double boxvol_rdndtd;
};
typedef struct userCAGeom * userCAGeomP;


struct userCADefMS
{
	unsigned char defmstype;
	uint32_t ngrx;			//GIA grain grid
	uint32_t ngry;
	uint32_t ngrz;
	uint32_t ngrxy;
	uint32_t ngrxyz;
	double u_xrd;		//coordinate of the GIA/grain grid origin to the local box origin absolute size in (m)
	double v_ynd;
	double w_ztd;
};
typedef struct userCADefMS * userCADefMSP;


struct nucData
{
	unsigned char gbnucleation;
	unsigned char csrnucleation;
	unsigned char clustnucleation;
	unsigned char tincubmodel;

	//##MK::add further parameter here
	uint32_t defaultnucdensity;
	double cluster_nclust;
	double cluster_lambda;
	double cluster_rvesize;
	double cluster_a;
	double cluster_b;
	double cluster_c;
	double gbnucleation_dens2num;
	double gbnucleation_drho2dens;
	double gbnucleation_scatter;
	double tincub_rayleigh_sigma;
};
typedef struct nucData * nucDataP;


struct jobticket
{
	int jobid;
	double estimatedMemory;
	double loadfactor;
};
typedef struct jobticket * jobticketP;


struct agrain
{
	double vol;
	double tincub;
	uint32_t ideal;
	int rank;
};
typedef struct agrain * agrainP;


struct profilingData
{
	long JobID;
	long Iterations;
	double tend;
	double MPIWTimeSpendDefMS;
	double MPIWTimeSpendGBDetection;
	double MPIWTimeSpendNucleation;
	double MPIWTimeSpendGrowthSim;
};
typedef struct profilingData * profilingDataP;


typedef struct
{
	long JobID;
	long Iterations;
	double tend;
	double MPIWTimeSpendDefMS;
	double MPIWTimeSpendGBDetection;
	double MPIWTimeSpendNucleation;
	double MPIWTimeSpendGrowthSim;
} MPI_IO_CAProfilingInfoData;

//##MK::further optimization necessary
typedef struct
{
	long nx;
	long ny;
	long nz;
	long nxyz;
	long nboundarycells;
	long ndgrseeds;
	long nrxgrseeds;
	long ndefmicrotexture;
	long nnucleimicrotexture;
	double storedenergy;
} MPI_IO_CAPhysicsInfoData;


typedef struct
{
	double finalvol;
	double tincub;
	uint32_t ideal;
} MPI_IO_FinalGSDInfoData;


typedef struct
{
	double finalvol;
	double tincub;
	int ideal;
	int rank;
} MPI_IO_FGSDComplete;


struct loginfo_ca_physics
{
	uint32_t nx;
	uint32_t ny;
	uint32_t nz;
	uint32_t nxyz;
	uint32_t nboundarycells; //how many cells where identified close to a boundary
	uint32_t ndgrseeds;		//how many deformed grain seeds where placed
	uint32_t nrxgrseeds;	//how many nuclei were place
	uint32_t ndefmicrotexture;	//into how many disjoint components was the deformation structure split
	uint32_t nnucleimicrotexture;
	uint32_t* defmicrotexture;	//how many cells
	uint32_t* nucleimicrotexture; //how many nuclei found in each texture class
	double storedenergy;
};
typedef struct loginfo_ca_physics * loginfo_ca_physicsP;


//the hostgrain switching function as defined by M K{\"u}hbach
struct hpf
{
	int dx;		//+- di in either direction
	int dy;
	int dz;
	int nx;		//2*di + 1, the size of the 3D window on which the function is calculated
	int ny;
	int nz;
	int nxy;
	int nxyz;
};
typedef struct hpf * hpfP;


struct defgseed
{
	uint32_t ensdefgpoolid;
	uint32_t mydefgpoolid;
	uint32_t location;			//linear implicitly addressed position along x, stacked in y, xy slabs stacked in z
};
typedef struct defgseed * defgseedP;


struct point
{
	double x;
	double y;
	double z;
};
typedef struct point * pointP;


//##MK::uint32_t Boundary
struct cellsBndFast
{
	//utilizes the fact that the dimensions of the automaton are constraint to uint32_t
	uint32_t location;
	uint32_t seedid;
};
typedef struct cellsBndFast * cellsBndFastP;


struct bndFaceFast
{
	double disori;					//disori angle of the two grains
	cellsBndFast* voxelizedface;	//the voxel pair that make up the boundary ALWAYS + 1 to have zero entries denoting positions at which nucleation criteria have already been analyzed
	uint_fast64_t id;				//a unique hashtag (MAXMIN) of disjoint IDs
	uint_fast32_t gposUp;			//the two grains adjoining the face
	uint_fast32_t gposDown;
	uint32_t nvoxeltotal;
	uint32_t nextfreevoxelslot;
};
typedef struct bndFaceFast * bndFaceFastP;


struct bndColumnFast
{
	bndFaceFastP thefaces;			//the datastructure carrying metadata and references to voxel of the face
	uint32_t len;					//actual number of faces
	uint32_t maxLen;				//preallocated buffer for faces
};
typedef struct bndColumnFast * bndColumnFastP;



struct loginfo_grainevo_ca
{
	double localtime;			//the local time when the timestamp was taken
	double localX;				//the local recrystallized fraction as the timestamp was taken
	uint32_t locallogcnt;
	uint32_t localtstep;				//the local integration step when the log was taken
	uint32_t localSv;				//the local number of cells active in the whole domain
	uint32_t nlocaldata;				//how many grains are tracked here
	uint32_t* localdatabucket;		//pointing to an array of nlocaldata elements where the volume of the grains is stored
};
typedef struct loginfo_grainevo_ca * loginfo_grainevo_caP;


struct loginfo_fragmentation_radar
{
	uint32_t localtstep;		//at which step was the radar value taken
	uint32_t blocklength;		//how many cells consecutive in myRXFront are summarized over
	uint32_t nblocks;			//how many blocks are considered, each block is represented by a pixel the color coding represents either //white, block never utilized, black to red RGB=(255,0,0) 1.0 of all cells in the block that are utilized are ACTIVE, red at most zero cells are active
	unsigned char* rgba;
};
typedef struct loginfo_fragmentation_radar * loginfo_fragmentation_radarP;


struct loginfo_rxfrontstats_ca
{
	double localtime;
	double localX;
	double localmemory;
	double localPmax;
	uint32_t localstep;
	uint32_t localSv;

	//track caching and fragmentation state of RXfront list my
	uint32_t ntotalRXFront;
	uint32_t nextSlotNeverActiveRXFront;				//int references still juvenile/INACTIVE end of the RXFront that has been cached
	uint32_t nCurrentlyActive;						//how many cells are currently in the interval [0, nextSlotNeverActiveRXFront) ACTIVE?
	//FullRXList
	uint32_t ntotalFullRXList;
	uint32_t nextSlotToFullRX;						//pointing in FullRXList which cell should next infect and then be switched off thus claimed as INACTIVE
	//RecyclingList
	uint32_t ntotalRecyclingList;
	uint32_t nextSlotThatBecomesRecycled;			//int id referencing which entry in myRXFront can be reused (as having been flagged INACTIVE) to keep the RXFront as filled/unfragmented and compact/short as possible
	uint32_t firstNotRecycledYet;
	//###add thread-based collection of cell fragmentation and utilization
};
typedef struct loginfo_rxfrontstats_ca * loginfo_memory_caP;


struct rediscr_window
{
	double tensmin;
	double tensmax;
	uint32_t nslots;									//rediscretize in so many disjoint slots [tensmin;..;tensmax]
	unsigned char strategy;						//REDISCR_TIME_EQUIDISTANT, linear in n steps between [tensmin; tensmax], REDISCR_TIME_LOGDISTANT
};
typedef struct rediscr_window * rediscr_windowP;


class caHdl;
typedef caHdl * caHdlP;

class ensembleHdl : public io, public mathMethods
{
	//friend class caHdl ##DEBUG
	//ensembleHdl is a friend of the caHdl so can for IO purpose access output from the caHdl

public:
	ensembleHdl();
	~ensembleHdl();

	//prototypes
	uint32_t check_disjunctness ( double * bunge );	//tests for disjunctness of orientations that are introduced into the orientation pool
	uint32_t get_closest_standardlage ( double * quat );

	void init_ensprng( bool SetDissimilarSeedsForAll ); //true - seed is -1*myRank, false - all the same DEFAULT_PRNG_SEED
	void init_mpidatatypes( void );
	void init_parameter ( void );
	void init_idealorientations( void );
	void init_processing ( void );
	void init_zenerdrag ( void );
	void init_defgpool (void );
	void init_rxgpool ( void );
	inline double calculateWorkload( uint32_t ncells, double averagefillperstep );
	inline double estimateMemory( uint32_t ncells, double frontcelloverhead, uint32_t avnuclei, uint32_t logsteps );
	void init_distributeWorkOnRanks( void );	//does the work distribution at process level, then caHdl takes over
	void init_distributeWorkOnStartedThreads( void );	//starts threads and then further suggest tieing each workPiece to that local threadid
	void SIMULATE_myCAs ( void );			//handles ##MK currently sequential the simulation of all myCAs

	//internally MPI synchronized postprocessing routines, MASTER writing output
	void postprocess_write_mycas( void );
	void postprocess_write_mycas_mpifast( void );
	void postprocess_initrediscrtimes( void );		//set the timeslots where local data are analyzed locally
	void postprocess_init ( void );					//ensembleHdl is friend to caHdl so it can access output data

	void postprocess_rediscr_kinetics( void );
	void postprocess_rediscr_macrotexture( void );
	void postprocess_rediscr_finalgrainsizedistribution( void );
	bool postprocess_rediscr_finalgrainsizedistribution_mpifast( void );
	void destroy_myCAs( void );


	int simid;
	vector<double> UserDefLogPoint_X_CellListDefragment;
	vector<double> UserDefLogPoint_X_Output;
	vector<double> UserDefLogPoint_MS_Rendering;
	vector<uint32_t> UserDefLogPoint_WhichCAtoOutput;

	//my ensemble of CAs
	vector<loginfo_ca_physics> myCAPhysics;	//keeps track about how many nuclei and other quantities in local domains
	vector<profilingData> myCAProfiler;	//keeps information how long individual simulations took
	vector<uint32_t> WhichCAonWhichRank;		//###identifies for each automaton which rank currently work on which other automaton vector[ca] dereferences rank
	vector<uint32_t> myIDs;
	vector<uint32_t> WhichMyCAonWhichThread;	//fixes the WhichMyCA-th entry that is job myIDs[WhichMyCA-th] in the world ID scope local CA to a certain thread - THUS THEN THREADS ARE BOUND TO COREs and THEIR MEMORY UTILIZING A ccNUMA ARCHITECTURE memory is threadlocal
	vector< vector<caHdlP> > myCAs;		//###the local CAs I am working on if each thread works exclusively on one entire automaton

	vector<ori> worldoripool;			//all polyxx know all orientations
	vector<ideal> standardlagen;		//the world famous texture components to discuss in aluminum
	vector<defg> worlddefgpool;
	vector<rxg> worldrxgpool;

	//processing which is the same for all cas in all ensembles
	vector<double> ttime;				//in (s)
	vector<double> ttemperature;		//in (K)
	vector<double> dispersoidtime;		//in (s)
	vector<double> dispersoidfr;		//zenerfac * f/r in (J/m^2 * 1/m)


	double* ensRediscrTime;
	struct rediscr_window ensRediscrWindow;

	//physical properties
	struct userCAGeom ensCAGeometry;
	struct physData ensPhysData;
	struct dragData ensDragData;
	struct recoveryData ensRecoveryModel;
	struct nucData ensNucleationModel;

	//##MK::no recovery, no drag

	//integrator accuracy
	//##DEBUG
	//MPI process relevant information in the parallel environment
	double maxfillperstep;				// how much partial volume of a cell does the fastest boundary sweep in one timestep
	double ensMemGuard;					//how much memory is currently consumed in the process in MB
	double initialRelCellCaching;		//how much of the total number of cells is reserved to cache information
	double transientRelCellRecaching;	//how much of the already existent number of cells is reserve to carry additional cells
	double XMAX;
	double TMAX;
	uint32_t NTSTEPSMAX;
	//double Pmax;

	int myRank;							//my MPI ID in the MPI_COMM_WORLD
	int nRanks;							//total number of MPI processes that work in the world
	uint32_t nworldCAs;						//so many in the whole world
	uint32_t nensCAs;						//so many in my ensemble that is part of the world

	//exiting, IO and file access at MPI process level
	unsigned char defmsmethod;
	unsigned char mobilitymodel;
	unsigned char outopt_rendermethod;
	unsigned char outopt_rendercolormodel;
	unsigned char outopt_logboundaries;
	unsigned char outopt_localrenderboundaries;
	unsigned char outopt_rxfront;
	unsigned char outopt_singlegrainevo;

	bool ensembleSuccess;
	bool outopt_hpf3d2donmaster;


	randomClass ensembleprng;
	FILE * score_input;

	//MPIDatatypes
	MPI_Datatype MPI_IO_CAProfilingInfoData_Type;
	MPI_Datatype MPI_IO_CAPhysicsInfoData_Type;
	MPI_Datatype MPI_IO_FinalGSDInfoData_Type;
	MPI_Datatype MPI_IO_FGSDComplete_Type;

	//performance counters
	double prof_t0;
	double prof_tstartpp;
	double prof_tend;
};
typedef class ensembleHdl * ensembleHdlP;


class caHdl : public randomClass, public mathMethods //, not public io because IO is only done by the ensembleHdl each object of caHdl is associated with
{
	friend class ensembleHdl;

public:
	caHdl();
	~caHdl();

private:
	ensembleHdlP myensHdl;					//to allow the caHdl access to the public available elements of his "group/ensemble" leader at the process level
	int ccnumatid;
	int jobid;								//a unique ID that identifies which automaton this caHdl object represents, its request is restricted as the jobid is private, only the ensembleHdl can ask in his myCAs[threadid][0-(myCAs[threadid].size-1)] which this ID references!

	void cleanupMemoryUsedForGrowthMachine( void );
	uint32_t ca_get_closest_standardlage( double * quat );
	uint32_t ca_check_disjunctness( double * bunge );

	void nes_networkgrowthmodel_vacancycorediff( void );
	void nes_networkgrowthmodel_solutedrag( void );
	void nes_michalakpaxton( void );
	void characterize_nucsite_environment( void );//computing the HPF distance-type function

	void solve_INITIALIZATION ( void );
	void solve_SYNTHETIZED_DEFORMEDSTRUCTURE_AUTOMATON ( void );
	void solve_SYNTHETIZED_DEFORMEDSTRUCTURE_CUBOIDBLOCKS ( void );

	void solve_SYNTHETIZE_DEFORMATIONSTRUCTURE( void );
	void solve_DETECT_GRAINBOUNDARIES( void );
	long gbfacenucnumbermodel_haase( double facearea, double rhodifference, double disori);
	long gbfacenucnumbermodel_phdmk( double facearea, double rhodifference, double disori);
	void solve_nucmodeling_gbnuc_physicalmodelFast( void );
	void solve_GRAINBOUNDARYNUCLEATION( void );
	void solve_nucmodeling_csrenforced( void );
	void solve_nucmodeling_csrpickrandomly( void );
	void init_myMaternUniverse( void );
	void pickrandomly_myMaternUniverse(void );
	void solve_nucmodeling_ellipsoidalcluster_pickrandomly( void );

	void solve_NUCLEATIONMODELING ( void );
	void solve_REPLACE_CA_STRUCTURE( void );
	void solve_RXGROWTH ( void );
	void log_ca_physics( loginfo_ca_physicsP  container );

	//upon initialization data in ensembleHdl serve as a library for all possible grain orientations and mobility weights
	void init_cahdlprng( void );
	void init_parameter ( void );			//local copies from ensembleHdl values in the caHdl
	void init_processing ( void );
	void init_zenerdrag ( void );
	//void init_defgpool (void );			//obsolete part of construct deformation microstructure, not all deformed grain orientations are necessary at the local scale
	//void init_rxgpool ( void );			//obsolete part during nucleation modelling ,not all those as well
	//##MK::of course there is always a tradeoff when each CA itself becomes very diverse and large but that is not the modeling paradigm!

//MK::there are two CAs, one that is utilized to synthetize the deformation structure
	//inline void id2xyz_cpfem ( int id, short * coor);
	//inline int xyz2id_cpfem ( short * coor );
	inline uint32_t get_damask_idperiodic ( int oid, int dnx, int dny, int dnz );
	defgseedP init_discreteWindowFromMasterCSR( uint32_t * howmanyseeds );
	void seed_deformedgrains_poisson( void );
	void grow_deformedgrains_voxelize_poisson( void );
	void grow_deformedgrains_voxelize_damask( void );

	//dependencies for the simpler synthetization automaton
	void grow_deformedgrains_append_to_fullseed_list( uint32_t seedfrontid );
	void grow_deformedgrains_append_to_recyc_list( uint32_t seedfrontid );
	uint32_t grow_deformedgrains_NextFreeSlotInSeedFrontAppendOrRecycle( bool how );
	void grow_deformedgrains_infect_JuvenileNucleation( uint32_t seedid, uint32_t WhichCellToInfect, double frac0 ); //here divine juvenile infection
	void grow_deformedgrains_infect_OutofSeedFront_BoundsCheck( uint32_t defgid, uint32_t seedfrontid, short dx, short dy, short dz, unsigned char direction, double carryOver );
	void grow_deformedgrains_infect_OutofSeedFront_NoBoundsCheck( uint32_t defgid, uint32_t seedfrontid, short dx, short dy, short dz, unsigned char direction, double carryOver );

	void grow_deformedgrains_initialization( void );
	void grow_deformedgrains_init_mySeedFront( void );
	void grow_deformedgrains_placeAllSeeds( void );
	void grow_deformedgrains_defragmentation( void );
	void grow_deformedgrains_calcGrowthStep( void );
	void grow_deformedgrains_updateFullCells( void );
	void grow_deformedgrains_cleanCellManagementDefMS( void );


	//upon synthetizing GIA blocks
	void pick_deformedgrains( void );		//selects from ensembleHdl GIA defgpool deformed grains
	void voxelize_defms( void );

	void determine_initial_deformationtexture_mydefgpoolids(void);
	void determine_initial_nucleitexture_myrxgpoolids( void );
	void init_myRXFront ( void );
	void append_to_fullrx_list( uint32_t rxfrontid );
	void append_to_recycling_list( uint32_t rxfrontid );
	uint32_t get_NextFreeSlotInRXFrontAppendOrRecycle( bool how );
	void defragment_myRXFront( void );			//performs a linear scan on the elements in the range [0; nextSlotNeverActiveRXFront) all INACTIVE cell slots are overwritten from the rightmost ACTIVE, such nextSlotNeverActiveRXFront is reduced by n such successful swops, such the total length for which to test ACTIVITY during calcGrowthStep is reduced as well as cache misses
	void defragmentation( void );
	void sim_myCA_sitesaturatedNucleation( void );
	inline double pdf_rayleigh( double realtime );
	inline double cdf_rayleigh( double realtime );
	void sim_myCA_timedependentNucleation( void );
	void sim_myCA_calcGrowthStep( void );
	void sim_myCA_updateFullRX ( void ); //initiates the infection of new cells, subsequently switches the infector cell inactive to enable recycling of memory
	void sim_myCA_infect_JuvenileNucleation( uint32_t rxgpoolid, uint32_t WhichCelltoInfect, double rxfrac0 ); //here divide juvenile infection
	void sim_myCA_infect_OutofRXFront_BoundsCheck( uint32_t rxgpoolid, uint32_t rxfrontid, short dx, short dy, short dz, unsigned char direction, double carryOver ); //###currently no overshoot! 
	void sim_myCA_infect_OutofRXFront_NoBoundsCheck( uint32_t rxgpoolid, uint32_t rxfrontid, short dx, short dy, short dz, unsigned char direction, double carryOver ); //###currently no overshoot! 


	//handle processing during the simulation
	void update_system_fixedtime( double time );
	void update_temperature( void );				//to the time set in this->t
	void update_atomisticproperties( void );		//G and b and temperature dependence
	void update_intrinsicmobilities( void );		//Grain boundary properties
	void update_microchemistry( void );
	void update_deformedsubstructure(void );

	//accessors to physical properties
	double calc_mobilityweight ( uint32_t rgp, uint32_t dgp );
	double calc_mobilityweight_sebaldgottstein ( uint32_t rgpoolid, uint32_t dgpoolid );
	double calc_mobilityweight_rollettholm ( uint32_t rgp, uint32_t dgp );
	inline double get_currentintrinsicmobility( double Pvalue );
	inline double get_rho ( uint32_t defgid );
	inline double get_zener ( void );

	//helper functions for managing the dynamic adaptive timestepping at the CA scale at any time
	double get_dtmax_instantslope_cells( void ); //double when );
	double get_dtmax_instantslope_temp ( void ); //double when );
	//exchange when recovery and Zener drag is simulated
	inline double get_dtmax_instantslope_rho ( void ); //double when );
	inline double get_dtmax_instantslope_zener ( void ); //double when );
	inline double get_dtmax_instantslope_nucleation( void );
	double get_dtmax_minimum ( void ); //double when );

	void log_colorize_myoripool_ipfz( void );
	void log_initialization( void );
	void log_rxfrontstats( void );						//traces status information on the current state of the automaton
	void log_grainevo( void );							//keeps track of what the individual grains did
	void log_OUTPUT( void );

	double get_interpCellCount( uint32_t localid, double when );

	//allow a particular node to write status information into a plain file
	void write_hpf3d( void );							//export local parent grain switching probability in 3D
	void write_hpf2d( int wz );							//export in 2D
	void write_rxfrontstats( void );					//outputs simulation metadata from the myca-th of my automata in a file
	void write_grainevolution( void );					//outputs grain-resolved metadata from the myca-th of my automata in a file
	void write_zsection_coloring_grainids( void );		//output rendering of the microstructure grainids
	void write_zsection_coloring_ipfz( void );			//output rendering of the microstructure ipfs colored
	void write_voxeldata_coloring_grainids( void );		//##DEBUG::store implicit array with grain IDs that can be read with paraview or avizo
	void write_voxeldata_coloring_ipfz( void );
#ifdef FRAGMENTATION_RADAR
	void write_fragmentation_radar( void );
#endif

	//deformed grain boundary tracking functionalities
	void initializeBoundariesFast( void );
	void addVoxelAtBoundary( uint32_t seed0, uint32_t seedtt, uint32_t site );
	uint32_t whichSeedAtLocationFast( uint32_t seedid0, int ix , int iy, int iz, short dx, short dy, short dz );
	void detectGrainBoundaryCellsFast( void );
	double calculateBoundaryDisoriFast( uint_fast32_t seedup, uint_fast32_t seeddown );
	void trimGBMemoryAndCalculateDisoriFast( void );
	void visualizeGBTopologyFast( void );
	void emptyBoundariesFast( void );

	vector<double> defragRXFront_atthisX;
	vector<double> output_atthisX;						//copied from ensembleHdl when should the automaton log local information?
	vector<double> rendering_atthisX;
	vector<loginfo_rxfrontstats_ca> myrxfrontstatus;	//keeps track of the state of the myRXFront array counters fragmentation and occupancy
	vector<loginfo_grainevo_ca> mygrainevolution;
	vector<loginfo_fragmentation_radar> myrxfrontfragmentation;

	//vector<ideal> standardlagen;

	vector<ori> myoripool;
	vector<cadefg> mydefgpool;
	vector<carxg> myrxgpool;					//directly references orientation and incubation time
	vector<point> mypointprocess;

	//local copies from ensembleHdl for faster access
	struct userCAGeom myCAGeometry;				//copy from the ensemble Hdl
	struct physData myPhysData;
	struct dragData myDragData;
	struct recoveryData myRecoveryModel;
	struct userCADefMS myCADefMS;				//local dimensions of GIA deformed grain grid topology
	struct nucData myNucleationModel;


//3DCA TO SYNTHETIZE DEFORMATION STRUCTURE
	defgseedP tmpdefgseeds;						//a list of deformed grains that is sampled from the input data set and placed at particular positions
	uint32_t ndefgseeds;

	//MK::the deformation synthesis and the subsequent RX CA operate both on mycellgrid but
	//the individual automata have own management operations as complexity is different that can be utilized without function pointing or complex ifs

	uint32_t* myFullSeedList;						//bookkeeps entries that should be considered during updateFull that should be considered
	uint32_t* myRecycList;							//int reference to currently cells that can be utilized again
	defcellP mySeedFront;						//the list of active defcells stored as structs

	//SeedFrontManagement at thread level, ##EACH THREAD HAS OWN CA
	//SeedFront
	uint32_t ntotalSeedFront;
	uint32_t nextSlotNeverActiveSeedFront;			//int references still juvenile/INACTIVE end of the SeedFront that has been cached
	uint32_t nCurrentlyActiveSeedFront;				//how many cells are currently in the interval [0, nextSlotNeverActiveSeedFront) ACTIVE?
	//FullSeedList
	uint32_t ntotalFullSeedList;
	uint32_t nextSlotToFullSeed;						//pointing in FullSeedList which cell should next infect and then be switched off thus claimed as INACTIVE
	//RecyclingList
	uint32_t ntotalRecycList;
	uint32_t nextSlotThatBecomesRecyc;				//int id referencing which entry in myRXFront can be reused (as having been flagged INACTIVE) to keep the RXFront as filled/unfragmented and compact/short as possible
	uint32_t firstNotRecycYet;						//##MK:: in the interval [0;nextToRecycle) all references in RecycleList have already been utilized


//3DCA TO PERFORM RX SIMULATION the following constructs are classic arrays to avoid STL overhead and in particular the invocation of the copy constructor
	//##DEBUGuint32_t* mydefggrid;							//IDs that dereference from defgpool
	uint32_t* mycellgrid;							//the voxelized version of mydefggrid - this is where the cells infect IDS are stored in range ##MK [ 0; .. mydefgpool .. ;(mydefgpool.size+myrxgpool.size) )
	uint32_t* myFullRXList;							//bookkeeps entries that should be considered during updateFullRX that should be considered
	uint32_t* myRecyclingList;						//int reference to currently cells that can be utilized again
	cellP myRXFront;							//the list of active cells stored as structs
	bool* myrxgpoolboundarycontact;				//has the nucleus tried to infect a cell outside the physical domain?, e.g. made use of local subdomain periodicity

	//RXFrontManagement at thread level, ##EACH THREAD HAS OWN CA
	//RXFront
	uint32_t ntotalRXFront;
	uint32_t nextSlotNeverActiveRXFront;				//int references still juvenile/INACTIVE end of the RXFront that has been cached
	uint32_t nCurrentlyActive;						//how many cells are currently in the interval [0, nextSlotNeverActiveRXFront) ACTIVE?
	//FullRXList
	uint32_t ntotalFullRXList;
	uint32_t nextSlotToFullRX;						//pointing in FullRXList which cell should next infect and then be switched off thus claimed as INACTIVE
	//RecyclingList
	uint32_t ntotalRecyclingList;
	uint32_t nextSlotThatBecomesRecycled;						//int id referencing which entry in myRXFront can be reused (as having been flagged INACTIVE) to keep the RXFront as filled/unfragmented and compact/short as possible
	uint32_t firstNotRecycledYet;					//##MK:: in the interval [0;nextToRecycle) all references in RecycleList have already been utilized
	//track fragmentation behavior

	vector<double> ttime;						//in (s)
	vector<double> ttemperature;				//in (K)
	double tprocessingend;

	vector<double> zenertime;					//in (s)
	vector<double> zenerdispersion;				//zenerfac * f/r in (J/m^2 * 1/m)

	//##MK::currently not utilized variables with which to track local MS evolution
	//vector<double> RhoMax;					//local recovery history
	//vector<uint32_t> RhoMaxID;						
	//vector<double> currZenerDrag;
	//vector<double> deltatime;

	//##MK::with a fast container cellsBndFastP
	bndColumnFastP boundaryBucketFast;
	uint32_t gbCountFast;
	uint32_t boundaryBucketSizeFast;
	uint32_t SvDeformed;
	uint32_t DefMicrotextureClasses;
	uint32_t NucleiMicrotextureClasses;
	uint32_t* DefMicrotexture;
	uint32_t* NucleiMicrotexture;


	double myMemGuard;							//how much memory is currently consumed in MB
	double StoredEnergy;
	double CurrentG;
	double CurrentBurgersVector;
	double CurrentTemperature;
	double t;									//current time, X(t), and temperature
	double X;									//current RX fraction
	double Xcells;
	double dXstep;								//accumulated change in rxVolume accounting for overshoot as well
	double dt;									//time increment
	double tsimend;								//maximum physical real time simulated
	double XMAX;
	double TMAX;
	uint32_t stepsimend;
	uint32_t Sv;										//number of cells in the front
	uint32_t step;
	uint32_t NTSTEPSMAX;

	double myMobilityWeightMax;					//maximum mobility weight detected in this automaton, if for instance locally there is no GS or only a factor of the maximum value, the local adaptive timestepping can profit from this information
	double myrhomax;							//maximum dislocation density in the system
	double mypzmin;
	double maxfillperstep;
	double initialRelCellCaching;				//what to do when too few memory in the cells
	double transientRelCellRecaching;

	//updated intrinsic mobilities
	double mLAGB;
	double mHAGB;
	double mGS;
	double mRHHAGB;
	double Gbhalfsq;
	double _cellsize;

	//double _kT;
	//double zenerfac;							//fgeo * gammaGB 

	//performance caching variables and things used to clarify the code
	uint32_t nmyoripool;
	uint32_t nmydefgpool;
	uint32_t nmyrxgpool;
	uint32_t nmynuclei;

	uint32_t loginfo_rxfrontstats_cnt;
	uint32_t loginfo_grainevo_cnt;
	uint32_t loginfo_rendering_cnt;
	uint32_t loginfo_defrag_cnt;

	//correlation functions characterizing the local environment
	struct hpf myhpfwindow;
	int* myhpfvalues;

	//Parallel environment relevant information
	int myensRank;
	int nRanks;

	unsigned char mobilitymodel;
	unsigned char nucleationmodel_status;
	unsigned char outopt_localrenderhow;
	unsigned char outopt_localrendercolor;
	unsigned char outopt_localrenderboundaries;
	unsigned char outopt_logboundaries;
	unsigned char outopt_localrxfront;
	unsigned char outopt_localsinglegrainevo;


	bool mySuccess;
	bool renderingForThisCA;
	bool outopt_localhpf3d2d;

	randomClass localprng;
	//file interaction with respect to input data should be performed by ensembleHdl at processLevel exclusively
};


#endif

