//SCORE automaton developed by M K{\"u}hbach, 2014/2015, for questions and details contact markus.kuehbach@rwth-aachen.de

#ifndef __SCORE_DEFINITIONS_H_INCLUDED__
#define	__SCORE_DEFINITIONS_H_INCLUDED__


#include <mpi.h>


//maximum dimensions and hard limits owing to uint32_t based datastructures and algorithms
//MK::DO APPLY GREATEST CARE WHEN CHANGING THESE WELL BALANCED AND HARD LIMITATIONS OF THE CODE!
#define UINT32T_MAX					4294967294 //2^32-1 otherwise wrap around but owing to MARKER_TO_IDENTIFY_EXPELLED_SITES better minus another one

#define CA_DIMENSIONS_MAXIMUM		1600	//to stay within UINT32_MAX range
#define CA_DIMENSIONS_MINIMUM		50		//so that nxyz / nxy is not resulting in wrap around and inside/outside bounds check not run in wraparound
#define CA_ALLOCATION_MAXIMUM		4096000000
#define CA_GRAINIDS_MAXIMUM			4096000000 //not more than 1600^3 disjoint IDs, anyway, discretizing so many grains in 1625^3 cells is way too inaccurate!
//##MK::MIND UINT32T ARITHMETICS CHOOSE CARE WHEN CHANGING VALUES HERE!
#define NO_GRAIN_ASSIGNED			4096000001
#define I_DONT_KNOW_YET				4096000001
#define NOT_ASSIGNED_YET			4096000001
#define CELL_IS_A_PARTICLE			4096001000
#define CURRENTLY_INFECTED			4096002000
#define FULLY_RECRYSTALLIZED		4096003000 //defgid run from 0 to ndefgpool - 1
#define NUCLEUS_ALREADY_PLACED		4096004000
#define SITE_ALREADY_CONSUMED		4096006000

//physical quantities
#define MICRON2METER				1.0e-6
#define SCALING_DEFAULT_DISDENS		1.0e14 //m^2
//unsigned char RGBA color model constants
#define SMALLEST_TIME_LOGARITMIZE	(1.0e-9)
#define UCHAR_RANGE_MIN				0
#define UCHAR_RANGE_MAX				255
#define INTEGER_RANGE_MAX			2147483647
#define SHORT_MAX					32767
#define RED							0
#define GREEN						1
#define BLUE						2
#define ALPHA						3
#define RGB							3
//RGB MAX MAX MAX is white whereas MIN MIN MIN is black


#define MAXFILLPERSTEP				(0.1)

//implemenation relevant defines memory guard
#define BYTE2MEGABYTE				1048576


//integrator accuracy limits to assure reasonable numerical solutions during transient annealing
#define VERYSMALLTIMESTEP			(1e-6)	//s - do not use smaller timesteps than this
#define INITIAL_DELTAT				(1e-6) //micro s
#define INFINITE					(1e36) //s

#define FGEOFACE					(1.0000) //1.0000 //1.0449		//1.0449 // 1 //1
#define FGEOEDGE					(0.7957) //1.1301 //0.8024 //1.1347 // 0.7957 //0.707106781186547
#define FGEODIAG					(0.6494) //1.1000 //0.5821 //1.0082 // 0.6494 //0.577350269189626 //###MK20121211, 26NN


//discretization to time integration scheme
#define SMALL_HEATRATE				(0.000277777777777777) //1 K/ 3600s
#define SMALL_RECOVERYRATE			(2.77e10) // 10^14 1/m^2/3600s
#define SMALL_DRAGGINGRATE			(0.027777777) // 100Pa/3600s
#define SMALL_HEAT					(1.0) //K
#define SMALL_RHO					(1.0e11) //1/m^2 ~ 0.01 rho min
#define SMALL_ZENERFORCE			(100.0) //Pa
#define SMALL_NUMBEROFNUCLEI		(1)
#define SMALL_NUCLEATIONRATE		(2.77e-4) //1 nuc/3600s
#define MINIMUM_RAYLEIGH_SIGMA		(1e-3) //x/sigma^2*exp(-0.5x^2/sigma^2)
//#define SMALL_DISTANCE			this is the dcell
#define SMALL_VELOCITY				(2.7e-13) //1nm/3600s


//physical constants
#define kboltzman					(1.3806488e-23)
#define TOFFSET						(273.15) //degree Celsius into Kelvin
#define echarge						(1.602176565e-19) //Coulomb
#define KUHLMANNFINALRHO			(1e10)



//definition of ideal texture component / ideal / standardlagen
#define RANDOM_ORIENTATION			-1
#define MAX_FCC						(1.099) //MacKenzie..., 62.8/180*_PI_
#define RESOLUTION_SO3GRID			(0.01745329251994330) // 1.0 degree raster
#define MAXDISORI_TO_40DEG111		(0.174532925199433) //10 / 180 * _PI_;
#define MAXDISORI_LAGB2HAGB			(0.261799387799149) //15°
#define NOT_WITHIN_GRIDRESU			-1
#define MAXIMUM_DISJOINT_ORIS		4096000000
#define CATEGORIZED_AS_RANDOM		4096000001
#define UNKNOWN_ORIENTATION			4096000001
#define INVALID_NUCSITE				4096000001



//default values in the construction phase
#define DEFAULT_SIMID				10
#define DEFAULT_CELLSIZE			(1e-6)
#define DEFAULT_SHEARMODULUS		(2.7e10)		//RT pure aluminium
#define DEFAULT_BURGERSVECTOR		(2.86e-10)	//RT pure aluminium
#define DEFAULT_PUREALU_C			(1.0)
#define DEFAULT_PUREALU_TMELT		(933.0)

#define DEFAULT_SMALL_NUMBER		(1e-6)
#define DEFAULT_LAGB_HACT			(1.30)		//eV
#define DEFAULT_HAGB_HACT			(1.20)
#define DEFAULT_GS_HACT				(1.10)
#define DEFAULT_DEFGSIZE			(1e-4)		//m
#define DEFAULT_MAXFILLPERSTEP		(0.1)
#define DEFAULT_MINFILLIN			(0.01)		//smaller is unnecessary fine sampling
#define DEFAULT_MAXFILLIN			(0.2)		//larger too coarse behavior

#define DEFAULT_RELCELLCACHING		(0.10)
#define DEFAULT_TRANSRELCELL		(1.0)
#define DEFAULT_PRNG_SEED			-46356
#define DEFAULT_DELTATIME			(1e-6)		//s
#define DEFAULT_XMAX				(1.0)
#define DEFAULT_TMAX				(600000.0)	//s
#define DEFAULT_NMAX				2000000
#define DEFAULT_NUCDENSITY_CSR		1
#define DEFAULT_CELLSIZE_MIN		(1e-7)		//smaller violates mesoscale behavior of atomic defects!
#define DEFAULT_CELLSIZE_MAX		(1e-5)		//larger cells not expected because of technical grain size barely orders of magnitude of 10micron
#define DEFAULT_MINRELCELL			(0.05)
#define DEFAULT_MAXRELCELL			(0.25)
#define DEFAULT_MINTRANSRELCELL		(0.5)
#define DEFAULT_MAXTRANSRELCELL		(2.0)		//according to GNU
#define DEFAULT_MIN_GRAINDISCR		10			// in cells
#define DEFAULT_PMAX				(0.0)		// assuming that their are fast boundaries in the system
#define MAXATTEMPTS_NUCSITE			100000
#define RHOMAX_WELLANNEALED			(1e8)

//limiting values defining which potential simulation setups can be realized
#define DEFAULT_MAX_NCAS			1048576		//corresponding to 32768 process each of which hosting 32 caHdl
#define DEFAULT_MAX_DEFORMEDGRAINS	65536
#define DEFAULT_MAX_RXGRAINS		1048576
#define DEFAULT_MAX_INPUTDATALEN	4096
#define DEFAULT_MAX_EDGELENGTH		1200		//MK::internal address calculation proceeds as int
#define MIN_DISCRETIZATION_GRAIN	14137		//4/3_PI_15^3 in cells


//nucleation modeling
#define NO_GBNUCLEATION				0x01
#define GB_ALLBOUNDARIES			0x02
#define GB_ONLYATHAGB				0x03
#define NO_CSRNUCLEATION			0x04
#define CSR_ENFORCED				0x05
#define CSR_PICKRANDOMLY			0x06
#define NO_CLUSTNUCLEATION			0x07
#define CLUSTNUCLEATION_PICKRANDOMLY		0x08
#define TINCUB_SITESATURATION		0x09
#define TINCUB_TIMEDEPENDENT		0x10


#define NUCSITES_STILL_SOME_FREE	0x01
#define NUCSITES_ALL_EXHAUSTED		0x02


//dispersoid drag
#define DISPERSOIDDRAG_NO			0x01
#define DISPERSOIDDRAG_CONSTANT		0x02
#define DISPERSOIDDRAG_TIMEDEP		0x03

#define MOBILITYMODEL_SEBALDGOTTSTEIN	0x01
#define MOBILITYMODEL_ROLLETTHOLM		0x02

#define RENDERING_MSNO				0x01
#define RENDERING_MS2D				0x02
#define RENDERING_MS3D				0x03

#define RENDERING_BOUNDARIES_NO		0x01
#define RENDERING_BOUNDARIES_YES	0x02

#define RENDERING_COLOR_GRAINID		0x01
#define RENDERING_COLOR_IPFZ		0x02

#define OUTPUT_LOGBND_NO			0x01
#define	OUTPUT_LOGBND_YES			0x02

#define OUTPUT_RXFRONTSTATS_NO		0x01
#define OUTPUT_RXFRONTSTATS_YES		0x02

#define OUTPUT_SINGLEGRAIN_NO		0x01
#define OUTPUT_SINGLEGRAIN_YES		0x02


//deformation structures
#define CUBOID_DEFMS				0x01
#define POISSONVORONOI_DEFMS		0x02
#define CPFEMDAMASK_DEFMS			0x03

#define MASTER_CSR_NPOINTS			1000000
#define XMI							0
#define XMX							1
#define YMI							2
#define YMX							3
#define ZMI							4
#define ZMX							5


//parallelization
#define MASTER						0

//defaults for the load calculation analysis
#define DEFAULT_MEMPERCELL_DEF		4			//int for each cell
#define DEFAULT_MEMPERCELL_RX		24			//active cells
#define DEFAULT_MEMPERGRAIN_LOG		8			//double for cell volume
#define DEFAULT_MEMPERCA_OVERHEAD	20480		//20MB/CA


//defines for deformation microstructure construction
#define DEFAULT_NTSTEPSMAX_DEFSYNTH	1000000


//grain boundary tracking management
#define MYHASH(a,b)					( (((uint_fast64_t) max(a,b)) << 32) | ((uint_fast64_t) min(a,b)) )
#define STDLEN						16	//has to be positive
#define BNDBIN_CACHE_MULTIPLIER		2
#define REPEAT_MAX_DRAWING			1000
#define GBDISORI_NOT_REQUIRED_YET	(0.0)
#define DEFAULT_VOXEL2FACE_REFERENCES	1000
#define EMPTY						0
#define THEFIRSTONE					0
#define MARKER_TO_IDENTIFY_EXPELLED_SITES 1 //MK::MUST NOT BE ZERO!!!!

//grain boundary tracking fast
#define DEFAULT_ALLOC_VOXELPAIRS	(0.25)
#define DEFAULT_REALLOC_VOXELPAIRS	(0.05)
#define DEFAULT_NSTANDARDLAGEN		30


//defines for CELL MANAGEMENT
#define NOTHING_TO_RECYCLE			0
#define	ACTIVE						true
#define INACTIVE					false
#define	CELLAPPEND					1 //everything not 0 is true
#define CELLRECYCLE					0 //is false by definition
#define NO_INFECTION				(0.0) //recrystallized fraction
#define ALMOST_FULLY_INFECTED		(1.0)

#define	PLOTBOUNDARIES_DEFORMED		1

#define INFECTOR_DEEP_IN_THE_CUBE			0x88
#define INFECTOR_CLOSE_TO_THE_DOMAINWALL	0x22

//by analyzing the microstructural path an estimate that in a front tracking approach on the fraction of 0.05-0.15 of all voxel
//are in one timestep participating in the transformation, however on average only maxfillperstep of them become fully transformed in that timestep
//so it is the objective of the code to make use of this fact that most ACTIVE cells are not able to induce a transformation of neighboring cells
#define FULLRXCACHING_ATT_FACTOR	(1.5)
#define FULLRECYCLING_ATT_FACTOR	(3.0)
#define MINLENGTH_ACTIVELIST		1000
#define MINLENGTH_RECYCLIST			1000
#define MINLENGTH_FULLSEEDLIST		1000


//raw data plotting constants
#define	RAW_PARTICLE				0
#define RAW_INFECTED				1


//rediscretization of ensemble information
#define REDISCR_TIMESLOTS_DEFAULT	1000
#define REDISCR_TIMESLOTS_MAX		10000
#define REDISCR_TIMESLOTS_MIN		100
#define REDISCR_TIME_EQUIDISTANT	0
#define REDISCR_TIME_LOGDISTANT		1
#define REDISCR_DTMIN				(1.0e-6)


//defines associated with the dataholder of cellFast structs
#define _MOBFACUNIT					0xFFFF
#define MOBFACUNIT					(1.0/_MOBFACUNIT)


//recovery modeling
#define RECOVERY_NO					0x01
#define RECOVERY_NES_VACCOREDIFF	0x02
#define RECOVERY_NES_SOLUTEDRAG		0x03
#define RECOVERY_NES_MICHALAKPAXTON	0x04


//correlation functions characterizing growth environment
#define DEFAULT_HPF_WINDOWSIZE		50 //cells


//##MK::first trial for Christian Haases problem
#define CPFEM_NGRAINSX				20
#define CPFEM_NGRAINSY				20
#define CPFEM_NGRAINSZ				20
#define CPFEM_NGRAINSXY				( (CPFEM_NGRAINSX) * (CPFEM_NGRAINSY) )
#define CPFEM_NGRAINSXYZ			( (CPFEM_NGRAINSX) * (CPFEM_NGRAINSY) * (CPFEM_NGRAINSZ) )
#define CPFEM_GRAINSIZEX			(59.0e-6)
#define CPFEM_GRAINSIZEY			(14.0e-6)
#define CPFEM_GRAINSIZEZ			(14.0e-6)
#define	SUBSEQ_VOXELIZATION			0


//MK::grain boundary nucleation number model
#define NO_NUCLEI				0
#define MINIMUM_DRHO			(1.0e1)
#define SCALING_LAMBDA			(10.0)
#define POISSON_CUMSUM_CUTOFF	19
#define POISSON_CUMSUM_TABLE	20 //cutoff+1
#define SCALE_NUCDENSITY		(1.0e-4) //1/micron^2
#define	LAGB_TO_HAGB_TRANS		0.2618	//15/180.0*_PI_

//MK::clustering nucleation model
#define DEFAULT_KMAX_POISSRND			5000 
#define DEFAULT_MASTERMATERN_MAXRADIUS	(0.1)


//lodepng
#define REDCHAN						0
#define GREENCHAN					1
#define BLUECHAN					2
#define ALPHACHAN					3
#define	DEFAULT_ZSECTIONING_ZPOS	(0.5)
#define BLACK						0
#define WHITE						255


//I/O relevant
#define MAXIMUM_NUMBER_OF_NUCLEI_MPI_BUFFERED	((10) * (1000000)) //*(8+8+4+4) byte internal buffer necessary on master!


#endif
