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


//MK::the SCORE --- statistical cellular operator for recrystallization --- modeling approach rests on the assumption
//that it is sufficient to study a statistical ensemble of many periodic independent volume domains (the so-called)
//solitary units rather than to discretize extremely large representative volume elements with only one domain
//this enables to reduce the total number of cells/voxels per simulation domain and hence allows to utilize
//leaner datatype with lower impact on caches, as such all domains geometry is unsigned int with occasional utilization
//of the integer datatype and accepting signed to unsigned conversions where necessary and by assuring well-limited
//automaton size to be safe against range wraparound issues...

#ifndef __SCORE_DEFS_FUNCTIONALITY_H_INCLUDED__
#define	__SCORE_DEFS_FUNCTIONALITY_H_INCLUDED__

//MK::here all default and constant values for the functionality of the implementation should be placed!


//versioning
#define VERSION_MAJOR				1
#define VERSION_MINOR				2
#define VERSION_REVISION			1
#define VERSION_BUILD				0
#define IS_BETAVERSION_YES


//STL functionalities
//...add here...
#include <vector>
#include <algorithm>
#include <iomanip>
#include <map>

//third party stuff


//#include <mpi.h>


//maximum dimensions and hard limits owing to uint32_t based datastructures and algorithms
//MK::DO APPLY GREATEST CARE WHEN CHANGING THESE WELL BALANCED DEFAULTS AS WELL AS HARD LIMITATIONS OF THE CODE!
#define CA_DIMENSIONS_MAXIMUM				1600		//to stay within UINT32_MAX range, remaining headroom of the uint range can be utilized for all sorts of flags
#define CA_DIMENSIONS_MINIMUM				50			//so that nxyz / nxy is not resulting in wrap around and inside/outside bounds check not run in wraparound
#define MINIMUM_REGION_SIZE					10			//smallest dimension of a thread region

#define CA_ALLOCATION_MAXIMUM				4096000000
#define CA_GRAINIDS_MAXIMUM					4096000000	//not more than 1600^3 disjoint IDs, anyway, who would be as stupid as attempting discretizing aS many grains with 1600^3 cells?
#define NO_GRAIN_ASSIGNED					4096000001	//MK::MUST NOT BE ZERO! otherwise one has to change the indexing and bookkeeping scheme for deformed and rxed grains over mycellgrid
#define I_DONT_KNOW_YET						4096000001
#define NOT_ASSIGNED_YET					4096000001
#define INVALID_CELLASSIGNMENT				4096000003
#define INVALID_COORDINATE					32769
#define CELL_IS_A_PARTICLE					4096001000
#define CURRENTLY_INFECTED					4096002000
#define FULLY_RECRYSTALLIZED				4096003000	//defgid run from [0 to ndefgpool - 1]
#define NUCLEUS_ALREADY_PLACED				4096004000
#define SITE_ALREADY_CONSUMED				4096006000
#define NUCLEUS_NOT_PLACABLE				(-1.0)
#define INVALID_ADDRESS						4294967290
#define UINT32T_MAX							4294967294	//2^32-1 otherwise wrap around but owing to MARKER_TO_IDENTIFY_EXPELLED_SITES better minus another one
#define DIEHL_RXG_OFFSET					1000000000

//caregions communication across their boundaries in a Moore environment
#define NUMBER_OF_NEIGHBORS					26
#define MYSELF								1
#define THIS_REGION							0
#define UNKNOWN_NEIGHBOR					42949672
#define UNKNOWN_POSITION					42949673

#define IDS_AND_LIMITS						7 //ID, xmi,xmx,ymi,ymx,zmi,zmx
#define THE_IDS								0
#define THE_XMIN							1
#define THE_XMAX							2
#define THE_YMIN							3
#define THE_YMAX							4
#define THE_ZMIN							5
#define THE_ZMAX							6

//halo interaction
#define HALO_NOT_EXISTENT					1048577 //MK::it is very unlikely to have as many threads per MPI process on machines in the next years
#define NO_HALO_CELL_INFECTED				0
#define UNKNOWN_ORIGIN						32769
#define UNKNOWN_HALO						42
#define HALOCELL_NEVER_VISITED				0
#define HALOCELL_ALREADY_INFECTED			1

//omp solveRXGROWTH statii
#define OMP_OPERATIVE						0x01
#define OMP_COMPLETED						0x02

//datatype range limits
#define UCHAR_RANGE_MIN						0
#define UCHAR_RANGE_MAX						255
#define INTEGER_RANGE_MAX					2147483647
#define SHORT_MAX							32767

//unsigned char RGBA color model constants
//RGB MAX MAX MAX is white whereas MIN MIN MIN is black
#define RED									0
#define GREEN								1
#define BLUE								2
#define ALPHA								3
#define RGB									3
#define RGB_MAX								255
#define RGB_MIN								0

//implementation relevant defines memory guard
#define BYTE2MEGABYTE						1048576


//definition of ideal texture component / ideal / standardlagen
#define RANDOM_ORIENTATION					-1
#define NOT_WITHIN_GRIDRESU					-1
#define MAXIMUM_DISJOINT_ORIS				4096000000
#define CATEGORIZED_AS_RANDOM				4096000001
#define UNKNOWN_ORIENTATION					4096000001
#define INVALID_NUCSITE						4096000001

//limiting values defining which potential simulation setups can be realized
#define DEFAULT_MAX_NCAS					1048576		//corresponding to 32768 process each of which hosting 32 caHdl
#define DEFAULT_MAX_DEFORMEDGRAINS			65536
#define DEFAULT_MAX_DEFORMEDGRAINS_EBSD		1048576		//##MK::should be sufficient
#define DEFAULT_MAX_RXGRAINS				1048576
#define DEFAULT_MAX_INPUTDATALEN			4096
#define MIN_DISCRETIZATION_GRAIN			14137		//4/3_PI_10^3  14137		//4/3_PI_15^3 in cells

//nucleation modeling
#define GBNUCLEATION_NO						0x01
#define GBNUCLEATION_YES_PHYS_LAGBHAGB		0x02
#define GBNUCLEATION_YES_PHYS_ONLYHAGB		0x03
#define GBNUCLEATION_YES_PICKRND			0x04
#define CSRNUCLEATION_NO					0x05
#define CSRNUCLEATION_YES_ENFORCED			0x06
#define CSRNUCLEATION_YES_PICKRND			0x07
#define CLUSTERNUCLEATION_NO				0x08
#define CLUSTERNUCLEATION_YES_PICKRND		0x09
#define TINCUB_SITESATURATION				0x10
#define TINCUB_TIMEDEPENDENT				0x11
#define NUCSITES_STILL_SOME_FREE			0x12
#define NUCSITES_ALL_EXHAUSTED				0x13

//further nucleation models
#define CSRNUCLEATION_YES_DIEHL				0x22
#define CSRNUCLEATION_DIEHL_RANDOMSO3		0x23
#define CSRNUCLEATION_DIEHL_RANDOMDEFORI	0x24
#define CSRNUCLEATION_DIEHL_SCATTEREXISTENT	0x25

//recovery modeling
#define RECOVERY_NO							0x01
#define RECOVERY_NES_VACCOREDIFF			0x02
#define RECOVERY_NES_SOLUTEDRAG				0x03
#define RECOVERY_NES_MICHALAKPAXTON			0x04

//dispersoid drag
#define DISPERSOIDDRAG_NO					0x01
#define DISPERSOIDDRAG_CONSTANT				0x02
#define DISPERSOIDDRAG_TIMEDEP				0x03

#define MOBILITYMODEL_SEBALDGOTTSTEIN		0x01
#define MOBILITYMODEL_ROLLETTHOLM			0x02

#define RENDERING_MSNO						0x01
#define RENDERING_MS2D						0x02
#define RENDERING_MS3D						0x03

#define RENDERING_BOUNDARIES_NO				0x01
#define RENDERING_BOUNDARIES_YES			0x02

#define RENDERING_COLOR_UNKNOWN				0x00
#define RENDERING_COLOR_GRAINID				0x01
#define RENDERING_COLOR_IPFZ				0x02
#define RENDERING_COLOR_RHO					0x03

#define RENDERING_FILEFORMAT_RAW			0x01
#define RENDERING_FILEFORMAT_HDF5			0x02

#define OUTPUT_LOGBND_NO					0x01
#define	OUTPUT_LOGBND_YES					0x02

#define OUTPUT_RXFRONTSTATS_NO				0x01
#define OUTPUT_RXFRONTSTATS_YES				0x02

#define OUTPUT_THREADPROFILING_NO			0x01
#define OUTPUT_THREADPROFILING_YES			0x02

#define OUTPUT_SINGLEGRAIN_NO				0x01
#define OUTPUT_SINGLEGRAIN_ASCII			0x02
#define OUTPUT_SINGLEGRAIN_BINARY			0x03

#define OUTPUT_DAMASK_GEOMETRYFILE_NO		0x00
#define OUTPUT_DAMASK_GEOMETRYFILE_YES		0x01

#define OUTPUT_SEMEBSD_NO					0x00
#define OUTPUT_SEMEBSD_YES					0x01

#define	DEFAULT_SEMEBSD_MODE				0x88

#define DATA_NOT_LOGGED						2147483648

//junction skeletonization
#define OUTPUT_JUNCTIONS_GBTHREADCOLORING	0x01
#define OUTPUT_JUNCTIONS_GBFACES			0x02
#define OUTPUT_JUNCTIONS_TRIJUNCTIONS		0x03
#define OUTPUT_JUNCTIONS_HIGHERORDER		0x04

#define VONNEUMANN_DISJOINT_SEEDS			6
#define INVALID_SEED						4294967293


//percolation analysis constants
#define PERCOLATION_ANALYZE_NO				0x00
#define PERCOLATION_ANALYZE_YES_NOSZDISTR	0x01
#define PERCOLATION_ANALYZE_YES_SZDISTR		0x02

//deformation structures
#define CUBOID_DEFMS						0x01
#define POISSONVORONOI_DEFMS				0x02
#define CPFEMDAMASK_DEFMS					0x03
#define SEMEBSD_2D3D_COLUMNSINGLE			0x04
#define SEMEBSD_2D3D_COLUMNSTACK			0x05
#define SEMEBSD_2D3D_COLUMNSHIFT			0x06
#define MASTER_CSR_NPOINTS					1000000
#define XMI									0
#define XMX									1
#define YMI									2
#define YMX									3
#define ZMI									4
#define ZMX									5


/*
//defaults for the load calculation analysis
#define DEFAULT_MEMPERCELL_DEF				4			//int for each cell
#define DEFAULT_MEMPERCELL_RX				24			//active cells
#define DEFAULT_MEMPERGRAIN_LOG				8			//double for cell volume
#define DEFAULT_MEMPERCA_OVERHEAD			20480		//20MB/CA
*/

//defines for deformation microstructure construction
#define DEFAULT_NTSTEPSMAX_DEFSYNTH			1000000


//grain boundary tracking management
#define MYHASH(a,b)							( (((uint_fast64_t) max(a,b)) << 32) | ((uint_fast64_t) min(a,b)) )
#define STDLEN								16	//has to be positive
#define DEFAULT_ALLOCATE_TJ					((1024)*(1024))
#define DEFAULT_ALLOCATE_HJ					((32)*(1024))
#define BNDBIN_CACHE_MULTIPLIER				2
#define REPEAT_MAX_DRAWING					1000
#define GBDISORI_NOT_REQUIRED_YET			(0.0)
#define DEFAULT_VOXEL2FACE_REFERENCES			1000
#define EMPTY								0
#define THEFIRSTONE							0
#define MARKER_TO_IDENTIFY_EXPELLED_SITES	1 //MK::MUST NOT BE ZERO!!!!
#define POSITION_UNKNOWN					4096333333
#define DEFAULT_GBPLOTTING_CACHESIZE		((1024)*(1024)*(1024)) //1GB by default
//grain boundary tracking fast
#define DEFAULT_ALLOC_VOXELPAIRS			(0.25)
#define DEFAULT_REALLOC_VOXELPAIRS			(0.05)
#define DEFAULT_NSTANDARDLAGEN				30
#define GBNAME_UNKNOWN						4097666666

//defines for CELL MANAGEMENT
#define NOTHING_TO_RECYCLE					0
#define	ACTIVE								true
#define INACTIVE							false
#define	CELLAPPEND							1 //everything not 0 is true
#define CELLRECYCLE							0 //is false by definition
#define NO_INFECTION						(0.0) //recrystallized fraction
#define ALMOST_FULLY_INFECTED				(1.0)

#define	PLOTBOUNDARIES_DEFORMED				1

#define INFECTOR_DEEP_IN_THE_CUBE			0x88
#define INFECTOR_CLOSE_TO_THE_DOMAINWALL	0x22

//by analyzing the microstructural path an estimate that in a front tracking approach on the fraction of 0.05-0.15 of all voxel
//are in one timestep participating in the transformation, however on average only maxfillperstep of them become fully transformed in that timestep
//so it is the objective of the code to make use of this fact that most ACTIVE cells are not able to induce a transformation of neighboring cells
#define MINLENGTH_ACTIVELIST				1000
#define MINLENGTH_RECYCLIST					1000
#define MINLENGTH_FULLSEEDLIST				1000


//raw data plotting constants
#define	RAW_PARTICLE						0
#define RAW_INFECTED						1


//rediscretization of ensemble information
#define REDISCR_TIME_EQUIDISTANT			0
#define REDISCR_TIME_LOGDISTANT				1

//defragmentation
#define L1_CACHE_LENGTH						((32)*(1024)/(8)/(32))

//lodepng
#define REDCHAN								0
#define GREENCHAN							1
#define BLUECHAN							2
#define ALPHACHAN							3
#define	DEFAULT_ZSECTIONING_ZPOS			(0.5)
#define BLACK								0
#define WHITE								255

//plotting mode
#define COLORIZE_RX_LEAVE_DEFORMED_BLACK	0x01
#define COLORIZE_DEF_LEAVE_RX_BLACK			0x02
#define COLORIZE_RX_RHOGREYSCALE_DEFORMED	0x03 //lowest dislocation density is white, black, highest dislocation density

//IO relevant
#define MAX_NUMBER_OF_NUCLEI_MPI_BUFFERED	((10) * (1000000)) //*(8+8+4+4) byte internal buffer necessary on master!
#define	CHARBUFFER_SIZE						256
#define THREAD_PROFILING_NPROPS_SYNTH		14
#define THREAD_PROFILING_NPROPS_GROWTH		13

//parallelization
//MPI
#define MASTER								0
//OpenMP
#define DEBUG_UTILIZE_ONLY_ONE_THREAD		0
#define ALL_THREADS_WORK_ONONECA			1


//HDF5 library
#define HDF5_MAXIMUM_SINGLE_WRITE			((1600)*(1024)*(1024)) //1.2GB
#define MPIIO_MAXIMUM_SINGLE_WRITE			((50)*(1024)*(1024)) //50MB

#endif
