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
#ifndef __SCORE_DEFS_PHYSICS_H_INCLUDED__
#define	__SCORE_DEFS_PHYSICS_H_INCLUDED__

//MK::here all default values of physical relevance should be placed!


//physical constants
#define kboltzman							(1.3806488e-23)
#define TOFFSET								(273.15)				//degree Celsius into Kelvin
#define echarge								(1.602176565e-19)		//Coulomb
#define KUHLMANNFINALRHO					(1e10)
#define pi									(3.1415926535897932384626433832795)


//conversion of SI units for specific order of magnitude community standards
#define MICRON2METER						(1.0e-6)
#define SCALING_DEFAULT_DISDENS				(1.0e14)				//m^2


//NUMERICAL ACCURACY
//integrator accuracy limits to assure reasonable numerical solutions during transient annealing
#define RHO_DEFORMED_MATERIAL_MAX			(1.0e18)				//m^-2
#define RHO_RECRYSTALLIZED_MATERIAL			(1.0e11)				//m^-2
#define MAXFILLPERSTEP						(0.1)
#define SMALLEST_TIME_LOGARITMIZE			(1.0e-9)
#define VERYSMALLTIMESTEP					(1e-6)					//s - do not use smaller timesteps than this
#define INITIAL_DELTAT						(1e-6)					//micro s
#define INFINITE							(1e36)					//s

#define FGEOFACE							(1.0000)
#define FGEOEDGE							(0.7957)
#define FGEODIAG							(0.6494)

//thresholds
#define RXFRACTION_THRESHOLD				(0.50)	//at which infection state do we consider a cell as recrystallized?


//discretization to time integration scheme
#define SMALL_HEATRATE						(0.000277777777777777)	//1 K/ 3600s
#define SMALL_RECOVERYRATE					(2.77e10)				// 10^14 1/m^2/3600s
#define SMALL_DRAGGINGRATE					(0.027777777)			// 100Pa/3600s
#define SMALL_HEAT							(1.0)					//K
#define SMALL_RHO							(1.0e11)				//1/m^2 ~ 0.01 rho min
#define SMALL_ZENERFORCE					(100.0)					//Pa
#define SMALL_NUMBEROFNUCLEI				(1)
#define SMALL_NUCLEATIONRATE				(2.77e-4)				//1 nuc/3600s
#define MINIMUM_RAYLEIGH_SIGMA				(1e-3)					//x/sigma^2*exp(-0.5x^2/sigma^2)
//#define SMALL_DISTANCE									this is the dcell
#define SMALL_VELOCITY						(2.7e-13)				//1nm/3600s


//Zener drag model defaults
#define DEFAULT_ZENER_ALPHA					((3.0)/(2.0))			//classical geometrical constant from Zener/Smith' model
#define DEFAULT_ZENER_GAMMA					(0.324)					//J/m^2 grain boundary energy Aluminium according to Murr et. al.


//definition of ideal texture component / ideal / standardlagen
#define MAX_FCC								(1.099)					//MacKenzie..., 62.8/180*_PI_
#define RESOLUTION_SO3GRID					((0.5)/(180.0)*(pi))	// ori degree raster
#define MAXDISORI_TO_40DEG111				((10.0)/(180.0)*(pi))
#define MAXDISORI_LAGB2HAGB					((15.0)/(180.0)*(pi))
#define DEFAULT_SCATTER						(10.0)
#define BUNGE_MIN							((0.0)/(180.0)*(pi))
#define BUNGE_MAX							((360.0)/(180.0)*(pi))


//default values in the construction phase
#define DEFAULT_SIMID						0
#define DEFAULT_CELLSIZE					(1e-6)
#define DEFAULT_SHEARMODULUS				(2.7e10)		//RT pure aluminium
#define DEFAULT_BURGERSVECTOR				(2.86e-10)	//RT pure aluminium
#define DEFAULT_PUREALU_C					(1.0)
#define DEFAULT_PUREALU_TMELT				(933.0)

#define DEFAULT_SMALL_NUMBER				(1e-6)
#define DEFAULT_LAGB_HACT					(1.30)		//eV
#define DEFAULT_HAGB_HACT					(1.20)
#define DEFAULT_GS_HACT						(1.10)
#define DEFAULT_DEFGSIZE					(1e-4)		//m
#define DEFAULT_MAXFILLPERSTEP				(0.1)
#define DEFAULT_MINFILLIN					(0.01)		//smaller is unnecessary fine sampling
#define DEFAULT_MAXFILLIN					(0.2)		//larger too coarse behavior

#define DEFAULT_RELCELLCACHING				(0.10)
#define DEFAULT_TRANSRELCELL				(1.0)
#define DEFAULT_PRNG_SEED					-46356
#define DEFAULT_DELTATIME					(1e-6)		//s
#define DEFAULT_XMAX						(1.0)
#define DEFAULT_TMAX						(360000.0)	//s
#define DEFAULT_NMAX						1000000
#define DEFAULT_NUCDENSITY_CSR				1
#define DEFAULT_CELLSIZE_MIN				(1e-7)		//smaller violates mesoscale behavior of atomic defects!
#define DEFAULT_CELLSIZE_MAX				(1e-5)		//larger cells not expected because of technical grain size barely orders of magnitude of 10micron
#define DEFAULT_MINRELCELL					(0.05)
#define DEFAULT_MAXRELCELL					(0.25)
#define DEFAULT_MINTRANSRELCELL				(0.5)
#define DEFAULT_MAXTRANSRELCELL				(2.0)		//according to GNU
#define DEFAULT_MIN_GRAINDISCR				10			// in cells
#define DEFAULT_PMAX						(0.0)		// assuming that their are fast boundaries in the system
#define MAXATTEMPTS_NUCSITE					100000
#define RHOMAX_WELLANNEALED					(1e8)


//grain boundary tracking fast
#define DEFAULT_ALLOC_VOXELPAIRS			(0.25)
#define DEFAULT_REALLOC_VOXELPAIRS			(0.05)
#define DEFAULT_NSTANDARDLAGEN				30


//by analyzing the microstructural path an estimate that in a front tracking approach on the fraction of 0.05-0.15 of all voxel
//are in one timestep participating in the transformation, however on average only maxfillperstep of them become fully transformed in that timestep
//so it is the objective of the code to make use of this fact that most ACTIVE cells are not able to induce a transformation of neighboring cells
#define FULLRXCACHING_ATT_FACTOR			(1.5)
#define FULLRECYCLING_ATT_FACTOR			(3.0)


//rediscretization of ensemble information
#define REDISCR_TIMESLOTS_DEFAULT			1000
#define REDISCR_TIMESLOTS_MAX				10000
#define REDISCR_TIMESLOTS_MIN				100
#define REDISCR_DTMIN						(1.0e-6)


//##MK::first trial for Christian Haases problem
#define CPFEM_NGRAINSX						20
#define CPFEM_NGRAINSY						20
#define CPFEM_NGRAINSZ						20
#define CPFEM_NGRAINSXY						( (CPFEM_NGRAINSX) * (CPFEM_NGRAINSY) )
#define CPFEM_NGRAINSXYZ					( (CPFEM_NGRAINSX) * (CPFEM_NGRAINSY) * (CPFEM_NGRAINSZ) )
#define CPFEM_GRAINSIZEX					(59.0e-6)
#define CPFEM_GRAINSIZEY					(14.0e-6)
#define CPFEM_GRAINSIZEZ					(14.0e-6)
#define	SUBSEQ_VOXELIZATION					0


//MK::grain boundary nucleation number model
#define NO_NUCLEI							0
#define MINIMUM_DRHO						(1.0e1)
#define SCALING_LAMBDA						(10.0)
#define POISSON_CUMSUM_CUTOFF				19
#define POISSON_CUMSUM_TABLE				20				//cutoff+1
#define SCALE_NUCDENSITY					(1.0e-4)		//1/micron^2
#define	LAGB_TO_HAGB_TRANS					((15.0)/(180.0)*(_PI_))

//MK::clustering nucleation model
#define DEFAULT_KMAX_POISSRND				5000 
#define DEFAULT_MASTERMATERN_MAXRAD			(0.1)

//MK::Diehl/Kuehbach conversion of Kernel Average Misorientation in Dislocation density
#define DEFAULT_EBSDSTEPSIZE_MINIMUM		(0.05e-6)			//spatial resolution limit SEM/EBSD meter

#endif
