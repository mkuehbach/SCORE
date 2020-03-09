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


#ifndef __SCORE_MPIIOTYPES_H_INCLUDED__
#define __SCORE_MPIIOTYPES_H_INCLUDED__

#include <mpi.h>

typedef struct
{
	long JobID;
	long Iterations;
	double tend;
	double MPIWTimeInitialization;
	double MPIWTimeSpendDefMS;
	double MPIWTimeSpendGBDetection;
	double MPIWTimeSpendNucleation;
	double MPIWTimeSpendGrowthSim;
	double MPIWTimeSpendFinalIO;
} MPI_IO_CAProfilingInfoData;

//##MK::further optimization necessary
typedef struct
{
	unsigned int nx;
	unsigned int ny;
	unsigned int nz;
	unsigned int nxyz;

	unsigned int regx;		//threaded execution
	unsigned int regy;
	unsigned int regz;
	unsigned int nboundarycells;

	unsigned int ndgrseeds;
	unsigned int nrxgrseeds;
	unsigned int ndefmicrotexture;
	unsigned int nnucleimicrotexture;

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

#endif
