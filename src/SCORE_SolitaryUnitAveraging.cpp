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


void ensembleHdl::postprocess_initrediscrtimes( void )
{
	//acceptable time frame investigated
	//double dtworld = ensRediscretization.tensmax - ensRediscretization.tensmin;


	//construct rediscretized time-annealing scheme
	ensRediscretization.ensRediscrTime = NULL;
	ensRediscretization.ensRediscrTime = new double[ensRediscretization.nslots];
	QUICKASSERT( ensRediscretization.ensRediscrTime != NULL );
	ensMemGuard = ensMemGuard + (ensRediscretization.nslots * sizeof(double));


	if ( ensRediscretization.strategy == REDISCR_TIME_EQUIDISTANT ) {
		double ttmin = ensRediscretization.tensmin;
		double ttmax = ensRediscretization.tensmax;
		uint32_t nmax = ensRediscretization.nslots;

		double dt = (ttmax - ttmin ) / ((double) (nmax - 1));
		QUICKASSERT ( dt > REDISCR_DTMIN );
		
		for ( uint32_t n = 0; n < nmax; ++n ) { //points to end of the interval
			ensRediscretization.ensRediscrTime[n] = ttmin + ((double) (n + 1) * dt);

			//cout << "->RediscrTime;myRank;n;value\t\t" << myRank << "\t" << n << "\t" << ensRediscretization.ensRediscrTime[n] << endl;
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


void ensembleHdl::postprocess_write_mycas_mpislow( void )
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
		log_profiling_file << "MyRank;JobID;tend;Iterations;MPIWTimeInitialization;MPIWTimeDefMS;MPIWTimeGBDetection;MPIWTimeNucleation;MPIWTimeGrowth;MPIWTimeFinalIO;MPIWtimeSinceMasterExitedMPIInit;CAnx;CAny;CAnz;CAnxyz;ThreadingRegionsX;ThreadingRegionsY;ThreadingRegionsZ;CASvDeformed;CAStoredEnergy/in1E+14;CAnDefSeeds;CAnRXG;CAnIdeal;IdealComponentsDeformation;IdealComponentsNucleiIdentified\n";

		for ( uint32_t myca = 0; myca < thisrank_nids; myca++ ) {
			prof_tio = MPI_Wtime();

			log_profiling_file << myRank << ";" << myCAProfiler[myca].JobID << ";" << myCAProfiler[myca].tend << ";" << myCAProfiler[myca].Iterations << ";" << setprecision(6) << myCAProfiler[myca].MPIWTimeInitialization << ";" << setprecision(6) << myCAProfiler[myca].MPIWTimeSpendDefMS << ";" << setprecision(6) << myCAProfiler[myca].MPIWTimeSpendGBDetection << ";" << setprecision(6) << myCAProfiler[myca].MPIWTimeSpendNucleation << ";" << setprecision(6) << myCAProfiler[myca].MPIWTimeSpendGrowthSim << ";" << setprecision(6) << myCAProfiler[myca].MPIWTimeSpendFinalIO;
			log_profiling_file << ";" << setprecision(6) << (prof_tio - prof_t0) << ";" << myCAPhysics[myca].nx << ";" << myCAPhysics[myca].ny << ";" << myCAPhysics[myca].nz << ";" << myCAPhysics[myca].nxyz << ";" << myCAPhysics[myca].regx << ";" << myCAPhysics[myca].regy << ";" << myCAPhysics[myca].regz << ";" << myCAPhysics[myca].nboundarycells << ";" << setprecision(8) << myCAPhysics[myca].storedenergy << ";" << myCAPhysics[myca].ndgrseeds << ";" << myCAPhysics[myca].nrxgrseeds << ";" << myCAPhysics[myca].ndefmicrotexture;

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
				bucket[i].MPIWTimeInitialization = this->myCAProfiler[i].MPIWTimeInitialization;
				bucket[i].MPIWTimeSpendDefMS = this->myCAProfiler[i].MPIWTimeSpendDefMS;
				bucket[i].MPIWTimeSpendGBDetection = this->myCAProfiler[i].MPIWTimeSpendGBDetection;
				bucket[i].MPIWTimeSpendNucleation = this->myCAProfiler[i].MPIWTimeSpendNucleation;
				bucket[i].MPIWTimeSpendGrowthSim = this->myCAProfiler[i].MPIWTimeSpendGrowthSim;
				bucket[i].MPIWTimeSpendFinalIO = this->myCAProfiler[i].MPIWTimeSpendFinalIO;

				//physical information
				physbucket[i].nx = this->myCAPhysics[i].nx;
				physbucket[i].ny = this->myCAPhysics[i].ny;
				physbucket[i].nz = this->myCAPhysics[i].nz;
				physbucket[i].nxyz = this->myCAPhysics[i].nxyz;
				physbucket[i].regx = this->myCAPhysics[i].regx;
				physbucket[i].regy = this->myCAPhysics[i].regy;
				physbucket[i].regz = this->myCAPhysics[i].regz;
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

			MPI_Send( nuctexture, thisrank_nids*nidealcomponents, MPI_LONG, MASTER, SENDMESRXTEX, MPI_COMM_WORLD );

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

				log_profiling_file << r << ";" << bucket[rca].JobID << ";" << bucket[rca].tend << ";" << bucket[rca].Iterations << ";" << setprecision(6) << bucket[rca].MPIWTimeInitialization << ";" << setprecision(6) << bucket[rca].MPIWTimeSpendDefMS << ";" << setprecision(6) << bucket[rca].MPIWTimeSpendGBDetection << ";" << setprecision(6) << bucket[rca].MPIWTimeSpendNucleation << ";" << setprecision(6) << bucket[rca].MPIWTimeSpendGrowthSim << ";" << setprecision(6) << bucket[rca].MPIWTimeSpendFinalIO;
				log_profiling_file << ";" << setprecision(6) << (prof_tio - prof_t0) << ";" << physbucket[rca].nx << ";" << physbucket[rca].ny << ";" << physbucket[rca].nz << ";" << physbucket[rca].nxyz << ";" << physbucket[rca].regx << ";" << physbucket[rca].regy << ";" << physbucket[rca].regz << ";" << physbucket[rca].nboundarycells << ";" << physbucket[rca].storedenergy << ";" << physbucket[rca].ndgrseeds << ";" << physbucket[rca].nrxgrseeds << ";" << physbucket[rca].ndefmicrotexture;
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
			bucket[i].MPIWTimeInitialization = this->myCAProfiler[i].MPIWTimeInitialization;
			bucket[i].MPIWTimeSpendDefMS = this->myCAProfiler[i].MPIWTimeSpendDefMS;
			bucket[i].MPIWTimeSpendGBDetection = this->myCAProfiler[i].MPIWTimeSpendGBDetection;
			bucket[i].MPIWTimeSpendNucleation = this->myCAProfiler[i].MPIWTimeSpendNucleation;
			bucket[i].MPIWTimeSpendGrowthSim = this->myCAProfiler[i].MPIWTimeSpendGrowthSim;
			bucket[i].MPIWTimeSpendFinalIO = this->myCAProfiler[i].MPIWTimeSpendFinalIO;

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
		log_profiling_file << "MyRank;JobID;tend;Iterations;MPIWTimeInitialization;MPIWTimeDefMS;MPIWTimeGBDetection;MPIWTimeNucleation;MPIWTimeGrowth;MPIWtimeSinceMasterExitedMPIInit;CAnx;CAny;CAnz;CAnxyz;CASvDeformed;CAStoredEnergy/in1E+14;CAnDefSeeds;CAnRXG;CAnIdeal;IdealComponentsDeformation;IdealComponentsNucleiIdentified\n";

		for ( uint32_t myca = 0; myca < thisrank_nids; myca++ ) {
			prof_tio = MPI_Wtime();

			log_profiling_file << myRank << ";" << myCAProfiler[myca].JobID << ";" << myCAProfiler[myca].tend << ";" << myCAProfiler[myca].Iterations << ";" << setprecision(6) << myCAProfiler[myca].MPIWTimeInitialization << ";" << setprecision(6) << myCAProfiler[myca].MPIWTimeSpendDefMS << ";" << setprecision(6) << myCAProfiler[myca].MPIWTimeSpendGBDetection << ";" << setprecision(6) << myCAProfiler[myca].MPIWTimeSpendNucleation << ";" << setprecision(6) << myCAProfiler[myca].MPIWTimeSpendGrowthSim << ";" << setprecision(6) << myCAProfiler[myca].MPIWTimeSpendFinalIO;
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

				log_profiling_file << r << ";" << bucket_from_all[nids_per_rank_cumul[r]+rca].JobID << ";" << bucket_from_all[nids_per_rank_cumul[r]+rca].tend << ";" << bucket_from_all[nids_per_rank_cumul[r]+rca].Iterations << ";" << setprecision(6) << bucket_from_all[nids_per_rank_cumul[r]+rca].MPIWTimeInitialization << ";" << setprecision(6) << bucket_from_all[nids_per_rank_cumul[r]+rca].MPIWTimeSpendDefMS << ";" << setprecision(6) << bucket_from_all[nids_per_rank_cumul[r]+rca].MPIWTimeSpendGBDetection << ";" << setprecision(6) << bucket_from_all[nids_per_rank_cumul[r]+rca].MPIWTimeSpendNucleation << ";" << setprecision(6) << bucket_from_all[nids_per_rank_cumul[r]+rca].MPIWTimeSpendGrowthSim << ";" << setprecision(6) << bucket_from_all[nids_per_rank_cumul[r]+rca].MPIWTimeSpendFinalIO;
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
	ensRediscretization.tensmin = 0.0;
	ensRediscretization.tensmax = worldmax_tsimend; //##MK::nslots was = 100
	ensRediscretization.strategy = REDISCR_TIME_EQUIDISTANT;

	postprocess_initrediscrtimes();
}


//##MK::the postprocessing functions
void ensembleHdl::postprocess_rediscr_kinetics( void )
{
	//thread based interpolation possible //##pragma but pending
	//##MK::strategy cache-friendly when, as here, each grain is analyzed along all times
	uint32_t nmax = ensRediscretization.nslots;

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

			double theca_nboxvol = theca->myCAGeometry.nboxvol_rdtdnd;

			//MK::how much volume sample for the myCAs at the ensemble level?
			ensemble_myCAs_nboxvol = ensemble_myCAs_nboxvol + theca_nboxvol;

			//reset temporary collector for the theca to 0.0
			for ( uint32_t n = 0; n < nmax; ++n ) { 
				rediscr_myCA_allgrains_vol[n] = 0.0; 
			}

			//how much volume still left deformed?
			uint32_t theca_ndefg = theca->mydefgpool.size();

			//##MK::trivial parallel but potentially false sharing, multithreading on the level of each grain
			for ( uint32_t dg = 0; dg < theca_ndefg; ++dg ) { //##MK::can be OpenMP parallelized as well
				for ( uint32_t n = 0; n < nmax; ++n ) {
					//interpolate volume of the dg grain in the mydefgpool of theca at that rediscretized time slot
					rediscr_myCA_allgrains_vol[n] = rediscr_myCA_allgrains_vol[n] + theca->get_interpCellCount( dg, ensRediscretization.ensRediscrTime[n] );
				}
			}

			//pipe this information to the local ensembleLevel
			//utilizing that always Sum(Vdef) + Sum(Vrxg) == Vnboxvol_rdtdnd holds valid
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
		log_kinetics_file << "Step;RediscrTime/s;SumRediscrVolume--" << (world_allCAs_nboxvol) << "--;X;ln(t);ln(-ln(1-X))\n";

		double _vref, rxfrac;
		QUICKASSERT ( world_allCAs_nboxvol > 0.0 );
		_vref = ((1.0) / world_allCAs_nboxvol);

		for ( uint32_t n = 0; n < nmax; ++n ) {
			rxfrac = rediscr_world_vol[n] * _vref;

			log_kinetics_file << n << ";" << ensRediscretization.ensRediscrTime[n] << ";" << rediscr_world_vol[n] << ";" << rxfrac << ";";

			//##MK::move this to a different cache from which one can later read simply
			if ( ensRediscretization.ensRediscrTime[n] > SMALLEST_TIME_LOGARITMIZE )
				log_kinetics_file << log(ensRediscretization.ensRediscrTime[n]);
			log_kinetics_file << ";";

			if ( (1.0 - rxfrac) > SMALLEST_TIME_LOGARITMIZE )
				log_kinetics_file << log(-1.0 * log( 1.0 - rxfrac ));
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
	uint32_t nmax = ensRediscretization.nslots;
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

			double theca_nboxvol = theca->myCAGeometry.nboxvol_rdtdnd;
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
			//##MK::in principle possible to parallelize but then threadlocal buffer necessary for rediscr_....
			for ( uint32_t dg = 0; dg < theca_ndefg; ++dg ) { //##MK::OpenMP tasking possible, but then each threads needs local copy for closestideal variants and then reduction
				//categorize according to standardlagen, catch RANDOM as the first element
				nideal = theca->myoripool[theca->mydefgpool[dg].caori].closestideal;

				//to which ideal orientation is dg belonging?
				for ( uint32_t n = 0; n < nmax; ++n ) {
					//interpolate volume of the dg grain in the mydefgpool of theca at that rediscretized time slot
					rediscr_myCA_allgrains_oridepvol[(nstandardlagen*n)+nideal] += theca->get_interpCellCount( dg, ensRediscretization.ensRediscrTime[n] );
				}
			}

			for ( uint32_t rxg = 0; rxg < theca_nrxg; ++rxg ) {
				//categorize according to standardlagen, catch random as the first element
				nideal = theca->myoripool[theca->myrxgpool[rxg].caori].closestideal;

				idx = theca_ndefg + rxg;

				for ( uint32_t n = 0; n < nmax; ++n ) {
					rediscr_myCA_allgrains_oridepvol[(nstandardlagen*n)+nideal] += theca->get_interpCellCount( idx, ensRediscretization.ensRediscrTime[n] );
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
			log_texture_file << ";" << id << "-" << RAD2DEG(standardlagen[id].bunge1) << "-" << RAD2DEG(standardlagen[id].bunge2) << "-" << RAD2DEG(standardlagen[id].bunge3);
		}
		log_texture_file << ";SumRediscrVolume=" << (world_allCAs_nboxvol) << endl;

		//pipe information to file
		for ( uint32_t n = 0; n < nmax; ++n ) {
			log_texture_file << n << ";" << ensRediscretization.ensRediscrTime[n];
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
	//uint32_t nmax = ensRediscretization.nslots;
	double now = ensRediscretization.tensmax;

	//each rank collects volume from all his nuclei at ensRediscretization.tensmax in all his automata,
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
		//##MK::initialize parallel region, grain log data are stored threadlocal to the master!
		for ( uint32_t tid = 0; tid < myCAs.size(); tid++ ) {
			for ( uint32_t ca = 0; ca < myCAs[tid].size(); ca++) {
				caHdlP theca = myCAs[tid][ca];

				uint32_t ndefg = theca->mydefgpool.size();
				uint32_t nrxg = theca->myrxgpool.size();
				uint32_t ndisjoint_grains = ndefg + nrxg;

				uint32_t nideal;
				//##MK::trivial parallel reading and linear interpolation of grain volume
				for ( uint32_t rxg = ndefg; rxg < ndisjoint_grains; rxg++) {
					nideal = theca->myoripool[theca->myrxgpool[rxg-ndefg].caori].closestideal;
					rediscr_ens_allnuclei[nid].finalvol = theca->get_interpCellCount( rxg, now );
					rediscr_ens_allnuclei[nid].tincub = theca->myrxgpool[rxg-ndefg].tincub;
					rediscr_ens_allnuclei[nid].ideal = nideal;
					nid++;
				} //because each grain was nucleated
			}
		}
	}
	else if ( this->ensNucleationModel.tincubmodel == TINCUB_TIMEDEPENDENT ) {
		//##MK::initialize parallel region, grain log data are stored threadlocal to the master!
		for ( uint32_t tid = 0; tid < myCAs.size(); tid++ ) {
			for ( uint32_t ca = 0; ca < myCAs[tid].size(); ca++) {
				caHdlP theca = myCAs[tid][ca];

				uint32_t ndefg = theca->mydefgpool.size();
				uint32_t nrxg = theca->myrxgpool.size();
				uint32_t ndisjoint_grains = ndefg + nrxg;

				uint32_t nideal;
				//##MK::trivial parallel reading and linear interpolation of grain volume
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
	//##MK::it is very likely that ens_nnuclei is populated but not assured so strictly speaking a barrier is necessary here
	MPI_Barrier( MPI_COMM_WORLD ); //its wasnt before

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
	double finalgsdmin = std::numeric_limits<double>:: max();
	double finalgsdmax = std::numeric_limits<double>:: lowest();

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

			if ( volnorm > SMALLEST_TIME_LOGARITMIZE )	log_finalgsd_file << ";" << setprecision(12) << (log(volnorm)) << ";" << setprecision(12) << probability << "\n";
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
	//uint32_t nmax = ensRediscretization.nslots;
	double now = ensRediscretization.tensmax;

	//each rank collects volume from all his nuclei at ensRediscretization.tensmax in all his automata,
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
	if ( ensemble_all_nnuclei > MAX_NUMBER_OF_NUCLEI_MPI_BUFFERED ) { return false; }

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


double caHdl::get_interpCellCount( uint32_t localid, double when )
{
	//##MK::further optimization

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

	//in the earlier version scanning was performed on mygrainevolution[i].localtime with large cache elements and hence much more cacheloads until finding i
	uint32_t nmygrainevo = mygrainevolution.size();
	//while ( mygrainevolution[i].localtime < when && i < nmygrainevo ) { //framing the interval ([i-1].localtime, [i].localtime)
	while ( mygrainevolutionlocaltime[i] < when && i < nmygrainevo ) {
		i++;
	}

	//MK::implicitly assuming the microstructure unchanged
	//MK::lift load from branch predictor by planning for the most likely case
	if ( i > 0 && i < nmygrainevo ) {
		//###disjunctness of rediscretization scheme has been assured by ensRediscretization
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

	return interpCellCount;
}

#endif