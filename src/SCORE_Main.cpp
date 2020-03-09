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


#include "SCORE_Kernel.h"

using namespace std;

//MK::dont forget to initialize proper environment by setting OMP_NUM_THREADS and 
//utilize process + thread pinning via $export OMP_NUM_THREADS=<n> and $export KMP_AFFINITY=verbose,scatter
//execute program via $mpiexec -n <numproc> scorehybrid_rhcs <ID> <INPUTFILE> optional:<EBSDMAPPING> <PRNGSEED>
#define	ID						1
#define INPUTFILE				2
#define EBSDMAPPING				3
#define PRNGSEED				4

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>




void die( const char* message, bool inmpi ) {
	cout << "ERROR::" << message << endl;
	if ( inmpi == true ) 
		MPI_Finalize();
}

void versioninfo( int myrank, int nranks ) {
		cout << "SCORE starting rank " << myrank << " within " << nranks << endl;
		cout << "SCORE v" << VERSION_MAJOR << "." << VERSION_MINOR << "." << VERSION_REVISION;
#ifdef IS_BETAVERSION_YES
		cout << "beta" << endl;
#else
		cout << endl;
#endif
		cout << "SCORE the coordinate system is right-handed with x pointing right and ||RD, y pointing inwards and ||TD, and z pointing upwards and ||ND" << endl;
#ifdef USE_QUATLIB_DEGRAEF
		cout << "SCORE the mathlibrary is the Rowenhorst/Rollett/Konijnenberg/deGraef et. al. library" << endl;
#endif
#ifdef USE_QUATLIB_DAMASK
		cout << "SCORE the mathlibrary is the DAMASK library" << endl;
#endif
#ifdef USE_QUATLIB_IMM
		cout << "SCORE the mathlibrary is the IMM library" << endl;
#endif
}

bool checkworkingdirectory( string dirname ) {
	//create an organizational folder structure specific wf simulation if not already existing
	int status = 0;
	struct stat st = {0};
	if (stat(dirname.c_str(), &st) == -1) { //if not existent attempt to generate simulation folder 
		status = mkdir( dirname.c_str(), 0777);
		//perror("mkdir");
		if ( stat(dirname.c_str(), &st) == -1) { return false; } //still not there go out
	}

	string pngs = dirname + "/png"; //microstructure snapshots
	if (stat(pngs.c_str(), &st) == -1 ) {
		status = mkdir( pngs.c_str(), 0777 );
		if ( stat(pngs.c_str(), &st) == -1) { return false; }
	}

	string cmodel1 = pngs + "/rxipfz_dgblack";
	if (stat(cmodel1.c_str(), &st) == -1 ) {
		status = mkdir( cmodel1.c_str(), 0777 );
		if ( stat(cmodel1.c_str(), &st) == -1) { return false; }
	}

	string cmodel2 = pngs + "/dgipfz_rxblack";
	if (stat(cmodel2.c_str(), &st) == -1 ) {
		status = mkdir( cmodel2.c_str(), 0777 );
		if ( stat(cmodel2.c_str(), &st) == -1) { return false; }
	}

	string cmodel3 = pngs + "/dgrhogry_rxipfz";
	if (stat(cmodel3.c_str(), &st) == -1 ) {
		status = mkdir( cmodel3.c_str(), 0777 );
		if ( stat(cmodel3.c_str(), &st) == -1) { return false; }
	}

	string ang = dirname + "/ang"; //virtual SEM/EBSD
	if (stat(ang.c_str(), &st) == -1 ) {
		status = mkdir( ang.c_str(), 0777 );
		if ( stat(ang.c_str(), &st) == -1) { return false; }
	}

	string res = dirname + "/res"; //simulation results
	if (stat(res.c_str(), &st) == -1 ) {
		status = mkdir( res.c_str(), 0777 );
		if ( stat(res.c_str(), &st) == -1) { return false; }
	}

	string prof = dirname + "/prof"; //profiling
	if (stat(prof.c_str(), &st) == -1 ) {
		status = mkdir( prof.c_str(), 0777 );
		if ( stat(prof.c_str(), &st) == -1) { return false; }
	}

	//##MK::add further folders to organize data storage
	return true;
}


int main(int argc, char** argv) {
	//set environment and working folders probe if input files are there and load them
	ensembleHdlP myensemble = new ensembleHdl;
	myensemble->simid = strtoul(argv[ID], NULL, 0 );
	string dirnm = "SCORE_" + std::to_string(myensemble->simid);
	myensemble->simwrkdir = dirnm;
	if ( checkworkingdirectory( dirnm ) == false ) { die("Unable to generate working directories!", false ); }
	cout << "SCORE initialized working directory " << dirnm.c_str() << " ..." << endl;

	if( !myensemble->open(READ, argv[INPUTFILE], &(myensemble->score_input)) ) { die("Reading input file!", false); }
	cout << "SCORE loaded input data " << argv[INPUTFILE] << " ..." << endl;

	if ( argc > (1+2) ) {
		if ( myensemble->init_ebsdmap_read( argv[EBSDMAPPING] ) == false ) {
			die("EBSDMapping input desired but failed!", true ); delete myensemble; return 0;
		} //if readmap didnt fail data are now loaded!
		else {
			cout << "SCORE loaded SEM/EBSD map " << argv[EBSDMAPPING] << " ..." << endl;
		}
	}

	//MK::Init MPI capable of hybrid MPI/OpenMP even though SCORE currently at no point calls MPI within a parallel region
	int supportlevel_desired = MPI_THREAD_FUNNELED;
	int supportlevel_provided;
	MPI_Init_thread(&argc, &argv, supportlevel_desired, &supportlevel_provided);
	if ( supportlevel_provided < supportlevel_desired ) { die("Insufficient threading capabilities of MPI library!", false); }
	else {
		MPI_Comm_size(MPI_COMM_WORLD, &myensemble->nRanks);
		MPI_Comm_rank(MPI_COMM_WORLD, &myensemble->myRank);
		string helloworld = "Firing up MPI rank " + std::to_string(myensemble->myRank) + " over MPI_COMM_WORLD with " + std::to_string(myensemble->nRanks) + " in total\n"; cout << helloworld;
		if ( myensemble->myRank == MASTER ) {
			versioninfo( myensemble->myRank, myensemble->nRanks );
		}

		//preparation phase
		myensemble->prof_t0 = MPI_Wtime();
		if ( argc > (1+3) )
			if ( atoi(argv[PRNGSEED]) > 0 )
				myensemble->init_ensprng( atoi(argv[PRNGSEED]) );
			else
				myensemble->init_ensprng( true );
		else
			myensemble->init_ensprng( true );

		myensemble->init_mpidatatypes();
		if ( myensemble->init_parameter() == false ) { die("Invalid parameterization and/or input data!", true ); delete myensemble; return 0; }

		myensemble->init_distributeWorkOnRanks();
		myensemble->init_distributeWorkOnThreads();

		//no Barrier necessary because each process can evolve his own SU queue independently from all other ranks
		myensemble->SIMULATE_myCAs();

		//MPIBarrier necessary before gathering all results in the solitary unit postProcessing
		MPI_Barrier( MPI_COMM_WORLD );

		myensemble->prof_tstart = MPI_Wtime();
		//if ( myensemble->experimentInput == false ) //##MK::SU analysis currently only done for summary
			myensemble->PERFORM_SOLITARYUNIT_MODELING_ANALYSIS();

		myensemble->destroy_myCAs();
		myensemble->prof_tend = MPI_Wtime();

		if ( myensemble->myRank == MASTER ) {
			cout << "Simulation finished with some new insights into pRX, skal!\n" << endl;
			cout << "Profiling..." << endl;
			cout << "\t\tMASTER Time between MPI startup up to begin of Postprocessing:" << setprecision(6) << (myensemble->prof_tstart - myensemble->prof_t0) << endl;
			cout << "\t\tMASTER Time for Postprocessing and clean-up " << setprecision(6) << (myensemble->prof_tend - myensemble->prof_tstart) << endl;
			cout << "\t\tMASTER Time Total " << setprecision(6) << (myensemble->prof_tend - myensemble->prof_t0) << endl;
		}
		MPI_Finalize();
	}

	delete myensemble; return 0;
}

/*
	caHdlP aca = NULL;
	aca = new caHdl;
	aca->test_lodepng();
	delete aca;
	return 0;

	//hacky overwrite for debugging
	//if ( argc == 3 ) { //programname and two files
		cout << "Performing only debugging!" << endl;
		ensembleHdlP debugger = new ensembleHdl;
		//debugger->ipf_colormap( atof(argv[1]), atof(argv[1]), atof(argv[1]), atoi(argv[2]) );
		//debugger->ang_file_coloring( argv[1], atoi(argv[2]) );
		//debugger->ang_file_coloring( argv[1], atoi(argv[2]) );
		//debugger->compare_two_ang_files( argv[1], argv[2] );
		delete debugger;
		return 0;
	//}

//proof that global instance of ensembleHdl all feed from the same copy of a globally defined mathLibrary
//therefore the processes will individually generate the same sequence of orientations
//i.e. given the same starting structure they will distribute the same random orientations
	ensembleHdlP debugger = new ensembleHdl;
	int sldes = MPI_THREAD_FUNNELED;
	int slpro;
	MPI_Init_thread(&argc, &argv, sldes, &slpro);
	MPI_Comm_size(MPI_COMM_WORLD, &debugger->nRanks);
	MPI_Comm_rank(MPI_COMM_WORLD, &debugger->myRank);
	unsigned int n = atoi(argv[1]); //generate local chian of random numbers
	MPI_Barrier( MPI_COMM_WORLD );
	double* masterchain = new double[4*n];
	if ( debugger->myRank == MASTER ) {
		debugger->ensmath.setprng( -debugger->myRank*100000 - 1 );
		for ( unsigned int i = 0; i < n; i++ ) {
			double qrnd[4] = {1.0, 0.0, 0.0, 0.0};
			//debugger->rnd_quat_shoemake(qrnd); //utilizing global namespace instance of mathMethods
			debugger->ensmath.rnd_quat_shoemake(qrnd);
			//masterchain[i] = debugger->r.MersenneTwister();
			masterchain[4*i+0] = qrnd[0];
			masterchain[4*i+1] = qrnd[1];
			masterchain[4*i+2] = qrnd[2];
			masterchain[4*i+3] = qrnd[3];
		}
	}
	MPI_Barrier( MPI_COMM_WORLD );
	MPI_Bcast( masterchain, 4*n, MPI_DOUBLE, MASTER, MPI_COMM_WORLD );
	MPI_Barrier( MPI_COMM_WORLD );
	for ( int r = MASTER; r < debugger->nRanks; r++ ) {
		double mydiff = 0.0;
		if ( debugger->myRank == r && debugger->myRank != MASTER) {
			debugger->ensmath.setprng( -debugger->myRank*100000 - 1 );
			for ( unsigned int i = 0; i < n; i++ ) {
				//mydiff = mydiff + fabs(debugger->r.MersenneTwister() - masterchain[i]);
				double qrnd[4] = {1.0, 0.0, 0.0, 0.0};
				//debugger->rnd_quat_shoemake(qrnd);
				debugger->ensmath.rnd_quat_shoemake(qrnd);
				mydiff += fabs(qrnd[0]-masterchain[4*i+0])+fabs(qrnd[1]-masterchain[4*i+1])+fabs(qrnd[2]-masterchain[4*i+2])+fabs(qrnd[3]-masterchain[4*i+3]);
			}
			cout << "Rank" << r << "\t\t" << setprecision(32) << mydiff << endl;
		}
		MPI_Barrier( MPI_COMM_WORLD );
	}
	MPI_Finalize();
	delete [] masterchain; masterchain = NULL;
	delete debugger; debugger = NULL;
	return 0;



*/
