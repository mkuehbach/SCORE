//MK::SCORE is a software developed by Markus K{\"u}hbach in 2014/2015 with the Institute of Physical Metallurgy and Metal Science, RWTH Aachen University
//for questions and comments please contact markus.kuehbach at rwth-aachen.de, last major revision 20150920
//##MK::marks comments that where the code is in rather quick and dirty,i.e. parts of the code which have been validated operative as is but offer certainly room for further improvements

#include "SCORE_Kernel.h"


//user-interaction passed arguments
#define INPUTFILE				1
#define	ID						2


using namespace std;


int main(int argc, char** argv) {

	int simid = atoi(argv[ID]);
	if ( simid < 1 && simid >= INTEGER_RANGE_MAX ) { cout << "ERROR::Invalid simulation ID to identify files." << endl; return 0; }

	ensembleHdlP myensemble = new ensembleHdl;
	myensemble->simid = simid;

	if( !myensemble->open(READ, argv[INPUTFILE], &(myensemble->score_input)) ) { cout << "ERROR::Reading input file." << endl; }

MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &myensemble->nRanks);
	MPI_Comm_rank(MPI_COMM_WORLD, &myensemble->myRank);
	myensemble->prof_t0 = MPI_Wtime();

//preparation phase
	myensemble->init_ensprng( true ); //true - seed is (-1*myRank)-1, false - all the same default seed
	myensemble->init_mpidatatypes();
	myensemble->init_parameter();
	myensemble->init_processing();
	myensemble->init_defgpool();
	myensemble->init_rxgpool();

	//no Barrier necessary when distributeWorkOnRanksOnRanks separates work without iteratively packing/work partitioning
	myensemble->init_distributeWorkOnRanks();
	myensemble->init_distributeWorkOnStartedThreads();
	//no Barrier necessary because now each process can evolve independently the associated automata

	myensemble->SIMULATE_myCAs();

	//MPIBarrier necessary,before gathering all results again at the local process level
	MPI_Barrier( MPI_COMM_WORLD );
	myensemble->prof_tstartpp = MPI_Wtime();

	myensemble->postprocess_write_mycas(); //MK::either explicit sequential retrieval or gatherv
	//myensemble->postprocess_write_mycas_mpifast();
//cout << myensemble->myRank << " beyond the Barrier!" << endl;

	myensemble->postprocess_init();
	myensemble->postprocess_rediscr_kinetics();
	myensemble->postprocess_rediscr_macrotexture();

	//try to output the final grain size distribution by collecting in parallel from all nodes at once, if it fails, sequential collecting which is much slower
	/*bool gsdfast = false;
	gsdfast = myensemble->postprocess_rediscr_finalgrainsizedistribution_mpifast();
	if ( gsdfast == false ) {
		myensemble->postprocess_rediscr_finalgrainsizedistribution();
		//##MK::implement more efficient hierarchical strategy here, if further performance gains are necessary...
	}*/
	//experience showed that explicit MPI blocking MPI comm was more efficient that collective routines
	//##MK::if in doubt use this function, otherwise use the collective version ...mpifast with the source above
	myensemble->postprocess_rediscr_finalgrainsizedistribution();

	myensemble->destroy_myCAs();

	myensemble->prof_tend = MPI_Wtime();
	if ( myensemble->myRank == MASTER ) {
		cout << "Simulation finished with some new insights into pRX, skal." << endl;
		cout << "MASTER exiting MPI_init prior to starting postprocessing seconds=" << setprecision(6) << (myensemble->prof_tstartpp - myensemble->prof_t0) << endl;
		cout << "MASTER starting postprocessing and prior to MPI_Finalize seconds=" << setprecision(6) << (myensemble->prof_tend - myensemble->prof_tstartpp) << endl;
		cout << "MASTER total time seconds=" << setprecision(6) << (myensemble->prof_tend - myensemble->prof_t0) << endl;
	}
MPI_Finalize();

	fclose(myensemble->score_input);
	delete myensemble;

	return 0;
}


	/*uint32_t nlim = 4294967295;
	uint32_t n = (uint32_t) atof(argv[1]);
	uint32_t nbx = 1625;
	uint32_t nby = 1625;
	uint32_t nbz = 1625;
	uint32_t nbxy = nbx * nby;
	uint32_t nbxyz = nbx * nby * nbz;
	cout << "Input was = " << atoi(argv[1]) << " uint32_t = " << n << " nlim = " << nlim << " nbxy = " << nbxy << " nbxyz = " << nbxyz << endl;
	//for ( uint32_t n = 0; n < nlim; ++n ) {
		uint32_t iz = n / nbxy;
		uint32_t rem = n - (nbxy * iz);
		uint32_t iy = rem / nbx;
		uint32_t ix = rem - (nbx * iy);
	cout << "ix\t\tiy\t\tiz----rem\t\t\t" << ix << "\t\t" << iy << "\t\t" << iz << "----" << rem << endl;
	cout << "sizeof(uint32_t) = " << sizeof(uint32_t) << endl;
	return 0;*/

	/*if (myensemble->myRank == MASTER) { 
		cout << "SCORe - Statistical Cellular Operator Ensemble Model for Recrystallization." << endl;
		cout << "defcell=" << sizeof(defcell) << endl;
		cout << "cell=" << sizeof(cell) << endl;
		cout << "rxg=" << sizeof(rxg) << endl;
		cout << "defg=" << sizeof(defg) << endl;
		cout << "ori=" << sizeof(ori) << endl;
		cout << "qsymmfcc=" << sizeof(qsymm_fcc) << endl;
		cout << "cadefg=" << sizeof(cadefg) << endl;
		cout << "carxg=" << sizeof(carxg) << endl;
		cout << "ideal=" << sizeof(ideal) << endl;
		cout << "loginfophys=" << sizeof(loginfo_ca_physics) << endl;
		cout << "defgseed=" << sizeof(defgseed) << endl;
		cout << "cellsBndFast=" << sizeof(cellsBndFast) << endl;
		cout << "bndFaceFast=" << sizeof(bndFaceFast) << endl;
		cout << "bndColumnFast=" << sizeof(bndColumnFast) << endl;
		cout << "loginfogrevo=" << sizeof(loginfo_grainevo_ca) << endl;
		cout << "loginfo_rxfront=" << sizeof(loginfo_rxfrontstats_ca) << endl;
		cout << "ensHdl=" << sizeof(ensembleHdl) << endl;
		cout << "caHdl=" << sizeof(caHdl) << endl;
	}*/
