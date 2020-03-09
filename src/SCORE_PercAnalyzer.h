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

#ifndef __SCORE_PERCANALYZER_H_INCLUDED__
#define __SCORE_PERCANALYZER_H_INCLUDED__

#include "SCORE_Defs_Functionality.h"
#include "SCORE_Defs_Physics.h"
#include "SCORE_MPIIOTypes.h"
#include "SCORE_Io.h"
#include "SCORE_Math.h" //thus implicitly including SCORE_Random as well


//MK::Cluster IDs of 0 mark non-recrystallized microstructure
#define UNKNOWN_ID							(0)
#define STDOUT_Z_PROMPT						(20) //to control when to prompt progress in scanning grid layers each -th layer stdout
#define BILLION								1000000000L
#define CLUSTERINDEXING_OFFSET				1		//first cluster is number 1


using namespace std;

/*
class timer
{
public:	
	timer(){};
	~timer(){};
	
	vector<long long unsigned int> howlong;
	vector<string> forwhat;
	
	void profiling_start(){
		clock_gettime( CLOCK_MONOTONIC, &ts );		
	};
	void profiling_end(){
		clock_gettime( CLOCK_MONOTONIC, &te );
	};
	void profiling_diff_mem( const char* where ){ 
		diff = BILLION * (te.tv_sec - ts.tv_sec) + te.tv_nsec - ts.tv_nsec; //in ns
		howlong.push_back( (long long unsigned int) diff );
		forwhat.push_back( where );
	};
	void profiling_diff_io( const char* where ){
		diff = BILLION * (te.tv_sec - ts.tv_sec) + te.tv_nsec - ts.tv_nsec;
		cout << "\t\t" << where << " " << (long long unsigned int) diff << " ns" << endl;
		howlong.push_back( (long long unsigned int) diff );
		forwhat.push_back( where );
	};
	void profiling_log() {
		if ( howlong.size() != forwhat.size() ) {
			cout << "Inconsistency in logging data!" << endl;
			return;
		}
		cout << endl;
		cout << "Profiling in milliseconds (higher accuracy possible...)" << endl;
		cout << "-------------------------------------------" << endl;
		for (unsigned int evts = 0; evts < forwhat.size(); evts++ ) {
			cout << forwhat.at(evts) << "\t\t\t" << ( howlong.at(evts) / 1000 / 1000 ) << endl;
		}
	};
	
	struct timespec ts,te;
	uint64_t diff;
};
*/


class UF
{
	unsigned int* id;
	unsigned int* sz;

public:
	// Create an empty union find data structure with N isolated sets.
	UF(unsigned int N) {
		id = NULL;
		sz = NULL;
		id = new unsigned int[N];
		sz = new unsigned int[N];
		for( unsigned int i=0; i<N; i++ ) {
			id[i] = i;
			sz[i] = 1;
		}
	}
	~UF() {
		delete [] id;
		delete [] sz;
	}
	
	//Return the id of component corresponding to object p.
	unsigned int initialAssgn(unsigned int i) {
		return i;
	}
	unsigned int find_root(unsigned int p)	{
		unsigned int root = p;
		while (root != id[root])
			root = id[root]; //traversal operation because for a root node i == id[i] holds
		while (p != root) { //path compression
			unsigned int newp = id[p];
			id[p] = root;
			p = newp;
		}
		return root;
	}

	//Replace sets containing x and y with their union, MK::was originally void
	unsigned int merge(unsigned int x, unsigned int y) {
		unsigned int i = find_root(x);
		unsigned int j = find_root(y);
		if ( i == j ) return i; //MK::return simply one...
	
		// make smaller root point to larger one
		if (sz[i] < sz[j])	{ //eq class j larger than eq class i
			id[i] = j; 
			sz[j] += sz[i];
			return j;
		}
		else { //eq class j smaller or equal to eq class i
			id[j] = i; 
			sz[i] += sz[j];
			return i;
		}
	}
};


struct percAnalysis
{
	unsigned int nClusterTotal; // = 0;
	unsigned int LargestClusterCnt; // = 0;
	unsigned int LargestClusterID; // = 0;
	double PL; // = 0.0;
	percAnalysis() : nClusterTotal(0), LargestClusterCnt(0), LargestClusterID(0), PL(0.0) {}
};
typedef percAnalysis * percAnalysisP;


struct loginfo_perc
{
	double ProfInitializing;
	double ProfHoshenKopeling;
	double ProfCompactifying;
	double ProfCheckLabeling;
	double ProfCharacterizing;
	double ProfCheckPercolating;
	double ProfPercTotal;
	unsigned int NCluster;
	unsigned int LargestCluster;
	unsigned int Percolation;
	loginfo_perc() : ProfInitializing(0.0), ProfHoshenKopeling(-1.0), ProfCompactifying(-2.0), ProfCheckLabeling(-3.0), ProfCharacterizing(-4.0), ProfCheckPercolating(-5.0), ProfPercTotal(-6.0), NCluster(4), LargestCluster(5), Percolation(666) {}
};
typedef loginfo_perc * loginfo_percP;



struct compactLabel
{
	unsigned int oky;
	unsigned int nky;
};
typedef compactLabel * compactLabelP;


class percAnalyzer
{
	unsigned int* id;
	UF* finder;

	unsigned int N; // = 0;
	unsigned int NN; // = SQR(N);
	unsigned int NNN; // = CUBE(N);
	//int dummyseed; // = 100;
	//double dummyoccupancy; // = 0.5;

	struct percAnalysis results;
	vector<unsigned int>* cnt;

public:	
	//class timer horatiocaine;

	void initialize( vector<bool>* threshold, unsigned int cube_edge_length );
	//void initialize( void );
	//void printStatus( bool io, const char* fname );
	void hk3d_core_nonperiodic( unsigned int x, unsigned int y, unsigned int z );
	//void hk3d_core_periodic( unsigned int x, unsigned int y, unsigned int z );
	void hoshen_kopelman( void );
	void compactify( void );
	void checkLabeling( void );
	//void characterize( bool io, const char* fname );
	bool percolates( void );
	void determine_clustersize_distr( void );
	void handover_distribution( unsigned int* dist );
	void reset( void );

	percAnalyzer(){
		id = NULL;
		finder = NULL;
		cnt = NULL;
		setGridsize( 0 );
		//setDummyOccupancy( 0.0 );
		//setDummySeed( 0 );
	}
	~percAnalyzer(){
		delete [] id;
		id = NULL;
		delete finder;
		finder = NULL;
		delete cnt;
		cnt = NULL;
		finder = NULL;
		setGridsize( 0 );
	}

	void setGridsize( unsigned int n ) { N = n; NN = SQR(N); NNN = CUBE(N); }
	//void setDummyOccupancy( double p ) { dummyoccupancy = p; }
	//void setDummySeed( int sd ){ dummyseed = sd; }
	unsigned int getNCluster( void ) {
		return results.nClusterTotal;
	}
	unsigned int getClusterSize( void ) {
		return results.LargestClusterCnt;
	}

};


#endif
