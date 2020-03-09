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


/*
	Implements functionalities to detect whether the microstructure has a percolating network of RX grains
	utilizes the Hoshen-Kopelman algorithm in 3D and implements a pathcompressed weighted union/find algorithm based
	on Kartik Kukreja's implementation http://kartikkukreja.wordpress.com
	Markus Kühbach, m.kuehbach (at) mpie.de

	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "SCORE_PercAnalyzer.h"
////#include <iostream>
//#include <random>
////#include <algorithm>
////#include <fstream>
////#include <iomanip>
////#include <stdlib.h>
//#include <stdint.h>	//uint64 def
//#include <time.h>	//clock_gettime
////#include <string>

using namespace std;

/*
//###these macros are also defined in SCORE
#define MIN(X,Y)				(((X) < (Y)) ? (X) : (Y))
#define MAX(X,Y) 				(((X) > (Y)) ? (X) : (Y))
#define SQR(X)					((X)*(X))
#define CUBE(X)					((X)*(X)*(X))
*/





/*
void percAnalyzer::initialize( void )
{
cout << "SCORE PercolationAnalyzer 3D, SC site perc N = " << (unsigned int) N << "^3" << " p = " << dummyoccupancy << " seed = " << dummyseed << endl;
cout << "Initializing..." << endl;
	//this->horatiocaine.profiling_start();

	//get memory for the HK analysis and the union/find
	id = new unsigned int[NNN];
	finder = new UF(NNN);		//##MK::consider optimize as this is a very large preallocation
	
cout << "Generating random network structure..." << endl;
	//generate a dummy structure instead of passing a SCORE grid
	std::mt19937 generator ( dummyseed );
	std::uniform_real_distribution<double> dis(0.0, 1.0);
	
	for (unsigned int c = 0; c < NNN; c++ ) { //id field 0 (non-rxed), 1 (rxed)
		id[c] = 0;
		if ( dis(generator) < dummyoccupancy ) { id[c] = 1; }
	}

	//this->horatiocaine.profiling_end();
	//this->horatiocaine.profiling_diff_mem( "Initializing" );
}
*/


void percAnalyzer::initialize( vector<bool>* inputdata, unsigned int cube_edge_length )
{
	setGridsize( cube_edge_length );

	cout << "\t\tSCORE PercolationAnalyzer 3D, SC site perc N = " << (unsigned int) N << "^3" << endl; //<< " p = " << dummyoccupancy << " seed = " << dummyseed << endl;
	cout << "\t\t\t\tInitializing..." << endl;
	//this->horatiocaine.profiling_start();

	//get memory for the HK analysis and the union/find
	id = new unsigned int[NNN];
	finder = new UF(NNN);		//##MK::consider optimize as this is a very large preallocation

	//cout << "Generating random network structure..." << endl;
	cout << "\t\t\t\tLoading binary input structure..." << endl;

	for (unsigned int c = 0; c < NNN; c++ ) { //id field 0 (non-rxed), 1 (rxed)
		id[c] = 0;
		if ( inputdata->at(c) == true ) 
			id[c] = 1;
	}

	//this->horatiocaine.profiling_end();
	//this->horatiocaine.profiling_diff_mem( "Initializing" );
}


/*
void percAnalyzer::printStatus( bool io, const char* fname )
{
cout << "I/O writing " << fname << endl;
	if ( io == true ) {
		FILE* pFile;
		pFile = fopen( fname, "wb" );
		fwrite( this->id, sizeof(unsigned int), NNN, pFile);
		fclose(pFile);
	}
	else {
		//###implement stdout print only
	}
}
*/


void percAnalyzer::hk3d_core_nonperiodic( unsigned int x, unsigned int y, unsigned int z )
{
	if ( id[x+y*N+z*NN] == 1 ) { //only RXed (1) to analyze for percolation, periodic boundaries
		//the following three lines DO NOT implement PBC correctly because they refer to labels not yet assigned
		//unsigned int front = ( z == 0 ? id[x+y*N+(N-1)*NN] : id[x+y*N+(z-1)*NN] ); //catches also the wraparound of the unsigned int
		//unsigned int bottom = ( y == 0 ? id[x+(N-1)*N+z*NN] : id[x+(y-1)*N+z*NN] );
		//unsigned int left = ( x == 0 ? id[(N-1)+y*N+z*NN] : id[(x-1)+y*N+z*NN] );
		
		unsigned int front = ( z == 0 ? 0 : id[x+y*N+(z-1)*NN] );
		unsigned int bottom = ( y == 0 ? 0 : id[x+(y-1)*N+z*NN] );
		unsigned int left = ( x == 0 ? 0 : id[(x-1)+y*N+z*NN] );					
		
		bool bl = false;
		bool bb = false;
		bool bf = false;
		if ( left > 0 ) bl = true; //> 0, i.e. cell exists surplus has already (at least initially been labelled) once
		if ( bottom > 0 ) bb = true;
		if ( front > 0 ) bf = true;
		
		//##MK::consider in further optimization that the likelihood for these cases is a function of the RX fraction to be able to test less ifs
		//##MK::utilize a checksum like bl*100+bb*10+bf*1 in a switch command maybe more efficient...
		//case 000 all neighbors are not-RXed
		if ( bl == false && bb == false && bf == false ) {
			id[x+y*N+z*NN] = finder->initialAssgn( x+y*N+z*NN ); //because UF assumes initially already NNN disjoint labels!
			return; //continue;
		}
		
		if ( bl == true && bb == false && bf == false ) { //100
			id[x+y*N+z*NN] = left;
			return; //continue;
		}
		
		if ( bl == false && bb == true && bf == false ) { //010
			id[x+y*N+z*NN] = bottom;
			return; //continue;
		}

		if ( bl == false && bb == false && bf == true ) { //001
			id[x+y*N+z*NN] = front;
			return; //continue;
		}
		
		if ( bl == true && bb == true && bf == false ) { //110
			id[x+y*N+z*NN] = this->finder->merge( left, bottom );
			return; //continue;
		}
		if ( bl == true && bb == false && bf == true ) { //101
			id[x+y*N+z*NN] = this->finder->merge( left, front );
			return; //continue;
		}
		if ( bl == false && bb == true && bf == true ) { //011
			id[x+y*N+z*NN] = this->finder->merge( bottom, front );
			return; //continue;
		}
		if ( bl == true && bb == true && bf == true ) {
			//unsigned int cand1 = this->finder->merge( left, bottom );
			//unsigned int cand2 = this->finder->merge( left, front );
			unsigned int cand3 = this->finder->merge( bottom, front );
			//id[x+y*N+z*NN] = MIN(cand2, cand3); //##MKmay be tricky
			id[x+y*N+z*NN] = cand3; //because call for cand1 and cand2 modify tree representation
			//return; //continue;
		}
	} //next cell
}


/*void percAnalyzer::hk3d_core_periodic( unsigned int x, unsigned int y, unsigned int z )
{
	//##MK::id[..] > UNKNOWN_ID...
	if ( id[x+y*N+z*NN] == 1 ) { //only RXed (1) to analyze for percolation, periodic boundaries
		//MK::these lines implement PBC only correctly if HK ran already over the entire domain! see also hoshen_kopelman function for further details!
		unsigned int front = ( z == 0 ? id[x+y*N+(N-1)*NN] : id[x+y*N+(z-1)*NN] ); //catches also the wraparound of the unsigned int
		unsigned int bottom = ( y == 0 ? id[x+(N-1)*N+z*NN] : id[x+(y-1)*N+z*NN] );
		unsigned int left = ( x == 0 ? id[(N-1)+y*N+z*NN] : id[(x-1)+y*N+z*NN] );
				
		bool bl = false;
		bool bb = false;
		bool bf = false;
		if ( left > 0 ) bl = true; //> 0, i.e. cell exists surplus has already (at least initially been labelled) once
		if ( bottom > 0 ) bb = true;
		if ( front > 0 ) bf = true;
		
		//##MK::consider in further optimization that the likelihood for these cases is a function of the RX fraction to be able to test less ifs
		//##MK::utilize a checksum like bl*100+bb*10+bf*1 in a switch command maybe more efficient...
		//case 000 all neighbors are not-RXed
		if ( bl == false && bb == false && bf == false ) {
			id[x+y*N+z*NN] = finder->initialAssgn( x+y*N+z*NN ); //because UF assumes initially already NNN disjoint labels!
			return; //continue;
		}
		
		if ( bl == true && bb == false && bf == false ) { //100
			id[x+y*N+z*NN] = left;
			return; //continue;
		}
		
		if ( bl == false && bb == true && bf == false ) { //010
			id[x+y*N+z*NN] = bottom;
			return; //continue;
		}

		if ( bl == false && bb == false && bf == true ) { //001
			id[x+y*N+z*NN] = front;
			return; //continue;
		}
		
		if ( bl == true && bb == true && bf == false ) { //110
			id[x+y*N+z*NN] = this->finder->merge( left, bottom );
			return; //continue;
		}
		if ( bl == true && bb == false && bf == true ) { //101
			id[x+y*N+z*NN] = this->finder->merge( left, front );
			return; //continue;
		}
		if ( bl == false && bb == true && bf == true ) { //011
			id[x+y*N+z*NN] = this->finder->merge( bottom, front );
			return; //continue;
		}
		if ( bl == true && bb == true && bf == true ) {
			unsigned int cand1 = this->finder->merge( left, bottom );
			unsigned int cand2 = this->finder->merge( left, front );
			unsigned int cand3 = this->finder->merge( bottom, front );
			//id[x+y*N+z*NN] = MIN(cand2, cand3); //##MKmay be tricky
			id[x+y*N+z*NN] = cand3; //because call for cand1 and cand2 modify tree representation
			//return; //continue;
		}
	} //next cell
}*/


void percAnalyzer::hoshen_kopelman( void )
{
cout << "\t\t\t\tLabeling cluster via Hoshen-Kopelman..." << endl;
	//this->horatiocaine.profiling_start();

	for( unsigned int zz = 0; zz < N; zz++ ) {
//if ( zz % STDOUT_Z_PROMPT == 0 ) cout << "\t\t" << zz << endl;
		for( unsigned int yy = 0; yy < N; yy++ ) {
			for( unsigned int xx = 0; xx < N; xx++ ) {
				hk3d_core_nonperiodic( xx, yy, zz );
			}
		}
	}
/*
	//this labeling is not aware of the desired periodic boundary conditions. However according to the suggestion by
	//N. Weik et. al. "Graphene with vacancies: supernumerary zero modes" arXiv:1603.00212v3 we can correct for this by 
	//sweeping subsequently again over the left (-x), front (-y), and the ottom (-z) face 
	//starting with the left face
	unsigned int xx = 0;
	for ( unsigned int zz = 0; zz < N; zz++ ) {
		for (unsigned int yy = 0; yy < N; yy++ ) {
			hk3d_core_periodic( xx, yy, zz );
		}
	}
	//now the front face
	unsigned int yy = 0;
	for ( unsigned int zz = 0; zz < N; zz++ ) {
		for (unsigned int xx = 0; xx < N; xx++ ) {
			hk3d_core_periodic( xx, yy, zz );
		}
	}
	//lastly the bottom face
	unsigned int zz = 0;
	for ( unsigned int yy = 0; yy < N; yy++ ) {
		for (unsigned int xx = 0; xx < N; xx++ ) {
			hk3d_core_periodic( xx, yy, zz );
		}
	}	
*/	
	//this->horatiocaine.profiling_end();
	//this->horatiocaine.profiling_diff_mem( "HKLabeling" );
}


void percAnalyzer::compactify( void )
{
cout << "\t\t\t\tCompactifying..." << endl;
	//this->horatiocaine.profiling_start();

	vector<compactLabel>* eqsc = NULL;
	eqsc = new vector<compactLabel>;
	QUICKASSERT( eqsc != NULL );
	unsigned int nnew = 0;
	unsigned int croot = 0;
	for ( unsigned int c = 0; c < NNN; c++ ) {
		if ( id[c] > 0 ) {
			croot = this->finder->find_root( id[c] );
			unsigned int e = 0;
			for ( e = 0; e < eqsc->size(); e++ ) {
				if ( eqsc->at(e).oky == croot ) {
					break;
				}
			}
			if ( e == eqsc->size() ) { //new key
				nnew++;
				struct compactLabel al;
				al.oky = croot;
				al.nky = nnew;
				eqsc->push_back( al );
				id[c] = al.nky;
				continue;
			}
			//old key
			id[c] = eqsc->at(e).nky;
		}
	}
	delete eqsc;

	/*unsigned int* eqsc = NULL;
	eqsc = new unsigned int[NNN]; //##MK::maybe unnecessary large buffer
	QUICKASSERT( eqsc != NULL );
	unsigned int nnew = 0;
	for ( unsigned int c = 0; c < NNN; c++ ) { eqsc[c] = 0; }
	for ( unsigned int z = 0; z < N; z++ ) {
		for ( unsigned int y = 0; y < N; y++ ) {
			for ( unsigned int x = 0; x < N; x++ ) {
				if ( id[x+y*N+z*NN] > 0 ) {
					unsigned int cc = this->finder->find_root( id[x+y*N+z*NN] ); //MK::finds the root?
					if ( eqsc[cc] == 0 ) { //if a compact label has not been assigned, do so now
						nnew++;
						eqsc[cc] = nnew;
					}
					id[x+y*N+z*NN] = eqsc[cc];
				}
			}
		}
	}
	delete [] eqsc;*/


	results.nClusterTotal = nnew;

	//this->horatiocaine.profiling_end();
	//this->horatiocaine.profiling_diff_mem( "Compactifying" );
cout << "\t\t\t\tTotal number of clusters = " << results.nClusterTotal << endl;
}

void percAnalyzer::checkLabeling( void )
{
cout << "\t\t\t\tChecking the correct labeling with the compactified IDs..." << endl;
	bool labelingerror = false;
	unsigned int LE,RI,FR,RE,BO,TO;
	for ( unsigned int z = 0; z < N; z++ ) {
//if ( z % STDOUT_Z_PROMPT == 0 ) cout << "\t\t" << z << endl;
		for (unsigned int y = 0; y < N; y++ ) {
			for ( unsigned int x = 0; x < N; x++ ) {
				if ( id[x+y*N+z*NN] > 0 ) {
					//original Tobin Fricke i-->y, j --> x
					/*
					LE = ( x == 0 	? id[(N-1)+y*N+z*NN]	: id[(x-1)+y*N+z*NN] );
					RI = ( x == N-1	? id[0+y*N+z*NN] 		: id[(x+1)+y*N+z*NN] );
					FR = ( y == 0	? id[x+(N-1)*N+z*NN]	: id[x+(y-1)*N+z*NN] );
					RE = ( y == N-1	? id[x+0*N+z*NN]		: id[x+(y+1)*N+z*NN] );
					BO = ( z == 0	? id[x+y*N+(N-1)*NN]	: id[x+y*N+(z-1)*NN] );
					TO = ( z == N-1	? id[x+y*N+0*NN]		: id[x+y*N+(z+1)*NN] );
					*/
					//non-periodic domain
					LE = ( x == 0 	? 0	: id[(x-1)+y*N+z*NN] );
					RI = ( x == N-1	? 0	: id[(x+1)+y*N+z*NN] );
					FR = ( y == 0	? 0	: id[x+(y-1)*N+z*NN] );
					RE = ( y == N-1	? 0	: id[x+(y+1)*N+z*NN] );
					BO = ( z == 0	? 0	: id[x+y*N+(z-1)*NN] );
					TO = ( z == N-1	? 0	: id[x+y*N+(z+1)*NN] );
							
					unsigned int cand = id[x+y*N+z*NN]; //von Neumann nearest neighbors must have the same label if they are not 0!
					if ( LE != 0 && LE != cand ) 	labelingerror = true; //cout << "LE error " << x << ";" << y << ";" << z << endl;
					if ( RI != 0 && RI != cand ) 	labelingerror = true; //cout << "RI error " << x << ";" << y << ";" << z << endl;
					if ( FR != 0 && FR != cand ) 	labelingerror = true; //cout << "FR error " << x << ";" << y << ";" << z << endl;
					if ( RE != 0 && RE != cand ) 	labelingerror = true; //cout << "RE error " << x << ";" << y << ";" << z << endl;
					if ( BO != 0 && BO != cand ) 	labelingerror = true; //cout << "BO error " << x << ";" << y << ";" << z << endl;
					if ( TO != 0 && TO != cand ) 	labelingerror = true; //cout << "TO error " << x << ";" << y << ";" << z << endl;

					QUICKASSERT( labelingerror == false );
				} //label consistency for cluster cell x,y,z checked
			}
		}
	}
}


/*void percAnalyzer::characterize( bool io, const char* fname )
{
cout << "Characterizing..." << endl;
	//this->horatiocaine.profiling_start();

	//count cluster size
	unsigned int* cnt = NULL;
	cnt = new unsigned int[1+results.nClusterTotal];
	for ( unsigned int i = 0; i < (1+results.nClusterTotal); i++ ) { cnt[i] = 0; }

	unsigned int cand = 0;
	for ( unsigned int c = 0; c < NNN; c++ ) {
			cand = id[c];
			if ( cand > 0 ) { cnt[cand] += 1; } //prevents counting of non-RXed volume into a large cluster
	}
	
	//MK::DEBUG analysis
	results.LargestClusterCnt = 0;
	results.LargestClusterID = 0;	
	results.PL = 0.0;
	
	for ( unsigned int cluster = 1; cluster <= results.nClusterTotal; cluster++ ) { 
		if ( cnt[cluster] > results.LargestClusterCnt ) { 
			results.LargestClusterCnt = cnt[cluster];
			results.LargestClusterID = cluster;
		}
	}
	
	if ( results.LargestClusterCnt > 0 && NNN > 0 ) {
		results.PL = (double) results.LargestClusterCnt / (double) NNN;
	}
	
	delete [] cnt;

	//this->horatiocaine.profiling_end();
	//this->horatiocaine.profiling_diff_mem( "Characterizing" );

cout << "\t\tLargest cluster is cluster " << results.LargestClusterID << " with " << results.LargestClusterCnt << " cells." << endl;
cout << "\t\tThis cluster occupies a volume fraction of P_L = " << setprecision(8) << results.PL << " of the entire domain." << endl;

	if ( io == true && results.LargestClusterID > UNKNOWN_ID ) {
		unsigned int* tmp = NULL;
		tmp = new unsigned int[NNN];
		
		for ( unsigned int c = 0; c < NNN; c++ ) { //threshold grid to extract the cells belonging only to the largest cluster
			tmp[c] = (id[c] == results.LargestClusterID ? results.LargestClusterID : UNKNOWN_ID);
			//tmp[c] = UNKNOWN_ID;
			//if ( id[c] == results.LargestClusterID ) { tmp[c] = results.LargestClusterID; }
		}
		
		FILE* pFile;
		pFile = fopen( fname, "wb" );
		fwrite( tmp, sizeof(unsigned int), NNN, pFile);
		fclose(pFile);
	
		delete [] tmp;
	}
}*/


bool percAnalyzer::percolates( void )
{
cout << "\t\t\t\tAnalyzing for possible percolating cluster..." << endl;
	//this->horatiocaine.profiling_start();

	bool networkstatus = false;

	//test three pairs of opposite sides (left/right, front/rear, bottom/top) for percolation
	bool* b1 = NULL;
	bool* b2 = NULL;
	b1 = new bool[1+results.nClusterTotal];
	b2 = new bool[1+results.nClusterTotal];
	for (unsigned int cluster = 0; cluster < (1+results.nClusterTotal); cluster++ ) { b1[cluster] = false; b2[cluster] = false; }
	
	//for an x,y,z ordered implicit array scanning of the YZ planes is cache inefficient these are the left -x and the right +x sides so test them last!
//check for percolation along the x direction by scanning left -x and right +x face
	unsigned int cand = 0;	
	unsigned int x = 0;
	for (unsigned int z = 0; z < N; z++ ) {
		for (unsigned int y = 0; y < N; y++ ) {
			cand = id[x+y*N+z*NN];
			if ( cand > 0 ) { b1[cand] = true; }
		}
	}
	cand = 0;
	x = N-1;
	for (unsigned int z = 0; z < N; z++ ) {
		for (unsigned int y = 0; y < N; y++ ) {
			cand = id[x+y*N+z*NN];
			if ( cand > 0 ) { b2[cand] = true; }
		}
	}
	
	for ( unsigned int cluster = 0; cluster < (1+results.nClusterTotal); cluster++ ) {
		if ( b1[cluster] == true && b2[cluster] == true ) {
			networkstatus = true;
cout << "\t\t\t\tPercolation along X for cluster " << cluster << endl;
		}
	}
	for (unsigned int cluster = 0; cluster < (1+results.nClusterTotal); cluster++ ) { b1[cluster] = false; b2[cluster] = false; }
	
//check for percolation along the y direction by scanning front -y and rear +y face
	cand = 0;
	unsigned int y = 0;
	unsigned int offset = 0;
	for (unsigned int z = 0; z < N; z++ ) {
		offset = y*N+z*NN;
		for (unsigned int x = 0; x < N; x++ ) {
			cand = id[x+offset];
			if ( cand > 0 ) { b1[cand] = true; }
		}
	}
	cand = 0;
	y = N-1;
	offset = 0;
	for (unsigned int z = 0; z < N; z++ ) {
		offset = y*N+z*NN;
		for (unsigned int x = 0; x < N; x++ ) {
			cand = id[x+offset];
			if ( cand > 0 ) { b2[cand] = true; }
		}
	}

	for ( unsigned int cluster = 0; cluster < (1+results.nClusterTotal); cluster++ ) {
		if ( b1[cluster] == true && b2[cluster] == true ) {
			networkstatus = true;
cout << "\t\t\t\tPercolation along Y for cluster " << cluster << endl;
		}
	}
	for (unsigned int cluster = 0; cluster < (1+results.nClusterTotal); cluster++ ) { b1[cluster] = false; b2[cluster] = false; }
	
//check for percolation along the z direction by scanning bottom -z and top +z face
	cand = 0;
	unsigned int z = 0;
	offset = 0;
	for (unsigned int y = 0; y < N; y++) {
		offset = y*N+z*NN;
		for (unsigned int x = 0; x < N; x++) {
			cand = id[x+offset];
			if ( cand > 0 ) { b1[cand] = true; }
		}
	}
	cand = 0;
	z = N-1;
	offset = 0;
	for (unsigned int y = 0; y < N; y++) {
		offset = y*N+z*NN;
		for (unsigned int x = 0; x < N; x++) {
			cand = id[x+offset];
			if ( cand > 0 ) { b2[cand] = true; }
		}
	}

	for ( unsigned int cluster = 0; cluster < (1+results.nClusterTotal); cluster++ ) {
		if ( b1[cluster] == true && b2[cluster] == true ) {
			networkstatus = true;
cout << "\t\t\t\tPercolation along Z for cluster " << cluster << endl;
		}
	}
	
	delete [] b1;
	delete [] b2;
	
	if ( networkstatus == true ) 
		cout << "\t\t\t\tThe structure percolates!" << endl;
	else 
		cout << "\t\t\t\tThe structure DOES NOT percolate!" << endl;

	//this->horatiocaine.profiling_end();
	//this->horatiocaine.profiling_diff_mem( "PercAnalysis" );

	return networkstatus;
}




void percAnalyzer::determine_clustersize_distr( void )
{
	//determine cluster size distribution first
	cout << "\t\t\t\tCharacterizing cluster size distribution..." << endl;

	cnt = new vector<unsigned int>;
	QUICKASSERT( cnt != NULL );
	cnt->assign( results.nClusterTotal, 0 );

	//indices + 1 of cnt are the cluster names, the content cnt[i] their size
	unsigned int cand = 0;
	for ( unsigned int c = 0; c < NNN; c++ ) {
		cand = id[c];
		if ( cand > 0 )
			cnt->at(cand-1) += 1;
	}

	//sort cluster by their size ascendingly and utilize this to find the largest cluster
	std::sort( cnt->begin(), cnt->end() );

	results.LargestClusterCnt = 0;
	results.LargestClusterID = 0;
	results.PL = 0.0;

	if ( results.nClusterTotal > 0 ) {
		results.LargestClusterCnt = cnt->at(cnt->size()-1);
		results.LargestClusterID = cnt->at(cnt->size()-1) + 1;
	}

	if ( NNN > 0 ) 
		results.PL = (double) results.LargestClusterCnt / (double) NNN;

	//add output routine

	//this->horatiocaine.profiling_end();
	//this->horatiocaine.profiling_diff_mem( "Characterizing" );

cout << "\t\t\t\tLargest cluster is cluster " << results.LargestClusterID << " with " << results.LargestClusterCnt << " cells." << endl;
cout << "\t\t\t\tThis cluster occupies a volume fraction of P_L = " << setprecision(8) << results.PL << " of the entire domain." << endl;
}


void percAnalyzer::handover_distribution( unsigned int* dist )
{
	//copy cluster size distribution into preallocated array dist
	for ( unsigned int c = 0; c < cnt->size(); c++ ) {
		dist[c] = cnt->at(c);
	}
}


void percAnalyzer::reset( void ) 
{
	delete [] id;
	id = NULL;
	delete finder;
	finder = NULL;
	setGridsize( 0 );
}


//default main to execute the HK algorithm as a standalone program
/*
#define LENGTH		(1)
#define OCC			(2)
#define SEED		(3)

int main(int argc, char** argv)
{
	//threshold fully RXed voxels from SCORE simulation grid at time t
	percAnalyzer* hk = new percAnalyzer;
	hk->setGridsize( strtoul(argv[LENGTH], NULL, 0 ) );
	hk->setDummyOccupancy( strtold(argv[OCC], NULL ) );
	hk->setDummySeed( atoi(argv[SEED]) );

	hk->initialize();
	//hk->printStatus( true, "SCOREPercAnalyzer3D_0_Init.raw" );

//2D version in HKCXX3D code
//####change UF structure to unsigned int! and in the output function also to unsigned int!
//####check why structure always percolates, furthermore, check coordinate change in DAMASK and SCORE!

	hk->hoshen_kopelman();
//	hk->printStatus( true, "SCOREPercAnalyzer3D_1_BeforeCompact.raw" );

	hk->compactify();
//	hk->printStatus( true, "SCOREPercAnalyzer3D_2_AfterCompact.raw" );

	hk->checkLabeling();

	//checkLabels( lab );
//	hk->characterize( true, "SCOREPercAnalyzer3D_3_LargestCluster.raw");

	bool result = hk->percolates();
	
//hk->horatiocaine.profiling_log();
	delete hk;
	
	//utilize result further... here however end
	return 0;
}
*/


/*
inline bool SortMiMxAsc(const mimx &m1, const mimx &m2)
{
	if (m1.mi < m2.mi)
		return true;
	else if (m1.mi > m2.mi)
		return false;

	//not returned yet: okay m1.mi == m2.mi
	if (m1.mx < m2.mx)
		return true;
	else if (m1.mx > m2.mx)
		return false;

	//##MK::mi and mx for both the same, leave as is
	return false;
}
*/