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

void caHdl::ompshar_init_regions( int threadid )
{
	//MK::IS CALLED FROM WITHIN PARALLEL REGION!

	//thread-local creation of caregionMemHdl in parallel and local first touch
	caregionMemHdlP tlr = new caregionMemHdl;

	//constructor implicitly called and making first touch to assure threadlocal memory utilization
	tlr->mycaHdl = this;
	this->regions[threadid] = tlr;

	this->regions[threadid]->mythreadid = threadid;

	//#pragma omp critical { cout << "t=" << omp_get_thread_num() << ",tid=" << this->regions[threadid]->mythreadid << endl; } //impossible

	//##MK::init further caregion variables
}


void caHdl::ompshar_ca2regions( void )
{
	//MK:: IS CALLED FROM WITHIN PARALLEL REGION
	//MK::static partitioning of structure in omp_get_num_threads() = npz * npy non-overlapping regions that comprise the whole domain
	//sectioned along y = td, z = nd for cache locality, npz >= npy \in \mathbb{N} to make implicit xy 2D array along y as long as possible

	//categorize first whether OMP_NUM_THREADS is prime
	uint32_t nregions = omp_get_num_threads();

	uint32_t npy = localmath.smallestYFactor( nregions );
	uint32_t npz = nregions / npy;
	QUICKASSERT ( (npy * npz) == nregions );
	QUICKASSERT ( npy <= npz );
	if ( this->myensHdl->myRank == MASTER && omp_get_thread_num() == MASTER ) { cout << "The nodegrid is X/Y/Z = " << 1 << "__" << npy << "__" << npz << endl; }

	//partition geometry
	uint32_t nbx = this->myCAGeometry.nboxedge_rd; //[0, nbi) only are included in global domain!
	uint32_t nby = this->myCAGeometry.nboxedge_td;
	uint32_t nbz = this->myCAGeometry.nboxedge_nd;
	//cells per region
	uint32_t nry = nby / npy;
	uint32_t nrz = nbz / npz;
	//each domain is required an inner and outer region, i.e. at least 3 cells thick but that is inefficient, as then Area/Volume --> 1
	QUICKASSERT ( nry > MINIMUM_REGION_SIZE && nrz > MINIMUM_REGION_SIZE );

	//determine geometrical mapping of SU domain to memory regions and region limits
	uint32_t* limy = new uint32_t[2*npy]; //min,max
	uint32_t* limz = new uint32_t[2*npz];

	for ( uint32_t y = 0; y < npy; y++ ) { //static partitioning
		limy[(2*y)+0] = (y*nry) + 0;		limy[(2*y)+1] = (y*nry) + nry - 1;
	}
	if ( limy[(2*(npy-1))+1] < (nby-1) )	{ limy[(2*(npy-1))+1] = nby - 1; } //last region too small? enlarge along y!
	if ( limy[(2*(npy-1))+1] >= nby )		{ limy[(2*(npy-1))+1] = nby - 1; } //last maximum protruding beyond domain? chop to size!

	for ( uint32_t z = 0; z < npz; z++ ) { //vice versa along x
		limz[(2*z)+0] = (z*nrz) + 0;		limz[(2*z)+1] = (z*nrz) + nrz - 1;
	}
	if ( limz[(2*(npz-1))+1] < (nbz-1) )	{ limz[(2*(npz-1))+1] = nbz - 1; }
	if ( limz[(2*(npz-1))+1] >= nbz )		{ limz[(2*(npz-1))+1] = nbz - 1; }

	//this orphaned function is executed in parallel, i.e. each thread may pick only once! and writes to disjoint positions in shared variable regions
	//exemplifying a 12 = 3*4 partitioning with npy = 3 and npz = 4
	//z=3::9,10,11
	//z=2::6,7,8
	//z=1::3,4,5
	//z=0::0,1,2

	uint32_t ti = 0;
	for ( uint32_t z = 0; z < npz; z++ ) { 
		//MK::here we define the internal organization of the regions buffer, see example for npy = 3 and npz = 4 as an example then the z coordinates will read as
		//id = 0 1 2 3 4 5 6 7 8 9 10 11
		//z = 0 0 0 1 1 1 2 2 2 3 3 3
		//y = 0 1 2 0 1 2 0 1 2 0 1 2
		//x = 0 0 0 0 0 0 0 0 0 0 0 0
		for ( uint32_t y = 0; y < npy; y++ ) { //domain grid hence, stacking in xyslabs in +z

			if ( omp_get_thread_num() == ti ) { //define local domain and allocate memory thread locally
				regions[ti]->myGeom.nreg_rdmin = 0;
				regions[ti]->myGeom.nreg_rdmax = nbx - 1; //-1 because min,max are for the regions ALWAYS inclusive!
				regions[ti]->myGeom.nreg_tdmin = limy[(2*y)+0];
				regions[ti]->myGeom.nreg_tdmax = limy[(2*y)+1];
				regions[ti]->myGeom.nreg_ndmin = limz[(2*z)+0];
				regions[ti]->myGeom.nreg_ndmax = limz[(2*z)+1];

				regions[ti]->myGeom.nreg_rd = regions[ti]->myGeom.nreg_rdmax - regions[ti]->myGeom.nreg_rdmin + 1;
				regions[ti]->myGeom.nreg_td = regions[ti]->myGeom.nreg_tdmax - regions[ti]->myGeom.nreg_tdmin + 1;
				regions[ti]->myGeom.nreg_nd = regions[ti]->myGeom.nreg_ndmax - regions[ti]->myGeom.nreg_ndmin + 1;
				regions[ti]->myGeom.nregarea_rdtd = regions[ti]->myGeom.nreg_rd * regions[ti]->myGeom.nreg_td;
				regions[ti]->myGeom.nregvol_rdtdnd = regions[ti]->myGeom.nreg_rd * regions[ti]->myGeom.nreg_td * regions[ti]->myGeom.nreg_nd;

				regions[ti]->myGeom.nedge_global_rd = nbx;
				regions[ti]->myGeom.nedge_global_td = nby;
				regions[ti]->myGeom.nedge_global_nd = nbz;

				regions[ti]->myGeom.cellsize = myCAGeometry.cellsize;

				regions[ti]->thePartitioning.nreg_rdx = 1;
				regions[ti]->thePartitioning.nreg_tdy = npy;
				regions[ti]->thePartitioning.nreg_ndz = npz;

				//now the thread may request for local memory to represent the domain
				regions[ti]->ompshar_init_mycellgrid();

				//MK::an openmp barrier is not necessary as all know the same partitioning
				regions[ti]->ompshar_init_neighbortopology( npy, npz );

				//MK::openmp barrier not necessary because regions[ti] knows already its dimensions
				//so it can setup its own region
				regions[ti]->ompshar_init_haloregions( npy, npz );
			}
			ti++;
		}
	}

	delete [] limy; limy = NULL;
	delete [] limz; limz = NULL;
}


void caHdl::ompshar_init_regionlimits( int threadid )
{
	//CALLED FROM WITHIN PARALLEL REGION
	//read into a threadlocal copy the limits of all myneighbors
	uint32_t nbortid;
	for (uint32_t nb = 0; nb < (NUMBER_OF_NEIGHBORS + MYSELF ); nb++ ) {
		nbortid = regions[threadid]->mynborregions[(nb*IDS_AND_LIMITS)+THE_IDS];

		regions[threadid]->mynborregions[(nb*IDS_AND_LIMITS)+THE_XMIN] = regions[nbortid]->myGeom.nreg_rdmin;
		regions[threadid]->mynborregions[(nb*IDS_AND_LIMITS)+THE_XMAX] = regions[nbortid]->myGeom.nreg_rdmax;
		regions[threadid]->mynborregions[(nb*IDS_AND_LIMITS)+THE_YMIN] = regions[nbortid]->myGeom.nreg_tdmin;
		regions[threadid]->mynborregions[(nb*IDS_AND_LIMITS)+THE_YMAX] = regions[nbortid]->myGeom.nreg_tdmax;
		regions[threadid]->mynborregions[(nb*IDS_AND_LIMITS)+THE_ZMIN] = regions[nbortid]->myGeom.nreg_ndmin;
		regions[threadid]->mynborregions[(nb*IDS_AND_LIMITS)+THE_ZMAX] = regions[nbortid]->myGeom.nreg_ndmax;
	}
}


uint32_t caHdl::exists_nborhalo( int regid, short dx, short dy, short dz, int npx, int npy, int npz, uint32_t cand )
{
	//CALLED FROM WITHIN PARALLEL REGION
	//finds the list position that a particular halo has in the list of a nboring memory region
	//the function assumes it sits inside the region regid and moves in dx,dy,dz direction to find its if in this neighbor 
	//a halo has been created,i.e thisneighbor->myhalos[cand]->theHaloCells != NULL && extendxyz > 0 function returns coordinates of regid
	int targetreg = regid;

	int ix = 0;						//always, because no partitioning along x direction to maintain long cache lines
	int iy = regid % npy;
	int iz = regid / npy;

	//get new coordinate
	ix += dx;
	iy += dy;
	iz += dz;

	//periodic boundary conditions
	if ( ix < 0 )		ix += npx;
	if ( ix >= npx )	ix -= npx;
	if ( iy < 0 )		iy += npy;
	if ( iy >= npy )	iy -= npy;
	if ( iz < 0 )		iz += npz;
	if ( iz >= npz )	iz -= npz;

	targetreg = ( (iz * npy) + iy ); 

	if ( this->regions[targetreg]->myhalos[cand].extendxyz > 0 && this->regions[targetreg]->myhalos[cand].theHaloCells != NULL ) {
		return cand;
	}

	return HALO_NOT_EXISTENT;
}


void caHdl::ompshar_link_haloregions( int threadid )
{
	//CALLED FROM WITHIN PARALLEL REGION
	//connects regions[threadid] physically to the neighboring regions by assigning which entry in the neighboring regions haloregions
	//is the one that is neighboring regions[threadid]. This ??might be required to read out quickly the haloregion
	//of a neighbor to synchronize with the regions[thread] ids own mycellgrid 
	//####
	regions[threadid]->mynborhalos = NULL;
	regions[threadid]->mynborhalos = new uint32_t[(MYSELF+NUMBER_OF_NEIGHBORS)]; //IDs of halos in list of neighboring regions, do not reference locally!
	QUICKASSERT ( regions[threadid]->mynborhalos != NULL );

	uint32_t party = regions[threadid]->thePartitioning.nreg_tdy;
	uint32_t partz = regions[threadid]->thePartitioning.nreg_ndz;

	//self
	regions[threadid]->mynborhalos[0] = HALO_NOT_EXISTENT; //0,0,0

	//six faces, exists_nborhalo checks whether there has been a halo created at the adjacent site, opposite cand
	regions[threadid]->mynborhalos[1] = exists_nborhalo( threadid, -1, 0, 0, 1, party, partz, 2 );
	regions[threadid]->mynborhalos[2] = exists_nborhalo( threadid, +1, 0, 0, 1, party, partz, 1 );
	regions[threadid]->mynborhalos[3] = exists_nborhalo( threadid, 0, -1, 0, 1, party, partz, 4 );
	regions[threadid]->mynborhalos[4] = exists_nborhalo( threadid, 0, +1, 0, 1, party, partz, 3 );
	regions[threadid]->mynborhalos[5] = exists_nborhalo( threadid, 0, 0, -1, 1, party, partz, 6 );
	regions[threadid]->mynborhalos[6] = exists_nborhalo( threadid, 0, 0, +1, 1, party, partz, 5 );

	//twelve edges
	regions[threadid]->mynborhalos[7] = exists_nborhalo( threadid, -1, 0, -1, 1, party, partz, 10 );
	regions[threadid]->mynborhalos[8] = exists_nborhalo( threadid, +1, 0, -1, 1, party, partz, 9 );
	regions[threadid]->mynborhalos[9] = exists_nborhalo( threadid, -1, 0, +1, 1, party, partz, 8 );
	regions[threadid]->mynborhalos[10] = exists_nborhalo( threadid, +1, 0, +1, 1, party, partz, 7 );

	regions[threadid]->mynborhalos[11] = exists_nborhalo( threadid, -1, -1, 0, 1, party, partz, 14 );
	regions[threadid]->mynborhalos[12] = exists_nborhalo( threadid, -1, +1, 0, 1, party, partz, 13 );
	regions[threadid]->mynborhalos[13] = exists_nborhalo( threadid, +1, -1, 0, 1, party, partz, 12 );
	regions[threadid]->mynborhalos[14] = exists_nborhalo( threadid, +1, +1, 0, 1, party, partz, 11 );

	regions[threadid]->mynborhalos[15] = exists_nborhalo( threadid, 0, -1, -1, 1, party, partz, 18 );
	regions[threadid]->mynborhalos[16] = exists_nborhalo( threadid, 0, -1, +1, 1, party, partz, 17 );
	regions[threadid]->mynborhalos[17] = exists_nborhalo( threadid, 0, +1, -1, 1, party, partz, 16 );
	regions[threadid]->mynborhalos[18] = exists_nborhalo( threadid, 0, +1, +1, 1, party, partz, 15 );

	//eight corners, eight corners you have ...
	regions[threadid]->mynborhalos[19] = exists_nborhalo( threadid, -1, -1, -1, 1, party, partz, 26 );
	regions[threadid]->mynborhalos[20] = exists_nborhalo( threadid, +1, -1, -1, 1, party, partz, 25 );
	regions[threadid]->mynborhalos[21] = exists_nborhalo( threadid, -1, +1, -1, 1, party, partz, 24 );
	regions[threadid]->mynborhalos[22] = exists_nborhalo( threadid, +1, +1, -1, 1, party, partz, 23 );
	regions[threadid]->mynborhalos[23] = exists_nborhalo( threadid, -1, -1, +1, 1, party, partz, 22 );
	regions[threadid]->mynborhalos[24] = exists_nborhalo( threadid, +1, -1, +1, 1, party, partz, 21 );
	regions[threadid]->mynborhalos[25] = exists_nborhalo( threadid, -1, +1, +1, 1, party, partz, 20 );
	regions[threadid]->mynborhalos[26] = exists_nborhalo( threadid, +1, +1, +1, 1, party, partz, 19 );

/*
#pragma omp critical
{
	cout << threadid << "--LINKING-HALOS;0;1;...;26--"; 
	for (uint32_t h = 0; h < (MYSELF+NUMBER_OF_NEIGHBORS); h++) { cout << ";" << regions[threadid]->mynborhalos[h];}
	cout << endl;
}
*/
}

#endif