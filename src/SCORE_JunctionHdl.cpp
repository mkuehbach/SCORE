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

juncHdl::juncHdl()
{
	mycaHdl = NULL;
	mymemreg = NULL;
	mythreadid = MASTER;

	fjuncSkeleton = NULL;
	tjuncSkeleton = NULL;
	hjuncSkeleton = NULL;

	fjuncSkeletonSize = 0;
	tjuncSkeletonSize = 0;
	hjuncSkeletonSize = 0;

	tjuncnextfreeslot = 0;
	hjuncnextfreeslot = 0;

	juncMemGuard = 0.0;
}


juncHdl::~juncHdl()
{
	//MK::YOU MUST NOT DELETE mycaHdl and mymemreg - THEY ARE BACKREFERENCES ONLY!

	//data structure which makes accessible individual grain boundary faces
	if ( this->fjuncSkeleton != NULL ) {
		for ( uint32_t i = 0; i < this->fjuncSkeletonSize; i++ ) {
			if ( this->fjuncSkeleton[i].maxLen > 0 ) {
				for ( uint32_t f = 0; f < this->fjuncSkeleton[i].maxLen; f++ ) {
					delete [] this->fjuncSkeleton[i].thefaces[f].voxelizedface;
				}
				//nothing to delete for on the interval [len+1,maxLen) because these point to zero
				delete [] this->fjuncSkeleton[i].thefaces;
			}
		}
	}
	delete [] this->fjuncSkeleton; this->fjuncSkeleton = NULL; this->fjuncSkeletonSize = 0;

	//two data structures which do not distinguish the locality of the junction-adjacent voxel
	//delete checks for deallocating NULL pointer either way...
	delete [] this->tjuncSkeleton; this->tjuncSkeleton = NULL; this->tjuncSkeletonSize = 0;
	delete [] this->hjuncSkeleton; this->hjuncSkeleton = NULL; this->hjuncSkeletonSize = 0;
}


//MK::functions to read state of mycellgrid in remote threads, individual functions for each direction
//of the von Neumann stencil to detect boundaries to reduce formal complexity of periodic boundary checks
//MK::at the moment however no tracking of periodic domain boundaries, therefore _left x is never 0
//determines ID field of neighbors in an implicit, topology of regions identified from 0, 1, 2, omp_get_num_threads() - 1 that are
//arranged on a party * partz region grid with party <= partz
//addressing is first in z at fixed y then in y+1 along z, y+2 along z and so forth...
//so party = 2 and partz = 4 generates: 0, 1,2,3 with y = 0 and 4,5,6,7 with y = 1
//regions are periodic so in that example

uint32_t juncHdl::read_cellstate_mx( int caller, uint32_t x, uint32_t y, uint32_t z )
{
	//CALLED FROM WITHIN PARALLEL REGION
	//reads mycellgrid state at global location x, y, z which is guaranteed to be in the left (-x) memory region neighboring the caller's or outside the domain
	QUICKASSERT( (x-1) >= 0 );
	//OMP region partitioning does not partition along the x axis
	return INVALID_SEED;
}


uint32_t juncHdl::read_cellstate_px( int caller, uint32_t x, uint32_t y, uint32_t z )
{
	//CALLED FROM WITHIN PARALLEL REGION
	//reads mycellgrid state at global location x, y, z which is guaranteed to be the right (+x) memory region neighbor of caller
	QUICKASSERT( (x+1) < this->mycaHdl->myCAGeometry.nboxedge_rd );
	return INVALID_SEED;
}


uint32_t juncHdl::read_cellstate_py( int caller, uint32_t x, uint32_t y, uint32_t z )
{
	//CALLED FROM WITHIN PARALLEL REGION
	//reads mycellgrid state at global location x, y, z which is guaranteed to be the top (+y) memory region neighbor of caller
	QUICKASSERT( (y+1) < this->mycaHdl->myCAGeometry.nboxedge_td );
	uint32_t npy = this->mycaHdl->regions[caller]->thePartitioning.nreg_tdy;
	//int ix = 0;//get coordinates of regid
	int iy = caller % npy;
	int iz = caller / npy;

	iy += +1;		//rear neighbor
	if ( iy >= this->mycaHdl->regions[caller]->thePartitioning.nreg_tdy )	//simplified periodic boundary conditions
		return INVALID_SEED;

	int rear = ((iz * npy)+iy); //get logical neighbor and his geometry
	uint32_t xmi = this->mycaHdl->regions[rear]->myGeom.nreg_rdmin;
	uint32_t ymi = this->mycaHdl->regions[rear]->myGeom.nreg_tdmin;
	uint32_t zmi = this->mycaHdl->regions[rear]->myGeom.nreg_ndmin;

#ifdef DEBUG
	uint32_t xmx = this->mycaHdl->regions[rear]->myGeom.nreg_rdmax;
	uint32_t ymx = this->mycaHdl->regions[rear]->myGeom.nreg_tdmax;
	uint32_t zmx = this->mycaHdl->regions[rear]->myGeom.nreg_ndmax;
	QUICKASSERT( x >= xmi );
	QUICKASSERT( x <= xmx );
	QUICKASSERT( y >= ymi );
	QUICKASSERT( y <= ymx );
	QUICKASSERT( z >= zmi );
	QUICKASSERT( z <= zmx );
#endif

	uint32_t xx = this->mycaHdl->regions[rear]->myGeom.nreg_rd;
	uint32_t xxyy = this->mycaHdl->regions[rear]->myGeom.nregarea_rdtd;
	uint32_t cxyz = ((z - zmi) * xxyy) + ((y - ymi) * xx) + (x - xmi);

	return this->mycaHdl->regions[rear]->mycellgrid[cxyz];
}


uint32_t juncHdl::read_cellstate_my( int caller, uint32_t x, uint32_t y, uint32_t z )
{
	//CALLED FROM WITHIN PARALLEL REGION
	//reads mycellgrid state at global location x, y, z which is guaranteed to be the bottom (-y) memory region neighbor of caller
	QUICKASSERT( (y-1) >= 0 );

	uint32_t npy = this->mycaHdl->regions[caller]->thePartitioning.nreg_tdy;
	//int ix = 0;//get coordinates of regid
	int iy = caller % npy;
	int iz = caller / npy;

	iy += -1;		//front neighbor
	if ( iy < 0 )	//simplified, do not track boundary with periodic boundary conditions
		return INVALID_SEED;

	int front = ((iz * npy)+iy); //get logical neighbor and his geometry
	uint32_t xmi = this->mycaHdl->regions[front]->myGeom.nreg_rdmin;
	uint32_t ymi = this->mycaHdl->regions[front]->myGeom.nreg_tdmin;
	uint32_t zmi = this->mycaHdl->regions[front]->myGeom.nreg_ndmin;

#ifdef DEBUG
	uint32_t xmx = this->mycaHdl->regions[front]->myGeom.nreg_rdmax;
	uint32_t ymx = this->mycaHdl->regions[front]->myGeom.nreg_tdmax;
	uint32_t zmx = this->mycaHdl->regions[front]->myGeom.nreg_ndmax;
	QUICKASSERT( x >= xmi );
	QUICKASSERT( x <= xmx );
	QUICKASSERT( y >= ymi );
	QUICKASSERT( y <= ymx );
	QUICKASSERT( z >= zmi );
	QUICKASSERT( z <= zmx );
#endif

	uint32_t xx = this->mycaHdl->regions[front]->myGeom.nreg_rd;
	uint32_t xxyy = this->mycaHdl->regions[front]->myGeom.nregarea_rdtd;
	uint32_t cxyz = ((z - zmi) * xxyy) + ((y - ymi) * xx) + (x - xmi);

	return this->mycaHdl->regions[front]->mycellgrid[cxyz];
}


uint32_t juncHdl::read_cellstate_mz( int caller, uint32_t x, uint32_t y, uint32_t z )
{
	//CALLED FROM WITHIN PARALLEL REGION
	//reads mycellgrid state at global location x, y, z which is guaranteed to be the front (-z) memory region neighbor of caller
	QUICKASSERT( (z-1) >= 0 );
	uint32_t npy = this->mycaHdl->regions[caller]->thePartitioning.nreg_tdy;
	//int ix = 0;//get coordinates of regid
	int iy = caller % npy;
	int iz = caller / npy;

	iz += -1;		//bottom neighbor
	if ( iz < 0 )	//simplified periodic boundary conditions
		return INVALID_SEED;

	int bottom = ((iz * npy)+iy); //get logical neighbor and his geometry
	uint32_t xmi = this->mycaHdl->regions[bottom]->myGeom.nreg_rdmin;
	uint32_t ymi = this->mycaHdl->regions[bottom]->myGeom.nreg_tdmin;
	uint32_t zmi = this->mycaHdl->regions[bottom]->myGeom.nreg_ndmin;

#ifdef DEBUG
	uint32_t xmx = this->mycaHdl->regions[bottom]->myGeom.nreg_rdmax;
	uint32_t ymx = this->mycaHdl->regions[bottom]->myGeom.nreg_tdmax;
	uint32_t zmx = this->mycaHdl->regions[bottom]->myGeom.nreg_ndmax;
	QUICKASSERT( x >= xmi );
	QUICKASSERT( x <= xmx );
	QUICKASSERT( y >= ymi );
	QUICKASSERT( y <= ymx );
	QUICKASSERT( z >= zmi );
	QUICKASSERT( z <= zmx );
#endif

	uint32_t xx = this->mycaHdl->regions[bottom]->myGeom.nreg_rd;
	uint32_t xxyy = this->mycaHdl->regions[bottom]->myGeom.nregarea_rdtd;
	uint32_t cxyz = ((z - zmi) * xxyy) + ((y - ymi) * xx) + (x - xmi);

	return this->mycaHdl->regions[bottom]->mycellgrid[cxyz];
}


uint32_t juncHdl::read_cellstate_pz( int caller, uint32_t x, uint32_t y, uint32_t z )
{
	//CALLED FROM WITHIN PARALLEL REGION
	//reads mycellgrid state at global location x, y, z which is guaranteed to be the rear (+z) memory region neighbor of caller
	QUICKASSERT( (z+1) < this->mycaHdl->myCAGeometry.nboxedge_nd );
	uint32_t npy = this->mycaHdl->regions[caller]->thePartitioning.nreg_tdy;
	//int ix = 0;//get coordinates of regid
	int iy = caller % npy;
	int iz = caller / npy;

	iz += +1;		//rear neighbor
	if ( iz >= this->mycaHdl->regions[caller]->thePartitioning.nreg_ndz )	//simplified periodic boundary conditions
		return INVALID_SEED;

	int top = ((iz * npy)+iy); //get logical neighbor and his geometry
	uint32_t xmi = this->mycaHdl->regions[top]->myGeom.nreg_rdmin;
	uint32_t ymi = this->mycaHdl->regions[top]->myGeom.nreg_tdmin;
	uint32_t zmi = this->mycaHdl->regions[top]->myGeom.nreg_ndmin;

#ifdef DEBUG
	uint32_t xmx = this->mycaHdl->regions[top]->myGeom.nreg_rdmax;
	uint32_t ymx = this->mycaHdl->regions[top]->myGeom.nreg_tdmax;
	uint32_t zmx = this->mycaHdl->regions[top]->myGeom.nreg_ndmax;
	QUICKASSERT( x >= xmi );
	QUICKASSERT( x <= xmx );
	QUICKASSERT( y >= ymi );
	QUICKASSERT( y <= ymx );
	QUICKASSERT( z >= zmi );
	QUICKASSERT( z <= zmx );
#endif

	uint32_t xx = this->mycaHdl->regions[top]->myGeom.nreg_rd;
	uint32_t xxyy = this->mycaHdl->regions[top]->myGeom.nregarea_rdtd;
	uint32_t cxyz = ((z - zmi) * xxyy) + ((y - ymi) * xx) + (x - xmi);

	return this->mycaHdl->regions[top]->mycellgrid[cxyz];
}


void juncHdl::addVoxelGB( uint32_t c, uint32_t a, uint32_t b )
{
	uint32_t pmax = MAX(a,b);
	uint32_t pmin = MIN(a,b); //pmax == pmin is impossible because detection algorithm does not generate boundaries inside bulk grain volume
	uint_fast64_t hashtag = MYHASH(pmax, pmin); //get unique hash for identifying the face, higher IDs

	//generate a set of grain boundary faces for key grain pmax
	if ( this->fjuncSkeleton[pmax].maxLen == 0 ) { //not yet any boundary found for this grain?
		this->fjuncSkeleton[pmax].len = 0;
		this->fjuncSkeleton[pmax].maxLen = STDLEN;
		this->fjuncSkeleton[pmax].thefaces = NULL;
		this->fjuncSkeleton[pmax].thefaces = new bndFaceFast[STDLEN];
		QUICKASSERT( this->fjuncSkeleton[pmax].thefaces != NULL );
		this->juncMemGuard = this->juncMemGuard + ( STDLEN * sizeof(bndFaceFast) );

		//initialize for safety
		bndFaceFastP fc =  this->fjuncSkeleton[pmax].thefaces;
		for ( uint32_t af = 0; af <  this->fjuncSkeleton[pmax].maxLen; af++ ) {
			fc[af].voxelizedface = NULL;
			fc[af].nvoxeltotal = EMPTY;
			fc[af].nextfreevoxelslot = THEFIRSTONE;
		}
	}

	uint32_t i;
	//attempt to add the voxel to an existent face, so search the face
	for ( i = 0; i < this->fjuncSkeleton[pmax].len; i++ ) {
		if( this->fjuncSkeleton[pmax].thefaces[i].id == hashtag ) {

			//check whether enough memory to store the cell
			bndFaceFastP aface = &(this->fjuncSkeleton[pmax].thefaces[i]);
			if ( aface->nextfreevoxelslot < aface->nvoxeltotal ) {
				aface->voxelizedface[aface->nextfreevoxelslot].location = MARKER_TO_IDENTIFY_EXPELLED_SITES + c;
				//MK::technically c could be 0, namely coordinate 0,0,0 of the CA because c encodes a global implicit coordinate
				//however the marking resets all locations to +1, i.e. to values > 0
				//thereby we can easily flag cells we have already marked as nucleation sites during the subsequent choosing of GB nucleation sites by flagging them 0
				aface->voxelizedface[aface->nextfreevoxelslot].seedid = a; //the seed at that location
				aface->nextfreevoxelslot++;
				return;
			}

			//no there was not enough memory preallocated in this face buffer, so reallocate
			size_t oldSize = aface->nvoxeltotal;
			size_t newSize = aface->nvoxeltotal + DEFAULT_VOXEL2FACE_REFERENCES;
			cellsBndFastP newblock = NULL;
			newblock = new cellsBndFast[newSize];
			QUICKASSERT( newblock != NULL );
			this->juncMemGuard = this->juncMemGuard + ( newSize * sizeof(cellsBndFast) );

			//copy the old voxelpairs first
			cellsBndFastP oldblock = aface->voxelizedface;
			for ( uint32_t vxp = 0; vxp < oldSize; vxp++ ) {
				newblock[vxp].location = oldblock[vxp].location;
				newblock[vxp].seedid = oldblock[vxp].seedid;
			}
			//free the old block
			delete [] this->fjuncSkeleton[pmax].thefaces[i].voxelizedface;
			this->juncMemGuard = this->juncMemGuard - (oldSize * sizeof(cellsBndFast) );

			//link in the new block
			this->fjuncSkeleton[pmax].thefaces[i].voxelizedface = newblock;
			this->fjuncSkeleton[pmax].thefaces[i].nvoxeltotal = newSize;
			this->fjuncSkeleton[pmax].thefaces[i].nextfreevoxelslot++;

			//now take the first memory location of the newblock to store
			newblock[oldSize].location = MARKER_TO_IDENTIFY_EXPELLED_SITES + c;
			newblock[oldSize].seedid = a;
			return;
		}
	}

	//mhh we are initializing a completely new face
	//do we have still memory to open a bucket to fill cell positions?
	if( i < this->fjuncSkeleton[pmax].maxLen ) {
		this->fjuncSkeleton[pmax].thefaces[i].gposUp = pmax;
		this->fjuncSkeleton[pmax].thefaces[i].gposDown = pmin;
		this->fjuncSkeleton[pmax].thefaces[i].id = hashtag;
		this->fjuncSkeleton[pmax].thefaces[i].disori = GBDISORI_NOT_REQUIRED_YET;
		this->fjuncSkeleton[pmax].thefaces[i].gbID = GBNAME_UNKNOWN;
		this->fjuncSkeleton[pmax].thefaces[i].voxelizedface = NULL;
		this->fjuncSkeleton[pmax].thefaces[i].voxelizedface = new cellsBndFast[DEFAULT_VOXEL2FACE_REFERENCES];
		QUICKASSERT( this->fjuncSkeleton[pmax].thefaces[i].voxelizedface != NULL );
		this->juncMemGuard = this->juncMemGuard + (DEFAULT_VOXEL2FACE_REFERENCES * sizeof(cellsBndFast));
		this->fjuncSkeleton[pmax].thefaces[i].nvoxeltotal = DEFAULT_VOXEL2FACE_REFERENCES;
		this->fjuncSkeleton[pmax].thefaces[i].nextfreevoxelslot = THEFIRSTONE;

		this->fjuncSkeleton[pmax].len++;
		this->fjuncSkeleton[pmax].thefaces[i].voxelizedface[0].location = MARKER_TO_IDENTIFY_EXPELLED_SITES + c;
		this->fjuncSkeleton[pmax].thefaces[i].voxelizedface[0].seedid = a;
		this->fjuncSkeleton[pmax].thefaces[i].nextfreevoxelslot++;

		return;
	}

	QUICKASSERT ( this->fjuncSkeleton[pmax].len == this->fjuncSkeleton[pmax].maxLen );

	//not we do not even have place for a face, so we add further faces into the collection of pmax
	size_t oSize = this->fjuncSkeleton[pmax].maxLen;
	size_t nSize = oSize + STDLEN;

	bndFaceFastP newfaces = NULL;
	newfaces = new bndFaceFast[nSize];
	QUICKASSERT( newfaces != NULL );
	this->juncMemGuard = this->juncMemGuard + ( nSize * sizeof(bndFaceFast) );

	//copy the old face data first
	bndFaceFastP ofaces = this->fjuncSkeleton[pmax].thefaces;
	for ( uint32_t oldf = 0; oldf < oSize; oldf++ ) {
		newfaces[oldf].gposUp = ofaces[oldf].gposUp;
		newfaces[oldf].gposDown = ofaces[oldf].gposDown;
		newfaces[oldf].id = ofaces[oldf].id;
		newfaces[oldf].disori = ofaces[oldf].disori;
		newfaces[oldf].gbID = ofaces[oldf].gbID;
		newfaces[oldf].voxelizedface = ofaces[oldf].voxelizedface; //carry over the pointer to the cell container building this face to avoid leaks
		newfaces[oldf].nvoxeltotal = ofaces[oldf].nvoxeltotal;
		newfaces[oldf].nextfreevoxelslot = ofaces[oldf].nextfreevoxelslot;
	}
	//initialize the additional potential face buckets
	for ( uint32_t newf = oSize; newf < nSize; newf++ ) {
		newfaces[newf].voxelizedface = NULL;
		newfaces[newf].nvoxeltotal = EMPTY;
		newfaces[newf].nextfreevoxelslot = THEFIRSTONE;
	}
	//delete old buckets and link in new ones
	delete [] this->fjuncSkeleton[pmax].thefaces; //the memory that was referred to in thefaces[..].voxelizedface was copied over
	this->juncMemGuard = this->juncMemGuard - (oSize * sizeof(bndFaceFast));
	this->fjuncSkeleton[pmax].thefaces = newfaces;
	this->fjuncSkeleton[pmax].maxLen = nSize;

	this->fjuncSkeleton[pmax].len++;

	//finally make up new face
	newfaces[i].gposUp = pmax;
	newfaces[i].gposDown = pmin;
	newfaces[i].id = hashtag;
	newfaces[i].disori = GBDISORI_NOT_REQUIRED_YET;
	newfaces[i].gbID = GBNAME_UNKNOWN;
	newfaces[i].voxelizedface = NULL;
	newfaces[i].voxelizedface = new cellsBndFast[DEFAULT_VOXEL2FACE_REFERENCES];
	QUICKASSERT( newfaces[i].voxelizedface != NULL );
	this->juncMemGuard = this->juncMemGuard + (DEFAULT_VOXEL2FACE_REFERENCES * sizeof(cellsBndFast));

	newfaces[i].nvoxeltotal = DEFAULT_VOXEL2FACE_REFERENCES;
	newfaces[i].nextfreevoxelslot = THEFIRSTONE;

	//add our cell to this new face
	newfaces[i].voxelizedface[0].location = MARKER_TO_IDENTIFY_EXPELLED_SITES + c;
	newfaces[i].voxelizedface[0].seedid = a;
	newfaces[i].nextfreevoxelslot++;
}


void juncHdl::addVoxelTR( uint32_t c )
{
	if ( this->tjuncSkeletonSize != 0 ) { //most likely we desire adding coordinates to an already allocated array
		if ( this->tjuncnextfreeslot < this->tjuncSkeletonSize ) {
			//write into preallocated memory
			this->tjuncSkeleton[this->tjuncnextfreeslot].location = MARKER_TO_IDENTIFY_EXPELLED_SITES + c;
			this->tjuncnextfreeslot++;
			return;
		}

		//end of preallocated buffer reached, reallocation necessary
		QUICKASSERT( this->tjuncnextfreeslot == this->tjuncSkeletonSize );
		size_t oldSize = this->tjuncSkeletonSize;
		size_t newSize = this->tjuncSkeletonSize + DEFAULT_ALLOCATE_TJ;
		cellsTJ* newBlock = NULL;
		newBlock = new cellsTJ[newSize];
		QUICKASSERT( newBlock != NULL );
		this->juncMemGuard = this->juncMemGuard + (newSize * sizeof(cellsTJ));

		//copy over old block
		for ( size_t i = 0; i < oldSize; i++ ) {
			newBlock[i].location = this->tjuncSkeleton[i].location;
		}
		//delete old block and repointer
		delete [] this->tjuncSkeleton;
		this->juncMemGuard = this->juncMemGuard - (oldSize * sizeof(cellsTJ));
		this->tjuncSkeleton = newBlock;
		this->tjuncSkeletonSize = this->tjuncSkeletonSize + DEFAULT_ALLOCATE_TJ;

		//use new memory now
		this->tjuncSkeleton[this->tjuncnextfreeslot].location = MARKER_TO_IDENTIFY_EXPELLED_SITES + c;
		this->tjuncnextfreeslot++;
		return;
	}

	//first time allocation ever
	QUICKASSERT( this->tjuncSkeleton == NULL );

	cellsTJ* newBlock = NULL;
	newBlock = new cellsTJ[DEFAULT_ALLOCATE_TJ];
	QUICKASSERT( newBlock != NULL );
	this->juncMemGuard = this->juncMemGuard + (DEFAULT_ALLOCATE_TJ * sizeof(cellsTJ));

	this->tjuncSkeleton = newBlock;
	this->tjuncSkeletonSize = this->tjuncSkeletonSize + DEFAULT_ALLOCATE_TJ;

	this->tjuncSkeleton[this->tjuncnextfreeslot].location = MARKER_TO_IDENTIFY_EXPELLED_SITES + c;
	this->tjuncnextfreeslot++;
}


void juncHdl::addVoxelHT( uint32_t c )
{
	if ( this->hjuncSkeletonSize != 0 ) { //most likely we desire adding coordinates to already allocated array
		if ( this->hjuncnextfreeslot < this->hjuncSkeletonSize ) {
			//write into existent
			this->hjuncSkeleton[this->hjuncnextfreeslot].location = MARKER_TO_IDENTIFY_EXPELLED_SITES + c;
			this->hjuncnextfreeslot++;
			return;
		}

		//end of preallocated buffer reached, reallocation necessary
		QUICKASSERT( this->hjuncnextfreeslot == this->hjuncSkeletonSize );
		size_t oldSize = this->hjuncSkeletonSize;
		size_t newSize = this->hjuncSkeletonSize + DEFAULT_ALLOCATE_HJ;
		cellsHJ* newBlock = NULL;
		newBlock = new cellsHJ[newSize];
		QUICKASSERT( newBlock != NULL );
		this->juncMemGuard = this->juncMemGuard + (newSize * sizeof(cellsHJ));

		//copy over old block
		for ( size_t i = 0; i < oldSize; i++ ) {
			newBlock[i].location = this->hjuncSkeleton[i].location;
		}
		//delete old block and repointer
		delete [] this->hjuncSkeleton;
		this->juncMemGuard = this->juncMemGuard - (oldSize * sizeof(cellsHJ));
		this->hjuncSkeleton = newBlock;
		this->hjuncSkeletonSize = this->hjuncSkeletonSize + DEFAULT_ALLOCATE_HJ;

		//use new memory now
		this->hjuncSkeleton[this->hjuncnextfreeslot].location = MARKER_TO_IDENTIFY_EXPELLED_SITES + c;
		this->hjuncnextfreeslot++;
		return;
	}

	//first time allocation ever
	QUICKASSERT( this->hjuncSkeleton == NULL );

	cellsHJ* newBlock = NULL;
	newBlock = new cellsHJ[DEFAULT_ALLOCATE_HJ];
	QUICKASSERT( newBlock != NULL );
	this->juncMemGuard = this->juncMemGuard + (DEFAULT_ALLOCATE_HJ * sizeof(cellsHJ));

	this->hjuncSkeleton = newBlock;
	this->hjuncSkeletonSize = this->hjuncSkeletonSize + DEFAULT_ALLOCATE_HJ;

	this->hjuncSkeleton[this->hjuncnextfreeslot].location = MARKER_TO_IDENTIFY_EXPELLED_SITES + c;
	this->hjuncnextfreeslot++;
}


void juncHdl::updateFaces( uint32_t p, vector<organizeBnd>* luu )
{
	//CALLED FROM WITHIN PARALLEL REGION
	//resets memory for storing faces in this->fjuncSkeleton to nf
	//add further faces into that collection
	if ( this->fjuncSkeleton[p].maxLen != luu->size() ) { //either too many buckets were already allocated or too few, either way set accordingly...
		size_t oldSize = this->fjuncSkeleton[p].maxLen;
		size_t usedSize = this->fjuncSkeleton[p].len; //always len <= maxLen holds
		size_t newSize = luu->size();

		bndFaceFastP newfaces = NULL; //get memory for these faces
		newfaces = new bndFaceFast[newSize];
		QUICKASSERT( newfaces != NULL );
		this->juncMemGuard = this->juncMemGuard + ( newSize * sizeof(bndFaceFast) );

		//copy old utilized portion over first
		bndFaceFastP oldfaces = this->fjuncSkeleton[p].thefaces;
		for ( uint32_t oldf = 0; oldf < usedSize; oldf++ ) { //may not execute at all if calling thread has no boundary segments with p as high index found
			newfaces[oldf].gposUp = oldfaces[oldf].gposUp;
			newfaces[oldf].gposDown = oldfaces[oldf].gposDown;
			newfaces[oldf].id = oldfaces[oldf].id;
			QUICKASSERT( luu->at(oldf).key == newfaces[oldf].gposDown );
			newfaces[oldf].disori = oldfaces[oldf].disori;
			newfaces[oldf].gbID = oldfaces[oldf].gbID;
			newfaces[oldf].voxelizedface = oldfaces[oldf].voxelizedface; //carry over the pointer to the cell container build this face to avoid leaks
			newfaces[oldf].nvoxeltotal = oldfaces[oldf].nvoxeltotal;
			newfaces[oldf].nextfreevoxelslot = oldfaces[oldf].nextfreevoxelslot;
		}
		//initialize additional faces
		for ( uint32_t newf = usedSize; newf < newSize; newf++ ) {
			newfaces[newf].gposUp = p;
			newfaces[newf].gposDown = luu->at(newf).key;
			newfaces[newf].id = MYHASH( p, luu->at(newf).key );
			newfaces[newf].disori = GBDISORI_NOT_REQUIRED_YET;
			newfaces[newf].gbID = GBNAME_UNKNOWN;

			newfaces[newf].voxelizedface = NULL;
			newfaces[newf].nvoxeltotal = EMPTY;
			newfaces[newf].nextfreevoxelslot = THEFIRSTONE;
			this->fjuncSkeleton[p].len++;
		}

		//delete old buckets and link in new ones
		delete [] this->fjuncSkeleton[p].thefaces; //the memory that was referred to in thefaces[..].voxelizedface was copied over
		this->juncMemGuard = this->juncMemGuard - (oldSize * sizeof(bndFaceFast));
		this->fjuncSkeleton[p].thefaces = newfaces;
		this->fjuncSkeleton[p].maxLen = newSize;
		QUICKASSERT( this->fjuncSkeleton[p].len == this->fjuncSkeleton[p].maxLen );
	}
}


void juncHdl::registerMissingFaces( uint32_t p, vector<organizeBnd>* luu ) 
{
	//CALLED FROM WITHIN PARALLEL REGION
	uint32_t len = this->fjuncSkeleton[p].len;
	uint32_t maxlen = this->fjuncSkeleton[p].maxLen;
	QUICKASSERT( len == luu->size() );
	QUICKASSERT( maxlen == luu->size() );
	QUICKASSERT( maxlen != 0 );

	//check consistency with luu
	for ( uint32_t f = 0; f < maxlen; f++ ) {
		QUICKASSERT( f == luu->at(f).poslocal );
		QUICKASSERT( this->fjuncSkeleton[p].thefaces[f].gposUp == p );
		QUICKASSERT( this->fjuncSkeleton[p].thefaces[f].gposDown == luu->at(f).key );
		QUICKASSERT( this->fjuncSkeleton[p].thefaces[f].id == MYHASH(p, luu->at(f).key) );
		QUICKASSERT( this->fjuncSkeleton[p].thefaces[f].disori == GBDISORI_NOT_REQUIRED_YET );
		QUICKASSERT( this->fjuncSkeleton[p].thefaces[f].gbID == GBNAME_UNKNOWN );
	}
}


void juncHdl::updateFaceCellBuffers( uint32_t p, vector<organizeBnd>* luu )
{
	//CALLED FROM WITHIN PARALLEL REGION
	QUICKASSERT( this->fjuncSkeleton[p].maxLen == luu->size() ); //checks correct length of threadlocal face container

	for ( uint32_t f = 0; f < this->fjuncSkeleton[p].maxLen; f++ ) { 
		//now as all faces are registered we may require further memory to store individual cell positions because the calling thread may have already found
		//all cells, only a few, none for each face i.e. pmax/pmin pairing
		QUICKASSERT( this->fjuncSkeleton[p].thefaces[f].gposDown == luu->at(f).key );

		uint32_t oldSize = this->fjuncSkeleton[p].thefaces[f].nvoxeltotal; //either EMPTY, i.e. 0 or > EMPTY
		uint32_t usedSize = this->fjuncSkeleton[p].thefaces[f].nextfreevoxelslot;
		uint32_t newSize = luu->at(f).cnt;
		QUICKASSERT( newSize != 0 );

		//most likely oldSize != newSize
		if ( oldSize != newSize ) {
			cellsBndFast* ac = NULL;
			ac = new cellsBndFast[newSize];
			this->juncMemGuard = this->juncMemGuard + ( newSize * sizeof(cellsBndFast));

			for ( uint32_t i = 0; i < usedSize; i++ ) {
				ac[i].location = this->fjuncSkeleton[p].thefaces[f].voxelizedface[i].location;
				ac[i].seedid = this->fjuncSkeleton[p].thefaces[f].voxelizedface[i].seedid;
			}

			delete [] this->fjuncSkeleton[p].thefaces[f].voxelizedface;
			this->juncMemGuard = this->juncMemGuard - (oldSize * sizeof(cellsBndFast));

			this->fjuncSkeleton[p].thefaces[f].voxelizedface = ac;
			this->fjuncSkeleton[p].thefaces[f].nvoxeltotal = newSize;
			this->fjuncSkeleton[p].thefaces[f].nextfreevoxelslot = usedSize;
			QUICKASSERT( this->fjuncSkeleton[p].thefaces[f].nvoxeltotal == luu->at(f).cnt );
			//section [nextfreevoxelslot,nvoxeltotal) is uninitialized, i.e. random values
		}
	}
}


void juncHdl::addRemoteFaceCells( uint32_t p, uint32_t targetthread, vector<organizeBnd>* luu ) 
{
//now collect pieces of information from other threads BUT THE SAME SEED ID!
	uint32_t nthreads = omp_get_num_threads(); //int 2 uint32 safe because >=1 and not in the near future >INT32_RANGE

	for ( uint32_t t = MASTER; t < nthreads; t++ ) {
		if ( t != targetthread ) { //querying of remote cell positions necessary
			for ( unsigned int f = 0; f < luu->size(); f++ ) {
				juncHdlP remotejunc = &(this->mycaHdl->regions[t]->myjunctions);

				for ( uint32_t cand = 0; cand < remotejunc->fjuncSkeleton[p].len; cand++ ) { //wont be executed if len == 0
					if ( luu->at(f).key == remotejunc->fjuncSkeleton[p].thefaces[cand].gposDown ) { //##MK::micoropt potential when checking key != 
						//copy all cells from this face
						QUICKASSERT( p == remotejunc->fjuncSkeleton[p].thefaces[cand].gposUp );
						QUICKASSERT( &(this->fjuncSkeleton[p].thefaces[luu->at(f).poslocal].voxelizedface) == &(this->mycaHdl->regions[targetthread]->myjunctions.fjuncSkeleton[p].thefaces[luu->at(f).poslocal].voxelizedface) );

						cellsBndFast* target = this->fjuncSkeleton[p].thefaces[luu->at(f).poslocal].voxelizedface;
						uint32_t offset = this->fjuncSkeleton[p].thefaces[luu->at(f).poslocal].nextfreevoxelslot;
						cellsBndFast* source = remotejunc->fjuncSkeleton[p].thefaces[cand].voxelizedface;
						uint32_t nsource = remotejunc->fjuncSkeleton[p].thefaces[cand].nextfreevoxelslot;
						for ( uint32_t i = 0; i < nsource; i++ ) {
							target[offset+i].location = source[i].location;
							target[offset+i].seedid = source[i].seedid;
						}
						this->fjuncSkeleton[p].thefaces[luu->at(f).poslocal].nextfreevoxelslot = this->fjuncSkeleton[p].thefaces[luu->at(f).poslocal].nextfreevoxelslot + nsource;
						break; //because we can consolidate the next face
					}
				}
			}
		} //evaluate potential data from next thread
	}
}


void juncHdl::ompshar_init_junctions( void )
{
	//CALLED FROM WITHIN PARALLEL REGIONS
	this->fjuncSkeletonSize = this->mycaHdl->ndefgseeds;
	this->fjuncSkeleton = NULL;
	this->fjuncSkeleton = new bndColumnFast[this->fjuncSkeletonSize];
	QUICKASSERT ( this->fjuncSkeleton != NULL );
	this->juncMemGuard = this->juncMemGuard + (this->fjuncSkeletonSize * sizeof(bndColumnFast));

	for( uint32_t pmax = 0; pmax < this->fjuncSkeletonSize; pmax++ ) { //no faces identified
		this->fjuncSkeleton[pmax].len = 0;
		this->fjuncSkeleton[pmax].maxLen = 0;
		this->fjuncSkeleton[pmax].thefaces = NULL;
	}
}


void juncHdl::ompshar_detect_junctions( void )
{
	//MK::scanning for adjacencies of a voxel/cells to differently IDed grains in von Neumann neighborhood
	//MK::the old method was slightly inaccurate as it always scanned for boundaries along a non-randomized probing order
	//MK::a more accurate approach would probe all neighbors, collect which and how many of them are different and
	//with this information skeletonize the network into pure faces (where there are only two neighbors), tri junctions,
	//and higher order junctions, i.e. quadruple points and alike as in the discrete picture of a CA octuple points are possible
	caHdlP mca = this->mycaHdl;
	uint32_t* thecells = this->mymemreg->mycellgrid;

	int xmi = this->mymemreg->myGeom.nreg_rdmin; //local region volume in global coordinates uint32_t --> int no problem because rvalue > 0 && <= CA_MAXIMUM_DIMENSIONS << INT32_RANGE
	int xmx = this->mymemreg->myGeom.nreg_rdmax;
	int ymi = this->mymemreg->myGeom.nreg_tdmin;
	int ymx = this->mymemreg->myGeom.nreg_tdmax;
	int zmi = this->mymemreg->myGeom.nreg_ndmin;
	int zmx = this->mymemreg->myGeom.nreg_ndmax;
	int xx = this->mymemreg->myGeom.nreg_rd; //this is positive and within INT32 range because the domain is also at least in z CA_MINIMUM_DIMENSIONS in length
	int xxyy = this->mymemreg->myGeom.nregarea_rdtd;

	int gxmx = mca->myCAGeometry.nboxedge_rd; //SU global only [0,gimx-1] is still inside the domain
	int gymx = mca->myCAGeometry.nboxedge_td;
	int gzmx = mca->myCAGeometry.nboxedge_nd;
	int gxy = mca->myCAGeometry.nboxarea_rdtd;

	uint32_t threeseeds[3] = { INVALID_SEED, INVALID_SEED, INVALID_SEED };
	uint32_t character = 0;
	uint32_t currseed, testseed;

	for ( uint32_t z = zmi; z <= zmx; z++ ) { //imi and imx are global coordinates!
		//##MK::further optimization could eliminate periodic boundary checks by filtering e.g. if ( z != zmi && z != zmx ) { //so z inside the domain else { //border }
		for ( uint32_t y = ymi; y <= ymx; y++ ) {
			for ( uint32_t x = xmi; x <= xmx; x++ ) {
				//find at least three additional disjoint and valid keys (seed of deformed grain) over the von Neumann neighbors
				//then we know the cell at global coordinate x,y,z, is adjacent to a higher order junction (>3)
				//first of all get seed at x,y,z
				currseed = INVALID_SEED;
				currseed = thecells[(x-xmi)+(y-ymi)*xx+(z-zmi)*xxyy]; //global to local coordinate inside the memory region
				if ( currseed != INVALID_SEED ) { //assure cell is no CELL_IS_PARTICLE
					//find at most three disjoint and valid seeds disjoint to currseed
					threeseeds[0] = INVALID_SEED;
					threeseeds[1] = INVALID_SEED;
					threeseeds[2] = INVALID_SEED;
					character = 0;

					//start search for at most three disjoint seeds
					testseed = INVALID_SEED; //-x
					//if ( (int) (x-1) > -1 ) { testseed = thecells[(x-1-xmi)+(y-ymi)*xx+(z-zmi)*xxyy]; } without pay for the costly explicit cast and add
					if ( x != 0 ) { //not at the global domain boundary, THEREFORE AT THE MOMENT NO PERIODIC TRACKING OF GB BOUNDARIES!
						if ( x > xmi ) //also not at the local domain boundary
							testseed = thecells[(x-1-xmi)+(y-ymi)*xx+(z-zmi)*xxyy];
						//else //data are in another threads memory region, ##MK::but we do not partition regions along X
						//	testseed = read_cellstate_mx( this->mymemreg->mythreadid, (x-1), y, z );
					} //else x == 0 is outside SU cube but we do not track boundaries with PBC therefore leave testseed as INVALID
					//returns either INVALID_SEED or integer > -1
					//at most::character == 0,threeseeds [000] --> 0, [000] with 1 -valid, 0 - false, furthermore testseed != currseed --> INVALID or valid and disjoint
					if ( testseed != currseed && testseed != threeseeds[0] ) { //i.e. when different to currseed and INVALID_SEED
						//MK::as initially character == 0 and all entries of threeseeds were set to INVALID we do not get here if testseed were INVALID
						threeseeds[character] = testseed; //at most [000] --> [100]
						character++; //at most::1, i.e. we have at most found one disjoint seed
					}
					testseed = INVALID_SEED; //+x
					if ( (x+1) < gxmx ) {
						if ( x < xmx ) 
							testseed = thecells[(x+1-xmi)+(y-ymi)*xx+(z-zmi)*xxyy];
						//else, ##MK::same as above
						//	testseed = read_cellstate_px( this->mymemreg->mythreadid, (x+1), y, z );
					}
					//1, [100]
					if ( testseed != currseed && testseed != threeseeds[0] && testseed != threeseeds[1] ) {
						//so if testseed is different from currseed different to the possibly already found threeseeds[0] and not INVALID
						threeseeds[character] = testseed; //at most [100] --> [110]
						character++; //2
					}
					testseed = INVALID_SEED; //-y
					if ( y != 0  ) {
						if ( y > ymi ) //staying in myself
							testseed = thecells[(x-xmi)+(y-1-ymi)*xx+(z-zmi)*xxyy];
						else
							testseed = read_cellstate_my( this->mymemreg->mythreadid, x, (y-1), z );
					}
					//by now character can be at most 2, [110]
					if ( testseed != currseed && testseed != threeseeds[0] && testseed != threeseeds[1] && testseed != threeseeds[2] ) {
						threeseeds[character] = testseed; //at most [110] --> [111]
						character++; //3
					}
					//by now we may already have found 3 disjoint seeds additional to currentseed so we only have to continue if not
					if ( character < 3 ) { //3, [111]
						//if this executes, character is < 3 and for sure [110] so access to threeseeds still valid and including a test for INVALID_SEED with the last condition
						testseed = INVALID_SEED; //+y
						if ( (y+1) < gymx ) { 
							if ( y < ymx )
								testseed = thecells[(x-xmi)+(y+1-ymi)*xx+(z-zmi)*xxyy];
							else
								testseed = read_cellstate_py( this->mymemreg->mythreadid, x, (y+1), z );
						}
						if ( testseed != currseed && testseed != threeseeds[0] && testseed != threeseeds[1] && testseed != threeseeds[2] ) {
							//even though < 3 we have to test against all threeseeds because by now character could be 2 so we have to accept
							//only the testseed if it is different to currseed, the two first of the threeseeds and additionally also not an invalid key
							threeseeds[character] = testseed;
							character++;
						}
					}
					if ( character < 3 ) {
						testseed = INVALID_SEED; //-z
						if ( z != 0  ) { 
							if ( z > zmi )
								testseed = thecells[(x-xmi)+(y-ymi)*xx+(z-1-zmi)*xxyy];
							else 
								testseed = read_cellstate_mz( this->mymemreg->mythreadid, x, y, (z-1) );
						}
						if ( testseed != currseed && testseed != threeseeds[0] && testseed != threeseeds[1] && testseed != threeseeds[2] ) {
							threeseeds[character] = testseed;
							character++;
						}
					}
					if ( character < 3 ) {
						testseed = INVALID_SEED; //+z
						if ( (z+1) < gzmx ) {
							if ( z < zmx ) //as zmx is inclusive z+1 can at most be zmx so we can still access within our local data
								testseed = thecells[(x-xmi)+(y-ymi)*xx+(z+1-zmi)*xxyy];
							else
								testseed = read_cellstate_pz( this->mymemreg->mythreadid, x, y, (z+1) );
						}
						if ( testseed != currseed && testseed != threeseeds[0] && testseed != threeseeds[1] && testseed != threeseeds[2] ) {
							threeseeds[character] = testseed;
							character++;
						}
					}
					//categorize character of the voxel, mutally exclusive
					if ( character == 1 ) {
						this->addVoxelGB( (x+y*gxmx+z*gxy), currseed, threeseeds[0] ); //we require a global coordinate for the boundary bucket!
						continue;
					}
					if ( character == 2 ) {
						this->addVoxelTR( (x+y*gxmx+z*gxy) ); //##MK::add monitoring of the IDs if desired
						continue;
					}
					if ( character == 3 ) {
						this->addVoxelHT( (x+y*gxmx+z*gxy) );
						continue;
					}
				} //done with the valid voxel
			}
		}
	}

	#pragma omp critical
	{
		cout << "\t\tThread ID " << omp_get_thread_num() << " detected junctions!" << endl;
	}
}

/*switch ( character ) //interpret voxel character
{
	case 0 : //none of the neighbors valid or all bulk so nothing to do
		continue;
	case 1 : //only one valid but disjoint neighbor
		addVoxelGB( x, y, z, currseed, dsjseed1 );
		continue;
	case 2 : //two valid and disjoint
		addVoxelTR( x, y, z, currseed ); //##MK::add monitoring of the IDs
		continue;
	case 3 : //three valid and disjoint
		addVoxelHI( x, y, z, currseed );
		continue;
}*/

void juncHdl::ompshar_report_junctions( void ) 
{
	//##MK::output all my junctions via pmax;pmi;poslocal;cnt
}


void juncHdl::ompshar_consolidate_junctions( void )
{
	//MK::now each memory region maintains an array of faces and other junctions, thereby grain boundary faces may rest in several
	//thread regions at the same time, how can one choose randomly from these? we have several options with their own pros and cons...
	//1) copy all cells to one thread and let it do the work: -copy over, -NUMA issues due to memory size, -NUMA issues during placement, +no sync during rx
	//2) send list of candidates, thread evaluate: -WE MAY LOOSE SITES BECAUSE THEY ARE OCCUPIED -- REQUIRES FEEDBACK
	//3) concatenate the lists only virtually and have the threads carry an concatenated list for each face
	//4) have the threads fire random infections from their lists, however ever the list may be very short thus we may run in finite size issues...
	//5) copy all cells of a grain to one thread

	//MK::the current compromise is to collect in thread omp_get_thread_num all seed indices for which PMAX % OMP_NUM_THREADS == omp_get_thread_num()
	//this may still incur load imbalances, because not all grain boundaries have the same count of cells and not each thread has the same local interfacial area density
	//however also partitioning work based on the local position may incur problems because all grain boundary interfaces may cluster in one memory region
	//with the effect of having all work accruing in one thread('s memory), therefore this compromise

	//CALLED FROM WITHIN A PARALLEL REGION on previously synchronized fjuncSkeleton which all threads can access via shared memory
	//scan all seed ids

	uint32_t mytid = omp_get_thread_num(); //int 2 uint32 safe because thread ids > -1 and << INT32_RANGE
	uint32_t nthreads = omp_get_num_threads(); //int 2 uint32 safe because >=1 and not in the near future >INT32_RANGE
	uint32_t target = MASTER;

	caHdlP myca = this->mycaHdl; //in all threads this points to the same mycaHdl class object instance
	vector<organizeBnd>* lu = NULL;

	//##MK::what happens if nthreads == 1 ?, then set size becomes essentially a trim only
	for ( uint32_t pmax = 0; pmax < this->fjuncSkeletonSize; pmax++ ) {
		target = workPartitioning(pmax, nthreads);
		//MK::static work decomposition which assures that threads do not invalidate each others data
		if ( target == mytid ) { //thread scans which boundary faces exists and how many cells are necessary in total per pmax,pmin hashtag
			//collect all faces with id max pmax from other threads and myself
			lu = new vector<organizeBnd>;
			QUICKASSERT( lu != NULL );
			uint32_t posloc = 0;
			//first add my own pairs pmax/pmin
			for ( uint32_t f = 0; f < myca->regions[target]->myjunctions.fjuncSkeleton[pmax].len; f++ ) { //may add nothing if len = 0
				//by construction of the faces pmi ID are mutually exclusive
				struct organizeBnd akey; 
				akey.key = myca->regions[target]->myjunctions.fjuncSkeleton[pmax].thefaces[f].gposDown;
				akey.cnt = myca->regions[target]->myjunctions.fjuncSkeleton[pmax].thefaces[f].nextfreevoxelslot;
				akey.poslocal = f;
				lu->push_back( akey );
			}
			posloc = myca->regions[target]->myjunctions.fjuncSkeleton[pmax].len;
			//sweep through all other threads to find further faces logically assigned to pmax or additional cells that require association with any of my pairs pmax/pmin
			for ( uint32_t t = MASTER; t < nthreads; t++ ) { 
				if ( t != target ) { //MK::query remote information
					for ( uint32_t f = 0; f < myca->regions[t]->myjunctions.fjuncSkeleton[pmax].len; f++ ) { //does not execute when len == 0
						uint32_t pmi = myca->regions[t]->myjunctions.fjuncSkeleton[pmax].thefaces[f].gposDown;
						uint32_t nc = myca->regions[t]->myjunctions.fjuncSkeleton[pmax].thefaces[f].nextfreevoxelslot;
						QUICKASSERT( nc > 0 );
						bool NewFace = true; //assume its a hastag unknown to the target thread
						for ( unsigned int k = 0; k < lu->size(); k++ ) {
							if ( lu->at(k).key == pmi ) {
								NewFace = false;
								lu->at(k).cnt += nc;
								break;
							}
						}
						if ( NewFace == true ) {
							struct organizeBnd akey;
							akey.key = pmi;
							akey.cnt = nc;
							akey.poslocal = posloc;
							lu->push_back( akey );
							posloc++;
						} //boundary pairs are ordered
					}
				} //inspect all pairs on thread t
			}
			//by now if lu->size() > 0 we know that the first entries refer to the exact order of pmi ids in the target's thread container | all which follow belong to other ids
			//only continue if there are at all faces for the grain, this has to be tested even though for sure all grains have physically faces, logically a face 
			//is associated always with the grain with the higher ID pmax, therefore grain 1 has for sure logically no boundaries as all neighbors (2,3,4,5) have higher indices
			if ( lu->size() == 0 ) { //nothing found neither in my database nor in those of the other threads...
				delete lu;
				lu = NULL;
			}
			else { //so if there are faces across some threads (not necessarily me!)
				//first logically organize the faces, as the local thread may already have identified some of them...fjuncSkeleton[pmax].thefaces but not all, if any!
				//that reset in a different region of the solitary unit domain, thus we need to prepare a face container of sufficient length
				this->updateFaces( pmax, lu );
				this->registerMissingFaces( pmax, lu );
				this->updateFaceCellBuffers( pmax, lu );
				this->addRemoteFaceCells( pmax, mytid, lu ); //heavily querying remote information

				delete lu; //##MK::microoptim potential one could swap and reutilize the buffer...
				lu = NULL;
			}
		} //thread has consolidated all pieces of information for seed id pmax
		//else memory needs to be freed after all threads have synchronized via the juncHdl::ompshar_cleanup_junctions function
	}

	#pragma omp critical
	{
		cout << "\t\tThread ID " << omp_get_thread_num() << " consolidated junctions!" << endl;
	}
	//##MK::currently we do not consolidate t and h junctions, this may be added in the future
}


void juncHdl::ompshar_calc_junction_properties( void )
{
	//calculates properties of junctions
	//CALLED FROM WITHIN PARALLEL REGION
	uint32_t mytid = omp_get_thread_num();
	uint32_t nthreads = omp_get_num_threads();
	caHdlP myca = this->mycaHdl;

	for ( uint32_t pmax = 0; pmax < this->fjuncSkeletonSize; pmax++ ) {
		//NEEDS TO BE THE SAME STATIC WORK PARTITIONING AS IN CONSOLIDATE_JUNCTIONS
		if ( workPartitioning(pmax, nthreads) == mytid ) { //if OMP_NUM_THREADS=1 concerns all faces
			for ( uint32_t f = 0; f < myca->regions[mytid]->myjunctions.fjuncSkeleton[pmax].len; f++ ) {
				uint32_t pup = myca->regions[mytid]->myjunctions.fjuncSkeleton[pmax].thefaces[f].gposUp;
				uint32_t pdw = myca->regions[mytid]->myjunctions.fjuncSkeleton[pmax].thefaces[f].gposDown;

				myca->regions[mytid]->myjunctions.fjuncSkeleton[pmax].thefaces[f].disori = this->mycaHdl->calculateBoundaryDisoriFast( pup, pdw );
			}
		}
	}
}


void juncHdl::ompshar_cleanup_foreign_junctions( void )
{
	//CALLED FROM WITHIN PARALLEL REGION
	//thread cleans all memory to boundaries which they do not took care of as these are stored in the memory of other threads
	uint32_t mytid = omp_get_thread_num(); //int 2 uint32 possible because thread ids > -1 and << INT32_RANGE
	uint32_t nthreads = omp_get_num_threads();
	caHdlP myca = this->mycaHdl;

	for ( uint32_t pmax = 0; pmax < this->fjuncSkeletonSize; pmax++ ) {
		//NEEDS TO BE THE SAME WORK PARTITIONING AS IN CONSOLIDATE JUNCTIONS
		if ( workPartitioning(pmax, nthreads) != mytid ) { //if OMP_NUM_THREADS=1 not executed
			if ( myca->regions[mytid]->myjunctions.fjuncSkeleton[pmax].maxLen > 0 ) {
				for ( uint32_t f = 0; f < myca->regions[mytid]->myjunctions.fjuncSkeleton[pmax].maxLen; f++ ) {
					delete [] myca->regions[mytid]->myjunctions.fjuncSkeleton[pmax].thefaces[f].voxelizedface;
				}
				myca->regions[mytid]->myjunctions.fjuncSkeleton[pmax].len = 0;
				myca->regions[mytid]->myjunctions.fjuncSkeleton[pmax].maxLen = 0;
				delete [] myca->regions[mytid]->myjunctions.fjuncSkeleton[pmax].thefaces;
				myca->regions[mytid]->myjunctions.fjuncSkeleton[pmax].thefaces = NULL;
			}
		}
	}
}

#endif