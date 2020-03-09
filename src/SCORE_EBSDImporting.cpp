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

void caHdl::ebsd2polycrystal( void )
{
	//determine CA volume geometry and lay EBSD map at x,y, z=0 then "extend" structure along z
	//MK::concept of M. Kuehbach/M. Diehl: extend 2D grains along z as sections of diameter ? Mean DefSizePoissonDiameter
	//how many deformed grains to synthesize in total?
	vector<defgseed>* tmp2defgseeds = NULL;
	tmp2defgseeds = new vector<defgseed>;
	//char scatteringmode = REFERENCE_TO_ROTATED_REFERENCE;

	for ( size_t g = 0; g < myensHdl->expGrains.size(); g++ ) {
		//uint32_t nz = myCAGeometry.nboxedge_td;
		double bunge_ref[3] = { myensHdl->expGrains.at(g).bunge1, myensHdl->expGrains.at(g).bunge2, myensHdl->expGrains.at(g).bunge3 };
		double quaternion_ref[4] = { myensHdl->expGrains.at(g).q0, myensHdl->expGrains.at(g).q1, myensHdl->expGrains.at(g).q2, myensHdl->expGrains.at(g).q3 };

		//baryxy are in EBSD steps, translate in pixel
		//double bx = myensHdl->expGrains.at(g).baryx;
		//double by = myensHdl->expGrains.at(g).baryy;
		double xx = floor(myensHdl->expGrains.at(g).baryx * MICRON2METER / ebsdstepsize); //myCAGeometry.boxedge_rd;
		double yy = floor(myensHdl->expGrains.at(g).baryy * MICRON2METER / ebsdstepsize); //makemyCAGeometry.boxedge_td;
		uint32_t xy[2] = { xx, yy }; //range of xx and yy should be fine...
		QUICKASSERT( xy[0] >= 0 && xy[0] < myCAGeometry.nboxedge_rd );
		QUICKASSERT( xy[1] >= 0 && xy[1] < myCAGeometry.nboxedge_td );

		//grains inherit properties from myensHdl->expGrains.at(g) with scatter
		for ( double zi = 0.0 + (0.5 * myPhysData.defgmean_poisson); zi < myCAGeometry.boxedge_nd; zi += myPhysData.defgmean_poisson ) { //sweep along positive z
			struct cadefg dgr;

			//dgr.worlddefgid is not of relevance for EBSD-based defms synthesis
			//dgr.cellcount already set upon construction of dgr

			//hand over measured orientation into simulation
			//no scatter?
			//dgr.caori = ca_check_disjunctness_core( quaternion_ref, true, bunge_ref[0], bunge_ref[1], bunge_ref[2] );

			//with scatter?
			/*double qderived[4] = {1.0, 0.0, 0.0, 0.0};
			scatter_oriquat( quaternion_ref, DEG2RAD(5.0), qderived, scatteringmode );
			dgr.caori = ca_check_disjunctness_core( qderived, false, 0.0, 0.0, 0.0 ); //internal quat2euler*/
			dgr.caori = ca_check_disjunctness_core( quaternion_ref, true, bunge_ref[0], bunge_ref[1], bunge_ref[2] );

			dgr.rho0 = kam2rho( myensHdl->expGrains.at(g).kam ); //##MK::Kernel Average Misorientation into Dislocation Density
			dgr.rho = dgr.rho0;
			//if (dgr.rho0 >= this->myrhomax) { myrhomax = dgr.rho0; }	//identify highest dislocation density in the SU for customizing an adaptive iteration scheme
			mydefgpool.push_back ( dgr );

			//add grain to maintenance structure which keeps track a reduction rho to speed up the explicit Euler forward integration
			struct dgopt anopt;
			anopt.rho = dgr.rho0;
			anopt.id = mydefgpool.size() - 1;
			anopt.cnt = 0;
			mydgoptimize.push_back( anopt );

			//now as we have a deformed grain we only need to store it for placement
			//where to place the seed for this grain?
			uint32_t z = (uint32_t) ( ((double) (zi/myCAGeometry.boxedge_nd)) * ((double) myCAGeometry.nboxedge_nd) );
			QUICKASSERT( z >= 0 && z < myCAGeometry.nboxedge_nd );
			uint32_t cxyz = (xy[0] + xy[1]*myCAGeometry.nboxedge_rd + z * myCAGeometry.nboxarea_rdtd);

			struct defgseed aseed;
			//aseed.ensdefgpoolid is irrelevant when we construct from EBSD data
			aseed.mydefgpoolid =  (mydefgpool.size() - 1);
			aseed.location = cxyz;

cout << "\t\tAdd defgrain seed kam/rho0/x/y/z/location = " << myensHdl->expGrains.at(g).kam << ";" << dgr.rho0 << ";" << xy[0] << ";" << xy[1] << ";" << z << ";" << aseed.location << endl;

			tmp2defgseeds->push_back( aseed );
		}
	}

	//transfer seeds into tmpdefgseeds
	//##MK::consider to replace tmpdefgseeds in the future

	ndefgseeds = tmp2defgseeds->size();

	QUICKASSERT ( ndefgseeds < CA_ALLOCATION_MAXIMUM );
	tmpdefgseeds = NULL;
	tmpdefgseeds = new struct defgseed[ndefgseeds];
	QUICKASSERT ( tmpdefgseeds != NULL );
	myMemGuard = myMemGuard + (ndefgseeds * sizeof(struct defgseed));
	for ( uint32_t sd = 0; sd < ndefgseeds; sd++ ) { //##MK::was ++sd
		//tmpdefgseeds[sd].ensdefgpoolid = tmp2defgseeds->at(sd).ensdefgpoolid;
		tmpdefgseeds[sd].mydefgpoolid = tmp2defgseeds->at(sd).mydefgpoolid;
		tmpdefgseeds[sd].location = tmp2defgseeds->at(sd).location;
	}

	delete tmp2defgseeds;
}

/*
//debug
//check also against bunge_ref
//double check_b[3] = { myoripool[dgr.caori].bunge1, myoripool[dgr.caori].bunge2, myoripool[dgr.caori].bunge3 }; //which orientation in oripool
//double check_q[4] = { myoripool[dgr.caori].q0, myoripool[dgr.caori].q1, myoripool[dgr.caori].q2, myoripool[dgr.caori].q3 }; //which quaternion in oripool
//double check_bq2b[3] = { -4.0*_PI_,-4.0*_PI_,-4.0*_PI_}; //is this quaternion when back to Euler bringing us to checkb?
//quat2euler( check_q, check_bq2b );
//double check_qinv[4] = { check_q[0], -1.0*check_q[1], -1.0*check_q[2], -1.0*check_q[3] }; //get us the inverse
//double check_binv[3] = {-4.0*_PI_,-4.0*_PI_,-4.0*_PI_};
//quat2euler( check_qinv, check_binv );
*/


void caHdl::ebsd2polycrystal2( void )
{
	//determine CA volume geometry and lay EBSD map at x,y, z=0 then "extend" structure along z
	//MK::concept of M. Kuehbach/M. Diehl: extend 2D mapping along z potentially into sections of diameter ? Mean DefSizePoissonDiameter

	vector<defgseed>* tmp2defgseeds = NULL;
	tmp2defgseeds = new vector<defgseed>;
	//char scatteringmode = REFERENCE_TO_ROTATED_REFERENCE;

	//grains inherit properties from myensHdl->expGrains.at(g) with scatter
	double columnheight = myCAGeometry.boxedge_nd;
	for ( double zi = 0.0; zi < myCAGeometry.boxedge_nd; zi += columnheight ) { //build only one column
		//entire map at different section enables columnar structure of successively larger misorientation to measured orientation
		for ( size_t g = 0; g < myensHdl->expGrains.size(); g++ ) { //expGrains is organized with gid ascending contiguous in [1,N]
			//get measured mean orientation from the TSL/OIM-Mtex data
			double bunge_ref[3] = { myensHdl->expGrains.at(g).bunge1, myensHdl->expGrains.at(g).bunge2, myensHdl->expGrains.at(g).bunge3 };
			double quaternion_ref[4] = { myensHdl->expGrains.at(g).q0, myensHdl->expGrains.at(g).q1, myensHdl->expGrains.at(g).q2, myensHdl->expGrains.at(g).q3 };

			//###############dgr.worlddefgid is not of relevance for EBSD-based defms synthesis
			//################dgr.cellcount already set upon construction of dgr
			struct cadefg dgr;
			//register deformed grain in local list of deformed grains
			/*double qderived[4] = {1.0, 0.0, 0.0, 0.0};
			scatter_oriquat( quaternion_ref, DEG2RAD(5.0), qderived, scatteringmode );
			dgr.caori = ca_check_disjunctness_core( qderived, false, 0.0, 0.0, 0.0 ); //internal quat2euler*/
			dgr.caori = ca_check_disjunctness_core( quaternion_ref, true, bunge_ref[0], bunge_ref[1], bunge_ref[2] );
			dgr.rho0 = kam2rho( myensHdl->expGrains.at(g).kam ); //##MK::Kernel Average Misorientation into Dislocation Density
			dgr.rho = dgr.rho0;
			//if (dgr.rho0 >= this->myrhomax) { myrhomax = dgr.rho0; }	//identify highest dislocation density in the SU for customizing an adaptive iteration scheme
			mydefgpool.push_back ( dgr );

			//add grain to maintenance structure which keeps track of a reduction rho to speed up the explicit Euler forward integration
			struct dgopt anopt;
			anopt.rho = dgr.rho0;
			anopt.id = mydefgpool.size() - 1;
			anopt.cnt = 0;
			mydgoptimize.push_back( anopt );

			//register this deformed grain as a with tmpdefgseeds for the REPLACE_CA_STRUCTURE function to work properly
			struct defgseed aseed;
			aseed.mydefgpoolid =  (mydefgpool.size() - 1);
			//aseed.ensdefgpoolid and aseed.location are irrelevant when we construct from EBSD data
			tmp2defgseeds->push_back( aseed );

			//ID transformation seed2mydefgpool? mydefgpool[tmp2defgseeds[seedid].mydefgpoolid]...
			//ID transformation EBSD pixel expGrain id to seed: 
			//this->myensHdl->expGrains[this->myensHdl->expData[px].gid-1].gid, this->myensHdl->expGrains[0].gid
			//connection to tmp2defgseeds->at()
			//########
		} //all grains from the mapping
	} //next layer of extension along z

	//transfer seeds from temporary into tmpdefgseeds, this workaround is necessary because total number of seeds is unknown at subfunction entry ##MK::consider to replace tmpdefgseeds in the future
	ndefgseeds = tmp2defgseeds->size(); //MK::total number of seeds is an integer multiple of expGrains.size()
	QUICKASSERT ( ndefgseeds < CA_ALLOCATION_MAXIMUM );
	tmpdefgseeds = NULL;
	tmpdefgseeds = new struct defgseed[ndefgseeds];
	QUICKASSERT ( tmpdefgseeds != NULL );
	myMemGuard = myMemGuard + (ndefgseeds * sizeof(struct defgseed));
	for ( uint32_t sd = 0; sd < ndefgseeds; sd++ ) { //##MK::was ++sd
		//tmpdefgseeds[sd].ensdefgpoolid = tmp2defgseeds->at(sd).ensdefgpoolid;
		//tmpdefgseeds[sd].location = tmp2defgseeds->at(sd).location;
		tmpdefgseeds[sd].mydefgpoolid = tmp2defgseeds->at(sd).mydefgpoolid;
	}
	delete tmp2defgseeds; tmp2defgseeds = NULL;

	//set explicitly initial state of every cell in every region
	//##MK::mind dimensions of global CA coordinates whose organization with +x lines along +y stacked slabs along +z is the same as expData, i.e. +x lines then +y lines downwards
	uint32_t column_block_id = 0;
	double height = 0.0;
	uint32_t px = 0;
	uint32_t rxmi, rxmx, rymi, rymx, rzmi, rzmx, extentx, extenty;
	//no domain partitioning into regions along x, but along y and z
	uint32_t npy = this->regions[0]->thePartitioning.nreg_tdy; //MK::fixed array access with 0 works because there is always at least one region and all regions know the global partitioning
	uint32_t npz = this->regions[0]->thePartitioning.nreg_ndz;
	uint32_t ry, rz;
	uint32_t ebsdgid, seedid;
	uint32_t totalexpgrains = myensHdl->expGrains.size();
	//organization of regions is such that first all regions at constant lowest z in direction +y, then next group along +y at const z+1, ...

	for ( uint32_t z = 0; z < myCAGeometry.nboxedge_nd; z++ ) {
		//in which xy slab of regions is this particular CA layer xy, z included?
		for ( uint32_t rr = 0; rr < (1*npy*npz); rr = rr + npy ) { //mind that implicit region id 0 unfolds into ry=0,rz=0, if npy e.g. =3 then next zlayer in 3
			if ( z >= this->regions[rr]->myGeom.nreg_ndmin && z <= this->regions[rr]->myGeom.nreg_ndmax ) { //works because linear cluster of regions along y has same dimensions
				rz = rr / npy; break;
			}
		}

		//read 2d EBSD mapping recurrently when building up columns along z
		px = 0;
		for ( uint32_t y = 0; y < myCAGeometry.nboxedge_td; y++ ) {
			//in which slab yregion is the particular coordinate line y?
			for ( uint32_t rr = 0+(rz*npy); rr < npy+(rz*npy); rr++ ) { //mind that caHdl->regions is organized as contiguous groups of regions along +y stacked in +z
				if ( y >= this->regions[rr]->myGeom.nreg_tdmin && y <= this->regions[rr]->myGeom.nreg_tdmax ) {
					ry = rr-(rz*npy); break;
				}
			}
			//no region partitioning along x so now we know in which particular region to place the xline assignments
			uint32_t here = ry+(rz*npy); //MK::mentality --- a pixel from the EBSD mapping becomes mapped on global CA volume that is partitioned into regions
			uint32_t* thecells = this->regions[here]->mycellgrid;

			//x,y,z are SU global coordinates, we have to translate into local coordinates...
			rxmi = this->regions[here]->myGeom.nreg_rdmin;
			rxmx = this->regions[here]->myGeom.nreg_rdmax;
			rymi = this->regions[here]->myGeom.nreg_tdmin;
			rymx = this->regions[here]->myGeom.nreg_tdmax; //local region cell coordinate inclusive limits
			rzmi = this->regions[here]->myGeom.nreg_ndmin;
			rzmx = this->regions[here]->myGeom.nreg_ndmax;
			extentx = rxmx - rxmi + 1;
			QUICKASSERT( extentx == myCAGeometry.nboxedge_rd );
			extenty = rymx - rymi + 1;

			uint32_t yzoffset = (y-rymi)*extentx + (z-rzmi)*extentx*extenty;
			for ( uint32_t x = 0; x < extentx; x++ ) { //setup region partitioning more efficiently
				ebsdgid = myensHdl->expData[px].gid; //the grain at this position, expGrains is contiguous with IDs on [1,N]

				seedid = (ebsdgid - 1) + (column_block_id * totalexpgrains); //MK::-1 substracting on ebsdgid necessary because expGrains.gid start at 1 but seedid are positions on linear array!

				//make the CA cell to represent a voxel of a particular seeded deformed grain
				thecells[x+yzoffset] = seedid;
				//translate next pixel to CA volume
				px++;
			} //xline completed
		} //ylines done

		height += this->myCAGeometry.cellsize;
		cout << "ZLayer " << column_block_id << " up to height " << height << " build, assigned a total of " << px << " cells." << endl;
		if ( height >= columnheight ) { //next block
			cout << "Beginning new block layer..." << endl;
			height = 0.0;
			column_block_id++;
		}
	} //all zlayers were build
}



#endif