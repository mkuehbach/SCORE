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


uint32_t caHdl::solve_nucmodeling_countgbcells( vector<organizeGBNuc>* luu )
{
	//boundary bucket containers have exact counts and array lengths!
	uint32_t nthreads = this->regions.size();
	uint32_t tid = MASTER;
	uint32_t nc = 0;
	for ( uint32_t pmax = 0; pmax < this->regions[MASTER]->myjunctions.fjuncSkeletonSize; pmax++ ) {
		tid = workPartitioning(pmax, nthreads);
		for ( uint32_t f = 0; f < this->regions[tid]->myjunctions.fjuncSkeleton[pmax].len; f++ ) {
			QUICKASSERT(  this->regions[tid]->myjunctions.fjuncSkeleton[pmax].thefaces[f].nvoxeltotal > 0 );

			struct organizeGBNuc cnt;
			cnt.pmax = pmax;
			cnt.thr = tid;
			cnt.fid = f;
			cnt.istart = nc;
			cnt.iend = nc + this->regions[tid]->myjunctions.fjuncSkeleton[pmax].thefaces[f].nvoxeltotal - 1; //inclusive bounds
			nc = nc + this->regions[tid]->myjunctions.fjuncSkeleton[pmax].thefaces[f].nvoxeltotal;
			luu->push_back( cnt );
		}
	}
	return nc;
}


struct organizeGBNuc caHdl::find_gbface( uint32_t idx, vector<organizeGBNuc>* luu )
{
	struct organizeGBNuc foo;
	//##MK::next line should not be necessary because of default constructor for the organizeGBNuc object...
	foo.pmax = UINT32T_MAX;	foo.thr = UINT32T_MAX;	foo.fid = UINT32T_MAX;	foo.istart = UINT32T_MAX;	foo.iend = UINT32T_MAX;

	for ( unsigned int i = 0; i < luu->size(); i++ ) {
		if ( idx < luu->at(i).istart ) continue;
		if ( idx > luu->at(i).iend ) continue;

		//foo = luu->at(i);
		foo.pmax = luu->at(i).pmax;
		foo.thr = luu->at(i).thr;
		foo.fid = luu->at(i).fid;
		foo.istart = luu->at(i).istart;
		foo.iend = luu->at(i).iend;
		return foo;
	}

	//nothing found, return foo with error flags to check against
	return foo;
}


void caHdl::solve_nucmodeling_gbnuc_pickrandomly( void )
{
	//this grain boundary nucleation model spreads a total of NucleiDesired randomly on the available boundary sites
	uint32_t NucleiDesired = this->myNucleationModel.defaultnucdensity;
	uint32_t NucleiPlaced = 0;
	uint32_t NucleiUnplaceable = 0;

	//MK::it is errorneous to simply pick a grain and a boundary at random and place a nucleus at a cell because 
	//the individual boundary face cell counts are different, therefore first we have to find the size of the boundaries

	//count boundary cell structure
	vector<organizeGBNuc>* lu = NULL;
	lu = new vector<organizeGBNuc>;
	QUICKASSERT( lu != NULL );

	uint32_t nBndCellsTotal = solve_nucmodeling_countgbcells( lu );
	double nBndCellsDbl = (double) nBndCellsTotal;
	bool found = false;

	for ( uint32_t nuc = 0; nuc < NucleiDesired; nuc++ ) { //pick cells at random
		found = false;
		uint32_t rndcell = nBndCellsTotal;
		uint32_t lc = 0;
		uint32_t repeat = 0;

		while ( found == false && repeat < REPEAT_MAX_DRAWING ) {
			rndcell = (localprng.MersenneTwister() * nBndCellsDbl);
			if ( rndcell >= nBndCellsTotal ) {
				repeat++;
				continue;
			}

			//we have a site, but where is it?
			organizeGBNuc whereisthis;
			whereisthis = find_gbface( rndcell, lu );
			QUICKASSERT( whereisthis.pmax != UINT32T_MAX );

			//MK::we remember that all cellsBndFast arrays maintain not directly the location but location + 1 to enable
			//the flagging of cells as already occupied, so lets check for this
			lc = rndcell - whereisthis.istart; 
			if ( this->regions[whereisthis.thr]->myjunctions.fjuncSkeleton[whereisthis.pmax].thefaces[whereisthis.fid].voxelizedface[lc].location == EMPTY ) {
				repeat++;
				continue;
			}

			//there is a site and it is not occupied, then claim it
			uint32_t nucsite = this->regions[whereisthis.thr]->myjunctions.fjuncSkeleton[whereisthis.pmax].thefaces[whereisthis.fid].voxelizedface[lc].location - MARKER_TO_IDENTIFY_EXPELLED_SITES;
			QUICKASSERT( nucsite >= 0 && nucsite < this->myCAGeometry.nboxvol_rdtdnd );

			//mark site as occupied
			this->regions[whereisthis.thr]->myjunctions.fjuncSkeleton[whereisthis.pmax].thefaces[whereisthis.fid].voxelizedface[lc].location = EMPTY;

			//pick a random grain from the ensemble input data
			uint32_t cand = ( (double) myensHdl->worldrxgpool.size() * localprng.MersenneTwister() );

			//get orientation from this candidate
			double nucbunge[3];
			uint32_t ensoriid = myensHdl->worldrxgpool[cand].ori;
			nucbunge[0] = myensHdl->worldoripool[ensoriid].bunge1;
			nucbunge[1] = myensHdl->worldoripool[ensoriid].bunge2;
			nucbunge[2] = myensHdl->worldoripool[ensoriid].bunge3;

			double nucq[4];
			nucq[0] = myensHdl->worldoripool[ensoriid].q0;
			nucq[1] = myensHdl->worldoripool[ensoriid].q1;
			nucq[2] = myensHdl->worldoripool[ensoriid].q2;
			nucq[3] = myensHdl->worldoripool[ensoriid].q3;

			//pick a random orientation from the available ones with properties inheriting from the two neighboring pairs

			//register orientation locally the orientation
			uint32_t noid = ca_check_disjunctness_core( nucq, true, nucbunge[0], nucbunge[1], nucbunge[2] );

			//take a random one MK::mind registering is different!

			//rnd_oriquat_shoemake( nucq );
			//uint32_t noid = ca_check_disjunctness_core( nucq, false, 0.0, 0.0, 0.0 );

			//##MK::or implement other models or more complex models or twinning models if there are better experimental data justifying there usage

			//register nucleus in the SU
			struct carxg newNucleus;

			newNucleus.caori = noid;
			newNucleus.cellcount = 0;
			newNucleus.nucsite = nucsite;
			newNucleus.startingsite = nucsite;
			newNucleus.tincub = 0.0; //##MK::incubation time model is site-saturated at the moment

			myrxgpool.push_back ( newNucleus );
			NucleiPlaced++;
			found = true;
#ifdef REPORTSTYLE_DEVELOPER
			std::cout << "myRank;jobid;i;caori;cellcnt;nucsite;startingsite;tincub\t\t" << this->myensHdl->myRank << ";" << this->jobid << ";" << nuc << ";" << myrxgpool[myrxgpool.size()-1].caori << ";" << myrxgpool[myrxgpool.size()-1].cellcount << ";" << myrxgpool[myrxgpool.size()-1].nucsite << ";" << myrxgpool[myrxgpool.size()-1].startingsite << ";" << myrxgpool[myrxgpool.size()-1].tincub << std::endl;
#endif
		}
		if ( repeat >= REPEAT_MAX_DRAWING ) {
			std::cout << "WARNING::GBNucleus " << nuc << " not placeable!" << std::endl;
			NucleiUnplaceable++;
		}
	}

	delete lu; lu = NULL;
	std::cout << this->myensHdl->myRank << " GBNuclei Desired/Placed/Unplaceable " << NucleiDesired << ";" << NucleiPlaced << ";" << NucleiUnplaceable << std::endl;
}


void caHdl::solve_nucmodeling_csr_enforced( void )
{
	//MK::place exactly myNucleationModel.defaultnucdensity in each solitary unit domain (THIS IS THE OLD SCORE NUCLEATION MODEL utilzed in the K\"uhbach et. al. 2015 paper)
	uint32_t nuc = 0;
	uint32_t maxattempts = 0;
	uint32_t ncellstotal = myCAGeometry.nboxvol_rdtdnd;
	ensembleHdlP ens = this->myensHdl;
	uint32_t nworldrxgpool = ens->worldrxgpool.size();

	do {
		//pick a side at random
		uint32_t luckySite = localprng.MersenneTwister() * ncellstotal;

		//assume it to be still unoccupied
		bool AlreadyOccupied = false;
		for ( uint32_t n = 0; n < nuc; n++ ) {
			if ( luckySite == myrxgpool[n].nucsite ) {
				AlreadyOccupied = true; //MK::in older version a mistake was made namely AlreadyOccupied == true !
				break;
			}
		}

		if ( AlreadyOccupied == false ) { //obviously juvenile place
			//nucleus orientation, ##MK::pick nucleus at random from worldrxgpool
			uint32_t luckyNucleus = localprng.MersenneTwister() * nworldrxgpool;

			//contrary to the deformed grains in mydefgpool, each nucleus requires entry in myrxgpool as the incubation time can be different and furthermore also the site and single-nucleus resolved kinetics should be stored
			struct carxg crxgr;

/*
			//get orientation from predefined list
			uint32_t ensoriid = ens->worldrxgpool[luckyNucleus].ori;

			double ensbunge[3];
			ensbunge[0] = ens->worldoripool[ensoriid].bunge1;
			ensbunge[1] = ens->worldoripool[ensoriid].bunge2;
			ensbunge[2] = ens->worldoripool[ensoriid].bunge3;
			double ensq[4];
			ensq[0] = ens->worldoripool[ensoriid].q0;
			ensq[1] = ens->worldoripool[ensoriid].q1;
			ensq[2] = ens->worldoripool[ensoriid].q2;
			ensq[3] = ens->worldoripool[ensoriid].q3;

			crxgr.caori = ca_check_disjunctness_core( ensq, true, ensbunge[0], ensbunge[1], ensbunge[2] );
*/

			//get a random rotation from SO3 and interpret as orientation, mind different registration of Euler angles!
			double ensq[4] = {1.0, 0.0, 0.0, 0.0};
			localmath.rnd_quat_shoemake( ensq );
			crxgr.caori = ca_check_disjunctness_core( ensq, false, 0.0, 0.0, 0.0 );

			crxgr.cellcount = 0;
			crxgr.nucsite = luckySite; //common locations for debug purposes 4060300; //32240600 (401^3 center); //4060300 (201^3 center); //luckySite; //(200^3 center = 4020100); //(100^3 center=505050); //mluckySite;
			crxgr.startingsite = luckySite;
			crxgr.tincub = 0; //ens->worldrxgpool[ensoriid].tincub;

			myrxgpool.push_back ( crxgr );

#ifdef REPORTSTYLE_DEVELOPER
			cout << "jobid;nuc;att;nucsite;luckyNucleusWorldRXGPool\t\t" << this->jobid << ";" << nuc << ";" << maxattempts << ";" << myrxgpool[myrxgpool.size()-1].startingsite << ";" << luckyNucleus << endl;
#endif
			nuc++;
			continue;
		}

		//damn, collision - AlreadyOccupied == true so attempt a new place
		maxattempts++;
	} 
	while (nuc < myNucleationModel.defaultnucdensity && maxattempts < MAXATTEMPTS_NUCSITE );

	//add characterization of resulting nucleus placement here...

#ifdef DETAILED_PROMPTS
	if ( myensHdl->myRank == MASTER ) { cout << this->jobid << "\t\tNucleation completed." << endl; }
#endif
}


void caHdl::solve_nucmodeling_csr_pickrandomly( void )
{
	//MK::in a large CSR point pattern with a unit density what is the size of the observation window in generalized coordinates?
	//with myNucleationModel.defaultnucdensity as the limiting unit density for very many draws of owins
	//MK::the difference of this model to CSR_ENFORCE is that when the observation volume, i.e. the SU is small defaultnucdensity is only the expectation value
	//of the entire ensemble of SU. The exact number of nuclei in a specific domain then follows a distribution with scatter
	double scale = pow( (((double) myNucleationModel.defaultnucdensity) / ((double) MASTER_CSR_NPOINTS)), (1.0/3.0) );

	double owin[3];
	owin[0] = scale;	if ( owin[0] > 1.0 ) { cout << "ERROR::MASTER point process is too small!" << endl; return; }
	owin[1] = scale;	if ( owin[1] > 1.0 ) { cout << "ERROR::MASTER point process is too small!" << endl; return; }
	owin[2] = scale;	if ( owin[2] > 1.0 ) { cout << "ERROR::MASTER point process is too small!" << endl; return; }

	//populate a CSR point pattern on a unit cube from which to take a small sample
	pointP pp3 = NULL;
	pp3 = new point[MASTER_CSR_NPOINTS];
	QUICKASSERT ( pp3 != NULL );

	for ( uint32_t p = 0; p < MASTER_CSR_NPOINTS; ++p ) {
		pp3[p].x = localprng.MersenneTwister(); //[0,1]^3
		pp3[p].y = localprng.MersenneTwister();
		pp3[p].z = localprng.MersenneTwister();
	}

	//sample an owin from pp3 on the unit cube
	double owinoo[6];
	bool inside = false;

	while (inside == false ) {
		owinoo[XMI] = localprng.MersenneTwister(); //lower bound but guaranteed positive
		owinoo[XMX] = owinoo[XMI] + owin[0];
		owinoo[YMI] = localprng.MersenneTwister();
		owinoo[YMX] = owinoo[YMI] + owin[1];
		owinoo[ZMI] = localprng.MersenneTwister();
		owinoo[ZMX] = owinoo[ZMI] + owin[2];

		if ( owinoo[XMX] < 1.0 && owinoo[YMX] < 1.0 && owinoo[ZMX] < 1.0) //all coordinates positive owinoo cannot protrude in negative direction
			inside = true;
	}

	//extract all points from this owinoo window and sample via (i - owinoo[iMI])/owin[i] * ni into an automaton coordinate
	double dx, dy, dz;
	uint32_t ix, iy, iz, luckysite;
	uint32_t nbx = myCAGeometry.nboxedge_rd; //avoid the cache collision when calling this in the for loop to select the points
	uint32_t nby = myCAGeometry.nboxedge_td;
	uint32_t nbz = myCAGeometry.nboxedge_nd;
	bool valid = false;

	for ( uint32_t p = 0; p < MASTER_CSR_NPOINTS; ++p ) {
		if ( pp3[p].x < owinoo[XMI] ) continue; //most likely case
		if ( pp3[p].x > owinoo[XMX] ) continue;
		if ( pp3[p].y < owinoo[YMI] ) continue;
		if ( pp3[p].y > owinoo[YMX] ) continue;
		if ( pp3[p].z < owinoo[ZMI] ) continue;
		if ( pp3[p].z > owinoo[ZMX] ) continue;

		//p is inside owinoo so, discretize location, dx, dy, dz > 0 < 1
		dx = (pp3[p].x - owinoo[XMI]) / (owinoo[XMX] - owinoo[XMI]);
		dy = (pp3[p].y - owinoo[YMI]) / (owinoo[YMX] - owinoo[YMI]);
		dz = (pp3[p].z - owinoo[ZMI]) / (owinoo[ZMX] - owinoo[ZMI]);
		ix = dx * nbx;
		iy = dy * nby;
		iz = dz * nbz;
		luckysite = ix + (iy * nbx) + (iz * nbx * nby);

//cout << "csrpickrnd;luckysite = " << luckysite << "\t\t" << ix << ";" << iy << ";" << iz << endl;

		//disprove the assumption that a nucleus was already placed within the grid resolution cellsize, if so discard
		valid = true;
		for ( uint32_t s = 0; s < myrxgpool.size(); s++ ) {
			if ( luckysite == myrxgpool[s].nucsite ) { valid = false; break; }
		}

		//still true, okay nothing placed at the site yet so plant a nucleus
		if ( valid == true ) {

			//assign orientation and a nucleus
			uint32_t luckyNucleus = localprng.MersenneTwister() * this->myensHdl->worldrxgpool.size();

//cout << "p/luckyNucleus/worldrxgpool = " << p << "\t\t" << luckyNucleus << "\t\t" << this->myensHdl->worldrxgpool.size() << endl;

			//contrary to the deformed grains in mydefgpool, each nucleus requires entry in myrxgpool as then incubation time can be different and furthermore also the site and single-nucleus resolved kinetics should be stored
			struct carxg crxgr;

			//get orientation
			uint32_t ensoriid = myensHdl->worldrxgpool[luckyNucleus].ori;

			double ensbunge[3];
			ensbunge[0] = myensHdl->worldoripool[ensoriid].bunge1;
			ensbunge[1] = myensHdl->worldoripool[ensoriid].bunge2;
			ensbunge[2] = myensHdl->worldoripool[ensoriid].bunge3;
			double ensq[4];
			ensq[0] = myensHdl->worldoripool[ensoriid].q0;
			ensq[1] = myensHdl->worldoripool[ensoriid].q1;
			ensq[2] = myensHdl->worldoripool[ensoriid].q2;
			ensq[3] = myensHdl->worldoripool[ensoriid].q3;

			//orientation recategorization utilizes that not all possible orientations are used in the local CA grid thus a more efficient hashtable can be generated however at additional memory costs...
			crxgr.caori = ca_check_disjunctness_core( ensq, true, ensbunge[0], ensbunge[1], ensbunge[2] );
			crxgr.cellcount = 0;
			crxgr.nucsite = luckysite;
			crxgr.startingsite = luckysite;
			crxgr.tincub = myensHdl->worldrxgpool[ensoriid].tincub;

			myrxgpool.push_back ( crxgr );
		} //next point
	}

	delete [] pp3;
	pp3 = NULL;
}


void caHdl::solve_nucmodeling_csrdiehl( void )
{
	if ( myensHdl->experimentInput == false ) { cout << "ERROR::This nucleation model is design only for utilization with a SEM/EBSD mapping!" << endl; return; }
	//MK::random subwindow of a CSR point process (without enforcing fixed number of nuclei), orientation inheriting from parent grains + scatter
	//MK::for simulation with myensHdl->experimentInput == true no recrystallized grains were defined a priori hence we have to
	//MK::add nuclei into myrxgpool and add their orientation as well-
	double scale = pow( (((double) myNucleationModel.defaultnucdensity) / ((double) MASTER_CSR_NPOINTS)), (1.0/3.0) );

	double owin[3];
	owin[0] = scale;	if ( owin[0] > 1.0 ) { cout << "ERROR::MASTER point process is too small!" << endl; return; }
	owin[1] = scale;	if ( owin[1] > 1.0 ) { cout << "ERROR::MASTER point process is too small!" << endl; return; }
	owin[2] = scale;	if ( owin[2] > 1.0 ) { cout << "ERROR::MASTER point process is too small!" << endl; return; }

	//populate a CSR point pattern on a unit cube from which to take a small sample
	pointP pp3 = NULL;
	pp3 = new point[MASTER_CSR_NPOINTS];
	QUICKASSERT ( pp3 != NULL );

	for ( uint32_t p = 0; p < MASTER_CSR_NPOINTS; ++p ) {
		pp3[p].x = localprng.MersenneTwister(); //[0,1]^3
		pp3[p].y = localprng.MersenneTwister();
		pp3[p].z = localprng.MersenneTwister();
	}

	//sample an owin from pp3 on the unit cube
	double owinoo[6];
	bool inside = false;

	while (inside == false ) {
		owinoo[XMI] = localprng.MersenneTwister(); //lower bound but guaranteed positive
		owinoo[XMX] = owinoo[XMI] + owin[0];
		owinoo[YMI] = localprng.MersenneTwister();
		owinoo[YMX] = owinoo[YMI] + owin[1];
		owinoo[ZMI] = localprng.MersenneTwister();
		owinoo[ZMX] = owinoo[ZMI] + owin[2];

		if ( owinoo[XMX] < 1.0 && owinoo[YMX] < 1.0 && owinoo[ZMX] < 1.0) //all coordinates positive owinoo cannot protrude in negative direction
			inside = true;
	}

	//extract all points from this owinoo window and sample via (i - owinoo[iMI])/owin[i] * ni into an automaton coordinate
	double dx, dy, dz;
	uint32_t ix, iy, iz, luckysite;
	uint32_t nbx = myCAGeometry.nboxedge_rd; //avoid the cache collision when calling this in the for loop to select the points
	uint32_t nby = myCAGeometry.nboxedge_td;
	uint32_t nbz = myCAGeometry.nboxedge_nd;
	bool valid = false;
	char scatteringmode = REFERENCE_TO_ROTATED_REFERENCE;

	//place all nuclei inside owinoo into the simulation domain
	for ( uint32_t p = 0; p < MASTER_CSR_NPOINTS; ++p ) {
		if ( pp3[p].x < owinoo[XMI] ) continue; //most likely case for small volume of owin in comparison to unit cube
		if ( pp3[p].x > owinoo[XMX] ) continue;
		if ( pp3[p].y < owinoo[YMI] ) continue;
		if ( pp3[p].y > owinoo[YMX] ) continue;
		if ( pp3[p].z < owinoo[ZMI] ) continue;
		if ( pp3[p].z > owinoo[ZMX] ) continue;

		//p is inside owinoo so, discretize location, dx, dy, dz > 0 < 1
		dx = (pp3[p].x - owinoo[XMI]) / (owinoo[XMX] - owinoo[XMI]);
		dy = (pp3[p].y - owinoo[YMI]) / (owinoo[YMX] - owinoo[YMI]);
		dz = (pp3[p].z - owinoo[ZMI]) / (owinoo[ZMX] - owinoo[ZMI]);
		ix = dx * nbx; //global coordinates
		iy = dy * nby;
		iz = dz * nbz;
		luckysite = ix + (iy * nbx) + (iz * nbx * nby);

//cout << "csrpickrnd;luckysite = " << luckysite << "\t\t" << ix << ";" << iy << ";" << iz << endl;

		//disprove the assumption that a nucleus was already placed within the grid resolution cellsize, if so discard
		valid = true;
		for ( uint32_t s = 0; s < myrxgpool.size(); s++ ) {
			if ( luckysite == myrxgpool[s].nucsite ) { valid = false; break; }
		}

		//still true, okay nothing placed at the site yet so plant a nucleus
		if ( valid == true ) {
			double qderived[4] = {1.0, 0.0, 0.0, 0.0 };

			//define the orientation of the nucleus
			if ( myNucleationModel.csrnucleation == CSRNUCLEATION_DIEHL_RANDOMSO3 ) {
				//##a random orientation from the SO3
				localmath.rnd_quat_shoemake( qderived );
				cout << "RANDOMSO3 q0/q1/q2/q3 " << qderived[0] << ";" << qderived[1] << ";" << qderived[2] << ";" << qderived[3] << "\n";
			}
			else if (  myNucleationModel.csrnucleation == CSRNUCLEATION_DIEHL_RANDOMDEFORI ) {
				//##the orientation from a random of the orientations of the deformed grains
				double pickdefgr = floor(localprng.MersenneTwister() * static_cast<double>(ndefgseeds));
				uint32_t thisone = (static_cast<uint32_t>(pickdefgr) < ndefgseeds) ? static_cast<uint32_t>(pickdefgr) : ndefgseeds-1;
				uint32_t thisid = tmpdefgseeds[thisone].mydefgpoolid;
				uint32_t oid = mydefgpool[thisid].caori;
				qderived[0] = myoripool[oid].q0;
				qderived[1] = myoripool[oid].q1;
				qderived[2] = myoripool[oid].q2;
				qderived[3] = myoripool[oid].q3;
				//##MK::no scatter is added
				cout << "RANDOMDEFORI q0/q1/q2/q3 " << qderived[0] << ";" << qderived[1] << ";" << qderived[2] << ";" << qderived[3] << "\n";
			}
			else { // myNucleationModel.csrnucleation == CSRNUCLEATION_DIEHL_SCATTEREXISTENT
				//inspect parent grain at the site, as nucleation functions are executed before solve_REPLACE_CA_STRUCTURE() we have to
				//access deformed grains via the tmpdefgseeds ID references
				uint32_t parentid = tmpdefgseeds[getDeformedGrainID( ix, iy, iz )].mydefgpoolid;
				uint32_t oid = mydefgpool[parentid].caori;
				double qparent_ref[4] = { myoripool[oid].q0, myoripool[oid].q1, myoripool[oid].q2, myoripool[oid].q3 };

				//nuclei with scatter to parent
				localmath.scatter_oriquat( qparent_ref, DEG2RAD(5.0), qderived, scatteringmode );

				cout << "SCATTEREXISTENT q0/q1/q2/q3 " << qderived[0] << ";" << qderived[1] << ";" << qderived[2] << ";" << qderived[3] << " disori angle " << RAD2DEG(localmath.disori_angle_oriquat_cubic(qparent_ref, qderived)) << "\n";
			}

			//contrary to the deformed grains in mydefgpool, each nucleus requires entry in myrxgpool as then incubation time can be different and furthermore also the site and single-nucleus resolved kinetics should be stored
			struct carxg rxg;
			rxg.caori = ca_check_disjunctness_core( qderived, false, 0.0, 0.0, 0.0 ); //internal quat2euler
			rxg.cellcount = 0;
			rxg.nucsite = luckysite;
			rxg.startingsite = luckysite;
			rxg.tincub = 0.0; //##assume site-saturation

			myrxgpool.push_back ( rxg );

			cout << "NucleationDiehl placing id/ix/iy/iz/bunge1/bunge2/bunge3 = " << (myrxgpool.size() - 1) << "\t\t" << ix << ";" << iy << ";" << iz << "\tBunge123\t" << myoripool[rxg.caori].bunge1 << ";" << myoripool[rxg.caori].bunge2 << ";" << myoripool[rxg.caori].bunge3 << endl;
		} //next point
	}

	delete [] pp3;
	pp3 = NULL;
}


void caHdl::init_myMaternUniverse( void )
{
	//a very large Matern cluster process on a cubic RVE of size rvesize in micron is squeezed on a unit cube in \mathbb{R}^3
	//point process is on cubic RVE with myNucleationModel.cluster_rvesize^3 micron each cluster radius is cluster_extenda * cellsize
	double scalingradius = myNucleationModel.cluster_a * myCAGeometry.cellsize / myNucleationModel.cluster_rvesize;
	QUICKASSERT ( scalingradius < DEFAULT_MASTERMATERN_MAXRAD );

	mypointprocess.reserve ( (long) (myNucleationModel.cluster_nclust * myNucleationModel.cluster_lambda) );

cout << this->jobid << " -- GenerateMaternClusterUniverse NClust/MeanDensityPerClust/ScalingRVESize/RadiusClust/ExtendABC = " << myNucleationModel.cluster_nclust << ";" << myNucleationModel.cluster_lambda << ";" << myNucleationModel.cluster_rvesize << ";" << scalingradius << ";" << myNucleationModel.cluster_a << ";" << myNucleationModel.cluster_b << ";" << myNucleationModel.cluster_c << endl;

	double nx, ny, nz, len;
	//##MK::distortion currently not implemented
	//##MK::assumption is that number of clusters is large (at least 1000!) to sample random 3D poisson point process

	double cx, cy, cz;
	long npc;
	for (unsigned long c = 0; c < myNucleationModel.cluster_nclust; c++) {
		//construction of Matern: define first a cluster center that itself does not count as a point
		cx = localprng.MersenneTwister();
		cy = localprng.MersenneTwister();
		cz = localprng.MersenneTwister();

		//how many points to locate in one cluster, evaluate poisson distribution, naive cumulant method
		long k = 0;
		double lambda = myNucleationModel.cluster_lambda;
		double unifrnd = localprng.MersenneTwister();
		double P = exp(-lambda);
		double sum = P;

		npc = 0; //handles case sum >= unifrnd
		if ( sum < unifrnd ) {
			for ( k = 1; k < DEFAULT_KMAX_POISSRND; k++) {
				P *= lambda / (double) k;
				sum += P;
				if ( sum >= unifrnd ) break;
			}
		}
		npc = k;


		//populate the process
		for (unsigned long p = 0; p < npc; p++) {
			//random rotor instead of predictor/branch/reject
			nx = 2*localprng.MersenneTwister() - 1.0; //(-1.0,1.0) a random 3D vector pointing into space...
			ny = 2*localprng.MersenneTwister() - 1.0;
			nz = 2*localprng.MersenneTwister() - 1.0;
			len = ( 1.0 / pow( ( SQR(nx)+SQR(ny)+SQR(nz) ) , 0.5) ) * (localprng.MersenneTwister() * scalingradius); //normalized component value * scale(1.0 / ... ) normalizes define scaling of random rotor inside a sphere of maximum radius mmatr
			//scale direction vector
			nx *= len;
			ny *= len;
			nz *= len;
			//mount scaled vector pointing in nxyz direction vector into the cluster center
			nx += cx;
			ny += cy;
			nz += cz;

			//apply PBC in the universe
			if ( nx < 0.0 ) nx += 1.0;
			if ( nx > 1.0 ) nx -= 1.0;
			if ( ny < 0.0 ) ny += 1.0;
			if ( ny > 1.0 ) ny -= 1.0;
			if ( nz < 0.0 ) nz += 1.0;
			if ( nz > 1.0 ) nz -= 1.0;

			struct point ap;
			ap.x = nx;
			ap.y = ny;
			ap.z = nz;

//cout << "\t\t" << c << "\t\t" << npc << "\t\t" << p << "\t\t" << lambda << ";" << unifrnd << ";" << P << ";" << sum << ";" << cx << ";" << cy << ";" << cz << "--" << len << ";" << nx << ";" << ny << ";" << nz << endl;

			mypointprocess.push_back(ap);
		} //populate the cluster
	} //for all clusters
}


void caHdl::pickrandomly_myMaternUniverse( void )
{
	//the solitary unit domain has a specific size (edge length in microns) that is a fraction of myNucleationModel.clust_rvesize
	double owinx = myCAGeometry.boxedge_rd / myNucleationModel.cluster_rvesize;
	double owiny = myCAGeometry.boxedge_td / myNucleationModel.cluster_rvesize;
	double owinz = myCAGeometry.boxedge_nd / myNucleationModel.cluster_rvesize;

	double owinoo[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	bool found = false;

	//MK::to center the domain on a particular point view PPANALYSIS source code...
	//cut myself a window from mypointprocess
	found = false;
	while ( found == false ) {
		owinoo[XMI] = localprng.MersenneTwister();
		owinoo[XMX] = owinoo[XMI] + owinx;
		owinoo[YMI] = localprng.MersenneTwister();
		owinoo[YMX] = owinoo[YMI] + owiny;
		owinoo[ZMI] = localprng.MersenneTwister();
		owinoo[ZMX] = owinoo[ZMI] + owinz;

		//inside?
		if (owinoo[XMI] < 0.0) 	continue;
		if (owinoo[XMX] > 1.0)	continue;
		if (owinoo[YMI] < 0.0) 	continue;
		if (owinoo[YMX] > 1.0)	continue;
		if (owinoo[ZMI] < 0.0) 	continue;
		if (owinoo[ZMX] > 1.0)	continue;

		found = true;
	}

//cout << "Owinxyz = " << owinx << ";" << owiny << ";" << owinz << endl;
cout << "The window = " << owinoo[XMI] << ";" << owinoo[XMX] << ";" << owinoo[YMI] << ";" << owinoo[YMX] << ";" << owinoo[ZMI] << ";" << owinoo[ZMX] << endl;

	QUICKASSERT ( found == true );

	//now test which points from mypointprocess are in this window and plan them for the simulation
	double afftransx = (1.0 / (owinoo[XMX] - owinoo[XMI])) * myCAGeometry.boxedge_rd;
	double afftransy = (1.0 / (owinoo[YMX] - owinoo[YMI])) * myCAGeometry.boxedge_td;
	double afftransz = (1.0 / (owinoo[ZMX] - owinoo[ZMI])) * myCAGeometry.boxedge_nd;

	struct point sp;
	uint32_t nrdtd = myCAGeometry.nboxarea_rdtd;
	uint32_t nrd = myCAGeometry.nboxedge_rd;
	//uint32_t nnd = myCAGeometry.nboxedge_nd;
	//uint32_t ntd = myCAGeometry.nboxedge_td;
	uint32_t ix, iy, iz;
	//double box_rd = myCAGeometry.boxedge_rd;
	//double box_nd = myCAGeometry.boxedge_nd;
	//double box_td = myCAGeometry.boxedge_td;
	double cellsize = myCAGeometry.cellsize;

	uint32_t nworldrxgpool = this->myensHdl->worldrxgpool.size();

	long n = mypointprocess.size();
	bool occupied = false;

	for (unsigned long p = 0; p < n; p++) { //all nuclei
		if (mypointprocess[p].x < owinoo[XMI]) continue;
		if (mypointprocess[p].x > owinoo[XMX]) continue;
		if (mypointprocess[p].y < owinoo[YMI]) continue;
		if (mypointprocess[p].y > owinoo[YMX]) continue;
		if (mypointprocess[p].z < owinoo[ZMI]) continue;
		if (mypointprocess[p].z > owinoo[ZMX]) continue;

		//location inside the CA domain MK:: mypointprocess[p].i >= owinoo[iMX]
		sp.x = (mypointprocess[p].x - owinoo[XMI]) * afftransx;
		sp.y = (mypointprocess[p].y - owinoo[YMI]) * afftransy;
		sp.z = (mypointprocess[p].z - owinoo[ZMI]) * afftransz;

		//map on discrete CA grid
		ix = sp.x / cellsize;
		iy = sp.y / cellsize;
		iz = sp.z / cellsize;

		uint32_t luckySite = ix + (iy * nrd) + (iz * nrdtd);

#ifdef REPORTSTYLE_DEVELOPER
	cout << "\t\tJobID;p;luckySite;ix;iy;iz\t\t" << this->jobid << "\t\t" << p << "\t\t" << luckySite << ";" << ix << ";" << iy << ";" << iz << endl; // << "--" << sp.x << ";" << sp.y << ";" << sp.z <<  endl;
#endif

		//site still unoccupied?
		occupied = false;
		//uint32_t nucSite = INVALID_NUCSITE;

		//disprove expectation that the site is still free...
		for ( uint32_t i = 0; i < myrxgpool.size(); i++ ) {
			if ( luckySite == myrxgpool[i].nucsite ) {
				occupied = true;
				break;
			}
		}

		//placement only if site still unoccupied
		if ( occupied == false ) {
			//###currently, pick nucleus at random from worldrxgpool 
			uint32_t luckyNucleus = localprng.MersenneTwister() * nworldrxgpool;

			//contrary to the deformed grains in mydefgpool, each nucleus requires entry in myrxgpool as then incubation time can be different and furthermore also the site and single-nucleus resolved kinetics should be stored
			struct carxg crxgr;

			//get orientation
			uint32_t ensoriid = myensHdl->worldrxgpool[luckyNucleus].ori;

			double ensbunge[3];
			ensbunge[0] = myensHdl->worldoripool[ensoriid].bunge1;
			ensbunge[1] = myensHdl->worldoripool[ensoriid].bunge2;
			ensbunge[2] = myensHdl->worldoripool[ensoriid].bunge3;
			double ensq[4];
			ensq[0] = myensHdl->worldoripool[ensoriid].q0;
			ensq[1] = myensHdl->worldoripool[ensoriid].q1;
			ensq[2] = myensHdl->worldoripool[ensoriid].q2;
			ensq[3] = myensHdl->worldoripool[ensoriid].q3;

			//orientation recategorization utilizes that not all possible orientations are used in the local CA grid thus a more efficient hashtable can be generated however at additional memory costs...
			crxgr.caori = ca_check_disjunctness_core( ensq, true, ensbunge[0], ensbunge[1], ensbunge[2] );
			crxgr.cellcount = 0;
			crxgr.nucsite = luckySite;
			crxgr.startingsite = luckySite;
			crxgr.tincub = myensHdl->worldrxgpool[ensoriid].tincub;

			myrxgpool.push_back ( crxgr );

#ifdef REPORTSTYLE_DEVELOPER
	cout << "\t\tJobID;p;spx;spy;spz;ix;iy;iz;luckySite;nucsite\t\t" << this->jobid << "\t\t" << p << "\t\t" << sp.x << ";" << sp.y << ";" << sp.z << "--" << ix << ";" << iy << ";" << iz << "--" << myrxgpool[myrxgpool.size()-1].nucsite << ";" << myrxgpool[myrxgpool.size()-1].startingsite << endl;
#endif

		} //a point was found and planned for nucleation

	} //test next point

	cout << "Clustering process performed with so many entries = " << myrxgpool.size() << endl;
}


void caHdl::solve_nucmodeling_ellipsoidalcluster_pickrandomly( void )
{
	
cout << "Pick cluster randomly" << endl;

	//generate a large Matern cluster process on the unit cube [0,1]^3
	init_myMaternUniverse();

	//pick a random window from this process and map all points in the automaton domain
	pickrandomly_myMaternUniverse();

//#ifdef REPORTSTYLE_DEVELOPERcout << "\t\tcluster;nuc;jobid;att;nucsite;luckyNucleusWorldRXGPool\t\t" << c << ";" << k << ";" << this->jobid << ";" << ntrials << ";" << myrxgpool[myrxgpool.size()-1].startingsite << ";" << luckyNucleus << endl;#endif

	//characterize the probabilities in the CA that nuclei are able to grow in their parent deformed grain in the surrounding x,y,z minding periodic boundary conditions

#ifdef DETAILED_PROMPTS
	if ( myensHdl->myRank == MASTER ) { cout << this->jobid << "\t\tMatern type ellipsoid nucleation completed." << endl; }
#endif
}


#endif
