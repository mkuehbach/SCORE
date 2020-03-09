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

void caHdl::ompshar_sim_myCA_sitesaturatedNucleation( vector<rxgcommit>* sync )
{
	//this function expects the mycellgrid positions to be marked with mydefgpool ids not tmpdefg seedids!
	//CALLED FROM WITHIN PARALLEL REGIONS BECAUSE OF THAT JUVENILE NUCLEATION INTO MYCELLGRID
	//TAKE CARE THAT CONCURRENT WRITES ARE DISALLOWED
	uint32_t nucleateWhere, nucleateWhereRegion;
	uint32_t card = myCAGeometry.nboxedge_rd;
	uint32_t cardtd = myCAGeometry.nboxarea_rdtd;
	uint32_t z, rem, y, x;

	uint32_t tid = omp_get_thread_num();

	uint32_t xmi = regions[tid]->myGeom.nreg_rdmin;
	uint32_t xmx = regions[tid]->myGeom.nreg_rdmax;
	uint32_t ymi = regions[tid]->myGeom.nreg_tdmin;
	uint32_t ymx = regions[tid]->myGeom.nreg_tdmax;
	uint32_t zmi = regions[tid]->myGeom.nreg_ndmin;
	uint32_t zmx = regions[tid]->myGeom.nreg_ndmax;
	uint32_t xx = regions[tid]->myGeom.nreg_rd;
	uint32_t xxyy = regions[tid]->myGeom.nregarea_rdtd;
	bool toggle = false;

	//scan and place cooperatively the nuclei
	//MK::be aware that it is necessary to distinguish whether to place the nucleus in the inside or the border!

	for ( uint32_t nuc = 0; nuc < myrxgpool.size(); nuc++) {
		nucleateWhere = myrxgpool[nuc].nucsite; //MK::only thread-safe if JuvenileNucleation does not modify myrxgpool

		//scan all coordinates
		z = nucleateWhere / cardtd;
		rem = nucleateWhere - (cardtd * z);
		y = rem / card;
		x = rem - (card * y);

		//check who is with me, mycellgrid of all regions form a nonoverlapping exactly the domain
		//HERE THE THREADS FILTER OUT ONLY THOSE NUCLEI THAT ARE LOCATED IN THEIR OWN CELLGRID
		//AS THE CELLGRIDS ARE NONOVERLAPPING AND FILLING COMPLETELY THE DOMAIN ALL THREADS EXECUTE
		//IN PARALLEL THE SAME CODE BUT NEVER INTERFERE...
		
		if ( x < xmi ) continue; //for performance utilize that most likely the nucleus is outside my own region
		if ( x > xmx ) continue;
		if ( y < ymi ) continue;
		if ( y > ymx ) continue;
		if ( z < zmi ) continue;
		if ( z > zmx ) continue;

		//inside!, interpret global coordinate from CA domain into local coordinate in my CAregion
		nucleateWhereRegion = ((z - zmi) * xxyy) + ((y - ymi) * xx) + (x - xmi);

		//inside or outer shell? most likely inside
		toggle = false; //assume most likely case, nucleus is in the shell
		if ( x > xmi && x < xmx && y > ymi && y < ymx && z > zmi && z < zmx ) {
			toggle = true;
		}

		if ( ompshar_sim_myCA_infect_JuvenileNucleation( nuc, x, y, z, nucleateWhereRegion, toggle, ALMOST_FULLY_INFECTED ) == true ) {
			struct rxgcommit success;
			success.nrxgid = nuc;
			success.nucsite = SITE_ALREADY_CONSUMED; //nucleateWhere;
			success.tincub = t;
			sync->push_back( success);
		}
	}

	if ( myensHdl->myRank == MASTER ) {
		#pragma omp critical
		{
			cout << "\t\tmyRank " << myensHdl->myRank << "; ThreadID " << tid << "; JobID " << this->jobid << "\n";
			cout << "\t\tmydefgpool.size() " << mydefgpool.size() << " myrxgpool.size() " << myrxgpool.size() << " nucleate site-saturated." << endl;
		}
	}
}


bool caHdl::ompshar_sim_myCA_infect_JuvenileNucleation( uint32_t rxgpoolid, uint32_t xglo, uint32_t yglo, uint32_t zglo, uint32_t infectWhereRegion, bool inside, double rxfrac0 )
{
	//CALLED FROM WITHIN PARALLEL REGION!
	//works because the calling thread is only assigned one unique region
	caregionMemHdlP me = this->regions[omp_get_thread_num()];

	//get out here, cell is already infected, but in case of site-saturated nucleation no cell can be infected already
	if ( (me->mycellgrid[infectWhereRegion] >= me->reg_nmydefgpool) || (me->mycellgrid[infectWhereRegion] == CURRENTLY_INFECTED) ) {
		//##MK::synchronization like the following two lines violates thread-safety!
		//##MK::myrxgpool[rxgpoolid].nucsite = SITE_ALREADY_CONSUMED;	//not critical as no thread accesses the same elements of myrxg than another do never overlap! but false sharing... flush not necessary!
		//##MK::myrxgpool[rxgpoolid].tincub = NUCLEUS_NOT_PLACABLE; //MK::startingsite varible must not be changed!
		return false;
	}

	//obviously cell is free, hehe, go for it!
	uint32_t mydgid = me->mycellgrid[infectWhereRegion];
	//nucleation sites disjoint, mark cell as infected to prevent multiple reinfections, do right here to localize memory access
	me->mycellgrid[infectWhereRegion] = CURRENTLY_INFECTED;

	double PP = calc_mobilityweight( rxgpoolid, mydgid );

	uint32_t freeplace = UINT32T_MAX;
	if ( inside == true ) {
		freeplace = me->omp_getNextFreeSlotInRXFrontInside( CELLAPPEND );
		QUICKASSERT ( freeplace != INVALID_ADDRESS);

		me->myRXFrontInside[freeplace].activity = ACTIVE;
		me->myRXFrontInside[freeplace].infector = 26;
		me->myRXFrontInside[freeplace].ix = (short) xglo; //explicit cast to lower precision short allowed because 0 < iglo <= CA_MAXIMUM_DIMENSIONS << SHORT_RANGE
		me->myRXFrontInside[freeplace].iy = (short) yglo;
		me->myRXFrontInside[freeplace].iz = (short) zglo;
		me->myRXFrontInside[freeplace].rxFrac = rxfrac0;
		me->myRXFrontInside[freeplace].mydefgid = mydgid;	//cell carries consumed orientation along
		me->myRXFrontInside[freeplace].myrxgid = rxgpoolid;
		me->myRXFrontInside[freeplace].P = PP;
	}
	else {
		freeplace = me->omp_getNextFreeSlotInRXFrontBorder( CELLAPPEND );
		QUICKASSERT (freeplace != INVALID_ADDRESS);

		me->myRXFrontBorder[freeplace].activity = ACTIVE;
		me->myRXFrontBorder[freeplace].infector = 26;
		me->myRXFrontBorder[freeplace].ix = (short) xglo;
		me->myRXFrontBorder[freeplace].iy = (short) yglo;
		me->myRXFrontBorder[freeplace].iz = (short) zglo;
		me->myRXFrontBorder[freeplace].rxFrac = rxfrac0;
		me->myRXFrontBorder[freeplace].mydefgid = mydgid;	//cell carries consumed orientation along
		me->myRXFrontBorder[freeplace].myrxgid = rxgpoolid;
		me->myRXFrontBorder[freeplace].P = PP;
	}

	if ( PP >= me->reg_myMobilityWeightMax ) me->reg_myMobilityWeightMax = PP;
	me->reg_nmynuclei++;


	//MK::get_NextFreeSlotInRXFront() ASSURES freeplace <= nextSlotNeverActiveRXFront < ntotalRXFront!
	
//#pragma omp critical{cout << "Thr/Nuc = " << omp_get_thread_num() << "\t\t" << rxgpoolid << endl;}
#ifdef REPORTSTYLE_CELLCYCLES
	cout << this->jobid << "\t\tjuvenile infection;rxgpoolid;WhichCelltoInfect;x;y;z;rxfrac0;freeplace" << rxgpoolid << ";" << WhichCelltoInfect << ";" << x << ";" << y << ";" << z << ";" << rxfrac0 << ";" << freeplace << endl;
#endif
	return true;
}


#endif
