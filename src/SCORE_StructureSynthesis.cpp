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

#ifndef SCORE_KERNEL_H
#define SCORE_KERNEL_H

#include "SCORE_Kernel.h"

void caHdl::ompshar_deformedgrains_infect_JuvenileNucleation( uint32_t seedid, uint32_t xglo, uint32_t yglo, uint32_t zglo, uint32_t infectWhereRegion, bool inside, double frac0 )
{
	//CALLED FROM WITHIN PARALLEL REGION!
	//works because the calling thread is only assigned one unique region
	caregionMemHdlP me = regions[omp_get_thread_num()];

	//get out here, there is already a seed at this location
	if ( (me->mycellgrid[infectWhereRegion] != NOT_ASSIGNED_YET) || (me->mycellgrid[infectWhereRegion] == CURRENTLY_INFECTED) ) {
		return;
	}

	//obviously cell is free, hehe, go for it!
	//place a seed at this location: in order to prevent multiple reinfections, do right here to localize memory access
	me->mycellgrid[infectWhereRegion] = CURRENTLY_INFECTED;

	uint32_t freeplace = INVALID_ADDRESS;
	if ( inside == true ) {
		freeplace = me->omp_getNextFreeSlotInSeedFrontInside( CELLAPPEND );
		QUICKASSERT ( freeplace != INVALID_ADDRESS);

		me->mySeedFrontInside[freeplace].activity = ACTIVE;
		me->mySeedFrontInside[freeplace].infector = 26;
		me->mySeedFrontInside[freeplace].ix = (short) xglo; //explicit cast to lower precision short allowed because 0 < iglo <= CA_MAXIMUM_DIMENSIONS << SHORT_RANGE
		me->mySeedFrontInside[freeplace].iy = (short) yglo;
		me->mySeedFrontInside[freeplace].iz = (short) zglo;
		me->mySeedFrontInside[freeplace].frac = frac0;
		me->mySeedFrontInside[freeplace].mydefgseedid = seedid; //cell carries reference to a seed along which in a REPLACE_CA_STRUCTURE BECOMES REPLACED BY DEFORMED GRAIN
	}
	else {
		freeplace = me->omp_getNextFreeSlotInSeedFrontBorder( CELLAPPEND );
		QUICKASSERT (freeplace != INVALID_ADDRESS);

		me->mySeedFrontBorder[freeplace].activity = ACTIVE;
		me->mySeedFrontBorder[freeplace].infector = 26;
		me->mySeedFrontBorder[freeplace].ix = (short) xglo;
		me->mySeedFrontBorder[freeplace].iy = (short) yglo;
		me->mySeedFrontBorder[freeplace].iz = (short) zglo;
		me->mySeedFrontBorder[freeplace].frac = frac0;
		me->mySeedFrontBorder[freeplace].mydefgseedid = seedid;
	}

	//########me->reg_nmynuclei++;

	//MK::get_NextFreeSlotIn..Front() ASSURES freeplace <= nextSlotNeverActive..Front < ntotal..Front!

#ifdef REPORTSTYLE_CELLCYCLES
	cout << this->jobid << "\t\tjuvenile infection;seedid;WhichCelltoInfect;x;y;z;rxfrac0;freeplace" << seedid << ";" << WhichCelltoInfect << ";" << x << ";" << y << ";" << z << ";" << rxfrac0 << ";" << freeplace << endl;
#endif
}


void caHdl::ompshar_placeAllSeeds( void )
{
	//this function expects the mycellgrid positions to be marked with mydefgpool ids not tmpdefg seedids!
	//CALLED FROM WITHIN PARALLEL REGIONS BECAUSE OF THAT JUVENILE NUCLEATION INTO MYCELLGRID
	//TAKE CARE THAT CONCURRENT WRITES ARE DISALLOWED
	uint32_t nucleateWhere, nucleateWhereRegion;
	uint32_t card = myCAGeometry.nboxedge_rd;
	uint32_t cardtd = myCAGeometry.nboxarea_rdtd;
	uint32_t z, rem, y, x;

	uint32_t tid = omp_get_thread_num();

	uint32_t xmi = regions[tid]->myGeom.nreg_rdmin; //inclusive...
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
	for ( uint32_t sd = 0; sd < ndefgseeds; ++sd) {
		nucleateWhere = tmpdefgseeds[sd].location;

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

		//place in the inside or in the border, i.e. outer shell of the local region? most likely inside
		toggle = false; //assume most likely case, nucleus is in the shell
		if ( x > xmi && x < xmx && y > ymi && y < ymx && z > zmi && z < zmx ) {
			toggle = true;
		}

		this->ompshar_deformedgrains_infect_JuvenileNucleation( sd, x, y, z, nucleateWhereRegion, toggle, ALMOST_FULLY_INFECTED );
	}

	if ( myensHdl->myRank == MASTER ) {
		cout << "\t\tmyRank " << myensHdl->myRank << "; ThreadID " << tid << "; JobID " << this->jobid << " ndefgseeds = " << ndefgseeds << " seeds placed." << endl;
	}
}


void caHdl::grow_deformedgrains_voxelize_omp( void )
{
	//grows synthetic deformed grains from the seed.locations with unit speed
	//function is an adaption of solve_RXGROWTH() but made leaner for computational efficiency
	double myomptimer = 0.0;
	unsigned char status = OMP_OPERATIVE;
	uint32_t integrationstep = 0;
	t = 0.0;
	X = 0.0;
	Xcells = 0.0;
	dXstep = 0.0;
	Sv = 0;
	step = 0;
	//only the fully synthesized structure is of interest, so therefore we stripped I/O here, inject here if though desired...
	//also no physics involved in a tessellation with unit speed, so dt also arbitrary
	dt = 1.0;

	//mechanism to enable the writing of HDF5 files (though sequentially) from within a threaded region
	omp_init_lock(&(this->h5lock)); //only the calling thread who sets a lock is allowed to unlock it again, so outside a parallel region locker is the master thread

	//the halo cells have already been initialized, we can utilize them both for this CA and AFTER RESETING again for the main CA
	#pragma omp parallel private( integrationstep, status, myomptimer )
	{
		/*initialize HDF5file
		#pragma omp master
		{
			if ( renderingForThisCA == true && outopt_localrenderhow == RENDERING_MS3D ) {
				if ( outopt_localrenderfileformat == RENDERING_FILEFORMAT_HDF5 && (outopt_localrendercolor == RENDERING_COLOR_GRAINID || outopt_localrendercolor == RENDERING_COLOR_IPFZ || outopt_localrendercolor == RENDERING_COLOR_RHO) ) {
					string h5fname = "SCORE." + std::to_string( this->myensHdl->simid ) + "." + "VoxelData.h5";
					bool h5status = true;
					omp_set_lock(&h5lock);
					h5status = hdf5_write_coordinates( h5fname, "/coordinates" );
					omp_unset_lock(&h5lock);
				}
			}
		}
		#pragma omp barrier*/

		integrationstep = 0;
		status = OMP_OPERATIVE;
		QUICKASSERT ( regions.size() == omp_get_num_threads() );

		regions[omp_get_thread_num()]->ompshar_init_mySeedFrontInside();
		regions[omp_get_thread_num()]->ompshar_init_mySeedFrontBorder();

		//MK::a barrier is necessary!, no placing of nuclei before RX front cell management structures have been initialized
		#pragma omp barrier

		//add nuclei into the simulation, per se always site saturated
		ompshar_placeAllSeeds();

		//now we require a barrier as we cannot start before all seeds are set surplus barrier implicitly flushes
		#pragma omp barrier

		#pragma omp master
		{
			for (uint32_t r = 0; r < regions.size(); r++)
				nmynuclei = nmynuclei + regions[r]->reg_nmynuclei;
		}
		//barrier necessary as omp master has no implicit one and otherwise the threads would start iterating already
		#pragma omp barrier


		//MAIN LOOP needed transformation as do/while not allowed in parallel region do {//write the correct value into the loop
		for ( integrationstep = 0; integrationstep < DEFAULT_NTSTEPSMAX_DEFSYNTH; integrationstep++ ) {
			status = OMP_OPERATIVE;
			//each thread reads globally synchronized variable and checks whether or not there is still something to do, because we are not allowed to break out uncontrolled of a parallel for construct
			if ( Xcells >= myCAGeometry.nboxvol_rdtdnd ) 
				status = OMP_COMPLETED;

			if ( status == OMP_OPERATIVE ) {
				//fill ACTIVE cells while keeping track of INACTIVE cells in myRecyclingList, totally local and independent of others
				regions[omp_get_thread_num()]->ompshar_voxelize_growthStepInside();
				regions[omp_get_thread_num()]->ompshar_voxelize_growthStepBorder();

				//let the fully synthesized ACTIVE cells now infect their neighbors
				//MK::inside is independent within the cellregion of one thread requiring no synchronization, as inside cells cannot in Moore neighborhood protrude beyond local mycellgrids at most only into local border
				//in this function cells on mySeedFront become activated again, but during updating in one thread all interesting positions NOT ADDRESSES! are already known in the FullSeed Lists
				//so it cannot happen that reactivated cells infect out of cycle...
				regions[omp_get_thread_num()]->ompshar_voxelize_updateFullInside();

				//possible infections outside this region are stored in the halo of the region or rejected directly by the halos
				//the halo hence not only acts as the synchronization element but also to contribute to an overall earlier rejection of attempts to infect
				//cells which have already been infected or transformed
				//MK::the halo is reutilized for both the VOXELIZING of the initial structure but also the solve_RXGROWTH CA itself!
				//MK::struct halocell t.myrxgid encodes either a reference to a seed tmpdefgseeds or a reference to a recrystallizing grain
				regions[omp_get_thread_num()]->ompshar_clear_halo_bookkeeping();
				regions[omp_get_thread_num()]->ompshar_voxelize_updateFullBorder();
				#pragma omp barrier
				//MK::necessary! because all threads should finish populating the halo regions before threads read out the halos in parallel to synchronize
				//they are guided by neighboring regions' haloref list which stores which positions were newly infected by the neighboring thread in each halo
				regions[omp_get_thread_num()]->ompshar_synchronize_haloregions_def();

				if ( outopt_localthreadprof == OUTPUT_THREADPROFILING_YES )
					regions[omp_get_thread_num()]->omp_log_synthfrontstats();

				#pragma omp barrier //##MK::probably necessary...
			} //growth step completed...
			myomptimer = 0.0;

			//...now all synchronization work on the master as writing ot global variables in memory shared by all threads
			#pragma omp master
			{
				if ( status == OMP_OPERATIVE ) {
					dXstep = 0.0;
					Sv = 0;
					for ( uint32_t tid = 0; tid < regions.size(); tid++ ) { //accumulate local results
						dXstep = dXstep + ( regions[tid]->reg_dXstepInside + regions[tid]->reg_dXstepBorder );
						Sv = Sv + (regions[tid]->reg_SvInside + regions[tid]->reg_SvBorder);
						Xcells = Xcells + ( (double) (regions[tid]->reg_nUpdatedCellsInside + regions[tid]->reg_nUpdatedCellsBorder) ); 
					}
					X = (Xcells / ((double) myCAGeometry.nboxvol_rdtdnd));

					//####unnecessarysynchronize_myrxgpool_counts();synchronize_mydefpool_counts();

					//neither physics updating nor time increment updating necessary as we continue to migrate with unit speed
					t += dt;
					step++;

//#ifdef REPORTSTYLE_USER
					if ( this->myensRank == MASTER ) 
						cout << "\tstep;t;X;Sv;dt\t\t" << this->step << "\t\t" << this->t << "\t\t" << this->X  << "\t\t" << this->Sv << "\t\t" << this->dt << endl;
//#endif
				}
			}
			#pragma omp barrier //necessary as otherwise threads would continue processing the next time step with wrong values for t and step!

			regions[omp_get_thread_num()]->profiling_synthmachine[regions[omp_get_thread_num()]->profiling_synthmachine.size()-1].tSeqOverhead = omp_get_wtime() - myomptimer;
		}
		//OMP_COMPLETED

		#pragma omp master
		{
			tsimend = t;
			stepsimend = step;

			if ( renderingForThisCA == true ) {
				if ( outopt_localrenderhow == RENDERING_MS3D ) {
					if ( outopt_localrenderfileformat == RENDERING_FILEFORMAT_RAW ) {
						if ( outopt_localrendercolor == RENDERING_COLOR_GRAINID )
							ompcrit_write_voxeldata_coloring_grainids( "SYNTH.3DGrainID" );
						else
							cout << "ERROR::3D output in RAW format currently supports only grain id coloring!" << endl;
					}
					if ( outopt_localrenderfileformat == RENDERING_FILEFORMAT_HDF5 ) {
						if ( outopt_localrendercolor == RENDERING_COLOR_GRAINID )
							ompcrit_write_voxeldata_h5( RENDERING_COLOR_GRAINID, "/GROWTH3DGrainID", false, "dummy", 0.0, 0.0, 0 );
						else if ( outopt_localrendercolor == RENDERING_COLOR_IPFZ ) {
							cout << "ERROR::3D output in HDF5 format currently not supported for the synthesizes structure in combination with IPFZ!" << endl;
							//MK::the current implementation of the ompcrit_write_voxeldata_h5 function with the option RENDERING_COLOR_IPFZ writes out raw files
							//as the HDF5 library since definately 1.8.12 does not provide in all cases the H5T_STD_U8LE type any longer
							//the current way the core function mpiio_... is implemented in SCORE must not be called before the CA structure was REPLACED...
							//as prior to this CA_REPLACEMENT function the ids in the local mycellgrids refer to tmpdefgseed ids, not mydefgpool IDs or myrxgpool ids!
						}
						else if ( outopt_localrendercolor == RENDERING_COLOR_RHO )
							ompcrit_write_voxeldata_h5( RENDERING_COLOR_RHO, "/GROWTH3DGrainID", false , "dummy", 0.0, 0.0, 0.0 );
						else
							cout << "ERROR::3D output in HDF5 format currently supports only either grain id or IPFZ coloring!" << endl;
					}
				}
			}
		}
		#pragma omp barrier

		regions[omp_get_thread_num()]->omp_voxelize_memoryCleanup();
		regions[omp_get_thread_num()]->omp_resetHaloRegions(); //MK::HALO REGIONS MUST NOT BE CLEANEDUP AS THEY ARE REUTILIZED

		#pragma omp master
		{
			regions[omp_get_thread_num()]->omp_resetInternalCounters();
		}

	} //implicit barrier and flush parallel execution done, go out when completely recrystallized, recrystallized to use defined value, when too many integration steps, when processing out
	omp_destroy_lock(&(this->h5lock));

	cout << myensHdl->myRank << " myRank,jobid = " << this->jobid << " local CA deformation structure was synthetized with " << ndefgseeds  << endl;

	if ( outopt_localthreadprof == OUTPUT_THREADPROFILING_YES )
		write_ThreadProfilingSummarySynthMachine();

	//##MK::output deformation structure//write_voxeldata_coloring_grainids();
}

void caHdl::grow_deformedgrains_voxelize_columns( void )
{
	//stretches an SEM/EBSD map along the 3d dimension such that the grain along z
	//this->myensHdl->expData
	//####
}


void caHdl::determine_polycrystalline_geometry( void )
{
	//synthesize polycrystalline aggregate of cuboids displace the lower left front origin of this aggregate randomly 
	//to the lower left front origin of the SU coordinate system to enable random sampling of the domain
	myCADefMS.u_xrd = -1.0 * localprng.MersenneTwister() * myPhysData.defgmean_rd;
	myCADefMS.v_ytd = -1.0 * localprng.MersenneTwister() * myPhysData.defgmean_td;
	myCADefMS.w_znd = -1.0 * localprng.MersenneTwister() * myPhysData.defgmean_nd;

	//therewith, determine required size of the 3D grain grid, which encloses the SU domain completely
	uint32_t nx = ( (myCAGeometry.boxedge_rd - myCADefMS.u_xrd) / myPhysData.defgmean_rd ); //will be zero if defg+u larger than box, so always add one grain
	nx++; 
	uint32_t ny = ( (myCAGeometry.boxedge_td - myCADefMS.v_ytd) / myPhysData.defgmean_td );
	ny++;
	uint32_t nz = ( (myCAGeometry.boxedge_nd - myCADefMS.w_znd) / myPhysData.defgmean_nd );
	nz++;

	myCADefMS.ngrx = nx;
	myCADefMS.ngry = ny;
	myCADefMS.ngrz = nz;
	myCADefMS.ngrxy = nx * ny;
	myCADefMS.ngrxyz = nx * ny * nz;

	//MK::allocate memory for the topology of the deformed cuboid-shaped grains, as grains from the input can be used multiple times but
	//the boundary detection requires grains to be of disjoint IDs it is absolutely necessary to keep a list of the seeds to which any grain are referenced
	ndefgseeds = myCADefMS.ngrxyz;

	QUICKASSERT ( ndefgseeds < CA_ALLOCATION_MAXIMUM );
	tmpdefgseeds = NULL;
	tmpdefgseeds = new struct defgseed[ndefgseeds];
	QUICKASSERT ( tmpdefgseeds != NULL );
	myMemGuard = myMemGuard + (ndefgseeds * sizeof(struct defgseed));
	for ( uint32_t sd = 0; sd < ndefgseeds; ++sd ) {
		tmpdefgseeds[sd].ensdefgpoolid = NOT_ASSIGNED_YET;
		tmpdefgseeds[sd].mydefgpoolid = I_DONT_KNOW_YET;
		tmpdefgseeds[sd].location = NOT_ASSIGNED_YET;
	}
}


void caHdl::pick_deformedgrains( void )
{
	//MK::everything is possible, GIA-aggregate packing, MODF optimization or simply random picking
	//##MK::here we start with sampling by picking randomly from the worlddefgpool grains resulting in an uncorrelated MODF
	ensembleHdlP ens = this->myensHdl;
	uint32_t nensdgpool = ens->worlddefgpool.size();
	uint32_t luckyGrain;

	for ( uint32_t dg = 0; dg < ndefgseeds; dg++ ) {
		luckyGrain = localprng.MersenneTwister() * nensdgpool; //MK::lucky grain refers to IDs from the ensembleHdl defgpool

		//MK::ASSIGNMENT IS WITH WORLDDEFGPOOL-ID, mind this is in contrast to mydefggrid where IDs refer to this->mydefgpool !
		tmpdefgseeds[dg].ensdefgpoolid = luckyGrain; //MK::to avoid carrying always the full pool of grains per SU, we cache locally in caHdl::mydefgpool

		//append only to mydefgpool if we have not picked this grain before therefore scan only the known ones
		bool DoIKnowLuckyGrain = false;
		for ( uint32_t tg = 0; tg < dg; tg++) {
			if ( luckyGrain == tmpdefgseeds[tg].ensdefgpoolid ) {
				DoIKnowLuckyGrain = true;
				break;
			}
		}

//cout << "2-dg;luckyGrain;DoIKnowLuckyGrain\t\t" << dg << ";" << luckyGrain << ";" << DoIKnowLuckyGrain << endl;

		if ( DoIKnowLuckyGrain == false ) { //most likely this condition is met when ensHdl->defgpool is large
			struct cadefg dgr;

			uint32_t ensoriid = ens->worlddefgpool[luckyGrain].ori; //pull orientation from global pool
			double ensq[4] = { ens->worldoripool[ensoriid].q0, ens->worldoripool[ensoriid].q1, ens->worldoripool[ensoriid].q2, ens->worldoripool[ensoriid].q3 };
			//orientation recategorization with known Bunge values in local caHdl to utilize that not all possible orientations are used in the local CA grid thus a more efficient hashtable can be generated however at additional memory costs...
			dgr.worlddefgid = luckyGrain;
			dgr.caori = ca_check_disjunctness_core( ensq, true, ens->worldoripool[ensoriid].bunge1, ens->worldoripool[ensoriid].bunge2, ens->worldoripool[ensoriid].bunge3 ); //checking against truely the closest avoids drifting orientation assignment
			dgr.cellcount = 0;
			dgr.rho0 = ens->worlddefgpool[luckyGrain].rho0;
			dgr.rho = ens->worlddefgpool[luckyGrain].rho;
			//if (dgr.rho0 >= this->myrhomax) { myrhomax = dgr.rho0; }//identify highest dislocation density in the SU for customizing an adaptive iteration scheme

			mydefgpool.push_back ( dgr );

			//add grain to maintenance structure which keeps track a reduction rho to speed up the explicit integration
			struct dgopt anopt;
			anopt.rho = dgr.rho0;
			anopt.id = mydefgpool.size() - 1;
			anopt.cnt = 0;
			mydgoptimize.push_back( anopt );

			//as I dont know this grain assign the copy of the grain local to caHdl to the pool
			tmpdefgseeds[dg].mydefgpoolid = (mydefgpool.size() - 1);

//cout << "3-dg;luckyGrain;UNKNOWN" << dg << ";" << luckyGrain << ";" << DoIKnowLuckyGrain << ";" << tmpdefgseeds[dg].ensdefgpoolid << "__" << tmpdefgseeds[dg].mydefgpoolid << endl;
			continue;
		}

		//else DoIKnowLuckyGrain == true, i.e. I picked the grain already a while ago well then find it in 
		//mylist utilizing that worlddefgid are >= 0 and unique
		for ( uint32_t testg = 0; testg < mydefgpool.size(); ++testg) {
			if ( luckyGrain == mydefgpool[testg].worlddefgid ) {
				tmpdefgseeds[dg].mydefgpoolid = testg;
//cout << "4-dg;luckyGrain;KNOWN" << dg << ";" << luckyGrain << ";" << DoIKnowLuckyGrain << ";" << "testgrain=" << testg << "\t\t" << tmpdefgseeds[dg].ensdefgpoolid << "__" << tmpdefgseeds[dg].mydefgpoolid << endl;
				break;
			}
		}
	}

	//##calculate local MODF if desired here
#ifdef REPORTSTYLE_DEVELOPER
	cout << this->jobid << "\t\t" << this->myensHdl->myRank << "\t\t" << myMemGuard << "\t\t" << mydefgpool.size() << "\t\t" << "\t\tjobid;myensHdlmyRank;myMemGuard (Byte);mydefgpool.size() successfully packed GIA grains." << endl;
#endif
}


void caHdl::ompshar_voxelize_defms( void )
{
	//if ( this->myensHdl->myRank == MASTER && omp_get_thread_num() == MASTER ) cout << "\t\tVoxelizing SU domain." << endl;
	//CALLED FROM WITHIN PARALLEL REGION, function call so variables private because on thread local function stack!

	//MK::a memory marking intensive function + number crunching
	//in which each voxel a reference to a deformed grain is assigned depending on in which cuboid grain the cell center ix,iy,iz is
	//implicit storage of 3D array: align positive x's successively, stack x lines positive in y and stack these xy slices in positive z direction

	caregionMemHdlP reg = this->regions[omp_get_thread_num()];

	uint32_t xmin = reg->myGeom.nreg_rdmin;
	uint32_t xmax = reg->myGeom.nreg_rdmax;
	uint32_t ymin = reg->myGeom.nreg_tdmin;
	uint32_t ymax = reg->myGeom.nreg_tdmax;
	uint32_t zmin = reg->myGeom.nreg_ndmin;
	uint32_t zmax = reg->myGeom.nreg_ndmax;

	uint32_t ngrx = myCADefMS.ngrx;
	uint32_t ngrxy = myCADefMS.ngrxy;

	double cs = reg->myGeom.cellsize;
	double cs2 = 0.5 * cs;
	double ux = myCADefMS.u_xrd;	//already absolute displacement in m of SU coordinate system to cuboid deformed grain grid!
	double uy = myCADefMS.v_ytd;
	double uz = myCADefMS.w_znd;

	double dgrd = myPhysData.defgmean_rd;
	double dgtd = myPhysData.defgmean_td;
	double dgnd = myPhysData.defgmean_nd;

	//now everything is between the stack and two linear arrays...
	uint32_t crxyz = 0;
	uint32_t gz, gz_ngrxy, gy, gy_ngrx, gx, gxyz;

//cout << "ngrxyz\t\t" << ngrx << ";" << ngry << ";" << ngrz << ";" << ngrxy << ";" << ngrxyz << endl;
	//bool validassignment = true;
	for (uint32_t iz = zmin; iz <= zmax; ++iz) {
		gz = ((( (double) iz * cs) + cs2 - uz) / dgnd); ///####testing!
		gz_ngrxy = gz * ngrxy;

		for (uint32_t iy = ymin; iy <= ymax; ++iy) {
			gy = ((( (double) iy * cs) + cs2 - uy) / dgtd);
			gy_ngrx = gy * ngrx;

			for (uint32_t ix = xmin; ix <= xmax; ++ix) {
				gx = ((( (double) ix * cs) + cs2 - ux) / dgrd);
				gxyz = gx + gy_ngrx + gz_ngrxy;

				//if ( gxyz >= ndefgseeds) 
				//	validassignment = false;

				reg->mycellgrid[crxyz] = gxyz;
//#pragma omp critical {cout << "iz;iy;ix;gz;gy;gx;gz_ngrxy;gy_ngrx;gxyz;tmpdefgsseeds;crxyz;mycellgridcrxyz\t\t" << iz << ";" << iy << ";" << ix << ";" << gz << ";" << gy << ";" << gx << ";" << gz_ngrxy << ";" << gy_ngrx << ";" << gxyz << ";" << tmpdefgseeds[gxyz].mydefgpoolid << ";" << crxyz << ";" << reg->mycellgrid[crxyz] << endl;}

				crxyz++;
			}
		} //stack xlines
	} //stack xy slices

	//if ( crxyz != (reg->myGeom.nreg_rd * reg->myGeom.nreg_nd * reg->myGeom.nreg_td) )
	//	validassignment = false;
}

#endif