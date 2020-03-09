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


void caHdl::write_rxfrontstats( void )
{
	//I/O interaction
	stringstream log_rxfrontstats_fname;
	ofstream log_rxfrontstats_file;
	
	log_rxfrontstats_fname << "SCORE." << myensHdl->simid << ".MyRank." << myensHdl->myRank << ".JobID." << this->jobid << ".RXFrontStats.csv";
#ifdef REPORTSTYLE_DEVELOPER
	cout << "File " << log_rxfrontstats_fname.str().c_str() << " is opened now" << endl;
#endif
	log_rxfrontstats_file.open ( log_rxfrontstats_fname.str().c_str() );

	//header
	log_rxfrontstats_file << "Step;Time/s;X;Sv;myMemory/Byte;currPmax;ntotalRXFront;nextSlotNeverActiveRXFront;nCurrentlyActive;ntotalFullRXList;nextSlotToFullRX;ntotalRecyclingList;nextSlotThatBecomesRecycled;firstNotRecycledYet" << endl;

	for ( uint32_t s = 0; s < myrxfrontstatus.size(); ++s ) {
		log_rxfrontstats_file << myrxfrontstatus[s].localstep << ";" << myrxfrontstatus[s].localtime << ";" << myrxfrontstatus[s].localX << ";" << myrxfrontstatus[s].localSv << ";" << myrxfrontstatus[s].localmemory << ";" << myrxfrontstatus[s].localPmax << ";" << myrxfrontstatus[s].ntotalRXFront << ";" << myrxfrontstatus[s].nextSlotNeverActiveRXFront << ";" << myrxfrontstatus[s].nCurrentlyActive << ";" << myrxfrontstatus[s].ntotalFullRXList << ";" << myrxfrontstatus[s].nextSlotToFullRX << ";" << myrxfrontstatus[s].ntotalRecyclingList << ";" << myrxfrontstatus[s].nextSlotThatBecomesRecycled << ";" << myrxfrontstatus[s].firstNotRecycledYet << endl;
	}

	log_rxfrontstats_file.flush();
	log_rxfrontstats_file.close();
}


void caHdl::write_ThreadProfilingSummarySynthMachine( void )
{
	//I/O interaction ASCII
	stringstream log_tprof_fname;
	ofstream log_tprof_file;
	
	log_tprof_fname << "SCORE." << myensHdl->simid << ".MyRank." << myensHdl->myRank << ".JobID." << this->jobid << ".ThreadProfilingSynthesis.csv";
#ifdef REPORTSTYLE_DEVELOPER
	cout << "File " << log_tprof_fname.str().c_str() << " is opened now" << endl;
#endif
	log_tprof_file.open ( log_tprof_fname.str().c_str() );

	//MK::many cache-misses and much loop overhead but then directly all results ready to be imported in OriginLAB...

	uint32_t nreg = this->regions.size();

	//header
	log_tprof_file << "IntegrationStep;X;";
	for ( uint32_t props = 0; props < THREAD_PROFILING_NPROPS_SYNTH; props++ ) { 
		for ( uint32_t r = 0; r < nreg; r++ ) log_tprof_file << ";" << "Thread-" << r;
	}
	log_tprof_file << endl;

	//property values for all threads next to one another rather than first all properties of the first, the second,...
	log_tprof_file << "IntegrationStep;X;";
	for ( uint32_t r = 0; r < nreg; r++ ) log_tprof_file << ";CalcGrowthInside";
	for ( uint32_t r = 0; r < nreg; r++ ) log_tprof_file << ";CalcGrowthBorder";
	for ( uint32_t r = 0; r < nreg; r++ ) log_tprof_file << ";UpdateFullInside";
	for ( uint32_t r = 0; r < nreg; r++ ) log_tprof_file << ";UpdateFullBorder";
	for ( uint32_t r = 0; r < nreg; r++ ) log_tprof_file << ";SyncHalos";
	for ( uint32_t r = 0; r < nreg; r++ ) log_tprof_file << ";SeqOverhead";
	for ( uint32_t r = 0; r < nreg; r++ ) log_tprof_file << ";reg_SvInside";
	for ( uint32_t r = 0; r < nreg; r++ ) log_tprof_file << ";reg_SvBorder";
	for ( uint32_t r = 0; r < nreg; r++ ) log_tprof_file << ";ntotalInside";
	for ( uint32_t r = 0; r < nreg; r++ ) log_tprof_file << ";ntotalBorder";
	for ( uint32_t r = 0; r < nreg; r++ ) log_tprof_file << ";firstNeverActiveInside";
	for ( uint32_t r = 0; r < nreg; r++ ) log_tprof_file << ";firstNeverActiveBorder";
	//for ( uint32_t r = 0; r < nreg; r++ ) log_tprof_file << ";nowActiveInside";
	//for ( uint32_t r = 0; r < nreg; r++ ) log_tprof_file << ";nowActiveBorder";
	for ( uint32_t r = 0; r < nreg; r++ ) log_tprof_file << ";fragmentationInside";
	for ( uint32_t r = 0; r < nreg; r++ ) log_tprof_file << ";fragmentationBorder";
	log_tprof_file << endl;

	//how many steps and consistent dataset?
	uint32_t nlogpoints = this->regions[MASTER]->profiling_synthmachine.size();
	for ( uint32_t r = 0; r < nreg; r++ ) { QUICKASSERT( this->regions[r]->profiling_synthmachine.size() == nlogpoints ); } //##DEBUG

	cout << "ThreadProfiling with " << nlogpoints << endl << endl;

	for ( uint32_t lp = 0; lp < nlogpoints; lp++ ) {
		log_tprof_file << this->regions[MASTER]->profiling_synthmachine[lp].localstep << ";" << setprecision(16) << this->regions[MASTER]->profiling_synthmachine[lp].localX << ";";

		//for all threads make profiling data matrix
		for ( uint32_t r = 0; r < nreg; r++ ) log_tprof_file << ";" << setprecision(8) << this->regions[r]->profiling_synthmachine[lp].tCalcGrowthInside;
		for ( uint32_t r = 0; r < nreg; r++ ) log_tprof_file << ";" << setprecision(8) << this->regions[r]->profiling_synthmachine[lp].tCalcGrowthBorder;
		for ( uint32_t r = 0; r < nreg; r++ ) log_tprof_file << ";" << setprecision(8) << this->regions[r]->profiling_synthmachine[lp].tUpdateInside;
		for ( uint32_t r = 0; r < nreg; r++ ) log_tprof_file << ";" << setprecision(8) << this->regions[r]->profiling_synthmachine[lp].tUpdateBorder;
		for ( uint32_t r = 0; r < nreg; r++ ) log_tprof_file << ";" << setprecision(8) << this->regions[r]->profiling_synthmachine[lp].tSyncHalos;
		for ( uint32_t r = 0; r < nreg; r++ ) log_tprof_file << ";" << setprecision(8) << this->regions[r]->profiling_synthmachine[lp].tSeqOverhead;
		for ( uint32_t r = 0; r < nreg; r++ ) log_tprof_file << ";" << setprecision(8) << this->regions[r]->profiling_synthmachine[lp].regionSvInside;
		for ( uint32_t r = 0; r < nreg; r++ ) log_tprof_file << ";" << setprecision(8) << this->regions[r]->profiling_synthmachine[lp].regionSvBorder;
		for ( uint32_t r = 0; r < nreg; r++ ) log_tprof_file << ";" << setprecision(8) << this->regions[r]->profiling_synthmachine[lp].ntotalSeedInside;
		for ( uint32_t r = 0; r < nreg; r++ ) log_tprof_file << ";" << setprecision(8) << this->regions[r]->profiling_synthmachine[lp].ntotalSeedBorder;
		for ( uint32_t r = 0; r < nreg; r++ ) log_tprof_file << ";" << setprecision(8) << this->regions[r]->profiling_synthmachine[lp].firstNeverActiveSeedInside;
		for ( uint32_t r = 0; r < nreg; r++ ) log_tprof_file << ";" << setprecision(8) << this->regions[r]->profiling_synthmachine[lp].firstNeverActiveSeedBorder;
		//for ( uint32_t r = 0; r < nreg; r++ ) log_tprof_file << ";" << setprecision(8) << this->regions[r]->profiling_synthmachine[lp].nowActiveSeedInside;
		//for ( uint32_t r = 0; r < nreg; r++ ) log_tprof_file << ";" << setprecision(8) << this->regions[r]->profiling_synthmachine[lp].nowActiveSeedBorder;
		for ( uint32_t r = 0; r < nreg; r++ ) {
			double fraginside = 0.0;
			if ( this->regions[r]->profiling_synthmachine[lp].firstNeverActiveSeedInside > 0 )
				fraginside = (double) this->regions[r]->profiling_synthmachine[lp].regionSvInside / (double) this->regions[r]->profiling_synthmachine[lp].firstNeverActiveSeedInside;
			log_tprof_file << ";" << setprecision(8) << (1.0 - fraginside);
		}
		for ( uint32_t r = 0; r < nreg; r++ ) {
			double fragborder = 0.0;
			if ( this->regions[r]->profiling_synthmachine[lp].firstNeverActiveSeedBorder > 0 )
				fragborder = (double) this->regions[r]->profiling_synthmachine[lp].regionSvBorder / (double) this->regions[r]->profiling_synthmachine[lp].firstNeverActiveSeedBorder;
			log_tprof_file << ";" << setprecision(8) << (1.0 - fragborder);
		}
		log_tprof_file << "\n";
	}


	log_tprof_file.flush();
	log_tprof_file.close();
}


void caHdl::write_ThreadProfilingSummaryGrowthMachine( void )
{
	//I/O interaction ASCII
	stringstream log_tprof_fname;
	ofstream log_tprof_file;
	
	log_tprof_fname << "SCORE." << myensHdl->simid << ".MyRank." << myensHdl->myRank << ".JobID." << this->jobid << ".ThreadProfilingGrowth.csv";
#ifdef REPORTSTYLE_DEVELOPER
	cout << "File " << log_tprof_fname.str().c_str() << " is opened now" << endl;
#endif
	log_tprof_file.open ( log_tprof_fname.str().c_str() );

	//MK::many cache-misses and much loop overhead but then directly all results ready to be imported in OriginLAB...

	uint32_t nreg = this->regions.size();

	//header
	log_tprof_file << "Step;Time/s;X;";
	for ( uint32_t props = 0; props < THREAD_PROFILING_NPROPS_GROWTH; props++ ) { 
		for ( uint32_t r = 0; r < nreg; r++ ) log_tprof_file << ";" << "Thread-" << r;
	}
	log_tprof_file << endl;

	//property values for all threads next to one another rather than first all properties of the first, the second,...
	log_tprof_file << "Step;Time/s;X;";
	for ( uint32_t r = 0; r < nreg; r++ ) log_tprof_file << ";CalcGrowthInside";
	for ( uint32_t r = 0; r < nreg; r++ ) log_tprof_file << ";CalcGrowthBorder";
	for ( uint32_t r = 0; r < nreg; r++ ) log_tprof_file << ";UpdFullRXInside";
	for ( uint32_t r = 0; r < nreg; r++ ) log_tprof_file << ";UpdFullRXBorder";
	for ( uint32_t r = 0; r < nreg; r++ ) log_tprof_file << ";SyncHalos";
	for ( uint32_t r = 0; r < nreg; r++ ) log_tprof_file << ";Defragmentation";
	for ( uint32_t r = 0; r < nreg; r++ ) log_tprof_file << ";SeqOverhead";
	for ( uint32_t r = 0; r < nreg; r++ ) log_tprof_file << ";reg_SvInside";
	for ( uint32_t r = 0; r < nreg; r++ ) log_tprof_file << ";reg_SvBorder";
	for ( uint32_t r = 0; r < nreg; r++ ) log_tprof_file << ";nextSlotNeverActiveRXFrontInside";
	for ( uint32_t r = 0; r < nreg; r++ ) log_tprof_file << ";nextSlotNeverActiveRXFrontBorder";
	for ( uint32_t r = 0; r < nreg; r++ ) log_tprof_file << ";fragmentationDegreeInside";
	for ( uint32_t r = 0; r < nreg; r++ ) log_tprof_file << ";fragmentationDegreeBorder";
	log_tprof_file << endl;

	//how many steps and consistent dataset?
	uint32_t nlogpoints = this->regions[MASTER]->profiling_growthmachine.size();
	for ( uint32_t r = 0; r < nreg; r++ ) { QUICKASSERT( this->regions[r]->profiling_growthmachine.size() == nlogpoints ); } //##DEBUG

	cout << "ThreadProfiling with " << nlogpoints << endl << endl;

	for ( uint32_t lp = 0; lp < nlogpoints; lp++ ) {
		log_tprof_file << this->regions[MASTER]->profiling_growthmachine[lp].localstep << ";" << setprecision(16) << this->regions[MASTER]->profiling_growthmachine[lp].localtime << ";" << setprecision(16) << this->regions[MASTER]->profiling_growthmachine[lp].localX << ";";

		//for all threads make profiling data matrix
		for ( uint32_t r = 0; r < nreg; r++ ) log_tprof_file << ";" << setprecision(8) << this->regions[r]->profiling_growthmachine[lp].tCalcGrowthInside;
		for ( uint32_t r = 0; r < nreg; r++ ) log_tprof_file << ";" << setprecision(8) << this->regions[r]->profiling_growthmachine[lp].tCalcGrowthBorder;
		for ( uint32_t r = 0; r < nreg; r++ ) log_tprof_file << ";" << setprecision(8) << this->regions[r]->profiling_growthmachine[lp].tUpdateInside;
		for ( uint32_t r = 0; r < nreg; r++ ) log_tprof_file << ";" << setprecision(8) << this->regions[r]->profiling_growthmachine[lp].tUpdateBorder;
		for ( uint32_t r = 0; r < nreg; r++ ) log_tprof_file << ";" << setprecision(8) << this->regions[r]->profiling_growthmachine[lp].tSyncHalos;
		for ( uint32_t r = 0; r < nreg; r++ ) log_tprof_file << ";" << setprecision(8) << this->regions[r]->profiling_growthmachine[lp].tDefragment;
		for ( uint32_t r = 0; r < nreg; r++ ) log_tprof_file << ";" << setprecision(8) << this->regions[r]->profiling_growthmachine[lp].tSeqOverhead;
		for ( uint32_t r = 0; r < nreg; r++ ) log_tprof_file << ";" << setprecision(8) << this->regions[r]->profiling_growthmachine[lp].regionSvInside;
		for ( uint32_t r = 0; r < nreg; r++ ) log_tprof_file << ";" << setprecision(8) << this->regions[r]->profiling_growthmachine[lp].regionSvBorder;
		for ( uint32_t r = 0; r < nreg; r++ ) log_tprof_file << ";" << setprecision(8) << this->regions[r]->profiling_growthmachine[lp].nextSlotNeverActiveRXFrontInside;
		for ( uint32_t r = 0; r < nreg; r++ ) log_tprof_file << ";" << setprecision(8) << this->regions[r]->profiling_growthmachine[lp].nextSlotNeverActiveRXFrontBorder;
		for ( uint32_t r = 0; r < nreg; r++ ) {
			double fraginside = 0.0;
			if ( this->regions[r]->profiling_growthmachine[lp].nextSlotNeverActiveRXFrontInside > 0 )
				fraginside = (double) this->regions[r]->profiling_growthmachine[lp].regionSvInside / (double) this->regions[r]->profiling_growthmachine[lp].nextSlotNeverActiveRXFrontInside;
			log_tprof_file << ";" << setprecision(8) << (1.0 - fraginside);
		}
		for ( uint32_t r = 0; r < nreg; r++ ) {
			double fragborder = 0.0;
			if ( this->regions[r]->profiling_growthmachine[lp].nextSlotNeverActiveRXFrontBorder > 0 )
				fragborder = (double) this->regions[r]->profiling_growthmachine[lp].regionSvBorder / (double) this->regions[r]->profiling_growthmachine[lp].nextSlotNeverActiveRXFrontBorder;
			log_tprof_file << ";" << setprecision(8) << (1.0 - fragborder);
		}		log_tprof_file << "\n";
	}


	log_tprof_file.flush();
	log_tprof_file.close();
}


void caHdl::write_grainevolution_ascii( void )
{
	stringstream log_grainevo_fname;
	ofstream log_grainevo_file;
	
	log_grainevo_fname << "SCORE." << myensHdl->simid << ".MyRank." << myensHdl->myRank << ".JobID." << this->jobid << ".SingleGrainData.csv";
#ifdef REPORTSTYLE_DEVELOPER
	cout << "File " << log_grainevo_fname.str().c_str() << " is opened now" << endl;
#endif
	log_grainevo_file.open ( log_grainevo_fname.str().c_str() );

	//header
	log_grainevo_file << "GrainID/Step;;;;;;;;";
	for ( uint32_t s = 0; s < loginfo_grainevo_cnt; s++ ) { log_grainevo_file << ";" << output_atthisX[s]; }
	log_grainevo_file << endl;
	log_grainevo_file << "Step;;;;;;;;";
	for ( uint32_t s = 0; s < loginfo_grainevo_cnt; s++ ) { log_grainevo_file << ";" << mygrainevolution[s].localtstep; }
	log_grainevo_file << endl;
	log_grainevo_file << "Time;;;;;;;;";
	for ( uint32_t s = 0; s < loginfo_grainevo_cnt; s++ ) { log_grainevo_file << ";" << mygrainevolution[s].localtime; }
	log_grainevo_file << endl;
	log_grainevo_file << "X;;;;;;;;";
	for ( uint32_t s = 0; s < loginfo_grainevo_cnt; s++ ) { log_grainevo_file << ";" << mygrainevolution[s].localX; }
	log_grainevo_file << endl;
	log_grainevo_file << "Sv;;;;phi1;PHI;phi2;IdealOri;";
	for ( uint32_t s = 0; s < loginfo_grainevo_cnt; s++ ) { log_grainevo_file << ";" << mygrainevolution[s].localSv; }
	log_grainevo_file << endl;

	//each grain one row, each step a column, this matrix transpose is practical because most often the number of logsteps << number of grains
	//however this is indeed a performance hit
	uint32_t ndisjoint_defgrains = mydefgpool.size();
	uint32_t ndisjoint_grains = ndisjoint_defgrains + myrxgpool.size();
	for ( uint32_t dgr = 0; dgr < ndisjoint_defgrains; dgr++ ) {
		uint32_t oriid = mydefgpool[dgr].caori;

		log_grainevo_file << dgr << ";;;;;" << myoripool[oriid].bunge1 << ";" << myoripool[oriid].bunge2 << ";" << myoripool[oriid].bunge3 << ";" << myoripool[oriid].closestideal; 

		//yeah, it hits caches badly I know, as matrix transposes do ...
		for ( uint32_t s = 0; s < loginfo_grainevo_cnt; s++) {
			log_grainevo_file << ";" << this->mygrainevolution[s].localdatabucket[dgr];
		}
		log_grainevo_file << endl;
	}

	log_grainevo_file << ";ix;iy;iz;boundarycontact;phi1;PHI;phi2;IdealOri";
	for ( uint32_t s = 0; s < loginfo_grainevo_cnt; s++) { log_grainevo_file << ";"; }
	log_grainevo_file << endl;

	for ( uint32_t rgr = ndisjoint_defgrains; rgr < ndisjoint_grains; rgr++ ) {
		uint32_t i = rgr - ndisjoint_defgrains;

		uint32_t z = myrxgpool[i].startingsite / myCAGeometry.nboxarea_rdtd;
		uint32_t rem = myrxgpool[i].startingsite - (myCAGeometry.nboxarea_rdtd * z);
		uint32_t y = rem / myCAGeometry.nboxedge_rd;
		uint32_t x = rem - (myCAGeometry.nboxedge_rd * y);

		uint32_t oriid = myrxgpool[i].caori;

		log_grainevo_file << i << ";" << x << ";" << y << ";" << z << ";" << "" << ";" << myoripool[oriid].bunge1 << ";" << myoripool[oriid].bunge2 << ";" << myoripool[oriid].bunge3 << ";" << myoripool[oriid].closestideal; //optionally add << "-RXG";
		//yeah, it hits caches badly I know, as matrix transposes do ...
		for ( uint32_t s = 0; s < loginfo_grainevo_cnt; s++) {
			log_grainevo_file << ";" << mygrainevolution[s].localdatabucket[rgr];
		}
		log_grainevo_file << endl;
	}

	log_grainevo_file.flush();
	log_grainevo_file.close();
}


void caHdl::write_grainevolution_binary( void )
{
	uint32_t ndisjoint_defgrains = mydefgpool.size();
	uint32_t ndisjoint_grains = ndisjoint_defgrains + myrxgpool.size(); //rxg grains written after the deformed ones

	//plain MPI I/O binary file structure
	uint32_t tx_nrows = loginfo_grainevo_cnt;
	uint32_t tx_ncols = 2; //time and X data pairs

	//allocate buffers for MPI I/O
	double* txbuffer = NULL;
	txbuffer = new double[tx_nrows*tx_ncols];
	QUICKASSERT( txbuffer != NULL );
	//fill txbuffer first all time then all X
	for (uint32_t ttt = 0; ttt < tx_nrows; ttt++ ) { txbuffer[ttt] =  mygrainevolution[ttt].localtime; }
	for (uint32_t xxx = 0; xxx < tx_nrows; xxx++ ) { txbuffer[tx_nrows+xxx] =  mygrainevolution[xxx].localX; }

	uint32_t dg_nrows = loginfo_grainevo_cnt;
	uint32_t rg_nrows = loginfo_grainevo_cnt;
	uint32_t dg_ncols = ndisjoint_defgrains;
	uint32_t rg_ncols = ndisjoint_grains - ndisjoint_defgrains;

	//lean buffer!
	uint32_t* dgbuffer = NULL;
	dgbuffer = new uint32_t[1*dg_nrows];
	QUICKASSERT( dgbuffer != NULL );

	uint32_t* rgbuffer = NULL;
	rgbuffer = new uint32_t[1*rg_nrows];
	QUICKASSERT( rgbuffer != NULL );

	//MPI file ...
	MPI_File msFileHdlTX;
	MPI_File msFileHdlDef;
	MPI_File msFileHdlRxg;
	MPI_Status msFileStatusTX;
	MPI_Status msFileStatusDef;
	MPI_Status msFileStatusRxg;

	// create C-consistent file name for MPI I/O
	stringstream msfntx;
	stringstream msfndg;
	stringstream msfnrg;
	msfntx << "SCORE." << myensHdl->simid << ".JobID." <<  this->jobid << ".txnrows." << tx_nrows << ".txncols." << tx_ncols << ".TimeLocalX.bin";
	msfndg << "SCORE." << myensHdl->simid << ".JobID." <<  this->jobid << ".dgnrows." << dg_nrows << ".dgncols." << dg_ncols << ".GrowthDefGrains.bin";
	msfnrg << "SCORE." << myensHdl->simid << ".JobID." << this->jobid << ".rgnrows." << rg_nrows << ".rgncols." << rg_ncols << ".GrowthRXGrains.bin";
	int msfntxl = msfntx.str().size();
	int msfndgl = msfndg.str().size();
	int msfnrgl = msfnrg.str().size();
	char* msfnttxx = new char[msfntxl+1];
	char* msfnd = new char[msfndgl+1];
	char* msfnr = new char[msfnrgl+1];
	strcpy(msfnttxx, msfntx.str().c_str());
	strcpy(msfnd, msfndg.str().c_str());
	strcpy(msfnr, msfnrg.str().c_str());

	//I output time-x data
	MPI_File_open(MPI_COMM_SELF, msfnttxx, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &msFileHdlTX);
	__int64 txoffset = 0;
	MPI_File_seek( msFileHdlTX, txoffset, MPI_SEEK_SET );
	MPI_File_write_at(msFileHdlTX, txoffset, txbuffer, tx_nrows * tx_ncols, MPI_DOUBLE, &msFileStatusTX);
	MPI_File_close(&msFileHdlTX);

	//now push growth curves into file
	MPI_File_open(MPI_COMM_SELF, msfnd, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &msFileHdlDef);
	MPI_File_open(MPI_COMM_SELF, msfnr, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &msFileHdlRxg);
	__int64 dgoffset = 0;
	__int64 rgoffset = 0;
	MPI_File_seek( msFileHdlDef, dgoffset, MPI_SEEK_SET ); //##MK::might as fp is at 0 after opening anyway be obsolete...
	MPI_File_seek( msFileHdlRxg, rgoffset, MPI_SEEK_SET );

	for ( uint32_t dg = 0; dg < dg_ncols; dg++ ) {
		for ( uint32_t tdd = 0; tdd < dg_nrows; tdd++ ) {
			dgbuffer[tdd] = mygrainevolution[tdd].localdatabucket[dg];
		}
		dgoffset = dg * dg_nrows * 4; //unsigned int is 4byte in MPI I/O

		MPI_File_write_at(msFileHdlDef, dgoffset, dgbuffer, dg_nrows, MPI_UNSIGNED, &msFileStatusDef);
	}

	for ( uint32_t rg = 0; rg < rg_ncols; rg++ ) {
		for ( uint32_t trr = 0; trr < rg_nrows; trr++ ) {
			rgbuffer[trr] = mygrainevolution[trr].localdatabucket[dg_ncols+rg];
		}
		rgoffset = rg * rg_nrows * 4;

		MPI_File_write_at(msFileHdlRxg, rgoffset, rgbuffer, rg_nrows, MPI_UNSIGNED, &msFileStatusRxg);
	}


	MPI_File_close(&msFileHdlDef); //no Barrier as MPI_COMM_SELF
	MPI_File_close(&msFileHdlRxg);

	delete [] msfnttxx;		msfnttxx = NULL;
	delete [] msfnd;		msfnd = NULL;
	delete [] msfnr;		msfnr = NULL;
	delete [] txbuffer;		txbuffer = NULL;
	delete [] dgbuffer;		dgbuffer = NULL;
	delete [] rgbuffer;		rgbuffer = NULL;

	cout << myensHdl->myRank << " successful export of growth curves, tx nrows/ncols = " << tx_nrows << ";" << tx_ncols << endl;
	cout << myensHdl->myRank << " successful export of growth curves, dg nrows/ncols = " << dg_nrows << ";" << dg_ncols << endl;
	cout << myensHdl->myRank << " successful export of growth curves, rg nrows/ncols = " << rg_nrows << ";" << rg_ncols << endl;

}


void caHdl::write_PercolationProfilingSummary( void )
{
	stringstream log_tprof_fname;
	ofstream log_tprof_file;

	log_tprof_fname << "SCORE." << myensHdl->simid << ".MyRank." << myensHdl->myRank << ".JobID." << this->jobid << ".PercolationProfiling.csv";
	log_tprof_file.open ( log_tprof_fname.str().c_str() );

	//header
	log_tprof_file << "Step;Time/s;X;NCluster;LargestCluster;Percolation;;Initializing;HoshenKopelmanLabeling;Compactifying;CheckLabeling;Characterizing;CheckPercolating;Total(all in seconds)";
	log_tprof_file << endl;

	for ( unsigned int s = 0; s < this->myhk.size(); s++ ) {
		log_tprof_file << this->regions[MASTER]->profiling_growthmachine[s].localstep << ";" << this->regions[MASTER]->profiling_growthmachine[s].localtime << ";" << this->regions[MASTER]->profiling_growthmachine[s].localX << ";";
		log_tprof_file << myhk[s].NCluster << ";" << myhk[s].LargestCluster << ";" << myhk[s].Percolation << ";;" << myhk[s].ProfInitializing << ";" << myhk[s].ProfHoshenKopeling << ";" << myhk[s].ProfCompactifying;
		log_tprof_file << ";" << myhk[s].ProfCheckLabeling << ";" << myhk[s].ProfCharacterizing << ";" << myhk[s].ProfCheckPercolating << ";"  << myhk[s].ProfPercTotal << "\n";
	}
	log_tprof_file << endl;

	log_tprof_file.flush();
	log_tprof_file.close();
}


void caHdl::write_rxareaprofile_ascii( void )
{
	//I/O
	stringstream csv_fnm;
	csv_fnm << "SCORE.SimID." << this->myensHdl->simid << ".RXAreaFractionNDDepthProfile.csv";

	ofstream csv;
	csv.open ( csv_fnm.str().c_str() );
	csv << "t";
	for( size_t i = 0; i < this->myrxprofile.size(); ++i ) { csv << ";" << this->myrxprofile[i].localtime; } csv << "\n";
	csv << "Xvolume";
	for( size_t i = 0; i < this->myrxprofile.size(); ++i ) { csv << ";" << this->myrxprofile[i].localX; } csv << "\n";
	csv << "Step";
	for( size_t i = 0; i < this->myrxprofile.size(); ++i ) { csv << ";" << this->myrxprofile[i].localstep; } csv << "\n";
	csv << "ZCoordinate" << "\n";

	uint32_t imgz = myCAGeometry.nboxedge_nd;
	for( uint32_t z = 0; z < imgz; z++ ) {
		csv << z;
		for( size_t i = 0; i < this->myrxprofile.size(); ++i ) {
			if ( this->myrxprofile[i].localdatabucket != NULL ) {
				uint32_t rx = this->myrxprofile[i].localdatabucket[2*z+0];
				uint32_t dg = this->myrxprofile[i].localdatabucket[2*z+1];
				if ( (rx+dg) > 0 )
					csv << ";" << (static_cast<double>(rx) / static_cast<double>((rx+dg)));
				else
					csv << ";" << -1.0;
			}
			else
				csv << ";" << -1.0;
		}
		csv << "\n";
	}

	csv.flush();
	csv.close();
}




#define MPIIO_STRIPESIZE	((1024*1024)*10)

void caHdl::ompcrit_write_voxeldata_coloring_regions( void )
{
	//region topology
	uint32_t npx = this->regions[MASTER]->thePartitioning.nreg_rdx; //1 required
	if ( npx != 1 ) { cout << "ERROR::Plotting of region assignment implemented currently only for npx == 1!" << endl; return;}
	uint32_t npy = this->regions[MASTER]->thePartitioning.nreg_tdy;
	uint32_t npz = this->regions[MASTER]->thePartitioning.nreg_ndz;

	// create C-consistent file name for MPI I/O
	stringstream msFileName;
	msFileName << "SCORE." << myensHdl->simid << ".JobID." <<  this->jobid << ".Iter." << this->step << ".X.";
	uint32_t XX = (X * 1000.0);
	if ( XX < 10 )					msFileName << "000" << XX;
	if ( XX >= 10 && XX < 100 )		msFileName << "00" << XX;
	if ( XX >= 100 && XX < 1000 )	msFileName << "0" << XX;
	if ( XX >= 1000 )				msFileName << XX;
	msFileName << ".MS3DRegionID.raw";

	uint32_t msFileNameLength = msFileName.str().size();
	char* CmsFileName = new char[msFileNameLength+1];
	strcpy(CmsFileName, msFileName.str().c_str());

	//only the MASTER node writes in parallel! in write-only mode
	MPI_File msFileHdl;
	MPI_Status msFileStatus;
	MPI_File_open(MPI_COMM_SELF, CmsFileName, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &msFileHdl);
	//MK:: with 64-bit architectures, offsets calculations can be performed as longs 2^64 bit thus on petabyte files...
	//in mpi.h MPI_Offset is defined as an __int64 which is long long, thus we can jump much more than 2^32 directly when unsigned int would be utilized
	__int64 totalOffset = 0;
	MPI_File_seek( msFileHdl, totalOffset, MPI_SEEK_SET );

	//single automaton smaller than uint32_t range on shared memory in one process!
	//uint32_t nx = myCAGeometry.nboxedge_rd;
	uint32_t nxy = myCAGeometry.nboxarea_rdtd;
	//one temporary container of size myCAGeometry.nboxvol_rdtdnd might not be allocateable
	//uint32_t zbuffer = MPIIO_STRIPESIZE / 4 / nxy;
	//if ( zbuffer < 1 ) zbuffer = 1;

	uint32_t* rawdata = NULL;
	rawdata = new uint32_t[nxy];
	QUICKASSERT ( rawdata != NULL );

	//the automaton comprises npz regions along z, x is aligned, and y as well, so simply read x*y slab stack in npy direction and drop such a global z section into rawdata
	//##MK::at the moment to find safe but less efficient sequential working through xy layers along z
	uint32_t c, r;
	uint32_t zper, yper, xper;

	for (uint32_t zregion = 0; zregion < npz; zregion++ ) { //utilize that regions form regular stacked aggregate
		zper = this->regions[(zregion*npy)+0]->myGeom.nreg_nd; //thickness in z direction, for each layer collect regions stacked in y...

		for ( uint32_t zz = 0; zz < zper; zz++ ) { //layer-by-layer MPI I/O several MB per write
			c = 0;

			//collect data over regions
			for(uint32_t yregion = 0; yregion < npy; yregion++ ) {
				r = (zregion*npy)+yregion;
				xper = this->regions[r]->myGeom.nreg_rd; //exclusive global coordinate
				yper = this->regions[r]->myGeom.nreg_td; //##MK::optimizable...

				for ( uint32_t yy = 0; yy < yper; yy++ ) {
					for ( uint32_t xx = 0; xx < xper; xx++ ) {
						//corr = xx + (xper * yy) + (xper * yper * zz);

						//global coordinate
						if ( xx > 0 && xx < (xper-2) && yy > 0 && yy < (yper-2) && zz > 0 && zz < (zper-2) ) { //most likely case the inside volume of the memory region
							rawdata[c] = 1 * (r+1);
						}
						else {// in the border of the memory region
							rawdata[c] = 1000 * (r+1);
						}
						c++;
					} //scan through xline
				} //scan +y the xlines stacked in region r

				//sum over yper adds to size of automaton in y

			} //next region stacked upon last one in +y

			//push in file at once, runtime environment will optimize, further potential by nowait multiple write from the thread, but tricky...
			//MPI_File_write_at(msFileHdl, totalOffset, rawdata, nxyz, MPI_INT, &msFileStatus);
			MPI_File_write(msFileHdl, rawdata, nxy, MPI_UNSIGNED, &msFileStatus); //implicit advancement of fp

		} //next zslice in the zregion
	} //next zregion slab of regions stacked in +y

	delete [] rawdata;
	delete [] CmsFileName;

//cout << myensHdl->myRank << " successful MPI I/O in timestep = " << this->step << endl;

	MPI_File_close(&msFileHdl); //no Barrier as MPI_COMM_SELF
}


void caHdl::ompcrit_write_voxeldata_coloring_grainids( string postfix )
{
	//TARGET FORMAT IS UINT32 FOR PARAVIEW
	//CELL_IS_PARTICLE	0
	//CELL_IS_INFECTED	1
	//DEFORMED 			2+mydefgid //because there is defgid = 0
	//RECRYSTALLIZED	3+mydefgid.size()+myrxgid

	//region topology
	uint32_t npx = this->regions[MASTER]->thePartitioning.nreg_rdx; //1 required
	uint32_t npy = this->regions[MASTER]->thePartitioning.nreg_tdy;
	uint32_t npz = this->regions[MASTER]->thePartitioning.nreg_ndz;
	if ( npx != 1 ) { std::cout << "ERROR::No 3D output because npx != 1!" << std::endl; }

	MPI_File msFileHdl;
	MPI_Status msFileStatus;
	// create C-consistent file name for MPI I/O
	stringstream msFileName;
	msFileName << "SCORE." << myensHdl->simid << ".JobID." <<  this->jobid << ".Iter." << this->step << ".X.";
	uint32_t XX = (X * 1000.0);
	if ( XX < 10 )					msFileName << "000" << XX;
	if ( XX >= 10 && XX < 100 )		msFileName << "00" << XX;
	if ( XX >= 100 && XX < 1000 )	msFileName << "0" << XX;
	if ( XX >= 1000 )				msFileName << XX;
	msFileName << "." << postfix << ".raw";

	uint32_t msFileNameLength = msFileName.str().size();
	char* CmsFileName = new char[msFileNameLength+1];
	strcpy(CmsFileName, msFileName.str().c_str());

	//in mpi.h MPI_Offset is defined as an __int64 which is long long, thus we can jump much more than 2^32 directly when unsigned int would be utilized
	MPI_File_open( MPI_COMM_SELF, CmsFileName, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &msFileHdl );
	__int64 totalOffset = 0;
	MPI_File_seek( msFileHdl, totalOffset, MPI_SEEK_SET );

	//one temporary container of size myCAGeometry.nboxvol_rdtdnd may not be allocateable or at least significantly fragment memory
	uint32_t nxy = myCAGeometry.nboxarea_rdtd;
	uint32_t* rawdata = NULL;
	rawdata = new uint32_t[nxy];
	QUICKASSERT ( rawdata != NULL );

	//the automaton comprises npz regions along z, x is aligned, and y as well, so simply read x*y slab stack in npy direction and drop such a global z section into rawdata
	//##MK::at the moment a safe but less efficient than an optimized buffering strategy of working sequentially through xy layers along z is performed
	uint32_t c, r, corr, gridvalue;
	uint32_t yper, xper;
	for (uint32_t zregion = 0; zregion < npz; zregion++ ) { //utilize that the regions form regular-stacked 3D region aggregate in z direction
		for ( uint32_t zz = 0; zz < this->regions[(zregion*npy)+0]->myGeom.nreg_nd; zz++ ) { //layer-by-layer MPI I/O several MB per write
			c = 0;
			//collect data over regions in a fixed z layer advancing in +z layer stacking direction
			for(uint32_t yregion = 0; yregion < npy; yregion++ ) {
				r = (zregion*npy) + yregion;
				yper = this->regions[r]->myGeom.nreg_td; //nreg_td is an extension of the region box along the y direction counting in cells ##MK::old left-handed csys yregion*npz
				xper = this->regions[r]->myGeom.nreg_rd;

				for ( uint32_t yy = 0; yy < yper; yy++ ) {
					for ( uint32_t xx = 0; xx < xper; xx++ ) {
						corr = xx + (xper * yy) + (xper * yper * zz); //implicit coordinate in the region

						gridvalue = regions[r]->mycellgrid[corr];

						rawdata[c] = INVALID_CELLASSIGNMENT;		//4096000003

						if ( gridvalue == CELL_IS_A_PARTICLE ) {	//0
							rawdata[c] = 0; //RAW_PARTICLE;
							c++;
							continue;
						}
						if ( gridvalue == CURRENTLY_INFECTED ) {	//1
							rawdata[c] = 1; //RAW_INFECTED;
							c++;
							continue;
						}
						//other cases already excluded
						if ( gridvalue <= CA_GRAINIDS_MAXIMUM ) {	//<=4096000002
							rawdata[c] = 2 + gridvalue;
							c++;
							continue;
						}
					} //scan through xline
				} //scan +y the xlines stacked in region r
			} //next region stacked upon last one in +y

			//push layer to file at once, runtime environment will optimize, further potential by nowait multiple write from the thread, but tricky...
			MPI_File_write(msFileHdl, rawdata, nxy, MPI_UNSIGNED, &msFileStatus); //implicit advancement of fp
		} //next zslice in zregions
	} //next region z with regions on stacked top of one another in y

	delete [] rawdata;
	delete [] CmsFileName;

//std::cout << myensHdl->myRank << " successful MPI I/O in timestep = " << this->step << std::endl;
	MPI_File_close(&msFileHdl); //no Barrier as MPI_COMM_SELF
}


void caHdl::ompcrit_write_zsection_coloring_ipfz( double zpos, double xnow, char mode )
{
	//master executes this has to pull from neighboring threads memory
	//render a xy section as a PNG
	uint32_t imgx = myCAGeometry.nboxedge_rd;
	uint32_t imgy = myCAGeometry.nboxedge_td;
	uint32_t imgxy = imgx * imgy;
	uint32_t imgz = (uint32_t) ( zpos * ((double) myCAGeometry.nboxedge_nd));
	if ( imgz == myCAGeometry.nboxedge_nd ) //MK::global coordinates [0,nboxedge_i) along each direction!
		imgz--;

	//in which regions are the image data located? get region topology
	//uint32_t nregions = this->regions.size();
	uint32_t npx = this->regions[MASTER]->thePartitioning.nreg_rdx; //1 required, no domain split in x
	uint32_t npy = this->regions[MASTER]->thePartitioning.nreg_tdy;
	uint32_t npz = this->regions[MASTER]->thePartitioning.nreg_ndz;
	if ( npx != 1 ) { cout << "ERROR::Not outputting sections because npx != 1!" << endl; }

	//find the zlayer by scanning along the first npz region limits
	uint32_t zregion = 0;
	uint32_t zmi, zmx, zz;
	bool found = false;
	for ( zregion = 0; zregion < npz; zregion++ ) {
		zmi = this->regions[(zregion*npy)+0]->myGeom.nreg_ndmin;
		zmx = this->regions[(zregion*npy)+0]->myGeom.nreg_ndmax;

		if ( imgz >= zmi && imgz <= zmx ) {
			zz = imgz - zmi; //local grid coordinate z for all regions that include the z position imgz adnd participate in I/O
			found = true;
			break;
		}
	}

	//now we know in which region we are and therefore which thread is responsible

	//allocate memory to store the resulting image
	unsigned char* img_rgba = NULL;
	ostringstream fname;

	if ( found == true ) {
		img_rgba = new unsigned char [4*imgx*imgy];

		//mode == COLORIZE_RX_LEAVE_DEFORMED_BLACK prepare image by filling in black for CELL_IS_DEFORMED or CURRENTLY_INFECTED
		//mode == COLORIZE_DEF_LEAVE_RX_BLACK prepare image by filling in black for everything which has no grain ID
		uint32_t cx;
		for ( uint32_t y = 0; y < imgy; y++ ) {
			for ( uint32_t x = 0; x < imgx; x++ ) {
				cx = 4*(x + (y*imgx));
				img_rgba[cx + REDCHAN ] = BLACK;
				img_rgba[cx + GREENCHAN ] = BLACK;
				img_rgba[cx + BLUECHAN ] = BLACK;
				img_rgba[cx + ALPHACHAN ] = UCHAR_RANGE_MAX;
			}
		}

		//the automaton comprises npz regions along z, x is aligned, and y as well, so simply read x*y slab stack in npy direction and drop such a global z section into the image
		uint32_t xper, yper;
		uint32_t r,corr, gridvalue, oid;
		uint32_t ndg = this->mydefgpool.size();
		QUICKASSERT ( ndg == this->nmydefgpool );

		uint32_t c = 0;
		for( uint32_t yregion = 0; yregion < npy; yregion++ ) { //collect data over regions
			r = (zregion*npy)+yregion;
			xper = this->regions[r]->myGeom.nreg_rd; //MK::extent in the r-th region not global SU coordinate
			yper = this->regions[r]->myGeom.nreg_td;

			//colorizing, but leaving first CURRENTLY_INFECTED cells untouched
			for ( uint32_t yy = 0; yy < yper; yy++ ) {
				for ( uint32_t xx = 0; xx < xper; xx++ ) {
					corr = xx + (xper * yy) + (xper * yper * zz);

					gridvalue = regions[r]->mycellgrid[corr];

					//CELL_IS_A_PARTICLE, CURRENTLY_INFECTED
					if ( mode == COLORIZE_RX_LEAVE_DEFORMED_BLACK ) { //also leaving currently infected black
						if ( gridvalue <= CA_GRAINIDS_MAXIMUM ) { //mind that CURRENTLY_INFECTED > CA_GRAINIDS_MAXIMUM most likely
							if ( gridvalue >= ndg ) { //fully recrystallized because by definition CURRENTLY_INFECTED > CA_GRAINIDS_MAXIMUM
								oid = myrxgpool[gridvalue-ndg].caori;
								img_rgba[4*c + REDCHAN ] = myoripool[oid].RGB_R;
								img_rgba[4*c + GREENCHAN ] = myoripool[oid].RGB_G;
								img_rgba[4*c + BLUECHAN ] = myoripool[oid].RGB_B;
								img_rgba[4*c + ALPHACHAN ] = myoripool[oid].RGB_A;
							} //leave gridvalue < ndg BLACK
						}
						/*if ( gridvalue == CURRENTLY_INFECTED ) {
							//read cellStatii and update later for now leave BLACK!
						}*/
						//else leave definately BLACK
					}
					else if ( mode == COLORIZE_DEF_LEAVE_RX_BLACK ) {
						if ( gridvalue <= CA_GRAINIDS_MAXIMUM ) { //most likely
							if ( gridvalue < ndg ) { //thus implicitly smaller than CA_GRAINIDS_MAXIMUM, also leaving CURRENTLY_INFECTED black
								oid = mydefgpool[gridvalue].caori;
								img_rgba[4*c + REDCHAN ] = myoripool[oid].RGB_R;
								img_rgba[4*c + GREENCHAN ] = myoripool[oid].RGB_G;
								img_rgba[4*c + BLUECHAN ] = myoripool[oid].RGB_B;
								img_rgba[4*c + ALPHACHAN ] = myoripool[oid].RGB_A;
							} //leave RX BLACK
						}
						/*if ( gridvalue == CURRENTLY_INFECTED ) {
							//read cellStatii and update later for now leave BLACK!
						}*/
						//else leave definately BLACK
					}
					else { //COLORIZE_RX_RHOGREYSCALE_DEFORMED
						if ( gridvalue <= CA_GRAINIDS_MAXIMUM ) {
							if ( gridvalue < ndg ) { //deformed ndg < CURRENTLY_INFECTED therefore implicitly also testing therewith != CURRENTLY_INFECTED!
								unsigned char alpha = rho2grayscale( mydefgpool[gridvalue].rho );
								img_rgba[4*c + REDCHAN ] = alpha; //setting all three RGB to the same value renders grey
								img_rgba[4*c + GREENCHAN ] = alpha;
								img_rgba[4*c + BLUECHAN ] = alpha;
								img_rgba[4*c + ALPHACHAN ] = UCHAR_RANGE_MAX;
							}
							else { //>= ndg
								//else if ( gridvalue <= CA_GRAINIDS_MAXIMUM ) { //implicitly >= ndg and explicitly <= CA_GRAINIDS_MAXIMUM < CURRENTLY_INFECTED by definition! hence implicitly RX because && >= ndg
								/*img_rgba[4*c + REDCHAN ] = UCHAR_RANGE_MIN;
								img_rgba[4*c + GREENCHAN ] = UCHAR_RANGE_MAX; //MK::DEBUG mask RX grains as pure green
								img_rgba[4*c + BLUECHAN ] = UCHAR_RANGE_MIN;
								img_rgba[4*c + ALPHACHAN ] = UCHAR_RANGE_MAX;*/
								//#######changing masking and instead IPFZ colorize RXGrains as in COLORIZE_RX
								oid = myrxgpool[gridvalue-ndg].caori;
								img_rgba[4*c + REDCHAN ] = myoripool[oid].RGB_R;
								img_rgba[4*c + GREENCHAN ] = myoripool[oid].RGB_G;
								img_rgba[4*c + BLUECHAN ] = myoripool[oid].RGB_B;
								img_rgba[4*c + ALPHACHAN ] = myoripool[oid].RGB_A;
							}
						}
						/*else { colorize > CA_GRAINIDS_MAXIMUM,i.e. implicitly == CURRENTLY_INFECTED later, i.e. infected, thus leave initial RGB assignment and explicitly set ALPHACHAN again
							img_rgba[4*c + REDCHAN ] = UCHAR_RANGE_MAX; //mask currently infected cells as pure red
							img_rgba[4*c + GREENCHAN ] = UCHAR_RANGE_MIN;
							img_rgba[4*c + BLUECHAN ] = UCHAR_RANGE_MIN;
							img_rgba[4*c + ALPHACHAN] = UCHAR_RANGE_MAX;
						}*/
					}

					c++; //stitching of the image occurs implicitly as yregion is increased

				} //scan through xline in region r
			} //scan +y the xlines stacked in region r

			//recolorize CURRENTLY_INFECTED cells correctly
			if ( this->cellStatii.at(r) != NULL ) { //if a clarification of the active cell status was done...
				uint32_t ncandidates = this->cellStatii.at(r)->size();
				vector<cellstatus>* candidates = this->cellStatii.at(r);

				for ( uint32_t cand = 0; cand < ncandidates; cand++ ) {
					//coordinate transform starting with extracting global z variable
					uint32_t z = candidates->at(cand).ixyz / imgxy;
					//already now we know whether we are in the correct layer...
					if ( z != zpos ) { //most likely we're not...
						continue;
					}
					//we are in the correct layer, so recolorize all pixel at which initially we knew only the cellstate == CURRENTLY_INFECTED
					uint32_t rem = candidates->at(cand).ixyz - (imgxy * z);
					uint32_t y = rem / imgx;
					uint32_t x = rem - (imgx * y); //x,y,z are SU global coordinates and thus translate into PNG image coordinates
					cx = 4*(x + (y*imgx));
					gridvalue = candidates->at(cand).gid;

					if ( mode == COLORIZE_RX_LEAVE_DEFORMED_BLACK ) {
						if ( gridvalue >= ndg ) { //remind that we explicitly handle here only cells in a state == CRRENTLY_INFECTED but rxFraction >= RXFRACTION_THRESHOLD
							oid = myrxgpool[gridvalue-ndg].caori;
							img_rgba[cx + REDCHAN ] = myoripool[oid].RGB_R;
							img_rgba[cx + GREENCHAN ] = myoripool[oid].RGB_G;
							img_rgba[cx + BLUECHAN ] = myoripool[oid].RGB_B;
							img_rgba[cx + ALPHACHAN ] = myoripool[oid].RGB_A;
						} //else it is deformed so leave BLACK
					}
					else if ( mode == COLORIZE_DEF_LEAVE_RX_BLACK ) {
						if ( gridvalue < ndg ) {
							oid = mydefgpool[gridvalue].caori;
							img_rgba[cx + REDCHAN ] = myoripool[oid].RGB_R;
							img_rgba[cx + GREENCHAN ] = myoripool[oid].RGB_G;
							img_rgba[cx + BLUECHAN ] = myoripool[oid].RGB_B;
							img_rgba[cx + ALPHACHAN ] = myoripool[oid].RGB_A;
						} //else it is RX so leave BLACK
					}
					else { //COLORIZE_RX_RHOGREYSCALE_DEFORMED )
						if ( gridvalue < ndg ) { //deformed and == CURRENTLY_INFECTED
							unsigned char alpha = rho2grayscale( mydefgpool[gridvalue].rho );
							img_rgba[cx + REDCHAN ] = alpha; //setting all three RGB to the same value renders grey
							img_rgba[cx + GREENCHAN ] = alpha;
							img_rgba[cx + BLUECHAN ] = alpha;
							img_rgba[cx + ALPHACHAN ] = UCHAR_RANGE_MAX;
						}
						else {
							/*img_rgba[c + REDCHAN ] = UCHAR_RANGE_MIN;
							img_rgba[c + GREENCHAN ] = UCHAR_RANGE_MAX; //MK::DEBUG mask RX grains as pure green
							img_rgba[c + BLUECHAN ] = UCHAR_RANGE_MIN;
							img_rgba[c + ALPHACHAN ] = UCHAR_RANGE_MAX;*/
							//#######changing masking and instead IPFZ colorize RXGrains as in COLORIZE_RX
							oid = myrxgpool[gridvalue-ndg].caori;
							img_rgba[cx + REDCHAN ] = myoripool[oid].RGB_R;
							img_rgba[cx + GREENCHAN ] = myoripool[oid].RGB_G;
							img_rgba[cx + BLUECHAN ] = myoripool[oid].RGB_B;
							img_rgba[cx + ALPHACHAN ] = myoripool[oid].RGB_A;
						}
					}
				} //next cell currently in state == CURRENTLY_INFECTED to potentially recolorize
			} //done with potential recoloring of CURRENTLY_INFECTED cells
		} //next region stacked upon last one in +y
		//done with all regions along y who include the zsection at global SU coordinate imgz

		//get folder
		string wrkdir;
		if ( mode == COLORIZE_RX_LEAVE_DEFORMED_BLACK ) { wrkdir = this->myensHdl->simwrkdir + "/png/rxipfz_dgblack/"; }
		else if ( mode == COLORIZE_DEF_LEAVE_RX_BLACK ) { wrkdir = this->myensHdl->simwrkdir + "/png/dgipfz_rxblack/"; }
		else { wrkdir = this->myensHdl->simwrkdir + "/png/dgrhogry_rxipfz/"; }

		//get file name
		fname << wrkdir << "SCORE." << this->myensHdl->simid << ".JobID." << this->jobid << ".XYIter." << this->step << ".X.";
		//if ( mode != COLORIZE_DEF_LEAVE_RX_BLACK ) {
		uint32_t XX = (xnow * 1000.0);
		if ( XX < 10 )					fname << "000" << XX;
		if ( XX >= 10 && XX < 100 )		fname << "00" << XX;
		if ( XX >= 100 && XX < 1000 )	fname << "0" << XX;
		if ( XX >= 1000 )				fname << XX;
		//}		else { //COLORIZE_DEF_LEAVE_RX_BLACk	fname << "0000";		}

		fname << ".ZPOS." << imgz;
		if ( mode == COLORIZE_RX_LEAVE_DEFORMED_BLACK ) fname << ".rxIPFZdgBLACK.png";
		else if ( mode == COLORIZE_DEF_LEAVE_RX_BLACK ) fname << ".dgIPFZrxBLACK.png";
		//else fname << ".dgRHOGRYrxGREENinfRED.png"; //GREYSCALE_RHO_LEAVE_RX_WHITE
		else fname << ".dgRHOGRYrxIPFZ.png";

		//encode the PNG file by calling the lodePNG library
		lodepng::encode( fname.str().c_str(), img_rgba , imgx, imgy );

		delete [] img_rgba;
		img_rgba = NULL;

		cout << myensHdl->myRank << " rendered " << fname.str().c_str() << endl;
	}
	else {
		cout << myensHdl->myRank << " warns::RENDERING of zsection IPFZ was UNSUCCESSFUL! " << this->step << endl;
	}
}


void caHdl::ompcrit_write_zsection_looping_ipfz( double rxfraction, char thatmode ) {
	for ( uint32_t zr = 0; zr < rendering_atthisZPos.size(); zr++ ) {
		ompcrit_write_zsection_coloring_ipfz( rendering_atthisZPos.at(zr), rxfraction, thatmode );
	}
}


void caHdl::ompcrit_clarify_status_infectedcells( double rxthreshold ) {
	//MAY ONLY BE CALLED BY ONE THREAD!
	//MK::SCORE model is implemented such that cell carries a) index of deformed grain, b) flag index CURRENTLY_INFECTED showing infected, c) index of rxgrain
	//thus, once a cell is infected only the myRXFrontLists know the initial and future assignment, therefore we require a function
	//which collects for all such cells the current grain assignment
	//most naively one could for each pixel to print when hitting a CURRENTLY_INFECTED scan all lists in all memregions to find that particular cell carrying the information
	//however this is O(N^2) complexity, better, compile and cache the assignments once and thereafter operate over all these cells to 
	//modify within the image all those pixels at which an infection happened, resulting in O(N) complexity
	//if ( omp_get_thread_num() != MASTER ) 
	//	return;

	//first of all delete all old region entries, because the structure evolves, hence the assignments change
	uint32_t nreg = this->regions.size();
	for ( uint32_t r = 0; r < nreg; r++ ) {
		if ( cellStatii.at(r) != NULL ) {
			delete cellStatii[r];
			cellStatii[r] = NULL;
		}
	}

	//now interpret the assignments for each region; master analyzes regions sequentially could also be parallelized ##MK
	vector<cellstatus>* rstatii = NULL;
	struct cellstatus ac;
	uint32_t ndefg = this->mydefgpool.size();
	uint32_t nbx = myCAGeometry.nboxedge_rd; //avoid the cache collision when calling this in the for loop to select the points
	uint32_t nby = myCAGeometry.nboxedge_td;
	uint32_t nbxy = nbx * nby;
	uint32_t ncandidates = 0;
	cellP front = NULL;

	for ( uint32_t r = 0; r < nreg; r++ ) {
		rstatii = NULL;
		rstatii = new vector<cellstatus>;
		QUICKASSERT( rstatii != NULL );

		//analyze cells inside
		ncandidates = this->regions[r]->nextSlotNeverActiveRXFrontInside;
		front = this->regions[r]->myRXFrontInside;
		for ( uint32_t c = 0; c < ncandidates; c++ ) {
			if ( front[c].activity == ACTIVE ) { //most likely case for densely populated list
				ac.ixyz = front[c].ix + (front[c].iy * nbx) + (front[c].iz * nbxy); //ix, ... are short but positive! so short to uint32_t no problem
				if ( front[c].rxFrac < rxthreshold ) {
					ac.gid = front[c].mydefgid;
					rstatii->push_back( ac );
					continue;
				}
				//else, rxed
				ac.gid = ndefg + front[c].myrxgid;
				rstatii->push_back( ac );
			}
		}

		//analyze cells border
		ncandidates = this->regions[r]->nextSlotNeverActiveRXFrontBorder;
		front = this->regions[r]->myRXFrontBorder;
		for ( uint32_t c = 0; c < ncandidates; c++ ) {
			if ( front[c].activity == ACTIVE ) { //most likely case for densely populated list
				ac.ixyz = front[c].ix + (front[c].iy * nbx) + (front[c].iz * nbxy); //ix, ... are short but positive! so short to uint32_t no problem
				if ( front[c].rxFrac < rxthreshold ) {
					ac.gid = front[c].mydefgid;
					rstatii->push_back( ac );
					continue;
				}
				//else, rxed
				ac.gid = ndefg + front[c].myrxgid;
				rstatii->push_back( ac );
			}
		}
		
		//hand over pointer to the memory which stores these assignments
		this->cellStatii.at(r) = rstatii;
	} //next region
}


struct sembuffer
{
	double phi1;
	double Phi;
	double phi2;
	unsigned int gid;
	bool rxed;
	bool pad1;
	bool pad2;
	bool pad3;
	sembuffer() : phi1(-4.0*_PI_), Phi(-4.0*_PI_), phi2(-4.0*_PI_),
			gid(NO_GRAIN_ASSIGNED), rxed(false),
			pad1(false), pad2(false), pad3(false) {}
};

void caHdl::ompcrit_write_semebsd_core( double zpos, char mode ) {
	//CALLED FROM WITHIN PARALLEL REGION, but only master executes fishing in other threads memory
	//output ASCII texture file

	//omp_set_lock(&h5lock);

	uint32_t imgx = myCAGeometry.nboxedge_rd;
	uint32_t imgy = myCAGeometry.nboxedge_td;
	uint32_t imgxy = imgx*imgy;
	uint32_t imgz = (uint32_t) ( zpos * ((double) myCAGeometry.nboxedge_nd));
	if ( imgz == myCAGeometry.nboxedge_nd ) //MK::global coordinates [0,nboxedge_i) along each direction!
		imgz--;

	//strategy when piping sequentially point by point into the file we have to scan the list of cellstatii whenever the was in a CURRENTLY_INFECTED state
	//instead we could build a fixed size array of double triplets, flooding with all values we know and continuing with flooding all cellStatii, then writing
	//comparison shows that algorithm is O(n) in cell != INFECTED and O(n^2) in cell == INFECTED surplus O(n) in pixel for output, the O(n^2) we can reduce to a O(n)
	//we chose the latter option as it eventually will also redue cache debris because of one after another calling CA functions than stream operations at the cost of some memory
	//however even the largest CA 2d mapping will take only 3*8B*SQR(CA_DIMENSIONS_MAXIMUM) <= 60MB

	double* iobuffer = NULL;
	iobuffer = new double[4*imgxy]; //+x then stack xlines along +y
	QUICKASSERT( iobuffer != NULL );
	for ( uint32_t px = 0; px < imgxy; px++ ) {
		iobuffer[4*px+0] = -4.0*_PI_; //dummy value to detect error
		iobuffer[4*px+1] = -4.0*_PI_;
		iobuffer[4*px+2] = -4.0*_PI_;
		iobuffer[4*px+3] = 0; //##MK::for Martin Diehl we use Fortran indexing no grain ID should be zero
	}

	//in which regions are the image data located? get region topology
	//uint32_t nregions = this->regions.size();
	uint32_t npx = this->regions[MASTER]->thePartitioning.nreg_rdx; //1 required, no domain split in x
	uint32_t npy = this->regions[MASTER]->thePartitioning.nreg_tdy;
	uint32_t npz = this->regions[MASTER]->thePartitioning.nreg_ndz;
	if ( npx != 1 ) { cout << "ERROR::Not outputting sections because npx != 1!" << endl; }

	//find the zlayer by scanning along the first npz region limits
	uint32_t zregion = 0;
	uint32_t zmi, zmx, zz;
	bool found = false;
	for ( zregion = 0; zregion < npz; zregion++ ) {
		zmi = this->regions[(zregion*npy)+0]->myGeom.nreg_ndmin; //global coordinates
		zmx = this->regions[(zregion*npy)+0]->myGeom.nreg_ndmax;

		if ( imgz >= zmi && imgz <= zmx ) {
			zz = imgz - zmi; //local grid coordinate z for all regions that include the z position imgz adnd participate in I/O
			found = true;
			break;
		}
	}


	//###need to implement scanning if sufficient size grain and in case add DIEHL_RXG_OFFSET

	stringstream ang_fname;
	ang_fname << this->myensHdl->simwrkdir << "/ang/" << "SCORE." << this->myensHdl->simid << ".JobID." << this->jobid << ".XYIter." << this->step << ".X.";
	uint32_t XX = (X * 1000.0);
	if ( XX < 10 )					ang_fname << "000" << XX;
	if ( XX >= 10 && XX < 100 )		ang_fname << "00" << XX;
	if ( XX >= 100 && XX < 1000 )	ang_fname << "0" << XX;
	if ( XX >= 1000 )				ang_fname << XX;
	ang_fname << ".ZPOS." << imgz << ".ang";

	ofstream ang;
	ang.open ( ang_fname.str().c_str() );

	//##MK::if really desiring ANG file add specific header, ##MK::currently implementing only dummy
	ang << "//SCORE auto-generated artificial SEM/EBSD file!\n";
	ang << "//SCORE coordinate system x || RD, y || TD, z || ND, right-handed\n";
	ang << "\n\n\n";
	ang << "BungePhi1;BungePHI;BungePhi2;FortranIndexedGrainID(PosForRXed,NegForDef)\n";

	if ( found == true ) {
		//the automaton comprises npz regions along z, x is aligned, and y as well, so simply read x*y slab stack in npy direction and drop such a global z section into the image
		uint32_t xper, yper;
		uint32_t r,region_corr, su_corr, gridvalue, oid;
		uint32_t ndg = this->mydefgpool.size();
		QUICKASSERT ( ndg == this->nmydefgpool );

		uint32_t xmi,xmx,ymi,ymx;
		for( uint32_t yregion = 0; yregion < npy; yregion++ ) { //collect data over regions
			r = (zregion*npy)+yregion;
			xper = this->regions[r]->myGeom.nreg_rd; //MK::extent in the r-th region not global SU coordinate
			yper = this->regions[r]->myGeom.nreg_td;
			xmi = this->regions[r]->myGeom.nreg_rdmin; //global coordinates delimitng memory region inclusively
			xmx = this->regions[r]->myGeom.nreg_rdmax;
			ymi = this->regions[r]->myGeom.nreg_tdmin;
			ymx = this->regions[r]->myGeom.nreg_tdmax;

			for ( uint32_t yy = 0; yy < yper; yy++ ) {
				uint32_t regyzoff = (xper * yy) + (xper * yper * zz);
				uint32_t yzoff = (ymi + yy)*imgx;
				for ( uint32_t xx = 0; xx < xper; xx++ ) {
					region_corr = xx + regyzoff;
					su_corr = xmi + xx + yzoff;

					gridvalue = regions[r]->mycellgrid[region_corr];

					if ( gridvalue <= CA_GRAINIDS_MAXIMUM ) { //therefore also != CURENTLY_INFECTED most likely either rx or defg
						if ( gridvalue >= ndg ) {
							oid = myrxgpool[gridvalue-ndg].caori;
							iobuffer[4*su_corr+3] = +1.0 * (gridvalue - ndg + 1); //##MK::for MartinDiehl Fortran indexed ID of RX grain
						}
						else {
							oid = mydefgpool[gridvalue].caori;
							iobuffer[4*su_corr+3] = -1.0 * (gridvalue + 1); //##MK::for MartinDiehl Fortran indexed ID of deformed grain
						}
						iobuffer[4*su_corr+0] = myoripool[oid].bunge1;
						iobuffer[4*su_corr+1] = myoripool[oid].bunge2;
						iobuffer[4*su_corr+2] = myoripool[oid].bunge3;

						//ang << setprecision(8) << myoripool[oid].bunge1 << ";" << setprecision(8) << myoripool[oid].bunge2 << ";" << setprecision(8) << myoripool[oid].bunge3 << endl;
					}
					/*
					else { //find what is at this location, MK::requires scanning the list of active cells begin with Inside if nothing there move to border
						//will be filled in later...
						//currently naive solution, scan entire list of active cells
						//get global coordinate
						su_corr = (xmi + xx) + ((ymi + yy)*imgx) + (imgz*(imgx*imgy)); //utilizing that imgx and imgy are myCAGeometry.nedge_rd, ... nedge_td and imgz a global coordinate
						uint32_t ncandidates = this->cellStatii.at(r)->size();
						vector<cellstatus>* candidates = this->cellStatii.at(r);
						bool match = false;

						for ( uint32_t cand = 0; cand < ncandidates; cand++ ) {
							if ( candidates->at(cand).ixyz != su_corr ) { //most likely a miss
								continue;
							}
							match = true;
							gridvalue = candidates->at(cand).gid;
							if ( gridvalue >= ndg )
								oid = myrxgpool[gridvalue-ndg].caori;
							else
								oid = mydefgpool[gridvalue].caori;

							ang << setprecision(8) << myoripool[oid].bunge1 << ";" << setprecision(8) << myoripool[oid].bunge2 << ";" << setprecision(8) << myoripool[oid].bunge3 << endl;
							break; //found my match!
						}
						//handle nothing was found case
						if ( match == false ) { 
							ang << "----> no identifiable orientation in CURRENTLY_INFECTED cell!" << endl; //setprecision(8) << 0.0 << ";" << setprecision(8) << 0.0 << ";" << setprecision(8) << 0.0 << endl;
						}
					} //done with identifying the CURRENTLY_INFECTED cell
					*/

				} //scan through xline in region r
			} //scan +y the xlines stacked in region r

			//fill in remaining values inside this region
			uint32_t ncandidates = this->cellStatii.at(r)->size();
			vector<cellstatus>* candidates = this->cellStatii.at(r);
			for ( uint32_t cand = 0; cand < ncandidates; cand++ ) {
				//is the infected pixel in layer imgz?
				uint32_t zzz = candidates->at(cand).ixyz / (imgxy);
				if ( zzz != imgz ) { //most likely we're not...
					continue;
				}
				//in this plane, zzz == imgz, i.e. we are in the correct layer and require an orientation for a cell in state == CURRENTLY_INFECTED
				uint32_t rem = candidates->at(cand).ixyz - (imgxy * zzz);
				uint32_t yyy = rem / imgx;
				uint32_t xxx = rem - (imgx * yyy); //x,y,z are SU global coordinates and thus translate into PNG image coordinates
				uint32_t px = xxx + (yyy*imgx);

				gridvalue = candidates->at(cand).gid;
				if ( gridvalue >= ndg ) {
					oid = myrxgpool[gridvalue-ndg].caori;
					iobuffer[4*px+3] = +1.0 * (gridvalue - ndg + 1); //##MK::for MartinDiehl Fortran indexed ID of RX grain
				}
				else {
					oid = mydefgpool[gridvalue].caori;
					iobuffer[4*px+3] = -1.0 * (gridvalue + 1); //##MK::for MartinDiehl Fortran indexed ID of deformed grain
				}

				iobuffer[4*px+0] = myoripool[oid].bunge1;
				iobuffer[4*px+1] = myoripool[oid].bunge2;
				iobuffer[4*px+2] = myoripool[oid].bunge3;

				//ang << setprecision(8) << myoripool[oid].bunge1 << ";" << setprecision(8) << myoripool[oid].bunge2 << ";" << setprecision(8) << myoripool[oid].bunge3 << endl;
				//break; //found my match!
			} //done checking all candidates of the region for contribution CURRENTLY_INFECTED pixel to the image
		} //next region stacked upon last one in +y
		//done with all regions along y who include the zsection at global SU coordinate imgz
	}

	//pipe values into file
	ang << setprecision(8);
	for ( uint32_t px = 0; px < imgx*imgy; px++ ) {
		ang << iobuffer[4*px+0] << ";" << iobuffer[4*px+1] << ";" << iobuffer[4*px+2] << ";" << iobuffer[4*px+3] << endl;
	}
	ang.flush();
	ang.close();

	delete [] iobuffer; iobuffer = NULL;

	//omp_unset_lock(&h5lock);
}


void caHdl::ompcrit_write_semebsd_core( double zpos, double xnow, uint32_t threshold, char mode ) {
	//CALLED FROM WITHIN PARALLEL REGION, but only master executes fishing in other threads memory
	//output ASCII texture file

	//omp_set_lock(&h5lock);

	uint32_t imgx = myCAGeometry.nboxedge_rd;
	uint32_t imgy = myCAGeometry.nboxedge_td;
	uint32_t imgxy = imgx*imgy;
	uint32_t imgz = (uint32_t) ( zpos * ((double) myCAGeometry.nboxedge_nd));
	if ( imgz == myCAGeometry.nboxedge_nd ) //MK::global coordinates [0,nboxedge_i) along each direction!
		imgz--;

	//see comments ompcrit_write_semebd_core
	//in which regions are the image data located? get region topology
	//uint32_t nregions = this->regions.size();
	uint32_t npx = this->regions[MASTER]->thePartitioning.nreg_rdx; //1 required, no domain split in x
	uint32_t npy = this->regions[MASTER]->thePartitioning.nreg_tdy;
	uint32_t npz = this->regions[MASTER]->thePartitioning.nreg_ndz;
	if ( npx != 1 ) { cout << "ERROR::Not outputting sections because npx != 1!" << endl; }

	//find the zlayer by scanning along the first npz region limits
	uint32_t zregion = 0;
	uint32_t zmi, zmx, zz;
	bool found = false;
	for ( zregion = 0; zregion < npz; zregion++ ) {
		zmi = this->regions[(zregion*npy)+0]->myGeom.nreg_ndmin; //global coordinates
		zmx = this->regions[(zregion*npy)+0]->myGeom.nreg_ndmax;

		if ( imgz >= zmi && imgz <= zmx ) {
			zz = imgz - zmi; //local grid coordinate z for all regions that include the z position imgz adnd participate in I/O
			found = true;
			break;
		}
	}


	//###need to implement scanning if sufficient size grain and in case add DIEHL_RXG_OFFSET
	double* iobuffer = NULL;
	iobuffer = new double[4*imgxy]; //+x then stack xlines along +y
	QUICKASSERT( iobuffer != NULL );
	for ( uint32_t px = 0; px < imgxy; px++ ) {
		iobuffer[4*px+0] = -4.0*_PI_; //dummy value to detect error
		iobuffer[4*px+1] = -4.0*_PI_;
		iobuffer[4*px+2] = -4.0*_PI_;
		iobuffer[4*px+3] = 0; //##MK::for Martin Diehl we use Fortran indexing no grain ID should be zero
	}

	if ( found == true ) {
		map<uint32_t,uint32_t> rxg_pixel_count;

		//the automaton comprises npz regions along z, x is aligned, and y as well, so simply read x*y slab stack in npy direction and drop such a global z section into the image
		uint32_t xper, yper;
		uint32_t r,region_corr, su_corr, gridvalue, oid;
		uint32_t ndg = this->mydefgpool.size();
		QUICKASSERT ( ndg == this->nmydefgpool );

		uint32_t xmi,xmx,ymi,ymx;
		for( uint32_t yregion = 0; yregion < npy; yregion++ ) { //collect data over regions
			r = (zregion*npy)+yregion;
			xper = this->regions[r]->myGeom.nreg_rd; //MK::extent in the r-th region not global SU coordinate
			yper = this->regions[r]->myGeom.nreg_td;
			xmi = this->regions[r]->myGeom.nreg_rdmin; //global coordinates delimitng memory region inclusively
			xmx = this->regions[r]->myGeom.nreg_rdmax;
			ymi = this->regions[r]->myGeom.nreg_tdmin;
			ymx = this->regions[r]->myGeom.nreg_tdmax;

			for ( uint32_t yy = 0; yy < yper; yy++ ) {
				uint32_t regyzoff = (xper * yy) + (xper * yper * zz);
				uint32_t yzoff = (ymi + yy)*imgx;
				for ( uint32_t xx = 0; xx < xper; xx++ ) {
					region_corr = xx + regyzoff;
					su_corr = xmi + xx + yzoff;

					gridvalue = regions[r]->mycellgrid[region_corr];

					if ( gridvalue <= CA_GRAINIDS_MAXIMUM ) { //therefore also != CURENTLY_INFECTED most likely either rx or defg
						if ( gridvalue >= ndg ) {
							uint32_t rxgid = gridvalue;
							auto it = rxg_pixel_count.find( rxgid );
							if ( it != rxg_pixel_count.end() )
								it->second++;
							else
								rxg_pixel_count.insert( pair<uint32_t,uint32_t>(rxgid, 1) );
						}
						//else in the accounting stage
					}
				} //scan through xline in region r
			} //scan +y the xlines stacked in region r

			//fill in remaining values inside this region
			uint32_t ncandidates = this->cellStatii.at(r)->size();
			vector<cellstatus>* candidates = this->cellStatii.at(r);
			for ( uint32_t cand = 0; cand < ncandidates; cand++ ) {
				//is the infected pixel in layer imgz?
				uint32_t zzz = candidates->at(cand).ixyz / (imgxy);
				if ( zzz != imgz ) { //most likely we're not...
					continue;
				}
				//in this plane, zzz == imgz, i.e. we are in the correct layer and require an orientation for a cell in state == CURRENTLY_INFECTED
				uint32_t rem = candidates->at(cand).ixyz - (imgxy * zzz);
				uint32_t yyy = rem / imgx;
				uint32_t xxx = rem - (imgx * yyy); //x,y,z are SU global coordinates and thus translate into PNG image coordinates
				uint32_t px = xxx + (yyy*imgx);

				gridvalue = candidates->at(cand).gid;
				if ( gridvalue >= ndg ) {
					uint32_t rxgid = gridvalue;
					auto it = rxg_pixel_count.find( rxgid );
					if ( it != rxg_pixel_count.end() )
						it->second++;
					else
						rxg_pixel_count.insert( pair<uint32_t,uint32_t>(rxgid, 1) );
				}
			} //done checking all candidates of the region for contribution CURRENTLY_INFECTED pixel to the image
		} //next region stacked upon last one in +y
		//done with all regions along y who include the zsection at global SU coordinate imgz

		//the automaton comprises npz regions along z, x is aligned, and y as well, so simply read x*y slab stack in npy direction and drop such a global z section into the image
		for( uint32_t yregion = 0; yregion < npy; yregion++ ) { //collect data over regions
			r = (zregion*npy)+yregion;
			xper = this->regions[r]->myGeom.nreg_rd; //MK::extent in the r-th region not global SU coordinate
			yper = this->regions[r]->myGeom.nreg_td;
			xmi = this->regions[r]->myGeom.nreg_rdmin; //global coordinates delimitng memory region inclusively
			xmx = this->regions[r]->myGeom.nreg_rdmax;
			ymi = this->regions[r]->myGeom.nreg_tdmin;
			ymx = this->regions[r]->myGeom.nreg_tdmax;

			for ( uint32_t yy = 0; yy < yper; yy++ ) {
				uint32_t regyzoff = (xper * yy) + (xper * yper * zz);
				uint32_t yzoff = (ymi + yy)*imgx;
				for ( uint32_t xx = 0; xx < xper; xx++ ) {
					region_corr = xx + regyzoff;
					su_corr = xmi + xx + yzoff;

					gridvalue = regions[r]->mycellgrid[region_corr];

					if ( gridvalue <= CA_GRAINIDS_MAXIMUM ) { //therefore also != CURENTLY_INFECTED most likely either rx or defg
						if ( gridvalue >= ndg ) {
							auto it = rxg_pixel_count.find( gridvalue );
							if ( it->second >= threshold ) {
								oid = myrxgpool[gridvalue-ndg].caori;
								iobuffer[4*su_corr+3] = +1.0 * (gridvalue - ndg + 1); //##MK::for MartinDiehl Fortran indexed ID of RX grain
							}
							else { //below threshold so ori of recrystallized but marker set that something is special
								oid = myrxgpool[gridvalue-ndg].caori;
								iobuffer[4*su_corr+3] = +1.0 * (gridvalue - ndg + DIEHL_RXG_OFFSET + 1); //##MK::for MartinDiehl Fortran indexed ID of RX grain
							}
						}
						else {
							oid = mydefgpool[gridvalue].caori;
							iobuffer[4*su_corr+3] = -1.0 * (gridvalue + 1); //##MK::for MartinDiehl Fortran indexed ID of deformed grain
						}
						iobuffer[4*su_corr+0] = myoripool[oid].bunge1;
						iobuffer[4*su_corr+1] = myoripool[oid].bunge2;
						iobuffer[4*su_corr+2] = myoripool[oid].bunge3;

						//ang << setprecision(8) << myoripool[oid].bunge1 << ";" << setprecision(8) << myoripool[oid].bunge2 << ";" << setprecision(8) << myoripool[oid].bunge3 << endl;
					}
				} //scan through xline in region r
			} //scan +y the xlines stacked in region r

			//fill in remaining values inside this region
			uint32_t ncandidates = this->cellStatii.at(r)->size();
			vector<cellstatus>* candidates = this->cellStatii.at(r);
			for ( uint32_t cand = 0; cand < ncandidates; cand++ ) {
				//is the infected pixel in layer imgz?
				uint32_t zzz = candidates->at(cand).ixyz / (imgxy);
				if ( zzz != imgz ) { //most likely we're not...
					continue;
				}
				//in this plane, zzz == imgz, i.e. we are in the correct layer and require an orientation for a cell in state == CURRENTLY_INFECTED
				uint32_t rem = candidates->at(cand).ixyz - (imgxy * zzz);
				uint32_t yyy = rem / imgx;
				uint32_t xxx = rem - (imgx * yyy); //x,y,z are SU global coordinates and thus translate into PNG image coordinates
				uint32_t px = xxx + (yyy*imgx);

				gridvalue = candidates->at(cand).gid;
				if ( gridvalue >= ndg ) {
					auto it = rxg_pixel_count.find( gridvalue );
					if ( it->second >= threshold ) {
						oid = myrxgpool[gridvalue-ndg].caori;
						iobuffer[4*px+3] = +1.0 * (gridvalue - ndg + 1); //##MK::for MartinDiehl Fortran indexed ID of RX grain
					}
					else { //something is special see above
						oid = myrxgpool[gridvalue-ndg].caori;
						iobuffer[4*su_corr+3] = +1.0 * (gridvalue - ndg + DIEHL_RXG_OFFSET + 1); //##MK::for MartinDiehl Fortran indexed ID of RX grain
					}
				}
				else {
					oid = mydefgpool[gridvalue].caori;
					iobuffer[4*px+3] = -1.0 * (gridvalue + 1); //##MK::for MartinDiehl Fortran indexed ID of deformed grain
				}

				iobuffer[4*px+0] = myoripool[oid].bunge1;
				iobuffer[4*px+1] = myoripool[oid].bunge2;
				iobuffer[4*px+2] = myoripool[oid].bunge3;

				//ang << setprecision(8) << myoripool[oid].bunge1 << ";" << setprecision(8) << myoripool[oid].bunge2 << ";" << setprecision(8) << myoripool[oid].bunge3 << endl;
				//break; //found my match!
			} //done checking all candidates of the region for contribution CURRENTLY_INFECTED pixel to the image
		} //next region stacked upon last one in +y
		//done with all regions along y who include the zsection at global SU coordinate imgz
	}

	stringstream ang_fname;
	ang_fname << this->myensHdl->simwrkdir << "/ang/" << "SCORE." << this->myensHdl->simid << ".JobID." << this->jobid << ".XYIter." << this->step << ".X.";
	uint32_t XX = (xnow * 1000.0);
	if ( XX < 10 )					ang_fname << "000" << XX;
	if ( XX >= 10 && XX < 100 )		ang_fname << "00" << XX;
	if ( XX >= 100 && XX < 1000 )	ang_fname << "0" << XX;
	if ( XX >= 1000 )				ang_fname << XX;
	ang_fname << ".ZPOS." << imgz << ".ang";

	ofstream ang;
	ang.open ( ang_fname.str().c_str() );

	//##MK::if really desiring ANG file add specific header, ##MK::currently implementing only dummy
	ang << "//SCORE auto-generated artificial SEM/EBSD file!\n";
	ang << "//SCORE coordinate system x || RD, y || TD, z || ND, right-handed\n";
	ang << "\n\n\n";
	ang << "BungePhi1;BungePHI;BungePhi2;FortranIndexedGrainID(PosForRXed,NegForDef)\n";

	//pipe values into file
	ang << setprecision(8);
	for ( uint32_t px = 0; px < imgx*imgy; px++ ) {
		ang << iobuffer[4*px+0] << ";" << iobuffer[4*px+1] << ";" << iobuffer[4*px+2] << ";" << iobuffer[4*px+3] << endl;
	}
	ang.flush();
	ang.close();

	delete [] iobuffer; iobuffer = NULL;

	//omp_unset_lock(&h5lock);
}


double caHdl::ompcrit_probe_semebsd_core( double zpos, uint32_t threshold ) {
	//CALLED FROM WITHIN PARALLEL REGION, but only master executes fishing in other threads memory
	//this function probes the accumulated recrystallized area for all grains larger than critical threshold
	//output ASCII texture file

	uint32_t imgx = myCAGeometry.nboxedge_rd;
	uint32_t imgy = myCAGeometry.nboxedge_td;
	uint32_t imgxy = imgx*imgy;
	uint32_t imgz = (uint32_t) ( zpos * ((double) myCAGeometry.nboxedge_nd));
	if ( imgz == myCAGeometry.nboxedge_nd ) //MK::global coordinates [0,nboxedge_i) along each direction!
		imgz--;

	//see comments for ompcrit_write_semebsd_core()

	//in which regions are the image data located? get region topology
	//uint32_t nregions = this->regions.size();
	uint32_t npx = this->regions[MASTER]->thePartitioning.nreg_rdx; //1 required, no domain split in x
	uint32_t npy = this->regions[MASTER]->thePartitioning.nreg_tdy;
	uint32_t npz = this->regions[MASTER]->thePartitioning.nreg_ndz;
	if ( npx != 1 ) { cout << "ERROR::Not outputting sections because npx != 1!" << endl; }

	//find the zlayer by scanning along the first npz region limits
	uint32_t zregion = 0;
	uint32_t zmi, zmx, zz;
	bool found = false;
	for ( zregion = 0; zregion < npz; zregion++ ) {
		zmi = this->regions[(zregion*npy)+0]->myGeom.nreg_ndmin; //global coordinates
		zmx = this->regions[(zregion*npy)+0]->myGeom.nreg_ndmax;

		if ( imgz >= zmi && imgz <= zmx ) {
			zz = imgz - zmi; //local grid coordinate z for all regions that include the z position imgz adnd participate in I/O
			found = true;
			break;
		}
	}

	/*map<uint32_t,uint32_t> rxg_pixel_count;*/
	vector<uint32_t> rxg_px_count(myrxgpool.size(), 0);

	if ( found == true ) {
		//the automaton comprises npz regions along z, x is aligned, and y as well, so simply read x*y slab stack in npy direction and drop such a global z section into the image
		uint32_t xper, yper;
		uint32_t r,region_corr, su_corr, gridvalue, oid;
		uint32_t ndg = this->mydefgpool.size();
		QUICKASSERT ( ndg == this->nmydefgpool );

		uint32_t xmi,xmx,ymi,ymx;
		for( uint32_t yregion = 0; yregion < npy; yregion++ ) { //collect data over regions
			r = (zregion*npy)+yregion;
			xper = this->regions[r]->myGeom.nreg_rd; //MK::extent in the r-th region not global SU coordinate
			yper = this->regions[r]->myGeom.nreg_td;
			xmi = this->regions[r]->myGeom.nreg_rdmin; //global coordinates delimitng memory region inclusively
			xmx = this->regions[r]->myGeom.nreg_rdmax;
			ymi = this->regions[r]->myGeom.nreg_tdmin;
			ymx = this->regions[r]->myGeom.nreg_tdmax;

			for ( uint32_t yy = 0; yy < yper; yy++ ) {
				uint32_t regyzoff = (xper * yy) + (xper * yper * zz);
				uint32_t yzoff = (ymi + yy)*imgx;
				for ( uint32_t xx = 0; xx < xper; xx++ ) {
					region_corr = xx + regyzoff;
					su_corr = xmi + xx + yzoff;

					gridvalue = regions[r]->mycellgrid[region_corr];

					if ( gridvalue <= CA_GRAINIDS_MAXIMUM ) { //therefore also != CURENTLY_INFECTED most likely either rx or defg
						if ( gridvalue >= ndg ) { //account for pixel count of recrystallized grain area section
							/*uint32_t rxgid = gridvalue - ndg + 1;
							auto it = rxg_pixel_count.find( rxgid );
							if ( it != rxg_pixel_count.end() ) { //likely value exists
								it->second++;
							}
							else {
								rxg_pixel_count.insert( pair<uint32_t,uint32_t>(rxgid, 1) );
							}*/
							uint32_t rxgid = gridvalue - ndg;
							rxg_px_count.at(rxgid)++;
						}
					}
				} //scan through xline in region r
			} //scan +y the xlines stacked in region r

			//fill in remaining values inside this region
			uint32_t ncandidates = this->cellStatii.at(r)->size();
			vector<cellstatus>* candidates = this->cellStatii.at(r);
			for ( uint32_t cand = 0; cand < ncandidates; cand++ ) {
				//is the infected pixel in layer imgz?
				uint32_t zzz = candidates->at(cand).ixyz / (imgxy);
				if ( zzz != imgz ) { //most likely we're not...
					continue;
				}
				//in this plane, zzz == imgz, i.e. we are in the correct layer and require an orientation for a cell in state == CURRENTLY_INFECTED
				uint32_t rem = candidates->at(cand).ixyz - (imgxy * zzz);
				uint32_t yyy = rem / imgx;
				uint32_t xxx = rem - (imgx * yyy); //x,y,z are SU global coordinates and thus translate into PNG image coordinates
				uint32_t px = xxx + (yyy*imgx);

				gridvalue = candidates->at(cand).gid;
				if ( gridvalue >= ndg ) {
					/*uint32_t rxgid = gridvalue - ndg + 1;
					auto it = rxg_pixel_count.find( rxgid );
					if ( it != rxg_pixel_count.end() ) { //likely value exists
						it->second++;
					}
					else {
						rxg_pixel_count.insert( pair<uint32_t,uint32_t>(rxgid, 1) );
					}*/
					uint32_t rxgid = gridvalue - ndg;
					rxg_px_count.at(rxgid)++;
				}
			} //done checking all candidates of the region for contribution CURRENTLY_INFECTED pixel to the image
		} //next region stacked upon last one in +y
		//done with all regions along y who include the zsection at global SU coordinate imgz
	}

	uint32_t sum_rx_area = 0;
	/*for( auto it = rxg_pixel_count.begin(); it != rxg_pixel_count.end(); ++it ) {
		if ( it->second >= threshold ) { //if grain section has minimum size we account for it
			//thereby we correct bias in distribution of sectioned rxg grain areas as
			//we account now like EBSD people in 2d experiments do, they exclude grain (section) if they do not have
			//threshold pixel count
			sum_rx_area += it->second;
		}
	}*/
	for( auto it = rxg_px_count.begin(); it != rxg_px_count.end(); ++it ) {
		if ( *it >= threshold ) {
			sum_rx_area += *it;
		}
	}

	double rx_area_zz = static_cast<double>(sum_rx_area) / static_cast<double>(imgxy);
	return rx_area_zz;
}


double caHdl::ompcrit_probe_semebsd_core2( double zpos, uint32_t threshold ) {
	//CALLED FROM WITHIN PARALLEL REGION, but only master executes fishing in other threads memory
	//this function probes the accumulated recrystallized area for all grains larger than critical threshold
	//output ASCII texture file

	uint32_t imgx = myCAGeometry.nboxedge_rd;
	uint32_t imgy = myCAGeometry.nboxedge_td;
	uint32_t imgxy = imgx*imgy;
	uint32_t imgz = (uint32_t) ( zpos * ((double) myCAGeometry.nboxedge_nd));
	if ( imgz == myCAGeometry.nboxedge_nd ) //MK::global coordinates [0,nboxedge_i) along each direction!
		imgz--;

	//see comments for ompcrit_write_semebsd_core()

	//in which regions are the image data located? get region topology
	uint32_t npx = this->regions[MASTER]->thePartitioning.nreg_rdx; //1 required, no domain split in x
	uint32_t npy = this->regions[MASTER]->thePartitioning.nreg_tdy;
	uint32_t npz = this->regions[MASTER]->thePartitioning.nreg_ndz;
	if ( npx != 1 ) { cout << "ERROR::Not outputting sections because npx != 1!" << endl; }

	//find the zlayer by scanning along the first npz region limits
	uint32_t zregion = 0;
	uint32_t zmi, zmx, zz;
	bool found = false;
	for ( zregion = 0; zregion < npz; zregion++ ) {
		zmi = this->regions[(zregion*npy)+0]->myGeom.nreg_ndmin; //global coordinates inclusive
		zmx = this->regions[(zregion*npy)+0]->myGeom.nreg_ndmax;

		if ( imgz >= zmi && imgz <= zmx ) {
			zz = imgz - zmi; //local grid coordinate z for all regions that include the z position imgz adnd participate in I/O
			found = true;
			break;
		}
	}

	vector<uint32_t> rxg_px_count = vector<uint32_t>( myrxgpool.size(), 0);

	if ( found == true ) {
		//the automaton comprises npz regions along z, x is aligned, and y as well, so simply read x*y slab stack in npy direction and drop such a global z section into the image
		uint32_t nrxg = this->myrxgpool.size();
		uint32_t ndg = this->mydefgpool.size();

		for( uint32_t yregion = 0; yregion < npy; yregion++ ) {
			//collect precomputed single rxg resolved serial sectioning data from thread regions
			//from which region?
			uint32_t r = (zregion*npy)+yregion;
			//local zmin of this region all along +y aligned thread sub-domains are at the same height level
			uint32_t zroff = zz*nrxg;
			vector<uint32_t> & thisone = this->regions[r]->rdtd_serialsection;
			for( uint32_t rxg = 0; rxg < nrxg; ++rxg ) {
				rxg_px_count[rxg] += thisone.at(zroff+rxg);
			}
		} //next region stacked upon last one in +y
	} //done with all regions along y who include the zsection at global SU coordinate imgz

	uint32_t sum_rx_area = 0;
	for( auto it = rxg_px_count.begin(); it != rxg_px_count.end(); ++it ) {
		if ( *it >= threshold ) {
			sum_rx_area += *it;
		}
	}

	double rx_area_zz = static_cast<double>(sum_rx_area) / static_cast<double>(imgxy);
	return rx_area_zz;
}



void caHdl::ompcrit_write_semebsd_looping( char thatmode ) {

	//ompcrit_write_semebsd_core( 0.0, thatmode );
	for ( uint32_t zr = 0; zr < rendering_atthisZPos.size(); zr++ ) {
		ompcrit_write_semebsd_core( rendering_atthisZPos.at(zr), thatmode );
	}
}


void caHdl::test_lodepng( void ) {
	
	//allocate memory to store the resulting image
	uint32_t imgx = 256;
	uint32_t imgy = 256;

	unsigned char* img_rgba = NULL;
	img_rgba = new unsigned char [4*imgx*imgy];

	uint32_t c = 0;
	for ( uint32_t y = 0; y < imgy; y++ ) {
		for ( uint32_t x = 0; x < imgx; x++ ) {
			c = 4*(x + (y*imgx));
			img_rgba[c + REDCHAN ] = BLACK;
			img_rgba[c + GREENCHAN ] = BLACK;
			img_rgba[c + BLUECHAN ] = BLACK;
			img_rgba[c + ALPHACHAN ] = y; //y; //works because imgy < 256 and positive UCHAR_RANGE_MAX;

			//insert naively a 120 120 120 circle inside
			double xx = x;
			double yy = y;
			double rr = SQR(50.0);
			if ( (SQR(xx - 0.5*256.0) + SQR(yy + 20.0 - 0.5*256.0)) < rr ) {
				img_rgba[c + REDCHAN ] = (unsigned char) 255;
				img_rgba[c + GREENCHAN ] = (unsigned char) 120;
				img_rgba[c + BLUECHAN ] = (unsigned char) 0;
				img_rgba[c + ALPHACHAN ] = UCHAR_RANGE_MAX;
			}

		}
	}

	string fname;
	fname = "SCORE.TestingLodePNG.png";

	//encode the PNG file by calling the lodePNG library
	lodepng::encode( fname.c_str(), img_rgba , imgx, imgy );

	delete [] img_rgba; img_rgba = NULL;
}


#define DAMASK_GRAININDEXING_OFFSET			1		//MK::because in SCORE grains start at 0
#define DAMASK_TEXTUREINDEXING_OFFSET		1

void caHdl::write_damask_matconfig_ascii( string suffix )
{
	//CALLED FROM WITHIN A PARALLEL REGION AND MUST NOT BE EXECUTED BY OTHERS THAN THE MASTER THREAD
	//MK::one can change this then however the OMP constructs in the solve_RXGROWTH need to be adjusted accordingly!
	if ( omp_get_thread_num() != MASTER ) 
		return;

	omp_set_lock(&h5lock);

	stringstream damask_matconfig_fname;
	ofstream damask_matconfig;

	damask_matconfig_fname << "SCORE." << myensHdl->simid << ".MyRank." << myensHdl->myRank << ".JobID." << this->jobid << "." << suffix << ".DAMASKMaterial.config";
	damask_matconfig.open ( damask_matconfig_fname.str().c_str() );

	//header
	damask_matconfig << "//SCORE auto-generated DAMASK material.config file!\n";
	damask_matconfig << "//Mind that not all parts of the material.config file required are generated by this tool!!!\n";
	damask_matconfig << "//We assume at the moment only one crystallite section and only one phase section and volume fraction 1.0\n";
	damask_matconfig << "//End of comments for SCORE autogenerated DAMASK material.config file\n";
	damask_matconfig << "\n\n\n"; 

	//write microstructure section
	damask_matconfig << "#-------------------#\n";
	damask_matconfig << "<microstructure>\n";
	damask_matconfig << "#-------------------#\n";
	//write first the deformed grains
	unsigned int grainid, textureid;
	for ( unsigned int dg = 0; dg < this->mydefgpool.size(); dg++ ) {
		grainid = DAMASK_GRAININDEXING_OFFSET + dg;
		textureid = DAMASK_TEXTUREINDEXING_OFFSET + this->mydefgpool[dg].caori;
		damask_matconfig << "[Grain" << grainid << "]\n";
		damask_matconfig << "crystallite 1\n";
		damask_matconfig << "(constituent)  phase 1   texture " << textureid << "   fraction 1.0\n";
	}
	for ( unsigned int rxg = 0; rxg < this->myrxgpool.size(); rxg++ ) {
		grainid = DAMASK_GRAININDEXING_OFFSET + nmydefgpool + rxg;
		textureid = DAMASK_TEXTUREINDEXING_OFFSET + this->myrxgpool[rxg].caori;
		damask_matconfig << "[Grain" << grainid << "]\n";
		damask_matconfig << "crystallite 1\n";
		damask_matconfig << "(constituent)  phase 1   texture " << textureid << "   fraction 1.0\n";
	}
	damask_matconfig << endl << endl;
	damask_matconfig << "#-------------------#\n";
	damask_matconfig << "<texture>\n";
	damask_matconfig << "#-------------------#\n";

	for ( unsigned int tx = 0; tx < this->myoripool.size(); tx++ ) {
		damask_matconfig << "[Texture" << (DAMASK_TEXTUREINDEXING_OFFSET + tx) << "]\n";
		damask_matconfig << "(gauss)  phi1 " << RAD2DEG(myoripool[tx].bunge1) << "    Phi " << RAD2DEG(myoripool[tx].bunge2) << "    phi2 " << RAD2DEG(myoripool[tx].bunge3) << "   scatter 0.0   fraction 1.0\n";
	}
	damask_matconfig << endl << endl;


	damask_matconfig.flush();
	damask_matconfig.close();

	omp_unset_lock(&h5lock);
}


void caHdl::write_damask_geom_ascii( string suffix )
{
	//CALLED FROM WITHIN A PARALLEL REGION AND MUST NOT BE EXECUTED BY OTHERS THAN THE MASTER THREAD
	//MK::one can change this then however the OMP constructs in the solve_RXGROWTH need to be adjusted accordingly!
	if ( omp_get_thread_num() != MASTER ) 
		return;

	//region topology
	uint32_t npx = this->regions[MASTER]->thePartitioning.nreg_rdx; //1 required
	uint32_t npy = this->regions[MASTER]->thePartitioning.nreg_tdy;
	uint32_t npz = this->regions[MASTER]->thePartitioning.nreg_ndz;
	if ( npx != 1 ) { 
		std::cout << "ERROR::Invalid region topology based on current hybrid MPI/OMP data/memory model!" << std::endl;
		return;
	}

	omp_set_lock(&h5lock);

	stringstream damask_geom_fname;
	ofstream damask_geom;

	damask_geom_fname << "SCORE." << myensHdl->simid << ".MyRank." << myensHdl->myRank << ".JobID." << this->jobid << "." << suffix << ".DAMASK.geom";
	damask_geom.open ( damask_geom_fname.str().c_str() );

	//header
	damask_geom << "//SCOREv1.2 auto-generated DAMASK geometry file!\n";
	damask_geom << "//Now coordinate systems SCORE and DAMASK both right-handed +x || RD pointing right (right), +y || TD pointing inwards (rear), and +z ||ND upwards (top), while\n";
	damask_geom << "//End of comments for SCORE autogenerated DAMASK geometry file\n";
	damask_geom << "\n\n\n";

	//write DAMASK header
	damask_geom << "5\theader\n";
	damask_geom << "grid\ta " << this->myCAGeometry.nboxedge_rd << "\tb " << this->myCAGeometry.nboxedge_td << "\tc " << this->myCAGeometry.nboxedge_nd << "\n";
	damask_geom << "size\tx 1.000000\ty 1.000000\t z 1.000000\n";
	damask_geom << "origin\tx 0.000000\ty 0.000000\t z 0.000000\n";
	damask_geom << "microstructures " << (this->mydefgpool.size() + this->myrxgpool.size() ) << "\n";
	damask_geom << "homogenization\t1" << endl; 
	//MK::these make the output not portable but one can change here the structure accordingly...

	//write uncompressed block of microstructure data...
	//the automaton comprises npz regions along z, x is aligned, and y as well, so simply read x*y slab stack in npy direction and drop such a global z section into rawdata
	//##MK::we work sequentially through xy layers along z is performed
	uint32_t r, corr, gridvalue, yper, xper, tmp;

	for (uint32_t zregion = 0; zregion < npz; zregion++ ) { //utilize that the regions form a regular-stacked 3D region aggregate
		for ( uint32_t zz = 0; zz < this->regions[(zregion*npy)+0]->myGeom.nreg_nd; zz++ ) { //layer-by-layer MPI I/O several MB per write
			//collect data over regions in the zlayer
			for(uint32_t yregion = 0; yregion < npy; yregion++ ) {
				r = (zregion*npy)+yregion;
				yper = this->regions[r]->myGeom.nreg_td;
				xper = this->regions[r]->myGeom.nreg_rd;

				for ( uint32_t yy = 0; yy < yper; yy++ ) {
					for ( uint32_t xx = 0; xx < xper; xx++ ) { 
						//##MK::utilizing for the formatting of the geom file that currently npx == 1 such that xper == nboxedge_rd

						//translate into local coordinate in memory region
						corr = xx + (xper * yy) + (xper * yper * zz);
						gridvalue = regions[r]->mycellgrid[corr];

						if ( gridvalue != CELL_IS_A_PARTICLE ) {
							if ( gridvalue == CURRENTLY_INFECTED ) {

								tmp = characterize_cellstatus( r, xx, yy, zz );

								if ( tmp != NO_GRAIN_ASSIGNED ) {
									if ( xx < (xper - 1) ) 
										damask_geom << (DAMASK_GRAININDEXING_OFFSET + tmp) << " ";
									else
										damask_geom << (DAMASK_GRAININDEXING_OFFSET + tmp) << endl; //xx == xper - 1 finish that line
									continue;
								}
								//else NO_GRAIN_ASSIGNED )
								std::cout << "Found CELL_IS_A_PARTICLE, damask output not implemented to handle this currently, now stopping!" << std::endl;
								damask_geom.flush(); damask_geom.close(); return;
							}

							//already correctly determined either deformed or recrystallized grain i.e. indices either [0,mydefgpool.size) (deformed grain) or [mydefpool.size+0,mydefpool.size+0+myrxgpool.size) (rx grain)
							if ( xx < (xper - 1) ) 
								damask_geom << (DAMASK_GRAININDEXING_OFFSET + gridvalue) << " ";
							else
								damask_geom << (DAMASK_GRAININDEXING_OFFSET + gridvalue) << endl;
						}
						else {
							std::cout << "Found CELL_IS_A_PARTICLE, damask output is currently not implemented to handle this, now stopping!" << std::endl;
							damask_geom.flush(); damask_geom.close(); return;
						}
					} //scan through xline
				} //scan +y the xlines stacked in region r
			} //next region stacked upon last one in +y
		} //next zslice in zregions
	} //next zregion with regions as xy slab

	damask_geom.flush();
	damask_geom.close();

	omp_unset_lock(&h5lock);
}


void caHdl::write_junction_skeleton( unsigned char mode )
{
	double t0 = MPI_Wtime();
	MPI_File msFileHdl;
	MPI_Status msFileStatus;

	// create C-consistent file name for MPI I/O
	stringstream msFileName;
	msFileName << "SCORE." << myensHdl->simid << ".JobID." <<  this->jobid << ".Iter." << this->step;
	if ( mode == OUTPUT_JUNCTIONS_GBTHREADCOLORING )
		msFileName << ".GBThreadID.raw";
	if ( mode == OUTPUT_JUNCTIONS_GBFACES )
		msFileName << ".GBFaces.raw";
	if ( mode ==  OUTPUT_JUNCTIONS_TRIJUNCTIONS )
		msFileName << ".TripleJunc.raw";
	if ( mode == OUTPUT_JUNCTIONS_HIGHERORDER )
		msFileName << ".HigherJunc.raw";

	uint32_t msFileNameLength = msFileName.str().size();
	char* CmsFileName = new char[msFileNameLength+1];
	strcpy(CmsFileName, msFileName.str().c_str());

	//split analysis into as few pieces as possible because each piece analysis requires the scanning of all boundary cells
	//however only one preallocated array may exceed the maximum size of a single MPI_File_write op therefore need a compromise solution
	uint32_t* rawdata = NULL;
	long nxyz = myCAGeometry.nboxvol_rdtdnd; //promotion to long no problem on 64bit architecture
	long nportionDesired = ((double) DEFAULT_GBPLOTTING_CACHESIZE / (double) sizeof(uint32_t));
	if ( nportionDesired < 1 ) { std::cout << "WARNING::No plotting because data portion too small!" << std::endl; return; }
	long nportionDone = 0;
	long nportionNow = 0;
	//MK::datatype long necessary because nxyz can exceed the INT32_RANGE !

	//only the MASTER node writes in parallel! in write-only mode
	MPI_File_open(MPI_COMM_SELF, CmsFileName, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &msFileHdl);
	__int64 totalOffset = 0; //MK:: with 64-bit architectures, offsets calculations can be performed as longs 2^64 bit thus on petabyte files...	//in mpi.h MPI_Offset is defined as an __int64 which is long long, thus we can jump much more than 2^32 directly when unsigned int would be utilized
	MPI_File_seek( msFileHdl, totalOffset, MPI_SEEK_SET );


	while( nportionDone < nxyz ) { //first fill cache then write cached data to the file
		nportionNow = nportionDesired;
		if ( nxyz - nportionDone < nportionNow ) //limit buffer size when coming to end of the implicit array namely then process only remaining cells
			nportionNow = nxyz - nportionDone;

		rawdata = NULL;
		rawdata = new uint32_t[nportionNow];
		QUICKASSERT( rawdata != NULL );
		for ( long i = 0; i < nportionNow; i++ ) //set dummy value
			rawdata[i] = 0;

		uint32_t tid = MASTER;
		uint32_t nthreads = this->regions.size(); //MK::something like nthreads = omp_get_num_threads() gives wrong results because we are currently not in a parallel region

		if ( mode == OUTPUT_JUNCTIONS_GBTHREADCOLORING ) {
			for ( uint32_t pmax = 0; pmax < this->regions[MASTER]->myjunctions.fjuncSkeletonSize; pmax++ ) {
				tid = workPartitioning(pmax, nthreads); //AGAIN THER SAME STATIC WORK PARTITIONING
				for ( uint32_t f = 0; f < this->regions[tid]->myjunctions.fjuncSkeleton[pmax].len; f++ ) {
					uint32_t flen = this->regions[tid]->myjunctions.fjuncSkeleton[pmax].thefaces[f].nvoxeltotal;
					cellsBndFast* bucket = this->regions[tid]->myjunctions.fjuncSkeleton[pmax].thefaces[f].voxelizedface;
					for ( uint32_t c = 0; c < flen; c++ ) {
						//MK::nportionDone and nportionNow are implicit global coordinates of the entire solitary unit domain volume!
						if ( bucket[c].location >= nportionDone && bucket[c].location < (nportionDone + nportionNow) ) 
							rawdata[bucket[c].location - nportionDone] = 1 + tid;
					}
				}
			}
		}
		if ( mode == OUTPUT_JUNCTIONS_GBFACES ) {
			for ( uint32_t pmax = 0; pmax < this->regions[MASTER]->myjunctions.fjuncSkeletonSize; pmax++ ) { //fjuncSize the same on all threads after consolidation
				tid = workPartitioning(pmax, nthreads);
				for ( uint32_t f = 0; f < this->regions[tid]->myjunctions.fjuncSkeleton[pmax].len; f++ ) {
					uint32_t flen = this->regions[tid]->myjunctions.fjuncSkeleton[pmax].thefaces[f].nvoxeltotal;
					cellsBndFast* bucket = this->regions[tid]->myjunctions.fjuncSkeleton[pmax].thefaces[f].voxelizedface;
					for ( uint32_t c = 0; c < flen; c++ ) {
						if ( bucket[c].location >= nportionDone && bucket[c].location < (nportionDone + nportionNow) ) 
							rawdata[bucket[c].location - nportionDone] = 1 + this->regions[tid]->myjunctions.fjuncSkeleton[pmax].thefaces[f].gbID; //##MK::colud be optimize to reduce lenght of necessary dereferencing
					}
				}
			}
		}
		if ( mode == OUTPUT_JUNCTIONS_TRIJUNCTIONS ) {
			for ( uint32_t thr = MASTER; thr < nthreads; thr++ ) {
				//MK::tjunc and hjunc sizes are in most cases different across the threads!
				uint32_t ntj = this->regions[thr]->myjunctions.tjuncSkeletonSize;
				cellsTJ* bucket = this->regions[thr]->myjunctions.tjuncSkeleton;
				for ( uint32_t tj = 0; tj < ntj; tj++ ) { //TRI junctions not consolidated yet
					if ( bucket[tj].location >= nportionDone && bucket[tj].location < (nportionDone + nportionNow) ) 
							rawdata[bucket[tj].location - nportionDone] = 1 + thr;
				}
			}
		}
		if ( mode == OUTPUT_JUNCTIONS_HIGHERORDER ) {
			for ( uint32_t thr = MASTER; thr < nthreads; thr++ ) {
				uint32_t nhj = this->regions[thr]->myjunctions.hjuncSkeletonSize;
				cellsHJ* bucket = this->regions[thr]->myjunctions.hjuncSkeleton;
				for ( uint32_t hj = 0; hj < nhj; hj++ ) { //HIGHER ORDER junctions not consolidated yet
					if ( bucket[hj].location >= nportionDone && bucket[hj].location < (nportionDone + nportionNow) ) 
							rawdata[bucket[hj].location - nportionDone] = 1 + thr;
				}
			}
		}

		//push into file at once, runtime environment will optimize, further potential by nowait multiple write from the thread, but tricky...
		//even though nxyz may exceed the range of INT32 and therefore the nportionNow must be long after having set the value of nportionNow we are
		//guaranteed nportionNow < INT32_RANGE such that long -> int is safe!
		MPI_File_write( msFileHdl, rawdata, (int) nportionNow, MPI_UNSIGNED, &msFileStatus); //implicit advancement of fp
		nportionDone += nportionNow;

		delete [] rawdata;
		rawdata = NULL;
	}

	delete [] CmsFileName; CmsFileName = NULL;
	std::cout << myensHdl->myRank << " I/O " << mode << " in timestep " << this->step << " in " << ( (double) MPI_Wtime() - t0 ) << " seconds" << std::endl;

	MPI_File_close(&msFileHdl); //no Barrier as MPI_COMM_SELF
}


void caHdl::write_junction_log_gb( uint32_t thread, bool consolidated )
{
	//MK::if consolidated == true threads write only boundaries of grain pairs which they handle in local memory
	stringstream gblogfname;
	ofstream gblog;
	gblogfname << "SCORe." << myensHdl->simid << ".JobID." << this->jobid << ".Thread." << thread;
	if ( consolidated == true )
		gblogfname << ".ConsolidatedBoundaries.csv";
	else
		gblogfname << ".UnconsolidatedBoundaries.csv";
	gblog.open( gblogfname.str().c_str() );
	gblog << "ConsecutiveGBID;InternalPosition;CallingThread;ManagingThread;gPosUpID;gPosDownID;CellCnt;Disori;rho0gPosUp;rho0gPosDown" << endl;

	uint32_t thr = thread;
	uint32_t nthreads = this->regions.size();
	uint32_t npmax = this->regions[MASTER]->myjunctions.fjuncSkeletonSize;
	uint32_t nf = 0;

	for( uint32_t pmax = 0; pmax < npmax; pmax++ ) {
		if ( consolidated == true ) 
			thr = workPartitioning(pmax, nthreads); //accessing remote memory in shared memory
		else
			thr = thread;

		nf = this->regions[thr]->myjunctions.fjuncSkeleton[pmax].len;
		if ( nf != 0 ) {
			bndFaceFast* theface = NULL;
			for ( uint32_t f = 0; f < nf; f++) {
				theface = &(this->regions[thr]->myjunctions.fjuncSkeleton[pmax].thefaces[f]);
				gblog << theface->gbID << ";" << nf << ";" << thread << ";" << thr << ";" << theface->gposUp << ";" << theface->gposDown << ";" << theface->nextfreevoxelslot << ";" << theface->disori << ";" << mydefgpool[tmpdefgseeds[theface->gposUp].mydefgpoolid].rho0 << ";" << mydefgpool[tmpdefgseeds[theface->gposDown].mydefgpoolid].rho0 << "\n";
			}
		}
	}
	gblog << endl;
	gblog.flush();
	gblog.close();

	std::cout << "Grain boundary logfile written" << std::endl;
}


void caHdl::write_junction_log_tj( void )
{
	std::cout << "write_junction_log_tj remains to be implemented!" << std::endl;
}


void caHdl::write_junction_log_hj( void )
{
	std::cout << "write_junction_log_hj remains to be implemented!" << std::endl;
}


#endif
