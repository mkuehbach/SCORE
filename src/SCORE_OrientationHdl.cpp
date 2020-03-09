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


uint32_t ensembleHdl::ens_get_closest_standardlage( double * quat )
{
	double closest_disori = MAX_FCC;
	uint32_t closest_ideal = CATEGORIZED_AS_RANDOM;

	double qs[4] = {1.0, 0.0, 0.0, 0.0 };
	for ( uint32_t cand = 0; cand < standardlagen.size(); cand++) {
		qs[0] = standardlagen[cand].q0;
		qs[1] = standardlagen[cand].q1;
		qs[2] = standardlagen[cand].q2;
		qs[3] = standardlagen[cand].q3;

		double disori = ensmath.disori_angle_oriquat_cubic( quat, qs );

		//scatter within the range and closer as to all other candidates?
		if ( disori <= standardlagen[cand].scatter && disori <= closest_disori ) {
			closest_disori = disori;
			closest_ideal = cand;
		}
		//for all possible standardlagen candidates find always the closest match even when scatter range of two components overlaps such that macrotexture properly adds to 1.0
	}

	//MK::within the code orientations categorized as RANDOM are at the zero-th position
	if ( closest_ideal == CATEGORIZED_AS_RANDOM )
		return 0;
	//else
	return closest_ideal + 1;
}

uint32_t ensembleHdl::ens_check_disjunctness_core( double * q_ori, bool thesebunge, double phi1, double Phi, double phi2 )
{
	//check disorientation to all other components in worldoripool, adds to worldoripool if dissimilar enough
	//if thesebunge == true, the provided phi1, Phi, phi2 is utilized otherwise the function internally converts from quat2euler
	double qcand[4] = {1.0, 0.0, 0.0, 0.0 };
	double disori = MAX_FCC;
	double closestdisori = RESOLUTION_SO3GRID;
	uint32_t closestid = UNKNOWN_ORIENTATION;

	for (uint32_t cand = 0; cand < this->worldoripool.size(); cand++ ) {
		qcand[0] = worldoripool[cand].q0;
		qcand[1] = worldoripool[cand].q1;
		qcand[2] = worldoripool[cand].q2;
		qcand[3] = worldoripool[cand].q3;

		disori = ensmath.disori_angle_oriquat_cubic( q_ori, qcand );

		//closer than the SO3GRID resolution to any other orientation already known?
		if (disori <= closestdisori) {
			closestdisori = disori;
			closestid = cand;
		}
		//however still testing against all other to get the closest match
	}

	//most likely case
	if ( closestid == UNKNOWN_ORIENTATION ) {
		struct ori anori;

		//only once conversion of a new orientation into Bunge representation
		if ( thesebunge == true ) {
			anori.bunge1 = phi1;
			anori.bunge2 = Phi;
			anori.bunge3 = phi2;
		}
		else { //convert quaternion into Bunge internally
			double bunge_ori[3] = {0.0, 0.0, 0.0};
			ensmath.quat2euler( q_ori, bunge_ori );
			anori.bunge1 = bunge_ori[0];
			anori.bunge2 = bunge_ori[1];
			anori.bunge3 = bunge_ori[2];
		}

		anori.q0 = q_ori[0];
		anori.q1 = q_ori[1];
		anori.q2 = q_ori[2];
		anori.q3 = q_ori[3];

		anori.closestideal = ens_get_closest_standardlage( q_ori ); //RANDOM or one of our components

		anori.RGB_R = UCHAR_RANGE_MIN;
		anori.RGB_G = UCHAR_RANGE_MIN;
		anori.RGB_B = UCHAR_RANGE_MIN;
		anori.RGB_A = UCHAR_RANGE_MAX;

		worldoripool.push_back( anori );

#ifdef REPORTSTYLE_DEVELOPER
		if (myRank == MASTER) {
			cout << "ADD\t\t" << (worldoripool.size() - 1) << "\t\t" << anori.bunge1 << "\t\t" << anori.bunge2 << "\t\t" << anori.bunge3;
			cout << "\t\t" << anori.q0 << "\t\t" << anori.q1 << "\t\t" << anori.q2 << "\t\t" << anori.q3 << endl;
			if ( anori.closestideal == 0 ) { cout << "Node;" << this->myRank << ", I have categorized close to RANDOM." << endl; }
			else { cout << "Node;" << this->myRank << ", I have categorized close to listentry " << (anori.closestideal + 1) << endl; }
		}
		//#endif DETAILED_PROMPTS
#endif

		//make known the new orientation
		return (worldoripool.size() - 1);
	}

	//else present closest already existent orientation id
	//cout << "GET\t\t" << closestid << "\t\t" << closestdisori << endl;

	//one of the goodfellas...
	return closestid;
}

uint32_t ensembleHdl::ens_check_disjunctness_io( double * bunge_ori )
{
	//utilize only for reading in datasets, otherwise stay within quaternion representation so call ca_check_disjunctness_core with a quaternion
	QUICKASSERT ( worldoripool.size() < MAXIMUM_DISJOINT_ORIS );	//still accepting orientations?

	double quaternion_ori[4] = {1.0, 0.0, 0.0, 0.0};
	ensmath.euler2quat( bunge_ori, quaternion_ori );				//only once convert a Bunge Euler orientation into a quaternion

	uint32_t newid = ens_check_disjunctness_core( quaternion_ori, true, bunge_ori[0], bunge_ori[1], bunge_ori[2] );

	return newid;
}



uint32_t caHdl::ca_get_closest_standardlage( double  * quat )
{
	double closest_disori = MAX_FCC;
	uint32_t closest_ideal = CATEGORIZED_AS_RANDOM;
	ensembleHdlP ens = this->myensHdl;

	double qs[4] = {1.0, 0.0, 0.0, 0.0 };
	for ( uint32_t cand = 0; cand < ens->standardlagen.size(); cand++) {
		qs[0] = ens->standardlagen[cand].q0;
		qs[1] = ens->standardlagen[cand].q1;
		qs[2] = ens->standardlagen[cand].q2;
		qs[3] = ens->standardlagen[cand].q3;

		double disori = localmath.disori_angle_oriquat_cubic( quat, qs );

		//scatter within the range and closer as to all other candidates?
		if ( disori <= ens->standardlagen[cand].scatter && disori <= closest_disori ) {
			closest_disori = disori;
			closest_ideal = cand;
		}
		//for all possible standardlagen candidates find always the closest match even when scatter range of two components overlaps such that macrotexture properly adds to 1.0
	}

	//MK::within the code orientations categorized as RANDOM are at the zero-th position
	if ( closest_ideal == CATEGORIZED_AS_RANDOM )
		return 0;
	//else
	return closest_ideal + 1;
}


uint32_t caHdl::ca_check_disjunctness_core( double* q_ori, bool thesebunge, double phi1, double Phi, double phi2 )
{
/*
	//check disorientation to all other components in myoripool
	uint32_t closestid = UNKNOWN_ORIENTATION;
	double qcand[4] = {1.0, 0.0, 0.0, 0.0 };
	double disori = MAX_FCC;
	double closestdisori = RESOLUTION_SO3GRID;

	for (uint32_t cand = 0; cand < this->myoripool.size(); cand++ ) {
		qcand[0] = myoripool[cand].q0;
		qcand[1] = myoripool[cand].q1;
		qcand[2] = myoripool[cand].q2;
		qcand[3] = myoripool[cand].q3;

		disori = localmath.disori_angle_oriquat_cubic( q_ori, qcand );

		//closer than the SO3GRID resolution to any other orientation already known?
		if (disori <= closestdisori) {
			closestdisori = disori;
			closestid = cand;
		}
		//however still testing against all other to get the closest match
	}

	//most likely case
	if ( closestid == UNKNOWN_ORIENTATION ) {
*/
		struct ori anori;

		if ( thesebunge == true ) {
			anori.bunge1 = phi1;
			anori.bunge2 = Phi;
			anori.bunge3 = phi2;
		}
		else {
			//only once conversion of a new orientation into Bunge representation
			double bunge_ori[3] = {0.0, 0.0, 0.0};
			localmath.quat2euler( q_ori, bunge_ori );

			anori.bunge1 = bunge_ori[0];
			anori.bunge2 = bunge_ori[1];
			anori.bunge3 = bunge_ori[2];
		}

		anori.q0 = q_ori[0];
		anori.q1 = q_ori[1];
		anori.q2 = q_ori[2];
		anori.q3 = q_ori[3];

		anori.closestideal = ca_get_closest_standardlage( q_ori ); //RANDOM or one of our components

		anori.RGB_R = UCHAR_RANGE_MIN;
		anori.RGB_G = UCHAR_RANGE_MIN;
		anori.RGB_B = UCHAR_RANGE_MIN;
		anori.RGB_A = UCHAR_RANGE_MAX;


		myoripool.push_back( anori );

		unsigned int ii = myoripool.size()-1;
cout << "MyOriPool adding " << myoripool[ii].bunge1 << ";" << myoripool[ii].bunge2 << ";" << myoripool[ii].bunge3 << ";" << myoripool[ii].q0 << ";" << myoripool[ii].q1 << ";" << myoripool[ii].q2 << ";" << myoripool[ii].q3 << endl;

//cout << "ADD\t\t" << (myoripool.size() - 1) << "\t\t" << anori.bunge1 << "\t\t" << anori.bunge2 << "\t\t" << anori.bunge3;
//cout << "\t\t" << anori.q0 << "\t\t" << anori.q1 << "\t\t" << anori.q2 << "\t\t" << anori.q3 << endl;
//#ifdef DETAILED_PROMPTS
//if ( anori.closestideal == RANDOM_ORIENTATION ) { cout << "Node;" << this->myRank << ", I have categorized close to RANDOM." << endl; }
//else { cout << "Node;" << this->myRank << ", I have categorized close to listentry " << (anori.closestideal + 1) << endl; }
//#endif DETAILED_PROMPTS

		//make known the new orientation
		return (myoripool.size() - 1);
/*
	}

	//else present closest already existent orientation id
	//cout << "GET\t\t" << closestid << "\t\t" << closestdisori << endl;

	//one of the goodfellas...
cout << "MyOriPool the Euler triplet " << phi1 << ";" << Phi << ";" << phi2 << " is already known closestid/disori " << closestid << "\t\t" << closestdisori << endl;
cout << "MyOriPool instead the following orientation is taken " << myoripool[closestid].bunge1 << ";" << myoripool[closestid].bunge2 << ";" << myoripool[closestid].bunge3 << endl;
	return closestid;
*/
}


uint32_t caHdl::ca_check_disjunctness_io( double* bunge_ori )
{
	//utilize only for reading in datasets, otherwise stay within quaternion representation so call ca_check_disjunctness_core with a quaternion
	QUICKASSERT ( myoripool.size() < MAXIMUM_DISJOINT_ORIS );	//still accepting orientations?

	double quaternion_ori[4] = {1.0, 0.0, 0.0, 0.0};
	localmath.euler2quat( bunge_ori, quaternion_ori );				//only once convert a Bunge Euler orientation into a quaternion

	uint32_t newid = ca_check_disjunctness_core( quaternion_ori, true, bunge_ori[0], bunge_ori[1], bunge_ori[2] );

	return newid;
}


void caHdl::colorize_myoripool_ipfz( void )
{
//##MK::OpenMP parallelization potential as oriquat2sst_ipfz_rgb and subroutines are threadsafe
	for ( uint32_t oi = 0; oi < myoripool.size(); oi++) {
		double q[4] = {myoripool[oi].q0, myoripool[oi].q1, myoripool[oi].q2, myoripool[oi].q3 };
		double ax[3] = {0.0, 0.0, 1.0}; //z ||ND for instance

		unsigned char rgb[3] = { UCHAR_RANGE_MIN, UCHAR_RANGE_MIN, UCHAR_RANGE_MIN };
		double locationinsst[2] = {FAIL_MYMATH_NUMERIC, FAIL_MYMATH_NUMERIC};

		localmath.oriquat2sst_cubic_ipf_quat2pos( q, ax, locationinsst );
		localmath.oriquat2sst_cubic_ipf_pos2rgb( locationinsst, rgb, IPF_COLOR_STRETCH_R, IPF_COLOR_STRETCH_G, IPF_COLOR_STRETCH_B ); //threadsafe as well as all dependencies

		myoripool[oi].RGB_R = rgb[0];
		myoripool[oi].RGB_G = rgb[1];
		myoripool[oi].RGB_B = rgb[2];
		myoripool[oi].RGB_A = UCHAR_RANGE_MAX;
//cout << "Bunge2IPFZ=" << oi << "\t" << (int) rgb[RED] << ";" << (int) rgb[GREEN] << (int) rgb[BLUE] << endl;

cout << "oi/phi1/Phi/phi2/RGB = " << oi << "\t" << RAD2DEG(myoripool[oi].bunge1) << ";" << RAD2DEG(myoripool[oi].bunge2) << ";" << RAD2DEG(myoripool[oi].bunge3) << "\t\t" << (int) myoripool[oi].RGB_R << ";" << (int) myoripool[oi].RGB_G << ";" << (int) myoripool[oi].RGB_B << endl;
	}
}


#endif
