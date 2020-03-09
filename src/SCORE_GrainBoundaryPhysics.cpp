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

double caHdl::calc_mobilityweight ( uint32_t rgpoolid , uint32_t dgpoolid )
{
	if ( this->mobilitymodel == MOBILITYMODEL_ROLLETTHOLM ) {
		return ( this->calc_mobilityweight_rollettholm( rgpoolid, dgpoolid ) );
	}
	//else
	return ( this->calc_mobilityweight_sebaldgottstein( rgpoolid, dgpoolid ) );
}


double caHdl::calc_mobilityweight_sebaldgottstein ( uint32_t rgpid , uint32_t dgpid )
{
	//Sebald Gottstein model
	double q1[4], q2[4], theta;

	uint32_t rgori = myrxgpool[rgpid].caori;
	uint32_t dgori = mydefgpool[dgpid].caori;

	//now calculate on the fly...
	q1[0] = myoripool[rgori].q0;
	q1[1] = myoripool[rgori].q1;
	q1[2] = myoripool[rgori].q2;
	q1[3] = myoripool[rgori].q3;

	q2[0] = myoripool[dgori].q0;
	q2[1] = myoripool[dgori].q1;
	q2[2] = myoripool[dgori].q2;
	q2[3] = myoripool[dgori].q3;

	//P=weightedMob for all HAGB already defined
	double weightedMob = 0.0;

	theta = localmath.disori_angle_oriquat_cubic( q1, q2 );
	//assign categorical intrinsic boundary mobility to disorientation among two orientation
	//LAGB detected, overwrite default value
	if( theta <= MAXDISORI_LAGB2HAGB ) {
		weightedMob = -1.0;

		return weightedMob;
	}

	//obviously a HAGB but also eligible for GS? check proximity to 40deg111

	double qdis[4] = {1.0, 0.0, 0.0, 0.0};
	localmath.disori_fromoriquat_cubic_imm( q1, q2, qdis );
	//the quaternion that describes a 40deg111 misorientation
	double maxDev40_111 = MAXDISORI_TO_40DEG111;
	double _sqrt3 = 1 / sqrt( 3.0 );
	double oneNinth = 1.0 / 9.0;
	double m40_111[4] = { cos( oneNinth * _PI_ ), _sqrt3 * sin( oneNinth * _PI_ ), _sqrt3 * sin( oneNinth * _PI_ ), _sqrt3 * sin( oneNinth * _PI_ ) };

	double dev_40_111 = localmath.disori_angle_oriquat_cubic( qdis, m40_111 );

	if( dev_40_111 < maxDev40_111 ) {
		weightedMob = SQR( cos( 0.5 * _PI_ * dev_40_111 / maxDev40_111 ) );
	}

	return weightedMob;
}


double caHdl::calc_mobilityweight_rollettholm ( uint32_t rgpid, uint32_t dgpid )
{
	double q1[4], q2[4], theta;

	uint32_t rgori = myrxgpool[rgpid].caori;
	uint32_t dgori = mydefgpool[dgpid].caori;

	//now calculate on the fly...
	q1[0] = myoripool[rgori].q0;
	q1[1] = myoripool[rgori].q1;
	q1[2] = myoripool[rgori].q2;
	q1[3] = myoripool[rgori].q3;

	q2[0] = myoripool[dgori].q0;
	q2[1] = myoripool[dgori].q1;
	q2[2] = myoripool[dgori].q2;
	q2[3] = myoripool[dgori].q3;

	theta = localmath.disori_angle_oriquat_cubic( q1, q2 );

	//P value allows to calculate the mobility as often occurring in works by Rollett and Holm via P * mRHHAGB
	double P = 1.0 - (myPhysData.RH_LAGBHAGBcut * exp( -1.0 * myPhysData.RH_LAGBHAGBtrans * pow( (theta/MAXDISORI_LAGB2HAGB), myPhysData.RH_LAGBHAGBexponent ) ));

	return P;
}

#endif
