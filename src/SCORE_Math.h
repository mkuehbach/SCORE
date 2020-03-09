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


#ifndef __SCORE_MATH_H_INCLUDED__
#define __SCORE_MATH_H_INCLUDED__

//#include <stdlib.h>
#include "SCORE_Random.h"
//#include <vector> //#MK::?

//pure nature
#define _PI_					3.1415926535897932384626433832795

//angle conversions
#define DEG2RAD(deg)			((deg)/(180.0)*(_PI_))
#define RAD2DEG(rad)			((rad)/(_PI_)*(180.0))

//numerics
#define DOUBLE_ACCURACY 		(5.0e-16)				//2^-51 as determined with MATLAB eps(2)
#define MYMATH_STANDARD_SEED	-4256
#define MINIMUM_ANGULAR_SPREAD	(DEG2RAD(1.0))
#define MAXIMUM_MISORI_FCC		(DEG2RAD(62.8))
#define SYMMETRIES_IN_FCC		24

#define SQR(a)					((a)*(a))
#define CUBE(a)					((a)*(a)*(a))
#define MIN(X,Y)				(((X) < (Y)) ? (X) : (Y))
#define MAX(X,Y)				(((X) > (Y)) ? (X) : (Y))


//IPFZ Coloring
#define FAIL_MYMATH_NUMERIC		(-1.0)
#define EPS_PROXIMITY_IPFZ		(0.001)
#define IPF_COLOR_STRETCH_R		(0.5)
#define IPF_COLOR_STRETCH_G		(0.5)
#define IPF_COLOR_STRETCH_B		(0.5)
#define EPS_ENVIRONMENT			(1e-7)
#define RGBRANGE				255


//computation involving primes
#define NPRIMES				32


//Distributional numerical limits
#define RAYLEIGH_MINIMUM_E		(0.1)


//Define which version of the quaternion library
//see: Consistent representations of and conversions between 3D rotations
//D. Rowenhorst, A. D. Rollett, G. S. Rohrer, M. Groeber, M. Jackson, P. J. Konijnenberg, M. de Graef
//Modelling and Simulation in Materials Science and Engineering 2015, Vol 23, doi: 10.1088/0965-0393/23/8/083501
//Fortran-caused index shift though: 0<=>1, 1<=>2, 2<=>3, 3<=>4 but keeping the ijk definition w+ix+jy+kz, w<=>q0, x<=>q1, y<=>q2, z<=>q3
//the comment "DAMASK:ok", means the function is consistent with the one in DAMASK v2.0.1 under the assumption of P_EIJK=+1
//in every case we follow the rotation conventions by H.-J. Bunge, Texture Analysis in Materials Science 1982
//Bunge utilizes the 3-1-3 i.e. ZXZ convention..., gimbal lock occurs when \Phi = 0
//right-handed Cartesian reference frames, passive interpretation of rotations, rotation angle \omega \in [0,\pi]
//counterclockwise rotations are positive to and outbound axis
//the range of the Bunge-Euler angles is \varphi_1 \in [0,2\pi], \Phi \in [0, \pi], \varphi_2 \in [0,2\pi]
//we utilize the ijk quaternion definition
//Quaternions i^2=j^2=k^2=-1, ij=-ji=k, jk = -kj=i, ki=-ik=j, q := (q_0, \textbf{q}) = (q_0, q_1, q_2, q_3)
//all quaternions are meant as unit quaternions!

#define P_EIJK					(+1.0) 		//MK::FOR USE_QUATLIB_DEGRAEF do not change, unless reading above reference! 
#define QUAT2EULER_EPSILON		(1.0e-14)

//choose the respective quaternion library, because there are still differences in their definitions
//YOU MAY DEFINE ONLY ONE OF THOSE AT A TIME!
//#define USE_QUATLIB_DEGRAEF				//following Rowenhorst 2015
//#define USE_QUATLIB_DAMASK					//the implementations utilized in DAMASK v2.0.1
#define USE_QUATLIB_IMM					//utilized with the IMM

//scattering modes
#define REFERENCE_TO_ROTATED_REFERENCE			(0x01)
#define REFERENCE_TO_ROTATED_ARBITRARY_VARIANT	(0x02)

using namespace std;

class randomClass;

struct quat
{
	double q0;
	double q1;
	double q2;
	double q3;
	quat() : q0(1.0), q1(0.0), q2(0.0), q3(0.0) {}
};
typedef struct quat * quatP;


struct euler313
{
	double bunge1;
	double bunge2;
	double bunge3;
	euler313() : bunge1(0.0), bunge2(0.0), bunge3(0.0) {}
};
typedef struct euler313 * euler313P;


class mathMethods
{
public:
	mathMethods( void );
	~mathMethods( void );

	void setprng( long s ); //MK::necessary so that the different ensembleHdl can set explicitly disjoint seeds so that not always the same nucleus orientations are being generated

	//helper LB
	void bubbleSort ( double arr [ ], int size );
	void swap( double& x, double& y );

	//partition a domain
	uint32_t smallestYFactor( uint32_t yz );

	//user-defined modifications to the input data because of different coordinate choices for experiment data and SCORE internal data
	void ebsdori2score( double* in, double* out );

	//general math
	double fac( long x );
	double poissondistribution( double lambda, long k );
	double RayleighDistrPDF( double x, double sigma );
	double RayleighDistrCDF( double x, double sigma );
	double RayleighDistrE2Sigma( double E );

	//Fortran intrinsic translations for USE_QUATLIB_DAMASK
	double sign( double A, double B );
	double huge( double A );
	bool is_DAMASK_NaN( double a, double left, double right );
	double math_limit( double a, double left, double right );

	//quaternion library
	void multiply_quat( double* q1, double* q2, double *q1q2 );
	void multiply_quat_degraef( double *p, double *q, double* pq );

	void quat_on_vec3d( double *q, double *v, double *qv);
	void quat_on_vec3d_degraef( double *p, double *r, double* rp );
	void quat_on_vec3d_imm( double *p, double *r, double* rp );

	void euler2quat( double* bunge, double *q );
	void euler2quat_degraef( double *bunge, double* q );
	void euler2quat_damask( double *bunge, double* q );
	void euler2quat_imm( double *bunge, double* q );

	void quat2euler( double* q, double* bunge );
	void quat2euler_degraef( double *q, double* bunge );
	void quat2euler_damask( double *q, double* bunge );
	void quat2euler_imm( double *q, double* bunge );

	double disori_angle_oriquat_cubic( double* gA, double* gB );
	double disori_angle_oriquat_cubic_degraef( double* p, double* q );
	double disori_angle_oriquat_cubic_imm( double* p, double* q );

	void conjugate_quaternion( double *q );
	void signreversal_quaternion( double *q );
	void normalize_quaternion( double* q, double *uq );
	double rotationangle_quaternion( double* q );
	void rnd_quat_shoemake( double *q );
	void scatter_oriquat( double* qOld, double dismean, double *qNew, char mode );
	void scatter_oriquat_degraef( double* q, double thetamean, double *n, char mode );
	void scatter_oriquat_damask( double* q, double thetamean, double *n, char mode );
	void scatter_oriquat_imm( double* q, double thetamean, double *n, char mode );

	void oriquat2sst_cubic_ipf_quat2pos( double* q1, double* refaxis, double *xy );
	void oriquat2sst_cubic_ipf_pos2rgb( double *pos, unsigned char *rgb,  double ipf_stretch_r, double ipf_stretch_g, double ipf_stretch_b );

	//undocumented
	void disori_fromoriquat_cubic_imm( double* p, double* q, double* quat  );

	//an own PRNG
//private:
	randomClass r;
};


#endif