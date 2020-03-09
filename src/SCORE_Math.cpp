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

#include "SCORE_Math.h"
//#include "SCORE_Random.h"


//##LB: Unless indicated, all functions in radians.
//##MK: Quaternion algebra: q = -q define an equivalent rotation.


mathMethods::mathMethods( void )
{
	r.init( MYMATH_STANDARD_SEED );
}


mathMethods::~mathMethods( void )
{
}


void mathMethods::setprng( long s )
{
	this->r.init( s );
}


void mathMethods::bubbleSort ( double arr [ ], int size ) // Sort components of a quaternion in ASCENDING order
{
	int last = size - 2; 
	int isChanged = 1; 

	while ( last >= 0 && isChanged ) {
		isChanged = 0; 
		for ( int k = 0; k <= last; k++ )
			if ( arr[k] > arr[k+1] ) {
				swap ( arr[k], arr[k+1] );
				isChanged = 1; 
			}

		last--;
	}
}


void mathMethods::swap( double& x, double& y ) //Required for the sorting
{
	double temp;
	temp = x;
	x = y;
	y = temp;
}


uint32_t mathMethods::smallestYFactor( uint32_t yz )
{
	//factorizes the unsigned integer yz in such a way into two unsigned integer y and z that y <= z and y * z = yz
	uint32_t z = yz;
	uint32_t y = 1;
	uint32_t test = 0;

	uint32_t Primes[NPRIMES] = {1,2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127};
	for ( uint32_t cand = 0; cand < (NPRIMES-1); cand++ ) { //NPRIMES-1 to remain in range of Primes...
		test = yz / Primes[cand];
		if ( (test * Primes[cand]) == yz ) {
			z = test;
			y = Primes[cand];
		}
		if ( Primes[cand+1] >= z ) break;
	}

	return y;
}


double mathMethods::fac( long x )
{
	if ( x <= 0 )	return 1.0;

	long res = 1;
	for ( long i = 1; i <= x; i++ ) {
		res = res * i;
	}

	return (double) res;
}

double mathMethods::poissondistribution( double lambda, long k )
{
	if ( k < 0 ) return 0.0;

	double res = pow( lambda, k );
	res = res * exp( - lambda );
	res = res / fac(k);

	return res;
}

double mathMethods::RayleighDistrPDF( double x, double sigma )
{
	//evalutes PDF of a Rayleigh distribution with parameter sigma
	double res =  (x < 0.0) ? 0.0 : ( x/SQR(sigma) * exp(-0.5*SQR(x/sigma)) );
	return res;
}

double mathMethods::RayleighDistrCDF( double x, double sigma )
{
	//evaluates CDF of a Rayleigh distribution with parameter sigma
	double res = (x < 0.0) ? 0.0 : (1.0 - exp(-0.5 * SQR(x/sigma)));
	return res;
}

double mathMethods::RayleighDistrE2Sigma( double E )
{
	//evaluates the sigma parameter of a Rayleigh Distribution whose Mean/Expectation Value should be E
	double res = (E < RAYLEIGH_MINIMUM_E) ? RAYLEIGH_MINIMUM_E / sqrt(0.5*_PI_) : E / sqrt(0.5*_PI_);
	return res;
}


void mathMethods::ebsdori2score( double* in, double* out )
{
	//here we implement necessary rotations of coordinate systems and alike
	//in and out are double[3]
	//##MK::accept as they are
	out[0] = in[0];
	out[1] = in[1];
	out[2] = in[2];
}

//some specific C++ equivalences to the DAMASK v2.0.1 Fortran intrinsics calls sign, huge, merge and math_limit
double mathMethods::sign(double A, double B) {
	//C++ equivalent to fortran intrinsic sign(A,B)
	//https://gcc.gnu.org/onlinedocs/gfortran/SIGN.html
	//The kind of the return value is that of A and B. If B\ge 0 then the result is ABS(A), else it is -ABS(A)
	if ( B >= 0.0 )	return fabs(A);
	else 			return (-1.0*fabs(A));
}

double mathMethods::huge(double A) {
	//C++ equivalent to fortran intrinsic sign(A,B)
	//HUGE(X) returns the largest number that is not an infinity in the model of the type of X. 
	//cannot assumed to be equivalent to the DAMASK function if a is not pReal kind 8
	//no type checking here ###
	return 1.8e+308; //std::numeric_limits<double>:: max();
}

//Select values from two arrays according to a logical mask. The result is equal to TSOURCE if MASK is .TRUE., or equal to FSOURCE if it is .FALSE.. 
bool mathMethods::is_DAMASK_NaN( double a, double left, double right ) {
	if ( left > right ) return true;
	else				return false;
}

double mathMethods::math_limit(double a, double left, double right) {
	//QUICKASSERT( left < right ); //see comment in quat2euler_damask in case of gimbal lock
	//C++ equivalent to Fortran DAMASK math_limit(a, left, right) 
	//when present(left) == true && present(right) == true
	//will return NaN if left > right
		/*use prec, only:  DAMASK_NaN
		implicit none
		real(pReal), intent(in) :: a
		real(pReal), intent(in), optional :: left, right
		math_limit = min ( max (merge(left, -huge(a), present(left)), a), merge(right, huge(a), present(right)) )
		if (present(left) .and. present(right)) math_limit = merge (DAMASK_NaN,math_limit, left>right)
		end function math_limit*/
	//merge(left, -huge(a), present(left) == true )
	double arg1 = left;
	//merge(right, huge(a), present(right) == true )
	double arg2 = right;
	double math_limit = MIN( MAX(arg1, a), arg2 );
	//merge( DAMASK_NaN, math_limit, because of left= -1 > right=1 --> false, implicit false, hence
	return math_limit;	
}


void mathMethods::multiply_quat( double* q1, double* q2, double *q1q2 ) {
#ifdef USE_QUATLIB_DEGRAEF
	multiply_quat_degraef( q1, q2, q1q2 );
#endif
#ifdef USE_QUATLIB_DAMASK
	multiply_quat_degraef( q1, q2, q1q2 ); //as it is constitent with deGraef
#endif
#ifdef USE_QUATLIB_IMM
	multiply_quat_degraef( q1, q2, q1q2 ); //as it is constitent with deGraef
#endif
}

void mathMethods::multiply_quat_degraef( double *p, double *q, double* pq ) {
	//Hamilton product of two quaternions, Eq.5 in the Rowenhorst2015
	pq[0] = p[0]*q[0] - p[1]*q[1] -   p[2]*q[2] - p[3]*q[3];
	
	pq[1] = q[0]*p[1] + p[0]*q[1] + (+P_EIJK*p[2]*q[3] - P_EIJK*p[3]*q[2]);
	pq[2] = q[0]*p[2] + p[0]*q[2] + (-P_EIJK*p[1]*q[3] + P_EIJK*p[3]*q[1]);
	pq[3] = q[0]*p[3] + p[0]*q[3] + (+P_EIJK*p[1]*q[2] - P_EIJK*p[2]*q[1]);
	//DAMASK:ok p<=>A, q<=>B
	//IMM:ok q<=>p, p<=>q
}


void mathMethods::quat_on_vec3d( double *q, double *v, double *qv) {
#ifdef USE_QUATLIB_DEGRAEF
	quat_on_vec3d_degraef( q, v, qv );
#endif
#ifdef USE_QUATLIB_DAMASK
	quat_on_vec3d_degraef( q, v, qv ); //DAMASK does not utilize such operations
#endif
#ifdef USE_QUATLIB_IMM
	quat_on_vec3d_imm( q, v, qv ); //##MK::only for backcompatibility, however inconsistent, if not faulty!
#endif
}

void mathMethods::quat_on_vec3d_degraef( double *p, double *r, double* rp ) {
	//rotate vector passively by quaternion p, Eq. 24
	double pre1 = SQR(p[0]) - ( SQR(p[1])+SQR(p[2])+SQR(p[3]) ); //p0^2-\||p\||^2
	double pre2 = 2.0*(p[1]*r[0] + p[2]*r[1] + p[3]*r[2]); //2(pr)
	double pre3 = 2.0*P_EIJK*p[0]; //2P*p0(p x r)
	rp[0] = pre1*r[0] + pre2*p[1] + pre3*(+p[2]*r[3] - p[3]*r[2]);
	rp[1] = pre1*r[1] + pre2*p[2] + pre3*(-p[1]*r[3] + p[3]*r[1]);
	rp[2] = pre1*r[2] + pre2*p[3] + pre3*(+p[1]*r[2] - p[2]*r[1]);	
}

void mathMethods::quat_on_vec3d_imm( double *p, double *r, double* rp ) {
	//MK::only for back-compatibility, however this is an actively interpreted rotation and therefore inconsistent if not incorrect...
	//get rotation matrix components
	double r11 = SQR(p[0]) + SQR(p[1]) - SQR(p[2]) - SQR(p[3]);
	double r12 = 2*(p[1]*p[2] - p[0]*p[3]);
	double r13 = 2*(p[1]*p[3] + p[0]*p[2]);

	double r21 = 2*(p[1]*p[2] + p[0]*p[3]);
	double r22 = SQR(p[0]) - SQR(p[1]) + SQR(p[2]) - SQR(p[3]);
	double r23 = 2*(p[2]*p[3] - p[0]*p[1]);

	double r31 = 2*(p[1]*p[3] - p[0]*p[2]);
	double r32 = 2*(p[2]*p[3] + p[0]*p[1]);
	double r33 = SQR(p[0]) - SQR(p[1]) - SQR(p[2]) + SQR(p[3]);

	//x,y,z vector3d v
	rp[0] = (r11 * r[0]) + (r12 * r[1]) + (r13 * r[2]);
	rp[1] = (r21 * r[0]) + (r22 * r[1]) + (r23 * r[2]);
	rp[2] = (r31 * r[0]) + (r32 * r[1]) + (r33 * r[2]);	
}


void mathMethods::euler2quat( double* bunge, double *q ) {
#ifdef USE_QUATLIB_DEGRAEF
	euler2quat_degraef( bunge, q );
#endif
#ifdef USE_QUATLIB_DAMASK
	euler2quat_damask( bunge, q );
#endif
#ifdef USE_QUATLIB_IMM
	euler2quat_imm( bunge, q );
#endif
}

void mathMethods::euler2quat_degraef( double *bunge, double* q ) {
	//Bunge passive Euler angles (radians) into passive orientation unit quaternion
	double sigma = 0.5*(bunge[0]+bunge[2]);
	double delta = 0.5*(bunge[0]-bunge[2]);
	double c = cos(0.5*bunge[1]);
	double s = sin(0.5*bunge[1]);
	
	q[0] = c*cos(sigma);
	q[1] = -1.0*P_EIJK*s*cos(delta);
	q[2] = -1.0*P_EIJK*s*sin(delta);
	q[3] = -1.0*P_EIJK*c*sin(sigma);
	
	//up to this point consistent with DAMASK only for P_EIJK=+1, however in DAMASK no enforced flipping into upper hemisphere
	
	//sign reversal to flip to upper hemisphere if q[0] is negative
	if ( q[0] < 0.0 )
		signreversal_quaternion( q );
		
	//MK::is not obviously consistent with DAMASK at least because upper hemisphere restriction is not considered	
}

void mathMethods::euler2quat_damask( double *bunge, double* q ) {
	//Bunge passive Euler angles (radians) into passive orientation unit quaternion
	//consistent with DAMASK v2.0.1
	double sigma = 0.5*(bunge[0]+bunge[2]);
	double delta = 0.5*(bunge[0]-bunge[2]);
	
	double c = cos(0.5*bunge[1]);
	double s = sin(0.5*bunge[1]);
	
	q[0] =        c * cos(sigma);
	q[1] = -1.0 * s * cos(delta); //the -1.0 in the vector part convert to passive by implementing the conjugate
	q[2] = -1.0 * s * sin(delta);
	q[3] = -1.0 * c * sin(sigma);
}

void mathMethods::euler2quat_imm( double *bunge, double* q ) {
	//OK, 20130326MK convention: Bunge ZXZ, which represents the (3,1,3) case analyzed in: Diebel 2006 
	//Representing Attitude: Euler Angles, Unit Quaternions, and Rotation Vectors, James Diebel, Stanford University, Stanford, California 94301-9010, Email: diebel@stanford.edu
	//mind utilization of addition theorems as validated with www.wolframalpha.com
	//cos(a+b) = c(a+b) = cacb-sasb
	//cos(a-b) = c(a-b) = cacb+sasb
	//sin(a+b) = s(a+b) = sacb+casb
	//sin(a-b) = s(a-b) = sacb-casb
	//Diebel PSI<=>phi1, theta<=>Phi, PHI<=>phi2
	
	double sigma = (bunge[0]+bunge[2])/2;
	double delta = (bunge[0]-bunge[2])/2;
	
	double c = cos(bunge[1]/2);
	double s = sin(bunge[1]/2);

	//strictly the second component should be p2-p1 however cos(a-b)=cos(b-a) as the addition theorems show...
	//strictly the third component should also be p2-p1 however -sin(a-b)=sin(b-a)
	q[0] = c*cos(sigma);
	q[1] = s*cos(delta);
	q[2] = s*sin(delta);
	q[3] = c*sin(sigma); //applying sin, cos addition theorems
	//MK::obviously this the active interpretation of the rotation quaternion as the conjugate_quaternion ( q ) <=> euler2quat_damask( bunge )
}

void mathMethods::quat2euler( double* q, double *bunge ) {
#ifdef USE_QUATLIB_DEGRAEF
	quat2euler_degraef( q, bunge );
#endif
#ifdef USE_QUATLIB_DAMASK
	quat2euler_damask( q, bunge ); //MK::it is not obvious how this is the same like _degraef, however it is consistent with a passive convention but does for sure not consider antipodal symmetry...
#endif
#ifdef USE_QUATLIB_IMM
	quat2euler_imm( q, bunge ); //MK::only for backward compatibility, division by sP strictly speaking incorrect, however sampling large Euler space in all directions with 0.1 does not find any case where the division by sP or not makes any substantial numerical difference
#endif
}

void mathMethods::quat2euler_degraef( double *q, double* bunge ) {
	//passive orientation unit quaternion in Bunge Euler orientation passive(radians), gimbal lock can occur if \Phi = 0 ...
	double q03 = SQR(q[0]) + SQR(q[3]);
	double q12 = SQR(q[1]) + SQR(q[2]);
	double varkappa = sqrt(q03*q12);
	//q03 and q12 positive so sqrt is positive, however sign bit for small number may be set therefore fabs
	
	if ( fabs(varkappa) > QUAT2EULER_EPSILON ) { //varkappa =/= most likely case, no gimbal lock
		double _varkappa = 1.0/varkappa;
		bunge[0] = atan2( (q[1]*q[3] - P_EIJK*q[0]*q[2]) * _varkappa, (-1.0*P_EIJK*q[0]*q[1] - q[2]*q[3]) * _varkappa );
		bunge[1] = atan2( 2*varkappa, q03 - q12 );
		bunge[2] = atan2( (P_EIJK*q[0]*q[2] + q[1]*q[3]) * _varkappa, (q[2]*q[3] - P_EIJK*q[0]*q[1]) * _varkappa );
		return;
	}
	//not returned? so varkappa numerically == 0
	//else {
		if ( fabs(q12) < QUAT2EULER_EPSILON ) {
			bunge[0] = atan2( -2.0*P_EIJK*q[0]*q[3] , (SQR(q[0]) - SQR(q[3])) );
			bunge[1] = 0.0;
			bunge[2] = 0.0;
		}
		else if ( fabs(q03) < QUAT2EULER_EPSILON ) {
			bunge[0] = atan2( 2.0*q[1]*q[2], (SQR(q[1]) - SQR(q[2])) );
			bunge[1] = _PI_;
			bunge[2] = 0.0;
		}
		else {
			bunge[0] = 0.0; //##MK::std::cout << "ERROR::Quaternion library!" << endl; //##??
			bunge[1] = 0.0;
			bunge[2] = 0.0;
		}
	//}
	//MK::is not obviously consistent with DAMASK at least because upper hemisphere restriction is not considered
}

void mathMethods::quat2euler_damask( double *q, double* bunge ) {
	//passive orientation unit quaternion in Bunge Euler orientation passive(radians), gimbal lock can occur if \Phi = 0 ...
	
	//passive-->active
	q[1] *= -1.0;
	q[2] *= -1.0;
	q[3] *= -1.0;
	
	//since formulas are defined for active rotations
	bunge[1] = acos( 1.0 - 2.0*( SQR(q[1])+SQR(q[2]) ) );
	
	if ( fabs(bunge[1]) < 1.0e-6 ) { //inefficient because least likely
		//the math_limit call has argument left := -1.0 and right := 1.0 so left > right = false, therefore math_limit(q[0],-1.0,1.0) returns not DAMASK_NaN but q[0]
		bunge[0] = sign( 2.0*acos(math_limit(q[0],-1.0,1.0)), q[3] ); //not a hemisphere flipping because angle is allowed in either +-
		bunge[2] = 0.0;
	}
	else {
		bunge[0] = atan2( +1.0*q[0]*q[2] + q[1]*q[3], q[0]*q[1] - q[2]*q[3] );
		bunge[2] = atan2( -1.0*q[0]*q[2] + q[1]*q[3], q[0]*q[1] + q[2]*q[3]);
	}
	
	//ensure correct range
	if ( bunge[0] < 0.0 ) 	bunge[0] += 2.0*_PI_;
	if ( bunge[1] < 0.0 )	bunge[1] += _PI_;
	if ( bunge[2] < 0.0 )	bunge[2] += 2.0*_PI_;
}

void mathMethods::quat2euler_imm( double *q, double* bunge ) {
	//##MK::only provided for back compatibility
	//convention: Bunge, ZXZ, equal to case (3,1,3) as analyzed in Diebel, James, 2006:
	//Representing Attitude: Euler Angles, Unit Quaternions and Rotation Vectors
	//Gimbal lock situation analyzed following the line of Melcher et. al., Conversion of EBSD data by a quaternion based algorithm....
	//TECHNISCHE MECHANIK, 30, 4, (2010), 401 - 413
	double PHI, sP, phi1, phi2;

	double cosPHI = SQR(q[3]) - SQR(q[2]) - SQR(q[1]) + SQR(q[0]);
	double y0 =	+2*q[1]*q[3]	-2*q[0]*q[2]; //as  following Diebel equation 434
	double x0 =	+2*q[2]*q[3]	+2*q[0]*q[1];
	double y1 =	+2*q[1]*q[3]	+2*q[0]*q[2];
	double x1 = -2*q[2]*q[3]	+2*q[0]*q[1]; //but atan2(y,x) yields a domain error when both SQR(x) and SQR(y) <= DOUBLE_ACCURACY zero!

	//this approach works only for properly normalized unit quaternions...

	//acos(x) has goes to +pi for the argument x approaching from right to -1.0 that is when q[3] and q[0] are numerically zero
	if( cosPHI > (1.0 - DOUBLE_ACCURACY) ) cosPHI = 1.0;
	if( cosPHI < (-1.0 + DOUBLE_ACCURACY) ) cosPHI = -1.0;

	//now application of acos function is safe to use... 
	PHI = acos( cosPHI );

	//special case: PHI=0.0, q[1]=q[2]=0 -->cosPHI=1
	//special case:and PHI=_PI_ q[0]=q[3]=0 -->cosPHI=-1
	//in both of which equation 434 would cause atan2(0,0) nevertheless then sin(PHI) is a singularity free indicator of gimbal lock occurs
	sP = sin(PHI);

	//double tmp1 = atan2( y1, x1 );
	//double tmp2 = atan2( y0, x0 );

	//according to WolframAlpha atan2(y/a,x/a)-atan2(y,x) = 0 für a, x, y positive
	//according to the same ref atan2(-y/a,-x/a)-atan2(y/a,x/a)=-pi
	if( SQR(sP) > DOUBLE_ACCURACY ) {
		phi2 = atan2( y0 / sP, x0 / sP ); //##MK::atan2( y0, x0 )
		phi1 = atan2( y1 / sP, x1 / sP ); //##MK::atan2( y1, x1 )
	}
	else {
		//gimbal lock in the case PHI=0.0 first rotation lets ND' || ND0 and next as well so as if at all additive rotation about phi1+phi2 about ND0, choice for either phi1 or phi2 is arbitrary
		phi1 = atan2( 2*(q[1]*q[2] + q[0]*q[3]), (SQR(q[0]) + SQR(q[1]) - SQR(q[2]) - SQR(q[3])) );
		phi2 = 0.0; //arbitrary choice, Rollett for instance partitions equally phi1 and phi2 as in the case PHI = 0 but:: atan2(a12,a11) for phi1 is for instance inconsistent with Diebel

		//more safe but equally heuristical is to make use of Melchers acos(q[0]^2 - q[3]^2) = phi1+phi2 and partition randomly the angular contribution in the case PHI=0 
		//or acos(q[1]^2 - q[2]^2) = phi1 - phi2
	}

	//it always holds the under m-3m that the Euler space is symmetric to 0 <= phi1 <= 2*_PI_ , 0 <= PHI <= _PI_, 0 <= phi2 <= 2*_PI_
	if (phi1 < (0.0 - DOUBLE_ACCURACY) )
		phi1 += 2 * _PI_;
	if (phi2 < (0.0 - DOUBLE_ACCURACY) )
		phi2 += 2 * _PI_;

	bunge[0] = phi1;
	bunge[1] = PHI;
	bunge[2] = phi2; //following the notation order used in Diebel, James 2006	
}


double mathMethods::disori_angle_oriquat_cubic( double* gA, double* gB ) {
	//minimum angle between two passive orientation! quaternions for a cubic crystal
	//following H. Grimmer, Disorientations and Coincidence Rotations for Cubic Lattices, Acta Crystallographica, 1974, 30, 685-688
	//\Delta g = gB * gA^-1 <=> gA gB^-1 because of additional switching symmetry for misorientation quaternions
	//we work in the crystal frame
	//for obtaining the axis one has to test all possible combinations and crystallographically equivalent gA and gB
	//the resulting misorientation with absolute minimum angle and axis in the SST is though not in general unique!
	//see Equations 8 and 9 in above reference

#ifndef USE_QUATLIB_IMM //DEGRAEF or DAMASK
	return disori_angle_oriquat_cubic_degraef( gA, gB );
#else 
	return disori_angle_oriquat_cubic_imm( gA, gB );
#endif
	/*
	//check whether axis is in SST, i.e. a[1] >= a[2] >= a[3] >= 0
	if (a[3] < 0.0) continue;
	//not continued, well then a[3] >= 0.0
	if (a[2] < a[3]) continue;
	//still not continue well then a[2] >= a[3]
	if ( a[1] < a[2]) continue;
	//again not continued, a1>=a2>=a3
	*/
}

double mathMethods::disori_angle_oriquat_cubic_degraef( double* p, double* q ) {
	double qinv[4] = { q[0], q[1], q[2], q[3] };	//observe that the inverse of a unit quaternion is <=> its conjugate

	conjugate_quaternion( qinv );

	double pqinv[4] = { 1.0, 0.0, 0.0, 0.0 }; //Resulting quaternion, rotation of the two previous quaternions pq-1
	multiply_quat_degraef( p, qinv, pqinv );

	//now pqinv is just some misorientation of which there are 24*24*2 possible, we however
	//are interested in the disorientation i.e. a particular but unfortunately not for all cases unique
	//misorientation with the absolute minimal angle, due to Grimmer

	//compute all 24 components of the six characteristic expressions
	double a0 = pqinv[0];
	double a1 = pqinv[1];
	double a2 = pqinv[2];
	double a3 = pqinv[3];
	double f = 1.0 / sqrt( 2.0 );
	
	double expr[24];
	expr[0] = a0;					expr[1] = a1;					expr[2] = a2;					expr[3] = a3;					//Eq. 7a
	expr[4] = f*(a0+a1);			expr[5] = f*(a0-a1);			expr[6] = f*(a2+a3);			expr[7] = f*(a2-a3);			//Eq. 7b
	expr[8] = f*(a0+a2);			expr[9] = f*(a0-a2);			expr[10] = f*(a1+a3);			expr[11] = f*(a1-a3);			//Eq. 7c
	expr[12] = f*(a0+a3);			expr[13] = f*(a0-a3);			expr[14] = f*(a1+a2);			expr[15] = f*(a1-a2);			//Eq. 7d
	expr[16] = 0.5*(a0+a1+a2+a3);	expr[17] = 0.5*(a0+a1-a2-a3);	expr[18] = 0.5*(a0-a1+a2-a3); 	expr[19] = 0.5*(a0-a1-a2+a3); 	//Eq. 7e
	expr[20] = 0.5*(a0+a1+a2-a3);	expr[21] = 0.5*(a0+a1-a2+a3);	expr[22] = 0.5*(a0-a1+a2+a3);	expr[23] = 0.6*(a0-a1-a2-a3);	//Eq. 7f
	
	double omega = 0.0;
	for ( unsigned int i = 0; i < 24; i++ ) {
		double b0 = fabs(expr[i]); //acos(fabs(x)) is minimal for fabs(x) maximal
		if ( b0 >= omega )
			omega = b0;
	}

	//avoid singularity of acos function knowing omega is non negative
	if ( omega < (1.0 - DOUBLE_ACCURACY) ) //most likely case [0.0,1)
		return 2.0*acos(omega);

	//else
	return 0.0;
}

double mathMethods::disori_angle_oriquat_cubic_imm( double* p, double* q ) {
	//##MK::this implementation computes the disorientation as if p and q were defined actively rather than passively!
	//because, according to Bunge the misorientation between two passive orientation matrices
	//i.e. their equivalent unit quaternion representation p and q is defined as p * inv(q) <=> pq^-1
	//switching symmetry for misorientations assures further that pq^-1 <=> qp^-1 BUT HERE...

	double qm1[4];	//Inverse of quaternion q
	qm1[0] = q[0];
	qm1[1] = -1.0 * q[1];
	qm1[2] = -1.0 * q[2];
	qm1[3] = -1.0 * q[3];

	double r[4]; //Resulting quaternion, ...ONE COMPUTES q^-1p, which is 
	//because of the lacking commutativity of quaternion multiplication not pq^-1 it is however only because of 
	//switching symmetry equivalent to p^-1q

	multiply_quat_degraef( qm1, p, r ); 
	//MK:: was p post applied to qm1 but that was not consistent with MTex, in particular not with -(q^-1) = -q0,q1,q2,q3
	//the consistency with MTex was in fact not coincidentally, because, MTex defines all orientation math ACTIVELY,
	//consequently given in MTex two orientations o1 and o2 the following results inv(o1)*o2 <=> inv(o2)*o1 does exactly
	//that
	//##MK::the working hypothesis is that the specifying of MTex Bunge Euler angles pertains only to the choice of 
	//the axis setup i.e. that MTex defines the rotation order with phi_1, PHI, phi_2 via ZXZ however, it does not
	//render the MTex Euler angle numeric value triplets representing Bunge-Euler-angles in passive convention, they
	//are the corresponding transpose <=> inverse!

	//either way, according to Grimmer we now have to determine the smallest angle.
	double r0[6][4];    //There are 12 possible angles

	//Note: The notation r0 is due to the definition of the quaternion which lie
	//in the fundamental zone, this vector possesses the smallest angle, in such a way
	//that r0 is actually the scalar part of this quaternion

	double a = r[0];
	double b = r[1];
	double c = r[2];
	double d = r[3];

	double fac = 1.0 / sqrt( 2.0 ); //0.70710678;

	//six fundamental quaternions
	r0[0][0] =(a+b)*fac; 		r0[0][1]=(a-b)*fac; 		r0[0][2] = (c+d)*fac; 		r0[0][3] = (c-d)*fac;
	r0[1][0] =(a+c)*fac; 		r0[1][1]=(a-c)*fac;	 		r0[1][2] = (b+d)*fac; 		r0[1][3] = (b-d)*fac;
	r0[2][0] =(a+d)*fac; 		r0[2][1]=(a-d)*fac; 		r0[2][2] = (b+c)*fac; 		r0[2][3] = (b-c)*fac;
	r0[3][0] =(a+b+c+d)*0.5; 	r0[3][1]=(a+b-c-d)*0.5; 	r0[3][2] = (a-b+c-d)*0.5; 	r0[3][3] = (a-b-c+d)*0.5;
	r0[4][0] =(a+b+c-d)*0.5; 	r0[4][1]=(a+b-c+d)*0.5; 	r0[4][2] = (a-b+c+d)*0.5; 	r0[4][3] = (a-b-c-d)*0.5;
	r0[5][0] = a;				r0[5][1] = b;				r0[5][2] = c;				r0[5][3] = d;

	double omega = 0.0;
	for(int i = 0; i < 6; i++)
		for( int j=0; j < 4; j++ )
			if( fabs(r0[i][j]) > omega )
				omega = fabs(r0[i][j]);

	//avoid singularity of acos function knowing omega is non negative
	if ( omega < (1.0 - DOUBLE_ACCURACY) ) //most likely case [0.0,1)
		return 2.0*acos(omega);

	//else
	return 0.0;

	//omega = (1.0 - DOUBLE_ACCURACY);
	//old approach Luis:if( omega > 1.0 ) omega = (double) (int) omega;
	//return (2.0*acos(omega));
}



void mathMethods::conjugate_quaternion( double *q ) {
	q[1] *= -1.0;
	q[2] *= -1.0;
	q[3] *= -1.0;
}

void mathMethods::signreversal_quaternion( double *q ) {
	q[0] *= -1.0;
	q[1] *= -1.0;
	q[2] *= -1.0;
	q[3] *= -1.0;
}

void mathMethods::normalize_quaternion( double* q, double *uq ) {
	double qnorm = 1.0 / sqrt( SQR(q[0]) + SQR(q[1]) + SQR(q[2]) + SQR(q[3]) );
	uq[0] *= qnorm;
	uq[1] *= qnorm;
	uq[2] *= qnorm;
	uq[3] *= qnorm;
}

double mathMethods::rotationangle_quaternion( double* q ) {
	if ( fabs(q[0]) < (1.0 - DOUBLE_ACCURACY) ) { //most likely argument safely within range [-1, +1]
		return 2.0*acos(q[0]);
	} 
	else {
		if ( q[0] < 0.0 ) //arg == -1.00000000000000000000
			return _PI_;

		return 0.0; //arg == 1.000000000000 
	}
}


void mathMethods::rnd_quat_shoemake( double *q ) {
	//generates random unit quaternion from the SO3 with K. Shoemake's algorithm
	//without considering symmetry or assuming the unit quaternion to present a passive or active orientation of a crystal!
	//or be restricted to the Northern hemisphere!
	//see: Graphic Gems III (editor D. Kirk) CalTech pp124-134
	//mind quaternion order of Shoemake w  + i*v + j*x + k*z <=> 0 +  1 2 3 with, consistent only when P_EIJK=+1
	double X0 = r.MersenneTwister();
	double X1 = r.MersenneTwister();
	double X2 = r.MersenneTwister();

	double r1 = sqrt(1-X0);
	double r2 = sqrt(X0);
	double theta1 = 2 * _PI_ * X1;
	double theta2 = 2 * _PI_ * X2;

	q[0] = r1 * sin(theta1); //w
	q[1] = r1 * cos(theta1); //v
	q[2] = r2 * sin(theta2); //x
	q[3] = r2 * cos(theta2); //z

	//normalize
	double qnorm = 1.0 / sqrt( SQR(q[0]) + SQR(q[1]) + SQR(q[2]) + SQR(q[3]) );
	q[0] *= qnorm;
	q[1] *= qnorm;
	q[2] *= qnorm;
	q[3] *= qnorm;	
}


void mathMethods::scatter_oriquat( double* qOld, double dismean, double *qNew, char mode ) {
#ifdef USE_QUATLIB_DEGRAEF
	scatter_oriquat_degraef( qOld, dismean, qNew, mode );
#endif
#ifdef USE_QUATLIB_DAMASK
	scatter_oriquat_damask( qOld, dismean, qNew, mode );
#endif
#ifdef USE_QUATLIB_IMM
	scatter_oriquat_imm( qOld, dismean, qNew, mode );
#endif
}

void mathMethods::scatter_oriquat_degraef( double* q, double thetamean, double *n, char mode ) {
	//finds a random quaternion from SO3 and interprets it as a Bunge orientation quaternion, compares with q Bunge orientation quaterion
	//and accept/rejects probabilistically according to a certain distribution (MK::currently Rayleigh)
	double theta, prob, coinFlip;
	double qcand[4] = {1.0, 0.0, 0.0, 0.0}; //identity quaternion
	bool foundvalid = false;

	if ( mode == REFERENCE_TO_ROTATED_ARBITRARY_VARIANT ) {
		while ( foundvalid == false ) {
			rnd_quat_shoemake( qcand ); //just some unit quaternion on the SO3 that is interpreted is at most only a crystallographic variant of q

			theta = disori_angle_oriquat_cubic_degraef( q, qcand ); //how distant from q is it at all?

			prob = RayleighDistrCDF( theta, RayleighDistrE2Sigma( thetamean ) );
			coinFlip = r.MersenneTwister();

			if ( prob <= coinFlip )
				foundvalid = true;
		}
		//sign reversal brings to upper hemisphere, because we interpret this random SO3 variant of q as a scatter passive orientation
		//hence we tacitly assume that all crystallographic variants are equally likely!
		if ( qcand[0] < 0.0 )
			signreversal_quaternion( qcand );

		n[0] = qcand[0];
		n[1] = qcand[1];
		n[2] = qcand[2];
		n[3] = qcand[3];
	}
	else { //mode == REFERENCE_TO_ROTATED_REFERENCE
		//generate random q on SO3 with theta = ||2*acos(b0)|| Rayleigh-distributed with sigma = thetamean, postmultiply resulting rotation quaternion on orientation
		//in such case the resulting orientation has not necessarily the a disorientation of theta to the reference, it becomes
		//spread in a manner which is not independent of the position in orientation
		while ( foundvalid == false ) {
			rnd_quat_shoemake( qcand ); //just some unit quaternion on the SO3 that is interpreted is at most only a crystallographic variant of q
			theta = rotationangle_quaternion( qcand );
			prob = RayleighDistrCDF( theta, RayleighDistrE2Sigma( thetamean ) );
			coinFlip = r.MersenneTwister();
			if ( prob <= coinFlip )
				foundvalid = true;
		}

		double qrot[4] = {1.0, 0.0, 0.0, 0.0};
		multiply_quat_degraef( qcand, q, qrot );
		if ( qrot[0] < 0.0 )
			signreversal_quaternion( qrot );
		n[0] = qrot[0];
		n[1] = qrot[1];
		n[2] = qrot[2];
		n[3] = qrot[3];
	}
}

void mathMethods::scatter_oriquat_damask( double* q, double thetamean, double *n, char mode ) {
	//finds a random quaternion and interprets it as a Bunge orientation quaternion, compares with q Bunge orientation quaterion
	//and accept/rejects probabilistically according to a certain distribution (MK::currently Rayleigh)
	double theta, prob, coinFlip;
	double qcand[4] = {1.0, 0.0, 0.0, 0.0}; //identity quaternion
	bool foundvalid = false;

	if ( mode == REFERENCE_TO_ROTATED_ARBITRARY_VARIANT ) {
		while ( foundvalid == false ) {
			rnd_quat_shoemake( qcand ); //just some unit quaternion on the SO3 that is interpreted is at most only a crystallographic variant of q

			theta = disori_angle_oriquat_cubic_degraef( q, qcand ); //how distant from q is it at all?

			prob = RayleighDistrCDF( theta, RayleighDistrE2Sigma( thetamean ) );
			coinFlip = r.MersenneTwister();

			if ( prob <= coinFlip )
				foundvalid = true;
		}

		n[0] = qcand[0];
		n[1] = qcand[1];
		n[2] = qcand[2];
		n[3] = qcand[3];
	}
	else { //mode == REFERENCE_TO_ROTATED_REFERENCE
		//generate random q on SO3 with theta = ||2*acos(b0)|| Rayleigh-distributed with sigma = thetamean, postmultiply resulting rotation quaternion on orientation
		//in such case the resulting orientation has not necessarily the a disorientation of theta to the reference, it becomes
		//spread in a manner which is not independent of the position in orientation
		while ( foundvalid == false ) {
			rnd_quat_shoemake( qcand ); //just some unit quaternion on the SO3 that is interpreted is at most only a crystallographic variant of q
			theta = rotationangle_quaternion( qcand );
			prob = RayleighDistrCDF( theta, RayleighDistrE2Sigma( thetamean ) );
			coinFlip = r.MersenneTwister();
			if ( prob <= coinFlip )
				foundvalid = true;
		}

		double qrot[4] = {1.0, 0.0, 0.0, 0.0};
		multiply_quat_degraef( qcand, q, qrot );
		n[0] = qrot[0];
		n[1] = qrot[1];
		n[2] = qrot[2];
		n[3] = qrot[3];
	}
}

void mathMethods::scatter_oriquat_imm( double* q, double thetamean, double *n, char mode ) {
	//finds a random quaternion and interprets it as a Bunge orientation quaternion, compares with q Bunge orientation quaterion
	//and accept/rejects probabilistically according to a certain distribution (MK::currently Rayleigh)
	double theta, prob, coinFlip;
	double qcand[4] = {1.0, 0.0, 0.0, 0.0}; //identity quaternion
	bool foundvalid = false;

	if ( mode == REFERENCE_TO_ROTATED_ARBITRARY_VARIANT ) {
		while ( foundvalid == false ) {
			rnd_quat_shoemake( qcand );

			theta = disori_angle_oriquat_cubic_imm( q, qcand );

			prob = RayleighDistrCDF( theta, RayleighDistrE2Sigma( thetamean ) );
			coinFlip = r.MersenneTwister();

			if ( prob <= coinFlip )
				foundvalid = true;
		}

		n[0] = qcand[0];
		n[1] = qcand[1];
		n[2] = qcand[2];
		n[3] = qcand[3];
	}
	else { //mode == REFERENCE_TO_ROTATED_REFERENCE
		//generate random q on SO3 with theta = ||2*acos(b0)|| Rayleigh-distributed with sigma = thetamean, postmultiply resulting rotation quaternion on orientation
		//in such case the resulting orientation has not necessarily the a disorientation of theta to the reference, it becomes
		//spread in a manner which is not independent of the position in orientation
		while ( foundvalid == false ) {
			rnd_quat_shoemake( qcand ); //just some unit quaternion on the SO3 that is interpreted is at most only a crystallographic variant of q
			theta = rotationangle_quaternion( qcand );
			prob = RayleighDistrCDF( theta, RayleighDistrE2Sigma( thetamean ) );
			coinFlip = r.MersenneTwister();
			if ( prob <= coinFlip )
				foundvalid = true;
		}

		double qrot[4] = {1.0, 0.0, 0.0, 0.0};
		multiply_quat_degraef( qcand, q, qrot );
		n[0] = qrot[0];
		n[1] = qrot[1];
		n[2] = qrot[2];
		n[3] = qrot[3];
	}
}



//identify orientation by RGB scheme
void mathMethods::oriquat2sst_cubic_ipf_quat2pos( double* q1, double* refaxis, double *xy ) {
	//returns position of the refaxis in inverse pole figure of a crystal with the passively defined orientation that is represented by unit quaternion q
	//MK::ok, define crystallographic m-3m symmetry operators in quaternion representation
	double qsymm[SYMMETRIES_IN_FCC][4];

	double _sqrt2 = 1.0 / pow ( 2.0, 0.5);
	double half = 0.5;

	//<100>90*degrees four-fold symmetries
	qsymm[0][0] = 1.0;				qsymm[0][1] = 0.0;				qsymm[0][2] = 0.0;				qsymm[0][3] = 0.0;			/*identity*/
	qsymm[1][0] = _sqrt2;			qsymm[1][1] = _sqrt2;			qsymm[1][2] = 0.0;				qsymm[1][3] = 0.0;			/*[100]90*/
	qsymm[2][0] = 0.0;				qsymm[2][1] = 1.0;				qsymm[2][2] = 0.0;				qsymm[2][3] = 0.0;			/*[100]180*/
	qsymm[3][0] = -1.0 * _sqrt2;	qsymm[3][1] = _sqrt2;			qsymm[3][2] = 0.0;				qsymm[3][3] = 0.0;			/*[100]270*/
	qsymm[4][0] = _sqrt2;			qsymm[4][1] = 0.0;				qsymm[4][2] = _sqrt2;			qsymm[4][3] = 0.0;			/*[010]90*/ //MTex -sq2 0 -sq2 0
	qsymm[5][0] = 0.0;				qsymm[5][1] = 0.0;				qsymm[5][2] = 1.0;				qsymm[5][3] = 0.0;			/*[010]180*/ //##MK::be careful, this is consistent with Pomana it is inconsistent with Mtex stating 0,0,-1,0 but q = -q are the same quaternions
	qsymm[6][0] = -1.0 * _sqrt2;	qsymm[6][1] = 0.0;				qsymm[6][2] = _sqrt2;			qsymm[6][3] = 0.0;			/*[010]270*/
	qsymm[7][0] = _sqrt2;			qsymm[7][1] = 0.0;				qsymm[7][2] = 0.0;				qsymm[7][3] = _sqrt2;		/*[001]90*/
	qsymm[8][0] = 0.0;				qsymm[8][1] = 0.0;				qsymm[8][2] = 0.0;				qsymm[8][3] = 1.0;			/*[001]180*/
	qsymm[9][0] = -1.0 * _sqrt2;	qsymm[9][1] = 0.0;				qsymm[9][2] = 0.0;				qsymm[9][3] = _sqrt2;		/*[001]270*/

	//<110>180*degrees two-fold symmetries
	qsymm[10][0] = 0.0;				qsymm[10][1] = _sqrt2;			qsymm[10][2] = _sqrt2;			qsymm[10][3] = 0.0;			/*[110]180*/
	qsymm[11][0] = 0.0;				qsymm[11][1] = _sqrt2;			qsymm[11][2] = -1.0 * _sqrt2;	qsymm[11][3] = 0.0;			/*[1-10]180*/
	qsymm[12][0] = 0.0;				qsymm[12][1] = _sqrt2;			qsymm[12][2] = 0.0;				qsymm[12][3] = _sqrt2;		/*[101]180*/
	qsymm[13][0] = 0.0;				qsymm[13][1] = -1.0 * _sqrt2;	qsymm[13][2] = 0.0;				qsymm[13][3] = _sqrt2;		/*[-101]180*/ //mtex332 0 sq2 0 -sq2
	qsymm[14][0] = 0.0;				qsymm[14][1] = 0.0;				qsymm[14][2] = _sqrt2;			qsymm[14][3] = _sqrt2;		/*[011]180*/ //mtex 0 0 -sq -sq
	qsymm[15][0] = 0.0;				qsymm[15][1] = 0.0;				qsymm[15][2] = -1.0 * _sqrt2;	qsymm[15][3] = _sqrt2;		/*[110]180*/

	//<111>120*degrees, three-fold symmetries
	qsymm[16][0] = half;			qsymm[16][1] = half;			qsymm[16][2] = half;			qsymm[16][3] = half;		/*[111]120*/
	qsymm[17][0] = -1.0 * half;		qsymm[17][1] = half;			qsymm[17][2] = half;			qsymm[17][3] = half;		/*[111]240*/
	qsymm[18][0] = half;			qsymm[18][1] = half;			qsymm[18][2] = -1.0 * half;		qsymm[18][3] = half;		/*[1-11]240*/
	qsymm[19][0] = -1.0 * half; 	qsymm[19][1] = half;			qsymm[19][2] = -1.0 * half;		qsymm[19][3] = half;		/*[1-11]240*/
	qsymm[20][0] = half;			qsymm[20][1] = -1.0 * half;		qsymm[20][2] = half;			qsymm[20][3] = half;		/*[-111]120*/
	qsymm[21][0] = -1.0 * half;		qsymm[21][1] = -1.0 * half;		qsymm[21][2] = half;			qsymm[21][3] = half;		/*[-111]240*/ //mtex h h -h -h
	qsymm[22][0] = half;			qsymm[22][1] = -1.0 * half;		qsymm[22][2] = -1.0 * half;		qsymm[22][3] = half;		/*[-1-11]120*/ //Mtex332-h h h -h
	qsymm[23][0] = -1.0 * half;		qsymm[23][1] = -1.0 * half;		qsymm[23][2] = -1.0 * half;		qsymm[23][3] = half;		/*[-1-11]240*/
//cout << "All m-3m fcc crystal symmetry quaternion operators loaded successfully." << endl;

	double axis[3] = {refaxis[0], refaxis[1], refaxis[2]}; //either RD, TD, ND

	//project to standard triangle by calculating symmetric variants, first ndp
	double triposp[SYMMETRIES_IN_FCC][2];

	for (int s = 0; s < SYMMETRIES_IN_FCC; s++) {
		double qr[4];

		double r[4] = { qsymm[s][0], qsymm[s][1], qsymm[s][2], qsymm[s][3] };

		multiply_quat_degraef( q1, r, qr ); //##why reversed order?

		//ctranspose for unit quaternions, ie. q0, -q1, -q2, -q3  //##why?
		qr[1] *= -1.0;
		qr[2] *= -1.0;
		qr[3] *= -1.0;

		double hbar[3];
		quat_on_vec3d_imm( qr, axis, hbar );


		//##MK::antipodal consistency with MTex3.3.2
		if ( hbar[2] < 0.0 ) {
			hbar[0] *= -1.0;
			hbar[1] *= -1.0;
			hbar[2] *= -1.0;
		}

		//now get the azimuth angle theta and assure no throw out of range exceptions by the acos function, but before range limiting
		if( hbar[2] > (1.0 - DOUBLE_ACCURACY) ) hbar[2]= 1.0;
		if( hbar[2] < (-1.0 + DOUBLE_ACCURACY) ) hbar[2] = -1.0;

		double theta = acos( hbar[2] ); //result is [0.0 <= theta <= pi]

		//the fact that theta goes to pi can cause tan singularities if theta --> pi/2 periodicities
		double tt = tan(theta/2);
		if ( tt > (1.0 - DOUBLE_ACCURACY) ) { tt = 1.0; }

		double psi = 0.0; //capture discontinuity of the atan2 at 0,0
		if ( fabs(hbar[0]) > DOUBLE_ACCURACY && fabs(hbar[1]) > DOUBLE_ACCURACY ) { 
			psi = atan2( hbar[1], hbar[0] ); 
		}

		//psi is limited against pi
		double epps = 1e-10;
		if ( psi > ( _PI_ - epps ) ) { psi = _PI_; }
		if ( psi < ( (-1.0 * _PI_) + epps ) ) { psi = -1.0 * _PI_; }


		triposp[s][0] = tt * cos(psi);
		triposp[s][1] = tt * sin(psi);
	}

	//check which are in the standard triangle and closest to the center in the standard triangle x >= y > 0
	double minnorm = 10.0;
	double norm = 10.0;
	double xmin = FAIL_MYMATH_NUMERIC;
	double ymin = FAIL_MYMATH_NUMERIC;

	//double sstr = pow( 2, 0.5 ) - 1.0;
	//sstr = pow ( sstr + EPS_ENVIRONMENT, 2); //radius of the sst circle

	for	( int s = 0; s < SYMMETRIES_IN_FCC; s++) {
		double x2y2 = SQR(triposp[s][0]) + SQR(triposp[s][1]);
		norm = pow( x2y2, 0.5 );

		if ( norm <= minnorm && triposp[s][0] >= 0.0 && triposp[s][1] >= 0.0 && triposp[s][0] >= triposp[s][1] ) { //&& x2y2 <= sstr ) {
			xmin = triposp[s][0];
			ymin = triposp[s][1];
			minnorm = norm;
		}
	}

	xy[0] = xmin;
	xy[1] = ymin;
}


void mathMethods::oriquat2sst_cubic_ipf_pos2rgb( double *pos, unsigned char *rgb, double ipf_stretch_r, double ipf_stretch_g, double ipf_stretch_b ) {
	//returns position and RGB color value of the refaxis in inverse pole figure of a crystal with the passively defined orientation that is represented by unit quaternion q
	if ( pos[0] == FAIL_MYMATH_NUMERIC || pos[1] == FAIL_MYMATH_NUMERIC ) {
		//color in black, a color otherwise not utilized any mismatch and failure in the function!
		rgb[0] = 0;
		rgb[1] = 0;
		rgb[2] = 0;
		return;
	}

	//heuristic approach to catch numeric cases too close to the IPFZ coloring triangle vertices
	//mind stereographic projection formulas with colatitude angle alpha and in-projection plane anticlockwise location angle phi
	//the center four-fold axis projected alpha = 0, phi = 0
	//(x,y) = R*tan(alpha/2)*[cos(phi),sin(phi)] with R := 1.0 per definition, i.e.
	double xr = 0.0;							double yr = 0.0; //red, (001) fourfold poles 1*tan(0)*[cos(phi),sin(phi)] arbitrarness of phi
	double xg = pow(2.0, 0.5) - 1.0;			double yg = 0.0; //green (011) twofold poles 1*tan(45/2/180*pi)[cos(0),sin(0)]
	double xb = (0.5 * pow(3.0, 0.5)) - 0.5;	double yb = xb; //blue (111) threefold inversion poles
	//cos(alpha) = (0,0,1)*(1,1,1)/(1*sqrt(3)) --> 54.74...., phi = 45...
	double xo = pos[0];					double yo = pos[1];

	//heuristic catch to avoid running in on-the-vertex-location singularities 
	/*guideline implementation of the scaling approach Molodov
	scaling suggestion from K. Molodov HKL and TSL implement different white
	points and different strategies with respect to how to rescale the RGB
	mixing, source unknown so heuristic decision, no problem: because
	an IPF colorcoding for the m-3m, 1 symmetry is possible purpose is to show
	qualitative location, coding necessarily is unpractical for sharp ipfz
	fiber textures in particular if close to the vertices, then a power-law
	rescale might exaggerate contrast the stronger the closer to the vertices
	but the less distinct for orientations which are located near the white point*/

	if ( fabs(xo-xr) < EPS_PROXIMITY_IPFZ && fabs(yo-yr) < EPS_PROXIMITY_IPFZ ) {
		rgb[0] = 255; //assign pure RED
		rgb[1] = 0;
		rgb[2] = 0;
		return;
	}
	if ( fabs(xo-xg) < EPS_PROXIMITY_IPFZ && fabs(yo-yg) < EPS_PROXIMITY_IPFZ ) {
		rgb[0] = 0;
		rgb[1] = 255; //assign pure GREEN
		rgb[2] = 0;
		return;
	}
	if ( fabs(xo-xb) < EPS_PROXIMITY_IPFZ && fabs(yo-yb) < EPS_PROXIMITY_IPFZ ) {
		rgb[0] = 0;
		rgb[1] = 0;
		rgb[2] = 255; //assign pure BLUE
		return;
	}

	//identify RGB color in stereographic triangle
	double k1g = (xg-xo)*(yr-yb) - (yg-yo)*(xr-xb);
	//k1g can only be come 0 if: then the lines are parallel or coincident  xg-xo and yg-yo are either both zero or both terms the same because or naturally yr-yb = xr-xb = const  != 0catch k1g <= DOUBLE_ACCURACY)
	//which based on the shape of the triangle can be safely assumed

	double xgbar = ( (xg*yo - yg*xo)*(xr-xb) - (xg-xo)*(xr*yb - yr*xb) ) / k1g;
	double ygbar = ( (xg*yo - yg*xo)*(yr-yb) - (yg-yo)*(xr*yb - yr*xb) ) / k1g;
	double absggbar = fabs(pow( (SQR(xg-xgbar) + SQR(yg-ygbar)) , 0.5));
	double absogbar = fabs(pow( (SQR(xo-xgbar) + SQR(yo-ygbar)) , 0.5));

	double k1b = (xb-xo)*(yr-yg) - (yb-yo)*(xr-xg);
	double xbbar =( (xb*yo - yb*xo)*(xr-xg) - (xb-xo)*(xr*yg - yr*xg) ) / k1b;
	double ybbar =( (xb*yo - yb*xo)*(yr-yg) - (yb-yo)*(xr*yg - yr*xg) ) / k1b;
	double absbbbar = fabs(pow( (SQR(xb-xbbar) + SQR(yb-ybbar)), 0.5));
	double absobbar = fabs(pow( (SQR(xo-xbbar) + SQR(yo-ybbar)) , 0.5));

	double k1r = pow( (yo/xo), 2.0 );
	//x0 == zero already excluded by R vertex proximity check and the coordiante system choice x>=y
	double xrbar = pow( (2+k1r), 0.5);
	xrbar = xrbar - 1;
	xrbar /= (k1r + 1.0);
	double yrbar = yo/xo * xrbar;
	double absrrbar = fabs(pow( (SQR(xr-xrbar) + SQR(yr-yrbar)), 0.5));
	double absorbar = fabs(pow( (SQR(xo-xrbar) + SQR(yo-yrbar)), 0.5));

	//IPF color stretch approach
	double rrggbb[3];
	rrggbb[0] = pow( (absorbar/absrrbar), ipf_stretch_r );
	rrggbb[1] = pow( (absogbar/absggbar), ipf_stretch_g );
	rrggbb[2] = pow( (absobbar/absbbbar), ipf_stretch_b );

	double maxx = rrggbb[0];
	if ( rrggbb[1] > maxx ) maxx = rrggbb[1];
	if ( rrggbb[2] > maxx ) maxx = rrggbb[2];

	//K. Molodov got a better agrreement with the HKL colo though by anisotropically IPF_COLOR_STRETCHING via reverse engineering
	int rr = rrggbb[0] * (1.0 / maxx) * 255;
	int gg = rrggbb[1] * (1.0 / maxx) * 255;
	int bb = rrggbb[2] * (1.0 / maxx) * 255;

	//"RRGGBB=" << rr << ";" << gg << ";" << bb << endl;

	rgb[0] = rr;
	rgb[1] = gg;
	rgb[2] = bb;
}


//undocumented
void mathMethods::disori_fromoriquat_cubic_imm( double* p, double* q, double* quat  ) {
	//MK::ok
	double qm1[4];    //Inverse of quaternion q

	//Inverse of quaternion q is the same like the conjugate for unit quaternions
	qm1[0] = q[0];
	qm1[1] = -1.0 * q[1];
	qm1[2] = -1.0 * q[2];
	qm1[3] = -1.0 * q[3];

	double r[4]; //Resulting misorientation quaternion, m = pq-1

	multiply_quat_degraef( qm1, p, r ); //MK:: was p post applied to qm1 but that was not consistent with MTex, in particular not with -(q^-1) = -q0,q1,q2,q3


	//Now, we have to determine the smallest angle, following Grimmer H, Acta Cryst 1974, A30 pp685-688
	double r0[6][4];    //There are 12 possible angles

	//Note: The notation r0 is due to the definition of the quaternion which lie
	//in the fundamental zone, this vector possesses the smallest angle, in such a way
	//that r0 is actually the scalar part of this quaternion

	double a = r[0];
	double b = r[1];
	double c = r[2];
	double d = r[3];

	double fac = 1.0 / sqrt( 2.0 ); //0.70710678;

	//six fundamental quaternions
	r0[0][0] =(a+b)*fac; 		r0[0][1]=(a-b)*fac; 		r0[0][2] = (c+d)*fac; 		r0[0][3] = (c-d)*fac;
	r0[1][0] =(a+c)*fac; 		r0[1][1]=(a-c)*fac;	 		r0[1][2] = (b+d)*fac; 		r0[1][3] = (b-d)*fac;
	r0[2][0] =(a+d)*fac; 		r0[2][1]=(a-d)*fac; 		r0[2][2] = (b+c)*fac; 		r0[2][3] = (b-c)*fac;
	r0[3][0] =(a+b+c+d)*0.5; 	r0[3][1]=(a+b-c-d)*0.5; 	r0[3][2] = (a-b+c-d)*0.5; 	r0[3][3] = (a-b-c+d)*0.5;
	r0[4][0] =(a+b+c-d)*0.5; 	r0[4][1]=(a+b-c+d)*0.5; 	r0[4][2] = (a-b+c+d)*0.5; 	r0[4][3] = (a-b-c-d)*0.5;
	r0[5][0] = a;				r0[5][1] = b;				r0[5][2] = c;				r0[5][3] = d;

	double rq[4];

	int mi;
	double max=0.0;

	for( int i=0;i<6;i++ )						//Determing the quaternion with the maximal component and the component itself
		for( int j=0;j<4;j++ )
		{
			if( fabs(r0[i][j]) > max )
			{
				max=fabs(r0[i][j]);
				mi=i;
			}
		}

	rq[0] = fabs( r0[mi][0] );					//Disorientation requires all components positive
	rq[1] = fabs( r0[mi][1] );
	rq[2] = fabs( r0[mi][2] );
	rq[3] = fabs( r0[mi][3] );

	bubbleSort( rq,4 );						//Sorting into ascending order, because a desorientation in the SST
											//requires a quaternion with q0>=q1>=q2>=q3 which represents a minimal 
	quat[0] = rq[3];						//rotation angle and an axis fulfilling h>k>l
	quat[1] = rq[2];
	quat[2] = rq[1];
	quat[3] = rq[0];
	//additionally it is required that rq3 >= sum of all others and rq3 * (sqrt(2)-1) >= rq2
}

