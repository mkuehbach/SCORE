//MK::SCORE is a software developed by Markus K{\"u}hbach in 2014/2015 with the Institute of Physical Metallurgy and Metal Science, RWTH Aachen University
//for questions and comments please contact markus.kuehbach at rwth-aachen.de, last major revision 20150920
//##MK::marks comments that where the code is in rather quick and dirty,i.e. parts of the code which have been validated operative as is but offer certainly room for further improvements
//Quaternion based algorithms based on Grimmer, H, Acta Cryst, 1974, conventions according to Bunge and Diebel J, 2006
//Quaternion based algorithms originally implemented by L. A. Barrales-Mora, later revised and validated by M K{\"u}hbach

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


void mathMethods::sortInt( long arr [ ], long size ) // Sort integers
{
	long last = size - 2; 
	long isChanged = 1; 

	while ( last >= 0 && isChanged ) {
		isChanged = 0; 
		for ( long k = 0; k <= last; k++ )
			if ( arr[k] > arr[k+1] ) {
				swapInt ( arr[k], arr[k+1] );
				isChanged = 1; 
			}

		last--;
	}
 }


void mathMethods::swapInt( long& x, long& y ) //Required for the sorting
{
	long temp;
	temp = x;
	x = y;
	y = temp;
}
 

void mathMethods::swap( double& x, double& y ) //Required for the sorting
{
	double temp;
	temp = x;
	x = y;
	y = temp;
}


void mathMethods::sort(int n, double *ra)
{
	int l, j, ir, i;
	float rra;

	l = (n >> 1) + 1;
	ir = n;

	for ( ; ; )
	{
		if (l > 1)					/* still in hiring phase */
			rra = ra[--l];
		else						/* in retirement-and-promotion phase */
		{
			rra = ra[ir];			/* clear a space at end of array */
			ra[ir]=ra[1];			/* retire the top of the heap into it */
			if (--ir == 1) 			/* done with last promotion */
			{
				ra[1] = rra;
				return;
			}						/* end if */
		}							/* end else */
		i = l;						/* whether we are in the hiring phase */
		j = l << 1;					/* or promotion phase, we here set up */
		while ( j <= ir )
		{
			if ( (j < ir) && (ra[j] < ra[j + 1]) )
				++j;				/* compare to the better underling */
				if ( rra < ra[j] )	/* demote rra */
				{
					ra[i] = ra[j];
					j += (i = j);
				}
				else
					j = ir + 1;		/* this is rra's level; set j to */
		}							/* terminate the sift-down */
		ra[i] = rra;				/* put rra into its slot */
	}
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


///////////////////////////////////////////////////////////////////////////
//LB: Unless indicated, all functions in radians.
//MK: Quaternion algebra: q = -q define an equivalent rotation. for unit quaternions ||q|| = 1 so the inverse q^-1 = q*/||q||^2 simplifies to = q* with q* the conjugated q0 - (q1,q2,q3)
void mathMethods::randomOrientationShoemake( double * result )
{
	//K. Shoemake, Graphic Gems III (editor D. Kirk) CalTech pp124-134
	//##MK::mind quaternion order of Shoemake w  + i*v + j*x + k*z <=> 0 +  1 2 3 with 
	double q[4]={0,0,0,0};
	double qnorm;

	double X0 = r.leEcuyer();
	double X1 = r.leEcuyer();
	double X2 = r.leEcuyer();

	double r1 = sqrt(1-X0);
	double r2 = sqrt(X0);
	double theta1 = 2 * _PI_ * X1;
	double theta2 = 2 * _PI_ * X2;

	q[0] = r1 * sin(theta1); //w
	q[1] = r1 * cos(theta1); //v
	q[2] = r2 * sin(theta2); //x
	q[3] = r2 * cos(theta2); //z

	qnorm = sqrt( SQR(q[0]) + SQR(q[1]) + SQR(q[2]) + SQR(q[3]) );

	//normalize
	q[0] = q[0] / qnorm;
	q[1] = q[1] / qnorm;
	q[2] = q[2] / qnorm;
	q[3] = q[3] / qnorm;

	double angles[3] = {0.0};

	quaternion2Euler( q, angles );

	result[0] = angles[0];
	result[1] = angles[1];
	result[2] = angles[2];
}


void mathMethods::randomMisorientationShoemake( double theta, double* qr  )
{
	//##MK::practically possible, but maybe not mathematically sound
	double q[4] = {0.0, 0.0, 0.0, 0.0};

	double qcrit = cos( 0.5 * theta);
	if ( theta < MINIMUM_ANGULAR_SPREAD ) { //limit to assure the while loop to finish
		qcrit = cos( 0.5 * MINIMUM_ANGULAR_SPREAD );
	}

	while( q[0] < qcrit ) {

		double X0 = r.leEcuyer();
		double X1 = r.leEcuyer();
		double X2 = r.leEcuyer();

		double r1 = sqrt(1-X0);
		double r2 = sqrt(X0);
		double theta1 = 2 * _PI_ * X1;
		double theta2 = 2 * _PI_ * X2;

		q[0] = r1 * sin(theta1); //w
		q[1] = r1 * cos(theta1); //v
		q[2] = r2 * sin(theta2); //x
		q[3] = r2 * cos(theta2); //z

		double qnorm = sqrt( SQR(q[0]) + SQR(q[1]) + SQR(q[2]) + SQR(q[3]) );

		//normalize
		q[0] = q[0] / qnorm;
		q[1] = q[1] / qnorm;
		q[2] = q[2] / qnorm;
		q[3] = q[3] / qnorm;
		//definately this algorithm samples random on the SO(3) the resulting quaternion however is not necessarily a disorientation!
	}

	qr[0] = q[0];
	qr[1] = q[1];
	qr[2] = q[2];
	qr[3] = q[3];
}


void mathMethods::multiplyQuaternions( double *q, double* p, double* r )
{
	//MK::ok, mathematically multiplies quaternions q and p, active or passive rotation convention does not matter
	//verified via http://mathworld.wolfram.com/Quaternion.html as well as http://www.mathworks.de/de/help/aeroblks/quaternionmultiplication.html
	//equivalent to the multiplication given in Grimmer, H, 1974 Acta Cryst A30, 685-688
	//mind vector cross product notation qp = q0p0 - qbar*pbar + q0 pbar + p0 qbar + qbar cross pbar with bar vector quantities parts of the quaternion
	r[0] = + q[0] *	p[0]	- q[1] * p[1]	- q[2] *	p[2]	- q[3] *	p[3];
	r[1] = + q[1] *	p[0]	+ q[0] * p[1]	- q[3] *	p[2]	+ q[2] *	p[3];
	r[2] = + q[2] *	p[0]	+ q[3] * p[1]	+ q[0] *	p[2]	- q[1] *	p[3];
	r[3] = + q[3] *	p[0]	- q[2] * p[1]	+ q[1] *	p[2]	+ q[0] *	p[3];
}


void mathMethods::misorientationQuaternionCubic( double* p, double* q, double* quat  )
{
	//MK::ok
	double qm1[4];    //Inverse of quaternion q

	//Inverse of quaternion q is the same like the conjugate for unit quaternions
	qm1[0] = q[0];
	qm1[1] = -1.0 * q[1];
	qm1[2] = -1.0 * q[2];
	qm1[3] = -1.0 * q[3];

	double r[4]; //Resulting misorientation quaternion, m = pq-1

	multiplyQuaternions( qm1, p, r ); //MK:: was p post applied to qm1 but that was not consistent with MTex, in particular not with -(q^-1) = -q0,q1,q2,q3


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


double mathMethods::misorientationCubicQxQ( double q01, double q11, double q21, double q31, double q02, double q12, double q22, double q32 )
{
	//MK::ok
	double p[4] = {q01,q11,q21,q31};
	double q[4] = {q02,q12,q22,q32};

	double qm1[4];    //Inverse of quaternion q
	qm1[0] = q[0];
	qm1[1] = -1.0 * q[1];
	qm1[2] = -1.0 * q[2];
	qm1[3] = -1.0 * q[3];

	double r[4]; //Resulting quaternion, rotation of the two previous quaternions pq-1

	multiplyQuaternions( qm1, p, r ); //MK:: was p post applied to qm1 but that was not consistent with MTex, in particular not with -(q^-1) = -q0,q1,q2,q3

	//Now, we have to determine the smallest angle.

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

	for(int i=0;i<6;i++)
		for( int j=0;j<4;j++ )
			if( fabs(r0[i][j]) > omega )
				omega=fabs(r0[i][j]);


	if( omega > 1.0 ) //avoid singularity of acos function
		omega = (double) (int) omega;

	omega=2*acos(omega);
	//QUICKASSERT( omega <= 1.099 );
	return omega;
}


double mathMethods::misorientationCubic( double pa1, double Pa, double pa2, double pb1, double Pb, double pb2 )
{
	//MK::ok
	double oria[3] = { pa1, Pa, pa2 };
	double orib[3] = { pb1, Pb, pb2 };

	double p[4];
	double q[4];

	euler2quaternion( oria, p );
	euler2quaternion( orib, q );

	double qm1[4];

	//Inverse of quaternion q utilizing unity quaternions
	qm1[0] = q[0];
	qm1[1] = -1.0 * q[1];
	qm1[2] = -1.0 * q[2];
	qm1[3] = -1.0 * q[3];

	double r[4]; //Resulting quaternion, rotation of the two previous quaternions pq-1

	multiplyQuaternions( qm1, p, r ); //MK:: was p post applied to qm1 but that was not consistent with MTex, in particular not with -(q^-1) = -q0,q1,q2,q3
	///#######MK:::::ONLY THIS BRING AGREEMENT WITH MTEX!!!!but not p, qm1 !!!!!!!


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

	double omega = 0.0;

	for(int i=0;i<6;i++)
		for( int j=0;j<4;j++ )
			if( fabs(r0[i][j]) > omega ) //as the sign changes are arbitrary
				omega=fabs(r0[i][j]);

	if( omega > 1 ) //avoid singularity of acos function
		omega = (double) (int) omega;

	omega = 2 * acos(omega);
	//QUICKASSERT( omega <= MAXIMUM_MISORI_FCC );
	return omega;
}


void mathMethods::rotateOrientation( double *oriOri, double angle, double u, double v, double w, double *newOri ) //Angle in radians
{
	//MK::should not be used as not properly validated but with the additional comment is correct?
	double qori[4], qrot[4];
	double _norm = 1.0 / sqrt(SQR(u)+SQR(v)+SQR(w));

	euler2quaternion( oriOri, qori );

	double qref[4] = { cos( 0.5* angle ), u * _norm * sin( 0.5 * angle ), v * _norm * sin( 0.5 * angle ), w * _norm * sin( 0.5 * angle ) };

	//MK::ok, an orientation qref in triclinic space can be used to post rotated another rotation which is parameterized as a quaternion as here qori
	multiplyQuaternions( qref, qori, qrot );

	double euler[3] = {0};

	quaternion2Euler( qrot, euler );

	newOri[0] = euler[0];
	newOri[1] = euler[1];
	newOri[2] = euler[2];
}


void mathMethods::euler2quaternion( double * euler, double * q )
{
	//OK, 20130326MK convention: Bunge ZXZ, which represents the (3,1,3) case analyzed in: Diebel 2006 
	//Representing Attitude: Euler Angles, Unit Quaternions, and Rotation Vectors, James Diebel, Stanford University, Stanford, California 94301-9010, Email: diebel@stanford.edu
	//mind utilization of addition theorems as validated with www.wolframalpha.com 
	//cos(a+b) = c(a+b) = cacb-sasb
	//cos(a-b) = c(a-b) = cacb+sasb
	//sin(a+b) = s(a+b) = sacb+casb
	//sin(a-b) = s(a-b) = sacb-casb

	double p1 = euler[0]; //Diebel PSI
	double t  = euler[1]; //Diebel theta
	double p2 = euler[2]; //Diebel PHI

	double co1 = cos(t/2);
	double s1 = sin(t/2);

	double p[4] = {co1*cos((p1+p2)/2), s1*cos((p1-p2)/2), s1*sin((p1-p2)/2), co1*sin((p1+p2)/2)}; //applying sin, cos addition theorems

	q[0] = p[0];
	q[1] = p[1];
	q[2] = p[2];
	q[3] = p[3];
}


void mathMethods::quaternion2Euler( double * quat, double * euler )
{
	//convention: Bunge, ZXZ, equal to case (3,1,3) as analyzed in Diebel, James, 2006:
	//Representing Attitude: Euler Angles, Unit Quaternions and Rotation Vectors
	//Gimbal lock situation analyzed following the line of Melcher et. al., Conversion of EBSD data by a quaternion based algorithm....
	//TECHNISCHE MECHANIK, 30, 4, (2010), 401 - 413

	double q0 = quat[0];
	double q1 = quat[1];
	double q2 = quat[2];
	double q3 = quat[3];
	double PHI, sP, phi1, phi2;

	double cosPHI = SQR(q3) - SQR(q2) - SQR(q1) + SQR(q0);
	double y0 =	2*q1*q3	-	2*q0*q2; //as  following Diebel equation 434
	double x0 =	2*q2*q3	+	2*q0*q1;
	double y1 =	2*q1*q3	+	2*q0*q2;
	double x1 = -	2*q2*q3	+	2*q0*q1; //but atan2(y,x) yields a domain error when both SQR(x) and SQR(y) <= DOUBLE_ACCURACY zero!

	//this approach works only for properly normalized unit quaternions...

	//acos(x) has goes to +pi for the argument x approaching from right to -1.0 that is when q3 and q0 are numerically zero
	if( cosPHI > (1.0 - DOUBLE_ACCURACY) ) cosPHI = 1.0; 
	if( cosPHI < (-1.0 + DOUBLE_ACCURACY) ) cosPHI = -1.0;

	//now application of acos function is safe to use... 
	PHI = acos( cosPHI );


	//special case: PHI=0.0, q1=q2=0 -->cosPHI=1
	//special case:and PHI=_PI_ q0=q3=0 -->cosPHI=-1
	//in both of which equation 434 would cause atan2(0,0) nevertheless then sin(PHI) is a singularity free indicator of gimbal lock occurs
	sP = sin(PHI);


	if( SQR(sP) > DOUBLE_ACCURACY ) {
		phi2 = atan2( y0 / sP, x0 / sP ); //##atan2( y0, x0 )
		phi1 = atan2( y1 / sP, x1 / sP ); //##atan2( y1, x1 )
	}
	else {
		//gimbal lock in the case PHI=0.0 first rotation lets ND' || ND0 and next as well so as if at all additive rotation about phi1+phi2 about ND0, choice for either phi1 or phi2 is arbitrary
		phi1 = atan2( 2*(q1*q2 + q0*q3), (SQR(q0) + SQR(q1) - SQR(q2) - SQR(q3)) );
		phi2 = 0.0; //arbitrary choice, Rollett for instance partitions equally phi1 and phi2 as in the case PHI = 0 but:: atan2(a12,a11) for phi1 is for instance inconsistent with Diebel

		//more safe but equally heuristical is to make use of Melchers acos(q0^2 - q3^2) = phi1+phi2 and partition randomly the angular contribution in the case PHI=0 
		//or acos(q1^2 - q2^2) = phi1 - phi2
	}


	//it always holds the under m-3m that the Euler space is symmetric to 0 <= phi1 <= 2*_PI_ , 0 <= PHI <= _PI_, 0 <= phi2 <= 2*_PI_
	if (phi1 < (0.0 - DOUBLE_ACCURACY) )
		phi1 += 2 * _PI_;
	if (phi2 < (0.0 - DOUBLE_ACCURACY) )
		phi2 += 2 * _PI_;


	euler[2] = phi2; //following the notation order used in Diebel, James 2006
	euler[1] = PHI;
	euler[0] = phi1;
}


void mathMethods::newOrientationFromReference( double *bunge, double deviation, double *newOri )
{
	//MK::ok
	double qrndmisori[4];
	double qbunge[4], rotated[4];

	euler2quaternion( bunge, qbunge );

	//deviation is internally limited to MINIMUM_ANGULAR_SPREAD if called
	randomMisorientationShoemake( deviation, qrndmisori );

	//applying a rotation on another quaternion parameterized orientation
	multiplyQuaternions( qrndmisori , qbunge, rotated );


	double newEuler[3];
	quaternion2Euler( rotated, newEuler );

	newOri[0] = newEuler[0];
	newOri[1] = newEuler[1];
	newOri[2] = newEuler[2];
}


void mathMethods::specificallyDisorientednewOriFromReference( double *bunge, double sigma_rayl_dis2bunge, double *newOri ) //double sigma_rayl_dis2bunge_max, 
{
	//MK::the conventional newOrientationFromReference applies a random unit quaternion on the SO3 on the quaternion represenation
	//of orientation bunge, clearly and in accordance with MacKenzie one expects now that the distribution of disorientation of newOri to bunge
	//generated with such algorithm to follow similar a MacKenzie distribution, however for modeling nucleation spectra this is not desired
	//because the subgrains, at least in those alloys were they are of relevance for nucleation usually have a distinct orientation spread about
	//the average orientation of the grain

	//so it is in addition necessary to filter all possible generated newOris according to their expectation, otherwise
	//a nucleation model always generates practically highly disoriented nuclei

	//##MK::here we utilize an acceptance-rejection sampling
	double matrix[3];
	matrix[0] = bunge[0];
	matrix[1] = bunge[1];
	matrix[2] = bunge[2];

	double candidate[3];
	double theta, prob, coinFlip;

	bool foundvalid = false;
	//acceptance rejection assuming the subgrains to follow a Rayleigh disorientation distribution
	while ( foundvalid == false ) {
		newOrientationFromReference( matrix, MAXIMUM_MISORI_FCC, candidate );

		theta = this->misorientationCubic( matrix[0], matrix[1], matrix[2], candidate[0], candidate[1], candidate[2] );

		prob = 1.0 - exp ( - 0.5 * SQR(theta/sigma_rayl_dis2bunge) );

		//how likely is theta to occur in Rayleigh distribution?, evaluate pdf Rayleigh with sigma_rayl_dis2bunge
		//prob = theta / SQR(sigma_rayl_dis2bunge) * exp ( -0.5 * SQR((theta/sigma_rayl_dis2bunge)) );

		//to avoid the implementation of calculating Lambert's W function that appears when trying to scale the f(sigma) to its maximum value
		//prob = prob / sigma_rayl_dis2bunge_max;

		//acceptance rejectance scheme
		coinFlip = r.leEcuyer();

		if ( prob <= coinFlip )
			foundvalid = true;
	}

	//report accepted orientation
	newOri[0] = candidate[0];
	newOri[1] = candidate[1];
	newOri[2] = candidate[2];
}



void mathMethods::newOrientationFromReferenceQuat( double *qbunge, double deviation, double *newquat )
{
	//MK::ok, does the same as newOrientationFromReference but overloaded with quaternions
	double qrndmisori[4];

	//this is really a problem when deviation approaches zero
	randomMisorientationShoemake( deviation, qrndmisori );

	//applying a rotation on another quaternion parameterized orientation
	multiplyQuaternions( qrndmisori , qbunge, newquat );
}


void mathMethods::devtorefEuler2RGB ( double *bunge, double *ideal, double maxDev, unsigned char *rgb)
{
	//blue channel stretch from 0.0 to maxDev in radians, all other orientations white
	double dev = misorientationCubic ( bunge[0], bunge[1], bunge[2], ideal[0], ideal[1], ideal[2] );

	if (dev < 0.0 || dev > MAXIMUM_MISORI_FCC ) { //black as error mark
		rgb[0] = 0;
		rgb[1] = 0; 
		rgb[2] = 0;
		return;
	}

	//valid and in range
	if ( dev <= maxDev ) { //more than maxDev_ref larger than maxDev are white
		double colorscale = 255.0 * pow( (dev/maxDev), 1.0 );
		int col = colorscale;

		rgb[0] = col;
		rgb[1] = col;
		rgb[2] = 255;
		return;
	}
	//otherwise white
	rgb[0] = RGBRANGE;
	rgb[1] = RGBRANGE;
	rgb[2] = RGBRANGE;
}


void mathMethods::QuatOnVector3D( double * q, double* v, double* r )
{
	//in accordance with Spieß2009 and Morawiec Pospiech 1989
	double a0 = q[0];
	double a1 = q[1];
	double a2 = q[2];
	double a3 = q[3];

	//get rotation matrix components
	double r11 = SQR(a0) + SQR(a1) - SQR(a2) - SQR(a3);
	double r12 = 2*(a1*a2 - a0*a3);
	double r13 = 2*(a1*a3 + a0*a2);

	double r21 = 2*(a1*a2 + a0*a3);
	double r22 = SQR(a0) - SQR(a1) + SQR(a2) - SQR(a3);
	double r23 = 2*(a2*a3 - a0*a1);

	double r31 = 2*(a1*a3 - a0*a2);
	double r32 = 2*(a2*a3 + a0*a1);
	double r33 = SQR(a0) - SQR(a1) - SQR(a2) + SQR(a3);

	//x,y,z vector3d v
	r[0] = (r11 * v[0]) + (r12 * v[1]) + (r13 * v[2]);
	r[1] = (r21 * v[0]) + (r22 * v[1]) + (r23 * v[2]);
	r[2] = (r31 * v[0]) + (r32 * v[1]) + (r33 * v[2]);
}


void mathMethods::project2fundamentalregion_ipfz( double *qtest, double *xy )
{
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

	double nd[3] = {0.0, 0.0, 1.0}; //z-direction normal vector

	//project to standard triangle by calculating symmetric variants, first ndp
	double triposp[SYMMETRIES_IN_FCC][2];

	for (int s = 0; s < SYMMETRIES_IN_FCC; s++) {
		double qtestqsym[4];

		double qq[4] = { qsymm[s][0], qsymm[s][1], qsymm[s][2], qsymm[s][3] }; //might be superficial

		//QuatOnVector3D( qq, ndrotp, hp ); //no renormalization necessary because pure rotation not stretching
		multiplyQuaternions( qtest, qq, qtestqsym );

		//ctranspose for unit quaternions, ie. q0, -q1, -q2, -q3
		qtestqsym[1] *= -1.0;
		qtestqsym[2] *= -1.0;
		qtestqsym[3] *= -1.0;

		double hbar[3];
		QuatOnVector3D( qtestqsym, nd, hbar );


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


void mathMethods::bunge2ipfz( double phi1, double PHI, double phi2, unsigned char *rgb, double * pos )
{
	//project to fundamental region
	double bunge[3];
	bunge[0] = phi1;
	bunge[1] = PHI;
	bunge[2] = phi2;
	double qbunge[4];
	double position[2];

	euler2quaternion( bunge, qbunge );
//cout << "Euler2Quat=" << qbunge[0] << ";" << qbunge[1] << ";" << qbunge[2] << ";" << qbunge[3] << endl;

	project2fundamentalregion_ipfz( qbunge, position );
//cout << "Position=" << position[0] << ";" << position[1] << endl;

	if ( position[0] == FAIL_MYMATH_NUMERIC || position[1] == FAIL_MYMATH_NUMERIC ) {
		//color in black, a color otherwise not utilized any mismatch and failure in the function!
		rgb[0] = 0;
		rgb[1] = 0;
		rgb[2] = 0;

		pos[0] = position[0];
		pos[1] = position[1];
		return;
	}

	//heuristic approach to catch numeric cases too close to the IPFZ coloring triangle vertices 
	//"red"
	double xr = 0.0;							double yr = 0.0; //red
	double xg = pow(2.0, 0.5) - 1.0;			double yg = 0.0; //green
	double xb = (0.5 * pow(3.0, 0.5)) - 0.5;	double yb = xb; //blue

	double xo = position[0];					double yo = position[1];

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

		pos[0] = position[0];
		pos[1] = position[1];
		return;
	}
	if ( fabs(xo-xg) < EPS_PROXIMITY_IPFZ && fabs(yo-yg) < EPS_PROXIMITY_IPFZ ) {
		rgb[0] = 0;
		rgb[1] = 255; //assign pure GREEN
		rgb[2] = 0;

		pos[0] = position[0];
		pos[1] = position[1];
		return;
	}
	if ( fabs(xo-xb) < EPS_PROXIMITY_IPFZ && fabs(yo-yb) < EPS_PROXIMITY_IPFZ ) {
		rgb[0] = 0;
		rgb[1] = 0;
		rgb[2] = 255; //assign pure BLUE

		pos[0] = position[0];
		pos[1] = position[1];
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
	rrggbb[0] = pow( (absorbar/absrrbar), IPF_COLOR_STRETCH_R );
	rrggbb[1] = pow( (absogbar/absggbar), IPF_COLOR_STRETCH_G );
	rrggbb[2] = pow( (absobbar/absbbbar), IPF_COLOR_STRETCH_B );

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

	pos[0] = position[0];
	pos[1] = position[1];
}
