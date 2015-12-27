//MK::SCORE is a software developed by Markus K{\"u}hbach in 2014/2015 with the Institute of Physical Metallurgy and Metal Science, RWTH Aachen University
//for questions and comments please contact markus.kuehbach at rwth-aachen.de, last major revision 20150920
//##MK::marks comments that where the code is in rather quick and dirty,i.e. parts of the code which have been validated operative as is but offer certainly room for further improvements
//Quaternion based algorithms based on Grimmer, H, Acta Cryst, 1974, conventions according to Bunge and Diebel J, 2006
//Quaternion based algorithms originally implemented by L. A. Barrales-Mora, later revised and validated by M K{\"u}hbach

#ifndef __SCORE_MATH_H_INCLUDED__
#define __SCORE_MATH_H_INCLUDED__

//#include <stdlib.h>
#include "SCORE_Random.h"


//numerics
#define DOUBLE_ACCURACY 		5e-16				//2^-51 as determined with MATLAB eps(2)
#define MYMATH_STANDARD_SEED	-4256
#define MINIMUM_ANGULAR_SPREAD	(0.017453292)		//1.0/180.0_PI_
#define MAXIMUM_MISORI_FCC		(1.09606677)		//FCC 62.8/180.0*_PI_
#define SYMMETRIES_IN_FCC		24
#define _PI_					3.1415926535897932384626433832795
#define SQR(a)					((a)*(a))
#define CUBE(a)					((a)*(a)*(a))
#define MIN(X,Y)				(((X) < (Y)) ? (X) : (Y))
#define MAX(X,Y)				(((X) > (Y)) ? (X) : (Y))


//IPFZ Coloring
#define FAIL_MYMATH_NUMERIC		(-1.0)
#define EPS_PROXIMITY_IPFZ		(0.01)
#define IPF_COLOR_STRETCH_R		(0.5)
#define IPF_COLOR_STRETCH_G		(0.5)
#define IPF_COLOR_STRETCH_B		(0.5)
#define EPS_ENVIRONMENT			(1e-7)
#define RGBRANGE				255




using namespace std;

class randomClass;


class mathMethods
{
public:
	mathMethods( void );
	~mathMethods( void );

	void setprng( long s ); //MK::necessary so that the different ensembleHdl can set explicitly disjoint seeds so that not always the same nucleus orientations are being generated

	//helper LB
	void bubbleSort ( double arr [ ], int size );
	void sortInt( long arr [], long size );
	void swap ( double& x, double& y );
	void swapInt ( long& x, long& y );
	void sort(int n, double *ra);

	//general math
	double fac( long x );
	double poissondistribution( double lambda, long k );

	//quaternion algebra
	void multiplyQuaternions( double *q, double* p, double* r );
	void QuatOnVector3D( double *q, double* v, double* r );

	//convert orientations among parametrizations, such as Bunge (3,1,3) rotation convention and quaternion
	void euler2quaternion( double * euler, double * q );
	void quaternion2Euler( double * quat, double * euler );


	//calculate disorientation among two orientations in various parametrizations
	double misorientationCubic( double pa1, double Pa, double pa2, double pb1, double Pb, double pb2 );
	void misorientationQuaternionCubic( double* p, double* q, double* quat  );
	double misorientationCubicQxQ( double q01, double q11, double q21, double q31, double q02, double q12, double q22, double q32 );


	//generate new orientations with some scatter
	void randomOrientationShoemake( double * result );
	void randomMisorientationShoemake( double theta, double* qr );
	
	void rotateOrientation( double *oriOri, double angle, double u, double v, double w, double *newOri );
	void newOrientationFromReference( double *oriOri, double deviation, double *newOri );
	void specificallyDisorientednewOriFromReference( double *bunge, double sigma_rayl_dis2bunge,  double *newOri ); //double sigma_rayl_dis2bunge_max,
	void newOrientationFromReferenceQuat( double *qbunge, double deviation, double *newquat );

	//identify orientation by RGB scheme
	void devtorefEuler2RGB ( double *bunge, double *ideal, double maxDev, unsigned char *rgb); //blue channel stretch from 0.0 to maxDev in radian, all other orientations white
	void project2fundamentalregion_ipfz ( double *qtest, double *xy );
	void bunge2ipfz( double phi1, double PHI, double phi2, unsigned char *rgb, double * pos);

	//an own PRNG
//private:
	randomClass r;
};


#endif