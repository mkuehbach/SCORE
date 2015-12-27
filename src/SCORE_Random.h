//MK::SCORE is a software developed by Markus K{\"u}hbach in 2014/2015 with the Institute of Physical Metallurgy and Metal Science, RWTH Aachen University
//for questions and comments please contact markus.kuehbach at rwth-aachen.de, last major revision 20150920
//##MK::marks comments that where the code is in rather quick and dirty,i.e. parts of the code which have been validated operative as is but offer certainly room for further performance improvements
//PRNG are the original work from ParkMiller, L'Ecuyer and Matsumoto, and as implemented by L. A. Barrales-Mora from references (see functions), some explicit 32-bit datatypes where necessary

#ifndef __SCORE_RANDOM_H_INCLUDED__
#define __SCORE_RANDOM_H_INCLUDED__

#include <math.h>
#include <stdint.h>
#include <float.h>
#include <stdlib.h>
#include <time.h>

//Mersenne Twister Parameters
#define NMT 624
#define MMT 397
#define MATRIX_A 0x9908b0df
#define UPPER_MASK 0x80000000
#define LOWER_MASK 0x7fffffff

#define TEMPERING_MASK_B 0x9d2c5680
#define TEMPERING_MASK_C 0xefc60000
#define TEMPERING_SHIFT_U(y) (y >> 11)
#define TEMPERING_SHIFT_S(y) (y << 7)
#define TEMPERING_SHIFT_T(y) (y << 15)
#define TEMPERING_SHIFT_L(y) (y >> 18)

//defines constants relevant to the MLCG PRNG
//Park Miller parameter

#define APM (16807)
#define MPM (2147483647)
#define QPM (127773)
#define RPM (2836)
#define RANDMAX (1.0-DBL_EPSILON)

//L'Ecuyer Combined MLCG ran2 Numerical Recipes C

#define M1 (2147483563)
#define A1 (40014)
#define M2 (2147483399)
#define A2 (40692)
#define D1RE (53668)
#define D2RE (52774)
#define M1RE (12211)
#define M2RE (3791)

#define EPS1	0.001
#define EPS2	1.0e-8

#define DEFAULT_SEED	-46356


using namespace std;


class randomClass
{
public:
	randomClass( int32_t sd = DEFAULT_SEED );
	~randomClass( void ){};

	//setting the seeds
	void init( int32_t sd );

	//type of implemented and tested generators
	double parkMiller( void );
	double leEcuyer( void );
	double MersenneTwister( void );
	double sgenrand( void );

	//interval functions
	double randomInterval( double rmin, double rmax );
	void generateShuffledLongArray( long lmin, long lmax, long * result );

	long randomIntervalLongInclusive_parkMiller( long lmin, long lmax );
	long randomIntervalLongInclusive_leEcuyer( long lmin, long lmax );


private:
	int32_t seedOLE1;	//Seed_1 L'Ecuyer 
	int32_t seedOLE2;	//Seed_2 L'Ecuyer
	uint32_t seedMT;	//Seed Mersenne twister
	uint32_t mt[NMT];		//Array for Mersenne
	int32_t seedOPM;	//Seed Park-Miller
	int mti;
};
typedef randomClass *randomClassP;


#endif