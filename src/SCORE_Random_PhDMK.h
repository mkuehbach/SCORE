//PRNG as detailed in Numerical Recipes in C 2nd version and implemented by L. A. Barrales-Mora, some explicit 32-bit datatypes where necessary
//SCORE automaton developed by M K{\"u}hbach, 2014/2015, for questions and details contact markus.kuehbach@rwth-aachen.de

#ifndef __SCORE_RANDOM_H_INCLUDED__
#define __SCORE_RANDOM_H_INCLUDED__

#include <math.h>
#include <stdint.h>
#include <float.h>
#include <stdlib.h>
//#include <time.h>



//defines constants relevant to the MLCG PRNG
//Park Miller Bays Durham parameter
#define IA				16807
#define IM				2147483647
#define AM				(1.0/IM)
#define IQ				127773
#define IR				2836
#define NTAB			32
#define NDIV			(1+(IM-1)/NTAB)
#define EPS				DBL_EPSILON
#define RNMX			(1.0-EPS)

//L'Ecuyer Combined MLCG ran2 Numerical Recipes C
#define IM1R2			2147483563
#define IM2R2			2147483399
#define AMR2			(1.0/IM1R2)
#define IMM1R2			(IM1R2-1)
#define IA1R2			40014
#define IA2R2			40692
#define IQ1R2			53668
#define IQ2R2			52774
#define IR1R2			12211
#define IR2R2			3791
#define NTABR2			32
#define NDIVR2			(1+IMM1R2/NTAB)
#define EPSR2			DBL_EPSILON
#define RNMXR2			(1.0-EPS)


#define DEFAULT_SEED	-46356


using namespace std;


class randomClass
{
public:
	randomClass( void ){ seed = DEFAULT_SEED; seedr2 = DEFAULT_SEED; }
	~randomClass( void ){};

	//setting the seeds
	void init( long sd ) { seed = sd; seedr2 = sd; }
	//void init(long sd) { long sd = (long) time(NULL); seed = sd; seedr2 = sd; }

	//type of implemented and tested generators
	double parkMiller( void );
	double leEcuyer( void );

	//interval functions
	double randomInterval( double rmin, double rmax );
	void generateShuffledLongArray( long lmin, long lmax, long * result );

	long randomIntervalLongInclusive_parkMiller( long lmin, long lmax );
	long randomIntervalLongInclusive_leEcuyer( long lmin, long lmax );


private:
	long seed;
	long seedr2;
};
typedef randomClass *randomClassP;


#endif