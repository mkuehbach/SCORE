//MK::SCORE is a software developed by Markus K{\"u}hbach in 2014/2015 with the Institute of Physical Metallurgy and Metal Science, RWTH Aachen University
//for questions and comments please contact markus.kuehbach at rwth-aachen.de, last major revision 20150920
//##MK::marks comments that where the code is in rather quick and dirty,i.e. parts of the code which have been validated operative as is but offer certainly room for further improvements
//PRNG are the original work from ParkMiller, L'Ecuyer and Matsumoto, and as implemented by L. A. Barrales-Mora from references (see functions), some explicit 32-bit datatypes where necessary


#include "SCORE_Random.h"

randomClass::randomClass(int32_t sd) 
{
	//Constructor calls always init but init can be called outside to re-initialize the seeds if desired
	init( sd );
}

void randomClass::init( int32_t sd )
{
	int32_t seed = DEFAULT_SEED; //##MK::was uninitialized

	seed = sd;
	
	//if( sd == DEFAULT_SEED ) seed = (int32_t) time( NULL ); //must be commented out in an MPI world otherwise all processes get different seeds depending on when they become executed!

	mti = NMT+1;

	//Park Miller
	seedOPM = seed;
	seedMT = (uint32_t) seed;
	
	//L'Ecuyer combination of two MLG PRNGs
	seedOLE1 = seed;  		//MK::was seedMT
	seedOLE2 = sd-1;		//MK::was (uint32_t) time(NULL);

	if( seedOPM < 0 ) seedOPM += MPM;

	sgenrand();
}

long randomClass::randomIntervalLongInclusive_parkMiller( long lmin, long lmax ) 
{
	//MK::ok
	if( lmax < lmin ) return 0;
	lmax++; //making lmax inclusive

	double r = parkMiller();

	return (long) ( lmin + r * ( lmax - lmin ) );
}


long randomClass::randomIntervalLongInclusive_leEcuyer( long lmin, long lmax )
{
	//MK::ok
	if( lmax < lmin ) return 0;
	lmax++; //making lmax inclusive

	double r = leEcuyer();

	return lmin + ( long ) (r * ( lmax - lmin  ) );
}


double randomClass::randomInterval( double rmin, double rmax ) 
{
	//MK::ok
	if( rmax < rmin ) return 0.0;
	double r = parkMiller();
	return rmin + r * ( rmax - rmin );
}


void randomClass::generateShuffledLongArray( long lmin, long lmax, long * result ) 
{
	//MK::ok
	if( lmax <= lmin ) { result = NULL; return; }

	for( long i = 0; i < ( lmax - lmin ); i++ ) result[i] = i + lmin;

	for( long i = 0; i < ( lmax - lmin ); i++ )
	{
		//swop element at a random pos with the one at i successively
		long pos = i + randomIntervalLongInclusive_parkMiller( 0,  lmax - lmin - i - 1 );
		long temp = result[i];
		result[i] = result[pos];
		result[pos] = temp;
	}
}

double randomClass::parkMiller( void )
{
	int32_t lo, hi, test;

	hi = seedOPM / QPM;
	lo = seedOPM - QPM * hi;

	test = APM * lo - RPM * hi;

	if( test > 0 )
		seedOPM = test;
	else
		seedOPM = test + MPM;

	double randr = ((double) seedOPM) / MPM;

	if( randr > RANDMAX ) return RANDMAX; //MK::was rand, typo?

	return randr;
}

double randomClass::sgenrand( void )	//Initialization of Mersenne twister
{
	mt[0] = seedMT & 0xffffffff;
	for( mti=1; mti<NMT; mti++ )
		mt[mti] = (69069 * mt[mti-1]) & 0xffffffff;
}

double randomClass::MersenneTwister( void )
{
/* Saito M, Matsumoto M. SIMD-oriented fast Mersenne twister: a 128-bit
   pseudorandom number generator, Monte-Carlo and Quasi-Monte Carlo Methods
   2006;2:607 */

//Original algorithm interval [0,1]; here modified for [0,1)

	uint32_t y;
	static uint32_t mag01[2]={0x0, MATRIX_A};

	if( mti >= NMT )
	{
		int kk;

		if( mti == NMT+1 )
		{
			seedMT = (uint32_t) time(NULL);
			sgenrand();
		}

		for( kk=0; kk<NMT-MMT; kk++ )
		{
			y = ( mt[kk] & UPPER_MASK ) | ( mt[kk+1] & LOWER_MASK );
			mt[kk] = mt[kk+MMT] ^ ( y >> 1 ) ^ mag01[y & 0x1];
		}

		for(;kk<NMT-1;kk++)
		{
			y = ( mt[kk] & UPPER_MASK ) | ( mt[kk+1] & LOWER_MASK );
			mt[kk] = mt[kk+(MMT-NMT)] ^ (y >> 1) ^ mag01[y & 0x1];
		}

		y = ( mt[NMT-1] & UPPER_MASK ) | ( mt[0] & LOWER_MASK );
		mt[NMT-1] = mt[MMT-1] ^ (y >> 1) ^ mag01[y & 0x1];

		mti = 0;
	}

	y = mt[mti++];
	y ^= TEMPERING_SHIFT_U(y);
	y ^= TEMPERING_SHIFT_S(y) & TEMPERING_MASK_B;
	y ^= TEMPERING_SHIFT_T(y) & TEMPERING_MASK_C;
	y ^= TEMPERING_SHIFT_L(y);
	
	double randr = ( (double) y / (uint32_t)0xffffffff );
	
	if( randr > RANDMAX )	return RANDMAX;

	return randr;
}

double randomClass::leEcuyer( void )
{

/* L'Ecuyer P. Efficient and portable combined random number generators, Communications of the ACM 1998;31:742 */
// Interval (0,1)

	int32_t Z, k;

	while( seedOLE1 > M1-1 )
	{
		seedOLE1 = seedOLE1 / 2;
	}

	while( seedOLE2 > M2-1 )
	{
		seedOLE2 = seedOLE2 / 2;
	}

	k = seedOLE1 / D1RE;
	seedOLE1 = A1 * (seedOLE1 - k * D1RE) - k * M1RE;

	if( seedOLE1 < 0 )	seedOLE1 += M1;

	k = seedOLE2 / D2RE;
	seedOLE2 = A2 * (seedOLE2 - k * D2RE) - k * M2RE;

	if( seedOLE2 < 0 )	seedOLE2 += M2;

	Z = seedOLE1 - seedOLE2;

	if( Z < 1 ) Z += (M1-1);

	return Z * (1.0/M1);
	
}
