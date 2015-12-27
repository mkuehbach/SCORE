//PRNG as detailed in Numerical Recipes in C 2nd version and implemented by L. A. Barrales-Mora, some explicit 32-bit datatypes where necessary
//SCORE automaton developed by M K{\"u}hbach, 2014/2015, for questions and details contact markus.kuehbach@rwth-aachen.de

#include "SCORE_Random.h"


double randomClass::parkMiller( void ) 
{
	//MK::ok, seeds from -1 to -2^16 tested and acceptable Monobit and spectral pass rate
	//Random Number Generator Park-Miller-Bays-Durham "Numerical Recipes in C Second Edition (ran1), seed value zero is not allowed!
	//Ij+1 = ( a*Ij) % m with M = 2^31 - 1, a = 7^5 giving a period length of [0, m-1] with spectral property for triplets m^(1/k=3) = 1290 at most!
	//period is about 2.1e9
	int32_t j;
	int32_t k;
	static int32_t iy=0;
	static int32_t iv[NTAB];
	double temp;

	if( seed <= 0 || !iy )
	{
		if( -seed < 1 ) seed=1;
		else seed = -seed;
		for( j=NTAB+7;j>=0;j-- )
		{
			k = (int32_t) (seed/IQ);
			seed=(long) (IA*((int32_t)seed-k*IQ)-IR*k);
			if( seed < 0 ) seed += IM;
			if( j < NTAB ) iv[j] = (int32_t) seed;
		}
		iy=iv[0];
	}
	k=(int32_t)(seed/IQ);
	seed=(long)(IA*(seed-k*IQ)-IR*k);
	if( seed < 0 ) seed += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j] = (int32_t) seed;
	temp = AM*iy;

	if( (temp=AM*iy) > RNMX ) return RNMX;
	else return temp;
}


double randomClass::leEcuyer(void) 
{
	// L'Ecuyer combined as defined ran2 in Numerical Recipes in C Second Edition
	// combined period is at most 2.3e18
	//MK::ok, seeds can be negative
	int32_t j=0;
	int32_t k=0;
	static int32_t idum2 = 123456789;
	static int32_t iy = 0;
	static int32_t iv[NTABR2];
	double temp;

	if( seedr2 <= 0 )
	{
		if( -seedr2 < 1 ) seedr2 = 1;
		else seedr2 = -seedr2;

		idum2 = seedr2;

		for( j=NTAB+7; j>=0; j-- ) {
			k = (int32_t) (seedr2 / IQ1R2);
			seedr2 = (long) (IA1R2 * ( seedr2 - k * IQ1R2 ) - k * IR1R2);

			if( seedr2 < 0 ) seedr2 += IM1R2;
			if( j < NTABR2 ) iv[j] = (int32_t) seedr2;
		}
		iy=iv[0];
	}

	k = (int32_t) (seedr2/IQ1R2);
	seedr2 = (long) (IA1R2 * ( seedr2 - k * IQ1R2 ) - k * IR1R2);
	if( seedr2 < 0 ) seedr2 += IM1R2;
	k = (int32_t) (idum2 / IQ1R2);
	idum2 = IA2R2 * (idum2-k*IQ2R2)-k*IR2R2;
	if( idum2 < 0 ) idum2 += IM2R2;
	j = iy/NDIVR2;
	iy=iv[j]-idum2;
	iv[j] = (int32_t) seedr2;

	if( iy < 1 ) iy += IMM1R2;

	if( (temp = AMR2 * iy) > RNMXR2 ) return RNMXR2;
	else return temp;
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
