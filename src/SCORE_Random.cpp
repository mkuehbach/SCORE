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
//PRNG as detailed in Numerical Recipes in C 2nd version and implemented by L. A. Barrales-Mora, some explicit 32-bit datatypes where necessary
//credit to the authors of the additional PRNGs is given in the appropriate code sections

#include "SCORE_Random.h"

void randomClass::init( long sd )
{
	//initialize Park Miller
	if( sd )
		seed = sd;
	else
		seed = (long) time(NULL);

	iy = 0;
	iv = new long[NTAB];

	//initialize Marsaglia's generators
	jsr_shr3 = sd;
	jsr_r4 = sd; //MK::initialize with negative small long will cause an implicit uint32 unroll an hence many bits set

	kn = (uint32_t*) new uint32_t[R4LEN];
	fn = (float*) new float[R4LEN];
	wn = (float*) new float[R4LEN];

	r4_nor_setup();

	//initialize MersenneTwister
	seedMT = (uint32_t) sd;
	mt = (uint32_t*) new uint32_t[NMT];
	mti = NMT+1;
	sgenrand();
	warmupMT( 500000 );
}


//****************************************************************************80

void randomClass::r4_nor_setup ( void )

//****************************************************************************80
//
//  Purpose:
//
//    R4_NOR_SETUP sets data needed by R4_NOR.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 October 2013
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    George Marsaglia, Wai Wan Tsang,
//    The Ziggurat Method for Generating Random Variables,
//    Journal of Statistical Software,
//    Volume 5, Number 8, October 2000, seven pages.
//
//  Parameters:
//
//    Output, uint32_t KN[128], data needed by R4_NOR.
//
//    Output, float FN[128], WN[128], data needed by R4_NOR.
//
{
  int i;
  double dn = R4DN;
  double q;
  double tn = R4TN;

  q = R4VN / exp ( - 0.5 * dn * dn );

  kn[0] = ( uint32_t ) ( ( dn / q ) * R4M1 );
  kn[1] = 0;

  wn[0] = ( float ) ( q / R4M1 );
  wn[127] = ( float ) ( dn / R4M1 );

  fn[0] = 1.0;
  fn[127] = ( float ) ( exp ( - 0.5 * dn * dn ) );

  for ( i = 126; 1 <= i; i-- )
  {
    dn = sqrt ( - 2.0 * log ( R4VN / dn + exp ( - 0.5 * dn * dn ) ) );
    kn[i+1] = ( uint32_t ) ( ( dn / tn ) * R4M1 );
    tn = dn;
    fn[i] = ( float ) ( exp ( - 0.5 * dn * dn ) );
    wn[i] = ( float ) ( dn / R4M1 );
  }

  return;
}


//exemplary cutoff value to skew distribution and avoid negative values
#define R4NOR_LOWERCUTOFF		(0.0)

double randomClass::r4_nor( double mu, double sigma )
{
	//exemplary implementation of a generator to get normally distributed random numbers
	double rr = ( ((double) this->r4_nor()) * sigma ) + mu;
	if ( rr < R4NOR_LOWERCUTOFF ) 
		rr = R4NOR_LOWERCUTOFF;
	return rr;
}

//****************************************************************************80

float randomClass::r4_nor ( void )

//****************************************************************************80
//
//  Purpose:
//
//    R4_NOR returns a normally distributed single precision real value.
//
//  Discussion:
//
//    The value returned is generated from a distribution with mean 0 and
//    variance 1.
//
//    The underlying algorithm is the ziggurat method.
//
//    Before the first call to this function, the user must call R4_NOR_SETUP
//    to determine the values of KN, FN and WN.
//
//    Thanks to Chad Wagner, 21 July 2014, for noticing a bug of the form
//      if ( x * x <= y * y );   <-- Stray semicolon!
//      {
//        break;
//      }
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 July 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    George Marsaglia, Wai Wan Tsang,
//    The Ziggurat Method for Generating Random Variables,
//    Journal of Statistical Software,
//    Volume 5, Number 8, October 2000, seven pages.
//
//  Parameters:
//
//    Input/output, uint32_t &JSR, the seed.
//
//    Input, uint32_t KN[128], data computed by R4_NOR_SETUP.
//
//    Input, float FN[128], WN[128], data computed by R4_NOR_SETUP.
//
//    Output, float R4_NOR, a normally distributed random value.
//
{
  int hz;
  uint32_t iz;
  float value;
  float x;
  float y;

  hz = ( int ) shr3_seeded();
  iz = ( hz & 127 );

  if ( abs ( hz ) < kn[iz] ) //MK::was fabs
  {
    value = ( float ) ( hz ) * wn[iz];
  }
  else
  {
    for ( ; ; )
    {
      if ( iz == 0 )
      {
        for ( ; ; )
        {
          x = - 0.2904764 * log ( r4_uni() );
          y = - log ( r4_uni() );
          if ( x * x <= y + y )
          {
            break;
          }
        }

        if ( hz <= 0 )
        {
          value = - R4R - x;
        }
        else
        {
          value = + R4R + x;
        }
        break;
      }

      x = ( float ) ( hz ) * wn[iz];

      if ( fn[iz] + r4_uni() * ( fn[iz-1] - fn[iz] )
        < exp ( - 0.5 * x * x ) )
      {
        value = x;
        break;
      }

      hz = ( int ) shr3_seeded();
      iz = ( hz & 127 );

      if ( abs ( hz ) < kn[iz] ) //MK::was fabs
      {
        value = ( float ) ( hz ) * wn[iz];
        break;
      }
    }
  }

  return value;
}



//****************************************************************************80

float randomClass::r4_uni ( void )

//****************************************************************************80
//
//  Purpose:
//
//    R4_UNI returns a uniformly distributed real value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 October 2013
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    George Marsaglia, Wai Wan Tsang,
//    The Ziggurat Method for Generating Random Variables,
//    Journal of Statistical Software,
//    Volume 5, Number 8, October 2000, seven pages.
//
//  Parameters:
//
//    Input/output, uint32_t &JSR, the seed.
//
//    Output, float R4_UNI, a uniformly distributed random value in
//    the range [0,1].
//
{
  uint32_t jsr_input;
  float value;

  jsr_input = jsr_r4;

  jsr_r4 = ( jsr_r4 ^ ( jsr_r4 <<   13 ) );
  jsr_r4 = ( jsr_r4 ^ ( jsr_r4 >>   17 ) );
  jsr_r4 = ( jsr_r4 ^ ( jsr_r4 <<    5 ) );

  value = fmod ( 0.5
    + ( float ) ( jsr_input + jsr_r4 ) / 65536.0 / 65536.0, 1.0 );

  return value;
}


//****************************************************************************80

uint32_t randomClass::shr3_seeded ( void )

//****************************************************************************80
//
//#define SHR3 (jsr^=(jsr<<17), jsr^=(jsr>>13), jsr^=(jsr<<5))
//http://www.cse.yorku.ca/~oz/marsaglia-rng.html
//SHR3 is a 3-shift-register generator with period
//2^32-1. It uses y(n)=y(n-1)(I+L^17)(I+R^13)(I+L^5),
//with the y's viewed as binary vectors, L the 32x32
//binary matrix that shifts a vector left 1, and R its
//transpose. SHR3 seems to pass all except those
//related to the binary rank test, since 32 successive
//values, as binary vectors, must be linearly
//independent, while 32 successive truly random 32-bit
//integers, viewed as binary vectors, will be linearly
//independent only about 29% of the time.*/
//  Purpose:
//
//    SHR3_SEEDED evaluates the SHR3 generator for integers.
//
//  Discussion:
//
//    Thanks to Dirk Eddelbuettel for pointing out that this code needed to
//    use the uint32_t data type in order to execute properly in 64 bit mode,
//    03 October 2013.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 October 2013
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    George Marsaglia, Wai Wan Tsang,
//    The Ziggurat Method for Generating Random Variables,
//    Journal of Statistical Software,
//    Volume 5, Number 8, October 2000, seven pages.
//
//  Parameters:
//
//    Input/output, uint32_t &JSR, the seed, which is updated
//    on each call.
//
//    Output, uint32_t SHR3_SEEDED, the new value.
//
{
  uint32_t value;

  value = jsr_shr3;

  jsr_shr3 = ( jsr_shr3 ^ ( jsr_shr3 <<   13 ) );
  jsr_shr3 = ( jsr_shr3 ^ ( jsr_shr3 >>   17 ) );
  jsr_shr3 = ( jsr_shr3 ^ ( jsr_shr3 <<    5 ) );

  value = value + jsr_shr3;

  return value;
}



void randomClass::sgenrand( void )
{
	//Initialization of Mersenne twister
	mt[0] = seedMT & 0xffffffff;
	for( mti=1; mti < NMT; mti++ )
		mt[mti] = (69069 * mt[mti-1]) & 0xffffffff;
}


double randomClass::MersenneTwister( void )
{
/* Saito M, Matsumoto M. SIMD-oriented fast Mersenne twister: a 128-bit
   pseudorandom number generator, Monte-Carlo and Quasi-Monte Carlo Methods 2006; Springer 2008, 2:607
   DOI:10.1007/978-3-540-74496-2_36*/

//Original algorithm interval [0,1]; here modified for [0,1)

	uint32_t y;

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
