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

#ifndef __SCORE_RANDOM_H_INCLUDED__
#define __SCORE_RANDOM_H_INCLUDED__

#include <math.h>
#include <stdint.h>
#include <float.h>
#include <stdlib.h>
#include <time.h>
#include <vector>



//Park Miller constants
#define IA			16807
#define IM			2147483647
#define AM			(1.0/IM)
#define IQ			127773
#define IR			2836
#define NTAB		32
#define NDIV		(1+(IM-1)/NTAB)
#define EPS			(1.2e-7)
#define RNMX		(1.0-EPS)


//Marsaglia's generator constants
#define R4LEN		128
#define R4DN		(3.442619855899)
#define R4M1		(2147483648.0)
#define R4TN		(3.442619855899)
#define R4VN		(9.91256303526217e-03)
#define R4R			(3.442620)


//Mersenne Twister Parameters
#define NMT			624
#define MMT			397
#define MATRIX_A	0x9908b0df
#define UPPER_MASK	0x80000000
#define LOWER_MASK	0x7fffffff

#define TEMPERING_MASK_B 0x9d2c5680
#define TEMPERING_MASK_C 0xefc60000
#define TEMPERING_SHIFT_U(y) (y >> 11)
#define TEMPERING_SHIFT_S(y) (y << 7)
#define TEMPERING_SHIFT_T(y) (y << 15)
#define TEMPERING_SHIFT_L(y) (y >> 18)

#define RANDMAX (1.0-DBL_EPSILON)

#define DEFAULT_SEED	-46356


using namespace std;


class randomClass
{
public:
	randomClass(){
		iv = NULL;
		kn = NULL;
		fn = NULL;
		wn = NULL;
		mt = NULL;
	}
	~randomClass(){
		delete [] iv;
		delete [] kn;
		delete [] fn;
		delete [] wn;
		delete [] mt;
	};

	void init(long sd);
	void initPM(long sd) { seed = sd; }

	void initR4Uni( uint32_t sd ) { jsr_r4 = sd; }
	void initSHR3( uint32_t sd ) { jsr_shr3 = sd; }
	void initMT( uint32_t sd ) { seedMT = sd; sgenrand(); warmupMT( 500000 ); }

	void r4_nor_setup( void );
	double r4_nor( double mu, double sigma );
	float r4_nor( void );
	float r4_uni( void );
	uint32_t shr3_seeded ( void );

	double MersenneTwister( void );
	void sgenrand( void );
	void warmupMT( unsigned int n ) {
		double toss = 0.0;
		for ( unsigned int r = 0; r < n; r++ ) { toss = this->MersenneTwister(); }
	};

private:
		long seed;			//Park-Miller minimal generator
		long iy;
		long* iv;

		uint32_t jsr_shr3;	//seed for Marsaglia's shr3 generator
		uint32_t jsr_r4;	//seed for Marsaglia's r4_uni generator

		uint32_t* kn;		//Ziggurat state set for Marsaglia's Ziggurat standard normal distributed generator
		float* fn;
		float* wn;

		uint32_t seedMT;	//state set for the MersenneTwister
		uint32_t* mt;
		uint32_t mti;
		uint32_t mag01[2] = {0x0, MATRIX_A};
};
typedef randomClass *randomClassP;


#endif