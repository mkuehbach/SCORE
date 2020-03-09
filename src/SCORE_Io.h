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


#ifndef __SCORE_IO_H_INCLUDED__
#define __SCORE_IO_H_INCLUDED__


//STL mainly I/O related includes
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <string.h>


#include <cassert>
//macros relevant to Volkers C IO code only
#define BUFSIZE 1024
#define _STR(x) _VAL(x)
#define _VAL(x) #x
#define ASSERT(cond) if(!(cond)) std::cerr << "failed assertion:" << __FILE__ << " line " << _STR(__LINE__) << ":" << #cond << "\n";
#define TOLFP(x) ((1.0)-(0.99999/((double) x))) //tolerance
//#define ERRTXT(text) (text" : file: "__FILE__" line:"_STR(__LINE__)"\n")
#define stringsNotEqual(a,b) strcmp(a,b)
#define stringsEqual(a,b) !stringsNotEqual(a,b)


#define DEBUG
/*#ifdef DEBUG
	#include <cassert>
	#define QUICKASSERT(condition) if(!(condition)) { std::cerr << "Assertion failed at " << __FILE__ << ":" << __LINE__ << "\n inside " << __FUNCTION__ << "\nCondition: " << condition << "\n"; abort(); }
#else
	#define QUICKASSERT(condition) ((void)0)
#endif*/

//should debug tests should be compiled in
//#define DEBUG

#ifndef DEBUG
		#define QUICKASSERT(cond)       ((void)0)
#else
		#define QUICKASSERT(cond)       ASSERT(cond)
#endif


//compile in the desired output
//#define REPORTSTYLE_DEVELOPER	//detailed function intermediate results
//#define REPORTSTYLE_USER		//only least physical relevant output

//#define REPORTSTYLE_CELLCYCLES


#define READ			0
#define WRITE			1
#define APPEND			2
#define FAILURE			0
#define SUCCESS			1


class dataLine;
class dataBlock;
typedef dataLine * dataLineP;
typedef dataBlock * dataBlockP;

using namespace std;

void reportError( char * message); //int rank = 0);
void exitus (const char * s);


typedef struct
{
	union
	{
		float f;
		long i;
		char *s;
	} d;
	char type;

} univData;
typedef univData * univDataP;


class dataLine
{
public:
	dataLine( void );
	dataLine( int size );
	~dataLine( void );

	univDataP dat;
	long dataCount;

	dataLineP next;
	dataLineP prev;
};


class dataBlock
{
public:
	dataBlock( void );
	~dataBlock( void );
	long columnCount;
	long lineCount;
	
	dataLineP first;
	dataLineP last;

	char head[BUFSIZE];
	char name[128];
};

class io
{

public:
	io( void  );
	~io( void );

	short open( short rw, const char * filename, FILE ** file );
	short open( const char * filename, ofstream * file );
	short open( const char * filename, ifstream file );
	short write( char * message, FILE * file );
	short write( string message, ofstream * file );

	dataBlockP readDataBlock( const char *name, const char *fileName );
	dataBlockP readDataBlock( const char *name, FILE * file );

	char *newString(const char *s);
	double getReal( dataLineP line, long column );
	double getReal2( dataLineP line, long column );
	long getInt( dataLineP line, long column );
	char * getString( dataLineP line, long column );

	double geTReal( const char *s, dataBlockP db );   
	long geTInt( const char *s, dataBlockP db );
	char * geTString( const char *s, dataBlockP db );

};


#endif
