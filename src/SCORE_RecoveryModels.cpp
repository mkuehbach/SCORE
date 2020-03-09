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

#ifndef __SCORE_KERNEL_H_INCLUDED__
//#define __SCORE_KERNEL_H_INCLUDED__

#include "SCORE_Kernel.h"

void caHdl::nes_networkgrowthmodel_vacancycorediff( void )
{
	//MK::template code for a recovery according to Nes vacancy diffusion controlled themrally-activated glide
	//MK::Euler forward difference scheme with dr/dt ~ r(t+dt) - r(t) / dt
	//model parameter at this->t
	double Bx;
	double B = _PI_ * myRecoveryModel.NesC3 * myRecoveryModel.NesC5 * myRecoveryModel.VacancyDiffD0 * exp ( -1.0 * myRecoveryModel.VacancyDiffHact / ( kboltzman * CurrentTemperature ));
	B /= myRecoveryModel.NesKappa2;
	double exponent = SQR(CurrentBurgersVector) / (kboltzman * CurrentTemperature);
	double expterm = exp( myRecoveryModel.NesAlpha3 * myRecoveryModel.NesKappa2 * CurrentG * CurrentBurgersVector * exponent );

	uint32_t ndg = mydefgpool.size();
	double r_tdt, r_t;
	double currdt = this->dt;
	double currrhomax = RHOMAX_WELLANNEALED;

	for ( uint32_t dg = 0; dg < ndg; dg++ ) {
		r_t = 1.0 / pow( mydefgpool[dg].rho, 0.5 ); //>1.0

		Bx = (B / r_t) * expterm;

		r_tdt = (Bx * currdt) + r_t;

		mydefgpool[dg].rho = 1.0 / SQR(r_tdt);

		if ( mydefgpool[dg].rho > currrhomax ) 
			currrhomax = mydefgpool[dg].rho;
	}

	//this->myrhomax = currrhomax; //##MK::not necessary any longer because of update_myrhomax
}


void caHdl::nes_networkgrowthmodel_solutedrag( void )
{
	//#MK::not implemented yet
}


void caHdl::nes_michalakpaxton( void ) 
{
//set of parameter determined as in Michalak Paxton 1961, TransAIME, 221, 1961 p853
//and Nes, Recovery revisited
	//double C5 = 2.0/3.0;
	//double kb = 8.61733e-5;
	//double k3a3 = 0.1;
	//double Usd = 2.90;
	//double atau2 = 6.92878e-13;
	//double btau5 = 4.72858e-14;
	double f1 = 0.25;
	double f2 = 0.75;

	//double temp = CurrentG * CurrentBurgersVector;
	//double tau2 = pow( (atau2 * exp( ((C5 * Usd) - (k3a3 * temp * SQR(CurrentBurgersVector))) / (kb*CurrentTemperature) )), -1.0);
	//double tau5 = pow( (btau5 * exp( (C5 * Usd)/(kb*CurrentTemperature) )), -1.0);
	double tau2 = 0.3605; //min
	double tau5 = 5.2824; //min

	uint32_t ndg = mydefgpool.size();
	//double currrhomax = RHOMAX_WELLANNEALED;
	double currtime_minutes = this->t / 60.0;

	double factor = f1 * pow( (1 + ( currtime_minutes / tau2 )), -1.0 ) + f2 * pow( (1 + ( currtime_minutes / tau5 )), -0.5 );

	//"recovery modeling" in each grain
	for ( uint32_t dg = 0; dg < ndg; dg++ ) {
		mydefgpool[dg].rho = mydefgpool[dg].rho0 * SQR(factor);

		//if ( mydefgpool[dg].rho > currrhomax ) 
		//	currrhomax = mydefgpool[dg].rho;
	}

	//this->myrhomax = currrhomax; //##MK::not necessary any longer because of update_myrhomax

cout << "Michalak Paxton Recovery;tau2;tau5;factor;currrhomax\t\t" << this->t << ";" << this->X << ";" << tau2 << ";" << tau5 << ";" << factor << endl; //";" << currrhomax << endl;
}


void caHdl::update_mydgoptimize_rho( void ) {
	//ONLY CALL FROM ONE THREAD
	uint32_t ndg = mydgoptimize.size();
	for ( uint32_t dg = 0; dg < ndg; dg++ ) {
		uint32_t i = mydgoptimize[dg].id; //remote misses
		mydgoptimize[dg].rho = mydefgpool[i].rho; //mydefgpool becomes never sorted
	}
	std::sort( mydgoptimize.begin(), mydgoptimize.end(), SortDoubleDescending );
}

void caHdl::update_mydgoptimize_cnt( void ) {
	//ONLY CALL FROM ONE THREAD
	uint32_t ndg = mydgoptimize.size();
	for ( uint32_t dg = 0; dg < ndg; dg++ ) { //##MK::utilize that  array thins out
		uint32_t i = mydgoptimize[dg].id; //remote misses because order of ...[dg].id is not strictly increasing as is for mydefgpool
		mydgoptimize[dg].cnt = mydefgpool[i].cellcount; //mydefgpool becomes never sorted
	}
}

void caHdl::update_myrhomax( void ) {
	uint32_t ndg = mydgoptimize.size();
	double rmax = RHOMAX_WELLANNEALED;
	//scan control structure of how much volume with specific dislocation density
	//is still left, the array is sorted in descending order so the first values give the highest 
	uint32_t dg = 0;
	for ( dg = 0; dg < ndg; dg++ ) {
		if ( mydgoptimize[dg].cnt > 0 ) { //is there a 
			if ( mydgoptimize[dg].rho >= rmax ) {
				rmax = mydgoptimize[dg].rho;
				break;
			}
		}
	}

	if ( this->myensRank == MASTER ) std::cout << "Updated myrhomax/dg/cnt now " << rmax << "\t\t" << dg << "\t\t" << mydgoptimize[dg].cnt << std::endl;

	this->myrhomax = rmax;
}

unsigned char caHdl::rho2grayscale( double r ) {
	//utilizing the mapping defined in log_initrho_scalebar

	if ( r < (myrhomin0 + DOUBLE_ACCURACY) )
		return UCHAR_RANGE_MIN; //smallest, dark 0x00
	if ( r > (myrhomax0 - DOUBLE_ACCURACY) ) 
		return UCHAR_RANGE_MAX; //highest, white 0xFF
	
	//double sc = 255.0 * (1.0 - (log10(r) - log10(myrhomin0))/(log10(myrhomax0) - log10(myrhomin0)) );
	//MK::problem with this scaling is that bin widths increase exponentially towards the most decisive rho values
	//MK::given for instance a population with 10^14 and 10^16 a GB of a nucleus with fixed m will migrate 100x faster seeing 10^16 than 10^14
	//so we scale that we get most sensitive for the high values
	double sc = r/myrhomax0 * 255.0;
	return (unsigned char) sc; //smallest value to get is 1/255*myrhomax0
}


void caHdl::log_initrho_scalebar( void ) {
	//MK::mapping of dislocation density values on [RHO_ on solitary unit domain dependent grayscale bar via
	double rmax = RHO_DEFORMED_MATERIAL_MAX;
	if ( this->mydgoptimize.size() > 0 ) //if existent take actual value
		rmax = mydgoptimize.at(0).rho;

	double rmin = RHO_RECRYSTALLIZED_MATERIAL;
	if ( rmin < DOUBLE_ACCURACY ) 
		rmin = 0.0;

	myrhomax0 = rmax;
	myrhomin0 = rmin;

	//thus defining the interval [rmin,rmax]
	//rho2grayscale double sc = 255.0 * (1.0 - ((log10(rho)-log10(rmin)) / (log10(rmax)-log10(rmin)))) if sc > 255.0-DOUBLE-ACCURACY -->0xFF, if sc < DOUBLE_ACCURACY -->0x00, else (char) sc
	//grayscale2rho = 10^[(255.0 - (double) grayscale)/255.0 * (log10(rmax)-log10(rmin)) + log10(rmin)]
	//double diff = log10(rmax) - log10(rmin);
	//double offset = log10(rmin);
	cout << "Dislocation density grayscale bar values" << endl;

	//scale bars may differ between the CAs because they are adapted to the maximum dislocation content
	ofstream scalebar;
	string scalebar_fn = "SCORE." + std::to_string(this->myensHdl->simid) + "CA." + std::to_string(this->jobid) + ".DGRHOGRYScalebar.csv";
	scalebar.open ( scalebar_fn.c_str() );

	scalebar << "//SCORE auto-generated scalebar for colormodel COLORIZE_RX_RHOGREYSCALE_DEFORMED!\n";
	scalebar << "//SCORE rhomin = " << setprecision(8) << rmin << "\n";
	scalebar << "//SCORE rhomax = " << setprecision(8) << rmax << "\n";
	scalebar << "\n\n";
	//scalebar << "unsigned char value;lg(rho/m^-2);rho/m^-2\n";
	scalebar << "unsigned char value;<=rho/m^-2\n";


	for ( unsigned int c = 1; c < 256; c++ ) {
		//double cc = (((double) c) / 255.0 * diff) + offset;
		//scalebar << (int) c << ";" << setprecision(8) << cc << ";" << setprecision(8) << pow( 10.0, cc ) << "\n";
	
		double cc = (double) c / 255.0 * rmax;
		scalebar << (int) c << ";" << setprecision(8) << cc<< "\n";
	}

	scalebar.flush();
	scalebar.close();

	//render the visual colorbar bottom is light == low values, upper is dark == high values
	uint32_t imgx = 200;
	uint32_t imgy = 2550;
	unsigned char* img_rgba = NULL;
	img_rgba = new unsigned char [4*imgx*imgy];
	QUICKASSERT( img_rgba != NULL );
	uint32_t c = 0;
	double sc = 0;
	unsigned char uc = 0x00;
	for ( uint32_t y = 0; y < imgy; y++ ) {
		//sc = 255.0 - (0.1*((double) y));
		sc = 0.1*((double) y);
		QUICKASSERT( sc >= 0.0 && sc <= 255.0 );
		uc = (unsigned char) sc;

		for ( uint32_t x = 0; x < imgx; x++ ) {
			img_rgba[c+REDCHAN] = uc;
			img_rgba[c+GREENCHAN] = uc;
			img_rgba[c+BLUECHAN] = uc;
			img_rgba[c+ALPHACHAN] = UCHAR_RANGE_MAX;
			c = c+4;
		}
	}

	string legend_fn = "SCORE." + std::to_string(this->myensHdl->simid) + "CA." + std::to_string(this->jobid) + ".DGRHOGRYScalebar.png";
	lodepng::encode( legend_fn.c_str(), img_rgba , imgx, imgy );

	delete [] img_rgba;
	img_rgba = NULL;
}

#endif