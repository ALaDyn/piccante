/*******************************************************************************
This file is part of piccante.

piccante is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

piccante is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with piccante.  If not, see <http://www.gnu.org/licenses/>.
*******************************************************************************/

#define _USE_MATH_DEFINES

#include <mpi.h>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <malloc.h>
#include <cmath>
#include <iomanip>
#include <cstring>
#include <ctime>       /* time */
#if defined(_MSC_VER)
#include "gsl/gsl_rng.h" // gnu scientific linux per generatore di numeri casuali
#include "gsl/gsl_randist.h"
#else
#include <gsl/gsl_rng.h> // gnu scientific linux per generatore di numeri casuali
#include <gsl/gsl_randist.h>
#endif
#include <cstdarg> //Per chiamare funzioni con numero variabile di argomenti
#include <vector>

using namespace std;

#define DIMENSIONALITY 2


#include "access.h"
#include "commons.h"
#include "grid.h"
#include "structures.h"
#include "current.h"
#include "em_field.h"
#include "particle_species.h"
#include "output_manager.h"

#define NPROC_ALONG_Y 1
#define NPROC_ALONG_Z 1


#define DIRECTORY_OUTPUT "TEST"
#define RANDOM_NUMBER_GENERATOR_SEED 5489
#define FREQUENCY_STDOUT_STATUS 5

int main(int narg, char **args)
{
	GRID grid;
	EM_FIELD myfield;
	CURRENT current;
	std::vector<SPECIE*> species;
	vector<SPECIE*>::const_iterator spec_iterator;
	int istep;
	gsl_rng* rng = gsl_rng_alloc(gsl_rng_ranlxd1);

	//*******************************************INIZIO DEFINIZIONE GRIGLIA*******************************************************

	grid.setXrange(-30.0, 30.0);
	grid.setYrange(-20.0, 20.0);
	grid.setZrange(-1.0, +1.0);

	grid.setNCells(2000, 1500, 0);
	grid.setNProcsAlongY(NPROC_ALONG_Y);
	grid.setNProcsAlongZ(NPROC_ALONG_Z);

	grid.enableStretchedGrid();
	grid.setXandNxLeftStretchedGrid(-10.0, 500);
	grid.setYandNyLeftStretchedGrid(-5.0, 70);
	grid.setXandNxRightStretchedGrid(10.0, 500);
	grid.setYandNyRightStretchedGrid(5.0, 70);

	grid.setBoundaries(xOpen | yOpen | zPBC); //LUNGO Z c'è solo PBC al momento !
	grid.mpi_grid_initialize(&narg, args);
	grid.setCourantFactor(0.98);

	grid.setSimulationTime(15.0);

	grid.with_particles = YES;//NO;
	grid.with_current = YES;//YES;

	//grid.setStartMovingWindow(0);
	//grid.setBetaMovingWindow(BETA);
	//grid.setFrequencyMovingWindow(FREQUENCY);


	grid.setMasterProc(0);

	grid.finalize();

	srand((unsigned int)time(NULL));
	grid.initRNG(rng, RANDOM_NUMBER_GENERATOR_SEED);

	grid.visualDiag();

	//********************************************FINE DEFINIZIONE GRIGLIA********************************************************

	//*******************************************INIZIO DEFINIZIONE CAMPI*********************************************************
	myfield.allocate(&grid);
	myfield.setAllValuesToZero();

	laserPulse pulse1;
	pulse1.type = GAUSSIAN;                        //Opzioni : GAUSSIAN, PLANE_WAVE, COS2_PLANE_WAVE
	pulse1.polarization = P_POLARIZATION;
	pulse1.t_FWHM = 9.0;
	pulse1.waist = 3.0;
	pulse1.focus_position = 0.0;
	pulse1.laser_pulse_initial_position = -9.01;
	pulse1.lambda0 = 1.0;
	pulse1.normalized_amplitude = 5.0;
	pulse1.rotation = false;
	pulse1.angle = 2.0*M_PI*(-90.0 / 360.0);
	pulse1.rotation_center_along_x = 0.0;

	laserPulse pulse2;
	pulse2 = pulse1;
	pulse2.angle = 2.0*M_PI*(15.0 / 360.0);
	pulse2.laser_pulse_initial_position = -10.0;

	myfield.addPulse(&pulse1);
	//myfield.addPulse(&pulse2);

	myfield.boundary_conditions();
	//myfield.smooth_filter(10);   //DA NON USARE CON GRIGLIA STRETCHATA !

	current.allocate(&grid);
	current.setAllValuesToZero();
	//*******************************************FINE DEFINIZIONE CAMPI***********************************************************

	//*******************************************INIZIO DEFINIZIONE SPECIE*********************************************************
	PLASMA plasma1;
	plasma1.density_function = left_linear_ramp;      //Opzioni: box, left_linear_ramp, left_soft_ramp, left_grating
	plasma1.setXRangeBox(0.0, 0.5);                  //double (* distrib_function)(double x, double y, double z, PLASMAparams plist, double Z, double A)
	plasma1.setYRangeBox(-50.0, 50.0);                 //PLASMAparams: rminbox[3], rmaxbox[3], ramp_length, density_coefficient,
	plasma1.setZRangeBox(grid.rmin[2], grid.rmax[2]);
	plasma1.setRampLength(0.25);                       //ramp_min_density,void *additional_params
	plasma1.setDensityCoefficient(10.0);         // Per grating double g_depth = paramlist[0];double g_lambda = paramlist[1];
	plasma1.setRampMinDensity(0.001);                 //double g_phase = paramlist[2];
	double gratingParams[3] = { 0.2, sqrt(2), 0.0 };
	plasma1.setAdditionalParams(gratingParams);

	PLASMA plasma2;
	plasma2.density_function = box;
	plasma2.setMinBox(3.0, -20.0, grid.rmin[2]);
	plasma2.setMaxBox(4.0, 20.0, grid.rmax[2]);
	plasma2.setDensityCoefficient(0.5);

	SPECIE  electrons1(&grid);
	electrons1.plasma = plasma1;
	electrons1.setParticlesPerCellXYZ(5, 5, 1);       //Se < 1 il nPPC viene sostituito con 1
	electrons1.setName("ELE1");
	electrons1.type = ELECTRON;
	electrons1.creation();                            //electrons.isTestSpecies=true disabilita deposizione corrente.
	species.push_back(&electrons1);


	SPECIE electrons2(&grid);
	electrons2.plasma = plasma1;
	electrons2.setParticlesPerCellXYZ(100, 1, 1);
	electrons2.setName("ELE2");
	electrons2.type = ELECTRON;
	// electrons2.creation();                            //electrons.isTestSpecies=true disabilita deposizione corrente.
	// species.push_back(&electrons2);


	SPECIE ions1(&grid);
	ions1.plasma = plasma1;
	ions1.setParticlesPerCellXYZ(5, 5, 1);
	ions1.setName("ION1");
	ions1.type = ION;
	ions1.Z = 50.0;
	ions1.A = 100.0;
	ions1.creation();
	species.push_back(&ions1);


	SPECIE  ions2(&grid);
	ions2.plasma = plasma2;
	ions2.setParticlesPerCellXYZ(9, 9, 1);       //Se < 1 il nPPC viene sostituito con 1
	ions2.setName("ION2");
	ions2.type = ION;
	ions2.Z = 1.0;
	ions2.A = 1.0;
	//    ions2.creation();                            //electrons.isTestSpecies=true disabilita deposizione corrente.
	//    species.push_back(&ions2);


	// WATERBAG        : [P0] -P0< px < + P0, -P0< py < + P0, -P0< pz < + P0 uniforme
	// WATERBAG_3TEMP  : [P0_X,P0_Y,P0_Z] -P0_X< px < + P0_X, -P0_Y< py < + P0_Y, -P0_Z< pz < + P0_Z uniforme
	// UNIF_SPHERE     : [P0] -P0 < p < +P0 uniforme
	// SUPERGAUSSIAN   : [P0, ALPHA] f(p) = C*exp(-abs(p/P0)^(ALPHA))
	// MAXWELL         : [Ta] Maxwell alla Macchi
	// JUTTNER         : [a] f(p) = C*exp(-a*gamma(p)) [DA RISCRIVERE]
	//si può ancora usare electrons.add_momenta(1.0, 0.0, 0);

	tempDistrib distribution;
	distribution.setMaxwell(0.5);

	//electrons1.add_momenta(rng, 0.0, 0.0, 0.0, distribution);
	//   electrons2.add_momenta(rng,-1.0, 0.0, 0.0, distribution);
	//    ions1.add_momenta(rng,0.0, 0.0, 0.0, distribution);
	//    ions2.add_momenta(rng,0.0, 0.0, 0.0, distribution);
	//    //*******************************************FINE DEFINIZIONE CAMPI***********************************************************

	//*******************************************INIZIO DEFINIZIONE DIAGNOSTICHE**************************************************

	OUTPUT_MANAGER manager(&grid, &myfield, &current, species);
	//AGGIUNGERE CONTROLLI DIAGS
	//manager.addEMFieldBinaryFrom(0.0, 2.0);

	//manager.addSpecDensityBinaryFrom(electrons1.name, 0.0, 2.0);
	//manager.addSpecDensityBinaryFrom(ions1.name, 0.0, 2.0);

	//manager.addCurrentBinaryFrom(0.0, 2.0);

	//manager.addSpecPhaseSpaceBinaryFrom(electrons2.name, 0.0, 5.0);

	manager.addDiagFrom(0.0, 1.0);

	manager.initialize(DIRECTORY_OUTPUT);

	//manager.autoVisualDiag();



	//*******************************************FINE DEFINIZIONE DIAGNOSTICHE**************************************************

	//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ CICLO PRINCIPALE (NON MODIFICARE) @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	if (grid.myid == grid.master_proc){
		printf("----- START temporal cicle -----\n");
		fflush(stdout);
	}


	int Nstep = grid.getTotalNumberOfTimesteps();
	for (istep = 0; istep <= Nstep; istep++)
	{
		grid.istep = istep;

		grid.printTStepEvery(FREQUENCY_STDOUT_STATUS);


		manager.callDiags(istep);  /// deve tornare all'inizo del ciclo

		myfield.openBoundariesE_1();
		myfield.new_halfadvance_B();
		myfield.boundary_conditions();



		/*   for(spec_iterator=species.begin(); spec_iterator!=species.end(); spec_iterator++){
				(*spec_iterator)->density_deposition_standard(&current);
				}*/

		current.setAllValuesToZero();

		for (spec_iterator = species.begin(); spec_iterator != species.end(); spec_iterator++){
			(*spec_iterator)->current_deposition_standard(&current);
		}

		current.pbc();

		for (spec_iterator = species.begin(); spec_iterator != species.end(); spec_iterator++){
			(*spec_iterator)->position_parallel_pbc();
		}

		//myfield.boundary_conditions();

		myfield.openBoundariesB();
		myfield.new_advance_E(&current);
		//myfield.new_advance_E();

		myfield.boundary_conditions();

		myfield.openBoundariesE_2();
		myfield.new_halfadvance_B();

		myfield.boundary_conditions();

		for (spec_iterator = species.begin(); spec_iterator != species.end(); spec_iterator++){
			(*spec_iterator)->momenta_advance(&myfield);
		}

		grid.time += grid.dt;


		grid.move_window();

		myfield.move_window();
		for (spec_iterator = species.begin(); spec_iterator != species.end(); spec_iterator++){
			(*spec_iterator)->move_window();
		}
	}

	manager.close();
	MPI_Finalize();
	exit(1);

}
