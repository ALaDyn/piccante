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

#define DIMENSIONALITY 1

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

    grid.setXrange(-50.0, 0.0);
    grid.setYrange(-1.0, 1.0);
	grid.setZrange(-1.0, +1.0);

    grid.setNCells(2500, 0, 0);
	grid.setNProcsAlongY(NPROC_ALONG_Y);
	grid.setNProcsAlongZ(NPROC_ALONG_Z);

    //grid.enableStretchedGrid();
	//grid.setXandNxLeftStretchedGrid(-15.0,1000);
    //grid.setYandNyLeftStretchedGrid(-5.0, 70);
	//grid.setXandNxRightStretchedGrid(15.0,1000);
    //grid.setYandNyRightStretchedGrid(5.0, 70);

    grid.setBoundaries(xOpen | yOpen | zPBC); //LUNGO Z c'è solo PBC al momento !
	grid.mpi_grid_initialize(&narg, args);
	grid.setCourantFactor(0.98);

    grid.setSimulationTime(60.0);

    grid.with_particles = YES;//NO;
    grid.with_current = YES;//YES;

    grid.setStartMovingWindow(0);
    grid.setBetaMovingWindow(1.0);
    //grid.setFrequencyMovingWindow(FREQUENCY);

	grid.setMasterProc(0);	

	srand(time(NULL));
	grid.initRNG(rng, RANDOM_NUMBER_GENERATOR_SEED);

    grid.finalize();

	grid.visualDiag();

	//********************************************FINE DEFINIZIONE GRIGLIA********************************************************

	//*******************************************INIZIO DEFINIZIONE CAMPI*********************************************************
	myfield.allocate(&grid);
	myfield.setAllValuesToZero();

	laserPulse pulse1;
    pulse1.type = COS2_PLANE_WAVE;                        //Opzioni : GAUSSIAN, PLANE_WAVE, COS2_PLANE_WAVE
    pulse1.polarization = P_POLARIZATION;
    pulse1.t_FWHM = 5.0;
    pulse1.laser_pulse_initial_position = -6.0;
    pulse1.lambda0 = 1.0;
    pulse1.normalized_amplitude = 1.0;

    //pulse1.waist = 3.0;
    //pulse1.focus_position = 0.0;
    //pulse1.rotation = false;
    //pulse1.angle = 2.0*M_PI*(-90.0 / 360.0);
    //pulse1.rotation_center_along_x = 0.0;

    myfield.addPulse(&pulse1);

	myfield.boundary_conditions();

	current.allocate(&grid);
	current.setAllValuesToZero();
	//*******************************************FINE DEFINIZIONE CAMPI***********************************************************

	//*******************************************INIZIO DEFINIZIONE SPECIE*********************************************************
	PLASMA plasma1;
    plasma1.density_function = left_soft_ramp;      //Opzioni: box, left_linear_ramp, left_soft_ramp, left_grating
    plasma1.setXRangeBox(0.0,100.0);                  //double (* distrib_function)(double x, double y, double z, PLASMAparams plist, int Z, int A)
    plasma1.setYRangeBox(grid.rmin[1],grid.rmax[1]);                 //PLASMAparams: rminbox[3], rmaxbox[3], ramp_length, density_coefficient,
	plasma1.setZRangeBox(grid.rmin[2],grid.rmax[2]);
    plasma1.setRampLength(20.0);                       //ramp_min_density,void *additional_params
    plasma1.setDensityCoefficient(0.01);         // Per grating double g_depth = paramlist[0];double g_lambda = paramlist[1];
    plasma1.setRampMinDensity(0.0);                 //double g_phase = paramlist[2];


	SPECIE  electrons1(&grid);
	electrons1.plasma = plasma1;
    electrons1.setParticlesPerCellXYZ(100, 1, 1);       //Se < 1 il nPPC viene sostituito con 1
	electrons1.setName("ELE1");
	electrons1.type = ELECTRON;
	electrons1.creation();                            //electrons.isTestSpecies=true disabilita deposizione corrente.
	species.push_back(&electrons1);


	SPECIE ions1(&grid);
	ions1.plasma = plasma1;
    ions1.setParticlesPerCellXYZ(100, 1, 1);
	ions1.setName("ION1");
    ions1.type = ION;
    ions1.Z = 6.0;
    ions1.A = 12.0;
    ions1.creation();
    species.push_back(&ions1);

	tempDistrib distribution;
    distribution.setMaxwell(1.0);

    electrons1.add_momenta(rng, 0.0, 0.0, 0.0, distribution);
    ions1.add_momenta(rng,0.0, 0.0, 0.0, distribution);

	//    //*******************************************FINE DEFINIZIONE CAMPI***********************************************************

	//*******************************************INIZIO DEFINIZIONE DIAGNOSTICHE**************************************************

	OUTPUT_MANAGER manager(&grid, &myfield, &current, species);

    manager.addEMFieldBinaryFrom(0.0, 2.0);

    manager.addSpecDensityBinaryFrom(electrons1.name, 0.0, 2.0);
    manager.addSpecDensityBinaryFrom(ions1.name, 0.0, 2.0);

    manager.addCurrentBinaryFrom(0.0, 5.0);

    manager.addSpecPhaseSpaceBinaryFrom(electrons1.name, 0.0, 5.0);

	manager.addDiagFrom(0.0, 1.0);

	manager.initialize(DIRECTORY_OUTPUT);

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

		current.setAllValuesToZero();

		for (spec_iterator = species.begin(); spec_iterator != species.end(); spec_iterator++){
			(*spec_iterator)->current_deposition_standard(&current);
		}

		current.pbc();

		for (spec_iterator = species.begin(); spec_iterator != species.end(); spec_iterator++){
			(*spec_iterator)->position_parallel_pbc();
		}		

		myfield.openBoundariesB();
		myfield.new_advance_E(&current);

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
