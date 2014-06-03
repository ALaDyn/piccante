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

#define NPROC_ALONG_Y 2
#define NPROC_ALONG_Z 1

#define _RESTART_FROM_DUMP 1
#define _DO_RESTART false
#define DO_DUMP true
#define TIME_BTW_DUMP 50

#define DIRECTORY_OUTPUT "TEST"
#define DIRECTORY_DUMP "TEST"
#define RANDOM_NUMBER_GENERATOR_SEED 5489
#define FREQUENCY_STDOUT_STATUS 5

int main(int narg, char **args)
{
	GRID grid;
	EM_FIELD myfield;
	CURRENT current;
	std::vector<SPECIE*> species;
	vector<SPECIE*>::const_iterator spec_iterator;
    gsl_rng* rng = gsl_rng_alloc(gsl_rng_ranlxd1);

	//*******************************************INIZIO DEFINIZIONE GRIGLIA*******************************************************

    grid.setXrange(-29.0, 1.0);
    grid.setYrange(-16, 16);
	grid.setZrange(-1, +1);

    grid.setNCells(30*24, 80+48, 1);
	grid.setNProcsAlongY(NPROC_ALONG_Y);
	grid.setNProcsAlongZ(NPROC_ALONG_Z);

    grid.enableStretchedGrid();
    //grid.setXandNxLeftStretchedGrid(-20.0,1000);
    grid.setYandNyLeftStretchedGrid(-8.0,21);
    //grid.setXandNxRightStretchedGrid(20.0,1000);
    grid.setYandNyRightStretchedGrid(8.0,21);

    grid.setBoundaries(xOpen | yOpen | zPBC); //LUNGO Z c'è solo PBC al momento !
	grid.mpi_grid_initialize(&narg, args);
	grid.setCourantFactor(0.98);

    grid.setSimulationTime(50.0);

    grid.with_particles = YES;//NO;
    grid.with_current = YES;//YES;

    grid.setStartMovingWindow(0);
    grid.setBetaMovingWindow(1.0);
    grid.setFrequencyMovingWindow(20);

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
    pulse1.type = GAUSSIAN;                        //Opzioni : GAUSSIAN, PLANE_WAVE, COS2_PLANE_WAVE
    pulse1.polarization = P_POLARIZATION;
    pulse1.t_FWHM = 10.0;
    pulse1.laser_pulse_initial_position = -10.1;
    pulse1.lambda0 = 1.0;
    pulse1.normalized_amplitude = 3.0;
    pulse1.waist = 4.0;
    pulse1.focus_position = 0.0;
    pulse1.rotation = false;
    pulse1.angle = 2.0*M_PI*(-30.0 / 360.0);
    pulse1.rotation_center_along_x = 0.0;

    myfield.addPulse(&pulse1);

    laserPulse pulse2;
    pulse2 = pulse1;
    pulse2.angle = 2.0*M_PI*(30.0 / 360.0);

    //myfield.addPulse(&pulse2);

	myfield.boundary_conditions();

	current.allocate(&grid);
	current.setAllValuesToZero();
	//*******************************************FINE DEFINIZIONE CAMPI***********************************************************

	//*******************************************INIZIO DEFINIZIONE SPECIE*********************************************************
	PLASMA plasma1;
    plasma1.density_function = left_right_linear_ramp;      //Opzioni: box, left_linear_ramp, left_soft_ramp, left_grating
    plasma1.setXRangeBox(0.0,100.0);                  //double (* distrib_function)(double x, double y, double z, PLASMAparams plist, int Z, int A)
    plasma1.setYRangeBox(grid.rmin[1]*0.95,grid.rmax[1]*0.95);                 //PLASMAparams: rminbox[3], rmaxbox[3], ramp_length, density_coefficient,
	plasma1.setZRangeBox(grid.rmin[2],grid.rmax[2]);
    plasma1.setLeftRampLength(5.0);                       //ramp_min_density,void *additional_params
    plasma1.setRightRampLength(5.0);                       //ramp_min_density,void *additional_params
    plasma1.setLeftScaleLength(2.0);                       //ramp_min_density,void *additional_params
    plasma1.setRightScaleLength(2.0);                       //ramp_min_density,void *additional_params
    plasma1.setDensityCoefficient(0.01);         // Per grating double g_depth = paramlist[0];double g_lambda = paramlist[1];
    plasma1.setLeftRampMinDensity(0.0);                 //double g_phase = paramlist[2];
    plasma1.setRightRampMinDensity(0.0);                 //double g_phase = paramlist[2];

    PLASMA plasma2;
    plasma2.density_function = box;      //Opzioni: box, left_linear_ramp, left_soft_ramp, left_grating
    plasma2.setXRangeBox(0.8,0.85);                  //double (* distrib_function)(double x, double y, double z, PLASMAparams plist, int Z, int A)
    plasma2.setYRangeBox(grid.rmin[1],grid.rmax[1]);                 //PLASMAparams: rminbox[3], rmaxbox[3], ramp_length, density_coefficient,
    plasma2.setZRangeBox(grid.rmin[2],grid.rmax[2]);
    plasma2.setRampLength(0.5);                       //ramp_min_density,void *additional_params
    plasma2.setDensityCoefficient(10);         // Per grating double g_depth = paramlist[0];double g_lambda = paramlist[1];
    plasma2.setRampMinDensity(0.0);                 //double g_phase = paramlist[2];


	SPECIE  electrons1(&grid);
	electrons1.plasma = plasma1;
    electrons1.setParticlesPerCellXYZ(4, 4, 1);       //Se < 1 il nPPC viene sostituito con 1
	electrons1.setName("ELE1");
	electrons1.type = ELECTRON;
    electrons1.creation();                            //electrons.isTestSpecies=true disabilita deposizione corrente.
    species.push_back(&electrons1);


	SPECIE ions1(&grid);
	ions1.plasma = plasma1;
    ions1.setParticlesPerCellXYZ(4, 4, 1);
	ions1.setName("ION1");
    ions1.type = ION;
    ions1.Z = 6.0;
    ions1.A = 12.0;
    //ions1.creation();
    //species.push_back(&ions1);

    //   SPECIE  electrons2(&grid);
    //electrons2.plasma = plasma1;
    //electrons2.setParticlesPerCellXYZ(1, 4, 4);       //Se < 1 il nPPC viene sostituito con 1
    //electrons2.setName("ELE2");
    //electrons2.type = ELECTRON;
    //electrons2.creation();                            //electrons.isTestSpecies=true disabilita deposizione corrente.
    //species.push_back(&electrons2);


//    SPECIE ions2(&grid);
//    ions2.plasma = plasma2;
//    ions2.setParticlesPerCellXYZ(10, 10, 1);
//    ions2.setName("ION2");
//    ions2.type = ION;
//    ions2.Z = 1.0;
//    ions2.A = 1.0;
//    //ions2.creation();
//    //species.push_back(&ions2);

	tempDistrib distribution;
    distribution.setWaterbag(1.0e-10);

    //electrons1.add_momenta(rng,0.0,0.0,0.0,distribution);
    //ions1.add_momenta(rng,0.0, 0.0, 0.0, distribution);
    //    electrons2.add_momenta(rng,0.0,-5.0,0.0,distribution);
//    ions2.add_momenta(rng,0.0, 0.0, 0.0, distribution);


    for (spec_iterator = species.begin(); spec_iterator != species.end(); spec_iterator++){
        (*spec_iterator)->printParticleNumber();
    }

	//    //*******************************************FINE DEFINIZIONE CAMPI***********************************************************

	//*******************************************INIZIO DEFINIZIONE DIAGNOSTICHE**************************************************

	OUTPUT_MANAGER manager(&grid, &myfield, &current, species);

    emProbe *probe1=new emProbe;
    probe1->setPointCoordinate(30,0,0);
    probe1->setName("A");

    outDomain *plane1= new outDomain;
    plane1->setPointCoordinate(0,0,0);
    plane1->setFreeDimensions(1 ,0,0);
    plane1->setName("B");

    outDomain *plane2= new outDomain;
    plane2->setPointCoordinate(0,0,0);
    plane2->setFreeDimensions(1 ,0,0);
    plane2->setName("C");

    manager.addEFieldFrom(plane1, 0.0, 5.0);
    manager.addEFieldFrom(plane2, 0.0, 2.0);
    //manager.addEMFieldProbeFrom(probe1,0.0,0.1);
    //manager.addEMFieldPlaneFrom(plane1,0.0,1.0);

    manager.addSpecDensityBinaryFrom(electrons1.name, 0.0, 5.0);
    manager.addSpecDensityBinaryFrom(plane1,electrons1.name, 0.0, 2.0);
    //manager.addSpecDensityBinaryFrom(ions1.name, 0.0, 5.0);
    //manager.addSpecDensityBinaryFrom(electrons2.name, 0.0, 2.0);
    //manager.addSpecDensityBinaryFrom(ions2.name, 0.0, 2.0);

    //manager.addCurrentBinaryFrom(0.0, 5.0);

    //manager.addSpecPhaseSpaceBinaryFrom(electrons1.name, 10.0, 10.0);
    //manager.addSpecPhaseSpaceBinaryFrom(electrons2.name, 0.0, 2.0);
    //manager.addSpecPhaseSpaceBinaryFrom(ions1.name, 10.0, 10.0);
    //manager.addSpecPhaseSpaceBinaryFrom(ions2.name, 0.0, 5.0);

	manager.addDiagFrom(0.0, 1.0);

	manager.initialize(DIRECTORY_OUTPUT);

	//*******************************************FINE DEFINIZIONE DIAGNOSTICHE**************************************************

	//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ CICLO PRINCIPALE (NON MODIFICARE) @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	if (grid.myid == grid.master_proc){
		printf("----- START temporal cicle -----\n");
		fflush(stdout);
	}

	int Nstep = grid.getTotalNumberOfTimesteps();
    int dumpID=1, dumpEvery=40;
    grid.istep=0;
    if(DO_DUMP){
        dumpEvery= (int)TIME_BTW_DUMP/grid.dt;
    }
    if (_DO_RESTART){
        dumpID=_RESTART_FROM_DUMP;
        std::ifstream dumpFile;
        std::stringstream dumpName;
        dumpName << DIRECTORY_DUMP << "/DUMP_";
        dumpName<< std::setw(2)<< std::setfill('0') << std::fixed << dumpID << "_";
        dumpName<< std::setw(5)<< std::setfill('0') << std::fixed << grid.myid << ".bin";
        dumpFile.open(dumpName.str().c_str());

        grid.reloadDump(dumpFile);
        myfield.reloadDump(dumpFile);
        for (spec_iterator = species.begin(); spec_iterator != species.end(); spec_iterator++){
            (*spec_iterator)->reloadDump(dumpFile);
        }
        dumpFile.close();
        dumpID++;
        grid.istep++;
    }
    for (; grid.istep <= Nstep; grid.istep++)
	{

		grid.printTStepEvery(FREQUENCY_STDOUT_STATUS);


        manager.callDiags(grid.istep);  /// deve tornare all'inizo del ciclo

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
         if(DO_DUMP){
            if (grid.istep!=0 && !(grid.istep % (dumpEvery))) {
                std::ofstream dumpFile;
                std::stringstream dumpName;
                dumpName << DIRECTORY_OUTPUT << "/DUMP_";
                dumpName<< std::setw(2)<< std::setfill('0') << std::fixed << dumpID << "_";
                dumpName<< std::setw(5)<< std::setfill('0') << std::fixed << grid.myid << ".bin";
                dumpFile.open(dumpName.str().c_str());

                grid.dump(dumpFile);
                myfield.dump(dumpFile);
                for (spec_iterator = species.begin(); spec_iterator != species.end(); spec_iterator++){
                    (*spec_iterator)->dump(dumpFile);
                }
                dumpFile.close();
                dumpID++;
            }
        }
	}

	manager.close();
	MPI_Finalize();
	exit(1);

}
