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
#include <ctime>
#if defined(_MSC_VER)
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"
#else
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#endif
#include <cstdarg>
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
#include "utilities.h"
#define NPROC_ALONG_Y 2
#define NPROC_ALONG_Z 1

#define _RESTART_FROM_DUMP 1
#define _DO_RESTART false
#define DO_DUMP true
#define TIME_BTW_DUMP 10

#define DIRECTORY_OUTPUT "TEST"
#define DIRECTORY_DUMP "DUMP"
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

    //*******************************************BEGIN GRID DEFINITION*******************************************************

    grid.setXrange(-29.0, 1.0);
    grid.setYrange(-16, 16);
    grid.setZrange(-16, +16);

    grid.setNCells(30*24, 80+48, 80+48);
	grid.setNProcsAlongY(NPROC_ALONG_Y);
	grid.setNProcsAlongZ(NPROC_ALONG_Z);

    //grid.enableStretchedGrid();
    //grid.setXandNxLeftStretchedGrid(-20.0,1000);
    grid.setYandNyLeftStretchedGrid(-8.0,21);
    //grid.setXandNxRightStretchedGrid(20.0,1000);
    grid.setYandNyRightStretchedGrid(8.0,21);

    grid.setBoundaries(xOpen | yOpen | zPBC);
	grid.mpi_grid_initialize(&narg, args);
	grid.setCourantFactor(0.98);

    grid.setSimulationTime(50.0);

    grid.with_particles = YES;
    grid.with_current = YES;

    grid.setStartMovingWindow(0);
    grid.setBetaMovingWindow(1.0);
    grid.setFrequencyMovingWindow(20);

	grid.setMasterProc(0);	

	srand(time(NULL));
	grid.initRNG(rng, RANDOM_NUMBER_GENERATOR_SEED);

    grid.finalize();

	grid.visualDiag();

    //********************************************END GRID DEFINITION********************************************************

    //*******************************************BEGIN FIELD DEFINITION*********************************************************
	myfield.allocate(&grid);
	myfield.setAllValuesToZero();

	laserPulse pulse1;
    pulse1.setGaussianPulse(3.0, 10.0, 3.0);
    pulse1.setPPol();
    pulse1.setPulseInitialPosition(-10.1);
    pulse1. setFocusPosition( 0.0);
    //pulse1.setRotationAngleAndCenter( 2.0*M_PI*(-30.0 / 360.0, 0.0);

    myfield.addPulse(&pulse1);

	myfield.boundary_conditions();

	current.allocate(&grid);
	current.setAllValuesToZero();
    //*******************************************END FIELD DEFINITION***********************************************************

    //*******************************************BEGIN SPECIES DEFINITION*********************************************************
	PLASMA plasma1;
    plasma1.density_function = left_right_linear_ramp;
    plasma1.setXRangeBox(0.0,100.0);
    plasma1.setYRangeBox(grid.rmin[1]*0.95,grid.rmax[1]*0.95);
	plasma1.setZRangeBox(grid.rmin[2],grid.rmax[2]);
    plasma1.setLeftRampLength(5.0);
    plasma1.setRightRampLength(5.0);
    plasma1.setLeftScaleLength(2.0);
    plasma1.setRightScaleLength(2.0);
    plasma1.setDensityCoefficient(0.01);
    plasma1.setLeftRampMinDensity(0.0);
    plasma1.setRightRampMinDensity(0.0);

    SPECIE  electrons1(&grid);
	electrons1.plasma = plasma1;
    electrons1.setParticlesPerCellXYZ(4, 4, 1);
	electrons1.setName("ELE1");
	electrons1.type = ELECTRON;
    electrons1.creation();
    species.push_back(&electrons1);


	SPECIE ions1(&grid);
	ions1.plasma = plasma1;
    ions1.setParticlesPerCellXYZ(4, 4, 1);
	ions1.setName("ION1");
    ions1.type = ION;
    ions1.Z = 6.0;
    ions1.A = 12.0;
    ions1.creation();
    species.push_back(&ions1);


	tempDistrib distribution;
    distribution.setWaterbag(1.0e-10);

    electrons1.add_momenta(rng,0.0,0.0,0.0,distribution);
    ions1.add_momenta(rng,0.0, 0.0, 0.0, distribution);


    for (spec_iterator = species.begin(); spec_iterator != species.end(); spec_iterator++){
        (*spec_iterator)->printParticleNumber();
    }

     //*******************************************END SPECIED DEFINITION***********************************************************

    //*******************************************BEGIN DIAGNOSTICS DEFINITION**************************************************
	OUTPUT_MANAGER manager(&grid, &myfield, &current, species);

    emProbe *probe1=new emProbe;
    probe1->setPointCoordinate(0,0,0);
    probe1->setName("A");

    outDomain *plane1= new outDomain;
    plane1->setPointCoordinate(0,0,0);
    plane1->setFreeDimensions(1 ,0,0);
    plane1->setName("B");

    outDomain *plane2= new outDomain;
    plane2->setPointCoordinate(0,0,0);
    plane2->setFreeDimensions(1 ,0,1);
    plane2->setName("D");
    plane2->setXRange(-20,10);
    plane2->setYRange(-5,10);
    plane2->setZRange(-5,10);
    manager.addEBFieldFrom(0.0, 5.0);
    manager.addEBFieldFrom(plane2, 0.0, 5.0);
    //manager.addEBFieldProbeFrom(probe1,0.0,0.1);

    //manager.addSpeciesDensityFrom(electrons1.name, 0.0, 5.0);
    manager.addSpeciesDensityFrom(plane2,electrons1.name, 0.0, 2.0);
    //manager.addSpeciesDensityFrom(ions1.name, 0.0, 5.0);

    //manager.addCurrentFrom(0.0, 5.0);
    //manager.addCurrentFrom(plane1, 0.0, 5.0);

    //manager.addSpeciesPhaseSpaceFrom(electrons1.name, 10.0, 10.0);
    //manager.addSpeciesPhaseSpaceFrom(ions1.name, 10.0, 10.0);

	manager.addDiagFrom(0.0, 1.0);

	manager.initialize(DIRECTORY_OUTPUT);
    //*******************************************END DIAGNOSTICS DEFINITION**************************************************
    grid.setDumpPath(DIRECTORY_DUMP);
    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ MAIN CYCLE (DO NOT MODIFY!!) @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	if (grid.myid == grid.master_proc){
		printf("----- START temporal cicle -----\n");
		fflush(stdout);
	}

	int Nstep = grid.getTotalNumberOfTimesteps();
    int dumpID=1, dumpEvery;
    if(DO_DUMP){
        dumpEvery= (int)TIME_BTW_DUMP/grid.dt;
    }

    grid.istep=0;
    if (_DO_RESTART){
        dumpID=_RESTART_FROM_DUMP;
        restartFromDump(&dumpID, &grid, &myfield, species);
    }
    while(grid.istep <= Nstep)
    {

        grid.printTStepEvery(FREQUENCY_STDOUT_STATUS);

        manager.callDiags(grid.istep);

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

        moveWindow(&grid, &myfield, species);

        grid.istep++;
        if(DO_DUMP){
            if (grid.istep!=0 && !(grid.istep % (dumpEvery))) {
                dumpFilesForRestart(&dumpID, &grid, &myfield, species);
            }
        }
    }

    manager.close();
    MPI_Finalize();
    exit(1);

}
