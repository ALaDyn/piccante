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
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"
#else
#include <gsl/gsl_rng.h> 
#include <gsl/gsl_randist.h>
#endif
#include <cstdarg>
#include <vector>

using namespace std;

#define DIMENSIONALITY 3
#include "access.h"
#include "commons.h"
#include "grid.h"
#include "structures.h"
#include "current.h"
#include "em_field.h"
#include "particle_species.h"
//#include "diag_manager.h"
#include "output_manager.h"

#define NPROC_ALONG_Y 16
#define NPROC_ALONG_Z 8

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
    int istep;
    gsl_rng* rng = gsl_rng_alloc(gsl_rng_ranlxd1);

   //*******************************************BEGIN GRID DEFINITION*******************************************************

    grid.setXrange(-2.0,2.0);
    grid.setYrange(-2.0,2.0);
    grid.setZrange(-0.5 ,+0.5);

    grid.setNCells(400,400,100);
    grid.setNProcsAlongY(NPROC_ALONG_Y);
    grid.setNProcsAlongZ(NPROC_ALONG_Z);

    //grid.enableStretchedGrid();
    grid.setXandNxLeftStretchedGrid(-20.0,250);
    grid.setXandNxRightStretchedGrid(20.0,250);
    grid.setYandNyLeftStretchedGrid(-20.0,250);
    grid.setYandNyRightStretchedGrid(20.0,250);

    grid.setBoundaries(xPBC|yPBC|zPBC);
    grid.mpi_grid_initialize(&narg, args);
    grid.setCourantFactor(0.98);

    grid.setSimulationTime(100.05);

    grid.with_particles = YES;//NO;
    grid.with_current   = YES;//YES;
    //double start, beta_mw;	int frequency_of_shifts;
    //grid.setMovingWindow(start=0, beta_mw=0.0, frequency_of_shifts=10);

    grid.setMasterProc(0);

    grid.finalize();

    srand (time(NULL));
    grid.initRNG(rng,RANDOM_NUMBER_GENERATOR_SEED);

    grid.visualDiag();

 	//********************************************END GRID DEFINITION********************************************************

	//*******************************************BEGIN FIELD DEFINITION*********************************************************
    myfield.allocate(&grid);
    myfield.setAllValuesToZero();    

    laserPulse pulse1;
    pulse1.type = GAUSSIAN;                   
    pulse1.polarization = CIRCULAR_POLARIZATION;
    pulse1.t_FWHM = 5.0;
    pulse1.waist = 3.0;
    pulse1.focus_position = 0.0;
    pulse1.laser_pulse_initial_position = -15;
    pulse1.lambda0 = 1.0;
    pulse1.normalized_amplitude = 10.0;
    pulse1.rotation = false;
    pulse1.angle = 2.0*M_PI*(45.0/360.0);
    pulse1.rotation_center_along_x = 0.0;

    laserPulse pulse2;
    pulse2 = pulse1;
    pulse2.angle = 2.0*M_PI*(15.0/360.0);
    pulse2.laser_pulse_initial_position = -10.0;

    //myfield.addPulse(&pulse1);
    //myfield.addPulse(&pulse2);

    myfield.boundary_conditions();
    //myfield.smooth_filter(10);   

    current.allocate(&grid);
    current.setAllValuesToZero();

 	//*******************************************END FIELD DEFINITION***********************************************************

	//*******************************************BEGIN SPECIES DEFINITION*********************************************************
    PLASMA plasma1;
    plasma1.density_function = box;      
    plasma1.setMinBox(-10.0 , -10.0, grid.rmin[2]);  
    plasma1.setMaxBox(10.0, 10.0, grid.rmax[2]);    
    plasma1.setRampLength(0.2);                     
    plasma1.setDensityCoefficient(1.0);         
    plasma1.setRampMinDensity(0.001);           

    SPECIE  electrons1(&grid);
    electrons1.plasma = plasma1;
    electrons1.setParticlesPerCellXYZ(3, 3, 3);       
    electrons1.setName("ELE1");
    electrons1.type=ELECTRON;
    electrons1.creation();                            
    species.push_back(&electrons1);


    SPECIE electrons2(&grid);
    electrons2.plasma=plasma1;
    electrons2.setParticlesPerCellXYZ(3, 3, 3);
    electrons2.setName("ELE2");
    electrons2.type=ELECTRON;
    electrons2.creation();                            
    species.push_back(&electrons2);



    tempDistrib distribution;
    distribution.setMaxwell(1.0e-5);

    electrons1.add_momenta(rng,0.0, 0.0, -1.0, distribution);
    electrons2.add_momenta(rng,0.0, 0.0, 1.0, distribution);
  
   	//*******************************************END SPECIES DEFINITION***********************************************************

	//*******************************************BEGIN DIAG DEFINITION**************************************************

   OUTPUT_MANAGER manager(&grid, &myfield, &current, species);

    manager.addEFieldFrom(5.0,5.0);
		manager.addBFieldFrom(5.0,5.0);

    manager.addSpecDensityBinaryFrom("ELE1", 5.0,5.0);
    manager.addSpecDensityBinaryFrom("ELE2", 5.0,5.0);
   
    manager.addCurrentBinaryFrom(5.0,5.0);


    manager.addSpecPhaseSpaceBinaryFrom("ELE1", 5.0,5.0);
    manager.addSpecPhaseSpaceBinaryFrom("ELE2", 5.0,5.0);
   

    manager.addDiagFrom(0.0,2.0);
	
    manager.initialize(DIRECTORY_OUTPUT);

    manager.autoVisualDiag();



	//*******************************************END DIAG DEFINITION**************************************************

	//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ MAIN CYCLE (DO NOT MODIFY) @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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
