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


#define DIMENSIONALITY 1

#include "access.h"
#include "commons.h"
#include "grid.h"
#include "structures.h"
#include "current.h"
#include "em_field.h"
#include "particle_species.h"
#include "output_manager.h"
#include "utilities.h"

#define NPROC_ALONG_Y 1
#define NPROC_ALONG_Z 1

#define _RESTART_FROM_DUMP 1
#define _DO_RESTART false
#define DO_DUMP false
#define TIME_BTW_DUMP 10

#define DIRECTORY_OUTPUT "TEST"
#define DIRECTORY_DUMP "DUMP"
#define RANDOM_NUMBER_GENERATOR_SEED 5489
#define FREQUENCY_STDOUT_STATUS 5

#define _FACT 0.333333

#define RADIATION_FRICTION

//#define ESIRKEPOV

int main(int narg, char **args)
{
  GRID grid;
  EM_FIELD myfield;
  CURRENT current;
  std::vector<SPECIE*> species;
  std::vector<SPECIE*>::const_iterator spec_iterator;
  gsl_rng* rng = gsl_rng_alloc(gsl_rng_ranlxd1);

  //*******************************************BEGIN GRID DEFINITION*******************************************************

  grid.setXrange(-20, 20);
  grid.setYrange(-20, 20);
  grid.setZrange(-0.5, 0.5);

  grid.setNCells(1000, 1000, 1);
  grid.setNProcsAlongY(NPROC_ALONG_Y);
  grid.setNProcsAlongZ(NPROC_ALONG_Z);

  //grid.enableStretchedGrid();
  //grid.setXandNxLeftStretchedGrid(-20.0,1000);
  grid.setYandNyLeftStretchedGrid(-8.0, 21);
  //grid.setXandNxRightStretchedGrid(20.0,1000);
  grid.setYandNyRightStretchedGrid(8.0, 21);

  grid.setBoundaries( xOpen | yPBC | zPBC); //LUNGO Z c'è solo PBC al momento !
  grid.mpi_grid_initialize(&narg, args);
  grid.setCourantFactor(0.92);

  grid.setSimulationTime(2.0);

  grid.with_particles = YES;//NO;
  grid.with_current = YES;//YES;

  grid.setStartMovingWindow(0);
  //grid.setBetaMovingWindow(1.0);
  //grid.setFrequencyMovingWindow(20);

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
  pulse1.setCos2PlaneWave();
  pulse1.setWaist(4.0);
  pulse1.setDurationFWHM(5.0);
  pulse1.setNormalizedAmplitude(10.0);
  pulse1.setPPolarization();
  pulse1.setPulseInitialPosition(-5.0);
  pulse1.setFocusPosition(0.0);
  pulse1.setLambda(1.0);
  pulse1.setFocusPosition(0.0);
  // pulse1.setRotationAngleAndCenter(2.0*M_PI*(90.0 / 360.0), 0.0);

  myfield.addPulse(&pulse1);

  laserPulse pulse2;
  pulse2 = pulse1;
  pulse2.angle = 2.0*M_PI*(30.0 / 360.0);

  //myfield.addPulse(&pulse2);

  myfield.boundary_conditions();

  current.allocate(&grid);
  current.setAllValuesToZero();
  //*******************************************END FIELD DEFINITION***********************************************************

  //*******************************************BEGIN SPECIES DEFINITION*********************************************************
  PLASMA plasma1;
  plasma1.density_function = box;
  plasma1.setXRangeBox(0.0, 10.0);
  plasma1.setYRangeBox(-0.5, 0.5);
  plasma1.setZRangeBox(-0.5, 0.5);
  plasma1.setDensityCoefficient(20.0);

  SPECIE  electrons1(&grid);
  electrons1.plasma = plasma1;
  electrons1.setParticlesPerCellXYZ(1, 2, 3);
  electrons1.setName("ELE1");
  electrons1.type = ELECTRON;
  electrons1.creation();
  species.push_back(&electrons1);

  SPECIE  ions1(&grid);
  ions1.plasma = plasma1;
  ions1.setParticlesPerCellXYZ(1, 2, 3);
  ions1.setName("POS2");
  ions1.type = POSITRON;
  //ions1.creation();
  //species.push_back(&ions1);


  tempDistrib distribution;
  distribution.setWaterbag(1.0e-4);

  electrons1.add_momenta(rng, 0.0, 0.0, 0.0, distribution);
  ions1.add_momenta(rng, 0.0, 0.0, 0.0, distribution);

  for (spec_iterator = species.begin(); spec_iterator != species.end(); spec_iterator++){
    (*spec_iterator)->printParticleNumber();
  }

  //*******************************************END SPECIED DEFINITION***********************************************************

  //*******************************************BEGIN DIAGNOSTICS DEFINITION**************************************************
  OUTPUT_MANAGER manager(&grid, &myfield, &current, species);


  outDomain *domain1 = new outDomain;
  domain1->setPointCoordinate(0, 0, 0);
  domain1->setFreeDimensions(1, 1, 1);
  domain1->setName("SUBD");
  domain1->setXRange(-10, 6);
  domain1->setYRange(-5, 5);

  manager.addEBFieldFrom(0.0, 5.0);
  manager.addSpeciesDensityFrom(electrons1.name, 0.0, 5.0);
  manager.addSpeciesDensityFrom(ions1.name, 0.0, 5.0);
  manager.addCurrentFrom(0.0, 5.0);
  manager.addDiagFrom(0.0, 0.5);

  manager.initialize(DIRECTORY_OUTPUT);
  //*******************************************END DIAGNOSTICS DEFINITION**************************************************
  grid.setDumpPath(DIRECTORY_DUMP);

  //RADIATION FRICTION
  double lambda = 2.0*M_PI*9.67e-6;
  if (grid.myid == grid.master_proc){
#ifdef RADIATION_FRICTION
    std::cout << "RADIATION FRICTION: ON   " << "lambda: " << lambda << std::endl;
#else
    std::cout << "RADIATION FRICTION: OFF  " << "lambda: " << "---" << std::endl;
#endif
  }
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ MAIN CYCLE (DO NOT MODIFY!!) @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  if (grid.myid == grid.master_proc){
    printf("----- START temporal cicle -----\n");
    fflush(stdout);
  }

  int Nstep = grid.getTotalNumberOfTimesteps();
  int dumpID = 1, dumpEvery;
  if (DO_DUMP){
    dumpEvery = (int)(TIME_BTW_DUMP / grid.dt);
  }
  grid.istep = 0;
  if (_DO_RESTART){
    dumpID = _RESTART_FROM_DUMP;
    restartFromDump(&dumpID, &grid, &myfield, species);
  }
  while (grid.istep <= Nstep)
  {

    grid.printTStepEvery(FREQUENCY_STDOUT_STATUS);

    manager.callDiags(grid.istep);  /// deve tornare all'inizo del ciclo

    myfield.openBoundariesE_1();
    myfield.new_halfadvance_B();
    myfield.boundary_conditions();

    current.setAllValuesToZero();
    for (spec_iterator = species.begin(); spec_iterator != species.end(); spec_iterator++){
#ifdef ESIRKEPOV
      (*spec_iterator)->current_deposition(&current);
#else
      (*spec_iterator)->current_deposition_standard(&current);
#endif

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
#ifdef RADIATION_FRICTION            
      (*spec_iterator)->momenta_advance_with_friction(&myfield, lambda);
#else
      (*spec_iterator)->momenta_advance(&myfield);
#endif
    }

    //        if(grid.istep%FIELD_FILTER_FREQ==0){
    //            myfield.applyFilter(fltr_Ex, dir_x);
    //            myfield.boundary_conditions();
    //        }

    grid.time += grid.dt;

    moveWindow(&grid, &myfield, species);

    grid.istep++;
    if (DO_DUMP){
      if (grid.istep != 0 && !(grid.istep % (dumpEvery))) {
        dumpFilesForRestart(&dumpID, &grid, &myfield, species);
      }
    }
  }

  manager.close();
  MPI_Finalize();
  exit(0);

}
