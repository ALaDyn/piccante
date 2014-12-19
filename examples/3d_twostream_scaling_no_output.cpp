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


#define DIMENSIONALITY 3

#include "access.h"
#include "commons.h"
#include "grid.h"
#include "structures.h"
#include "current.h"
#include "em_field.h"
#include "particle_species.h"
#include "output_manager.h"
#include "utilities.h"

TODO

#define _RESTART_FROM_DUMP 1
#define _DO_RESTART false
#define DO_DUMP false
#define TIME_BTW_DUMP 0.1

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
  std::vector<SPECIE*>::const_iterator spec_iterator;
  int istep;
  gsl_rng* rng = gsl_rng_alloc(gsl_rng_ranlxd1);

  //*******************************************BEGIN GRID DEFINITION*******************************************************

  grid.setXrange(-2.0*Xfactor, +2.0*Xfactor);
  grid.setYrange(-2.0*Yfactor, +2.0*Yfactor);
  grid.setZrange(-0.5*Zfactor, +0.5*Zfactor);

  int Nxcell=(int)(Xfactor*512);
  int Nycell=(int)(Yfactor*512);
  int Nzcell=(int)(Zfactor*128);
  grid.setNCells(Nxcell, Nycell, Nzcell);
  grid.setNProcsAlongY(NPROC_ALONG_Y);
  grid.setNProcsAlongZ(NPROC_ALONG_Z);

  //grid.enableStretchedGrid();
  grid.setXandNxLeftStretchedGrid(-20.0, 250);
  grid.setXandNxRightStretchedGrid(20.0, 250);
  grid.setYandNyLeftStretchedGrid(-20.0, 250);
  grid.setYandNyRightStretchedGrid(20.0, 250);

  grid.setBoundaries(xPBC | yPBC | zPBC);
  grid.mpi_grid_initialize(&narg, args);
  grid.setCourantFactor(0.98);

  grid.setSimulationTime(0.15);

  grid.with_particles = YES;//NO;
  grid.with_current = YES;//YES;
  //double start, beta_mw;	int frequency_of_shifts;
  //grid.setMovingWindow(start=0, beta_mw=0.0, frequency_of_shifts=10);

  grid.setMasterProc(0);

  grid.finalize();

  srand(time(NULL));
  grid.initRNG(rng, RANDOM_NUMBER_GENERATOR_SEED);

  grid.visualDiag();

  //********************************************END GRID DEFINITION********************************************************

  //*******************************************BEGIN FIELD DEFINITION*********************************************************
  myfield.allocate(&grid);
  myfield.setAllValuesToZero();

  myfield.boundary_conditions();
  //myfield.smooth_filter(10);

  current.allocate(&grid);
  current.setAllValuesToZero();

  //*******************************************END FIELD DEFINITION***********************************************************

  //*******************************************BEGIN SPECIES DEFINITION*********************************************************
  PLASMA plasma1;
  plasma1.density_function = box;
  plasma1.setMinBox(-10.0, -10.0, grid.rmin[2]);
  plasma1.setMaxBox(10.0, 10.0, grid.rmax[2]);
  plasma1.setRampLength(0.2);
  plasma1.setDensityCoefficient(1.0);
  plasma1.setRampMinDensity(0.001);

  SPECIE  electrons1(&grid);
  electrons1.plasma = plasma1;
  electrons1.setParticlesPerCellXYZ(3, 3, 3);
  electrons1.setName("ELE1");
  electrons1.type = ELECTRON;
  electrons1.creation();
  species.push_back(&electrons1);


  SPECIE electrons2(&grid);
  electrons2.plasma = plasma1;
  electrons2.setParticlesPerCellXYZ(3, 3, 3);
  electrons2.setName("ELE2");
  electrons2.type = ELECTRON;
  electrons2.creation();
  species.push_back(&electrons2);


  tempDistrib distribution;
  distribution.setMaxwell(1.0e-5);

  electrons1.add_momenta(rng, 0.0, 0.0, -1.0, distribution);
  electrons2.add_momenta(rng, 0.0, 0.0, 1.0, distribution);

  for (spec_iterator = species.begin(); spec_iterator != species.end(); spec_iterator++){
    (*spec_iterator)->printParticleNumber();
  }
  //*******************************************END SPECIES DEFINITION***********************************************************

  //*******************************************BEGIN DIAG DEFINITION**************************************************
  OUTPUT_MANAGER manager(&grid, &myfield, &current, species);

  double startOutputA=0.0, freqOutputA=5.0;
  double startOutputB=0.0, freqOutputB=1.0;

  manager.addDiagFrom(startOutputB, freqOutputB);

  manager.addEFieldFrom(startOutputA, freqOutputA);
  manager.addBFieldFrom(startOutputA, freqOutputA);

  manager.addSpeciesDensityFrom("ELE1", startOutputA, freqOutputA);
  manager.addSpeciesDensityFrom("ELE2", startOutputA, freqOutputA);

  manager.addCurrentFrom(startOutputA, freqOutputA);

  manager.addSpeciesPhaseSpaceFrom("ELE1", startOutputA, freqOutputA);
  manager.addSpeciesPhaseSpaceFrom("ELE2", startOutputA, freqOutputA);


  manager.initialize(DIRECTORY_OUTPUT);
  //*******************************************END DIAG DEFINITION**************************************************

  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ MAIN CYCLE (DO NOT MODIFY) @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
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
    // manager.callDiags(grid.istep);

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
