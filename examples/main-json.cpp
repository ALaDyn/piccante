
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

#include "commons.h"
#include "grid.h"
#include "structures.h"
#include "current.h"
#include "em_field.h"
#include "particle_species.h"
#include "output_manager.h"
#include "utilities.h"


#define DEFAULT_DIMENSIONALITY 1



#define DIRECTORY_OUTPUT "TEST"
#define DIRECTORY_DUMP "DUMP"
#define RANDOM_NUMBER_GENERATOR_SEED 5489
#define FREQUENCY_STDOUT_STATUS 5
#include "rapidjson/document.h"     // rapidjson's DOM-style API
//#include "rapidjson/prettywriter.h" // for stringify JSON
//#include "rapidjson/filestream.h"   // wrapper of C stream for prettywriter as output


struct mySpecialParameters{
  static const int Nint=3, Ndouble=3, Nbool=3;
  int paramI[Nint];
  double paramD[Ndouble];
  bool paramB[Nbool];
//const char* nam="cacca";
//  static const char* namesI[]={"int1","int2","int3"};
//  static const char* NamesD[]={"double1", "double2", "double3"};
//  static const char* NamesB[]={"bool1","bool2", "bool3"};
};

int main(int narg, char **args)
{
  rapidjson::Document document;
  parseJsonInputFile(document,"inputPiccante.json");
  int dim = getDimensionalityFromJson(document, DEFAULT_DIMENSIONALITY);
  GRID grid(dim);
  EM_FIELD myfield;
  CURRENT current;
  std::vector<SPECIE*> species;
  std::vector<SPECIE*>::const_iterator spec_iterator;
  gsl_rng* rng = gsl_rng_alloc(gsl_rng_ranlxd1);

  //*******************************************BEGIN GRID DEFINITION*******************************************************

  setXrangeFromJson(document,&grid);
  setYrangeFromJson(document,&grid);
  setZrangeFromJson(document,&grid);

  setNCellsFromJson(document,&grid);
  setNprocsFromJson(document,&grid);

  setStretchedGridFromJson(document,&grid);

  grid.setBoundaries(xOpen | yOpen | zPBC);
  grid.mpi_grid_initialize(&narg, args);
  grid.setCourantFactor(0.98);
  setSimulationTimeFromJson(document,&grid);

  setMovingWindowFromJson(document,&grid);


  grid.setMasterProc(0);

  srand(time(NULL));
  grid.initRNG(rng, RANDOM_NUMBER_GENERATOR_SEED);

  grid.finalize();

  setDumpControlFromJson(document, &grid.dumpControl);
  grid.visualDiag();

  //********************************************END GRID DEFINITION********************************************************
  //******************** BEGIN TO READ OF user defined INPUT - PARAMETERS ****************************************
int myIntVariable=0;
  double myDoubleVariable=0;
  bool isThereSpecial=false;
  rapidjson::Value special;
  if(isThereSpecial=setValueFromJson(special,document,"special")){
    setIntFromJson( &myIntVariable, special, "variabile1");
    setDoubleFromJson( &myDoubleVariable, special, "variabile2");
   }
//********************  END READ OF "SPECIAL" (user defined) INPUT - PARAMETERS  ****************************************
  //*******************************************BEGIN FIELD DEFINITION*********************************************************
  myfield.allocate(&grid);
  myfield.setAllValuesToZero();

  laserPulse pulse1;
  pulse1.type = COS2_PLANE_WAVE;
  pulse1.setPPolarization();
  pulse1.setDurationFWHM(5.0);
  pulse1.setPulseInitialPosition(-5.5);
  pulse1.setNormalizedAmplitude(8.0);

  myfield.addPulse(&pulse1);

  myfield.boundary_conditions();

  current.allocate(&grid);
  current.setAllValuesToZero();
  //*******************************************END FIELD DEFINITION***********************************************************

  //*******************************************BEGIN SPECIES DEFINITION*********************************************************
  PLASMA plasma1;
  plasma1.density_function = left_soft_ramp;
  plasma1.setXRangeBox(0.0, 1.5);
  plasma1.setYRangeBox(grid.rmin[1], grid.rmax[1]);
  plasma1.setZRangeBox(grid.rmin[2], grid.rmax[2]);
  plasma1.setRampLength(0.5);
  plasma1.setDensityCoefficient(80);
  plasma1.setRampMinDensity(0.0);


  PLASMA plasma2;
  plasma2.density_function = box;
  plasma2.setXRangeBox(1.5, 1.505);
  plasma2.setYRangeBox(grid.rmin[1], grid.rmax[1]);
  plasma2.setZRangeBox(grid.rmin[2], grid.rmax[2]);
  plasma2.setRampLength(0.5);
  plasma2.setDensityCoefficient(10);
  plasma2.setRampMinDensity(0.0);


  SPECIE  electrons1(&grid);
  electrons1.plasma = plasma1;
  electrons1.setParticlesPerCellXYZ(300, 1, 1);
  electrons1.setName("ELE1");
  electrons1.type = ELECTRON;
  electrons1.creation();
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
  distribution.setMaxwell(1.0e-5);

  electrons1.add_momenta(rng, 0.0, 0.0, 0.0, distribution);
  ions1.add_momenta(rng, 0.0, 0.0, 0.0, distribution);

  for (spec_iterator = species.begin(); spec_iterator != species.end(); spec_iterator++){
    (*spec_iterator)->printParticleNumber();
  }

  //*******************************************END SPECIES DEFINITION***********************************************************

  //*******************************************BEGIN DIAG DEFINITION**************************************************

  OUTPUT_MANAGER manager(&grid, &myfield, &current, species);

  manager.addEFieldFrom(0.0, 2.0);
  manager.addBFieldFrom(0.0, 2.0);

  manager.addSpeciesDensityFrom(electrons1.name, 0.0, 2.0);
  manager.addSpeciesDensityFrom(ions1.name, 0.0, 2.0);

  manager.addCurrentFrom(0.0, 5.0);

  manager.addSpeciesPhaseSpaceFrom(electrons1.name, 0.0, 5.0);
  manager.addSpeciesPhaseSpaceFrom(ions1.name, 0.0, 5.0);

  manager.addDiagFrom(0.0, 1.0);

  manager.initialize(DIRECTORY_OUTPUT);

  //*******************************************END DIAG DEFINITION**************************************************
grid.setDumpPath(DIRECTORY_DUMP);
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ MAIN CYCLE (DO NOT MODIFY) @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  if (grid.myid == grid.master_proc){
    printf("----- START temporal cicle -----\n");
    fflush(stdout);
  }

  int dumpID = 1;
  grid.istep = 0;
  if (grid.dumpControl.doRestart){
    dumpID = grid.dumpControl.restartFromDump;
    std::cout << "restartID = " << dumpID << "\n";
    restartFromDump(&dumpID, &grid, &myfield, species);
  }
  int Nstep = grid.getTotalNumberOfTimesteps();
  while (grid.istep <= Nstep)
  {
#ifdef NO_ALLOCATION
    manager.close();
    MPI_Finalize();
    exit(0);
#endif

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
    if (grid.dumpControl.doDump){
      if (grid.istep != 0 && !(grid.istep % ((int)(grid.dumpControl.dumpEvery / grid.dt)))) {
        dumpFilesForRestart(&dumpID, &grid, &myfield, species);
      }
    }
  }

  manager.close();
  MPI_Finalize();
  exit(0);

}
