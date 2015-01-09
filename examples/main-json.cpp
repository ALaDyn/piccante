
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
#include <map>

#include "commons.h"
#include "grid.h"
#include "structures.h"
#include "current.h"
#include "em_field.h"
#include "particle_species.h"
#include "output_manager.h"
#include "utilities.h"
#include "jsonparser.h"

#define DEFAULT_DIMENSIONALITY 1


#define DIRECTORY_OUTPUT "TEST"
#define DIRECTORY_DUMP "DUMP"
#define RANDOM_NUMBER_GENERATOR_SEED 5489
#define FREQUENCY_STDOUT_STATUS 5
#include "rapidjson/document.h"     // rapidjson's DOM-style API


int main(int narg, char **args)
{
    const int masterProc=0;

  Json::Value root;
  jsonParser::parseJsonInputFile(root,"inputPiccante.json");
  int dim = jsonParser::getDimensionalityFromJson(root, DEFAULT_DIMENSIONALITY);

  GRID grid(dim);
  EM_FIELD myfield;
  CURRENT current;
  std::vector<SPECIE*> species;
  std::vector<SPECIE*>::const_iterator spec_iterator;
  gsl_rng* rng = gsl_rng_alloc(gsl_rng_ranlxd1);

  //*******************************************BEGIN GRID DEFINITION*******************************************************
  jsonParser::setXrangeFromJson(root,&grid);
  jsonParser::setYrangeFromJson(root,&grid);
  jsonParser::setZrangeFromJson(root,&grid);
  jsonParser::setNCellsFromJson(root,&grid);
  jsonParser::setNprocsFromJson(root,&grid);
  jsonParser::setStretchedGridFromJson(root,&grid);
  jsonParser::setBoundaryConditionsFromJson(root, &grid);

  grid.mpi_grid_initialize(&narg, args);
  if(grid.myid==masterProc)
      jsonParser::isThisJsonMaster=true;
  else
      jsonParser::isThisJsonMaster=false;

  grid.setCourantFactor(0.98);
  jsonParser::setSimulationTimeFromJson(root,&grid);
  jsonParser::setMovingWindowFromJson(root,&grid);

  grid.setMasterProc(masterProc);

  srand(time(NULL));
  grid.initRNG(rng, RANDOM_NUMBER_GENERATOR_SEED);

  grid.finalize();

  jsonParser::setDumpControlFromJson(root, &grid.dumpControl);
  grid.visualDiag();

  //********************************************END GRID DEFINITION********************************************************
  //******************** BEGIN TO READ OF user defined INPUT - PARAMETERS ****************************************
int myIntVariable=0;
  double myDoubleVariable=0;
  bool isThereSpecial=false;
  Json::Value special;
  if(isThereSpecial=jsonParser::setValueFromJson(special,root,"special")){
    jsonParser::setIntFromJson( &myIntVariable, special, "variabile1");
    jsonParser::setDoubleFromJson( &myDoubleVariable, special, "variabile2");
   }
//********************  END READ OF "SPECIAL" (user defined) INPUT - PARAMETERS  ****************************************
  //*******************************************BEGIN FIELD DEFINITION*********************************************************
  myfield.allocate(&grid);
  myfield.setAllValuesToZero();


  jsonParser::setLaserPulsesFromJson(root, &myfield);
  myfield.boundary_conditions();

  current.allocate(&grid);
  current.setAllValuesToZero();
  //*******************************************END FIELD DEFINITION***********************************************************

  //*******************************************BEGIN SPECIES DEFINITION*********************************************************

  std::map<std::string, PLASMA*> plasmas;
  jsonParser::setPlasmasFromJson(root, plasmas);

  jsonParser::setSpeciesFromJson(root, species, plasmas, &grid, rng);

//  PLASMA plasma1;
//  plasma1.density_function = box;
//  plasma1.setMinBox(-10.0, -10.0, grid.rmin[2]);
//  plasma1.setMaxBox(10.0, 10.0, grid.rmax[2]);
//  plasma1.setRampLength(0.2);
//  plasma1.setDensityCoefficient(1.0);
//  plasma1.setRampMinDensity(0.001);

//  SPECIE  electrons1(&grid);
//  electrons1.plasma = plasma1;
//  electrons1.setParticlesPerCellXYZ(3, 3, 3);
//  electrons1.setName("ELE1");
//  electrons1.type = ELECTRON;
//  electrons1.creation();
//  species.push_back(&electrons1);


//  SPECIE electrons2(&grid);
//  electrons2.plasma = plasma1;
//  electrons2.setParticlesPerCellXYZ(3, 3, 3);
//  electrons2.setName("ELE2");
//  electrons2.type = ELECTRON;
//  electrons2.creation();
//  species.push_back(&electrons2);


//  tempDistrib distribution;
//  distribution.setMaxwell(1.0e-5);

//  electrons1.add_momenta(rng, 0.0, 0.0, -1.0, distribution);
//  electrons2.add_momenta(rng, 0.0, 0.0, 1.0, distribution);

  for (spec_iterator = species.begin(); spec_iterator != species.end(); spec_iterator++){
    (*spec_iterator)->printParticleNumber();
  }
  //*******************************************END SPECIES DEFINITION***********************************************************

  //*******************************************BEGIN DIAG DEFINITION**************************************************
  OUTPUT_MANAGER manager(&grid, &myfield, &current, species);

  double startOutputA=0.2, freqOutputA=1.0;
  double startOutputB=0.2, freqOutputB=1.0;

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

  while (grid.istep <= grid.getTotalNumberOfTimesteps())
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
