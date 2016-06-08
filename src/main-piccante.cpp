/*   Copyright 2014-2016 - Andrea Sgattoni, Luca Fedeli, Stefano Sinigardi   */

/******************************************************************************
* This file is part of piccante.                                              *
*                                                                             *
* piccante is free software: you can redistribute it and/or modify            *
* it under the terms of the GNU General Public License as published by        *
* the Free Software Foundation, either version 3 of the License, or           *
* (at your option) any later version.                                         *
*                                                                             *
* piccante is distributed in the hope that it will be useful,                 *
* but WITHOUT ANY WARRANTY; without even the implied warranty of              *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the                *
* GNU General Public License for more details.                                *
*                                                                             *
* You should have received a copy of the GNU General Public License           *
* along with piccante. If not, see <http://www.gnu.org/licenses/>.            *
******************************************************************************/

#define _USE_MATH_DEFINES

#include <mpi.h>
#include "commons.h"
#include "grid.h"
#include "structures.h"
#include "current.h"
#include "em_field.h"
#include "particle_species.h"
#include "output_manager.h"
#include "utilities.h"
#include "jsonparser.h"

#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <iomanip>
#include <cstring>
#include <ctime>
#include <cstdarg>
#include <vector>
#include <map>

#define DEFAULT_DIMENSIONALITY 1


#define DIRECTORY_OUTPUT "OUTPUT"

#define RANDOM_NUMBER_GENERATOR_SEED 5489

#include "rapidjson/document.h"     // rapidjson's DOM-style API


int main(int narg, char **args)
{
  MPI_Init(&narg, &args);

  Json::Value root;
  std::string inputFileName = jsonParser::parseJsonInputFile(root, narg, args);
  int dim = jsonParser::getDimensionality(root, DEFAULT_DIMENSIONALITY);

  GRID grid(dim);
  EM_FIELD myfield;
  CURRENT current;
  std::vector<SPECIE*> species;
  std::vector<SPECIE*>::const_iterator spec_iterator;
  my_rng_generator mt_rng;
  //*******************************************BEGIN GRID DEFINITION*******************************************************
  jsonParser::setXrange(root, &grid);
  jsonParser::setYrange(root, &grid);
  jsonParser::setZrange(root, &grid);
  jsonParser::setNCells(root, &grid);
  jsonParser::setNprocs(root, &grid);
  jsonParser::setStretchedGrid(root, &grid);
  jsonParser::setBoundaryConditions(root, &grid);

  jsonParser::setRadiationFriction(root, &grid);
  jsonParser::setMasterProc(root, &grid);
  grid.mpi_grid_initialize(&narg, args);

  jsonParser::setCourantFactor(root, &grid);
  jsonParser::setSimulationTime(root, &grid);
  jsonParser::setMovingWindow(root, &grid);
  jsonParser::setFrequencyStdoutStatus(root, &grid);

  //srand(time(NULL));
  grid.initRNG(mt_rng, RANDOM_NUMBER_GENERATOR_SEED);

  grid.finalize();

  jsonParser::setDumpControl(root, &grid);
  grid.visualDiag();

  //********************************************END GRID DEFINITION********************************************************
   //******************** BEGIN TO READ OF user defined INPUT - PARAMETERS ****************************************
//  int myIntVariable = 0;
//  double myDoubleVariable = 0;
//  bool isThereSpecial = false;
//  Json::Value special= jsonParser::setValue(special, root, "special");
//  if (isThereSpecial) {
//    jsonParser::setInt(&myIntVariable, special, "variabile1");
//    jsonParser::setDouble(&myDoubleVariable, special, "variabile2");
//  }

   bool isThereSpecial = false;
   bool areThereSpheres = false;

  std::string fileSpheresName;
  Json::Value special;
  SPHERES myspheres;
  isThereSpecial = jsonParser::setValue(special, root, "special");
  if (isThereSpecial) {
      areThereSpheres = jsonParser::setString(&fileSpheresName, special, "spheresFile");
    if (areThereSpheres) {
      UTILITIES::readAndAllocateSpheres(myspheres, fileSpheresName, grid);
      UTILITIES::selectSpheres(myspheres, grid);
    }
  }
  std::map<std::string, PLASMA*>::iterator pIterator;


  //********************  END READ OF "SPECIAL" (user defined) INPUT - PARAMETERS  ****************************************

  //*******************************************BEGIN SPECIES DEFINITION*********************************************************

  std::map<std::string, PLASMA*> plasmas;
  jsonParser::setPlasmas(root, plasmas);
  if (areThereSpheres) {
    for (pIterator = plasmas.begin(); pIterator != plasmas.end(); pIterator++) {
      (pIterator)->second->params.spheres = &myspheres;
    }
  }
  jsonParser::setSpecies(root, species, plasmas, &grid, mt_rng);

  uint64_t totPartNum = 0;
  for (spec_iterator = species.begin(); spec_iterator != species.end(); spec_iterator++) {
    totPartNum += (*spec_iterator)->printParticleNumber();
  }
  if (grid.myid == grid.master_proc) {
    std::cout << "Total particle number: " << totPartNum << std::endl;
  }

  //*******************************************END SPECIES DEFINITION***********************************************************
  //*******************************************BEGIN FIELD DEFINITION*********************************************************
  myfield.allocate(&grid);
  myfield.setAllValuesToZero();
  current.allocate(&grid);

  //*******************************************    POISSON SOLVER    *********************************************************
  jsonParser::setPoissonSolver(root, &grid);

  if(grid.isWithPoisson()){
    bool withSign = true;
    std::cout << " evaluating density..." << std::endl;

    current.setAllValuesToZero();
    for (spec_iterator = species.begin(); spec_iterator != species.end(); spec_iterator++) {
      (*spec_iterator)->density_deposition_standard(&current, withSign);
    }
    current.pbc();
    std::cout << "   done... now into Poisson solver" << std::endl;
    myfield.poissonSolver(&current);
  }

  jsonParser::setLaserPulses(root, &myfield);
  myfield.boundary_conditions();


  current.setAllValuesToZero();
  //*******************************************END FIELD DEFINITION***********************************************************

  //*******************************************BEGIN DIAG DEFINITION**************************************************
  std::map<std::string, outDomain*> outDomains;
  OUTPUT_MANAGER manager(&grid, &myfield, &current, species);
  jsonParser::setDomains(root, outDomains);
  jsonParser::setOutputRequests(root, manager, outDomains, species);
  jsonParser::setOutputDirPath(root, manager);
  jsonParser::setOutputParameters(root, manager);

  manager.initialize();
  manager.copyInputFileInOutDir(inputFileName);
  //*******************************************END DIAG DEFINITION**************************************************

  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ MAIN CYCLE (DO NOT MODIFY) @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  if (grid.myid == grid.master_proc) {
    printf("----- START temporal cicle -----\n");
    fflush(stdout);
  }

  int dumpID = 1;
  grid.istep = 0;
  if (grid.dumpControl.doRestart) {
    dumpID = grid.dumpControl.restartFromDump;
    UTILITIES::restartFromDump(&dumpID, &grid, &myfield, species);
  }

  while (grid.istep <= grid.getTotalNumberOfTimesteps())
  {

    grid.printTStepAsPlanned();

    manager.callDiags(grid.istep);

    myfield.openBoundariesE_1();
    myfield.new_halfadvance_B();
    myfield.boundary_conditions();

    current.setAllValuesToZero();
    for (spec_iterator = species.begin(); spec_iterator != species.end(); spec_iterator++) {
      (*spec_iterator)->current_deposition_standard(&current);
    }
    current.pbc();

    for (spec_iterator = species.begin(); spec_iterator != species.end(); spec_iterator++) {
      (*spec_iterator)->position_parallel_pbc();
    }

    myfield.openBoundariesB();
    myfield.new_advance_E(&current);

    myfield.boundary_conditions();
    myfield.openBoundariesE_2();
    myfield.new_halfadvance_B();
    myfield.boundary_conditions();

    for (spec_iterator = species.begin(); spec_iterator != species.end(); spec_iterator++) {
      if (grid.isRadiationFrictionEnabled()) {
        (*spec_iterator)->momenta_advance_with_friction(&myfield, grid.getLambda0());
      }
      else {
        (*spec_iterator)->momenta_advance(&myfield);
      }
    }

    grid.time += grid.dt;

    UTILITIES::moveWindow(&grid, &myfield, species);

    grid.istep++;
    if (grid.dumpControl.doDump) {
      if (grid.istep != 0 && !(grid.istep % ((int)(grid.dumpControl.dumpEvery / grid.dt)))) {
        UTILITIES::dumpFilesForRestart(&dumpID, &grid, &myfield, species);
      }
    }
  }

  manager.close();
  MPI_Finalize();
  exit(0);

}
