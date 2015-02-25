

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


#define DIRECTORY_OUTPUT "OUTPUT"
#define DIRECTORY_DUMP "DUMP"
#define RANDOM_NUMBER_GENERATOR_SEED 5489
#define FREQUENCY_STDOUT_STATUS 5
#include "rapidjson/document.h"     // rapidjson's DOM-style API

void readAndAllocateSpheres(SPHERES &spheres, std::string filename){

  std::ifstream myfile;
  myfile.open(filename.c_str());
  if (!myfile.good()){
//    std::cout << "     disonore e tempesta! file \"" << filename << "\" does not exists" << std::endl;
  }
  else{
  //  std::cout << " file exists "<< filename << std::endl;

  }
  int Nsph=0;

  myfile.read((char*)&Nsph, sizeof(int));

  spheres.NSpheres = Nsph;
  //std::cout << "     Nsph = " << Nsph << std::endl;
  spheres.coords = new float[spheres.NSpheres*4];

  myfile.read((char*)&(spheres.fillingFactor), sizeof(float));
  //std::cout << "     fillingFactor = " << spheres.fillingFactor << std::endl;
  myfile.read((char*)spheres.rmin, sizeof(float)*3);
  myfile.read((char*)spheres.rmax, sizeof(float)*3);
  myfile.read((char*)spheres.coords, sizeof(float)*spheres.NSpheres*4);
  myfile.close();

}


int main(int narg, char **args)
{
  MPI_Init(&narg, &args);

  Json::Value root;
  jsonParser::parseJsonInputFile(root,narg, args);
  int dim = jsonParser::getDimensionality(root, DEFAULT_DIMENSIONALITY);

  GRID grid(dim);
  EM_FIELD myfield;
  CURRENT current;
  std::vector<SPECIE*> species;
  std::vector<SPECIE*>::const_iterator spec_iterator;
  gsl_rng* rng = gsl_rng_alloc(gsl_rng_ranlxd1);

  //*******************************************BEGIN GRID DEFINITION*******************************************************
  jsonParser::setXrange(root,&grid);
  jsonParser::setYrange(root,&grid);
  jsonParser::setZrange(root,&grid);
  jsonParser::setNCells(root,&grid);
  jsonParser::setNprocs(root,&grid);
  jsonParser::setStretchedGrid(root,&grid);
  jsonParser::setBoundaryConditions(root, &grid);

  jsonParser::setRadiationFriction(root, &grid);
  jsonParser::setMasterProc(root, &grid);

  grid.mpi_grid_initialize(&narg, args);

  jsonParser::setCourantFactor(root, &grid);
  jsonParser::setSimulationTime(root,&grid);
  jsonParser::setMovingWindow(root,&grid);

  srand(time(NULL));
  grid.initRNG(rng, RANDOM_NUMBER_GENERATOR_SEED);

  grid.finalize();

  jsonParser::setDumpControl(root, &grid);
  grid.visualDiag();

  //********************************************END GRID DEFINITION********************************************************
  //*******************************************BEGIN FIELD DEFINITION*********************************************************
  myfield.allocate(&grid);
  myfield.setAllValuesToZero();


  jsonParser::setLaserPulses(root, &myfield);
  myfield.boundary_conditions();

  current.allocate(&grid);
  current.setAllValuesToZero();
  //*******************************************END FIELD DEFINITION***********************************************************
  //******************** BEGIN TO READ OF user defined INPUT - PARAMETERS ****************************************
  bool isThereSpecial=false;
  bool areThereSpheres=false;

  std::string fileSpheresName;
  Json::Value special;
  SPHERES myspheres;
  if(isThereSpecial=jsonParser::setValue(special,root,"special")){
    if(areThereSpheres = jsonParser::setString( &fileSpheresName, special, "spheresFile")){
     readAndAllocateSpheres(myspheres, fileSpheresName);
    }
  }
  std::map<std::string, PLASMA*>::iterator pIterator;

  //********************  END READ OF "SPECIAL" (user defined) INPUT - PARAMETERS  ****************************************

  //*******************************************BEGIN SPECIES DEFINITION*********************************************************

  std::map<std::string, PLASMA*> plasmas;
  jsonParser::setPlasmas(root, plasmas);

  if(areThereSpheres){
    for(pIterator = plasmas.begin(); pIterator != plasmas.end(); pIterator++){
       (pIterator)->second->params.spheres = &myspheres;
     }
  }
 jsonParser::setSpecies(root, species, plasmas, &grid, rng);

  for (spec_iterator = species.begin(); spec_iterator != species.end(); spec_iterator++){
    (*spec_iterator)->printParticleNumber();
  }
  if(areThereSpheres){
    delete [] myspheres.coords;
  }
  //*******************************************END SPECIES DEFINITION***********************************************************

  //*******************************************BEGIN DIAG DEFINITION**************************************************
  std::map<std::string, outDomain*> outDomains;
  OUTPUT_MANAGER manager(&grid, &myfield, &current, species);
  jsonParser::setDomains(root, outDomains);
  jsonParser::setOutputRequests(root, manager, outDomains, species);
  jsonParser::setOutputDirPath(root,manager);

  manager.initialize();
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

    if(!(grid.istep%20))
      myfield.applyFilter(fltr_Ex|fltr_Ey, dir_x|dir_y);

    myfield.new_halfadvance_B();
    myfield.boundary_conditions();

    for (spec_iterator = species.begin(); spec_iterator != species.end(); spec_iterator++){
      if(grid.isRadiationFrictionEnabled()){
        (*spec_iterator)->momenta_advance_with_friction(&myfield, grid.getLambda0());
      }
      else{
        (*spec_iterator)->momenta_advance(&myfield);
      }
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
