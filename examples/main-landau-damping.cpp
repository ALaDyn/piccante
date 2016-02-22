/*   Copyright 2016 - Andrea Sgattoni, Luca Fedeli, Stefano Sinigardi        */

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
#ifdef _USE_FFTW_FILTER
#include<fftw3-mpi.h>
#endif

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

void readAndAllocateSpheres(SPHERES &spheres, std::string filename, GRID &grid) {

  if (grid.myid == grid.master_proc) {
    std::ifstream myfile;
    myfile.open(filename.c_str());
    if (!myfile.good()) {
      std::cout << "     spheres file: \"" << filename << "\" does not exists" << std::endl;
    }
    int Nsph = 0;

    myfile.read((char*)&Nsph, sizeof(int));

    spheres.NSpheres = Nsph;
    spheres.coords = new float[spheres.NSpheres * 4];

    myfile.read((char*)&(spheres.fillingFactor), sizeof(float));
    myfile.read((char*)spheres.rmin, sizeof(float) * 3);
    myfile.read((char*)spheres.rmax, sizeof(float) * 3);
    myfile.read((char*)spheres.coords, sizeof(float)*spheres.NSpheres * 4);
    myfile.close();
  }


  MPI_Bcast(&(spheres.NSpheres), 1, MPI_INT, grid.master_proc, MPI_COMM_WORLD);

  if (grid.myid != grid.master_proc)
    spheres.coords = new float[spheres.NSpheres * 4];

  MPI_Bcast(&(spheres.fillingFactor), 1, MPI_FLOAT, grid.master_proc, MPI_COMM_WORLD);
  MPI_Bcast(spheres.rmin, 3, MPI_FLOAT, grid.master_proc, MPI_COMM_WORLD);
  MPI_Bcast(spheres.rmax, 3, MPI_FLOAT, grid.master_proc, MPI_COMM_WORLD);
  MPI_Bcast(spheres.coords, 4 * spheres.NSpheres, MPI_FLOAT, grid.master_proc, MPI_COMM_WORLD);

}

void fromCoordsToSpheresCoords(double &x, double min, double max) {
  double box = max - min;
  if (x < min) {
    x += box * ((int)((max - x) / box));
  }
  if (x > max) {
    x -= box* ((int)((x - min) / box));
  }
}

bool isSphereInside(SPHERES& spheres, int index, GRID &grid) {
  float x = spheres.coords[index * 4];
  float y = spheres.coords[index * 4 + 1];
  float z = spheres.coords[index * 4 + 2];
  float r = spheres.coords[index * 4 + 3];

  double xl = grid.rminloc[0] - r;
  double yl = grid.rminloc[1] - r;
  double zl = grid.rminloc[2] - r;

  double xr = grid.rmaxloc[0] + r;
  double yr = grid.rmaxloc[1] + r;
  double zr = grid.rmaxloc[2] + r;

  fromCoordsToSpheresCoords(yl, spheres.rmin[1], spheres.rmax[1]);
  fromCoordsToSpheresCoords(zl, spheres.rmin[2], spheres.rmax[2]);
  fromCoordsToSpheresCoords(yr, spheres.rmin[1], spheres.rmax[1]);
  fromCoordsToSpheresCoords(zr, spheres.rmin[2], spheres.rmax[2]);

  bool chkX, chkY, chkZ;
  chkX = chkY = chkZ = true;

  if (grid.getDimensionality() >= 2) {
    if ((grid.rmaxloc[1] - grid.rminloc[1]) >= (spheres.rmax[1] - spheres.rmin[1]))
      chkY = true;
    else {
      if (yl <= yr)
        chkY = ((y >= yl) && (y <= yr));
      else
        chkY = ((y >= yr) || (y <= yl));
    }
  }

  if (grid.getDimensionality() == 3) {
    if ((grid.rmaxloc[2] - grid.rminloc[2]) >= (spheres.rmax[2] - spheres.rmin[2]))
      chkZ = true;
    else {

      if (zl <= zr)
        chkZ = ((z >= zl) && (z <= zr));
      else
        chkZ = ((z >= zr) || (z <= zl));
    }
  }

  return chkX && chkY && chkZ;

}

void swapSpheres(SPHERES &spheres, int i, int j) {
  float dummyCoords[4];
  dummyCoords[0] = spheres.coords[i * 4];
  dummyCoords[1] = spheres.coords[i * 4 + 1];
  dummyCoords[2] = spheres.coords[i * 4 + 2];
  dummyCoords[3] = spheres.coords[i * 4 + 3];

  spheres.coords[i * 4] = spheres.coords[j * 4];
  spheres.coords[i * 4 + 1] = spheres.coords[j * 4 + 1];
  spheres.coords[i * 4 + 2] = spheres.coords[j * 4 + 2];
  spheres.coords[i * 4 + 3] = spheres.coords[j * 4 + 3];

  spheres.coords[j * 4] = dummyCoords[0];
  spheres.coords[j * 4 + 1] = dummyCoords[1];
  spheres.coords[j * 4 + 2] = dummyCoords[2];
  spheres.coords[j * 4 + 3] = dummyCoords[3];
}

void selectSpheres(SPHERES &spheres, GRID &grid) {
  int counter = 0;
  for (int i = 0; i < spheres.NSpheres; i++) {
    if (isSphereInside(spheres, i, grid)) {
      swapSpheres(spheres, i, counter);
      counter++;
    }
  }
  spheres.NSpheres = counter;
  spheres.coords = (float*)realloc(spheres.coords, counter * 4 * sizeof(float));
}


int main(int narg, char **args)
{
  MPI_Init(&narg, &args);
#ifdef _USE_FFTW_FILTER
  fftw_mpi_init();
#endif
  Json::Value root;
  jsonParser::parseJsonInputFile(root, narg, args);
  int dim = jsonParser::getDimensionality(root, DEFAULT_DIMENSIONALITY);

  GRID grid(dim);
  EM_FIELD myfield;
  CURRENT current;
  std::vector<SPECIE*> species;
  std::vector<SPECIE*>::const_iterator spec_iterator;
  std::mt19937 mt_rng;

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

  srand(time(NULL));
  grid.initRNG(mt_rng, RANDOM_NUMBER_GENERATOR_SEED);

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
  bool isThereSpecial = false;
  bool areThereSpheres = false;

  std::string fileSpheresName;
  Json::Value special;
  SPHERES myspheres;
  if (isThereSpecial = jsonParser::setValue(special, root, "special")) {
    if (areThereSpheres = jsonParser::setString(&fileSpheresName, special, "spheresFile")) {
      readAndAllocateSpheres(myspheres, fileSpheresName, grid);
      selectSpheres(myspheres, grid);
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


  if (areThereSpheres) {
    delete[] myspheres.coords;
  }
  //*******************************************END SPECIES DEFINITION***********************************************************

  //*******************************************BEGIN DIAG DEFINITION**************************************************
  std::map<std::string, outDomain*> outDomains;
  OUTPUT_MANAGER manager(&grid, &myfield, &current, species);
  jsonParser::setDomains(root, outDomains);
  jsonParser::setOutputRequests(root, manager, outDomains, species);
  jsonParser::setOutputDirPath(root, manager);

  manager.initialize();
  //*******************************************END DIAG DEFINITION**************************************************
  grid.setDumpPath(DIRECTORY_DUMP);
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ MAIN CYCLE (DO NOT MODIFY) @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  if (grid.myid == grid.master_proc) {
    printf("----- START temporal cicle -----\n");
    fflush(stdout);
  }

  int dumpID = 1;
  grid.istep = 0;
  if (grid.dumpControl.doRestart) {
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

#ifdef _USE_FFTW_FILTER
    myfield.fftw_filter_Efield();
    myfield.boundary_conditions();
#endif

    myfield.openBoundariesE_2();

    if (!(grid.istep % 20)) {
      //myfield.applyFilter(fltr_Ex|fltr_Ey, dir_x|dir_y);
    }

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

    moveWindow(&grid, &myfield, species);

    grid.istep++;
    if (grid.dumpControl.doDump) {
      if (grid.istep != 0 && !(grid.istep % ((int)(grid.dumpControl.dumpEvery / grid.dt)))) {
        dumpFilesForRestart(&dumpID, &grid, &myfield, species);
      }
    }
  }

  manager.close();
  MPI_Finalize();
  exit(0);

}
