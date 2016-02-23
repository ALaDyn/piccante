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
#include <random>

#define DEFAULT_DIMENSIONALITY 1


#define DIRECTORY_OUTPUT "OUTPUT"

#define RANDOM_NUMBER_GENERATOR_SEED 5489
#define FREQUENCY_STDOUT_STATUS 5
#include "rapidjson/document.h"     // rapidjson's DOM-style API

void moveParticles(GRID* grid, SPECIE* specie, double amplitude,double lambda){
  int Npart=specie->Np;
  double kdx=2*M_PI/lambda;
  double density=specie->plasma.params.density_coefficient;
  double deltaX = amplitude;
  double deltaV = amplitude*2*M_PI*sqrt(density);
  double oldX;

  if(false){
    std::cout<< "sposto le particelle che sono" << Npart << std::endl;
    std::cout<< "density = " << density << std::endl;
    std::cout<< " ===================== "<< std::endl;
  }
  for(int n=0;n<Npart;n++){
    oldX = specie->r0(n);
    specie->r0(n) += deltaX*cos(kdx*oldX);
    //specie->u0(n) += deltaV*sin(kdx*oldX-2*M_PI*sqrt(density)*grid->dt*0.5);
  }

}

void deformEx(GRID* grid, EM_FIELD* field, double amplitude, double lambda){
  double kdx=2*M_PI/lambda;
  int Ngrid[3];
  Ngrid[0] = grid->NGridNodes[0];
  Ngrid[1] = grid->NGridNodes[1];
  Ngrid[2] = grid->NGridNodes[2];

  double x;

  if(false){
    std::cout<< "deformo Ex che ha " << Ngrid[0] << " punti" << std::endl;
    std::cout<< "Ny = " << Ngrid[1] << " " << std::endl;
    std::cout<< "Nz = " << Ngrid[2] << " " << std::endl;
    std::cout<< "amplitude = " << amplitude << " " << std::endl;
    std::cout<< "lambda = " << lambda << " " << std::endl;
  }
  for(int k=0; k<Ngrid[2]; k++){
    for(int j=0; j<Ngrid[1]; j++){
      for(int i=0; i<Ngrid[0]; i++){

        x = grid->rminloc[0] + grid->dr[0]*i;

        field->E0(i,j,k)+=amplitude*M_PI*cos(kdx*x);
      }
    }
  }
}

void poissonTest(GRID* grid, EM_FIELD* field, CURRENT* current){
  int Ngrid[3];
  Ngrid[0] = grid->NGridNodes[0];
  Ngrid[1] = grid->NGridNodes[1];
  Ngrid[2] = grid->NGridNodes[2];

  double x;
  std::ofstream density("density_0.txt");
  std::ofstream ExBefore("Ex-before.txt");
  std::ofstream ExAfter("Ex-after.txt");
  double totalSum=0;

  for(int k=0; k<Ngrid[2]; k++){
    for(int j=0; j<Ngrid[1]; j++){
      for(int i=0; i<Ngrid[0]; i++){
        x = grid->rmin[0] + grid->dr[0]*i;
        density << x << "   " << (1-current->density(i,j,k)) << std::endl;
        ExBefore << x << "   " << field->E0(i,j,k) << std::endl;
      }
    }
  }

  for(int k=0; k<Ngrid[2]; k++){
    for(int j=0; j<Ngrid[1]; j++){
      totalSum=0;
      int i=0;
      totalSum += field->E0(i,j,k) = 0;
      for(i=1; i<Ngrid[0]; i++){
        field->E0(i,j,k) = grid->dr[0]*grid->den_factor*(1-current->density(i,j,k)) + field->E0(i-1,j,k);
        totalSum += field->E0(i,j,k);
      }
      totalSum /= (Ngrid[0]-1);
      for(i=0; i<Ngrid[0]; i++){
        field->E0(i,j,k) -= totalSum;
      }
    }
  }

  for(int k=0; k<Ngrid[2]; k++){
    for(int j=0; j<Ngrid[1]; j++){
      for(int i=0; i<Ngrid[0]; i++){
        x = grid->rmin[0] + grid->dr[0]*i;
        ExAfter << x << "   " << field->E0(i,j,k) << std::endl;

      }
    }
  }
  density.close();
  ExBefore.close();
  ExAfter.close();
}

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
  bool isThereSpecial=false;
  bool isThereAmpli=false;
  bool isThereLambda=false;
  bool isWaveOK = false;
  double amplitude;
  double lambda;

  Json::Value special;
  if(isThereSpecial=jsonParser::setValue(special,root,"special")){
    isThereAmpli  = jsonParser::setDouble(&amplitude, special, "amplitude");
    isThereLambda = jsonParser::setDouble(&lambda,    special, "lambda");
    isWaveOK = isThereAmpli&&isThereLambda;
  }
  if(false){
    std::cout << "siamo sicuri sia tutto OK? " << " isThereSpecial=" << isThereSpecial << std::endl;
    std::cout << "siamo sicuri sia tutto OK? " << " amplitude=" << amplitude << std::endl;
    std::cout << "siamo sicuri sia tutto OK? " << " lambda=" << lambda << std::endl;
    std::cout << "siamo sicuri sia tutto OK? " << " isWaveOK=" << isWaveOK << std::endl;
  }

  //********************  END READ OF "SPECIAL" (user defined) INPUT - PARAMETERS  ****************************************

  //*******************************************BEGIN SPECIES DEFINITION*********************************************************

  std::map<std::string, PLASMA*> plasmas;
  jsonParser::setPlasmas(root, plasmas);
  jsonParser::setSpecies(root, species, plasmas, &grid, mt_rng);

  if(isWaveOK){
    int counter=0;
    for (spec_iterator = species.begin(); spec_iterator != species.end(); spec_iterator++){
      moveParticles(&grid,(*spec_iterator),amplitude,lambda);
      (*spec_iterator)->position_parallel_pbc();
      double ampliEx;
      ampliEx = 4*M_PI*(*spec_iterator)->chargeSign*(*spec_iterator)->Z*(*spec_iterator)->plasma.params.density_coefficient*amplitude;

      //current.eraseDensity();
      //(*spec_iterator)->density_deposition_standard(&current);
      //current.pbc();
      //poissonTest(&grid,&myfield,&current);
      //std::cout << "counter= " << counter<< "  ampliEx = " << ampliEx << "  lambda = " << lambda << std::endl;
      deformEx(&grid,&myfield,ampliEx,lambda);
      myfield.boundary_conditions();
      counter++;
    }
  }

  uint64_t totPartNum = 0;
  for (spec_iterator = species.begin(); spec_iterator != species.end(); spec_iterator++) {
    totPartNum += (*spec_iterator)->printParticleNumber();
  }
  if (grid.myid == grid.master_proc) {
    std::cout << "Total particle number: " << totPartNum << std::endl;
  }

  //*******************************************END SPECIES DEFINITION***********************************************************

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
    restartFromDump(&dumpID, &grid, &myfield, species);
  }

  while (grid.istep <= grid.getTotalNumberOfTimesteps())
  {

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
