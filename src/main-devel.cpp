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

void moveParticles(GRID* grid, SPECIE* specie, double amplitude,double lambda){
  int Npart = specie->Np;
  double kdx = 2*M_PI/lambda;
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
    specie->u0(n) += deltaV*sin(kdx*oldX+2*M_PI*sqrt(density)*grid->dt*0.5);
  }
}

void moveParticles(GRID* grid, SPECIE* specie, std::vector<KMODE> myKModes){
  int Npart=specie->Np;
  double density=specie->plasma.params.density_coefficient;
  double kk, phi, dr;
  double kx, dx, dVx, oldx;
  double ky, dy, dVy, oldy;
  double kz, dz, dVz, oldz;

  for(int n=0;n<Npart;n++){
    oldx = specie->r0(n);
    oldy = specie->r1(n);
    oldz = specie->r2(n);
    for(int m=0; m < myKModes.size(); m++){
      kx  = myKModes[m].k[0];
      ky  = myKModes[m].k[1];
      kz  = myKModes[m].k[2];
      kk  = sqrt(kx*kx + ky*ky + kz*kz);
      if(fabs(kk)>1e-2){
        dr  = myKModes[m].amplitude/kk;
        phi = myKModes[m].phase;

        dx  = dr*kx/kk;
        dy  = dr*ky/kk;
        dz  = dr*kz/kk;
        dVx  = dx*2*M_PI*sqrt(density);
        dVy  = dy*2*M_PI*sqrt(density);
        dVz  = dz*2*M_PI*sqrt(density);

        phi += kx*oldx + ky*oldy + kz*oldz;
        specie->r0(n) += dx*cos(phi);
        specie->r1(n) += dy*cos(phi);
        specie->r2(n) += dz*cos(phi);

        phi += 2*M_PI*sqrt(density)*grid->dt*0.5;
        specie->u0(n) += dVx*sin(phi);
        specie->u1(n) += dVy*sin(phi);
        specie->u2(n) += dVz*sin(phi);
      }
    }
  }
}

void deformEx(GRID* grid, EM_FIELD* field, double amplitude, double lambda){
  double kdx=2*M_PI/lambda;
  int Ngrid[3];

  Ngrid[0] = grid->Nloc[0];
  Ngrid[1] = grid->Nloc[1];
  Ngrid[2] = grid->Nloc[2];

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

  srand(time(NULL));
  grid.initRNG(mt_rng, RANDOM_NUMBER_GENERATOR_SEED);

  grid.finalize();

  jsonParser::setDumpControl(root, &grid);
  grid.visualDiag();

  //********************************************END GRID DEFINITION********************************************************
  //******************** BEGIN TO READ OF user defined INPUT - PARAMETERS ****************************************
  bool isThereSpecial=false;
  bool isThereAmpli=false;
  bool isThereLambda=false;
  bool isWaveOK = false;
  double amplitude;
  double lambda;
  std::stringstream messaggio;
  Json::Value special;
  isThereSpecial=jsonParser::setValue(special,root,"special");
  if(isThereSpecial){
    isThereAmpli  = jsonParser::setDouble(&amplitude, special, "amplitude");
    isThereLambda = jsonParser::setDouble(&lambda,    special, "lambda");
    isWaveOK = isThereAmpli&&isThereLambda;
  }
  if(true){
    messaggio << "siamo sicuri sia tutto OK? " << " isThereSpecial=" << isThereSpecial << std::endl;
    messaggio << "siamo sicuri sia tutto OK? " << " amplitude=" << amplitude << std::endl;
    messaggio << "siamo sicuri sia tutto OK? " << " lambda=" << lambda << std::endl;
    messaggio << "siamo sicuri sia tutto OK? " << " isWaveOK=" << isWaveOK << std::endl;
  }

  Json::Value langmuirSet;
  bool isThereLangmuirSet = jsonParser::setValue(langmuirSet, root, "langmuirSpectrum");
  std::vector<KMODE> myKModes;
  GRIDmodes gridModes;
  if (isThereLangmuirSet) {

    double amplitude=0.0;
    jsonParser::setDouble(&amplitude, langmuirSet, "amplitude");
    double centralK[3];
    centralK[0]=centralK[1]=centralK[2]=0;
    jsonParser::setDouble(&centralK[0], langmuirSet, "centralKx");
    jsonParser::setDouble(&centralK[1], langmuirSet, "centralKy");
    jsonParser::setDouble(&centralK[2], langmuirSet, "centralKz");
    double sigmaK[3];
    sigmaK[0]=sigmaK[1]=sigmaK[2]=0.0;
    jsonParser::setDouble(&sigmaK[0], langmuirSet, "sigmaKx");
    jsonParser::setDouble(&sigmaK[1], langmuirSet, "sigmaKy");
    jsonParser::setDouble(&sigmaK[2], langmuirSet, "sigmaKz");

    messaggio << "siamo sicuri sia tutto OK? " << " isThereLangmuirSet=" << isThereLangmuirSet << std::endl;
    messaggio << "siamo sicuri sia tutto OK? " << " amplitude=" << amplitude << std::endl;
    messaggio << "siamo sicuri sia tutto OK? " << " centralKx=" << centralK[0] << std::endl;
    messaggio << "siamo sicuri sia tutto OK? " << " sigmaKx=" << sigmaK[0] << std::endl;

    UTILITIES::allocateAccessibleKModes(gridModes, grid);
    UTILITIES::writeGridModes(gridModes, grid);

    UTILITIES::setKModesToBeInitialised(myKModes, gridModes, amplitude, centralK, sigmaK);

    //UTILITIES::setKModesToBeInitialised(myKModes, gridModes, amplitude, centralK, sigmaK);
    UTILITIES::exchangeKModesToBeInitialised(myKModes, grid);
    UTILITIES::writeKModesToBeInitialised(myKModes, grid);

  }
  grid.printMessage(messaggio.str());

  //********************  END READ OF "SPECIAL" (user defined) INPUT - PARAMETERS  ****************************************

  //*******************************************BEGIN SPECIES DEFINITION*********************************************************


  std::map<std::string, PLASMA*> plasmas;
  jsonParser::setPlasmas(root, plasmas);
  jsonParser::setSpecies(root, species, plasmas, &grid, mt_rng);

  //*******************************************  START LANGMUIR WAVE  *********************************************************
  if(isThereLangmuirSet){
    int counter=0;
    for (spec_iterator = species.begin(); spec_iterator != species.end(); spec_iterator++){
      moveParticles(&grid,(*spec_iterator),myKModes);
      (*spec_iterator)->position_parallel_pbc();
      (*spec_iterator)->position_parallel_pbc();
      counter++;
    }
  }
  if(isWaveOK){
    int counter=0;
    for (spec_iterator = species.begin(); spec_iterator != species.end(); spec_iterator++){
      moveParticles(&grid,(*spec_iterator),amplitude,lambda);
      (*spec_iterator)->position_parallel_pbc();
      (*spec_iterator)->position_parallel_pbc();
      counter++;
    }
  }
  //*******************************************   END  LANGMUIR WAVE  *********************************************************

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

  std::cout << messaggio.str();



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
