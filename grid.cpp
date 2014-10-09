/* Copyright 2014 - Andrea Sgattoni, Luca Fedeli, Stefano Sinigardi */

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

#include"grid.h"

GRID::GRID()
{
  std::time(&unix_time_start);  // get current time
  beta_mw = 0;
  time = 0;
  mark_mw = 0;
  master_proc = 0;
  istep = 0;
  with_particles = YES;
  with_current = YES;
  withMovingWindow = false;
  proc_totUniquePoints = NULL;
  cyclic[0] = cyclic[1] = cyclic[2] = 1;
  lambda0 = 1.0;   //set lenght of the normalization
  ref_den = 1.0; //= critical density
  den_factor = (2 * M_PI)*(2 * M_PI);
  dumpPath = "./";
  GRID::initializeStretchParameters();
  rnproc[1]=rnproc[2]=1;
}

GRID::~GRID(){
  delete[] proc_totUniquePoints;
  for (int c = 0; c < 3; c++){
    free(iStretchingDerivativeCorrection[c]);
    free(hStretchingDerivativeCorrection[c]);
    free(rproc_rmin[c]);
    free(rproc_rmax[c]);
    free(rproc_csimin[c]);
    free(rproc_csimax[c]);
    free(rproc_imin[c]);
    free(rproc_imax[c]);
    free(rproc_Nloc[c]);
    free(rproc_NuniquePointsloc[c]);
    free(cir[c]);
    free(chr[c]);
    free(cirloc[c]);
    free(chrloc[c]);
  }
}

void GRID::initializeStretchParameters(){
  flagStretched = false;
  flagStretchedAlong[0] = flagStretchedAlong[1] = flagStretchedAlong[2] = false;
  flagLeftStretchedAlong[0] = flagLeftStretchedAlong[1] = flagLeftStretchedAlong[2] = false;
  flagRightStretchedAlong[0] = flagRightStretchedAlong[1] = flagRightStretchedAlong[2] = false;
}

void GRID::setMasterProc(int idMasterProc){
  master_proc = idMasterProc;
}

void GRID::setXrange(double min, double max){
  rmin[0] = min;
  rmax[0] = max;
}
void GRID::setYrange(double min, double max){
  rmin[1] = min;
  rmax[1] = max;
  if (accesso.dimensions < 2){
    rmin[1] = -1;
    rmax[1] = +1;
  }
}
void GRID::setZrange(double min, double max){
  rmin[2] = min;
  rmax[2] = max;
  if (accesso.dimensions < 3){
    rmin[2] = -1;
    rmax[2] = +1;
  }
}

void GRID::setNCells(int xcells, int ycells, int zcells){
  NGridNodes[0] = 1 + xcells;
  NGridNodes[1] = 1 + ycells;
  NGridNodes[2] = 1 + zcells;
}

void GRID::setNProcsAlongY(int nprocy){
  rnproc[1] = nprocy;
  if (accesso.dimensions < 2)
    rnproc[1] = 1;
}
void GRID::setNProcsAlongZ(int nprocz){
  rnproc[2] = nprocz;
  if (accesso.dimensions < 3)
    rnproc[2] = 1;
}

void GRID::setCourantFactor(double courant_factor){
  switch (accesso.dimensions){
  case 1:
    dt = courant_factor*(1 / (sqrt(dri[0] * dri[0])));
    break;
  case 2:
    dt = courant_factor*(1 / (sqrt(dri[0] * dri[0] + dri[1] * dri[1])));
    break;
  case 3:
    dt = courant_factor*(1 / (sqrt(dri[0] * dri[0] + dri[1] * dri[1] + dri[2] * dri[2])));
    break;
  default:
    printf("WRONG definition of DIMENSIONALITY\n");
    exit(17);
    break;
  }
}
void GRID::setSimulationTime(double tot_time){
  totalTime = tot_time;
  totalNumberOfTimesteps = (int)(tot_time / dt) + 2;
}

double GRID::getTotalTime(){
  return totalTime;
}

void GRID::setMovingWindow(double start, double beta, int frequency_mw){
  withMovingWindow = true;
  shouldIMove = false;
  beta_mw = beta;
  t_start_moving_mw = start;
  frequency_mw_shifts = frequency_mw;
  mark_mw = 0;
  imove_mw = 0;
  fmove_mw = 0;
}

void GRID::setStartMovingWindow(double start){
  withMovingWindow = true;
  shouldIMove = false;
  beta_mw = 1.0;
  t_start_moving_mw = start;
  frequency_mw_shifts = 20;
  mark_mw = 0;
  imove_mw = 0;
  fmove_mw = 0;
}
void GRID::setBetaMovingWindow(double beta){
  beta_mw = beta;
}

void GRID::setFrequencyMovingWindow(int frequency_mw){
  frequency_mw_shifts = frequency_mw;
}

int GRID::getTotalNumberOfTimesteps(){
  return totalNumberOfTimesteps;
}
void GRID::move_window(){
  int cell_num;
  double buff, move;

  shouldIMove = false;
  imove_mw = 0;
  fmove_mw = 0;

  if (!withMovingWindow)
    return;

  if (time > t_start_moving_mw){

    if (!(istep % frequency_mw_shifts))
    {
      buff = (time - t_start_moving_mw)*beta_mw - mark_mw;
      cell_num = (int)(buff / dr[0]);
      move = cell_num*dr[0];
      mark_mw += move;
      fmove_mw = move;
      imove_mw = cell_num;

      if (imove_mw > 0){
        shouldIMove = true;

        rmin[0] += move;
        rmax[0] += move;
        rminloc[0] += move;
        rmaxloc[0] += move;

        csimax[0] += move;
        csimin[0] += move;
        csiminloc[0] += move;
        csimaxloc[0] += move;

        for (int pp = 0; pp < rnproc[0]; pp++){
          rproc_rmin[0][pp] += move;
          rproc_rmax[0][pp] += move;
        }
        for (int i = 0; i < NGridNodes[0]; i++){
          cir[0][i] += move;
          chr[0][i] += move;
        }
        for (int i = 0; i < Nloc[0]; i++){
          cirloc[0][i] += move;
          chrloc[0][i] += move;
        }
      }

    }
  }
}

void GRID::printTStepEvery(int every){
  int Nstep = totalNumberOfTimesteps;
  if (!(istep % (every)))
  {
    if (myid == master_proc){
      time_t timer;
      std::time(&timer);  /* get current time; same as: timer = time(NULL)  */

      struct tm * now = localtime(&timer);

      printf("%6i/%i  %f   %2.2i:%2.2i:%2.2i  (%2.2i/%2.2i/%4i)   %8i sec.\n", istep, Nstep, time, now->tm_hour, now->tm_min, now->tm_sec, now->tm_mday, (now->tm_mon + 1), (now->tm_year + 1900), (int)(timer - unix_time_start));
      fflush(stdout);
    }
  }
}

void GRID::initRNG(gsl_rng* rng, unsigned long int auxiliary_seed){
  //INIZIALIZZO IL GENERATORE DI NUMERI CASUALI
  //Seeding del generatore di numeri casuali.
  //Strategia: con il generatore Mersenne Twister di GSL (inizializzato allo stesso modo per tutti)
  //calcolo per il processo con rank i, i numeri casuali tutti diversi.
  //Tengo soltanto l'ultimo numero, che uso per inizializzare un altro generatore
  //di numeri casuali di tipo differente (RANLUX).
  //Questo generatore verrà passato come argomento in add_momenta

  gsl_rng* rng_aux;

  unsigned long int* seeds;

  rng_aux = gsl_rng_alloc(gsl_rng_mt19937);

  gsl_rng_set(rng_aux, auxiliary_seed);
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  seeds = new unsigned long int[myrank + 1];

  for (int i = 0; i <= myrank; i++){
    seeds[i] = gsl_rng_get(rng_aux);
    //Questo controllo potrebbe essere superfluo...
    for (int j = 0; j < i; j++){
      if (seeds[j] == seeds[i]){ i--; break; }
    }

  }
  gsl_rng_set(rng, seeds[myrank]);

  delete[] seeds;
}

void GRID::visualDiag(){
  int Nstep = totalNumberOfTimesteps;
  MPI_Barrier(MPI_COMM_WORLD);

  const char* s_pbc = "PBC";
  const char* s_open = "OPEN";
  const char* s_pml = "PML";
  const char* s_err = "ERRORE !";

  const char* s_x, *s_y, *s_z;

  switch (xBoundaryConditions){
  case _PBC:
    s_x = s_pbc;
    break;

  case _Open:
    s_x = s_open;
    break;

  case _PML:
    s_x = s_pml;
    break;

  default:
    s_x = s_err;
  }

  switch (yBoundaryConditions){
  case _PBC:
    s_y = s_pbc;
    break;

  case _Open:
    s_y = s_open;
    break;

  case _PML:
    s_y = s_pml;
    break;

  default:
    s_y = s_err;
  }

  switch (zBoundaryConditions){
  case _PBC:
    s_z = s_pbc;
    break;

  case _Open:
    s_z = s_open;
    break;

  case _PML:
    s_z = s_pml;
    break;

  default:
    s_z = s_err;
  }

  if (myid == master_proc){
    GRID::printLogo();
    GRID::printProcInformations();
    GRID::printGridProcessorInformation();
    printf("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
    printf("dt=%f\n", dt);
    printf("dx=%f\n", dr[0]);
    printf("dy=%f\n", dr[1]);
    printf("dz=%f\n", dr[2]);
    printf("Nstep=%i\n", Nstep);
    printf("Boundaries (X,Y,Z) : %s , %s , %s \n", s_x, s_y, s_z);
    printf("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
  }
}
void GRID::setGridDeltar(){
  if (flagStretched){
    GRID::setGridDeltarStretched();
  }
  else{
    GRID::setGridDeltarNormal();
  }
}

void GRID::setGridDeltarNormal(){
  switch (accesso.dimensions)
  {
    //ALERT GRIGLIA
  case 3:
    dr[0] = (rmax[0] - rmin[0]) / (NGridNodes[0] - 1);
    dr[1] = (rmax[1] - rmin[1]) / (NGridNodes[1] - 1);
    dr[2] = (rmax[2] - rmin[2]) / (NGridNodes[2] - 1);

    dri[0] = 1. / dr[0];
    dri[1] = 1. / dr[1];
    dri[2] = 1. / dr[2];
    break;
  case 2:
    dr[0] = (rmax[0] - rmin[0]) / (NGridNodes[0] - 1);
    dr[1] = (rmax[1] - rmin[1]) / (NGridNodes[1] - 1);
    NGridNodes[2] = Nloc[2] = 1;
    dr[2] = (rmax[2] - rmin[2]);

    dri[0] = 1. / dr[0];
    dri[1] = 1. / dr[1];
    dri[2] = 1. / dr[2];
    break;
  case 1:
    dr[0] = (rmax[0] - rmin[0]) / (NGridNodes[0] - 1);
    NGridNodes[1] = NGridNodes[2] = Nloc[1] = Nloc[2] = 1;
    dr[1] = rmax[1] - rmin[1];
    dr[2] = rmax[2] - rmin[2];

    dri[0] = 1. / dr[0];
    dri[1] = 1. / dr[1];
    dri[2] = 1. / dr[2];
    break;
  default:
    printf("WRONG definition of DIMENSIONALITY\n");
    exit(17);
    break;
  }
}
void GRID::setGridDeltarStretched(){
  switch (accesso.dimensions)
  {
  case 3:
    dr[0] = (rmaxUniformGrid[0] - rminUniformGrid[0]) / (NUniformGrid[0] - 1);
    dr[1] = (rmaxUniformGrid[1] - rminUniformGrid[1]) / (NUniformGrid[1] - 1);
    dr[2] = (rmaxUniformGrid[2] - rminUniformGrid[2]) / (NUniformGrid[2] - 1);

    dri[0] = 1. / dr[0];
    dri[1] = 1. / dr[1];
    dri[2] = 1. / dr[2];
    break;
  case 2:
    dr[0] = (rmaxUniformGrid[0] - rminUniformGrid[0]) / (NUniformGrid[0] - 1);
    dr[1] = (rmaxUniformGrid[1] - rminUniformGrid[1]) / (NUniformGrid[1] - 1);
    NGridNodes[2] = Nloc[2] = 1;
    dr[2] = (rmax[2] - rmin[2]);

    dri[0] = 1. / dr[0];
    dri[1] = 1. / dr[1];
    dri[2] = 1. / dr[2];
    break;
  case 1:
    dr[0] = (rmaxUniformGrid[0] - rminUniformGrid[0]) / (NUniformGrid[0] - 1);
    NGridNodes[1] = NGridNodes[2] = Nloc[1] = Nloc[2] = 1;
    dr[1] = rmax[1] - rmin[1];
    dr[2] = rmax[2] - rmin[2];

    dri[0] = 1. / dr[0];
    dri[1] = 1. / dr[1];
    dri[2] = 1. / dr[2];
    break;
  default:
    printf("WRONG definition of DIMENSIONALITY\n");
    exit(17);
    break;
  }
}
void GRID::printLogo()
{
  std::cout << common_logo;
  std::cout << common_versionLine << std::endl;
  fflush(stdout);
}
void GRID::printProcInformations()
{
  printf("========== %19s ==========\n", "GRID_INITIALIZATION");
  printf("master_proc=%i:\n", myid);
  double dtotcell = ((long int)uniquePoints[0])*((long int)uniquePoints[1])*uniquePoints[2];
  printf("Ncells   =%g :     ( %5i, %5i, %5i)\n", dtotcell, uniquePoints[0], uniquePoints[1], uniquePoints[2]);
  printf("Nprocs   =%6i :     ( %5i, %5i, %5i)\n", nproc, rnproc[0], rnproc[1], rnproc[2]);
  printf("Xrange = [ %6g : %6g ]\n", rmin[0], rmax[0]);
  printf("Yrange = [ %6g : %6g ]\n", rmin[1], rmax[1]);
  printf("Zrange = [ %6g : %6g ]\n", rmin[2], rmax[2]);

#ifdef _FLAG_DEBUG
  printf("single proc id and 3D coordinates:\n");
  printf("%6s = (%4s,%4s,%4s)\n", "id", "idx", "idy", "idz");

  for (int id = 0; id < nproc; id++){
    int rid[3];
    MPI_Cart_coords(cart_comm, id, 3, rid);
    printf("%6i = (%4i,%4i,%4i)\n", id, rid[0], rid[1], rid[2]);
  }
#endif
}
void GRID::checkProcNumber(){
  rnproc[0] = nproc / (rnproc[1] * rnproc[2]);
  if (nproc != (rnproc[0] * rnproc[1] * rnproc[2]))
  {
    printf("BADLY defined number of processors nproc=%i  rnprocx=%i  rnprocy=%i  rnprocz=%i!!!\n", nproc, rnproc[0], rnproc[1], rnproc[2]);
    exit(18);
  }
  if (rnproc[0] < 1 || NGridNodes[0] < 1)
  {
    printf("BADLY defined rnprocx=%i  or  Nx=%i   !!!\n", rnproc[0], NGridNodes[0]);
    exit(18);
  }
  if (rnproc[1] < 1)
  {
    printf("BADLY defined rnprocy=%i  or  Ny=%i   !!!\n", rnproc[1], NGridNodes[1]);
    exit(18);
  }
  if (rnproc[2] < 1)
  {
    printf("BADLY defined rnprocz=%i  or  Nz=%i   !!!\n", rnproc[2], NGridNodes[2]);
    exit(18);
  }

}
void GRID::printGridProcessorInformation(){
#ifdef _FLAG_DEBUG
  printf("==========         grid         ==========\n");
  printf("\t%4s: %5s = [ %6s : %6s ]\n", "id", "Nloc", "rmin", "rmax");

  int c = 0;
  printf("X:  #proc=%i\tNx=%i\n", rnproc[c], NGridNodes[c]);
  for (int pp = 0; pp < rnproc[c]; pp++)
    printf("%16i: %5i = [ %6g : %6g ]\n", pp, rproc_Nloc[c][pp], rproc_rmin[c][pp], rproc_rmax[c][pp]);

  c = 1;
  printf("Y:  #proc=%i\tNy=%i\n", rnproc[c], NGridNodes[c]);
  for (int pp = 0; pp < rnproc[c]; pp++)
    printf("%16i: %5i = [ %6g : %6g ]\n", pp, rproc_Nloc[c][pp], rproc_rmin[c][pp], rproc_rmax[c][pp]);

  c = 2;
  printf("Z:  #proc=%i\tNz=%i\n", rnproc[c], NGridNodes[c]);
  for (int pp = 0; pp < rnproc[c]; pp++)
    printf("%16i: %5i = [ %6g : %6g ]\n", pp, rproc_Nloc[c][pp], rproc_rmin[c][pp], rproc_rmax[c][pp]);

#endif
  printf("========== %20s ==========\n", "");
  if (flagStretched){
    printf("Stretched GRID!!!\n");
    printf("\t%4s: %5s = [ %6s : %6s ]\n", "id", "Nloc", "Ximin", "Ximax");
    int c = 0;
    if (flagStretchedAlong[c]){
      printf("Stretched Grid along X enabled\n");
      printf("%20s = %i \n", "NUniformGrid", NUniformGrid[c]);
      printf("%20s = %i\n", "NLeftStretcheGrid", NLeftStretcheGrid[c]);
      printf("%20s = %i\n", "NRightStretcheGrid", NRightStretcheGrid[c]);
      printf("%20s = %i\n", "NGridNodes", NGridNodes[c]);
      printf("%20s = [ %g : %g ]\n", "UniformGrid", rminUniformGrid[c], rmaxUniformGrid[c]);
      printf("%20s = %i\n", "#proc", rnproc[c]);
      printf("%20s = [ %6g : %6g ]\n", "alpha", leftAlphaStretch[c], rightAlphaStretch[c]);
      //			for (int pp = 0; pp < rnproc[c]; pp++)
      //				printf("%16i: %5i = [ %6g : %6g ]\n", pp, rproc_Nloc[c][pp], rproc_csimin[c][pp], rproc_csimax[c][pp]);

    }
    c = 1;
    if (flagStretchedAlong[c]){
      printf("Stretched Grid along Y enabled\n");
      printf("%20s = %i \n", "NUniformGrid", NUniformGrid[c]);
      printf("%20s = %i\n", "NLeftStretcheGrid", NLeftStretcheGrid[c]);
      printf("%20s = %i\n", "NRightStretcheGrid", NRightStretcheGrid[c]);
      printf("%20s = %i\n", "NGridNodes", NGridNodes[c]);
      printf("%20s = [ %g : %g ]\n", "UniformGrid", rminUniformGrid[c], rmaxUniformGrid[c]);
      printf("%20s = %i\n", "#proc", rnproc[c]);
      printf("%20s = [ %6g : %6g ]\n", "alpha", leftAlphaStretch[c], rightAlphaStretch[c]);
      //			for (int pp = 0; pp < rnproc[c]; pp++)
      //				printf("%16i: %5i = [ %6g : %6g ]\n", pp, rproc_Nloc[c][pp], rproc_csimin[c][pp], rproc_csimax[c][pp]);

    }
    c = 2;
    if (flagStretchedAlong[c]){
      printf("Stretched Grid along Z enabled\n");
      printf("%20s = %i \n", "NUniformGrid", NUniformGrid[c]);
      printf("%20s = %i\n", "NLeftStretcheGrid", NLeftStretcheGrid[c]);
      printf("%20s = %i\n", "NRightStretcheGrid", NRightStretcheGrid[c]);
      printf("%20s = %i\n", "NGridNodes", NGridNodes[c]);
      printf("%20s = [ %g : %g ]\n", "UniformGrid", rminUniformGrid[c], rmaxUniformGrid[c]);
      printf("%20s = %i\n", "#proc", rnproc[c]);
      printf("%20s = [ %6g : %6g ]\n", "alpha", leftAlphaStretch[c], rightAlphaStretch[c]);
      //			for (int pp = 0; pp < rnproc[c]; pp++)
      //				printf("%16i: %5i = [ %6g : %6g ]\n", pp, rproc_Nloc[c][pp], rproc_csimin[c][pp], rproc_csimax[c][pp]);

    }

  }

}
void GRID::allocateRProcQuantities(){

  for (int c = 0; c < 3; c++){
    rproc_rmin[c] = (double*)malloc(rnproc[c] * sizeof(double));
    rproc_rmax[c] = (double*)malloc(rnproc[c] * sizeof(double));
    rproc_csimin[c] = (double*)malloc(rnproc[c] * sizeof(double));
    rproc_csimax[c] = (double*)malloc(rnproc[c] * sizeof(double));
    rproc_imin[c] = (int*)malloc(rnproc[c] * sizeof(int));
    rproc_imax[c] = (int*)malloc(rnproc[c] * sizeof(int));
    rproc_Nloc[c] = (int*)malloc(rnproc[c] * sizeof(int));
    rproc_NuniquePointsloc[c] = (int*)malloc(rnproc[c] * sizeof(int));
  }
}
void GRID::computeRProcNloc(){

  int ibuffer1[3], ibuffer2[3];

  for (int c = 0; c < 3; c++) {
    ibuffer1[c] = (NGridNodes[c] - 1) / rnproc[c];  //integer number of grid points (-1)
    ibuffer2[c] = (NGridNodes[c] - 1) % rnproc[c];  //if number of points is not multiple of
  }
  //============= cycling through all the processor in the "c" direction
  //============= if number of points (N[c]-1) is not multiple of rnproc[c]
  //============= then I consider the rest... ibuffer2[c] = (N[c]-1)- ibuffer1[c]*rnproc[c]

  for (int c = 0; c < 3; c++) {
    for (int pp = 0; pp < rnproc[c]; pp++){
      if (pp < ibuffer2[c]){ // Nloc= N/nproc + 1 + rest
        rproc_Nloc[c][pp] = ibuffer1[c] + 2;
        rproc_NuniquePointsloc[c][pp] = rproc_Nloc[c][pp] - 1;
        if (rproc_NuniquePointsloc[c][pp] == 0)
          rproc_NuniquePointsloc[c][pp] = 1;
      }
      else{ // Nloc= N/nproc + 1
        rproc_Nloc[c][pp] = ibuffer1[c] + 1;
        rproc_NuniquePointsloc[c][pp] = rproc_Nloc[c][pp] - 1;
        if (rproc_NuniquePointsloc[c][pp] == 0)
          rproc_NuniquePointsloc[c][pp] = 1;
      }
    }
  }
}
void GRID::computeRProcNuniquePointsLoc(){

  for (int c = 0; c < 3; c++) {
    for (int pp = 0; pp < rnproc[c]; pp++){
      if (rproc_NuniquePointsloc[c][pp] == 0)
        rproc_NuniquePointsloc[c][pp] = 1;
    }
  }
}
void GRID::setRminRmax(){
  for (int c = 0; c < 3; c++){

    // =========== definition of local minima maxima for integer indexes and real extrems
    rproc_rmin[c][0] = rmin[c];
    for (int pp = 1; pp < rnproc[c]; pp++){
      double prevRmin = rproc_rmin[c][pp - 1];
      int prevNloc = rproc_Nloc[c][pp - 1];
      double prevRmax;
      prevRmax = prevRmin + (prevNloc - 1)*dr[c];
      rproc_rmax[c][pp - 1] = prevRmax;
      rproc_rmin[c][pp] = rproc_rmax[c][pp - 1];
    }
    rproc_rmax[c][rnproc[c] - 1] = rmax[c];
  }
}

void GRID::setRminRmaxStretched(){
  GRID::setCsiminCsimax();
  for (int c = 0; c < 3; c++){

    // =========== definition of local minima maxima for integer indexes and real extrems
    for (int pp = 0; pp < rnproc[c]; pp++){
      rproc_rmin[c][pp] = stretchGrid(rproc_csimin[c][pp], c);
      rproc_rmax[c][pp] = stretchGrid(rproc_csimax[c][pp], c);
    }
  }
}

void GRID::setCsiminCsimax(){
  for (int c = 0; c < 3; c++){

    // =========== definition of local minima maxima for integer indexes and real extrems
    rproc_csimin[c][0] = csimin[c];
    for (int pp = 1; pp < rnproc[c]; pp++){
      double prevRmin = rproc_csimin[c][pp - 1];
      int prevNloc = rproc_Nloc[c][pp - 1];
      double prevRmax;
      prevRmax = prevRmin + (prevNloc - 1)*dr[c];
      rproc_csimax[c][pp - 1] = prevRmax;
      rproc_csimin[c][pp] = rproc_csimax[c][pp - 1];
    }
    rproc_csimax[c][rnproc[c] - 1] = csimax[c];
  }
}

void GRID::setIminImax(){
  for (int c = 0; c < 3; c++){

    rproc_imin[c][0] = 0;
    for (int pp = 1; pp < rnproc[c]; pp++){
      int prevImin = rproc_imin[c][pp - 1];
      int prevNloc = rproc_Nloc[c][pp - 1];
      int prevImax;
      prevImax = prevImin + prevNloc - 1;
      rproc_imax[c][pp - 1] = prevImax;
      rproc_imin[c][pp] = rproc_imax[c][pp - 1];
    }
    rproc_imax[c][rnproc[c] - 1] = (NGridNodes[c] - 1);
  }
}
void GRID::checkProcNumberInitialization(){
  for (int c = 0; c < 3; c++){
    int check = 0;
    for (int pp = 0; pp < rnproc[c]; pp++)
      check += (rproc_Nloc[c][pp] - 1);

    if ((check != (NGridNodes[c] - 1))){
      printf("ERROR BAD PROC INITIALIZATION at dimension number %i  !!!\n %i - %i \n", c, check, NGridNodes[c] - 1);
      exit(18);
    }
  }
}
void GRID::setLocalExtrems(){
  for (int c = 0; c < 3; c++){
    rminloc[c] = rproc_rmin[c][rmyid[c]];
    rmaxloc[c] = rproc_rmax[c][rmyid[c]];
    Nloc[c] = rproc_Nloc[c][rmyid[c]];
    csiminloc[c] = rproc_csimin[c][rmyid[c]];
    csimaxloc[c] = rproc_csimax[c][rmyid[c]];
    if (c >= accesso.dimensions)
      uniquePointsloc[c] = Nloc[c];
    else
      uniquePointsloc[c] = Nloc[c] - 1;
  }
}
void GRID::mpi_grid_initialize(int *narg, char **args)
{
  MPI_Init(narg, &args);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  GRID::checkStretchedGridInitialization();
  if (flagStretched){
    for (int c = 0; c < accesso.dimensions; c++){
      GRID::checkStretchedGridNpointsAlong(c);
      GRID::computeNStretchedPointsAlong(c);
    }
  }
  GRID::setGridDeltar();

  if (flagStretched){
    for (int c = 0; c < accesso.dimensions; c++){
      GRID::checkStretchedGridExtensionAlong(c);
      GRID::computeAlphaStretchAlong(c);
      GRID::computeExtremsStretchedGrid(c);
    }
  }

  GRID::checkProcNumber();
  MPI_Cart_create(MPI_COMM_WORLD, 3, rnproc, cyclic, _REORDER_MPI_CART_PROCESSES, &cart_comm);
  MPI_Cart_coords(cart_comm, myid, 3, rmyid);

  GRID::allocateRProcQuantities();
  GRID::computeRProcNloc();
  GRID::computeRProcNuniquePointsLoc();
  GRID::setIminImax();
  if (flagStretched)
    GRID::setRminRmaxStretched();
  else
    GRID::setRminRmax();

  GRID::checkProcNumberInitialization();
  GRID::setLocalExtrems();
  GRID::computeDerivativeCorrection();
  GRID::computeTotUniquePoints();
}

void GRID::finalize()
{
  int c;
  for (c = 0; c < accesso.dimensions; c++)
  {
    cir[c] = (double*)malloc(NGridNodes[c] * sizeof(double));
    chr[c] = (double*)malloc(NGridNodes[c] * sizeof(double));
    cirloc[c] = (double*)malloc(Nloc[c] * sizeof(double));
    chrloc[c] = (double*)malloc(Nloc[c] * sizeof(double));
    if (flagStretchedAlong[c]){
      for (int i = 0; i < NGridNodes[c]; i++){
        cir[c][i] = stretchGrid((csimin[c] + dr[c] * i), c);
        chr[c][i] = stretchGrid((csimin[c] + dr[c] * (i + 0.5)), c);
      }
      for (int i = 0; i < Nloc[c]; i++){
        cirloc[c][i] = stretchGrid((csiminloc[c] + dr[c] * i), c);
        chrloc[c][i] = stretchGrid((csiminloc[c] + dr[c] * (i + 0.5)), c);
      }
    }
    else{
      for (int i = 0; i < NGridNodes[c]; i++){
        cir[c][i] = rmin[c] + dr[c] * i;
        chr[c][i] = cir[c][i] + 0.5*dr[c];
      }
      for (int i = 0; i < Nloc[c]; i++){
        cirloc[c][i] = rminloc[c] + dr[c] * i;
        chrloc[c][i] = cirloc[c][i] + 0.5*dr[c];
      }
    }
  }
  for (; c < 3; c++)
  {
    cir[c] = (double*)malloc(1 * sizeof(double));
    chr[c] = (double*)malloc(1 * sizeof(double));
    cirloc[c] = (double*)malloc(1 * sizeof(double));
    chrloc[c] = (double*)malloc(1 * sizeof(double));

    cir[c][0] = 0.5*(rmin[c] + rmax[c]);
    chr[c][0] = 0.5*(rmin[c] + rmax[c]);
    cirloc[c][0] = 0.5*(rmin[c] + rmax[c]);
    chrloc[c][0] = 0.5*(rmin[c] + rmax[c]);
  }
}

void GRID::computeTotUniquePoints()
{
  proc_totUniquePoints = new int[nproc];

  for (int rank = 0; rank < nproc; rank++){
    int rid[3];
    MPI_Cart_coords(cart_comm, rank, 3, rid);
    proc_totUniquePoints[rank] = rproc_NuniquePointsloc[0][rid[0]];
    proc_totUniquePoints[rank] *= rproc_NuniquePointsloc[1][rid[1]];
    proc_totUniquePoints[rank] *= rproc_NuniquePointsloc[2][rid[2]];
  }
  int c;
  for (c = 0; c < accesso.dimensions; c++)
    uniquePoints[c] = NGridNodes[c] - 1;
  for (; c < 3; c++)
    uniquePoints[c] = NGridNodes[c];

}
/* ********* STRETCHED GRID ************* */
bool GRID::isStretched(){
  return flagStretched;
}

bool GRID::isStretchedAlong(int c){
  return flagStretchedAlong[c];
}
bool GRID::isLeftStretchedAlong(int c){
  return flagLeftStretchedAlong[c];
}
bool GRID::isRightStretchedAlong(int c){
  return flagRightStretchedAlong[c];
}

void GRID::checkStretchedGridInitialization(){
  if (!flagStretched)
    return;
  if (withMovingWindow && isStretchedAlong(0)){
    if (myid == 0){
      std::cout << "ERROR: moving window is incompatible with stretching along x!" << std::endl;
      std::cout.flush();

    }
    emergencyStop();
  }
}

void GRID::enableStretchedGrid(){
  flagStretched = true;
  for (int c = 0; c < 3; c++){
    rminUniformGrid[c] = rmin[c];
    rmaxUniformGrid[c] = rmax[c];
    csimin[c] = rmin[c];
    csimax[c] = rmax[c];
    NRightStretcheGrid[c] = NLeftStretcheGrid[c] = 0;
    NUniformGrid[c] = NGridNodes[c];
    leftAlphaStretch[c] = 1e200;
    rightAlphaStretch[c] = 1e200;
  }
}

void GRID::setXandNxLeftStretchedGrid(double min, int N){
  if (!flagStretched)
    return;
  flagStretchedAlong[0] = true;
  flagLeftStretchedAlong[0] = true;
  rminUniformGrid[0] = min;
  if ((N >= NGridNodes[0]) || (N < 0)) {
    printf("ERROR NxLeftStretcheGrid cannot be set to %i \n", N);
    exit(19);
  }
  NLeftStretcheGrid[0] = N;
}
void GRID::setYandNyLeftStretchedGrid(double min, int N){
  if (!flagStretched)
    return;
  flagStretchedAlong[1] = true;
  flagLeftStretchedAlong[1] = true;
  rminUniformGrid[1] = min;
  if ((N >= NGridNodes[1]) || (N < 0)) {
    printf("ERROR NyLeftStretcheGrid cannot be set to %i \n", N);
    exit(19);
  }
  NLeftStretcheGrid[1] = N;
}
void GRID::setZandNzLeftStretchedGrid(double min, int N){
  if (!flagStretched)
    return;
  flagStretchedAlong[2] = true;
  flagLeftStretchedAlong[2] = true;
  rminUniformGrid[2] = min;
  if ((N >= NGridNodes[2]) || (N < 0)) {
    printf("ERROR NzLeftStretcheGrid cannot be set to %i \n", N);
    exit(19);
  }
  NLeftStretcheGrid[2] = N;
}
void GRID::setXandNxRightStretchedGrid(double max, int N){
  if (!flagStretched)
    return;
  flagStretchedAlong[0] = true;
  flagRightStretchedAlong[0] = true;
  rmaxUniformGrid[0] = max;
  if ((N >= NGridNodes[0]) || (N < 0)) {
    printf("ERROR NxRightStretcheGrid cannot be set to %i \n", N);
    exit(19);
  }
  NRightStretcheGrid[0] = N;
}
void GRID::setYandNyRightStretchedGrid(double max, int N){
  if (!flagStretched)
    return;
  flagStretchedAlong[1] = true;
  flagRightStretchedAlong[1] = true;
  rmaxUniformGrid[1] = max;
  if ((N >= NGridNodes[1]) || (N < 0)) {
    printf("ERROR NyRightStretcheGrid cannot be set to %i \n", N);
    exit(19);
  }
  NRightStretcheGrid[1] = N;
}
void GRID::setZandNzRightStretchedGrid(double max, int N){
  if (!flagStretched)
    return;
  flagStretchedAlong[2] = true;
  flagRightStretchedAlong[2] = true;
  rmaxUniformGrid[2] = max;
  if ((N >= NGridNodes[2]) || (N < 0)) {
    printf("ERROR NzRightStretcheGrid cannot be set to %i \n", N);
    exit(19);
  }
  NRightStretcheGrid[2] = N;
}

//void GRID::setNxLeftStretchedGrid(int N){
//    if(!flagStretched)
//        return;
//    if((N>NGridNodes[0])||(N<0)) {
//        printf("ERROR NxLeftStretcheGrid cannot be set to %i \n", N);
//        exit(0);
//    }
//    NLeftStretcheGrid[0]=N;
//}
//void GRID::setNyLeftStretchedGrid(int N){
//    if(!flagStretched)
//        return;
//    if((N>NGridNodes[1])||(N<0)) {
//        printf("ERROR NyLeftStretcheGrid cannot be set to %i \n", N);
//        exit(0);
//    }
//    NLeftStretcheGrid[1]=N;

//}
//void GRID::setNzLeftStretchedGrid(int N){
//    if(!flagStretched)
//        return;
//    if((N>NGridNodes[2])||(N<0)) {
//        printf("ERROR NzLeftStretcheGrid cannot be set to %i \n", N);
//        exit(0);
//    }
//    NLeftStretcheGrid[2]=N;
//}


//void GRID::setNxUniformGrid(int N){
//    if(!flagStretched)
//        return;
//    NUniformGrid[0]=N+1;
//}
//void GRID::setNyUniformGrid(int N){
//    if(!flagStretched)
//        return;
//    NUniformGrid[1]=N+1;
//}
//void GRID::setNzUniformGrid(int N){
//    if(!flagStretched)
//        return;
//    NUniformGrid[2]=N+1;
//}

void GRID::checkStretchedGridNpointsAlong(int c){
  if (!flagStretched)
    return;
  if (!flagStretchedAlong[c])
    return;

  if ((flagLeftStretchedAlong[c]) || (flagRightStretchedAlong[c])) {
    if ((NRightStretcheGrid[c] + NLeftStretcheGrid[c]) >= (NGridNodes[c])) {
      printf("ERROR!!! stretched Grid along %c enabled\n", c);
      printf("\tNLeftStretcheGrid[%i]= %i", c, NLeftStretcheGrid[c]);
      printf("\tNRightStretcheGrid[%i]= %i", c, NRightStretcheGrid[c]);
      printf("\tand NGridNodes[%i]=%i", c, NGridNodes[c]);
      exit(13);
    }
  }

  fflush(stdout);
}
void GRID::computeNStretchedPointsAlong(int c){
  //if ((flagLeftStretchedAlong[c]) && (flagRightStretchedAlong[c])) {
  NUniformGrid[c] = NGridNodes[c] - NRightStretcheGrid[c] - NLeftStretcheGrid[c];
  //}
  //	else if ((flagLeftStretchedAlong[c]) && (!flagRightStretchedAlong[c])) {
  //		NRightStretcheGrid[c] = 0;
  //		NUniformGrid[c] = NGridNodes[c] - NLeftStretcheGrid[c];
  //	}
  //	else if ((!flagLeftStretchedAlong[c]) && (flagRightStretchedAlong[c])) {
  //		NLeftStretcheGrid[c] = 0;
  //		NUniformGrid[c] = NGridNodes[c] - NRightStretcheGrid[c];
  //	}
  //	else{
  //		NRightStretcheGrid[c] = NLeftStretcheGrid[c] = 0;
  //		flagStretchedAlong[c] = false;
  //		NUniformGrid[c] = NGridNodes[c];
  //	}
  //    printf("Stretched Grid along %i enabled\n",c);
  //    printf("\tNUniformGrid[%i]=%i \n",c,NUniformGrid[c]);
  //    printf("\tNLeftStretcheGrid[%i]= %i\n", c,NLeftStretcheGrid[c]);
  //    printf("\tNRightStretcheGrid[%i]= %i\n", c,NRightStretcheGrid[c]);
  //    printf("\tand NGridNodes[%i]=%i\n", c, NGridNodes[c]);
  //    printf("\tUniformGrid[ %g : %g ]\n" ,rminUniformGrid[c],rminUniformGrid[c]);
}


void GRID::checkStretchedGridExtensionAlong(int c){
  if (!flagStretched)
    return;
  if (!flagStretchedAlong[c])
    return;
  //    printf("Stretched Grid along %i enabled\n",c);
  //    printf("\tNUniformGrid[%i]=%i \n",c,NUniformGrid[c]);
  if (flagLeftStretchedAlong[c]){
    double Dx = rminUniformGrid[c] - rmin[c];
    double Dxi = NLeftStretcheGrid[c] * dr[c];

    //        printf("LEFT stretched Grid along %i enabled\n",c);
    //        printf("\tleft layer physical    thickness = %g\n", Dx);
    //        printf("\tleft layer unstretched thickness = %g\n", Dxi);
    //        printf("\tNLeftStretcheGrid[%i]=%i \n",c,NLeftStretcheGrid[c]);

    if (Dxi >= Dx){
      printf("ERROR!!! LEFT stretched Grid along %i enabled\n", c);
      printf("\tleft layer physical    thickness = %g\n", Dx);
      printf("\tleft layer unstretched thickness = %g\n", Dxi);
      exit(14);
    }
  }
  if (flagRightStretchedAlong[c]) {
    double Dx = rmax[c] - rmaxUniformGrid[c];
    double Dxi = NRightStretcheGrid[c] * dr[c];

    //        printf("RIGHT stretched Grid along %i enabled\n",c);
    //        printf("\tright layer physical    thickness = %g\n", Dx);
    //        printf("\tright layer unstretched thickness = %g\n", Dxi);
    //        printf("\tNRightStretcheGrid[%i]=%i \n",c,NRightStretcheGrid[c]);
    if (Dxi >= Dx){
      printf("ERROR!!! RIGHT stretched Grid along %i enabled\n", c);
      printf("\tright layer physical    thickness = %g\n", Dx);
      printf("\tright layer unstretched thickness = %g\n", Dxi);
      exit(14);
    }
  }


}

double myHelpFunction(double Dx, double Dxi, double alpha) {
  return Dx - alpha*tan(Dxi / alpha);
}
double computeAlpha(double Dx, double Dxi){
  double alphaMin = 1.01*2.0*Dxi / M_PI; //1% more than minimum
  double alphaMax = 20 * alphaMin;
  double difference = 2 * Dx;
  double guess = 0.0;
  double result = 0.0;

  while (fabs(difference) > (Dx*0.00001)){
    guess = 0.5*(alphaMin + alphaMax);
    result = myHelpFunction(Dx, Dxi, guess);
    if (result < 0)
      alphaMin = guess;
    else
      alphaMax = guess;
    difference = result;
  }
  return guess;
}

void GRID::computeAlphaStretchLeft(int c){
  double Dx = rminUniformGrid[c] - rmin[c];
  double Dxi = NLeftStretcheGrid[c] * dr[c];
  leftAlphaStretch[c] = computeAlpha(Dx, Dxi);
}
void GRID::computeAlphaStretchRight(int c){
  double Dx = rmax[c] - rmaxUniformGrid[c];
  double Dxi = NRightStretcheGrid[c] * dr[c];
  rightAlphaStretch[c] = computeAlpha(Dx, Dxi);
}

void GRID::computeAlphaStretchAlong(int c){
  if (flagLeftStretchedAlong[c]) {
    GRID::computeAlphaStretchLeft(c);
  }
  if (flagRightStretchedAlong[c]){
    GRID::computeAlphaStretchRight(c);
  }
  //printf("alphastretch along %i are L:%g   R:%g\n", c, leftAlphaStretch[c],rightAlphaStretch[c]);
}

void GRID::computeExtremsStretchedGrid(int c){
  if (!flagStretchedAlong[c])
    return;

  //csimin[c] = unStretchGrid(rmin[c], c);
  csimin[c] = rminUniformGrid[c] - dr[c] * NLeftStretcheGrid[c];
  csimax[c] = csimin[c] + dr[c] * (NGridNodes[c] - 1);
  rmax[c] = stretchGrid(csimax[c], c);
  rmin[c] = stretchGrid(csimin[c], c);
  //    printf("\tdr= %g NGridNodes=%i \n",dr[c],NGridNodes[c]);
  //    printf("\tcsi=[ %g:%g] \n",csimin[c],csimax[c]);
  //    printf("\tr=[ %g:%g] \n",rmin[c],rmax[c]);



}


double GRID::stretchingFunction(double csi, int c)
{
  // x = f(csi)
  double alphaL = leftAlphaStretch[c];
  double alphaR = rightAlphaStretch[c];
  double x0neg = rminUniformGrid[c];
  double x0pos = rmaxUniformGrid[c];
  if (csi < x0neg)
    return x0neg + alphaL*tan((csi - x0neg) / alphaL);
  else if (csi < x0pos)
    return csi;
  else
    return x0pos + alphaR*tan((csi - x0pos) / alphaR);

}

double GRID::inverseStretchingFunction(double x, int c)
{
  // --> csi=f^{-1}(x)
  double alphaL = leftAlphaStretch[c];
  double alphaR = rightAlphaStretch[c];
  double x0neg = rminUniformGrid[c];
  double x0pos = rmaxUniformGrid[c];
  if (x < x0neg)
    return x0neg + alphaL*atan((x - x0neg) / alphaL);
  else if (x < x0pos)
    return x;
  else
    return x0pos + alphaR*atan((x - x0pos) / alphaR);
}

double GRID::derivativeStretchingFunction(double csi, int c)
{
  if (!flagStretched)
    return 1;
  if (!flagStretchedAlong[c])
    return 1;
  // --> 1/f'(csi)
  double alphaL = leftAlphaStretch[c];
  double alphaR = rightAlphaStretch[c];
  double x0neg = rminUniformGrid[c];
  double x0pos = rmaxUniformGrid[c];
  double buf;

  if (csi < x0neg){
    buf = cos((csi - x0neg) / alphaL);
    return 1.0 / (buf*buf);
  }
  else if (csi < x0pos)
    return 1.0;
  else{
    buf = cos((csi - x0pos) / alphaR);
    return 1.0 / (buf*buf);
  }
}

double GRID::stretchGrid(double xi_x, int c)
{
  if (!flagStretched)
    return xi_x;
  // from rescaled (xi_x) to physical (x=x0+f(xi_x)) coordinates

  if (flagStretchedAlong[c])
    return stretchingFunction(xi_x, c);
  else
    return xi_x;
}

double GRID::unStretchGrid(double x, int c)
{
  if (!flagStretched)
    return x;
  // from physical (x) to rescaled (xi_x=f^{-1}(x-x0)) coordinates

  if (flagStretchedAlong[c])
    return inverseStretchingFunction(x, c);
  else
    return x;
}


void GRID::computeDerivativeCorrection()
{
  double csi;

  for (int c = 0; c < 3; c++)
  {
    iStretchingDerivativeCorrection[c] = (double*)malloc((Nloc[c])*sizeof(double));
    hStretchingDerivativeCorrection[c] = (double*)malloc((Nloc[c])*sizeof(double));

    for (int i = 0; i < Nloc[c]; i++)
    {
      if (flagStretchedAlong[c]){
        csi = csiminloc[c] + i*dr[c];
        iStretchingDerivativeCorrection[c][i] = 1. / derivativeStretchingFunction(csi, c);
        hStretchingDerivativeCorrection[c][i] = 1. / derivativeStretchingFunction(csi + dr[c] * 0.5, c);
      }
      else{
        iStretchingDerivativeCorrection[c][i] = 1.0;
        hStretchingDerivativeCorrection[c][i] = 1.0;
      }
    }
  }

}

bool GRID::checkAssignBoundary(axisBoundaryConditions cond, axisBoundaryConditions* axisCond){
  bool isOk = false;
  if (*axisCond == _notAssigned){
    *axisCond = cond;
    isOk = true;
  }
  else{
    *axisCond = cond;
    isOk = false;
  }
  return isOk;
}

void GRID::setBoundaries(int flags){
  xBoundaryConditions = yBoundaryConditions = zBoundaryConditions = _notAssigned;
  bool global_isOk = true;

  if (flags & xPBC)
    global_isOk = global_isOk && checkAssignBoundary(_PBC, &xBoundaryConditions);
  if (flags & yPBC)
    global_isOk = global_isOk && checkAssignBoundary(_PBC, &yBoundaryConditions);
  if (flags & zPBC)
    global_isOk = global_isOk && checkAssignBoundary(_PBC, &zBoundaryConditions);

  if (flags & xOpen)
    global_isOk = global_isOk && checkAssignBoundary(_Open, &xBoundaryConditions);
  if (flags & yOpen)
    global_isOk = global_isOk && checkAssignBoundary(_Open, &yBoundaryConditions);
  if (flags & zOpen)
    global_isOk = global_isOk && checkAssignBoundary(_Open, &zBoundaryConditions);

  if (flags & xPML)
    global_isOk = global_isOk && checkAssignBoundary(_PML, &xBoundaryConditions);
  if (flags & yPML)
    global_isOk = global_isOk && checkAssignBoundary(_PML, &yBoundaryConditions);
  if (flags & zPML)
    global_isOk = global_isOk && checkAssignBoundary(_PML, &zBoundaryConditions);

  if (xBoundaryConditions == _notAssigned)
    xBoundaryConditions = _PBC;
  if (yBoundaryConditions == _notAssigned)
    yBoundaryConditions = _PBC;
  if (zBoundaryConditions == _notAssigned)
    zBoundaryConditions = _PBC;

  if (accesso.dimensions == 1){
    yBoundaryConditions = _PBC;
    zBoundaryConditions = _PBC;
  }
  else if (accesso.dimensions == 2){
    zBoundaryConditions = _PBC;
  }


  if (xBoundaryConditions != _PBC)
    cyclic[0] = 0;

  if (yBoundaryConditions != _PBC)
    cyclic[1] = 0;

  if (zBoundaryConditions != _PBC)
    cyclic[2] = 0;


  if (!global_isOk && myid == master_proc){
    std::cout << "Warning! Incompatible boundary conditions! PML overrides Open and PBC, Open overrides PBC." << std::endl;
  }

}

axisBoundaryConditions GRID::getXBoundaryConditions(){
  return xBoundaryConditions;
}

axisBoundaryConditions GRID::getYBoundaryConditions(){
  return yBoundaryConditions;
}

axisBoundaryConditions GRID::getZBoundaryConditions(){
  return zBoundaryConditions;
}

void GRID::emergencyStop(){
  MPI_Finalize();
  exit(17);
}

void GRID::dump(std::ofstream &ff){
  ff.write((char*)&istep, sizeof(int));
  ff.write((char*)&time, sizeof(double));
  ff.write((char*)&mark_mw, sizeof(double));

  if (!withMovingWindow)
    return;
  ff.write((char*)&rmin[0], sizeof(double));
  ff.write((char*)&rmax[0], sizeof(double));
  ff.write((char*)&rminloc[0], sizeof(double));
  ff.write((char*)&rmaxloc[0], sizeof(double));

  ff.write((char*)&csimax[0], sizeof(double));
  ff.write((char*)&csimin[0], sizeof(double));
  ff.write((char*)&csiminloc[0], sizeof(double));
  ff.write((char*)&csimaxloc[0], sizeof(double));

  for (int pp = 0; pp < rnproc[0]; pp++){
    ff.write((char*)&rproc_rmin[0][pp], sizeof(double));
    ff.write((char*)&rproc_rmax[0][pp], sizeof(double));
  }
  for (int i = 0; i < NGridNodes[0]; i++){
    ff.write((char*)&cir[0][i], sizeof(double));
    ff.write((char*)&chr[0][i], sizeof(double));
  }
  for (int i = 0; i < Nloc[0]; i++){
    ff.write((char*)&cirloc[0][i], sizeof(double));
    ff.write((char*)&chrloc[0][i], sizeof(double));
  }
}

void GRID::reloadDump(std::ifstream &ff){
  ff.read((char*)&istep, sizeof(int));
  ff.read((char*)&time, sizeof(double));
  ff.read((char*)&mark_mw, sizeof(double));

  if (!withMovingWindow)
    return;
  if (1){
    ff.read((char*)&rmin[0], sizeof(double));
    ff.read((char*)&rmax[0], sizeof(double));
    ff.read((char*)&rminloc[0], sizeof(double));
    ff.read((char*)&rmaxloc[0], sizeof(double));

    ff.read((char*)&csimax[0], sizeof(double));
    ff.read((char*)&csimin[0], sizeof(double));
    ff.read((char*)&csiminloc[0], sizeof(double));
    ff.read((char*)&csimaxloc[0], sizeof(double));

    for (int pp = 0; pp < rnproc[0]; pp++){
      ff.read((char*)&rproc_rmin[0][pp], sizeof(double));
      ff.read((char*)&rproc_rmax[0][pp], sizeof(double));
    }
    for (int i = 0; i < NGridNodes[0]; i++){
      ff.read((char*)&cir[0][i], sizeof(double));
      ff.read((char*)&chr[0][i], sizeof(double));
    }
    for (int i = 0; i < Nloc[0]; i++){
      ff.read((char*)&cirloc[0][i], sizeof(double));
      ff.read((char*)&chrloc[0][i], sizeof(double));
    }
  }
  else{

    rmin[0] += mark_mw;
    rmax[0] += mark_mw;
    rminloc[0] += mark_mw;
    rmaxloc[0] += mark_mw;

    csimax[0] += mark_mw;
    csimin[0] += mark_mw;
    csiminloc[0] += mark_mw;
    csimaxloc[0] += mark_mw;

    for (int pp = 0; pp < rnproc[0]; pp++){
      rproc_rmin[0][pp] += mark_mw;
      rproc_rmax[0][pp] += mark_mw;
    }
    for (int i = 0; i < NGridNodes[0]; i++){
      cir[0][i] += mark_mw;
      chr[0][i] += mark_mw;
    }
    for (int i = 0; i < Nloc[0]; i++){
      cirloc[0][i] += mark_mw;
      chrloc[0][i] += mark_mw;
    }
  }
}
double GRID::getMarkMW(){
  return mark_mw;
}

void GRID::setDumpPath(std::string _dumpDir){
#if defined (USE_BOOST)
  if (myid == master_proc){
    if (!boost::filesystem::exists(_dumpDir)){
      boost::filesystem::create_directories(_dumpDir);
    }
  }
#endif
  dumpPath = _dumpDir;
}
std::string GRID::composeDumpFileName(int dumpID){
  std::stringstream dumpName;
  dumpName << dumpPath << "/DUMP_";
  dumpName << std::setw(2) << std::setfill('0') << std::fixed << dumpID << "_";
  dumpName << std::setw(5) << std::setfill('0') << std::fixed << myid << ".bin";
  return dumpName.str();
}
