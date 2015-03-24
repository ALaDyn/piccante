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

#include "particle_species.h"
//#define OLD_ACCESS


SPECIE::SPECIE()
{
  Ncomp = 7;
  allocated = false;
  Z = A = 0;
  isTestSpecies = false;
  spectrum.values = NULL;
  energyExtremesFlag = false;
  lastParticle = 0;
  flagWithMarker = false;
}
SPECIE::SPECIE(GRID *grid)
{
  Ncomp = 7;
  allocated = false;
  Z = A = 0;
  mygrid = grid;
  isTestSpecies = false;
  spectrum.values = NULL;
  energyExtremesFlag = false;
  lastParticle = 0;
  flagWithMarker = false;
}
void SPECIE::allocate_species()
{
  if (mygrid->withParticles == NO)
    return;
#ifndef NO_ALLOCATION
#ifdef _ACC_SINGLE_POINTER
  pData = (double*)malloc((Np*Ncomp)*sizeof(double));
#else
  val = (double**)malloc(Ncomp*sizeof(double*));
  for (int c = 0; c < Ncomp; c++){
    val[c] = (double*)malloc(Np*sizeof(double));
  }
#endif

  valSize = Np;
  allocated = true;
#endif

}
SPECIE::~SPECIE(){
#ifdef _ACC_SINGLE_POINTER
  free(pData);
#else
  for (int c = 0; c < Ncomp; c++){
    free(val[c]);
  }
  free(val);
#endif
}

void SPECIE::erase()
{
  if (mygrid->withParticles == NO)
    return;

  if (!allocated){
    printf("ERROR: species not allocated!!!\n");
    exit(11);
  }
#ifdef _ACC_SINGLE_POINTER
  memset((void*)pData, 0, (Np*Ncomp)*sizeof(double));
#else
  for (int c = 0; c < Ncomp; c++){
    memset((void*)val[c], 0, Np*sizeof(double));
  }
#endif
}
void SPECIE::reallocate_species()
{
  if (mygrid->withParticles == NO)
    return;

  if (!allocated)
  {
#ifndef NO_ALLOCATION
    printf("\nERROR: species not allocated\n\n");
    exit(11);
#else
    return;
#endif
  }
  if (Np > valSize){
    valSize = Np + allocsize;
#ifdef _ACC_SINGLE_POINTER
    pData = (double *)realloc((void*)pData, valSize*Ncomp*sizeof(double));
#else
    for (int c = 0; c < Ncomp; c++){
      val[c] = (double *)realloc((void*)val[c], valSize*sizeof(double));
    }
#endif
  }
  else if (Np < (valSize - allocsize)){
    valSize = Np + allocsize;
#ifdef _ACC_SINGLE_POINTER
    pData = (double *)realloc((void*)pData, valSize*Ncomp*sizeof(double));
#else
    for (int c = 0; c < Ncomp; c++){
      val[c] = (double *)realloc((void*)val[c], valSize*sizeof(double));
    }
#endif
  }
  return;

}

SPECIE SPECIE::operator = (SPECIE &destro)
{
  if (!destro.allocated){ printf("---ERROR---\noperation not permitted\nSPECIE=SPECIE\nnot allocated\n"); exit(11); }

  Np = destro.Np;
  Ncomp = destro.Ncomp;
  type = destro.type;
  coupling = destro.coupling;
  Numerical2Physical_particles = destro.Numerical2Physical_particles;
  mygrid = destro.mygrid;
  chargeSign = destro.chargeSign;
  mass = destro.mass;
  plasma = destro.plasma;
  isTestSpecies = destro.isTestSpecies;
  for (int i = 0; i < 3; i++)
  {
    particlePerCellXYZ[i] = destro.particlePerCellXYZ[i];
    minima[i] = destro.minima[i];
    minima[3 + i] = destro.minima[3 + i];
    maxima[i] = destro.maxima[i];
    maxima[3 + i] = destro.maxima[3 + i];
  }
  if (!allocated)
  {
    allocate_species();
  }
  else reallocate_species();
#ifdef _ACC_SINGLE_POINTER
  memcpy((void*)pData, (void*)destro.pData, Np*Ncomp*sizeof(double));
#else
  for (int c = 0; c < Ncomp; c++){
    memcpy((void*)val[c], (void*)destro.val[c], Np*sizeof(double));
  }
#endif
  return *this;
}

//CREATION: create particles at the beginning of the simulation, evaluates the number of particles needed for the simulation
// allocate a slightly bigger array, effectively initialize the particles with SEVEN 7 doubles, x,y,z,ux,uy,uz,weight
// summing the weights of the particles in one celle gives ONE*the_electron_density
// if the plasma density is 0.1 times the reference density (i.e. the critical density) and 40 particles per cell are used
// each particle has w=0.0025

std::string SPECIE::getName(){
  return name;
}

void SPECIE::addMarker(){
  flagWithMarker = true;
}
void SPECIE::setTestSpecies(){
  isTestSpecies = true;
}

bool SPECIE::amIWithMarker(){
  return flagWithMarker;
}

void SPECIE::computeParticleMassChargeCoupling(){
  if (type == ELECTRON){
    coupling = -1.;
    mass = 1.0;
    Z = -1.0;
    chargeSign = -1.0;
  }
  if (type == POSITRON){
    coupling = 1.;
    mass = 1.0;
    Z = 1.0;
    chargeSign = 1.0;
  }
  if (type == ION){
    if (Z == 0 || A == 0)
    {
      printf("ERROR: Ion charge or mass NOT defined!\n");
      exit(11);
    }
    else{
      coupling = Z / (1.8362e3*A);
      mass = 1.8362e3*A;
    }
    chargeSign = 1.0;
  }
}
int SPECIE::getNumberOfParticlesWithin(double plasmarmin[3], double plasmarmax[3]){

  int counter = 0;
  double xloc, yloc, zloc;
  int Nx = mygrid->Nloc[0];
  int Ny = mygrid->Nloc[1];
  int Nz = mygrid->Nloc[2];

  for (int k = 0; k < Nz; k++)
    for (int j = 0; j < Ny; j++)
      for (int i = 0; i < Nx; i++)
      {
        xloc = mygrid->chrloc[0][i];
        yloc = mygrid->chrloc[1][j];
        zloc = mygrid->chrloc[2][k];

        if (xloc >= plasmarmin[0] && xloc <= plasmarmax[0])
          if (yloc >= plasmarmin[1] && yloc <= plasmarmax[1])
            if (zloc >= plasmarmin[2] && zloc <= plasmarmax[2]){
              if (plasma.density_function(xloc, yloc, zloc, plasma.params, Z, A) > 0)
                counter += particlePerCell;
            }
      }
  return counter;
}
int SPECIE::getNumberOfParticlesWithinFromFile1D(double plasmarmin[], double plasmarmax[], std::string name){

  std::ifstream fileDensity(name.c_str(), std::ifstream::in);
  int Nx_in;
  double xmin, xmax, dx_in;
  double temp;
  fileDensity >> Nx_in;
  fileDensity >> temp;
  fileDensity >> xmin;
  fileDensity >> xmax;
  printf("Benvenuto! leggero' il file \"%s\"\n", name.c_str());
  fflush(stdout);

  double *density, *velocity, *x_in, *uy, *uz;
  density = (double*)malloc(sizeof(double)*Nx_in);
  velocity = (double*)malloc(sizeof(double)*Nx_in);
  x_in = (double*)malloc(sizeof(double)*Nx_in);
  uy = (double*)malloc(sizeof(double)*Nx_in);
  uz = (double*)malloc(sizeof(double)*Nx_in);
  for (int i = 0; i < Nx_in; i++){
    fileDensity >> x_in[i];
    fileDensity >> density[i];
    fileDensity >> velocity[i];
    fileDensity >> uy[i];
    fileDensity >> uz[i];

  }
  fileDensity.close();
  xmin = x_in[0];
  xmax = x_in[Nx_in - 1];
  dx_in = x_in[1] - x_in[0];


  int counter = 0;
  double xloc, yloc, zloc;
  int Nx = mygrid->Nloc[0];
  int Ny = mygrid->Nloc[1];
  int Nz = mygrid->Nloc[2];

  double axh, wh[2], value;
  int ih, ihleft, ihright;

  for (int k = 0; k < Nz; k++)
    for (int j = 0; j < Ny; j++)
      for (int i = 0; i < Nx; i++)
      {
        xloc = mygrid->chrloc[0][i];
        yloc = mygrid->chrloc[1][j];
        zloc = mygrid->chrloc[2][k];
        ih = (int)((xloc - xmin) / dx_in);
        axh = (xloc - xmin) / dx_in - ih;
        wh[0] = 1 - axh;
        wh[1] = axh;
        ihleft = (ih + Nx_in - 1) % (Nx_in - 1);
        ihright = (ih + 1) % (Nx_in - 1);

        if (xloc >= plasmarmin[0] && xloc <= plasmarmax[0])
          if (yloc >= plasmarmin[1] && yloc <= plasmarmax[1])
            if (zloc >= plasmarmin[2] && zloc <= plasmarmax[2]){
              value = wh[0] * density[ihleft] + wh[1] * density[ihright];

              if (value > 0)
                counter += particlePerCell;
            }
      }

  free(density);
  free(velocity);
  free(x_in);
  free(uy);
  free(uz);
  return counter;
}
void SPECIE::createParticlesWithinFrom(double plasmarmin[3], double plasmarmax[3], int oldNumberOfParticles, uint64_t disp){
#ifdef NO_ALLOCATION
  if(!allocated){
    return;
  }
#endif
  int counter = oldNumberOfParticles;
  double xloc, yloc, zloc;
  int Nx = mygrid->Nloc[0];
  int Ny = mygrid->Nloc[1];
  int Nz = mygrid->Nloc[2];
  double dx = mygrid->dr[0];
  double dy = mygrid->dr[1];
  double dz = mygrid->dr[2];
  double dxp = dx / particlePerCellXYZ[0];
  double dyp = dy / particlePerCellXYZ[1];
  double dzp = dz / particlePerCellXYZ[2];
  double  weight;

  for (int k = 0; k < Nz; k++)
    for (int j = 0; j < Ny; j++)
      for (int i = 0; i < Nx; i++)
      {
        xloc = mygrid->chrloc[0][i];
        yloc = mygrid->chrloc[1][j];
        zloc = mygrid->chrloc[2][k];

        if (xloc >= plasmarmin[0] && xloc <= plasmarmax[0])
          if (yloc >= plasmarmin[1] && yloc <= plasmarmax[1])
            if (zloc >= plasmarmin[2] && zloc <= plasmarmax[2])
            {
              if (plasma.density_function(xloc, yloc, zloc, plasma.params, Z, A) > 0)
              {
                weight = plasma.density_function(xloc, yloc, zloc, plasma.params, Z, A) / particlePerCell;

                xloc -= 0.5*dx;
                yloc -= 0.5*dy;
                zloc -= 0.5*dz;
                for (int ip = 0; ip < particlePerCellXYZ[0]; ip++)
                  for (int jp = 0; jp < particlePerCellXYZ[1]; jp++)
                    for (int kp = 0; kp < particlePerCellXYZ[2]; kp++)
                    {
                      r0(counter) = xloc + dxp*(ip + 0.5);
                      r1(counter) = yloc + dyp*(jp + 0.5);
                      r2(counter) = zloc + dzp*(kp + 0.5);
                      u0(counter) = u1(counter) = u2(counter) = 0;
                      w(counter) = weight;
                      if(flagWithMarker)
                        marker(counter) = (counter + disp);
                      if (isTestSpecies)
                        w(counter) = (double)(counter + disp);
                      counter++;
                    }
              }
            }
      }
}

void SPECIE::createStretchedParticlesWithinFrom(double plasmarmin[3], double plasmarmax[3], int oldNumberOfParticles, uint64_t disp){
#ifdef NO_ALLOCATION
  if(!allocated){
    return;
  }
#endif
  int counter = oldNumberOfParticles;
  double xloc, yloc, zloc;
  double myx, myy, myz;
  double mydx, mydy, mydz;
  double mycsix, mycsiy, mycsiz;
  double csilocx, csilocy, csilocz;
  int Nx = mygrid->Nloc[0];
  int Ny = mygrid->Nloc[1];
  int Nz = mygrid->Nloc[2];
  double dx = mygrid->dr[0];
  double dy = mygrid->dr[1];
  double dz = mygrid->dr[2];
  double dxp = dx / particlePerCellXYZ[0];
  double dyp = dy / particlePerCellXYZ[1];
  double dzp = dz / particlePerCellXYZ[2];
  double  weight;

  for (int k = 0; k < Nz; k++){
    zloc = mygrid->chrloc[2][k];
    csilocz = mygrid->csiminloc[2] + dz*k;
    for (int j = 0; j < Ny; j++){
      yloc = mygrid->chrloc[1][j];
      csilocy = mygrid->csiminloc[1] + dy*j;
      for (int i = 0; i < Nx; i++){
        xloc = mygrid->chrloc[0][i];
        csilocx = mygrid->csiminloc[0] + dx*i;

        if (xloc >= plasmarmin[0] && xloc <= plasmarmax[0]){
          if (yloc >= plasmarmin[1] && yloc <= plasmarmax[1]){
            if (zloc >= plasmarmin[2] && zloc <= plasmarmax[2])
            {
              if (plasma.density_function(xloc, yloc, zloc, plasma.params, Z, A) > 0)
              {
                weight = plasma.density_function(xloc, yloc, zloc, plasma.params, Z, A) / particlePerCell;
                for (int kp = 0; kp < particlePerCellXYZ[2]; kp++){
                  mycsiz = csilocz + dzp*(kp + 0.5);
                  myz = mygrid->stretchGrid(mycsiz, 2);
                  mydz = mygrid->derivativeStretchingFunction(mycsiz, 2);

                  for (int jp = 0; jp < particlePerCellXYZ[1]; jp++){
                    mycsiy = csilocy + dyp*(jp + 0.5);
                    myy = mygrid->stretchGrid(mycsiy, 1);
                    mydy = mygrid->derivativeStretchingFunction(mycsiy, 1);

                    for (int ip = 0; ip < particlePerCellXYZ[0]; ip++){
                      mycsix = csilocx + dxp*(ip + 0.5);
                      myx = mygrid->stretchGrid(mycsix, 0);
                      mydx = mygrid->derivativeStretchingFunction(mycsix, 0);

                      r0(counter) = myx;
                      r1(counter) = myy;
                      r2(counter) = myz;
                      u0(counter) = u1(counter) = u2(counter) = 0;
                      w(counter) = weight*mydx*mydy*mydz;
                      if(flagWithMarker)
                        marker(counter) = (counter + disp);
                      if (isTestSpecies)
                        w(counter) = (double)(counter + disp);
                      counter++;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

void SPECIE::createParticlesWithinFromButFromFile1D(double plasmarmin[3], double plasmarmax[3], int oldNumberOfParticles, uint64_t disp, std::string name){
  int counter = oldNumberOfParticles;
  double xloc, yloc, zloc;
  int Nx = mygrid->Nloc[0];
  int Ny = mygrid->Nloc[1];
  int Nz = mygrid->Nloc[2];
  double dx = mygrid->dr[0];
  double dy = mygrid->dr[1];
  double dz = mygrid->dr[2];
  double dxp = dx / particlePerCellXYZ[0];
  double dyp = dy / particlePerCellXYZ[1];
  double dzp = dz / particlePerCellXYZ[2];
  double  weight;
  double myuy, myuz;

  std::ifstream fileDensity(name.c_str(), std::ifstream::in);
  int Nx_in;
  double xmin, xmax, dx_in;
  double temp;
  fileDensity >> Nx_in;
  fileDensity >> temp;
  fileDensity >> xmin >> xmax;

  double *density, *velocity, *x_in, *uy, *uz;
  density = (double*)malloc(sizeof(double)*Nx_in);
  velocity = (double*)malloc(sizeof(double)*Nx_in);
  x_in = (double*)malloc(sizeof(double)*Nx_in);
  uy = (double*)malloc(sizeof(double)*Nx_in);
  uz = (double*)malloc(sizeof(double)*Nx_in);
  for (int i = 0; i < Nx_in; i++){
    fileDensity >> x_in[i];
    fileDensity >> density[i];
    fileDensity >> velocity[i];
    fileDensity >> uy[i];
    fileDensity >> uz[i];
  }

  xmin = x_in[0];
  xmax = x_in[Nx_in - 1];
  dx_in = x_in[1] - x_in[0];

  double axh, wh[2], denValue, velValue;
  int ih, ihleft, ihright;

  for (int k = 0; k < Nz; k++)
    for (int j = 0; j < Ny; j++)
      for (int i = 0; i < Nx; i++)
      {
        xloc = mygrid->chrloc[0][i];
        yloc = mygrid->chrloc[1][j];
        zloc = mygrid->chrloc[2][k];

        ih = (int)((xloc - xmin) / dx_in);
        axh = (xloc - xmin) / dx_in - ih;
        wh[0] = 1 - axh;
        wh[1] = axh;
        ihleft = (ih + Nx_in - 1) % (Nx_in - 1);
        ihright = (ih + 1) % (Nx_in - 1);

        if (xloc >= plasmarmin[0] && xloc <= plasmarmax[0])
          if (yloc >= plasmarmin[1] && yloc <= plasmarmax[1])
            if (zloc >= plasmarmin[2] && zloc <= plasmarmax[2])
            {
              denValue = wh[0] * density[ihleft] + wh[1] * density[ihright];
              if (denValue > 0)
              {
                velValue = wh[0] * velocity[ihleft] + wh[1] * velocity[ihright];
                weight = denValue / particlePerCell;
                xloc -= 0.5*dx;
                yloc -= 0.5*dy;
                zloc -= 0.5*dz;
                myuy = wh[0] * uy[ihleft] + wh[1] * uy[ihright];
                myuz = wh[0] * uz[ihleft] + wh[1] * uz[ihright];
                for (int ip = 0; ip < particlePerCellXYZ[0]; ip++)
                  for (int jp = 0; jp < particlePerCellXYZ[1]; jp++)
                    for (int kp = 0; kp < particlePerCellXYZ[2]; kp++)
                    {
                      r0(counter) = xloc + dxp*(ip + 0.5);
                      r1(counter) = yloc + dyp*(jp + 0.5);
                      r2(counter) = zloc + dzp*(kp + 0.5);
                      u0(counter) = 0;//velValue/(sqrt(1-velValue*velValue));
                      u1(counter) = myuy;
                      u2(counter) = myuz;
                      w(counter) = weight;
                      marker(counter) = (counter + disp);
                      if (isTestSpecies)
                        w(counter) = (double)(counter + disp);
                      counter++;
                    }
              }
            }
      }
  free(density);
  free(velocity);
  free(x_in);
  free(uy);
  free(uz);
}

void SPECIE::createStretchedParticlesWithinFromButFromFile1D(double plasmarmin[3], double plasmarmax[3], int oldNumberOfParticles, uint64_t disp, std::string name){
  int counter = oldNumberOfParticles;
  double xloc, yloc, zloc;
  double myx, myy, myz;
  double mydx, mydy, mydz;
  double mycsix, mycsiy, mycsiz;
  double csilocx, csilocy, csilocz;
  int Nx = mygrid->Nloc[0];
  int Ny = mygrid->Nloc[1];
  int Nz = mygrid->Nloc[2];
  double dx = mygrid->dr[0];
  double dy = mygrid->dr[1];
  double dz = mygrid->dr[2];
  double dxp = dx / particlePerCellXYZ[0];
  double dyp = dy / particlePerCellXYZ[1];
  double dzp = dz / particlePerCellXYZ[2];
  double  weight;
  double myuy, myuz;

  std::ifstream fileDensity(name.c_str(), std::ifstream::in);
  int Nx_in;
  double xmin, xmax, dx_in;
  double temp;
  fileDensity >> Nx_in;
  fileDensity >> temp;
  fileDensity >> xmin >> xmax;

  double *density, *velocity, *x_in, *uy, *uz;
  density = (double*)malloc(sizeof(double)*Nx_in);
  velocity = (double*)malloc(sizeof(double)*Nx_in);
  x_in = (double*)malloc(sizeof(double)*Nx_in);
  uy = (double*)malloc(sizeof(double)*Nx_in);
  uz = (double*)malloc(sizeof(double)*Nx_in);
  for (int i = 0; i < Nx_in; i++){
    fileDensity >> x_in[i];
    fileDensity >> density[i];
    fileDensity >> velocity[i];
    fileDensity >> uy[i];
    fileDensity >> uz[i];
  }

  xmin = x_in[0];
  xmax = x_in[Nx_in - 1];
  dx_in = x_in[1] - x_in[0];

  double axh, wh[2], denValue, velValue;
  int ih, ihleft, ihright;

  for (int k = 0; k < Nz; k++){
    zloc = mygrid->chrloc[2][k];
    csilocz = mygrid->csiminloc[2] + dz*k;
    for (int j = 0; j < Ny; j++){
      yloc = mygrid->chrloc[1][j];
      csilocy = mygrid->csiminloc[1] + dy*j;
      for (int i = 0; i < Nx; i++){
        xloc = mygrid->chrloc[0][i];
        csilocx = mygrid->csiminloc[0] + dx*i;

        ih = (int)((xloc - xmin) / dx_in);
        axh = (xloc - xmin) / dx_in - ih;
        wh[0] = 1 - axh;
        wh[1] = axh;
        ihleft = (ih + Nx_in - 1) % (Nx_in - 1);
        ihright = (ih + 1) % (Nx_in - 1);

        if (xloc >= plasmarmin[0] && xloc <= plasmarmax[0]){
          if (yloc >= plasmarmin[1] && yloc <= plasmarmax[1]){
            if (zloc >= plasmarmin[2] && zloc <= plasmarmax[2])
            {
              denValue = wh[0] * density[ihleft] + wh[1] * density[ihright];

              if (denValue > 0)
              {
                velValue = wh[0] * velocity[ihleft] + wh[1] * velocity[ihright];
                weight = denValue / particlePerCell;
                myuy = wh[0] * uy[ihleft] + wh[1] * uy[ihright];
                myuz = wh[0] * uz[ihleft] + wh[1] * uz[ihright];
                for (int kp = 0; kp < particlePerCellXYZ[2]; kp++){
                  mycsiz = csilocz + dzp*(kp + 0.5);
                  myz = mygrid->stretchGrid(mycsiz, 2);
                  mydz = mygrid->derivativeStretchingFunction(mycsiz, 2);

                  for (int jp = 0; jp < particlePerCellXYZ[1]; jp++){
                    mycsiy = csilocy + dyp*(jp + 0.5);
                    myy = mygrid->stretchGrid(mycsiy, 1);
                    mydy = mygrid->derivativeStretchingFunction(mycsiy, 1);

                    for (int ip = 0; ip < particlePerCellXYZ[0]; ip++){
                      mycsix = csilocx + dxp*(ip + 0.5);
                      myx = mygrid->stretchGrid(mycsix, 0);
                      mydx = mygrid->derivativeStretchingFunction(mycsix, 0);

                      r0(counter) = myx;
                      r1(counter) = myy;
                      r2(counter) = myz;
                      u0(counter) = 0;//velValue/(sqrt(1-velValue*velValue));
                      u1(counter) = myuy;
                      u2(counter) = myuz;
                      w(counter) = weight*mydx*mydy*mydz;
                      if (isTestSpecies)
                        w(counter) = (double)(counter + disp);
                      counter++;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  free(density);
  free(velocity);
  free(x_in);
  free(uy);
  free(uz);
}

void SPECIE::setNumberOfParticlePerCell(){
  particlePerCell = 1;
  for (int c = 0; c < 3; c++){
    if (c < mygrid->getDimensionality()){
      particlePerCell *= particlePerCellXYZ[c];   //number of particles per cell(npc) total = npcx*npcy*npcz
    }
    else{
      particlePerCellXYZ[c] = 1;
    }
  }
}
void SPECIE::setLocalPlasmaMinimaAndMaxima(double *plasmarmin, double *plasmarmax){
  for (int c = 0; c < 3; c++){
    if (c < mygrid->getDimensionality()){
      plasmarmin[c] = MAX(plasma.params.rminbox[c], mygrid->rminloc[c]);
      plasmarmax[c] = MIN(plasma.params.rmaxbox[c], mygrid->rmaxloc[c]);
    }
    else{
      plasmarmin[c] = mygrid->rminloc[c];
      plasmarmax[c] = mygrid->rmaxloc[c];
    }
  }

}
uint64_t SPECIE::getSumNewParticlesOfAllPreviousProcessors(int number){
  int* NpartLoc = new int[mygrid->nproc];
  NpartLoc[mygrid->myid] = number;

  MPI_Allgather(MPI_IN_PLACE, 1, MPI_INT, NpartLoc, 1, MPI_INT, MPI_COMM_WORLD);
  uint64_t disp = lastParticle;
  for (int pp = 0; pp < mygrid->myid; pp++)
    disp += NpartLoc[pp];
  for (int pp = 0; pp < mygrid->nproc; pp++)
    lastParticle += (uint64_t)NpartLoc[pp];
  delete[] NpartLoc;
  return disp;
}

void SPECIE::creation()
{
  if (mygrid->withParticles == NO)
    return;

  if (flagWithMarker){
    Ncomp = 8;
  }

  computeParticleMassChargeCoupling();
  setNumberOfParticlePerCell();

  double plasmarmin[3], plasmarmax[3];
  setLocalPlasmaMinimaAndMaxima(plasmarmin, plasmarmax);
  Np = getNumberOfParticlesWithin(plasmarmin, plasmarmax);
  allocate_species();

  uint64_t disp = getSumNewParticlesOfAllPreviousProcessors(Np);

  if (mygrid->isStretched())
    createStretchedParticlesWithinFrom(plasmarmin, plasmarmax, 0, disp);
  else
    createParticlesWithinFrom(plasmarmin, plasmarmax, 0, disp);

}


void SPECIE::creationFromFile1D(std::string name){
  std::ifstream fileDensity(name.c_str(), std::ifstream::in);

  if (mygrid->withParticles == NO)
    return;

  int n_points;
  double plasmaXMin, plasmaXMax, plasmarmin[3], plasmarmax[3];
  double temp;
  fileDensity >> n_points;
  fileDensity >> temp;
  fileDensity >> plasmaXMin >> plasmaXMax;
  fileDensity.close();

  SPECIE::computeParticleMassChargeCoupling();
  particlePerCell = 1;
  for (int c = 0; c < 3; c++){
    if (c < 1){
      plasmarmin[c] = MAX(plasmaXMin, mygrid->rminloc[c]);
      plasmarmax[c] = MIN(plasmaXMax, mygrid->rmaxloc[c]);
      particlePerCell *= particlePerCellXYZ[c];   //number of particles per cell(npc) total = npcx*npcy*npcz
    }
    else{
      plasmarmin[c] = mygrid->rminloc[c];
      plasmarmax[c] = mygrid->rmaxloc[c];
      particlePerCellXYZ[c] = 1;
    }
  }

  Np = SPECIE::getNumberOfParticlesWithinFromFile1D(plasmarmin, plasmarmax, name);
  allocate_species();

  int* NpartLoc = new int[mygrid->nproc];
  NpartLoc[mygrid->myid] = Np;

  MPI_Allgather(MPI_IN_PLACE, 1, MPI_INT, NpartLoc, 1, MPI_INT, MPI_COMM_WORLD);
  uint64_t disp = 0;
  for (int pp = 0; pp < mygrid->myid; pp++)
    disp += NpartLoc[pp];
  for (int pp = 0; pp < mygrid->nproc; pp++)
    lastParticle += (uint64_t)NpartLoc[pp];


  if (mygrid->isStretched())
    SPECIE::createStretchedParticlesWithinFromButFromFile1D(plasmarmin, plasmarmax, 0, disp, name);
  else
    SPECIE::createParticlesWithinFromButFromFile1D(plasmarmin, plasmarmax, 0, disp, name);

  delete[] NpartLoc;

}
//CREATE PARTICLES IN THE NEW STRIPE OF DOMAIN "grown" from the window movement
// as "create()" but for a smaller reagion of the space
//reallocation occurs

void SPECIE::move_window()
{
  if (!mygrid->withParticles)
    return;
  if (!mygrid->shouldIMove)
    return;

  SPECIE::position_parallel_pbc();

  double plasmarmin[3], plasmarmax[3];
  setLocalPlasmaMinimaAndMaxima(plasmarmin, plasmarmax);
  plasmarmin[0] = mygrid->rmaxloc[0] - mygrid->fmove_mw;

  int newNumberOfParticles, oldNumberOfParticles = Np;

  if (mygrid->rmyid[0] == (mygrid->rnproc[0] - 1)){
    newNumberOfParticles = SPECIE::getNumberOfParticlesWithin(plasmarmin, plasmarmax);
  }
  else{
    newNumberOfParticles = 0;
  }
  Np += newNumberOfParticles;
  reallocate_species();

  uint64_t disp = getSumNewParticlesOfAllPreviousProcessors(newNumberOfParticles);

  if ((mygrid->rmyid[0] == (mygrid->rnproc[0] - 1)) && newNumberOfParticles > 0){

    if (mygrid->isStretched())
      createStretchedParticlesWithinFrom(plasmarmin, plasmarmax, oldNumberOfParticles, disp);
    else
      createParticlesWithinFrom(plasmarmin, plasmarmax, oldNumberOfParticles, disp);
  }

}
//void SPECIE::output_bin(ofstream &ff)
//{
//	if (mygrid->with_particles == NO)
//		return;

//	ff.write((char *)val, sizeof(double)*Ncomp*Np);
//}
void SPECIE::output(std::ofstream &ff)
{
  if (mygrid->withParticles == NO)
    return;

  for (int nn = 0; nn < Np; nn++)
  {
    for (int cc = 0; cc < Ncomp; cc++)
      ff << ru(cc, nn) << ", ";
    ff << "\n";
  }
}

void SPECIE::init_output_diag(std::ofstream &ff){
  if (mygrid->withParticles == NO)
    return;

  if (mygrid->myid == mygrid->master_proc){
    ff << std::setw(myNarrowWidth) << "#step" << " " << std::setw(myWidth) << "time" << " " << std::setw(myWidth) << "Etot";
    ff << " " << std::setw(myWidth) << "Px" << " " << std::setw(myWidth) << "Py" << " " << std::setw(myWidth) << "Pz" << std::endl;
  }
}
void SPECIE::output_diag(int istep, std::ofstream &ff){
  if (mygrid->withParticles == NO)
    return;

  //double extrema[14];
  computeKineticEnergyWExtrems();
  if (mygrid->myid == mygrid->master_proc){
    ff << std::setw(myWidth) << istep << " " << std::setw(myWidth) << mygrid->time << " " << std::setw(myWidth) << totalEnergy;
    for (int c = 0; c < 3; c++){
      ff << " " << std::setw(myWidth) << totalMomentum[c];
    }
    ff << std::endl;
  }
}
void SPECIE::init_output_extrems(std::ofstream &ff){
  if (mygrid->withParticles == NO)
    return;

  if (mygrid->myid == mygrid->master_proc){
    ff << std::setw(myNarrowWidth) << "#step" << " " << std::setw(myWidth) << "time"
       << " " << std::setw(myWidth) << "xmin" << " " << std::setw(myWidth) << "xmax"
       << " " << std::setw(myWidth) << "ymin" << " " << std::setw(myWidth) << "ymax"
       << " " << std::setw(myWidth) << "zmin" << " " << std::setw(myWidth) << "zmax"
       << " " << std::setw(myWidth) << "pxmin" << " " << std::setw(myWidth) << "pxmax"
       << " " << std::setw(myWidth) << "pymin" << " " << std::setw(myWidth) << "pymax"
       << " " << std::setw(myWidth) << "pzmin" << " " << std::setw(myWidth) << "pzmax"
       << " " << std::setw(myWidth) << "gammamin" << " " << std::setw(myWidth) << "gammamax" << std::endl;
  }
}
void SPECIE::output_extrems(int istep, std::ofstream &ff){
  if (mygrid->withParticles == NO)
    return;

  //double extrema[14];
  computeKineticEnergyWExtrems();
  if (mygrid->myid == mygrid->master_proc){
    ff << std::setw(myNarrowWidth) << istep << " " << std::setw(myWidth) << mygrid->time;
    for (int c = 0; c < 7; c++){
      ff << " " << std::setw(myWidth) << minima[c] << " " << std::setw(myWidth) << maxima[c];
    }
    ff << std::endl;
  }
}
void SPECIE::init_output_stat(std::ofstream &fdiag, std::ofstream &fextrem){
  if (mygrid->withParticles == NO)
    return;

  init_output_extrems(fextrem);
  init_output_diag(fdiag);
}
void SPECIE::output_stat(int istep, std::ofstream &fdiag, std::ofstream &fextrem, std::ofstream &fspectrum){
  if (mygrid->withParticles == NO)
    return;


  computeKineticEnergyWExtrems();

  if (mygrid->myid == mygrid->master_proc){
    fdiag << std::setw(myWidth) << istep << " " << std::setw(myWidth) << mygrid->time << " " << std::setw(myWidth) << totalEnergy;
    for (int c = 0; c < 3; c++){
      fdiag << " " << std::setw(myWidth) << totalMomentum[c];
    }
    fdiag << std::endl;

    fextrem << std::setw(myNarrowWidth) << istep << " " << std::setw(myWidth) << mygrid->time;
    for (int c = 0; c < 7; c++){
      fextrem << " " << std::setw(myWidth) << minima[c] << " " << std::setw(myWidth) << maxima[c];
    }
    fextrem << std::endl;

    fspectrum << "#" << std::setw(myWidth) << "Ebinmin";
    fspectrum << " " << std::setw(myWidth) << "Ebinmax";
    fspectrum << " " << std::setw(myWidth) << "value";
    fspectrum << std::endl;
    for (int ibin = 0; ibin < spectrum.Nbin; ibin++){
      fspectrum << " " << std::setw(myWidth) << (ibin*spectrum.Dk);
      fspectrum << " " << std::setw(myWidth) << ((ibin + 1)*spectrum.Dk);
      fspectrum << " " << std::setw(myWidth) << (spectrum.values[ibin]);
      fspectrum << std::endl;

    }
  }

}

void SPECIE::outputSpectrum(std::ofstream &fspectrum){
  computeKineticEnergyWExtrems();
  if (mygrid->myid == mygrid->master_proc){
    fspectrum << "#" << std::setw(myWidth) << "Ebinmin";
    fspectrum << " " << std::setw(myWidth) << "Ebinmax";
    fspectrum << " " << std::setw(myWidth) << "value";
    fspectrum << std::endl;
    for (int ibin = 0; ibin < spectrum.Nbin; ibin++){
      fspectrum << " " << std::setw(myWidth) << (ibin*spectrum.Dk);
      fspectrum << " " << std::setw(myWidth) << ((ibin + 1)*spectrum.Dk);
      fspectrum << " " << std::setw(myWidth) << (spectrum.values[ibin]);
      fspectrum << std::endl;

    }
  }
}


void SPECIE::position_advance()
{
  if (mygrid->withParticles == NO)
    return;

  double dt, gamma_i;
  int p;

  dt = mygrid->dt;

  switch (mygrid->getDimensionality())
  {
    case 3:
      //#pragma omp parallel for
      for (p = 0; p < Np; p++)
      {
        gamma_i = 1. / sqrt(1 + u0(p)*u0(p) + u1(p)*u1(p) + u2(p)*u2(p));
        r0(p) += dt*gamma_i*u0(p);
        r1(p) += dt*gamma_i*u1(p);
        r2(p) += dt*gamma_i*u2(p);

      }
      break;
    case 2:

      for (p = 0; p < Np; p++)
      {
        gamma_i = 1. / sqrt(1 + u0(p)*u0(p) + u1(p)*u1(p) + u2(p)*u2(p));
        r0(p) += dt*gamma_i*u0(p);
        r1(p) += dt*gamma_i*u1(p);

      }
      break;
    case 1:

      for (p = 0; p < Np; p++)
      {
        gamma_i = 1. / sqrt(1 + u0(p)*u0(p) + u1(p)*u1(p) + u2(p)*u2(p));
        r0(p) += dt*gamma_i*u0(p);

      }
      break;

  }
}
void SPECIE::position_pbc()
{
  if (mygrid->withParticles == NO)
    return;

  int p;

  /*
      cambiare in questo senso:
      se una particella è persa, la scambio con l'ultima particella attiva e diminiusco di uno il numero di particelle attive
      se anche questa è da buttare via, proseguon nella ricerca di partcielle buone con la penultima attiva fino a quando non trovo una buana
      da sostiuire a quella in esame
      così dicendo riduco Np_loc che è anche l'estremo del ciclo for

      */

  for (p = 0; p<Np; p++)
  {
    if (ru(0, p)> mygrid->rmaxloc[0])
      ru(0, p) -= (mygrid->rmaxloc[0] - mygrid->rminloc[0]);
    if (ru(0, p) < mygrid->rminloc[0])
      ru(0, p) += (mygrid->rmaxloc[0] - mygrid->rminloc[0]);

    if (ru(1, p) > mygrid->rmaxloc[1])
      ru(1, p) -= (mygrid->rmaxloc[1] - mygrid->rminloc[1]);
    if (ru(1, p) < mygrid->rminloc[1])
      ru(1, p) += (mygrid->rmaxloc[1] - mygrid->rminloc[1]);

    if (ru(2, p) > mygrid->rmaxloc[2])
      ru(2, p) -= (mygrid->rmaxloc[2] - mygrid->rminloc[2]);
    if (ru(2, p) < mygrid->rminloc[2])
      ru(2, p) += (mygrid->rmaxloc[2] - mygrid->rminloc[2]);
  }
}
void SPECIE::position_parallel_pbc()
{
  if (mygrid->withParticles == NO)
    return;

  int p, c;
  int nlost, nnew, nold;
  int ninright, ninleft, nright, nleft;
  ninright = ninleft = nright = nleft = 0;
  double *sendr_buffer = NULL, *sendl_buffer = NULL, *recv_buffer = NULL;
  MPI_Status status;
  int iright, ileft;
  /*
      si potrebbe cambiare in questo senso:
      se una particella è persa, la scambio con l'ultima particella attiva e diminiusco di uno il numero di particelle attive
      se anche questa è da buttare via, proseguon nella ricerca di partcielle buone con la penultima attiva fino a quando non trovo una buana
      da sostiuire a quella in esame
      così dicendo riduco Np_loc che è anche l'estremo del ciclo for

      */


  for (int direction = 0; direction < mygrid->getDimensionality(); direction++)
  {
    nlost = 0;
    ninright = ninleft = nright = nleft = 0;
    for (p = 0; p < Np; p++)
    {
      if (ru(direction, p) > mygrid->rmaxloc[direction])
      {
        nlost++;
        nright++;
        if (mygrid->rmyid[direction] == mygrid->rnproc[direction] - 1)
          ru(direction, p) -= (mygrid->rmax[direction] - mygrid->rmin[direction]);
        sendr_buffer = (double*)realloc(sendr_buffer, nright*Ncomp*sizeof(double));
        for (c = 0; c < Ncomp; c++)
          sendr_buffer[c + Ncomp*(nright - 1)] = ru(c, p);

      }
      else if (ru(direction, p) < mygrid->rminloc[direction])
      {
        nlost++;
        nleft++;
        if (mygrid->rmyid[direction] == 0)
          ru(direction, p) += (mygrid->rmax[direction] - mygrid->rmin[direction]);
        sendl_buffer = (double*)realloc(sendl_buffer, nleft*Ncomp*sizeof(double));
        for (c = 0; c < Ncomp; c++)
          sendl_buffer[c + Ncomp*(nleft - 1)] = ru(c, p);
      }
      else
      {
        for (c = 0; c < Ncomp; c++)
          ru(c, p - nlost) = ru(c, p);
      }
    }
    MPI_Cart_shift(mygrid->cart_comm, direction, 1, &ileft, &iright);
    // ====== send right receive from left
    ninleft = 0;
    MPI_Sendrecv(&nright, 1, MPI_INT, iright, 13,
                 &ninleft, 1, MPI_INT, ileft, 13,
                 MPI_COMM_WORLD, &status);
    nnew = ninleft;

    // ====== send left receive from right
    ninright = 0;
    MPI_Sendrecv(&nleft, 1, MPI_INT, ileft, 13,
                 &ninright, 1, MPI_INT, iright, 13,
                 MPI_COMM_WORLD, &status);
    nnew += ninright;
    recv_buffer = (double*)realloc(recv_buffer, nnew*Ncomp*sizeof(double));
    // ====== send right receive from left
    MPI_Sendrecv(sendr_buffer, nright*Ncomp, MPI_DOUBLE, iright, 13,
                 recv_buffer, ninleft*Ncomp, MPI_DOUBLE, ileft, 13,
                 MPI_COMM_WORLD, &status);
    MPI_Sendrecv(sendl_buffer, nleft*Ncomp, MPI_DOUBLE, ileft, 13,
                 (recv_buffer + ninleft*Ncomp), ninright*Ncomp, MPI_DOUBLE, iright, 13,
                 MPI_COMM_WORLD, &status);
    nold = Np;
    Np = Np - nlost + nnew;


    reallocate_species();

    for (int pp = 0; pp < nnew; pp++){
      for (c = 0; c < Ncomp; c++){
        ru(c, pp + nold - nlost) = recv_buffer[pp*Ncomp + c];
      }
    }
  }

  free(sendl_buffer);
  free(sendr_buffer);
  free(recv_buffer);
}
void SPECIE::position_obc()
{
  if (mygrid->withParticles == NO)
    return;

  int p, c;
  int nlost = 0;

  /*
      cambiare in questo senso:
      se una particella è persa, la scambio con l'ultima particella attiva e diminiusco di uno il numero di particelle attive
      se anche questa è da buttare via, proseguon nella ricerca di partcielle buone con la penultima attiva fino a quando non trovo una buana
      da sostiuire a quella in esame
      così dicendo riduco Np_loc che è anche l'estremo del ciclo for

      */
  for (p = 0; p < Np; p++)
  {
    for (c = 0; c < Ncomp; c++)
      ru(c, p - nlost) = ru(c, p);

    if ((ru(0, p) > mygrid->rmaxloc[0]) || (ru(0, p) < mygrid->rminloc[0]))
    {
      nlost++;
      continue;
    }

    if (ru(1, p) > mygrid->rmaxloc[1])
      ru(1, p) -= (mygrid->rmaxloc[1] - mygrid->rminloc[1]);
    if (ru(1, p) < mygrid->rminloc[1])
      ru(1, p) += (mygrid->rmaxloc[1] - mygrid->rminloc[1]);

    if (ru(2, p) > mygrid->rmaxloc[2])
      ru(2, p) -= (mygrid->rmaxloc[2] - mygrid->rminloc[2]);
    if (ru(2, p) < mygrid->rminloc[2])
      ru(2, p) += (mygrid->rmaxloc[2] - mygrid->rminloc[2]);



  }
  Np -= nlost;
  reallocate_species();
}


void SPECIE::momenta_advance(EM_FIELD *ebfield)
{

  energyExtremesFlag = false;

  if (mygrid->withParticles == NO)
    return;
  if (mygrid->isStretched()){
    SPECIE::momentaStretchedAdvance(ebfield);
    return;
  }
  double dt, gamma_i;
  int p, c;  // particle_int, component_int
  int i, j, k, i1, j1, k1, i2, j2, k2;
  //int indexMaxQuadraticShape[]={1,4};
  int hii[3], wii[3];           // half integer index,   whole integer index
  double hiw[3][3], wiw[3][3];  // half integer weight,  whole integer weight
  double rr, rh, rr2, rh2;          // local coordinate to integer grid point and to half integer,     local coordinate squared
  double dvol, xx[3];           // tensor_product,       absolute particle position
  double E[3], B[3];
  double u_plus[3], u_minus[3], u_prime[3], tee[3], ess[3], dummy;

  dt = mygrid->dt;

  double *myfield = ebfield->getDataPointer();
  int edge = mygrid->getEdge();
  int ebComp = ebfield->getNcomp();
  int N_grid[3];
  ebfield->writeN_grid(N_grid);
  switch (mygrid->getDimensionality())
  {

    case 3:
      for (p = 0; p < Np; p++)
      {
        //gamma_i=1./sqrt(1+u0(p)*u0(p)+u1(p)*u1(p)+u2(p)*u2(p));
        for (c = 0; c < 3; c++){
          xx[c] = pData[c + p*Ncomp];
          hiw[c][1] = wiw[c][1] = 1;
          hii[c] = wii[c] = 0;
        }
        for (c = 0; c < 3; c++)
        {
          rr = mygrid->dri[c] * (xx[c] - mygrid->rminloc[c]);
          rh = rr - 0.5;
          wii[c] = (int)floor(rr + 0.5); //whole integer int
          hii[c] = (int)floor(rr);     //half integer int
          rr -= wii[c];
          rh -= hii[c];
          rr2 = rr*rr;
          rh2 = rh*rh;

          wiw[c][1] = 0.75 - rr2;
          wiw[c][2] = 0.5*(0.25 + rr2 + rr);
          wiw[c][0] = 1. - wiw[c][1] - wiw[c][2];

          hiw[c][1] = 0.75 - rh2;
          hiw[c][2] = 0.5*(0.25 + rh2 + rh);
          hiw[c][0] = 1. - hiw[c][1] - hiw[c][2];
        }
        E[0] = E[1] = E[2] = B[0] = B[1] = B[2] = 0;

        for (k = 0; k < 3; k++)
        {
          k1 = k + wii[2] - 1;
          k2 = k + hii[2] - 1;
          for (j = 0; j < 3; j++)
          {
            j1 = j + wii[1] - 1;
            j2 = j + hii[1] - 1;
            for (i = 0; i < 3; i++)
            {
              i1 = i + wii[0] - 1;
              i2 = i + hii[0] - 1;
              double EX, EY, EZ;
              EX = myfield[my_indice(edge,1, 1, 0, i2, j1, k1, N_grid[0], N_grid[1], N_grid[2], ebComp)];
              EY = myfield[my_indice(edge,1, 1, 1, i1, j2, k1, N_grid[0], N_grid[1], N_grid[2], ebComp)];
              EZ = myfield[my_indice(edge,1, 1, 2, i1, j1, k2, N_grid[0], N_grid[1], N_grid[2], ebComp)];
              double BX, BY, BZ;
              BX = myfield[my_indice(edge,1, 1, 3, i1, j2, k2, N_grid[0], N_grid[1], N_grid[2], ebComp)];
              BY = myfield[my_indice(edge,1, 1, 4, i2, j1, k2, N_grid[0], N_grid[1], N_grid[2], ebComp)];
              BZ = myfield[my_indice(edge,1, 1, 5, i2, j2, k1, N_grid[0], N_grid[1], N_grid[2], ebComp)];

              dvol = hiw[0][i] * wiw[1][j] * wiw[2][k];
              E[0] += EX*dvol;  //Ex
              dvol = wiw[0][i] * hiw[1][j] * wiw[2][k];
              E[1] += EY*dvol;  //Ey
              dvol = wiw[0][i] * wiw[1][j] * hiw[2][k];
              E[2] += EZ*dvol;  //Ez

              dvol = wiw[0][i] * hiw[1][j] * hiw[2][k];
              B[0] += BX*dvol;  //Bx
              dvol = hiw[0][i] * wiw[1][j] * hiw[2][k];
              B[1] += BY*dvol;  //By
              dvol = hiw[0][i] * hiw[1][j] * wiw[2][k];
              B[2] += BZ*dvol;  //Bz
            }
          }
        }

        u_minus[0] = pData[3 + p*Ncomp] + 0.5*dt*coupling*E[0];
        u_minus[1] = pData[4 + p*Ncomp] + 0.5*dt*coupling*E[1];
        u_minus[2] = pData[5 + p*Ncomp] + 0.5*dt*coupling*E[2];

        gamma_i = 1. / sqrt(1 + u_minus[0] * u_minus[0] + u_minus[1] * u_minus[1] + u_minus[2] * u_minus[2]);

        tee[0] = 0.5*dt*coupling*B[0] * gamma_i;
        tee[1] = 0.5*dt*coupling*B[1] * gamma_i;
        tee[2] = 0.5*dt*coupling*B[2] * gamma_i;

        u_prime[0] = u_minus[0] + (u_minus[1] * tee[2] - u_minus[2] * tee[1]);
        u_prime[1] = u_minus[1] + (u_minus[2] * tee[0] - u_minus[0] * tee[2]);
        u_prime[2] = u_minus[2] + (u_minus[0] * tee[1] - u_minus[1] * tee[0]);

        dummy = 1 / (1 + tee[0] * tee[0] + tee[1] * tee[1] + tee[2] * tee[2]);

        ess[0] = 2 * dummy*tee[0];
        ess[1] = 2 * dummy*tee[1];
        ess[2] = 2 * dummy*tee[2];

        u_plus[0] = u_minus[0] + u_prime[1] * ess[2] - u_prime[2] * ess[1];
        u_plus[1] = u_minus[1] + u_prime[2] * ess[0] - u_prime[0] * ess[2];
        u_plus[2] = u_minus[2] + u_prime[0] * ess[1] - u_prime[1] * ess[0];

        pData[3 + p*Ncomp] = (u_plus[0] + 0.5*dt*coupling*E[0]);
        pData[4 + p*Ncomp] = (u_plus[1] + 0.5*dt*coupling*E[1]);
        pData[5 + p*Ncomp] = (u_plus[2] + 0.5*dt*coupling*E[2]);
      }
      break;

    case 2:
      for (p = 0; p < Np; p++)
      {
        //gamma_i=1./sqrt(1+u0(p)*u0(p)+u1(p)*u1(p)+u2(p)*u2(p));
        for (c = 0; c < 3; c++)
        {
          xx[c] = pData[c + p*Ncomp];
          hiw[c][1] = wiw[c][1] = 1;
          hii[c] = wii[c] = 0;
        }
        for (c = 0; c < 2; c++)
        {
          rr = mygrid->dri[c] * (xx[c] - mygrid->rminloc[c]);
          rh = rr - 0.5;
          wii[c] = (int)floor(rr + 0.5); //whole integer int
          hii[c] = (int)floor(rr);     //half integer int
          rr -= wii[c];
          rh -= hii[c];
          rr2 = rr*rr;
          rh2 = rh*rh;

          wiw[c][1] = 0.75 - rr2;
          wiw[c][2] = 0.5*(0.25 + rr2 + rr);
          wiw[c][0] = 1. - wiw[c][1] - wiw[c][2];

          hiw[c][1] = 0.75 - rh2;
          hiw[c][2] = 0.5*(0.25 + rh2 + rh);
          hiw[c][0] = 1. - hiw[c][1] - hiw[c][2];
        }
        E[0] = E[1] = E[2] = B[0] = B[1] = B[2] = 0;

        k1 = k2 = 0;
        for (j = 0; j < 3; j++)
        {
          j1 = j + wii[1] - 1;
          j2 = j + hii[1] - 1;
          for (i = 0; i < 3; i++){
            i1 = i + wii[0] - 1;
            i2 = i + hii[0] - 1;
#ifndef OLD_ACCESS
            double EX, EY, EZ;
            EX = myfield[my_indice(edge,1, 0, 0, i2, j1, k1, N_grid[0], N_grid[1], N_grid[2], ebComp)];
            EY = myfield[my_indice(edge,1, 0, 1, i1, j2, k1, N_grid[0], N_grid[1], N_grid[2], ebComp)];
            EZ = myfield[my_indice(edge,1, 0, 2, i1, j1, k2, N_grid[0], N_grid[1], N_grid[2], ebComp)];
            double BX, BY, BZ;
            BX = myfield[my_indice(edge,1, 0, 3, i1, j2, k2, N_grid[0], N_grid[1], N_grid[2], ebComp)];
            BY = myfield[my_indice(edge,1, 0, 4, i2, j1, k2, N_grid[0], N_grid[1], N_grid[2], ebComp)];
            BZ = myfield[my_indice(edge,1, 0, 5, i2, j2, k1, N_grid[0], N_grid[1], N_grid[2], ebComp)];

            dvol = hiw[0][i] * wiw[1][j];
            E[0] += EX*dvol;  //Ex
            dvol = wiw[0][i] * hiw[1][j];
            E[1] += EY*dvol;  //Ey
            dvol = wiw[0][i] * wiw[1][j];
            E[2] += EZ*dvol;  //Ez

            dvol = wiw[0][i] * hiw[1][j];
            B[0] += BX*dvol;  //Bx
            dvol = hiw[0][i] * wiw[1][j];
            B[1] += BY*dvol;  //By
            dvol = hiw[0][i] * hiw[1][j];
            B[2] += BZ*dvol;  //Bz
#else
            dvol = hiw[0][i] * wiw[1][j],
                E[0] += ebfield->E0(i2, j1, k1)*dvol;  //Ex
            dvol = wiw[0][i] * hiw[1][j],
                E[1] += ebfield->E1(i1, j2, k1)*dvol;  //Ey
            dvol = wiw[0][i] * wiw[1][j],
                E[2] += ebfield->E2(i1, j1, k2)*dvol;  //Ez

            dvol = wiw[0][i] * hiw[1][j],
                B[0] += ebfield->B0(i1, j2, k2)*dvol;  //Bx
            dvol = hiw[0][i] * wiw[1][j],
                B[1] += ebfield->B1(i2, j1, k2)*dvol;  //By
            dvol = hiw[0][i] * hiw[1][j],
                B[2] += ebfield->B2(i2, j2, k1)*dvol;  //Bz
#endif
          }
        }

        u_minus[0] = pData[3 + p*Ncomp] + 0.5*dt*coupling*E[0];
        u_minus[1] = pData[4 + p*Ncomp] + 0.5*dt*coupling*E[1];
        u_minus[2] = pData[5 + p*Ncomp] + 0.5*dt*coupling*E[2];

        gamma_i = 1. / sqrt(1 + u_minus[0] * u_minus[0] + u_minus[1] * u_minus[1] + u_minus[2] * u_minus[2]);

        tee[0] = 0.5*dt*coupling*B[0] * gamma_i;
        tee[1] = 0.5*dt*coupling*B[1] * gamma_i;
        tee[2] = 0.5*dt*coupling*B[2] * gamma_i;

        u_prime[0] = u_minus[0] + (u_minus[1] * tee[2] - u_minus[2] * tee[1]);
        u_prime[1] = u_minus[1] + (u_minus[2] * tee[0] - u_minus[0] * tee[2]);
        u_prime[2] = u_minus[2] + (u_minus[0] * tee[1] - u_minus[1] * tee[0]);

        dummy = 1 / (1 + tee[0] * tee[0] + tee[1] * tee[1] + tee[2] * tee[2]);

        ess[0] = 2 * dummy*tee[0];
        ess[1] = 2 * dummy*tee[1];
        ess[2] = 2 * dummy*tee[2];

        u_plus[0] = u_minus[0] + u_prime[1] * ess[2] - u_prime[2] * ess[1];
        u_plus[1] = u_minus[1] + u_prime[2] * ess[0] - u_prime[0] * ess[2];
        u_plus[2] = u_minus[2] + u_prime[0] * ess[1] - u_prime[1] * ess[0];

        pData[3 + p*Ncomp] = (u_plus[0] + 0.5*dt*coupling*E[0]);
        pData[4 + p*Ncomp] = (u_plus[1] + 0.5*dt*coupling*E[1]);
        pData[5 + p*Ncomp] = (u_plus[2] + 0.5*dt*coupling*E[2]);
      }
      break;

    case 1:
      for (p = 0; p < Np; p++)
      {
        //gamma_i=1./sqrt(1+u0(p)*u0(p)+u1(p)*u1(p)+u2(p)*u2(p));
        for (c = 0; c < 3; c++)
        {
          xx[c] = pData[c + p*Ncomp];
          hiw[c][1] = wiw[c][1] = 1;
          hii[c] = wii[c] = 0;
        }
        for (c = 0; c < 1; c++)
        {
          rr = mygrid->dri[c] * (xx[c] - mygrid->rminloc[c]);
          rh = rr - 0.5;
          wii[c] = (int)floor(rr + 0.5); //whole integer int
          hii[c] = (int)floor(rr);     //half integer int
          rr -= wii[c];
          rh -= hii[c];
          rr2 = rr*rr;
          rh2 = rh*rh;

          wiw[c][1] = 0.75 - rr2;
          wiw[c][2] = 0.5*(0.25 + rr2 + rr);
          wiw[c][0] = 1. - wiw[c][1] - wiw[c][2];

          hiw[c][1] = 0.75 - rh2;
          hiw[c][2] = 0.5*(0.25 + rh2 + rh);
          hiw[c][0] = 1. - hiw[c][1] - hiw[c][2];
        }
        E[0] = E[1] = E[2] = B[0] = B[1] = B[2] = 0;

        k1 = k2 = j1 = j2 = 0;
        for (i = 0; i < 3; i++)
        {
          i1 = i + wii[0] - 1;
          i2 = i + hii[0] - 1;
          double EX, EY, EZ;
          EX = myfield[my_indice(edge,0, 0, 0, i2, j1, k1, N_grid[0], N_grid[1], N_grid[2], ebComp)];
          EY = myfield[my_indice(edge,0, 0, 1, i1, j2, k1, N_grid[0], N_grid[1], N_grid[2], ebComp)];
          EZ = myfield[my_indice(edge,0, 0, 2, i1, j1, k2, N_grid[0], N_grid[1], N_grid[2], ebComp)];
          double BX, BY, BZ;
          BX = myfield[my_indice(edge,0, 0, 3, i1, j2, k2, N_grid[0], N_grid[1], N_grid[2], ebComp)];
          BY = myfield[my_indice(edge,0, 0, 4, i2, j1, k2, N_grid[0], N_grid[1], N_grid[2], ebComp)];
          BZ = myfield[my_indice(edge,0, 0, 5, i2, j2, k1, N_grid[0], N_grid[1], N_grid[2], ebComp)];

          dvol = hiw[0][i];
          E[0] += EX*dvol;  //Ex
          dvol = wiw[0][i];
          E[1] += EY*dvol;  //Ey
          dvol = wiw[0][i];
          E[2] += EZ*dvol;  //Ez

          dvol = wiw[0][i];
          B[0] += BX*dvol;  //Bx
          dvol = hiw[0][i];
          B[1] += BY*dvol;  //By
          dvol = hiw[0][i];
          B[2] += BZ*dvol;  //Bz
        }

        u_minus[0] = pData[3 + p*Ncomp] + 0.5*dt*coupling*E[0];
        u_minus[1] = pData[4 + p*Ncomp] + 0.5*dt*coupling*E[1];
        u_minus[2] = pData[5 + p*Ncomp] + 0.5*dt*coupling*E[2];

        gamma_i = 1. / sqrt(1 + u_minus[0] * u_minus[0] + u_minus[1] * u_minus[1] + u_minus[2] * u_minus[2]);

        tee[0] = 0.5*dt*coupling*B[0] * gamma_i;
        tee[1] = 0.5*dt*coupling*B[1] * gamma_i;
        tee[2] = 0.5*dt*coupling*B[2] * gamma_i;

        u_prime[0] = u_minus[0] + (u_minus[1] * tee[2] - u_minus[2] * tee[1]);
        u_prime[1] = u_minus[1] + (u_minus[2] * tee[0] - u_minus[0] * tee[2]);
        u_prime[2] = u_minus[2] + (u_minus[0] * tee[1] - u_minus[1] * tee[0]);

        dummy = 1 / (1 + tee[0] * tee[0] + tee[1] * tee[1] + tee[2] * tee[2]);

        ess[0] = 2 * dummy*tee[0];
        ess[1] = 2 * dummy*tee[1];
        ess[2] = 2 * dummy*tee[2];

        u_plus[0] = u_minus[0] + u_prime[1] * ess[2] - u_prime[2] * ess[1];
        u_plus[1] = u_minus[1] + u_prime[2] * ess[0] - u_prime[0] * ess[2];
        u_plus[2] = u_minus[2] + u_prime[0] * ess[1] - u_prime[1] * ess[0];

        pData[3 + p*Ncomp] = (u_plus[0] + 0.5*dt*coupling*E[0]);
        pData[4 + p*Ncomp] = (u_plus[1] + 0.5*dt*coupling*E[1]);
        pData[5 + p*Ncomp] = (u_plus[2] + 0.5*dt*coupling*E[2]);
      }
      break;
  }
}


void SPECIE::momenta_advance_with_friction(EM_FIELD *ebfield, double lambda)
{
  energyExtremesFlag = false;

  if (mygrid->withParticles == NO)
    return;
  if (mygrid->isStretched()){
    //SPECIE::momentaStretchedAdvance(ebfield);
    std::cout << "RR not yet implemented with stretched grid!" << std::endl;
    exit(13);
    return;
  }
  double dt, gamma_i;
  int p, c;  // particle_int, component_int
  int i, j, k, i1, j1, k1, i2, j2, k2;
  //int indexMaxQuadraticShape[]={1,4};
  int hii[3], wii[3];           // half integer index,   whole integer index
  double hiw[3][3], wiw[3][3];  // half integer weight,  whole integer weight
  double rr, rh, rr2, rh2;          // local coordinate to integer grid point and to half integer,     local coordinate squared
  double dvol, xx[3];           // tensor_product,       absolute particle position
  double E[3], B[3];
  double u_plus[3], u_minus[3], u_prime[3], tee[3], ess[3], dummy;
  double oldP[3];
  double pn[3]; double vn[3]; double fLorentz[3]; double fLorentz2; double vdotE2; double gamman;

  dt = mygrid->dt;

  double RRcoefficient = 4.0 / 3.0*M_PI*(CLASSICAL_ELECTRON_RADIUS / lambda);

  switch (mygrid->getDimensionality())
  {

    case 3:
      for (p = 0; p < Np; p++)
      {
        //gamma_i=1./sqrt(1+u0(p)*u0(p)+u1(p)*u1(p)+u2(p)*u2(p));
        for (c = 0; c < 3; c++)
        {
          xx[c] = ru(c, p);
          hiw[c][1] = wiw[c][1] = 1;
          hii[c] = wii[c] = 0;
        }
        for (c = 0; c < 3; c++)
        {
          rr = mygrid->dri[c] * (xx[c] - mygrid->rminloc[c]);
          rh = rr - 0.5;
          wii[c] = (int)floor(rr + 0.5); //whole integer int
          hii[c] = (int)floor(rr);     //half integer int
          rr -= wii[c];
          rh -= hii[c];
          rr2 = rr*rr;
          rh2 = rh*rh;

          wiw[c][1] = 0.75 - rr2;
          wiw[c][2] = 0.5*(0.25 + rr2 + rr);
          wiw[c][0] = 1. - wiw[c][1] - wiw[c][2];

          hiw[c][1] = 0.75 - rh2;
          hiw[c][2] = 0.5*(0.25 + rh2 + rh);
          hiw[c][0] = 1. - hiw[c][1] - hiw[c][2];
        }
        E[0] = E[1] = E[2] = B[0] = B[1] = B[2] = 0;

        for (k = 0; k < 3; k++)
        {
          k1 = k + wii[2] - 1;
          k2 = k + hii[2] - 1;
          for (j = 0; j < 3; j++)
          {
            j1 = j + wii[1] - 1;
            j2 = j + hii[1] - 1;
            for (i = 0; i < 3; i++)
            {
              i1 = i + wii[0] - 1;
              i2 = i + hii[0] - 1;
              dvol = hiw[0][i] * wiw[1][j] * wiw[2][k],
                  E[0] += ebfield->E0(i2, j1, k1)*dvol;  //Ex
              dvol = wiw[0][i] * hiw[1][j] * wiw[2][k],
                  E[1] += ebfield->E1(i1, j2, k1)*dvol;  //Ey
              dvol = wiw[0][i] * wiw[1][j] * hiw[2][k],
                  E[2] += ebfield->E2(i1, j1, k2)*dvol;  //Ez

              dvol = wiw[0][i] * hiw[1][j] * hiw[2][k],
                  B[0] += ebfield->B0(i1, j2, k2)*dvol;  //Bx
              dvol = hiw[0][i] * wiw[1][j] * hiw[2][k],
                  B[1] += ebfield->B1(i2, j1, k2)*dvol;  //By
              dvol = hiw[0][i] * hiw[1][j] * wiw[2][k],
                  B[2] += ebfield->B2(i2, j2, k1)*dvol;  //Bz
            }
          }
        }

        oldP[0] = ru(3, p);
        oldP[1] = ru(4, p);
        oldP[2] = ru(5, p);

        u_minus[0] = ru(3, p) + 0.5*dt*coupling*E[0];
        u_minus[1] = ru(4, p) + 0.5*dt*coupling*E[1];
        u_minus[2] = ru(5, p) + 0.5*dt*coupling*E[2];

        gamma_i = 1. / sqrt(1 + u_minus[0] * u_minus[0] + u_minus[1] * u_minus[1] + u_minus[2] * u_minus[2]);

        tee[0] = 0.5*dt*coupling*B[0] * gamma_i;
        tee[1] = 0.5*dt*coupling*B[1] * gamma_i;
        tee[2] = 0.5*dt*coupling*B[2] * gamma_i;

        u_prime[0] = u_minus[0] + (u_minus[1] * tee[2] - u_minus[2] * tee[1]);
        u_prime[1] = u_minus[1] + (u_minus[2] * tee[0] - u_minus[0] * tee[2]);
        u_prime[2] = u_minus[2] + (u_minus[0] * tee[1] - u_minus[1] * tee[0]);

        dummy = 1 / (1 + tee[0] * tee[0] + tee[1] * tee[1] + tee[2] * tee[2]);

        ess[0] = 2 * dummy*tee[0];
        ess[1] = 2 * dummy*tee[1];
        ess[2] = 2 * dummy*tee[2];

        u_plus[0] = u_minus[0] + u_prime[1] * ess[2] - u_prime[2] * ess[1];
        u_plus[1] = u_minus[1] + u_prime[2] * ess[0] - u_prime[0] * ess[2];
        u_plus[2] = u_minus[2] + u_prime[0] * ess[1] - u_prime[1] * ess[0];

        ru(3, p) = (u_plus[0] + 0.5*dt*coupling*E[0]);
        ru(4, p) = (u_plus[1] + 0.5*dt*coupling*E[1]);
        ru(5, p) = (u_plus[2] + 0.5*dt*coupling*E[2]);

        pn[0] = 0.5*(oldP[0] + ru(3, p));
        pn[1] = 0.5*(oldP[1] + ru(4, p));
        pn[2] = 0.5*(oldP[2] + ru(5, p));

        gamman = sqrt(1.0 + (pn[0] * pn[0] + pn[1] * pn[1] + pn[2] * pn[2]));

        vn[0] = pn[0] / gamman;
        vn[1] = pn[1] / gamman;
        vn[2] = pn[2] / gamman;

        fLorentz[0] = coupling*(E[0] + vn[1] * B[2] - vn[2] * B[1]);
        fLorentz[1] = coupling*(E[1] + vn[2] * B[0] - vn[0] * B[2]);
        fLorentz[2] = coupling*(E[2] + vn[0] * B[1] - vn[1] * B[0]);

        fLorentz2 = fLorentz[0] * fLorentz[0] + fLorentz[1] * fLorentz[1] + fLorentz[2] * fLorentz[2];

        vdotE2 = vn[0] * E[0] + vn[1] * E[1] + vn[2] * E[2];
        vdotE2 = vdotE2*vdotE2;

        ru(3, p) -= RRcoefficient*(gamman*gamman)*(fLorentz2 - vdotE2)*vn[0] * dt;
        ru(4, p) -= RRcoefficient*(gamman*gamman)*(fLorentz2 - vdotE2)*vn[1] * dt;
        ru(5, p) -= RRcoefficient*(gamman*gamman)*(fLorentz2 - vdotE2)*vn[2] * dt;


      }
      break;

    case 2:
      for (p = 0; p < Np; p++)
      {
        //gamma_i=1./sqrt(1+u0(p)*u0(p)+u1(p)*u1(p)+u2(p)*u2(p));
        for (c = 0; c < 3; c++)
        {
          xx[c] = ru(c, p);
          hiw[c][1] = wiw[c][1] = 1;
          hii[c] = wii[c] = 0;
        }
        for (c = 0; c < 2; c++)
        {
          rr = mygrid->dri[c] * (xx[c] - mygrid->rminloc[c]);
          rh = rr - 0.5;
          wii[c] = (int)floor(rr + 0.5); //whole integer int
          hii[c] = (int)floor(rr);     //half integer int
          rr -= wii[c];
          rh -= hii[c];
          rr2 = rr*rr;
          rh2 = rh*rh;

          wiw[c][1] = 0.75 - rr2;
          wiw[c][2] = 0.5*(0.25 + rr2 + rr);
          wiw[c][0] = 1. - wiw[c][1] - wiw[c][2];

          hiw[c][1] = 0.75 - rh2;
          hiw[c][2] = 0.5*(0.25 + rh2 + rh);
          hiw[c][0] = 1. - hiw[c][1] - hiw[c][2];
        }
        E[0] = E[1] = E[2] = B[0] = B[1] = B[2] = 0;

        k1 = k2 = 0;
        for (j = 0; j < 3; j++)
        {
          j1 = j + wii[1] - 1;
          j2 = j + hii[1] - 1;
          for (i = 0; i < 3; i++)
          {
            i1 = i + wii[0] - 1;
            i2 = i + hii[0] - 1;
            dvol = hiw[0][i] * wiw[1][j],
                E[0] += ebfield->E0(i2, j1, k1)*dvol;  //Ex
            dvol = wiw[0][i] * hiw[1][j],
                E[1] += ebfield->E1(i1, j2, k1)*dvol;  //Ey
            dvol = wiw[0][i] * wiw[1][j],
                E[2] += ebfield->E2(i1, j1, k2)*dvol;  //Ez

            dvol = wiw[0][i] * hiw[1][j],
                B[0] += ebfield->B0(i1, j2, k2)*dvol;  //Bx
            dvol = hiw[0][i] * wiw[1][j],
                B[1] += ebfield->B1(i2, j1, k2)*dvol;  //By
            dvol = hiw[0][i] * hiw[1][j],
                B[2] += ebfield->B2(i2, j2, k1)*dvol;  //Bz
          }
        }

        oldP[0] = ru(3, p);
        oldP[1] = ru(4, p);
        oldP[2] = ru(5, p);

        u_minus[0] = ru(3, p) + 0.5*dt*coupling*E[0];
        u_minus[1] = ru(4, p) + 0.5*dt*coupling*E[1];
        u_minus[2] = ru(5, p) + 0.5*dt*coupling*E[2];

        gamma_i = 1. / sqrt(1 + u_minus[0] * u_minus[0] + u_minus[1] * u_minus[1] + u_minus[2] * u_minus[2]);

        tee[0] = 0.5*dt*coupling*B[0] * gamma_i;
        tee[1] = 0.5*dt*coupling*B[1] * gamma_i;
        tee[2] = 0.5*dt*coupling*B[2] * gamma_i;

        u_prime[0] = u_minus[0] + (u_minus[1] * tee[2] - u_minus[2] * tee[1]);
        u_prime[1] = u_minus[1] + (u_minus[2] * tee[0] - u_minus[0] * tee[2]);
        u_prime[2] = u_minus[2] + (u_minus[0] * tee[1] - u_minus[1] * tee[0]);

        dummy = 1 / (1 + tee[0] * tee[0] + tee[1] * tee[1] + tee[2] * tee[2]);

        ess[0] = 2 * dummy*tee[0];
        ess[1] = 2 * dummy*tee[1];
        ess[2] = 2 * dummy*tee[2];

        u_plus[0] = u_minus[0] + u_prime[1] * ess[2] - u_prime[2] * ess[1];
        u_plus[1] = u_minus[1] + u_prime[2] * ess[0] - u_prime[0] * ess[2];
        u_plus[2] = u_minus[2] + u_prime[0] * ess[1] - u_prime[1] * ess[0];

        ru(3, p) = (u_plus[0] + 0.5*dt*coupling*E[0]);
        ru(4, p) = (u_plus[1] + 0.5*dt*coupling*E[1]);
        ru(5, p) = (u_plus[2] + 0.5*dt*coupling*E[2]);

        pn[0] = 0.5*(oldP[0] + ru(3, p));
        pn[1] = 0.5*(oldP[1] + ru(4, p));
        pn[2] = 0.5*(oldP[2] + ru(5, p));

        gamman = sqrt(1.0 + (pn[0] * pn[0] + pn[1] * pn[1] + pn[2] * pn[2]));

        vn[0] = pn[0] / gamman;
        vn[1] = pn[1] / gamman;
        vn[2] = pn[2] / gamman;

        fLorentz[0] = coupling*(E[0] + vn[1] * B[2] - vn[2] * B[1]);
        fLorentz[1] = coupling*(E[1] + vn[2] * B[0] - vn[0] * B[2]);
        fLorentz[2] = coupling*(E[2] + vn[0] * B[1] - vn[1] * B[0]);

        fLorentz2 = fLorentz[0] * fLorentz[0] + fLorentz[1] * fLorentz[1] + fLorentz[2] * fLorentz[2];

        vdotE2 = vn[0] * E[0] + vn[1] * E[1] + vn[2] * E[2];
        vdotE2 = vdotE2*vdotE2;

        ru(3, p) -= RRcoefficient*(gamman*gamman)*(fLorentz2 - vdotE2)*vn[0] * dt;
        ru(4, p) -= RRcoefficient*(gamman*gamman)*(fLorentz2 - vdotE2)*vn[1] * dt;
        ru(5, p) -= RRcoefficient*(gamman*gamman)*(fLorentz2 - vdotE2)*vn[2] * dt;

      }
      break;

    case 1:
      for (p = 0; p < Np; p++)
      {
        //gamma_i=1./sqrt(1+u0(p)*u0(p)+u1(p)*u1(p)+u2(p)*u2(p));
        for (c = 0; c < 3; c++)
        {
          xx[c] = ru(c, p);
          hiw[c][1] = wiw[c][1] = 1;
          hii[c] = wii[c] = 0;
        }
        for (c = 0; c < 1; c++)
        {
          rr = mygrid->dri[c] * (xx[c] - mygrid->rminloc[c]);
          rh = rr - 0.5;
          wii[c] = (int)floor(rr + 0.5); //whole integer int
          hii[c] = (int)floor(rr);     //half integer int
          rr -= wii[c];
          rh -= hii[c];
          rr2 = rr*rr;
          rh2 = rh*rh;

          wiw[c][1] = 0.75 - rr2;
          wiw[c][2] = 0.5*(0.25 + rr2 + rr);
          wiw[c][0] = 1. - wiw[c][1] - wiw[c][2];

          hiw[c][1] = 0.75 - rh2;
          hiw[c][2] = 0.5*(0.25 + rh2 + rh);
          hiw[c][0] = 1. - hiw[c][1] - hiw[c][2];
        }
        E[0] = E[1] = E[2] = B[0] = B[1] = B[2] = 0;

        k1 = k2 = j1 = j2 = 0;
        for (i = 0; i < 3; i++)
        {
          i1 = i + wii[0] - 1;
          i2 = i + hii[0] - 1;
          dvol = hiw[0][i],
              E[0] += ebfield->E0(i2, j1, k1)*dvol;  //Ex
          dvol = wiw[0][i],
              E[1] += ebfield->E1(i1, j2, k1)*dvol;  //Ey
          dvol = wiw[0][i],
              E[2] += ebfield->E2(i1, j1, k2)*dvol;  //Ez

          dvol = wiw[0][i],
              B[0] += ebfield->B0(i1, j2, k2)*dvol;  //Bx
          dvol = hiw[0][i],
              B[1] += ebfield->B1(i2, j1, k2)*dvol;  //By
          dvol = hiw[0][i],
              B[2] += ebfield->B2(i2, j2, k1)*dvol;  //Bz
        }

        oldP[0] = ru(3, p);
        oldP[1] = ru(4, p);
        oldP[2] = ru(5, p);

        u_minus[0] = ru(3, p) + 0.5*dt*coupling*E[0];
        u_minus[1] = ru(4, p) + 0.5*dt*coupling*E[1];
        u_minus[2] = ru(5, p) + 0.5*dt*coupling*E[2];

        gamma_i = 1. / sqrt(1 + u_minus[0] * u_minus[0] + u_minus[1] * u_minus[1] + u_minus[2] * u_minus[2]);

        tee[0] = 0.5*dt*coupling*B[0] * gamma_i;
        tee[1] = 0.5*dt*coupling*B[1] * gamma_i;
        tee[2] = 0.5*dt*coupling*B[2] * gamma_i;

        u_prime[0] = u_minus[0] + (u_minus[1] * tee[2] - u_minus[2] * tee[1]);
        u_prime[1] = u_minus[1] + (u_minus[2] * tee[0] - u_minus[0] * tee[2]);
        u_prime[2] = u_minus[2] + (u_minus[0] * tee[1] - u_minus[1] * tee[0]);

        dummy = 1 / (1 + tee[0] * tee[0] + tee[1] * tee[1] + tee[2] * tee[2]);

        ess[0] = 2 * dummy*tee[0];
        ess[1] = 2 * dummy*tee[1];
        ess[2] = 2 * dummy*tee[2];

        u_plus[0] = u_minus[0] + u_prime[1] * ess[2] - u_prime[2] * ess[1];
        u_plus[1] = u_minus[1] + u_prime[2] * ess[0] - u_prime[0] * ess[2];
        u_plus[2] = u_minus[2] + u_prime[0] * ess[1] - u_prime[1] * ess[0];

        ru(3, p) = (u_plus[0] + 0.5*dt*coupling*E[0]);
        ru(4, p) = (u_plus[1] + 0.5*dt*coupling*E[1]);
        ru(5, p) = (u_plus[2] + 0.5*dt*coupling*E[2]);

        pn[0] = 0.5*(oldP[0] + ru(3, p));
        pn[1] = 0.5*(oldP[1] + ru(4, p));
        pn[2] = 0.5*(oldP[2] + ru(5, p));

        gamman = sqrt(1.0 + (pn[0] * pn[0] + pn[1] * pn[1] + pn[2] * pn[2]));

        vn[0] = pn[0] / gamman;
        vn[1] = pn[1] / gamman;
        vn[2] = pn[2] / gamman;

        fLorentz[0] = coupling*(E[0] + vn[1] * B[2] - vn[2] * B[1]);
        fLorentz[1] = coupling*(E[1] + vn[2] * B[0] - vn[0] * B[2]);
        fLorentz[2] = coupling*(E[2] + vn[0] * B[1] - vn[1] * B[0]);

        fLorentz2 = fLorentz[0] * fLorentz[0] + fLorentz[1] * fLorentz[1] + fLorentz[2] * fLorentz[2];

        vdotE2 = vn[0] * E[0] + vn[1] * E[1] + vn[2] * E[2];
        vdotE2 = vdotE2*vdotE2;

        ru(3, p) -= RRcoefficient*(gamman*gamman)*(fLorentz2 - vdotE2)*vn[0] * dt;
        ru(4, p) -= RRcoefficient*(gamman*gamman)*(fLorentz2 - vdotE2)*vn[1] * dt;
        ru(5, p) -= RRcoefficient*(gamman*gamman)*(fLorentz2 - vdotE2)*vn[2] * dt;

      }
      break;
  }
}




void SPECIE::momentaStretchedAdvance(EM_FIELD *ebfield)
{
  energyExtremesFlag = false;
  if (mygrid->withParticles == NO)
    return;

  double dt, gamma_i;
  int p, c;  // particle_int, component_int
  int i, j, k, i1, j1, k1, i2, j2, k2;
  //int indexMaxQuadraticShape[]={1,4};
  int hii[3], wii[3];           // half integer index,   whole integer index
  double hiw[3][3], wiw[3][3];  // half integer weight,  whole integer weight
  double rr, rh, rr2, rh2;          // local coordinate to integer grid point and to half integer,     local coordinate squared
  double dvol, xx[3];           // tensor_product,       absolute particle position
  double E[3], B[3];
  double u_plus[3], u_minus[3], u_prime[3], tee[3], ess[3], dummy;
  double mycsi[3];

  dt = mygrid->dt;
  for (p = 0; p < Np; p++)
  {
    //gamma_i=1./sqrt(1+u0(p)*u0(p)+u1(p)*u1(p)+u2(p)*u2(p));
    for (c = 0; c < 3; c++)
    {
      xx[c] = ru(c, p);
      hiw[c][1] = wiw[c][1] = 1;
      hii[c] = wii[c] = 0;
    }
    for (c = 0; c < mygrid->getDimensionality(); c++)
    {
      mycsi[c] = mygrid->unStretchGrid(xx[c], c);
      rr = mygrid->dri[c] * (mycsi[c] - mygrid->csiminloc[c]);

      rh = rr - 0.5;
      wii[c] = (int)floor(rr + 0.5); //whole integer int
      hii[c] = (int)floor(rr);     //half integer int
      rr -= wii[c];
      rh -= hii[c];
      rr2 = rr*rr;
      rh2 = rh*rh;

      wiw[c][1] = 0.75 - rr2;
      wiw[c][2] = 0.5*(0.25 + rr2 + rr);
      wiw[c][0] = 1. - wiw[c][1] - wiw[c][2];

      hiw[c][1] = 0.75 - rh2;
      hiw[c][2] = 0.5*(0.25 + rh2 + rh);
      hiw[c][0] = 1. - hiw[c][1] - hiw[c][2];
    }
    E[0] = E[1] = E[2] = B[0] = B[1] = B[2] = 0;

    switch (mygrid->getDimensionality())
    {
      case 3:
        for (k = 0; k < 3; k++)
        {
          k1 = k + wii[2] - 1;
          k2 = k + hii[2] - 1;
          for (j = 0; j < 3; j++)
          {
            j1 = j + wii[1] - 1;
            j2 = j + hii[1] - 1;
            for (i = 0; i < 3; i++)
            {
              i1 = i + wii[0] - 1;
              i2 = i + hii[0] - 1;
              dvol = hiw[0][i] * wiw[1][j] * wiw[2][k],
                  E[0] += ebfield->E0(i2, j1, k1)*dvol;  //Ex
              dvol = wiw[0][i] * hiw[1][j] * wiw[2][k],
                  E[1] += ebfield->E1(i1, j2, k1)*dvol;  //Ey
              dvol = wiw[0][i] * wiw[1][j] * hiw[2][k],
                  E[2] += ebfield->E2(i1, j1, k2)*dvol;  //Ez

              dvol = wiw[0][i] * hiw[1][j] * hiw[2][k],
                  B[0] += ebfield->B0(i1, j2, k2)*dvol;  //Bx
              dvol = hiw[0][i] * wiw[1][j] * hiw[2][k],
                  B[1] += ebfield->B1(i2, j1, k2)*dvol;  //By
              dvol = hiw[0][i] * hiw[1][j] * wiw[2][k],
                  B[2] += ebfield->B2(i2, j2, k1)*dvol;  //Bz
            }
          }
        }
        break;

      case 2:
        k1 = k2 = 0;
        for (j = 0; j < 3; j++)
        {
          j1 = j + wii[1] - 1;
          j2 = j + hii[1] - 1;
          for (i = 0; i < 3; i++)
          {
            i1 = i + wii[0] - 1;
            i2 = i + hii[0] - 1;
            dvol = hiw[0][i] * wiw[1][j],
                E[0] += ebfield->E0(i2, j1, k1)*dvol;  //Ex
            dvol = wiw[0][i] * hiw[1][j],
                E[1] += ebfield->E1(i1, j2, k1)*dvol;  //Ey
            dvol = wiw[0][i] * wiw[1][j],
                E[2] += ebfield->E2(i1, j1, k2)*dvol;  //Ez

            dvol = wiw[0][i] * hiw[1][j],
                B[0] += ebfield->B0(i1, j2, k2)*dvol;  //Bx
            dvol = hiw[0][i] * wiw[1][j],
                B[1] += ebfield->B1(i2, j1, k2)*dvol;  //By
            dvol = hiw[0][i] * hiw[1][j],
                B[2] += ebfield->B2(i2, j2, k1)*dvol;  //Bz
          }
        }
        break;

      case 1:
        k1 = k2 = j1 = j2 = 0;
        for (i = 0; i < 3; i++)
        {
          i1 = i + wii[0] - 1;
          i2 = i + hii[0] - 1;
          dvol = hiw[0][i],
              E[0] += ebfield->E0(i2, j1, k1)*dvol;  //Ex
          dvol = wiw[0][i],
              E[1] += ebfield->E1(i1, j2, k1)*dvol;  //Ey
          dvol = wiw[0][i],
              E[2] += ebfield->E2(i1, j1, k2)*dvol;  //Ez

          dvol = wiw[0][i],
              B[0] += ebfield->B0(i1, j2, k2)*dvol;  //Bx
          dvol = hiw[0][i],
              B[1] += ebfield->B1(i2, j1, k2)*dvol;  //By
          dvol = hiw[0][i],
              B[2] += ebfield->B2(i2, j2, k1)*dvol;  //Bz
        }
        break;
    }

    u_minus[0] = ru(3, p) + 0.5*dt*coupling*E[0];
    u_minus[1] = ru(4, p) + 0.5*dt*coupling*E[1];
    u_minus[2] = ru(5, p) + 0.5*dt*coupling*E[2];

    gamma_i = 1. / sqrt(1 + u_minus[0] * u_minus[0] + u_minus[1] * u_minus[1] + u_minus[2] * u_minus[2]);

    tee[0] = 0.5*dt*coupling*B[0] * gamma_i;
    tee[1] = 0.5*dt*coupling*B[1] * gamma_i;
    tee[2] = 0.5*dt*coupling*B[2] * gamma_i;

    u_prime[0] = u_minus[0] + (u_minus[1] * tee[2] - u_minus[2] * tee[1]);
    u_prime[1] = u_minus[1] + (u_minus[2] * tee[0] - u_minus[0] * tee[2]);
    u_prime[2] = u_minus[2] + (u_minus[0] * tee[1] - u_minus[1] * tee[0]);

    dummy = 1 / (1 + tee[0] * tee[0] + tee[1] * tee[1] + tee[2] * tee[2]);

    ess[0] = 2 * dummy*tee[0];
    ess[1] = 2 * dummy*tee[1];
    ess[2] = 2 * dummy*tee[2];

    u_plus[0] = u_minus[0] + u_prime[1] * ess[2] - u_prime[2] * ess[1];
    u_plus[1] = u_minus[1] + u_prime[2] * ess[0] - u_prime[0] * ess[2];
    u_plus[2] = u_minus[2] + u_prime[0] * ess[1] - u_prime[1] * ess[0];

    ru(3, p) = (u_plus[0] + 0.5*dt*coupling*E[0]);
    ru(4, p) = (u_plus[1] + 0.5*dt*coupling*E[1]);
    ru(5, p) = (u_plus[2] + 0.5*dt*coupling*E[2]);

  }
}

void SPECIE::current_deposition(CURRENT *current)
{
  if (mygrid->withParticles == NO)
    return;
  if (mygrid->withCurrent == NO)
    return;

  double dt, gamma_i;
  int p, c;  // particle_int, component_int
  int i2, j2, k2, ti, tj, tk;
  //int indexMaxQuadraticShape[]={1,4};
  int ii1[3], ii2[3], di;           // half integer index,   whole integer index
  double w1[3][5], w2[3][5];  // half integer weight,  whole integer weight
  double r1, r2, r12, r22;          // local coordinate to integer grid point and to half integer,     local coordinate squared
  double xx1[3], xx2[3];           // tensor_product,       absolute particle position
  double s0x, s0y, s0z, dsx, dsy, dsz;
  double J[3][5][5][5], W[3][5][5][5], norm, vz, vy;

  dt = mygrid->dt;

  if (mygrid->getDimensionality() == 3)
    for (p = 0; p < Np; p++)
    {
      memset((void*)J, 0, 3 * 5 * 5 * 5 * sizeof(double));
      memset((void*)W, 0, 3 * 5 * 5 * 5 * sizeof(double));
      gamma_i = 1. / sqrt(1 + u0(p)*u0(p) + u1(p)*u1(p) + u2(p)*u2(p));

      xx1[0] = ru(0, p);
      ru(0, p) += dt*gamma_i*u0(p);
      xx2[0] = ru(0, p);

      xx1[1] = ru(1, p);
      ru(1, p) += dt*gamma_i*u1(p);
      xx2[1] = ru(1, p);

      xx1[2] = ru(2, p);
      ru(2, p) += dt*gamma_i*u2(p);
      xx2[2] = ru(2, p);

      for (c = 0; c < mygrid->getDimensionality(); c++)
      {
        r1 = mygrid->dri[c] * (xx1[c] - mygrid->rminloc[c]);
        r2 = mygrid->dri[c] * (xx2[c] - mygrid->rminloc[c]);
        ii1[c] = (int)floor(r1 + 0.5);
        ii2[c] = (int)floor(r2 + 0.5);
        r1 -= ii1[c];
        r2 -= ii2[c];
        r12 = r1*r1;
        r22 = r2*r2;
        di = ii2[c] - ii1[c];

        w1[c][4] = 0;
        w1[c][3] = 0.5*(0.25 + r12 + r1);
        w1[c][2] = 0.75 - r12;
        w1[c][1] = 1. - w1[c][3] - w1[c][2];
        w1[c][0] = 0;

        w2[c][(4 + di) % 5] = 0.;
        w2[c][3 + di] = 0.5*(0.25 + r22 + r2);
        w2[c][2 + di] = 0.75 - r22;
        w2[c][1 + di] = 1. - w2[c][2 + di] - w2[c][3 + di];
        w2[c][(0 + di) % 5] = 0;

      }
      norm = 1.;
      for (tk = 0; tk < 5; tk++)
      {
        k2 = tk + ii1[2] - 2;
        for (tj = 0; tj < 5; tj++)
        {
          j2 = tj + ii1[1] - 2;
          for (ti = 0; ti < 5; ti++)
          {
            i2 = ti + ii1[0] - 2;

            s0x = w1[0][ti];
            s0y = w1[1][tj];
            s0z = w1[2][tk];
            dsx = -w1[0][ti] + w2[0][ti];
            dsy = -w1[1][tj] + w2[1][tj];
            dsz = -w1[2][tk] + w2[2][tk];

            W[0][0][tj][tk] += norm*dsx*(s0y*s0z + 0.5*dsy*s0z + 0.5*s0y*dsz + UN_TERZO*dsy*dsz) / dt;
            W[1][ti][0][tk] += norm*dsy*(s0z*s0x + 0.5*dsz*s0x + 0.5*s0z*dsx + UN_TERZO*dsz*dsx) / dt;
            W[2][ti][tj][0] += norm*dsz*(s0x*s0y + 0.5*dsx*s0y + 0.5*s0x*dsy + UN_TERZO*dsx*dsy) / dt;

            J[0][ti][tj][tk] = -mygrid->dr[0] * W[0][0][tj][tk];
            J[1][ti][tj][tk] = -mygrid->dr[1] * W[1][ti][0][tk];
            J[2][ti][tj][tk] = -mygrid->dr[2] * W[2][ti][tj][0];
            current->Jx(i2, j2, k2) += chargeSign*w(p)*J[0][ti][tj][tk];
            current->Jy(i2, j2, k2) += chargeSign*w(p)*J[1][ti][tj][tk];
            current->Jz(i2, j2, k2) += chargeSign*w(p)*J[2][ti][tj][tk];

          }
        }
      }
    }
  if (mygrid->getDimensionality() == 2)
    for (p = 0; p < Np; p++)
    {
      memset((void*)J, 0, 3 * 5 * 5 * 5 * sizeof(double));
      memset((void*)W, 0, 3 * 5 * 5 * 5 * sizeof(double));
      gamma_i = 1. / sqrt(1 + u0(p)*u0(p) + u1(p)*u1(p) + u2(p)*u2(p));
      vz = gamma_i*u2(p);
      //ru(2, p) += dt*vz;

      xx1[0] = ru(0, p);
      ru(0, p) += dt*gamma_i*u0(p);
      xx2[0] = ru(0, p);

      xx1[1] = ru(1, p);
      ru(1, p) += dt*gamma_i*u1(p);
      xx2[1] = ru(1, p);

      for (c = 0; c < mygrid->getDimensionality(); c++)
      {
        r1 = mygrid->dri[c] * (xx1[c] - mygrid->rminloc[c]);
        r2 = mygrid->dri[c] * (xx2[c] - mygrid->rminloc[c]);
        ii1[c] = (int)floor(r1 + 0.5);
        ii2[c] = (int)floor(r2 + 0.5);
        //ii1[c]=rint(r1);
        //ii2[c]=rint(r2);
        //ii1[c]=(int)(r1+0.5);
        //ii2[c]=(int)(r2+0.5);
        r1 -= ii1[c];
        r2 -= ii2[c];
        r12 = r1*r1;
        r22 = r2*r2;
        di = ii2[c] - ii1[c];

        w1[c][4] = 0;
        w1[c][3] = 0.5*(0.25 + r12 + r1);
        w1[c][2] = 0.75 - r12;
        w1[c][1] = 1. - w1[c][3] - w1[c][2];
        w1[c][0] = 0;

        w2[c][(4 + di) % 5] = 0.;
        w2[c][3 + di] = 0.5*(0.25 + r22 + r2);
        w2[c][2 + di] = 0.75 - r22;
        w2[c][1 + di] = 1. - w2[c][2 + di] - w2[c][3 + di];
        w2[c][(0 + di) % 5] = 0;

      }
      norm = 1.;

      tk = k2 = 0;//tk+ii1[2]-2;
      for (tj = 0; tj < 5; tj++)
      {
        j2 = tj + ii1[1] - 2;
        for (ti = 0; ti < 5; ti++)
        {
          i2 = ti + ii1[0] - 2;

          s0x = w1[0][ti];
          s0y = w1[1][tj];

          dsx = -w1[0][ti] + w2[0][ti];
          dsy = -w1[1][tj] + w2[1][tj];


          W[0][0][tj][tk] += norm*dsx*(s0y + 0.5*dsy) / dt;
          W[1][ti][0][tk] += norm*dsy*(s0x + 0.5*dsx) / dt;
          W[2][ti][tj][0] = norm*vz*(s0x*s0y + 0.5*dsx*s0y + 0.5*s0x*dsy + UN_TERZO*dsx*dsy);

          J[0][ti][tj][tk] = -mygrid->dr[0] * W[0][0][tj][tk];
          J[1][ti][tj][tk] = -mygrid->dr[1] * W[1][ti][0][tk];
          J[2][ti][tj][tk] = W[2][ti][tj][0];
          current->Jx(i2, j2, k2) += chargeSign*w(p)*J[0][ti][tj][tk];
          current->Jy(i2, j2, k2) += chargeSign*w(p)*J[1][ti][tj][tk];
          current->Jz(i2, j2, k2) += chargeSign*w(p)*J[2][ti][tj][tk];

        }
      }

    }
  if (mygrid->getDimensionality() == 1)
    for (p = 0; p < Np; p++)
    {
      memset((void*)J, 0, 3 * 5 * 5 * 5 * sizeof(double));
      memset((void*)W, 0, 3 * 5 * 5 * 5 * sizeof(double));
      gamma_i = 1. / sqrt(1 + u0(p)*u0(p) + u1(p)*u1(p) + u2(p)*u2(p));
      vy = gamma_i*u1(p);
      vz = gamma_i*u2(p);
      //			ru(1, p) += dt*vy;
      //			ru(2, p) += dt*vz;
      for (c = 0; c < mygrid->getDimensionality(); c++)
      {
        xx1[c] = ru(c, p);
        ru(c, p) += dt*gamma_i*u0(p);
        xx2[c] = ru(c, p);
      }
      for (c = 0; c < mygrid->getDimensionality(); c++)
      {
        r1 = mygrid->dri[c] * (xx1[c] - mygrid->rminloc[c]);
        r2 = mygrid->dri[c] * (xx2[c] - mygrid->rminloc[c]);
        ii1[c] = (int)floor(r1 + 0.5);
        ii2[c] = (int)floor(r2 + 0.5);
        r1 -= ii1[c];
        r2 -= ii2[c];
        r12 = r1*r1;
        r22 = r2*r2;
        di = ii2[c] - ii1[c];

        w1[c][4] = 0;
        w1[c][3] = 0.5*(0.25 + r12 + r1);
        w1[c][2] = 0.75 - r12;
        w1[c][1] = 1. - w1[c][3] - w1[c][2];
        w1[c][0] = 0;

        w2[c][(4 + di) % 5] = 0.;
        w2[c][3 + di] = 0.5*(0.25 + r22 + r2);
        w2[c][2 + di] = 0.75 - r22;
        w2[c][1 + di] = 1. - w2[c][2 + di] - w2[c][3 + di];
        w2[c][(0 + di) % 5] = 0;

      }
      norm = 1.;

      tj = j2 = 0;//tk+ii1[2]-2;
      tk = k2 = 0;//tk+ii1[2]-2;

      for (ti = 0; ti < 5; ti++)
      {
        i2 = ti + ii1[0] - 2;

        s0x = w1[0][ti];
        dsx = w2[0][ti] - w1[0][ti];

        W[0][0][tj][tk] += norm*dsx / mygrid->dt;
        W[1][ti][0][tk] = norm*vy*(s0x + 0.5*dsx);
        W[2][ti][tj][0] = norm*vz*(s0x + 0.5*dsx);

        J[0][ti][tj][tk] = -mygrid->dr[0] * W[0][0][tj][tk];
        J[1][ti][tj][tk] = W[1][ti][0][tk];
        J[2][ti][tj][tk] = W[2][ti][tj][0];
        current->Jx(i2, j2, k2) += chargeSign*w(p)*J[0][ti][tj][tk];
        current->Jy(i2, j2, k2) += chargeSign*w(p)*J[1][ti][tj][tk];
        current->Jz(i2, j2, k2) += chargeSign*w(p)*J[2][ti][tj][tk];

      }


    }

}



void SPECIE::add_momenta(double uxin, double uyin, double uzin) //Aggiunge semplicemente un drift a tutta la specie
{

  if (!allocated){
    return;
  }
  int p;
  for (p = 0; p < Np; p++)
  {
    u0(p) += uxin;
    u1(p) += uyin;
    u2(p) += uzin;
  }
}

// WATERBAG        : [P0] -P0< px < + P0, -P0< py < + P0, -P0< pz < + P0 uniforme
// WATERBAG_3TEMP  : [P0_X,P0_Y,P0_Z] -P0_X< px < + P0_X, -P0_Y< py < + P0_Y, -P0_Z< pz < + P0_Z uniforme
// UNIF_SPHERE     : [P0] -P0 < p < +P0 uniforme
// SUPERGAUSSIAN   : [P0, ALPHA] f(p) = C*exp(-abs(p/P0)^(ALPHA))
// MAXWELL         : [Ta] Maxwell alla Macchi
// JUTTNER         : [a] f(p) = C*exp(-a*gamma(p)) [DA RISCRIVERE]

void SPECIE::computeLorentzMatrix(double ux, double uy, double uz, double *matr){
  double gm = sqrt(1.0 + ux*ux + uy*uy + uz*uz);
  double bx = -ux / gm;
  double by = -uy / gm;
  double bz = -uz / gm;
  double b = sqrt(bx*bx + by*by + bz*bz);

  matr[0 * 4 + 0] = gm;     matr[0 * 4 + 1] = -gm*bx;                 matr[0 * 4 + 2] = -gm*by;                 matr[0 * 4 + 3] = -gm*bz;
  matr[1 * 4 + 0] = -gm*bx; matr[1 * 4 + 1] = 1.0 + (gm - 1.0)*bx*bx / b / b; matr[1 * 4 + 2] = (gm - 1.0)*by*bx / b / b;     matr[1 * 4 + 3] = (gm - 1.0)*bz*bx / b / b;
  matr[2 * 4 + 0] = -gm*by; matr[2 * 4 + 1] = (gm - 1.0)*bx*by / b / b;     matr[2 * 4 + 2] = 1.0 + (gm - 1.0)*by*by / b / b; matr[2 * 4 + 3] = (gm - 1.0)*bz*by / b / b;
  matr[3 * 4 + 0] = -gm*bz; matr[3 * 4 + 1] = (gm - 1.0)*bx*bz / b / b;     matr[3 * 4 + 2] = (gm - 1.0)*by*bz / b / b;     matr[3 * 4 + 3] = 1.0 + (gm - 1.0)*bz*bz / b / b;

}

void SPECIE::callWaterbag(gsl_rng* ext_rng, double p0_x, double p0_y, double p0_z, double uxin, double uyin, double uzin){

  if (uxin*uxin + uyin*uyin + uzin*uzin < _VERY_SMALL_MOMENTUM*_VERY_SMALL_MOMENTUM){

    for (int p = 0; p < Np; p++)
    {
      u0(p) = uxin + p0_x*gsl_ran_flat(ext_rng, -1.0, 1.0);
      u1(p) = uyin + p0_y*gsl_ran_flat(ext_rng, -1.0, 1.0);
      u2(p) = uzin + p0_z*gsl_ran_flat(ext_rng, -1.0, 1.0);
    }
  }
  else{
    double L[16];
    computeLorentzMatrix(uxin, uyin, uzin, L);
    double Ett, u0t, u1t, u2t;
    for (int p = 0; p < Np; p++)
    {
      u0(p) = p0_x*gsl_ran_flat(ext_rng, -1.0, 1.0);
      u1(p) = p0_y*gsl_ran_flat(ext_rng, -1.0, 1.0);
      u2(p) = p0_z*gsl_ran_flat(ext_rng, -1.0, 1.0);
      Ett = sqrt(1.0 + u0(p)*u0(p) + u1(p)*u1(p) + u2(p)*u2(p));

      u0t = L[1 * 4 + 0] * Ett + L[1 * 4 + 1] * u0(p) + L[1 * 4 + 2] * u1(p) + L[1 * 4 + 3] * u2(p);
      u1t = L[2 * 4 + 0] * Ett + L[2 * 4 + 1] * u0(p) + L[2 * 4 + 2] * u1(p) + L[2 * 4 + 3] * u2(p);
      u2t = L[3 * 4 + 0] * Ett + L[3 * 4 + 1] * u0(p) + L[3 * 4 + 2] * u1(p) + L[3 * 4 + 3] * u2(p);

      u0(p) = u0t;
      u1(p) = u1t;
      u2(p) = u2t;
    }
  }

}

void SPECIE::callUnifSphere(gsl_rng* ext_rng, double p0, double uxin, double uyin, double uzin){
  double pmod;
  double phi;
  double cos_theta, sin_theta;

  if (uxin*uxin + uyin*uyin + uzin*uzin < _VERY_SMALL_MOMENTUM*_VERY_SMALL_MOMENTUM){
    for (int p = 0; p < Np; p++)
    {
      pmod = pow(gsl_ran_flat(ext_rng, 0.0, 1.0), 1. / 3.);
      phi = gsl_ran_flat(ext_rng, 0.0, 2.0*M_PI);
      cos_theta = gsl_ran_flat(ext_rng, -1.0, 1.0);
      sin_theta = sqrt(1.0 - cos_theta*cos_theta);
      u0(p) = uxin + p0*pmod*sin_theta*cos(phi);
      u1(p) = uyin + p0*pmod*sin_theta*sin(phi);
      u2(p) = uzin + p0*pmod*cos_theta;
    }
  }
  else{
    double L[16];
    computeLorentzMatrix(uxin, uyin, uzin, L);
    double Ett, u0t, u1t, u2t;
    for (int p = 0; p < Np; p++)
    {
      pmod = pow(gsl_ran_flat(ext_rng, 0.0, 1.0), 1. / 3.);
      phi = gsl_ran_flat(ext_rng, 0.0, 2.0*M_PI);
      cos_theta = gsl_ran_flat(ext_rng, -1.0, 1.0);
      sin_theta = sqrt(1.0 - cos_theta*cos_theta);
      u0(p) = p0*pmod*sin_theta*cos(phi);
      u1(p) = p0*pmod*sin_theta*sin(phi);
      u2(p) = p0*pmod*cos_theta;
      Ett = sqrt(1.0 + u0(p)*u0(p) + u1(p)*u1(p) + u2(p)*u2(p));

      u0t = L[1 * 4 + 0] * Ett + L[1 * 4 + 1] * u0(p) + L[1 * 4 + 2] * u1(p) + L[1 * 4 + 3] * u2(p);
      u1t = L[2 * 4 + 0] * Ett + L[2 * 4 + 1] * u0(p) + L[2 * 4 + 2] * u1(p) + L[2 * 4 + 3] * u2(p);
      u2t = L[3 * 4 + 0] * Ett + L[3 * 4 + 1] * u0(p) + L[3 * 4 + 2] * u1(p) + L[3 * 4 + 3] * u2(p);

      u0(p) = u0t;
      u1(p) = u1t;
      u2(p) = u2t;
    }
  }
}

void SPECIE::callSupergaussian(gsl_rng* ext_rng, double p0, double alpha, double uxin, double uyin, double uzin){
  if (uxin*uxin + uyin*uyin + uzin*uzin < _VERY_SMALL_MOMENTUM*_VERY_SMALL_MOMENTUM){
    for (int p = 0; p < Np; p++)
    {
      u0(p) = uxin + gsl_ran_exppow(ext_rng, p0, alpha);
      u1(p) = uyin + gsl_ran_exppow(ext_rng, p0, alpha);
      u2(p) = uzin + gsl_ran_exppow(ext_rng, p0, alpha);
    }
  }
  else{
    double L[16];
    computeLorentzMatrix(uxin, uyin, uzin, L);
    double Ett, u0t, u1t, u2t;
    for (int p = 0; p < Np; p++)
    {
      u0(p) = gsl_ran_exppow(ext_rng, p0, alpha);
      u1(p) = gsl_ran_exppow(ext_rng, p0, alpha);
      u2(p) = gsl_ran_exppow(ext_rng, p0, alpha);
      Ett = sqrt(1.0 + u0(p)*u0(p) + u1(p)*u1(p) + u2(p)*u2(p));

      u0t = L[1 * 4 + 0] * Ett + L[1 * 4 + 1] * u0(p) + L[1 * 4 + 2] * u1(p) + L[1 * 4 + 3] * u2(p);
      u1t = L[2 * 4 + 0] * Ett + L[2 * 4 + 1] * u0(p) + L[2 * 4 + 2] * u1(p) + L[2 * 4 + 3] * u2(p);
      u2t = L[3 * 4 + 0] * Ett + L[3 * 4 + 1] * u0(p) + L[3 * 4 + 2] * u1(p) + L[3 * 4 + 3] * u2(p);

      u0(p) = u0t;
      u1(p) = u1t;
      u2(p) = u2t;
    }

  }

}

void SPECIE::callMaxwell(gsl_rng* ext_rng, double Ta, double uxin, double uyin, double uzin){
  double ptot;
  double temp;
  double phi;
  double cos_theta, sin_theta;
  double theta;
  if (uxin*uxin + uyin*uyin + uzin*uzin < _VERY_SMALL_MOMENTUM*_VERY_SMALL_MOMENTUM){
    for (int p = 0; p < Np; p++)
    {

      temp = gsl_ran_exponential(ext_rng, Ta);
      ptot = sqrt((temp + 1)*(temp + 1) - 1 * 1);
      phi = gsl_ran_flat(ext_rng, 0.0, 2.0*M_PI);
      cos_theta = gsl_ran_flat(ext_rng, -1.0, 1.0);
      sin_theta = sqrt(1.0 - cos_theta*cos_theta);
      u0(p) = uxin + ptot*sin_theta*cos(phi);
      u1(p) = uyin + ptot*sin_theta*sin(phi);
      u2(p) = uzin + ptot*cos_theta;
    }
  }
  else{
    double L[16];
    computeLorentzMatrix(uxin, uyin, uzin, L);
    double Ett, u0t, u1t, u2t;
    for (int p = 0; p < Np; p++)
    {

      temp = gsl_ran_exponential(ext_rng, Ta);
      ptot = sqrt((temp + 1)*(temp + 1) - 1 * 1);
      phi = gsl_ran_flat(ext_rng, 0.0, 2.0*M_PI);
      cos_theta = gsl_ran_flat(ext_rng, -1.0, 1.0);
      sin_theta = sqrt(1.0 - cos_theta*cos_theta);
      u0(p) = ptot*sin_theta*cos(phi);
      u1(p) = ptot*sin_theta*sin(phi);
      u2(p) = ptot*cos_theta;
      Ett = sqrt(1.0 + u0(p)*u0(p) + u1(p)*u1(p) + u2(p)*u2(p));

      u0t = L[1 * 4 + 0] * Ett + L[1 * 4 + 1] * u0(p) + L[1 * 4 + 2] * u1(p) + L[1 * 4 + 3] * u2(p);
      u1t = L[2 * 4 + 0] * Ett + L[2 * 4 + 1] * u0(p) + L[2 * 4 + 2] * u1(p) + L[2 * 4 + 3] * u2(p);
      u2t = L[3 * 4 + 0] * Ett + L[3 * 4 + 1] * u0(p) + L[3 * 4 + 2] * u1(p) + L[3 * 4 + 3] * u2(p);

      u0(p) = u0t;
      u1(p) = u1t;
      u2(p) = u2t;
    }
  }
}
void SPECIE::callJuttner(gsl_rng* ext_rng, double Ta, double uxin, double uyin, double uzin){
  //DA DEFINIRE
}
double densityFunctionMaxwell(double px, double alpha, double temp){
  return exp(-(sqrt(alpha*alpha + px*px) - alpha) / temp);
}
void SPECIE::callSpecial(gsl_rng* ext_rng, double Ta){
  //double ptot, temp, cos_theta, sin_theta, segno, phi;
  double alpha;
  double auxPX, auxDF;
  for (int p = 0; p < Np; p++)
  {
    double uperp2 = u1(p)*u1(p) - u2(p)*u2(p);
    alpha = sqrt(1 + uperp2);
    while (1){
      auxDF = gsl_ran_flat(ext_rng, 0.0, 1.0);
      auxPX = gsl_ran_flat(ext_rng, -50 * sqrt(Ta), 50 * sqrt(Ta));
      if (densityFunctionMaxwell(auxPX, alpha, Ta) > auxDF)
        break;
    }
    u0(p) = auxPX;
  }
}

void SPECIE::add_momenta(gsl_rng* ext_rng, double uxin, double uyin, double uzin, tempDistrib distribution)
{
  if (mygrid->withParticles == NO)
    return;

  if (!allocated){
#ifndef NO_ALLOCATION
    std::cout << "Warning: species " << name << " is not allocated !" << std::endl;
#endif
    return;
  }

  if (!distribution.init){
    std::cout << "Warning: distribution function is not initialized !" << std::endl;
    return;
  }

  switch (distribution.type)
  {
    //TRASFORMARE LE ISTRUZIONI NEI DIVERSI CASI IN CHIAMATE A FUNZIONI PRIVATE
    case WATERBAG:
      callWaterbag(ext_rng,
                   distribution.p0,
                   distribution.p0,
                   distribution.p0,
                   uxin, uyin, uzin);
      break;

    case WATERBAG_3TEMP:

      callWaterbag(ext_rng,
                   distribution.p0_x,
                   distribution.p0_y,
                   distribution.p0_z,
                   uxin, uyin, uzin);
      break;

    case UNIF_SPHERE:
      callUnifSphere(ext_rng,
                     distribution.p0,
                     uxin, uyin, uzin);
      break;

    case SUPERGAUSSIAN:
      callSupergaussian(ext_rng,
                        distribution.p0,
                        distribution.alpha,
                        uxin, uyin, uzin);
      break;

    case MAXWELL:
      callMaxwell(ext_rng,
                  distribution.temp,
                  uxin, uyin, uzin);
      break;

    case JUTTNER:
      callJuttner(ext_rng,
                  distribution.a,
                  uxin, uyin, uzin);
      break;

    case SPECIAL:
      callSpecial(ext_rng,
                  distribution.a);
      break;

    default:
      break;


  }

}


void SPECIE::current_deposition_standard(CURRENT *current)
{
  if (mygrid->withParticles == NO)
    return;
  if (mygrid->withCurrent == NO)
    return;
  if (mygrid->isStretched()){
    SPECIE::currentStretchedDepositionStandard(current);
    return;
  }

  double dt, gamma_i;
  int p, c;  // particle_int, component_int
  int i, j, k, i1, j1, k1, i2, j2, k2;
  //int indexMaxQuadraticShape[]={1,4};
  int hii[3], wii[3];           // half integer index,   whole integer index
  double hiw[3][3], wiw[3][3];  // half integer weight,  whole integer weight
  double rr, rh, rr2, rh2;          // local coordinate to integer grid point and to half integer,     local coordinate squared
  double dvol, xx[3], vv[3];           // tensor_product,       absolute particle position

  dt = mygrid->dt;
  double* myCurrent = current->getDataPointer();
  int N_grid[3];
  current->writeN_grid(N_grid);
  int edge=mygrid->getEdge();
  if (!(mygrid->withCurrent == YES && (!isTestSpecies)))
  {
    for (p = 0; p < Np; p++)
    {
      double ux, uy, uz;
      ux = pData[pIndex( 3, p, Ncomp, Np)];
      uy = pData[pIndex( 4, p, Ncomp, Np)];
      uz = pData[pIndex( 5, p, Ncomp, Np)];
      gamma_i = 1. / sqrt(1 + ux*ux + uy*uy + uz*uz);

      for (c = 0; c < mygrid->getDimensionality(); c++)
      {
        vv[c] = gamma_i*pData[pIndex( c + 3, p, Ncomp, Np)];
        pData[pIndex( c, p, Ncomp, Np)] += dt*vv[c];

      }
    }
    return;
  }

  switch (mygrid->getDimensionality())
  {
    case 3:
      for (p = 0; p < Np; p++)
      {
        double ux, uy, uz;
        ux = pData[pIndex( 3, p, Ncomp, Np)];
        uy = pData[pIndex( 4, p, Ncomp, Np)];
        uz = pData[pIndex( 5, p, Ncomp, Np)];
        gamma_i = 1. / sqrt(1 + ux*ux + uy*uy + uz*uz);

        for (c = 0; c < 3; c++)
        {
          vv[c] = gamma_i*pData[pIndex( c + 3, p, Ncomp, Np)];
          hiw[c][1] = wiw[c][1] = 1;
          hiw[c][0] = wiw[c][0] = 0;
          hiw[c][2] = wiw[c][2] = 0;
          hii[c] = wii[c] = 0;
        }
        for (c = 0; c < 3; c++)
        {
          xx[c] = pData[pIndex( c, p, Ncomp, Np)] + 0.5*dt*vv[c];
          pData[pIndex( c, p, Ncomp, Np)] += dt*vv[c];

          rr = mygrid->dri[c] * (xx[c] - mygrid->rminloc[c]);
          rh = rr - 0.5;

          wii[c] = (int)floor(rr + 0.5); //whole integer int
          hii[c] = (int)floor(rr);     //half integer int
          rr -= wii[c];
          rh -= hii[c];
          rr2 = rr*rr;
          rh2 = rh*rh;

          wiw[c][1] = 0.75 - rr2;
          wiw[c][2] = 0.5*(0.25 + rr2 + rr);
          wiw[c][0] = 1. - wiw[c][1] - wiw[c][2];

          hiw[c][1] = 0.75 - rh2;
          hiw[c][2] = 0.5*(0.25 + rh2 + rh);
          hiw[c][0] = 1. - hiw[c][1] - hiw[c][2];
        }


        for (k = 0; k < 3; k++)
        {
          k1 = k + wii[2] - 1;
          k2 = k + hii[2] - 1;
          for (j = 0; j < 3; j++)
          {
            j1 = j + wii[1] - 1;
            j2 = j + hii[1] - 1;
            for (i = 0; i < 3; i++)
            {
              i1 = i + wii[0] - 1;
              i2 = i + hii[0] - 1;
#ifndef OLD_ACCESS

              double weight = pData[pIndex( 6, p, Ncomp, Np)];
              double *JX, *JY, *JZ;
              JX = &myCurrent[my_indice(edge,1, 1, 0, i2, j1, k1, N_grid[0], N_grid[1], N_grid[2], current->Ncomp)];
              JY = &myCurrent[my_indice(edge,1, 1, 1, i1, j2, k1, N_grid[0], N_grid[1], N_grid[2], current->Ncomp)];
              JZ = &myCurrent[my_indice(edge,1, 1, 2, i1, j1, k2, N_grid[0], N_grid[1], N_grid[2], current->Ncomp)];

              dvol = hiw[0][i] * wiw[1][j] * wiw[2][k];
              *JX += weight*dvol*vv[0] * chargeSign;
              dvol = wiw[0][i] * hiw[1][j] * wiw[2][k];
              *JY += weight*dvol*vv[1] * chargeSign;
              dvol = wiw[0][i] * wiw[1][j] * hiw[2][k];
              *JZ += weight*dvol*vv[2] * chargeSign;
#else
                dvol = hiw[0][i] * wiw[1][j] * wiw[2][k],
                    current->Jx(i2, j1, k1) += myweight*dvol*vv[0] * chargeSign;
                dvol = wiw[0][i] * hiw[1][j] * wiw[2][k],
                    current->Jy(i1, j2, k1) += myweight*dvol*vv[1] * chargeSign;
                dvol = wiw[0][i] * wiw[1][j] * hiw[2][k],
                    current->Jz(i1, j1, k2) += myweight*dvol*vv[2] * chargeSign;
 #endif


            }
          }
        }
      }
      break;

    case 2:
      for (p = 0; p < Np; p++)
      {
        double ux, uy, uz;
        ux = pData[pIndex( 3, p, Ncomp, Np)];
        uy = pData[pIndex( 4, p, Ncomp, Np)];
        uz = pData[pIndex( 5, p, Ncomp, Np)];


        gamma_i = 1. / sqrt(1 + ux*ux + uy*uy + uz*uz);
        for (c = 0; c < 3; c++)
        {
          vv[c] = gamma_i*pData[pIndex( c + 3, p, Ncomp, Np)];
          hiw[c][1] = wiw[c][1] = 1;
          hiw[c][0] = wiw[c][0] = 0;
          hiw[c][2] = wiw[c][2] = 0;
          hii[c] = wii[c] = 0;
        }
        for (c = 0; c < 2; c++)
        {
          xx[c] = pData[pIndex( c, p, Ncomp, Np)] + 0.5*dt*vv[c];
          pData[pIndex( c, p, Ncomp, Np)] += dt*vv[c];

          rr = mygrid->dri[c] * (xx[c] - mygrid->rminloc[c]);
          rh = rr - 0.5;

          wii[c] = (int)floor(rr + 0.5); //whole integer int
          hii[c] = (int)floor(rr);     //half integer int
          rr -= wii[c];
          rh -= hii[c];
          rr2 = rr*rr;
          rh2 = rh*rh;

          wiw[c][1] = 0.75 - rr2;
          wiw[c][2] = 0.5*(0.25 + rr2 + rr);
          wiw[c][0] = 1. - wiw[c][1] - wiw[c][2];

          hiw[c][1] = 0.75 - rh2;
          hiw[c][2] = 0.5*(0.25 + rh2 + rh);
          hiw[c][0] = 1. - hiw[c][1] - hiw[c][2];
        }


        k1 = k2 = 0;
        for (j = 0; j < 3; j++)
        {
          j1 = j + wii[1] - 1;
          j2 = j + hii[1] - 1;
          for (i = 0; i < 3; i++)
          {
            i1 = i + wii[0] - 1;
            i2 = i + hii[0] - 1;
#ifndef OLD_ACCESS
            double weight =  pData[pIndex( 6, p, Ncomp, Np)]; //
            double *JX, *JY, *JZ;
            JX = &myCurrent[my_indice(edge,1, 0, 0, i2, j1, k1, N_grid[0], N_grid[1], N_grid[2], current->Ncomp)];
            JY = &myCurrent[my_indice(edge,1, 0, 1, i1, j2, k1, N_grid[0], N_grid[1], N_grid[2], current->Ncomp)];
            JZ = &myCurrent[my_indice(edge,1, 0, 2, i1, j1, k2, N_grid[0], N_grid[1], N_grid[2], current->Ncomp)];

            dvol = hiw[0][i] * wiw[1][j];
            *JX += weight*dvol*vv[0] * chargeSign;
            dvol = wiw[0][i] * hiw[1][j];
            *JY += weight*dvol*vv[1] * chargeSign;
            dvol = wiw[0][i] * wiw[1][j];
            *JZ += weight*dvol*vv[2] * chargeSign;

#else
            dvol = hiw[0][i] * wiw[1][j],
                current->Jx(i2, j1, k1) += w(p)*dvol*vv[0] * chargeSign;
            dvol = wiw[0][i] * hiw[1][j],
                current->Jy(i1, j2, k1) += w(p)*dvol*vv[1] * chargeSign;
            dvol = wiw[0][i] * wiw[1][j],
                current->Jz(i1, j1, k2) += w(p)*dvol*vv[2] * chargeSign;
#endif
          }
        }
      }
      break;

    case 1:
      for (p = 0; p < Np; p++)
      {
        double ux, uy, uz;
        ux = pData[pIndex( 3, p, Ncomp, Np)];
        uy = pData[pIndex( 4, p, Ncomp, Np)];
        uz = pData[pIndex( 5, p, Ncomp, Np)];

        gamma_i = 1. / sqrt(1 + ux*ux + uy*uy + uz*uz);
        for (c = 0; c < 3; c++)
        {
          vv[c] = gamma_i*pData[pIndex( c + 3, p, Ncomp, Np)];
          hiw[c][1] = wiw[c][1] = 1;
          hiw[c][0] = wiw[c][0] = 0;
          hiw[c][2] = wiw[c][2] = 0;
          hii[c] = wii[c] = 0;
        }
        for (c = 0; c < 1; c++)
        {
          xx[c] = pData[pIndex( c, p, Ncomp, Np)] + 0.5*dt*vv[c];
          pData[pIndex( c, p, Ncomp, Np)] += dt*vv[c];


          rr = mygrid->dri[c] * (xx[c] - mygrid->rminloc[c]);
          rh = rr - 0.5;

          wii[c] = (int)floor(rr + 0.5); //whole integer int
          hii[c] = (int)floor(rr);     //half integer int
          rr -= wii[c];
          rh -= hii[c];
          rr2 = rr*rr;
          rh2 = rh*rh;

          wiw[c][1] = 0.75 - rr2;
          wiw[c][2] = 0.5*(0.25 + rr2 + rr);
          wiw[c][0] = 1. - wiw[c][1] - wiw[c][2];

          hiw[c][1] = 0.75 - rh2;
          hiw[c][2] = 0.5*(0.25 + rh2 + rh);
          hiw[c][0] = 1. - hiw[c][1] - hiw[c][2];
        }


        k1 = k2 = j1 = j2 = 0;
        for (i = 0; i < 3; i++)
        {
          i1 = i + wii[0] - 1;
          i2 = i + hii[0] - 1;
#ifndef OLD_ACCESS

          double weight = pData[pIndex( 6, p, Ncomp, Np)];
          double *JX, *JY, *JZ;
          JX = &myCurrent[my_indice(edge,0, 0, 0, i2, j1, k1, N_grid[0], N_grid[1], N_grid[2], current->Ncomp)];
          JY = &myCurrent[my_indice(edge,0, 0, 1, i1, j2, k1, N_grid[0], N_grid[1], N_grid[2], current->Ncomp)];
          JZ = &myCurrent[my_indice(edge,0, 0, 2, i1, j1, k2, N_grid[0], N_grid[1], N_grid[2], current->Ncomp)];

          dvol = hiw[0][i];
          *JX += weight*dvol*vv[0] * chargeSign;
          dvol = wiw[0][i];
          *JY += weight*dvol*vv[1] * chargeSign;
          dvol = wiw[0][i];
          *JZ += weight*dvol*vv[2] * chargeSign;
#else

            dvol = hiw[0][i],
                current->Jx(i2, j1, k1) += myweight*dvol*vv[0] * chargeSign;
            dvol = wiw[0][i],
                current->Jy(i1, j2, k1) += myweight*dvol*vv[1] * chargeSign;
            dvol = wiw[0][i],
                current->Jz(i1, j1, k2) += myweight*dvol*vv[2] * chargeSign;
#endif
        }
      }
      break;
  }
}

void SPECIE::debug_warning_particle_outside_boundaries(double x, double y, double z, int nump){
  if (x < mygrid->rminloc[0]){
    std::cout << "Particle " << nump << " of proc. " << mygrid->myid << " is at x = " << x << " (boundary is " << mygrid->rminloc[0] << ")" << std::endl;
    flush(std::cout);
  }

  if (x > mygrid->rmaxloc[0]){
    std::cout << "Particle " << nump << " of proc. " << mygrid->myid << " is at x = " << x << " (boundary is " << mygrid->rmaxloc[0] << ")" << std::endl;
    flush(std::cout);
  }


  if (y < mygrid->rminloc[1]){
    std::cout << "Particle " << nump << " of proc. " << mygrid->myid << " is at y = " << y << " (boundary is " << mygrid->rminloc[1] << ")" << std::endl;
    flush(std::cout);
  }

  if (y > mygrid->rmaxloc[1]){
    std::cout << "Particle " << nump << " of proc. " << mygrid->myid << " is at y = " << y << " (boundary is " << mygrid->rmaxloc[1] << ")" << std::endl;
    flush(std::cout);
  }


  if (z < mygrid->rminloc[2]){
    std::cout << "Particle " << nump << " of proc. " << mygrid->myid << " is at z = " << z << " (boundary is " << mygrid->rminloc[2] << ")" << std::endl;
    flush(std::cout);
  }

  if (z > mygrid->rmaxloc[2]){
    std::cout << "Particle " << nump << " of proc. " << mygrid->myid << " is at z = " << z << " (boundary is " << mygrid->rmaxloc[2] << ")" << std::endl;
    flush(std::cout);
  }



}

void SPECIE::currentStretchedDepositionStandard(CURRENT *current)
{

  if (mygrid->withParticles == NO)
    return;

  double dt, gamma_i;
  int p, c;  // particle_int, component_int
  int i, j, k, i1, j1, k1, i2, j2, k2;
  //int indexMaxQuadraticShape[]={1,4};
  int hii[3], wii[3];           // half integer index,   whole integer index
  double hiw[3][3], wiw[3][3];  // half integer weight,  whole integer weight
  double rr, rh, rr2, rh2;          // local coordinate to integer grid point and to half integer,     local coordinate squared
  double dvol, xx[3], vv[3];           // tensor_product,       absolute particle position
  double mydr[3], myweight;
  double mycsi[3];

  dt = mygrid->dt;
  double* myCurrent = current->getDataPointer();
  int N_grid[3];
  current->writeN_grid(N_grid);
  int edge=mygrid->getEdge();
  if (mygrid->withCurrent == YES && (!isTestSpecies))
  {
    for (p = 0; p < Np; p++)
    {
      //debug_warning_particle_outside_boundaries(r0(p), r1(p), r2(p), p);
      double ux, uy, uz;
      ux = pData[pIndex( 3, p, Ncomp, Np)];
      uy = pData[pIndex( 4, p, Ncomp, Np)];
      uz = pData[pIndex( 5, p, Ncomp, Np)];

      gamma_i = 1. / sqrt(1 + ux*ux + uy*uy + uz*uz);

      for (c = 0; c < 3; c++)
      {
        vv[c] = gamma_i*pData[pIndex( c + 3, p, Ncomp, Np)];
        hiw[c][1] = wiw[c][1] = 1;
        hiw[c][0] = wiw[c][0] = 0;
        hiw[c][2] = wiw[c][2] = 0;
        hii[c] = wii[c] = 0;
      }
      for (c = 0; c < mygrid->getDimensionality(); c++)
      {
        xx[c] = pData[pIndex( c, p, Ncomp, Np)] + 0.5*dt*vv[c];
        pData[pIndex( c, p, Ncomp, Np)] += dt*vv[c];

        mycsi[c] = mygrid->unStretchGrid(xx[c], c);
        mydr[c] = mygrid->derivativeStretchingFunction(mycsi[c], c);
        rr = mygrid->dri[c] * (mycsi[c] - mygrid->csiminloc[c]);
        rh = rr - 0.5;

        wii[c] = (int)floor(rr + 0.5); //whole integer int
        hii[c] = (int)floor(rr);     //half integer int
        rr -= wii[c];
        rh -= hii[c];
        rr2 = rr*rr;
        rh2 = rh*rh;

        wiw[c][1] = 0.75 - rr2;
        wiw[c][2] = 0.5*(0.25 + rr2 + rr);
        wiw[c][0] = 1. - wiw[c][1] - wiw[c][2];

        hiw[c][1] = 0.75 - rh2;
        hiw[c][2] = 0.5*(0.25 + rh2 + rh);
        hiw[c][0] = 1. - hiw[c][1] - hiw[c][2];
      }
      switch (mygrid->getDimensionality())
      {
        case 3:
          myweight = w(p) / (mydr[0] * mydr[1] * mydr[2]);

          for (k = 0; k < 3; k++)
          {
            k1 = k + wii[2] - 1;
            k2 = k + hii[2] - 1;
            for (j = 0; j < 3; j++)
            {
              j1 = j + wii[1] - 1;
              j2 = j + hii[1] - 1;
              for (i = 0; i < 3; i++)
              {
                i1 = i + wii[0] - 1;
                i2 = i + hii[0] - 1;
#ifndef OLD_ACCESS
                double weight = pData[pIndex( 6, p, Ncomp, Np)];
                double *JX, *JY, *JZ;
                JX = &myCurrent[my_indice(edge,1, 1, 0, i2, j1, k1, N_grid[0], N_grid[1], N_grid[2], current->Ncomp)];
                JY = &myCurrent[my_indice(edge,1, 1, 1, i1, j2, k1, N_grid[0], N_grid[1], N_grid[2], current->Ncomp)];
                JZ = &myCurrent[my_indice(edge,1, 1, 2, i1, j1, k2, N_grid[0], N_grid[1], N_grid[2], current->Ncomp)];

                dvol = hiw[0][i] * wiw[1][j] * wiw[2][k];
                *JX += weight*dvol*vv[0] * chargeSign;
                dvol = wiw[0][i] * hiw[1][j] * wiw[2][k];
                *JY += weight*dvol*vv[1] * chargeSign;
                dvol = wiw[0][i] * wiw[1][j] * hiw[2][k];
                *JZ += weight*dvol*vv[2] * chargeSign;

#else
                dvol = hiw[0][i] * wiw[1][j] * wiw[2][k],
                    current->Jx(i2, j1, k1) += myweight*dvol*vv[0] * chargeSign;
                dvol = wiw[0][i] * hiw[1][j] * wiw[2][k],
                    current->Jy(i1, j2, k1) += myweight*dvol*vv[1] * chargeSign;
                dvol = wiw[0][i] * wiw[1][j] * hiw[2][k],
                    current->Jz(i1, j1, k2) += myweight*dvol*vv[2] * chargeSign;
 #endif

              }
            }
          }
          break;

        case 2:
          myweight = w(p) / (mydr[0] * mydr[1]);

          k1 = k2 = 0;
          for (j = 0; j < 3; j++)
          {
            j1 = j + wii[1] - 1;
            j2 = j + hii[1] - 1;
            for (i = 0; i < 3; i++)
            {
              i1 = i + wii[0] - 1;
              i2 = i + hii[0] - 1;
#ifndef OLD_ACCESS
                double weight =  pData[6 + p*Ncomp];//pData[pIndex( 6, p, Ncomp, Np)]; //
                double *JX, *JY, *JZ;
                JX = &myCurrent[my_indice(edge,1, 0, 0, i2, j1, k1, N_grid[0], N_grid[1], N_grid[2], current->Ncomp)];
                JY = &myCurrent[my_indice(edge,1, 0, 1, i1, j2, k1, N_grid[0], N_grid[1], N_grid[2], current->Ncomp)];
                JZ = &myCurrent[my_indice(edge,1, 0, 2, i1, j1, k2, N_grid[0], N_grid[1], N_grid[2], current->Ncomp)];

                dvol = hiw[0][i] * wiw[1][j];
                *JX += weight*dvol*vv[0] * chargeSign;
                dvol = wiw[0][i] * hiw[1][j];
                *JY += weight*dvol*vv[1] * chargeSign;
                dvol = wiw[0][i] * wiw[1][j];
                *JZ += weight*dvol*vv[2] * chargeSign;
#else
              dvol = hiw[0][i] * wiw[1][j],
                  current->Jx(i2, j1, k1) += myweight*dvol*vv[0] * chargeSign;
              dvol = wiw[0][i] * hiw[1][j],
                  current->Jy(i1, j2, k1) += myweight*dvol*vv[1] * chargeSign;
              dvol = wiw[0][i] * wiw[1][j],
                  current->Jz(i1, j1, k2) += myweight*dvol*vv[2] * chargeSign;
#endif
            }
          }
          break;

        case 1:
          myweight = w(p) / mydr[0];

          k1 = k2 = j1 = j2 = 0;
          for (i = 0; i < 3; i++)
          {
            i1 = i + wii[0] - 1;
            i2 = i + hii[0] - 1;
#ifndef OLD_ACCESS
                double weight =  pData[pIndex( 6, p, Ncomp, Np)];
                double *JX, *JY, *JZ;
                JX = &myCurrent[my_indice(edge,0, 0, 0, i2, j1, k1, N_grid[0], N_grid[1], N_grid[2], current->Ncomp)];
                JY = &myCurrent[my_indice(edge,0, 0, 1, i1, j2, k1, N_grid[0], N_grid[1], N_grid[2], current->Ncomp)];
                JZ = &myCurrent[my_indice(edge,0, 0, 2, i1, j1, k2, N_grid[0], N_grid[1], N_grid[2], current->Ncomp)];

                dvol = hiw[0][i];
                *JX += weight*dvol*vv[0] * chargeSign;
                dvol = wiw[0][i];
                *JY += weight*dvol*vv[1] * chargeSign;
                dvol = wiw[0][i];
                *JZ += weight*dvol*vv[2] * chargeSign;
#else

            dvol = hiw[0][i],
                current->Jx(i2, j1, k1) += myweight*dvol*vv[0] * chargeSign;
            dvol = wiw[0][i],
                current->Jy(i1, j2, k1) += myweight*dvol*vv[1] * chargeSign;
            dvol = wiw[0][i],
                current->Jz(i1, j1, k2) += myweight*dvol*vv[2] * chargeSign;
#endif
          }
          break;
      }

    }
  }
  else
  {
    for (p = 0; p < Np; p++)
    {
      double ux, uy, uz;
      ux = pData[pIndex( 3, p, Ncomp, Np)];
      uy = pData[pIndex( 4, p, Ncomp, Np)];
      uz = pData[pIndex( 5, p, Ncomp, Np)];

      gamma_i = 1. / sqrt(1 + ux*ux + uy*uy + uz*uz);
      for (c = 0; c < mygrid->getDimensionality(); c++)
      {
        vv[c] = gamma_i*pData[pIndex( c + 3, p, Ncomp, Np)];
        pData[pIndex( c, p, Ncomp, Np)] += dt*vv[c];
      }
    }
  }


}
void SPECIE::density_deposition_standard(CURRENT *current)
{
  if (mygrid->withParticles == NO){
    return;
  }


  if (mygrid->isStretched()){
    SPECIE::densityStretchedDepositionStandard(current);
    return;
  }

  int p, c;  // particle_int, component_int
  int i, j, k, i1, j1, k1;
  //int indexMaxQuadraticShape[]={1,4};
  int wii[3];           // whole integer index
  double wiw[3][3];  // whole integer weight
  double rr, rr2;          // local coordinate to integer grid point,     local coordinate squared
  double dvol, xx[3];           // tensor_product,       absolute particle position

  if (mygrid->withParticles != YES)
    return;

  for (p = 0; p < Np; p++)
  {
    //debug_warning_particle_outside_boundaries(r0(p), r1(p), r2(p), p);
    for (c = 0; c < mygrid->getDimensionality(); c++)
    {
      xx[c] = ru(c, p);

      rr = mygrid->dri[c] * (xx[c] - mygrid->rminloc[c]);
      wii[c] = (int)floor(rr + 0.5); //whole integer int
      rr -= wii[c];
      rr2 = rr*rr;

      wiw[c][1] = 0.75 - rr2;
      wiw[c][2] = 0.5*(0.25 + rr2 + rr);
      wiw[c][0] = 1. - wiw[c][1] - wiw[c][2];
    }
    switch (mygrid->getDimensionality())
    {
      case 3:
        for (k = 0; k < 3; k++)
        {
          k1 = k + wii[2] - 1;
          for (j = 0; j < 3; j++)
          {
            j1 = j + wii[1] - 1;
            for (i = 0; i < 3; i++)
            {
              i1 = i + wii[0] - 1;

              dvol = wiw[0][i] * wiw[1][j] * wiw[2][k],
                  current->density(i1, j1, k1) += w(p)*dvol;
            }
          }
        }
        break;

      case 2:
        k1 = 0;
        for (j = 0; j < 3; j++)
        {
          j1 = j + wii[1] - 1;
          for (i = 0; i < 3; i++)
          {
            i1 = i + wii[0] - 1;
            dvol = wiw[0][i] * wiw[1][j],
                current->density(i1, j1, k1) += w(p)*dvol;
          }
        }
        break;

      case 1:
        k1 = j1 = 0;
        for (i = 0; i < 3; i++)
        {
          i1 = i + wii[0] - 1;
          dvol = wiw[0][i],
              current->density(i1, j1, k1) += w(p)*dvol;
        }
        break;
    }

  }
}
void SPECIE::densityStretchedDepositionStandard(CURRENT *current)
{
  if (mygrid->withParticles == NO)
    return;

  int p, c;  // particle_int, component_int
  int i, j, k, i1, j1, k1;
  //int indexMaxQuadraticShape[]={1,4};
  int wii[3];           // whole integer index
  double wiw[3][3];  // whole integer weight
  double rr, rr2;          // local coordinate to integer grid point,     local coordinate squared
  double dvol, xx[3];           // tensor_product,       absolute particle position
  double mydr[3], myweight;
  double mycsi[3];


  for (p = 0; p < Np; p++)
  {
    //debug_warning_particle_outside_boundaries(r0(p), r1(p), r2(p), p);
    for (c = 0; c < mygrid->getDimensionality(); c++)
    {
      xx[c] = ru(c, p);
      mycsi[c] = mygrid->unStretchGrid(xx[c], c);
      mydr[c] = mygrid->derivativeStretchingFunction(mycsi[c], c);
      rr = mygrid->dri[c] * (mycsi[c] - mygrid->csiminloc[c]);

      wii[c] = (int)floor(rr + 0.5); //whole integer int
      rr -= wii[c];
      rr2 = rr*rr;

      wiw[c][1] = 0.75 - rr2;
      wiw[c][2] = 0.5*(0.25 + rr2 + rr);
      wiw[c][0] = 1. - wiw[c][1] - wiw[c][2];
    }
    switch (mygrid->getDimensionality())
    {
      case 3:
        myweight = w(p) / (mydr[0] * mydr[1] * mydr[2]);
        for (k = 0; k < 3; k++)
        {
          k1 = k + wii[2] - 1;
          for (j = 0; j < 3; j++)
          {
            j1 = j + wii[1] - 1;
            for (i = 0; i < 3; i++)
            {
              i1 = i + wii[0] - 1;

              dvol = wiw[0][i] * wiw[1][j] * wiw[2][k],
                  current->density(i1, j1, k1) += myweight*dvol;
            }
          }
        }
        break;

      case 2:
        myweight = w(p) / (mydr[0] * mydr[1]);
        k1 = 0;
        for (j = 0; j < 3; j++)
        {
          j1 = j + wii[1] - 1;
          for (i = 0; i < 3; i++)
          {
            i1 = i + wii[0] - 1;
            dvol = wiw[0][i] * wiw[1][j],
                current->density(i1, j1, k1) += myweight*dvol;
          }
        }
        break;

      case 1:
        myweight = w(p) / mydr[0];
        k1 = j1 = 0;
        for (i = 0; i < 3; i++)
        {
          i1 = i + wii[0] - 1;
          dvol = wiw[0][i],
              current->density(i1, j1, k1) += myweight*dvol;
        }
        break;
    }

  }
}

void SPECIE::setParticlesPerCellXYZ(int numX, int numY, int numZ){
  numX = (numX <= 0) ? 1 : numX;
  numY = (numY <= 0) ? 1 : numY;
  numZ = (numZ <= 0) ? 1 : numZ;
  particlePerCellXYZ[0] = numX;
  particlePerCellXYZ[1] = numY;
  particlePerCellXYZ[2] = numZ;
}

void SPECIE::setName(std::string iname){
  name = iname;
}

double SPECIE::getKineticEnergy(){
  if (mygrid->withParticles == NO){
    return 0.0;
  }

  if (!allocated){
    return 0.0;
  }
  double energy = 0.0;
  for (int p = 0; p < Np; p++){
    energy += (sqrt(1.0 + (u0(p)*u0(p) + u1(p)*u1(p) + u2(p)*u2(p))) - 1.0)*w(p);
  }
  energy *= mygrid->dr[0] * mygrid->dr[1] * mygrid->dr[2] * mygrid->ref_den*M_PI / coupling*chargeSign;
  return energy;
}

//xmin, ymin, zmin, pxmin,pymin,pzmin,emin, xmax,ymax,zmax, pxmax,pymax,pzmax,emax
void SPECIE::computeKineticEnergyWExtrems(){

  if (mygrid->withParticles == NO){
    return;
  }
  if (!allocated){
    return;
  }

  if (energyExtremesFlag){
    return;
  }

  const double VERY_BIG_NUM_POS = 1.0e30;
  const double VERY_BIG_NUM_NEG = -1.0e30;

  for (int i = 0; i < 7; i++){
    minima[i] = VERY_BIG_NUM_POS;
    maxima[i] = VERY_BIG_NUM_NEG;
  }


  double energy = 0.0;
  totalEnergy = 0.0;
  totalMomentum[0] = 0;
  totalMomentum[1] = 0;
  totalMomentum[2] = 0;
  double gamma_minus_1 = 0;

  for (int p = 0; p < Np; p++){
    gamma_minus_1 = (sqrt(1.0 + (u0(p)*u0(p) + u1(p)*u1(p) + u2(p)*u2(p))) - 1.0);
    energy += gamma_minus_1*w(p);
    totalMomentum[0] += u0(p)*w(p);
    totalMomentum[1] += u1(p)*w(p);
    totalMomentum[2] += u2(p)*w(p);

    if (r0(p) <= minima[0])minima[0] = r0(p);
    if (r1(p) <= minima[1])minima[1] = r1(p);
    if (r2(p) <= minima[2])minima[2] = r2(p);

    if (r0(p) >= maxima[0])maxima[0] = r0(p);
    if (r1(p) >= maxima[1])maxima[1] = r1(p);
    if (r2(p) >= maxima[2])maxima[2] = r2(p);

    if (u0(p) <= minima[3])minima[3] = u0(p);
    if (u1(p) <= minima[4])minima[4] = u1(p);
    if (u2(p) <= minima[5])minima[5] = u2(p);

    if (u0(p) >= maxima[3])maxima[3] = u0(p);
    if (u1(p) >= maxima[4])maxima[4] = u1(p);
    if (u2(p) >= maxima[5])maxima[5] = u2(p);

    if (gamma_minus_1 <= minima[6])minima[6] = gamma_minus_1;
    if (gamma_minus_1 >= maxima[6])maxima[6] = gamma_minus_1;
  }
  energy *= mygrid->dr[0] * mygrid->dr[1] * mygrid->dr[2] * mygrid->ref_den*M_PI / coupling*chargeSign;
  totalMomentum[0] *= mass*mygrid->dr[0] * mygrid->dr[1] * mygrid->dr[2] * mygrid->ref_den*M_PI / coupling*chargeSign;
  totalMomentum[1] *= mass*mygrid->dr[0] * mygrid->dr[1] * mygrid->dr[2] * mygrid->ref_den*M_PI / coupling*chargeSign;
  totalMomentum[2] *= mass*mygrid->dr[0] * mygrid->dr[1] * mygrid->dr[2] * mygrid->ref_den*M_PI / coupling*chargeSign;


  MPI_Allreduce(&energy, &totalEnergy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, totalMomentum, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, minima, 7, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, maxima, 7, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  spectrum.Kmax = maxima[6];
  spectrum.Nbin = NBIN_SPECTRUM;
  spectrum.values = (double*)realloc((void*)spectrum.values, spectrum.Nbin*sizeof(double));
  spectrum.Dk = spectrum.Kmax / spectrum.Nbin;
  double Dki = 1 / spectrum.Dk;
  memset((void*)spectrum.values, 0, spectrum.Nbin*sizeof(double));
  for (int p = 0; p < Np; p++){
    int ibin;
    gamma_minus_1 = (sqrt(1.0 + (u0(p)*u0(p) + u1(p)*u1(p) + u2(p)*u2(p))) - 1.0);
    ibin = (int)(gamma_minus_1 / spectrum.Dk);
    if (ibin < 0)
      ibin = 0;
    if (ibin >= spectrum.Nbin)
      ibin = spectrum.Nbin - 1;
    spectrum.values[ibin] += Dki*w(p)*mygrid->ref_den;
  }
  MPI_Allreduce(MPI_IN_PLACE, spectrum.values, spectrum.Nbin, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);



  energyExtremesFlag = true;
}

void SPECIE::dump(std::ofstream &ff){
  ff.write((char*)&Np, sizeof(Np));

  for (int i = 0; i < Np; i++){
    for (int c = 0; c < Ncomp; c++){
      ff.write((char*)&ru(c, i), sizeof(double));
    }
  }
}
void SPECIE::dumpBigBuffer(std::ofstream &ff){
  ff.write((char*)&Np, sizeof(Np));
  ff.write((char*)pData, sizeof(double)*Np*Ncomp);
}
void SPECIE::debugDump(std::ofstream &ff){
  ff << this->name << Np << std::endl;
}

void SPECIE::reloadDump(std::ifstream &ff){
  ff.read((char*)&Np, sizeof(Np));
  SPECIE::reallocate_species();
  for (int i = 0; i < Np; i++){
    for (int c = 0; c < Ncomp; c++){
      ff.read((char*)&ru(c, i), sizeof(double));
    }
  }
}

void SPECIE::reloadBigBufferDump(std::ifstream &ff){
  ff.read((char*)&Np, sizeof(Np));
  SPECIE::reallocate_species();
  ff.read((char*)&ru(0, 0), sizeof(double)*Np*Ncomp);
}

bool SPECIE::areEnergyExtremesAvailable(){
  return energyExtremesFlag;
}


uint64_t SPECIE::printParticleNumber(){
  if (mygrid->myid != mygrid->master_proc)
    return lastParticle;
  std::cout << name << " has " << lastParticle << " particles." << std::endl;
  return lastParticle;
}


