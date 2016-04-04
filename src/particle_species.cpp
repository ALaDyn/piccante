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

#include "particle_species.h"
//#define OLD_ACCESS
//#define OLDPUSHER
//#define OLDCURRENT


SPECIE::SPECIE()
{
  Ncomp = 7;
  allocated = false;
  Z = A = 0;
  isTestSpecies = false;
  isFrozen = false;
  isQuiet = false;
  quietShuffle = 1;
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
  isFrozen = false;
  isQuiet = false;
  quietShuffle = 1;
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
  for (int c = 0; c < Ncomp; c++) {
    val[c] = (double*)malloc(Np*sizeof(double));
  }
#endif

  valSize = Np;
  allocated = true;
#endif

}
SPECIE::~SPECIE() {
#ifdef _ACC_SINGLE_POINTER
  free(pData);
#else
  for (int c = 0; c < Ncomp; c++) {
    free(val[c]);
  }
  free(val);
#endif
}

void SPECIE::erase()
{
  if (mygrid->withParticles == NO)
    return;

  if (!allocated) {
    printf("ERROR: species not allocated!!!\n");
    exit(11);
  }
#ifdef _ACC_SINGLE_POINTER
  memset((void*)pData, 0, (Np*Ncomp)*sizeof(double));
#else
  for (int c = 0; c < Ncomp; c++) {
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
  if (Np > valSize) {
    valSize = Np + allocsize;
#ifdef _ACC_SINGLE_POINTER
    pData = (double *)realloc((void*)pData, valSize*Ncomp*sizeof(double));
#else
    for (int c = 0; c < Ncomp; c++) {
      val[c] = (double *)realloc((void*)val[c], valSize*sizeof(double));
    }
#endif
  }
  else if (Np < (valSize - allocsize)) {
    valSize = Np + allocsize;
#ifdef _ACC_SINGLE_POINTER
    pData = (double *)realloc((void*)pData, valSize*Ncomp*sizeof(double));
#else
    for (int c = 0; c < Ncomp; c++) {
      val[c] = (double *)realloc((void*)val[c], valSize*sizeof(double));
    }
#endif
  }
  return;

}

SPECIE SPECIE::operator = (SPECIE &destro)
{
  if (!destro.allocated) { printf("---ERROR---\noperation not permitted\nSPECIE=SPECIE\nnot allocated\n"); exit(11); }

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
  isFrozen = destro.isFrozen;
  isQuiet = destro.isQuiet;
  quietShuffle = destro.quietShuffle;
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
  for (int c = 0; c < Ncomp; c++) {
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

std::string SPECIE::getName() {
  return name;
}

void SPECIE::addMarker() {
  flagWithMarker = true;
}
void SPECIE::setTestSpecies() {
  isTestSpecies = true;
}
void SPECIE::setFrozenSpecies() {
  isFrozen = true;
  isTestSpecies = true;
}
void SPECIE::setQuietStart(){
  isQuiet = true;
  std::cout<< "setQuietStart()!!!" << std::endl;
}
void SPECIE::setQuietShuffle(int Nshuffle){
  quietShuffle = Nshuffle;
}

bool SPECIE::amIWithMarker() {
  return flagWithMarker;
}

void SPECIE::computeParticleMassChargeCoupling() {
  if (type == ELECTRON) {
    coupling = -1.;
    mass = 1.0;
    Z = -1.0;
    chargeSign = -1.0;
  }
  if (type == POSITRON) {
    coupling = 1.;
    mass = 1.0;
    Z = 1.0;
    chargeSign = 1.0;
  }
  if (type == ION) {
    if (Z == 0 || A == 0)
    {
      printf("ERROR: Ion charge or mass NOT defined!\n");
      exit(11);
    }
    else {
      coupling = Z / (1.8362e3*A);
      mass = 1.8362e3*A;
    }
    chargeSign = 1.0;
  }
}
int SPECIE::getNumberOfParticlesWithin(double plasmarmin[3], double plasmarmax[3]) {

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
            if (zloc >= plasmarmin[2] && zloc <= plasmarmax[2]) {
              if (plasma.density_function(xloc, yloc, zloc, plasma.params, Z, A) > 0)
                counter += particlePerCell;
            }
      }
  return counter;
}
int SPECIE::getNumberOfParticlesWithinFromFile1D(double plasmarmin[], double plasmarmax[], std::string name) {

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
  for (int i = 0; i < Nx_in; i++) {
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
            if (zloc >= plasmarmin[2] && zloc <= plasmarmax[2]) {
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
void SPECIE::createParticlesWithinFrom(double plasmarmin[3], double plasmarmax[3], int oldNumberOfParticles, uint64_t disp) {
#ifdef NO_ALLOCATION
  if (!allocated) {
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
              if (plasma.density_function(xloc, yloc, zloc, plasma.params, Z, A) > 0){
                weight = plasma.density_function(xloc, yloc, zloc, plasma.params, Z, A) / particlePerCell;
                xloc -= 0.5*dx;
                yloc -= 0.5*dy;
                zloc -= 0.5*dz;


                  for (int ip = 0; ip < particlePerCellXYZ[0]; ip++)
                    for (int jp = 0; jp < particlePerCellXYZ[1]; jp++)
                      for (int kp = 0; kp < particlePerCellXYZ[2]; kp++){
                        r0(counter) = xloc + dxp*(ip + 0.5);
                        r1(counter) = yloc + dyp*(jp + 0.5);
                        r2(counter) = zloc + dzp*(kp + 0.5);
                        u0(counter) = u1(counter) = u2(counter) = 0;
                        w(counter) = weight;
                        if (flagWithMarker)
                          marker(counter) = (counter + disp);
                        if (isTestSpecies)
                          w(counter) = (double)(counter + disp);
                        counter++;
                      }

              }

            }
      }
}

void SPECIE::createStretchedParticlesWithinFrom(double plasmarmin[3], double plasmarmax[3], int oldNumberOfParticles, uint64_t disp) {
#ifdef NO_ALLOCATION
  if (!allocated) {
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

  for (int k = 0; k < Nz; k++) {
    zloc = mygrid->chrloc[2][k];
    csilocz = mygrid->csiminloc[2] + dz*k;
    for (int j = 0; j < Ny; j++) {
      yloc = mygrid->chrloc[1][j];
      csilocy = mygrid->csiminloc[1] + dy*j;
      for (int i = 0; i < Nx; i++) {
        xloc = mygrid->chrloc[0][i];
        csilocx = mygrid->csiminloc[0] + dx*i;

        if (xloc >= plasmarmin[0] && xloc <= plasmarmax[0]) {
          if (yloc >= plasmarmin[1] && yloc <= plasmarmax[1]) {
            if (zloc >= plasmarmin[2] && zloc <= plasmarmax[2])
            {
              if (plasma.density_function(xloc, yloc, zloc, plasma.params, Z, A) > 0)
              {
                weight = plasma.density_function(xloc, yloc, zloc, plasma.params, Z, A) / particlePerCell;
                for (int kp = 0; kp < particlePerCellXYZ[2]; kp++) {
                  mycsiz = csilocz + dzp*(kp + 0.5);
                  myz = mygrid->stretchGrid(mycsiz, 2);
                  mydz = mygrid->derivativeStretchingFunction(mycsiz, 2);

                  for (int jp = 0; jp < particlePerCellXYZ[1]; jp++) {
                    mycsiy = csilocy + dyp*(jp + 0.5);
                    myy = mygrid->stretchGrid(mycsiy, 1);
                    mydy = mygrid->derivativeStretchingFunction(mycsiy, 1);

                    for (int ip = 0; ip < particlePerCellXYZ[0]; ip++) {
                      mycsix = csilocx + dxp*(ip + 0.5);
                      myx = mygrid->stretchGrid(mycsix, 0);
                      mydx = mygrid->derivativeStretchingFunction(mycsix, 0);

                      r0(counter) = myx;
                      r1(counter) = myy;
                      r2(counter) = myz;
                      u0(counter) = u1(counter) = u2(counter) = 0;
                      w(counter) = weight*mydx*mydy*mydz;
                      if (flagWithMarker)
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

void SPECIE::createParticlesWithinFromButFromFile1D(double plasmarmin[3], double plasmarmax[3], int oldNumberOfParticles, uint64_t disp, std::string name) {
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
  for (int i = 0; i < Nx_in; i++) {
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

void SPECIE::createStretchedParticlesWithinFromButFromFile1D(double plasmarmin[3], double plasmarmax[3], int oldNumberOfParticles, uint64_t disp, std::string name) {
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
  for (int i = 0; i < Nx_in; i++) {
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

  for (int k = 0; k < Nz; k++) {
    zloc = mygrid->chrloc[2][k];
    csilocz = mygrid->csiminloc[2] + dz*k;
    for (int j = 0; j < Ny; j++) {
      yloc = mygrid->chrloc[1][j];
      csilocy = mygrid->csiminloc[1] + dy*j;
      for (int i = 0; i < Nx; i++) {
        xloc = mygrid->chrloc[0][i];
        csilocx = mygrid->csiminloc[0] + dx*i;

        ih = (int)((xloc - xmin) / dx_in);
        axh = (xloc - xmin) / dx_in - ih;
        wh[0] = 1 - axh;
        wh[1] = axh;
        ihleft = (ih + Nx_in - 1) % (Nx_in - 1);
        ihright = (ih + 1) % (Nx_in - 1);

        if (xloc >= plasmarmin[0] && xloc <= plasmarmax[0]) {
          if (yloc >= plasmarmin[1] && yloc <= plasmarmax[1]) {
            if (zloc >= plasmarmin[2] && zloc <= plasmarmax[2])
            {
              denValue = wh[0] * density[ihleft] + wh[1] * density[ihright];

              if (denValue > 0)
              {
                velValue = wh[0] * velocity[ihleft] + wh[1] * velocity[ihright];
                weight = denValue / particlePerCell;
                myuy = wh[0] * uy[ihleft] + wh[1] * uy[ihright];
                myuz = wh[0] * uz[ihleft] + wh[1] * uz[ihright];
                for (int kp = 0; kp < particlePerCellXYZ[2]; kp++) {
                  mycsiz = csilocz + dzp*(kp + 0.5);
                  myz = mygrid->stretchGrid(mycsiz, 2);
                  mydz = mygrid->derivativeStretchingFunction(mycsiz, 2);

                  for (int jp = 0; jp < particlePerCellXYZ[1]; jp++) {
                    mycsiy = csilocy + dyp*(jp + 0.5);
                    myy = mygrid->stretchGrid(mycsiy, 1);
                    mydy = mygrid->derivativeStretchingFunction(mycsiy, 1);

                    for (int ip = 0; ip < particlePerCellXYZ[0]; ip++) {
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

void SPECIE::setNumberOfParticlePerCell() {
  particlePerCell = 1;
  for (int c = 0; c < 3; c++) {
    if (c < mygrid->getDimensionality()) {
      particlePerCell *= particlePerCellXYZ[c];   //number of particles per cell(npc) total = npcx*npcy*npcz
    }
    else {
      particlePerCellXYZ[c] = 1;
    }
  }
}
void SPECIE::setLocalPlasmaMinimaAndMaxima(double *plasmarmin, double *plasmarmax) {
  for (int c = 0; c < 3; c++) {
    if (c < mygrid->getDimensionality()) {
      plasmarmin[c] = MAX(plasma.params.rminbox[c], mygrid->rminloc[c]);
      plasmarmax[c] = MIN(plasma.params.rmaxbox[c], mygrid->rmaxloc[c]);
    }
    else {
      plasmarmin[c] = mygrid->rminloc[c];
      plasmarmax[c] = mygrid->rmaxloc[c];
    }
  }

}
uint64_t SPECIE::getSumNewParticlesOfAllPreviousProcessors(int number) {
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

  if (flagWithMarker) {
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


void SPECIE::creationFromFile1D(std::string name) {
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
  for (int c = 0; c < 3; c++) {
    if (c < 1) {
      plasmarmin[c] = MAX(plasmaXMin, mygrid->rminloc[c]);
      plasmarmax[c] = MIN(plasmaXMax, mygrid->rmaxloc[c]);
      particlePerCell *= particlePerCellXYZ[c];   //number of particles per cell(npc) total = npcx*npcy*npcz
    }
    else {
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

  if (mygrid->rmyid[0] == (mygrid->rnproc[0] - 1)) {
    newNumberOfParticles = SPECIE::getNumberOfParticlesWithin(plasmarmin, plasmarmax);
  }
  else {
    newNumberOfParticles = 0;
  }
  Np += newNumberOfParticles;
  reallocate_species();

  uint64_t disp = getSumNewParticlesOfAllPreviousProcessors(newNumberOfParticles);

  if ((mygrid->rmyid[0] == (mygrid->rnproc[0] - 1)) && newNumberOfParticles > 0) {

    if (mygrid->isStretched())
      createStretchedParticlesWithinFrom(plasmarmin, plasmarmax, oldNumberOfParticles, disp);
    else
      createParticlesWithinFrom(plasmarmin, plasmarmax, oldNumberOfParticles, disp);
  }

}
//void SPECIE::output_bin(ofstream &ff)
//{
//  if (mygrid->with_particles == NO)
//    return;

//  ff.write((char *)val, sizeof(double)*Ncomp*Np);
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

//void SPECIE::init_output_diag(std::ofstream &ff){
//  if (mygrid->withParticles == NO)
//    return;

//  if (mygrid->myid == mygrid->master_proc){
//    ff << std::setw(myNarrowWidth) << "#step" << " " << std::setw(myWidth) << "time" << " " << std::setw(myWidth) << "Etot";
//    ff << " " << std::setw(myWidth) << "Px" << " " << std::setw(myWidth) << "Py" << " " << std::setw(myWidth) << "Pz" << std::endl;
//  }
//}
//void SPECIE::output_diag(int istep, std::ofstream &ff){
//  if (mygrid->withParticles == NO)
//    return;

//  //double extrema[14];
//  computeKineticEnergyWExtrems();
//  if (mygrid->myid == mygrid->master_proc){
//    ff << std::setw(myWidth) << istep << " " << std::setw(myWidth) << mygrid->time << " " << std::setw(myWidth) << totalEnergy;
//    for (int c = 0; c < 3; c++){
//      ff << " " << std::setw(myWidth) << totalMomentum[c];
//    }
//    ff << std::endl;
//  }
//}

void SPECIE::init_output_extrems(std::ofstream &ff) {
  if (mygrid->withParticles == NO)
    return;

  if (mygrid->myid == mygrid->master_proc) {
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
void SPECIE::output_extrems(int istep, std::ofstream &ff) {
  if (mygrid->withParticles == NO)
    return;

  //double extrema[14];
  computeKineticEnergyWExtrems();
  if (mygrid->myid == mygrid->master_proc) {
    ff << std::setw(myNarrowWidth) << istep << " " << std::setw(myWidth) << mygrid->time;
    for (int c = 0; c < 7; c++) {
      ff << " " << std::setw(myWidth) << minima[c] << " " << std::setw(myWidth) << maxima[c];
    }
    ff << std::endl;
  }
}

void SPECIE::outputSpectrum(std::ofstream &fspectrum) {
  computeKineticEnergyWExtrems();
  if (mygrid->myid == mygrid->master_proc) {
    fspectrum << "#" << std::setw(myWidth) << "Ebinmin";
    fspectrum << " " << std::setw(myWidth) << "Ebinmax";
    fspectrum << " " << std::setw(myWidth) << "value";
    fspectrum << std::endl;
    for (int ibin = 0; ibin < spectrum.Nbin; ibin++) {
      fspectrum << " " << std::setw(myWidth) << (ibin*spectrum.Dk);
      fspectrum << " " << std::setw(myWidth) << ((ibin + 1)*spectrum.Dk);
      fspectrum << " " << std::setw(myWidth) << (spectrum.values[ibin]);
      fspectrum << std::endl;

    }
  }
}

void SPECIE::position_parallel_pbc()
{
  if (mygrid->withParticles == NO||isFrozen)
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
    double Length = (mygrid->rmax[direction] - mygrid->rmin[direction]);
    for (p = 0; p < Np; p++)
    {
#ifdef _ACC_SINGLE_POINTER
      if (pData[direction + p*Ncomp] > mygrid->rmaxloc[direction])
      {
        nlost++;
        nright++;
        if (mygrid->rmyid[direction] == mygrid->rnproc[direction] - 1)
          pData[direction + p*Ncomp] -= Length;
        sendr_buffer = (double*)realloc(sendr_buffer, nright*Ncomp*sizeof(double));
        for (c = 0; c < Ncomp; c++)
          sendr_buffer[c + Ncomp*(nright - 1)] = pData[c + p*Ncomp];

      }

      else if (pData[direction + p*Ncomp] < mygrid->rminloc[direction])
      {
        nlost++;
        nleft++;
        if (mygrid->rmyid[direction] == 0)
          pData[direction + p*Ncomp] += Length;
        sendl_buffer = (double*)realloc(sendl_buffer, nleft*Ncomp*sizeof(double));
        for (c = 0; c < Ncomp; c++)
          sendl_buffer[c + Ncomp*(nleft - 1)] = pData[c + p*Ncomp];
      }
      else
      {
        for (c = 0; c < Ncomp; c++)
          pData[c + (p - nlost)*Ncomp] = pData[c + p*Ncomp];
      }
#else
      if (val[direction][p*Ncomp] > mygrid->rmaxloc[direction])
      {
        nlost++;
        nright++;
        if (mygrid->rmyid[direction] == mygrid->rnproc[direction] - 1)
          val[direction][p*Ncomp] -= Length;
        sendr_buffer = (double*)realloc(sendr_buffer, nright*Ncomp*sizeof(double));
        for (c = 0; c < Ncomp; c++)
          sendr_buffer[c + Ncomp*(nright - 1)] = val[c][p*Ncomp];

      }
      else if (val[direction][p*Ncomp] < mygrid->rminloc[direction])
      {
        nlost++;
        nleft++;
        if (mygrid->rmyid[direction] == 0)
          val[direction][p*Ncomp] += Length;
        sendl_buffer = (double*)realloc(sendl_buffer, nleft*Ncomp*sizeof(double));
        for (c = 0; c < Ncomp; c++)
          sendl_buffer[c + Ncomp*(nleft - 1)] = val[c][p*Ncomp];
      }
      else
      {
        for (c = 0; c < Ncomp; c++)
          val[c][(p - nlost)*Ncomp] = val[c][p*Ncomp];
      }
#endif
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

    #ifdef _ACC_SINGLE_POINTER
    for (int pp = 0; pp < nnew; pp++) {
      for (c = 0; c < Ncomp; c++) {
        pData[c + (pp + nold - nlost)*Ncomp] = recv_buffer[pp*Ncomp + c];
      }
    }
#else
      for (c = 0; c < Ncomp; c++) {
        for (int pp = 0; pp < nnew; pp++) {
        val[c][(pp + nold - nlost)*Ncomp] = recv_buffer[pp*Ncomp + c];
      }
    }
#endif
  }

  free(sendl_buffer);
  free(sendr_buffer);
  free(recv_buffer);
}
void SPECIE::position_obc()
{
  if (mygrid->withParticles == NO||isFrozen)
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

  if (mygrid->withParticles == NO||isFrozen)
    return;
  if (mygrid->isStretched()) {
    SPECIE::momentaStretchedAdvance(ebfield);
    return;
  }
  double dt, gamma_i;
  int p, c;  // particle_int, component_int
  int i, i1, j1, k1, i2, j2, k2;
#ifdef OLDPUSHER
  int j, k;
#endif
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
  int Nx, Ny, Nz;
  Nx = N_grid[0];
  Ny = N_grid[1];
  Nz = N_grid[2];
  switch (mygrid->getDimensionality())
  {

  case 3:
    for (p = 0; p < Np; p++)
    {
      //gamma_i=1./sqrt(1+u0(p)*u0(p)+u1(p)*u1(p)+u2(p)*u2(p));
      for (c = 0; c < 3; c++) {
#ifdef _ACC_SINGLE_POINTER
        xx[c] = pData[c + p*Ncomp];
#else
        xx[c] = val[c][p*Ncomp];
#endif
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



#ifdef OLDPUSHER
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
            EX = myfield[my_indice(edge, 1, 1, 0, i2, j1, k1, N_grid[0], N_grid[1], N_grid[2], ebComp)];
            EY = myfield[my_indice(edge, 1, 1, 1, i1, j2, k1, N_grid[0], N_grid[1], N_grid[2], ebComp)];
            EZ = myfield[my_indice(edge, 1, 1, 2, i1, j1, k2, N_grid[0], N_grid[1], N_grid[2], ebComp)];
            double BX, BY, BZ;
            BX = myfield[my_indice(edge, 1, 1, 3, i1, j2, k2, N_grid[0], N_grid[1], N_grid[2], ebComp)];
            BY = myfield[my_indice(edge, 1, 1, 4, i2, j1, k2, N_grid[0], N_grid[1], N_grid[2], ebComp)];
            BZ = myfield[my_indice(edge, 1, 1, 5, i2, j2, k1, N_grid[0], N_grid[1], N_grid[2], ebComp)];

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
#else
      {
        double EXLLL, EXCLL, EXRLL, EXLCL, EXCCL, EXRCL, EXLRL, EXCRL, EXRRL;
        double EYLLL, EYCLL, EYRLL, EYLCL, EYCCL, EYRCL, EYLRL, EYCRL, EYRRL;
        double EZLLL, EZCLL, EZRLL, EZLCL, EZCCL, EZRCL, EZLRL, EZCRL, EZRRL;

        double EXLLC, EXCLC, EXRLC, EXLCC, EXCCC, EXRCC, EXLRC, EXCRC, EXRRC;
        double EYLLC, EYCLC, EYRLC, EYLCC, EYCCC, EYRCC, EYLRC, EYCRC, EYRRC;
        double EZLLC, EZCLC, EZRLC, EZLCC, EZCCC, EZRCC, EZLRC, EZCRC, EZRRC;

        double EXLLR, EXCLR, EXRLR, EXLCR, EXCCR, EXRCR, EXLRR, EXCRR, EXRRR;
        double EYLLR, EYCLR, EYRLR, EYLCR, EYCCR, EYRCR, EYLRR, EYCRR, EYRRR;
        double EZLLR, EZCLR, EZRLR, EZLCR, EZCCR, EZRCR, EZLRR, EZCRR, EZRRR;

        double BXLLL, BXCLL, BXRLL, BXLCL, BXCCL, BXRCL, BXLRL, BXCRL, BXRRL;
        double BYLLL, BYCLL, BYRLL, BYLCL, BYCCL, BYRCL, BYLRL, BYCRL, BYRRL;
        double BZLLL, BZCLL, BZRLL, BZLCL, BZCCL, BZRCL, BZLRL, BZCRL, BZRRL;

        double BXLLC, BXCLC, BXRLC, BXLCC, BXCCC, BXRCC, BXLRC, BXCRC, BXRRC;
        double BYLLC, BYCLC, BYRLC, BYLCC, BYCCC, BYRCC, BYLRC, BYCRC, BYRRC;
        double BZLLC, BZCLC, BZRLC, BZLCC, BZCCC, BZRCC, BZLRC, BZCRC, BZRRC;

        double BXLLR, BXCLR, BXRLR, BXLCR, BXCCR, BXRCR, BXLRR, BXCRR, BXRRR;
        double BYLLR, BYCLR, BYRLR, BYLCR, BYCCR, BYRCR, BYLRR, BYCRR, BYRRR;
        double BZLLR, BZCLR, BZRLR, BZLCR, BZCCR, BZRCR, BZLRR, BZCRR, BZRRR;
        int iiL, iiC, iiR, hiL, hiC, hiR;
        iiL = wii[0] - 1;
        iiC = wii[0];
        iiR = wii[0] + 1;
        hiL = hii[0] - 1;
        hiC = hii[0];
        hiR = hii[0] + 1;

        int ijL, ijC, ijR, hjL, hjC, hjR;
        ijL = wii[1] - 1;
        ijC = wii[1];
        ijR = wii[1] + 1;
        hjL = hii[1] - 1;
        hjC = hii[1];
        hjR = hii[1] + 1;

        int ikL, ikC, ikR, hkL, hkC, hkR;
        ikL = wii[2] - 1;
        ikC = wii[2];
        ikR = wii[2] + 1;
        hkL = hii[2] - 1;
        hkC = hii[2];
        hkR = hii[2] + 1;

        {
          EXLLL = myfield[my_indice(edge, 1, 1, 0, hiL, ijL, ikL, Nx, Ny, Nz, ebComp)];
          EXCLL = myfield[my_indice(edge, 1, 1, 0, hiC, ijL, ikL, Nx, Ny, Nz, ebComp)];
          EXRLL = myfield[my_indice(edge, 1, 1, 0, hiR, ijL, ikL, Nx, Ny, Nz, ebComp)];
          EXLCL = myfield[my_indice(edge, 1, 1, 0, hiL, ijC, ikL, Nx, Ny, Nz, ebComp)];
          EXCCL = myfield[my_indice(edge, 1, 1, 0, hiC, ijC, ikL, Nx, Ny, Nz, ebComp)];
          EXRCL = myfield[my_indice(edge, 1, 1, 0, hiR, ijC, ikL, Nx, Ny, Nz, ebComp)];
          EXLRL = myfield[my_indice(edge, 1, 1, 0, hiL, ijR, ikL, Nx, Ny, Nz, ebComp)];
          EXCRL = myfield[my_indice(edge, 1, 1, 0, hiC, ijR, ikL, Nx, Ny, Nz, ebComp)];
          EXRRL = myfield[my_indice(edge, 1, 1, 0, hiR, ijR, ikL, Nx, Ny, Nz, ebComp)];

          EYLLL = myfield[my_indice(edge, 1, 1, 1, iiL, hjL, ikL, Nx, Ny, Nz, ebComp)];
          EYCLL = myfield[my_indice(edge, 1, 1, 1, iiC, hjL, ikL, Nx, Ny, Nz, ebComp)];
          EYRLL = myfield[my_indice(edge, 1, 1, 1, iiR, hjL, ikL, Nx, Ny, Nz, ebComp)];
          EYLCL = myfield[my_indice(edge, 1, 1, 1, iiL, hjC, ikL, Nx, Ny, Nz, ebComp)];
          EYCCL = myfield[my_indice(edge, 1, 1, 1, iiC, hjC, ikL, Nx, Ny, Nz, ebComp)];
          EYRCL = myfield[my_indice(edge, 1, 1, 1, iiR, hjC, ikL, Nx, Ny, Nz, ebComp)];
          EYLRL = myfield[my_indice(edge, 1, 1, 1, iiL, hjR, ikL, Nx, Ny, Nz, ebComp)];
          EYCRL = myfield[my_indice(edge, 1, 1, 1, iiC, hjR, ikL, Nx, Ny, Nz, ebComp)];
          EYRRL = myfield[my_indice(edge, 1, 1, 1, iiR, hjR, ikL, Nx, Ny, Nz, ebComp)];

          EZLLL = myfield[my_indice(edge, 1, 1, 2, iiL, ijL, hkL, Nx, Ny, Nz, ebComp)];
          EZCLL = myfield[my_indice(edge, 1, 1, 2, iiC, ijL, hkL, Nx, Ny, Nz, ebComp)];
          EZRLL = myfield[my_indice(edge, 1, 1, 2, iiR, ijL, hkL, Nx, Ny, Nz, ebComp)];
          EZLCL = myfield[my_indice(edge, 1, 1, 2, iiL, ijC, hkL, Nx, Ny, Nz, ebComp)];
          EZCCL = myfield[my_indice(edge, 1, 1, 2, iiC, ijC, hkL, Nx, Ny, Nz, ebComp)];
          EZRCL = myfield[my_indice(edge, 1, 1, 2, iiR, ijC, hkL, Nx, Ny, Nz, ebComp)];
          EZLRL = myfield[my_indice(edge, 1, 1, 2, iiL, ijR, hkL, Nx, Ny, Nz, ebComp)];
          EZCRL = myfield[my_indice(edge, 1, 1, 2, iiC, ijR, hkL, Nx, Ny, Nz, ebComp)];
          EZRRL = myfield[my_indice(edge, 1, 1, 2, iiR, ijR, hkL, Nx, Ny, Nz, ebComp)];
        }

        {
          EXLLC = myfield[my_indice(edge, 1, 1, 0, hiL, ijL, ikC, Nx, Ny, Nz, ebComp)];
          EXCLC = myfield[my_indice(edge, 1, 1, 0, hiC, ijL, ikC, Nx, Ny, Nz, ebComp)];
          EXRLC = myfield[my_indice(edge, 1, 1, 0, hiR, ijL, ikC, Nx, Ny, Nz, ebComp)];
          EXLCC = myfield[my_indice(edge, 1, 1, 0, hiL, ijC, ikC, Nx, Ny, Nz, ebComp)];
          EXCCC = myfield[my_indice(edge, 1, 1, 0, hiC, ijC, ikC, Nx, Ny, Nz, ebComp)];
          EXRCC = myfield[my_indice(edge, 1, 1, 0, hiR, ijC, ikC, Nx, Ny, Nz, ebComp)];
          EXLRC = myfield[my_indice(edge, 1, 1, 0, hiL, ijR, ikC, Nx, Ny, Nz, ebComp)];
          EXCRC = myfield[my_indice(edge, 1, 1, 0, hiC, ijR, ikC, Nx, Ny, Nz, ebComp)];
          EXRRC = myfield[my_indice(edge, 1, 1, 0, hiR, ijR, ikC, Nx, Ny, Nz, ebComp)];

          EYLLC = myfield[my_indice(edge, 1, 1, 1, iiL, hjL, ikC, Nx, Ny, Nz, ebComp)];
          EYCLC = myfield[my_indice(edge, 1, 1, 1, iiC, hjL, ikC, Nx, Ny, Nz, ebComp)];
          EYRLC = myfield[my_indice(edge, 1, 1, 1, iiR, hjL, ikC, Nx, Ny, Nz, ebComp)];
          EYLCC = myfield[my_indice(edge, 1, 1, 1, iiL, hjC, ikC, Nx, Ny, Nz, ebComp)];
          EYCCC = myfield[my_indice(edge, 1, 1, 1, iiC, hjC, ikC, Nx, Ny, Nz, ebComp)];
          EYRCC = myfield[my_indice(edge, 1, 1, 1, iiR, hjC, ikC, Nx, Ny, Nz, ebComp)];
          EYLRC = myfield[my_indice(edge, 1, 1, 1, iiL, hjR, ikC, Nx, Ny, Nz, ebComp)];
          EYCRC = myfield[my_indice(edge, 1, 1, 1, iiC, hjR, ikC, Nx, Ny, Nz, ebComp)];
          EYRRC = myfield[my_indice(edge, 1, 1, 1, iiR, hjR, ikC, Nx, Ny, Nz, ebComp)];

          EZLLC = myfield[my_indice(edge, 1, 1, 2, iiL, ijL, hkC, Nx, Ny, Nz, ebComp)];
          EZCLC = myfield[my_indice(edge, 1, 1, 2, iiC, ijL, hkC, Nx, Ny, Nz, ebComp)];
          EZRLC = myfield[my_indice(edge, 1, 1, 2, iiR, ijL, hkC, Nx, Ny, Nz, ebComp)];
          EZLCC = myfield[my_indice(edge, 1, 1, 2, iiL, ijC, hkC, Nx, Ny, Nz, ebComp)];
          EZCCC = myfield[my_indice(edge, 1, 1, 2, iiC, ijC, hkC, Nx, Ny, Nz, ebComp)];
          EZRCC = myfield[my_indice(edge, 1, 1, 2, iiR, ijC, hkC, Nx, Ny, Nz, ebComp)];
          EZLRC = myfield[my_indice(edge, 1, 1, 2, iiL, ijR, hkC, Nx, Ny, Nz, ebComp)];
          EZCRC = myfield[my_indice(edge, 1, 1, 2, iiC, ijR, hkC, Nx, Ny, Nz, ebComp)];
          EZRRC = myfield[my_indice(edge, 1, 1, 2, iiR, ijR, hkC, Nx, Ny, Nz, ebComp)];
        }

        {
          EXLLR = myfield[my_indice(edge, 1, 1, 0, hiL, ijL, ikR, Nx, Ny, Nz, ebComp)];
          EXCLR = myfield[my_indice(edge, 1, 1, 0, hiC, ijL, ikR, Nx, Ny, Nz, ebComp)];
          EXRLR = myfield[my_indice(edge, 1, 1, 0, hiR, ijL, ikR, Nx, Ny, Nz, ebComp)];
          EXLCR = myfield[my_indice(edge, 1, 1, 0, hiL, ijC, ikR, Nx, Ny, Nz, ebComp)];
          EXCCR = myfield[my_indice(edge, 1, 1, 0, hiC, ijC, ikR, Nx, Ny, Nz, ebComp)];
          EXRCR = myfield[my_indice(edge, 1, 1, 0, hiR, ijC, ikR, Nx, Ny, Nz, ebComp)];
          EXLRR = myfield[my_indice(edge, 1, 1, 0, hiL, ijR, ikR, Nx, Ny, Nz, ebComp)];
          EXCRR = myfield[my_indice(edge, 1, 1, 0, hiC, ijR, ikR, Nx, Ny, Nz, ebComp)];
          EXRRR = myfield[my_indice(edge, 1, 1, 0, hiR, ijR, ikR, Nx, Ny, Nz, ebComp)];

          EYLLR = myfield[my_indice(edge, 1, 1, 1, iiL, hjL, ikR, Nx, Ny, Nz, ebComp)];
          EYCLR = myfield[my_indice(edge, 1, 1, 1, iiC, hjL, ikR, Nx, Ny, Nz, ebComp)];
          EYRLR = myfield[my_indice(edge, 1, 1, 1, iiR, hjL, ikR, Nx, Ny, Nz, ebComp)];
          EYLCR = myfield[my_indice(edge, 1, 1, 1, iiL, hjC, ikR, Nx, Ny, Nz, ebComp)];
          EYCCR = myfield[my_indice(edge, 1, 1, 1, iiC, hjC, ikR, Nx, Ny, Nz, ebComp)];
          EYRCR = myfield[my_indice(edge, 1, 1, 1, iiR, hjC, ikR, Nx, Ny, Nz, ebComp)];
          EYLRR = myfield[my_indice(edge, 1, 1, 1, iiL, hjR, ikR, Nx, Ny, Nz, ebComp)];
          EYCRR = myfield[my_indice(edge, 1, 1, 1, iiC, hjR, ikR, Nx, Ny, Nz, ebComp)];
          EYRRR = myfield[my_indice(edge, 1, 1, 1, iiR, hjR, ikR, Nx, Ny, Nz, ebComp)];

          EZLLR = myfield[my_indice(edge, 1, 1, 2, iiL, ijL, hkR, Nx, Ny, Nz, ebComp)];
          EZCLR = myfield[my_indice(edge, 1, 1, 2, iiC, ijL, hkR, Nx, Ny, Nz, ebComp)];
          EZRLR = myfield[my_indice(edge, 1, 1, 2, iiR, ijL, hkR, Nx, Ny, Nz, ebComp)];
          EZLCR = myfield[my_indice(edge, 1, 1, 2, iiL, ijC, hkR, Nx, Ny, Nz, ebComp)];
          EZCCR = myfield[my_indice(edge, 1, 1, 2, iiC, ijC, hkR, Nx, Ny, Nz, ebComp)];
          EZRCR = myfield[my_indice(edge, 1, 1, 2, iiR, ijC, hkR, Nx, Ny, Nz, ebComp)];
          EZLRR = myfield[my_indice(edge, 1, 1, 2, iiL, ijR, hkR, Nx, Ny, Nz, ebComp)];
          EZCRR = myfield[my_indice(edge, 1, 1, 2, iiC, ijR, hkR, Nx, Ny, Nz, ebComp)];
          EZRRR = myfield[my_indice(edge, 1, 1, 2, iiR, ijR, hkR, Nx, Ny, Nz, ebComp)];
        }

        {
          BXLLL = myfield[my_indice(edge, 1, 1, 3, iiL, hjL, hkL, Nx, Ny, Nz, ebComp)];
          BXCLL = myfield[my_indice(edge, 1, 1, 3, iiC, hjL, hkL, Nx, Ny, Nz, ebComp)];
          BXRLL = myfield[my_indice(edge, 1, 1, 3, iiR, hjL, hkL, Nx, Ny, Nz, ebComp)];
          BXLCL = myfield[my_indice(edge, 1, 1, 3, iiL, hjC, hkL, Nx, Ny, Nz, ebComp)];
          BXCCL = myfield[my_indice(edge, 1, 1, 3, iiC, hjC, hkL, Nx, Ny, Nz, ebComp)];
          BXRCL = myfield[my_indice(edge, 1, 1, 3, iiR, hjC, hkL, Nx, Ny, Nz, ebComp)];
          BXLRL = myfield[my_indice(edge, 1, 1, 3, iiL, hjR, hkL, Nx, Ny, Nz, ebComp)];
          BXCRL = myfield[my_indice(edge, 1, 1, 3, iiC, hjR, hkL, Nx, Ny, Nz, ebComp)];
          BXRRL = myfield[my_indice(edge, 1, 1, 3, iiR, hjR, hkL, Nx, Ny, Nz, ebComp)];

          BYLLL = myfield[my_indice(edge, 1, 1, 4, hiL, ijL, hkL, Nx, Ny, Nz, ebComp)];
          BYCLL = myfield[my_indice(edge, 1, 1, 4, hiC, ijL, hkL, Nx, Ny, Nz, ebComp)];
          BYRLL = myfield[my_indice(edge, 1, 1, 4, hiR, ijL, hkL, Nx, Ny, Nz, ebComp)];
          BYLCL = myfield[my_indice(edge, 1, 1, 4, hiL, ijC, hkL, Nx, Ny, Nz, ebComp)];
          BYCCL = myfield[my_indice(edge, 1, 1, 4, hiC, ijC, hkL, Nx, Ny, Nz, ebComp)];
          BYRCL = myfield[my_indice(edge, 1, 1, 4, hiR, ijC, hkL, Nx, Ny, Nz, ebComp)];
          BYLRL = myfield[my_indice(edge, 1, 1, 4, hiL, ijR, hkL, Nx, Ny, Nz, ebComp)];
          BYCRL = myfield[my_indice(edge, 1, 1, 4, hiC, ijR, hkL, Nx, Ny, Nz, ebComp)];
          BYRRL = myfield[my_indice(edge, 1, 1, 4, hiR, ijR, hkL, Nx, Ny, Nz, ebComp)];

          BZLLL = myfield[my_indice(edge, 1, 1, 5, hiL, hjL, ikL, Nx, Ny, Nz, ebComp)];
          BZCLL = myfield[my_indice(edge, 1, 1, 5, hiC, hjL, ikL, Nx, Ny, Nz, ebComp)];
          BZRLL = myfield[my_indice(edge, 1, 1, 5, hiR, hjL, ikL, Nx, Ny, Nz, ebComp)];
          BZLCL = myfield[my_indice(edge, 1, 1, 5, hiL, hjC, ikL, Nx, Ny, Nz, ebComp)];
          BZCCL = myfield[my_indice(edge, 1, 1, 5, hiC, hjC, ikL, Nx, Ny, Nz, ebComp)];
          BZRCL = myfield[my_indice(edge, 1, 1, 5, hiR, hjC, ikL, Nx, Ny, Nz, ebComp)];
          BZLRL = myfield[my_indice(edge, 1, 1, 5, hiL, hjR, ikL, Nx, Ny, Nz, ebComp)];
          BZCRL = myfield[my_indice(edge, 1, 1, 5, hiC, hjR, ikL, Nx, Ny, Nz, ebComp)];
          BZRRL = myfield[my_indice(edge, 1, 1, 5, hiR, hjR, ikL, Nx, Ny, Nz, ebComp)];
        }

        {
          BXLLC = myfield[my_indice(edge, 1, 1, 3, iiL, hjL, hkC, Nx, Ny, Nz, ebComp)];
          BXCLC = myfield[my_indice(edge, 1, 1, 3, iiC, hjL, hkC, Nx, Ny, Nz, ebComp)];
          BXRLC = myfield[my_indice(edge, 1, 1, 3, iiR, hjL, hkC, Nx, Ny, Nz, ebComp)];
          BXLCC = myfield[my_indice(edge, 1, 1, 3, iiL, hjC, hkC, Nx, Ny, Nz, ebComp)];
          BXCCC = myfield[my_indice(edge, 1, 1, 3, iiC, hjC, hkC, Nx, Ny, Nz, ebComp)];
          BXRCC = myfield[my_indice(edge, 1, 1, 3, iiR, hjC, hkC, Nx, Ny, Nz, ebComp)];
          BXLRC = myfield[my_indice(edge, 1, 1, 3, iiL, hjR, hkC, Nx, Ny, Nz, ebComp)];
          BXCRC = myfield[my_indice(edge, 1, 1, 3, iiC, hjR, hkC, Nx, Ny, Nz, ebComp)];
          BXRRC = myfield[my_indice(edge, 1, 1, 3, iiR, hjR, hkC, Nx, Ny, Nz, ebComp)];

          BYLLC = myfield[my_indice(edge, 1, 1, 4, hiL, ijL, hkC, Nx, Ny, Nz, ebComp)];
          BYCLC = myfield[my_indice(edge, 1, 1, 4, hiC, ijL, hkC, Nx, Ny, Nz, ebComp)];
          BYRLC = myfield[my_indice(edge, 1, 1, 4, hiR, ijL, hkC, Nx, Ny, Nz, ebComp)];
          BYLCC = myfield[my_indice(edge, 1, 1, 4, hiL, ijC, hkC, Nx, Ny, Nz, ebComp)];
          BYCCC = myfield[my_indice(edge, 1, 1, 4, hiC, ijC, hkC, Nx, Ny, Nz, ebComp)];
          BYRCC = myfield[my_indice(edge, 1, 1, 4, hiR, ijC, hkC, Nx, Ny, Nz, ebComp)];
          BYLRC = myfield[my_indice(edge, 1, 1, 4, hiL, ijR, hkC, Nx, Ny, Nz, ebComp)];
          BYCRC = myfield[my_indice(edge, 1, 1, 4, hiC, ijR, hkC, Nx, Ny, Nz, ebComp)];
          BYRRC = myfield[my_indice(edge, 1, 1, 4, hiR, ijR, hkC, Nx, Ny, Nz, ebComp)];

          BZLLC = myfield[my_indice(edge, 1, 1, 5, hiL, hjL, ikC, Nx, Ny, Nz, ebComp)];
          BZCLC = myfield[my_indice(edge, 1, 1, 5, hiC, hjL, ikC, Nx, Ny, Nz, ebComp)];
          BZRLC = myfield[my_indice(edge, 1, 1, 5, hiR, hjL, ikC, Nx, Ny, Nz, ebComp)];
          BZLCC = myfield[my_indice(edge, 1, 1, 5, hiL, hjC, ikC, Nx, Ny, Nz, ebComp)];
          BZCCC = myfield[my_indice(edge, 1, 1, 5, hiC, hjC, ikC, Nx, Ny, Nz, ebComp)];
          BZRCC = myfield[my_indice(edge, 1, 1, 5, hiR, hjC, ikC, Nx, Ny, Nz, ebComp)];
          BZLRC = myfield[my_indice(edge, 1, 1, 5, hiL, hjR, ikC, Nx, Ny, Nz, ebComp)];
          BZCRC = myfield[my_indice(edge, 1, 1, 5, hiC, hjR, ikC, Nx, Ny, Nz, ebComp)];
          BZRRC = myfield[my_indice(edge, 1, 1, 5, hiR, hjR, ikC, Nx, Ny, Nz, ebComp)];
        }

        {
          BXLLR = myfield[my_indice(edge, 1, 1, 3, iiL, hjL, hkR, Nx, Ny, Nz, ebComp)];
          BXCLR = myfield[my_indice(edge, 1, 1, 3, iiC, hjL, hkR, Nx, Ny, Nz, ebComp)];
          BXRLR = myfield[my_indice(edge, 1, 1, 3, iiR, hjL, hkR, Nx, Ny, Nz, ebComp)];
          BXLCR = myfield[my_indice(edge, 1, 1, 3, iiL, hjC, hkR, Nx, Ny, Nz, ebComp)];
          BXCCR = myfield[my_indice(edge, 1, 1, 3, iiC, hjC, hkR, Nx, Ny, Nz, ebComp)];
          BXRCR = myfield[my_indice(edge, 1, 1, 3, iiR, hjC, hkR, Nx, Ny, Nz, ebComp)];
          BXLRR = myfield[my_indice(edge, 1, 1, 3, iiL, hjR, hkR, Nx, Ny, Nz, ebComp)];
          BXCRR = myfield[my_indice(edge, 1, 1, 3, iiC, hjR, hkR, Nx, Ny, Nz, ebComp)];
          BXRRR = myfield[my_indice(edge, 1, 1, 3, iiR, hjR, hkR, Nx, Ny, Nz, ebComp)];

          BYLLR = myfield[my_indice(edge, 1, 1, 4, hiL, ijL, hkR, Nx, Ny, Nz, ebComp)];
          BYCLR = myfield[my_indice(edge, 1, 1, 4, hiC, ijL, hkR, Nx, Ny, Nz, ebComp)];
          BYRLR = myfield[my_indice(edge, 1, 1, 4, hiR, ijL, hkR, Nx, Ny, Nz, ebComp)];
          BYLCR = myfield[my_indice(edge, 1, 1, 4, hiL, ijC, hkR, Nx, Ny, Nz, ebComp)];
          BYCCR = myfield[my_indice(edge, 1, 1, 4, hiC, ijC, hkR, Nx, Ny, Nz, ebComp)];
          BYRCR = myfield[my_indice(edge, 1, 1, 4, hiR, ijC, hkR, Nx, Ny, Nz, ebComp)];
          BYLRR = myfield[my_indice(edge, 1, 1, 4, hiL, ijR, hkR, Nx, Ny, Nz, ebComp)];
          BYCRR = myfield[my_indice(edge, 1, 1, 4, hiC, ijR, hkR, Nx, Ny, Nz, ebComp)];
          BYRRR = myfield[my_indice(edge, 1, 1, 4, hiR, ijR, hkR, Nx, Ny, Nz, ebComp)];

          BZLLR = myfield[my_indice(edge, 1, 1, 5, hiL, hjL, ikR, Nx, Ny, Nz, ebComp)];
          BZCLR = myfield[my_indice(edge, 1, 1, 5, hiC, hjL, ikR, Nx, Ny, Nz, ebComp)];
          BZRLR = myfield[my_indice(edge, 1, 1, 5, hiR, hjL, ikR, Nx, Ny, Nz, ebComp)];
          BZLCR = myfield[my_indice(edge, 1, 1, 5, hiL, hjC, ikR, Nx, Ny, Nz, ebComp)];
          BZCCR = myfield[my_indice(edge, 1, 1, 5, hiC, hjC, ikR, Nx, Ny, Nz, ebComp)];
          BZRCR = myfield[my_indice(edge, 1, 1, 5, hiR, hjC, ikR, Nx, Ny, Nz, ebComp)];
          BZLRR = myfield[my_indice(edge, 1, 1, 5, hiL, hjR, ikR, Nx, Ny, Nz, ebComp)];
          BZCRR = myfield[my_indice(edge, 1, 1, 5, hiC, hjR, ikR, Nx, Ny, Nz, ebComp)];
          BZRRR = myfield[my_indice(edge, 1, 1, 5, hiR, hjR, ikR, Nx, Ny, Nz, ebComp)];
        }

        {
          E[0] += wiw[2][0] * ((wiw[1][0] * (hiw[0][0] * EXLLL + hiw[0][1] * EXCLL + hiw[0][2] * EXRLL))
            + (wiw[1][1] * (hiw[0][0] * EXLCL + hiw[0][1] * EXCCL + hiw[0][2] * EXRCL))
            + (wiw[1][2] * (hiw[0][0] * EXLRL + hiw[0][1] * EXCRL + hiw[0][2] * EXRRL)))

            + wiw[2][1] * ((wiw[1][0] * (hiw[0][0] * EXLLC + hiw[0][1] * EXCLC + hiw[0][2] * EXRLC))
              + (wiw[1][1] * (hiw[0][0] * EXLCC + hiw[0][1] * EXCCC + hiw[0][2] * EXRCC))
              + (wiw[1][2] * (hiw[0][0] * EXLRC + hiw[0][1] * EXCRC + hiw[0][2] * EXRRC)))

            + wiw[2][2] * ((wiw[1][0] * (hiw[0][0] * EXLLR + hiw[0][1] * EXCLR + hiw[0][2] * EXRLR))
              + (wiw[1][1] * (hiw[0][0] * EXLCR + hiw[0][1] * EXCCR + hiw[0][2] * EXRCR))
              + (wiw[1][2] * (hiw[0][0] * EXLRR + hiw[0][1] * EXCRR + hiw[0][2] * EXRRR)));


          E[1] += wiw[2][0] * ((hiw[1][0] * (wiw[0][0] * EYLLL + wiw[0][1] * EYCLL + wiw[0][2] * EYRLL))
            + (hiw[1][1] * (wiw[0][0] * EYLCL + wiw[0][1] * EYCCL + wiw[0][2] * EYRCL))
            + (hiw[1][2] * (wiw[0][0] * EYLRL + wiw[0][1] * EYCRL + wiw[0][2] * EYRRL)))

            + wiw[2][1] * ((hiw[1][0] * (wiw[0][0] * EYLLC + wiw[0][1] * EYCLC + wiw[0][2] * EYRLC))
              + (hiw[1][1] * (wiw[0][0] * EYLCC + wiw[0][1] * EYCCC + wiw[0][2] * EYRCC))
              + (hiw[1][2] * (wiw[0][0] * EYLRC + wiw[0][1] * EYCRC + wiw[0][2] * EYRRC)))

            + wiw[2][2] * ((hiw[1][0] * (wiw[0][0] * EYLLR + wiw[0][1] * EYCLR + wiw[0][2] * EYRLR))
              + (hiw[1][1] * (wiw[0][0] * EYLCR + wiw[0][1] * EYCCR + wiw[0][2] * EYRCR))
              + (hiw[1][2] * (wiw[0][0] * EYLRR + wiw[0][1] * EYCRR + wiw[0][2] * EYRRR)));


          E[2] += hiw[2][0] * ((wiw[1][0] * (wiw[0][0] * EZLLL + wiw[0][1] * EZCLL + wiw[0][2] * EZRLL))
            + (wiw[1][1] * (wiw[0][0] * EZLCL + wiw[0][1] * EZCCL + wiw[0][2] * EZRCL))
            + (wiw[1][2] * (wiw[0][0] * EZLRL + wiw[0][1] * EZCRL + wiw[0][2] * EZRRL)))

            + hiw[2][1] * ((wiw[1][0] * (wiw[0][0] * EZLLC + wiw[0][1] * EZCLC + wiw[0][2] * EZRLC))
              + (wiw[1][1] * (wiw[0][0] * EZLCC + wiw[0][1] * EZCCC + wiw[0][2] * EZRCC))
              + (wiw[1][2] * (wiw[0][0] * EZLRC + wiw[0][1] * EZCRC + wiw[0][2] * EZRRC)))

            + hiw[2][2] * ((wiw[1][0] * (wiw[0][0] * EZLLR + wiw[0][1] * EZCLR + wiw[0][2] * EZRLR))
              + (wiw[1][1] * (wiw[0][0] * EZLCR + wiw[0][1] * EZCCR + wiw[0][2] * EZRCR))
              + (wiw[1][2] * (wiw[0][0] * EZLRR + wiw[0][1] * EZCRR + wiw[0][2] * EZRRR)));


          B[0] += hiw[2][0] * ((hiw[1][0] * (wiw[0][0] * BXLLL + wiw[0][1] * BXCLL + wiw[0][2] * BXRLL))
            + (hiw[1][1] * (wiw[0][0] * BXLCL + wiw[0][1] * BXCCL + wiw[0][2] * BXRCL))
            + (hiw[1][2] * (wiw[0][0] * BXLRL + wiw[0][1] * BXCRL + wiw[0][2] * BXRRL)))

            + hiw[2][1] * ((hiw[1][0] * (wiw[0][0] * BXLLC + wiw[0][1] * BXCLC + wiw[0][2] * BXRLC))
              + (hiw[1][1] * (wiw[0][0] * BXLCC + wiw[0][1] * BXCCC + wiw[0][2] * BXRCC))
              + (hiw[1][2] * (wiw[0][0] * BXLRC + wiw[0][1] * BXCRC + wiw[0][2] * BXRRC)))

            + hiw[2][2] * ((hiw[1][0] * (wiw[0][0] * BXLLR + wiw[0][1] * BXCLR + wiw[0][2] * BXRLR))
              + (hiw[1][1] * (wiw[0][0] * BXLCR + wiw[0][1] * BXCCR + wiw[0][2] * BXRCR))
              + (hiw[1][2] * (wiw[0][0] * BXLRR + wiw[0][1] * BXCRR + wiw[0][2] * BXRRR)));


          B[1] += hiw[2][0] * ((wiw[1][0] * (hiw[0][0] * BYLLL + hiw[0][1] * BYCLL + hiw[0][2] * BYRLL))
            + (wiw[1][1] * (hiw[0][0] * BYLCL + hiw[0][1] * BYCCL + hiw[0][2] * BYRCL))
            + (wiw[1][2] * (hiw[0][0] * BYLRL + hiw[0][1] * BYCRL + hiw[0][2] * BYRRL)))

            + hiw[2][1] * ((wiw[1][0] * (hiw[0][0] * BYLLC + hiw[0][1] * BYCLC + hiw[0][2] * BYRLC))
              + (wiw[1][1] * (hiw[0][0] * BYLCC + hiw[0][1] * BYCCC + hiw[0][2] * BYRCC))
              + (wiw[1][2] * (hiw[0][0] * BYLRC + hiw[0][1] * BYCRC + hiw[0][2] * BYRRC)))

            + hiw[2][2] * ((wiw[1][0] * (hiw[0][0] * BYLLR + hiw[0][1] * BYCLR + hiw[0][2] * BYRLR))
              + (wiw[1][1] * (hiw[0][0] * BYLCR + hiw[0][1] * BYCCR + hiw[0][2] * BYRCR))
              + (wiw[1][2] * (hiw[0][0] * BYLRR + hiw[0][1] * BYCRR + hiw[0][2] * BYRRR)));


          B[2] += wiw[2][0] * ((hiw[1][0] * (hiw[0][0] * BZLLL + hiw[0][1] * BZCLL + hiw[0][2] * BZRLL))
            + (hiw[1][1] * (hiw[0][0] * BZLCL + hiw[0][1] * BZCCL + hiw[0][2] * BZRCL))
            + (hiw[1][2] * (hiw[0][0] * BZLRL + hiw[0][1] * BZCRL + hiw[0][2] * BZRRL)))

            + wiw[2][1] * ((hiw[1][0] * (hiw[0][0] * BZLLC + hiw[0][1] * BZCLC + hiw[0][2] * BZRLC))
              + (hiw[1][1] * (hiw[0][0] * BZLCC + hiw[0][1] * BZCCC + hiw[0][2] * BZRCC))
              + (hiw[1][2] * (hiw[0][0] * BZLRC + hiw[0][1] * BZCRC + hiw[0][2] * BZRRC)))

            + wiw[2][2] * ((hiw[1][0] * (hiw[0][0] * BZLLR + hiw[0][1] * BZCLR + hiw[0][2] * BZRLR))
              + (hiw[1][1] * (hiw[0][0] * BZLCR + hiw[0][1] * BZCCR + hiw[0][2] * BZRCR))
              + (hiw[1][2] * (hiw[0][0] * BZLRR + hiw[0][1] * BZCRR + hiw[0][2] * BZRRR)));
        }

      }

#endif
#ifdef _ACC_SINGLE_POINTER
      u_minus[0] = pData[3 + p*Ncomp] + 0.5*dt*coupling*E[0];
      u_minus[1] = pData[4 + p*Ncomp] + 0.5*dt*coupling*E[1];
      u_minus[2] = pData[5 + p*Ncomp] + 0.5*dt*coupling*E[2];
#else
      u_minus[0] = val[3][p*Ncomp] + 0.5*dt*coupling*E[0];
      u_minus[1] = val[4][p*Ncomp] + 0.5*dt*coupling*E[1];
      u_minus[2] = val[5][p*Ncomp] + 0.5*dt*coupling*E[2];
#endif
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

#ifdef _ACC_SINGLE_POINTER
      pData[3 + p*Ncomp] = (u_plus[0] + 0.5*dt*coupling*E[0]);
      pData[4 + p*Ncomp] = (u_plus[1] + 0.5*dt*coupling*E[1]);
      pData[5 + p*Ncomp] = (u_plus[2] + 0.5*dt*coupling*E[2]);
#else
      val[3][p*Ncomp] = (u_plus[0] + 0.5*dt*coupling*E[0]);
      val[4][p*Ncomp] = (u_plus[1] + 0.5*dt*coupling*E[1]);
      val[5][p*Ncomp] = (u_plus[2] + 0.5*dt*coupling*E[2]);
#endif
    }
    break;

  case 2:
    for (p = 0; p < Np; p++)
    {
      //gamma_i=1./sqrt(1+u0(p)*u0(p)+u1(p)*u1(p)+u2(p)*u2(p));
      for (c = 0; c < 3; c++)
      {
#ifdef _ACC_SINGLE_POINTER
        xx[c] = pData[c + p*Ncomp];
#else
        xx[c] = val[c][p*Ncomp];
#endif
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

#ifdef OLDPUSHER
      k1 = k2 = 0;
      for (j = 0; j < 3; j++)
      {
        j1 = j + wii[1] - 1;
        j2 = j + hii[1] - 1;
        for (i = 0; i < 3; i++) {
          i1 = i + wii[0] - 1;
          i2 = i + hii[0] - 1;
#ifndef OLD_ACCESS
          double EX, EY, EZ;
          EX = myfield[my_indice(edge, 1, 0, 0, i2, j1, k1, N_grid[0], N_grid[1], N_grid[2], ebComp)];
          EY = myfield[my_indice(edge, 1, 0, 1, i1, j2, k1, N_grid[0], N_grid[1], N_grid[2], ebComp)];
          EZ = myfield[my_indice(edge, 1, 0, 2, i1, j1, k2, N_grid[0], N_grid[1], N_grid[2], ebComp)];
          double BX, BY, BZ;
          BX = myfield[my_indice(edge, 1, 0, 3, i1, j2, k2, N_grid[0], N_grid[1], N_grid[2], ebComp)];
          BY = myfield[my_indice(edge, 1, 0, 4, i2, j1, k2, N_grid[0], N_grid[1], N_grid[2], ebComp)];
          BZ = myfield[my_indice(edge, 1, 0, 5, i2, j2, k1, N_grid[0], N_grid[1], N_grid[2], ebComp)];

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
#else
      {

        double EXLLC, EXCLC, EXRLC, EXLCC, EXCCC, EXRCC, EXLRC, EXCRC, EXRRC;
        double EYLLC, EYCLC, EYRLC, EYLCC, EYCCC, EYRCC, EYLRC, EYCRC, EYRRC;
        double EZLLC, EZCLC, EZRLC, EZLCC, EZCCC, EZRCC, EZLRC, EZCRC, EZRRC;


        double BXLLC, BXCLC, BXRLC, BXLCC, BXCCC, BXRCC, BXLRC, BXCRC, BXRRC;
        double BYLLC, BYCLC, BYRLC, BYLCC, BYCCC, BYRCC, BYLRC, BYCRC, BYRRC;
        double BZLLC, BZCLC, BZRLC, BZLCC, BZCCC, BZRCC, BZLRC, BZCRC, BZRRC;

        int iiL, iiC, iiR, hiL, hiC, hiR;
        iiL = wii[0] - 1;
        iiC = wii[0];
        iiR = wii[0] + 1;
        hiL = hii[0] - 1;
        hiC = hii[0];
        hiR = hii[0] + 1;

        int ijL, ijC, ijR, hjL, hjC, hjR;
        ijL = wii[1] - 1;
        ijC = wii[1];
        ijR = wii[1] + 1;
        hjL = hii[1] - 1;
        hjC = hii[1];
        hjR = hii[1] + 1;

        int ikC, hkC;
        ikC = 0;
        hkC = 0;


        {
          EXLLC = myfield[my_indice(edge, 1, 0, 0, hiL, ijL, ikC, Nx, Ny, Nz, ebComp)];
          EXCLC = myfield[my_indice(edge, 1, 0, 0, hiC, ijL, ikC, Nx, Ny, Nz, ebComp)];
          EXRLC = myfield[my_indice(edge, 1, 0, 0, hiR, ijL, ikC, Nx, Ny, Nz, ebComp)];
          EXLCC = myfield[my_indice(edge, 1, 0, 0, hiL, ijC, ikC, Nx, Ny, Nz, ebComp)];
          EXCCC = myfield[my_indice(edge, 1, 0, 0, hiC, ijC, ikC, Nx, Ny, Nz, ebComp)];
          EXRCC = myfield[my_indice(edge, 1, 0, 0, hiR, ijC, ikC, Nx, Ny, Nz, ebComp)];
          EXLRC = myfield[my_indice(edge, 1, 0, 0, hiL, ijR, ikC, Nx, Ny, Nz, ebComp)];
          EXCRC = myfield[my_indice(edge, 1, 0, 0, hiC, ijR, ikC, Nx, Ny, Nz, ebComp)];
          EXRRC = myfield[my_indice(edge, 1, 0, 0, hiR, ijR, ikC, Nx, Ny, Nz, ebComp)];

          EYLLC = myfield[my_indice(edge, 1, 0, 1, iiL, hjL, ikC, Nx, Ny, Nz, ebComp)];
          EYCLC = myfield[my_indice(edge, 1, 0, 1, iiC, hjL, ikC, Nx, Ny, Nz, ebComp)];
          EYRLC = myfield[my_indice(edge, 1, 0, 1, iiR, hjL, ikC, Nx, Ny, Nz, ebComp)];
          EYLCC = myfield[my_indice(edge, 1, 0, 1, iiL, hjC, ikC, Nx, Ny, Nz, ebComp)];
          EYCCC = myfield[my_indice(edge, 1, 0, 1, iiC, hjC, ikC, Nx, Ny, Nz, ebComp)];
          EYRCC = myfield[my_indice(edge, 1, 0, 1, iiR, hjC, ikC, Nx, Ny, Nz, ebComp)];
          EYLRC = myfield[my_indice(edge, 1, 0, 1, iiL, hjR, ikC, Nx, Ny, Nz, ebComp)];
          EYCRC = myfield[my_indice(edge, 1, 0, 1, iiC, hjR, ikC, Nx, Ny, Nz, ebComp)];
          EYRRC = myfield[my_indice(edge, 1, 0, 1, iiR, hjR, ikC, Nx, Ny, Nz, ebComp)];

          EZLLC = myfield[my_indice(edge, 1, 0, 2, iiL, ijL, hkC, Nx, Ny, Nz, ebComp)];
          EZCLC = myfield[my_indice(edge, 1, 0, 2, iiC, ijL, hkC, Nx, Ny, Nz, ebComp)];
          EZRLC = myfield[my_indice(edge, 1, 0, 2, iiR, ijL, hkC, Nx, Ny, Nz, ebComp)];
          EZLCC = myfield[my_indice(edge, 1, 0, 2, iiL, ijC, hkC, Nx, Ny, Nz, ebComp)];
          EZCCC = myfield[my_indice(edge, 1, 0, 2, iiC, ijC, hkC, Nx, Ny, Nz, ebComp)];
          EZRCC = myfield[my_indice(edge, 1, 0, 2, iiR, ijC, hkC, Nx, Ny, Nz, ebComp)];
          EZLRC = myfield[my_indice(edge, 1, 0, 2, iiL, ijR, hkC, Nx, Ny, Nz, ebComp)];
          EZCRC = myfield[my_indice(edge, 1, 0, 2, iiC, ijR, hkC, Nx, Ny, Nz, ebComp)];
          EZRRC = myfield[my_indice(edge, 1, 0, 2, iiR, ijR, hkC, Nx, Ny, Nz, ebComp)];
        }


        {
          BXLLC = myfield[my_indice(edge, 1, 0, 3, iiL, hjL, hkC, Nx, Ny, Nz, ebComp)];
          BXCLC = myfield[my_indice(edge, 1, 0, 3, iiC, hjL, hkC, Nx, Ny, Nz, ebComp)];
          BXRLC = myfield[my_indice(edge, 1, 0, 3, iiR, hjL, hkC, Nx, Ny, Nz, ebComp)];
          BXLCC = myfield[my_indice(edge, 1, 0, 3, iiL, hjC, hkC, Nx, Ny, Nz, ebComp)];
          BXCCC = myfield[my_indice(edge, 1, 0, 3, iiC, hjC, hkC, Nx, Ny, Nz, ebComp)];
          BXRCC = myfield[my_indice(edge, 1, 0, 3, iiR, hjC, hkC, Nx, Ny, Nz, ebComp)];
          BXLRC = myfield[my_indice(edge, 1, 0, 3, iiL, hjR, hkC, Nx, Ny, Nz, ebComp)];
          BXCRC = myfield[my_indice(edge, 1, 0, 3, iiC, hjR, hkC, Nx, Ny, Nz, ebComp)];
          BXRRC = myfield[my_indice(edge, 1, 0, 3, iiR, hjR, hkC, Nx, Ny, Nz, ebComp)];

          BYLLC = myfield[my_indice(edge, 1, 0, 4, hiL, ijL, hkC, Nx, Ny, Nz, ebComp)];
          BYCLC = myfield[my_indice(edge, 1, 0, 4, hiC, ijL, hkC, Nx, Ny, Nz, ebComp)];
          BYRLC = myfield[my_indice(edge, 1, 0, 4, hiR, ijL, hkC, Nx, Ny, Nz, ebComp)];
          BYLCC = myfield[my_indice(edge, 1, 0, 4, hiL, ijC, hkC, Nx, Ny, Nz, ebComp)];
          BYCCC = myfield[my_indice(edge, 1, 0, 4, hiC, ijC, hkC, Nx, Ny, Nz, ebComp)];
          BYRCC = myfield[my_indice(edge, 1, 0, 4, hiR, ijC, hkC, Nx, Ny, Nz, ebComp)];
          BYLRC = myfield[my_indice(edge, 1, 0, 4, hiL, ijR, hkC, Nx, Ny, Nz, ebComp)];
          BYCRC = myfield[my_indice(edge, 1, 0, 4, hiC, ijR, hkC, Nx, Ny, Nz, ebComp)];
          BYRRC = myfield[my_indice(edge, 1, 0, 4, hiR, ijR, hkC, Nx, Ny, Nz, ebComp)];

          BZLLC = myfield[my_indice(edge, 1, 0, 5, hiL, hjL, ikC, Nx, Ny, Nz, ebComp)];
          BZCLC = myfield[my_indice(edge, 1, 0, 5, hiC, hjL, ikC, Nx, Ny, Nz, ebComp)];
          BZRLC = myfield[my_indice(edge, 1, 0, 5, hiR, hjL, ikC, Nx, Ny, Nz, ebComp)];
          BZLCC = myfield[my_indice(edge, 1, 0, 5, hiL, hjC, ikC, Nx, Ny, Nz, ebComp)];
          BZCCC = myfield[my_indice(edge, 1, 0, 5, hiC, hjC, ikC, Nx, Ny, Nz, ebComp)];
          BZRCC = myfield[my_indice(edge, 1, 0, 5, hiR, hjC, ikC, Nx, Ny, Nz, ebComp)];
          BZLRC = myfield[my_indice(edge, 1, 0, 5, hiL, hjR, ikC, Nx, Ny, Nz, ebComp)];
          BZCRC = myfield[my_indice(edge, 1, 0, 5, hiC, hjR, ikC, Nx, Ny, Nz, ebComp)];
          BZRRC = myfield[my_indice(edge, 1, 0, 5, hiR, hjR, ikC, Nx, Ny, Nz, ebComp)];
        }

        {
          E[0] += ((wiw[1][0] * (hiw[0][0] * EXLLC + hiw[0][1] * EXCLC + hiw[0][2] * EXRLC))
            + (wiw[1][1] * (hiw[0][0] * EXLCC + hiw[0][1] * EXCCC + hiw[0][2] * EXRCC))
            + (wiw[1][2] * (hiw[0][0] * EXLRC + hiw[0][1] * EXCRC + hiw[0][2] * EXRRC)));

          E[1] += ((hiw[1][0] * (wiw[0][0] * EYLLC + wiw[0][1] * EYCLC + wiw[0][2] * EYRLC))
            + (hiw[1][1] * (wiw[0][0] * EYLCC + wiw[0][1] * EYCCC + wiw[0][2] * EYRCC))
            + (hiw[1][2] * (wiw[0][0] * EYLRC + wiw[0][1] * EYCRC + wiw[0][2] * EYRRC)));

          E[2] += ((wiw[1][0] * (wiw[0][0] * EZLLC + wiw[0][1] * EZCLC + wiw[0][2] * EZRLC))
            + (wiw[1][1] * (wiw[0][0] * EZLCC + wiw[0][1] * EZCCC + wiw[0][2] * EZRCC))
            + (wiw[1][2] * (wiw[0][0] * EZLRC + wiw[0][1] * EZCRC + wiw[0][2] * EZRRC)));


          B[0] += ((hiw[1][0] * (wiw[0][0] * BXLLC + wiw[0][1] * BXCLC + wiw[0][2] * BXRLC))
            + (hiw[1][1] * (wiw[0][0] * BXLCC + wiw[0][1] * BXCCC + wiw[0][2] * BXRCC))
            + (hiw[1][2] * (wiw[0][0] * BXLRC + wiw[0][1] * BXCRC + wiw[0][2] * BXRRC)));

          B[1] += ((wiw[1][0] * (hiw[0][0] * BYLLC + hiw[0][1] * BYCLC + hiw[0][2] * BYRLC))
            + (wiw[1][1] * (hiw[0][0] * BYLCC + hiw[0][1] * BYCCC + hiw[0][2] * BYRCC))
            + (wiw[1][2] * (hiw[0][0] * BYLRC + hiw[0][1] * BYCRC + hiw[0][2] * BYRRC)));

          B[2] += ((hiw[1][0] * (hiw[0][0] * BZLLC + hiw[0][1] * BZCLC + hiw[0][2] * BZRLC))
            + (hiw[1][1] * (hiw[0][0] * BZLCC + hiw[0][1] * BZCCC + hiw[0][2] * BZRCC))
            + (hiw[1][2] * (hiw[0][0] * BZLRC + hiw[0][1] * BZCRC + hiw[0][2] * BZRRC)));
        }

      }

#endif
#ifdef _ACC_SINGLE_POINTER
      u_minus[0] = pData[3 + p*Ncomp] + 0.5*dt*coupling*E[0];
      u_minus[1] = pData[4 + p*Ncomp] + 0.5*dt*coupling*E[1];
      u_minus[2] = pData[5 + p*Ncomp] + 0.5*dt*coupling*E[2];
#else
      u_minus[0] = val[3][p*Ncomp] + 0.5*dt*coupling*E[0];
      u_minus[1] = val[4][p*Ncomp] + 0.5*dt*coupling*E[1];
      u_minus[2] = val[5][p*Ncomp] + 0.5*dt*coupling*E[2];
#endif

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

#ifdef _ACC_SINGLE_POINTER
      pData[3 + p*Ncomp] = (u_plus[0] + 0.5*dt*coupling*E[0]);
      pData[4 + p*Ncomp] = (u_plus[1] + 0.5*dt*coupling*E[1]);
      pData[5 + p*Ncomp] = (u_plus[2] + 0.5*dt*coupling*E[2]);
#else
      val[3][p*Ncomp] = (u_plus[0] + 0.5*dt*coupling*E[0]);
      val[4][p*Ncomp] = (u_plus[1] + 0.5*dt*coupling*E[1]);
      val[5][p*Ncomp] = (u_plus[2] + 0.5*dt*coupling*E[2]);
#endif
    }
    break;

  case 1:
    for (p = 0; p < Np; p++)
    {
      //gamma_i=1./sqrt(1+u0(p)*u0(p)+u1(p)*u1(p)+u2(p)*u2(p));
      for (c = 0; c < 3; c++)
      {
#ifdef _ACC_SINGLE_POINTER
        xx[c] = pData[c + p*Ncomp];
#else
        xx[c] = val[c][p*Ncomp];
#endif
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
        EX = myfield[my_indice(edge, 0, 0, 0, i2, j1, k1, N_grid[0], N_grid[1], N_grid[2], ebComp)];
        EY = myfield[my_indice(edge, 0, 0, 1, i1, j2, k1, N_grid[0], N_grid[1], N_grid[2], ebComp)];
        EZ = myfield[my_indice(edge, 0, 0, 2, i1, j1, k2, N_grid[0], N_grid[1], N_grid[2], ebComp)];
        double BX, BY, BZ;
        BX = myfield[my_indice(edge, 0, 0, 3, i1, j2, k2, N_grid[0], N_grid[1], N_grid[2], ebComp)];
        BY = myfield[my_indice(edge, 0, 0, 4, i2, j1, k2, N_grid[0], N_grid[1], N_grid[2], ebComp)];
        BZ = myfield[my_indice(edge, 0, 0, 5, i2, j2, k1, N_grid[0], N_grid[1], N_grid[2], ebComp)];

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

#ifdef _ACC_SINGLE_POINTER
      u_minus[0] = pData[3 + p*Ncomp] + 0.5*dt*coupling*E[0];
      u_minus[1] = pData[4 + p*Ncomp] + 0.5*dt*coupling*E[1];
      u_minus[2] = pData[5 + p*Ncomp] + 0.5*dt*coupling*E[2];
#else
      u_minus[0] = val[3][p*Ncomp] + 0.5*dt*coupling*E[0];
      u_minus[1] = val[4][p*Ncomp] + 0.5*dt*coupling*E[1];
      u_minus[2] = val[5][p*Ncomp] + 0.5*dt*coupling*E[2];
#endif

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

#ifdef _ACC_SINGLE_POINTER
      pData[3 + p*Ncomp] = (u_plus[0] + 0.5*dt*coupling*E[0]);
      pData[4 + p*Ncomp] = (u_plus[1] + 0.5*dt*coupling*E[1]);
      pData[5 + p*Ncomp] = (u_plus[2] + 0.5*dt*coupling*E[2]);
#else
      val[3][p*Ncomp] = (u_plus[0] + 0.5*dt*coupling*E[0]);
      val[4][p*Ncomp] = (u_plus[1] + 0.5*dt*coupling*E[1]);
      val[5][p*Ncomp] = (u_plus[2] + 0.5*dt*coupling*E[2]);
#endif
    }
    break;
  }
}


void SPECIE::momenta_advance_with_friction(EM_FIELD *ebfield, double lambda)
{
  energyExtremesFlag = false;

  if (mygrid->withParticles == NO||isFrozen)
    return;
  if (mygrid->isStretched()) {
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
  if (mygrid->withParticles == NO||isFrozen)
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
  if (mygrid->withParticles == NO||isFrozen)
    return;
  if (mygrid->withCurrent == NO||isTestSpecies)
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
      //      ru(1, p) += dt*vy;
      //      ru(2, p) += dt*vz;
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

  if (!allocated) {
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

void SPECIE::computeLorentzMatrix(double ux, double uy, double uz, double *matr) {
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

void SPECIE::callWaterbag(my_rng_generator& ext_rng, double p0_x, double p0_y, double p0_z, double uxin, double uyin, double uzin) {
  my_uniform_real_distribution dist(-1.0, 1.0);
  if (uxin*uxin + uyin*uyin + uzin*uzin < _VERY_SMALL_MOMENTUM*_VERY_SMALL_MOMENTUM) {

    for (int p = 0; p < Np; p++)
    {
      u0(p) = uxin + p0_x*dist(ext_rng);
      u1(p) = uyin + p0_y*dist(ext_rng);
      u2(p) = uzin + p0_z*dist(ext_rng);

    }
  }
  else {
    double L[16];
    computeLorentzMatrix(uxin, uyin, uzin, L);
    double Ett, u0t, u1t, u2t;
    for (int p = 0; p < Np; p++)
    {
      u0(p) = p0_x*dist(ext_rng);
      u1(p) = p0_y*dist(ext_rng);
      u2(p) = p0_z*dist(ext_rng);
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

void SPECIE::callUnifSphere(my_rng_generator& ext_rng, double p0, double uxin, double uyin, double uzin) {
  double pmod;
  double phi;
  double cos_theta, sin_theta;

  my_uniform_real_distribution distR(0, 1.0);
  my_uniform_real_distribution distA(-1.0, 1.0);
  my_uniform_real_distribution distB(0, 2.0*M_PI);

  if (uxin*uxin + uyin*uyin + uzin*uzin < _VERY_SMALL_MOMENTUM*_VERY_SMALL_MOMENTUM) {
    for (int p = 0; p < Np; p++)
    {
      pmod = pow(distR(ext_rng), 1. / 3.);
      phi = distB(ext_rng);
      cos_theta = distB(ext_rng);
      sin_theta = sqrt(1.0 - cos_theta*cos_theta);
      u0(p) = uxin + p0*pmod*sin_theta*cos(phi);
      u1(p) = uyin + p0*pmod*sin_theta*sin(phi);
      u2(p) = uzin + p0*pmod*cos_theta;
    }
  }
  else {
    double L[16];
    computeLorentzMatrix(uxin, uyin, uzin, L);
    double Ett, u0t, u1t, u2t;
    for (int p = 0; p < Np; p++)
    {
      pmod = pow(distR(ext_rng), 1. / 3.);
      phi = distB(ext_rng);
      cos_theta = distB(ext_rng);
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

void SPECIE::callSupergaussian(my_rng_generator& ext_rng, double p0, double alpha, double uxin, double uyin, double uzin) {

  my_exponential_distribution expDist(alpha);
  if (uxin*uxin + uyin*uyin + uzin*uzin < _VERY_SMALL_MOMENTUM*_VERY_SMALL_MOMENTUM) {
    for (int p = 0; p < Np; p++)
    {
      u0(p) = uxin + p0*expDist(ext_rng);
      u1(p) = uyin + p0*expDist(ext_rng);
      u2(p) = uzin + p0*expDist(ext_rng);
    }
  }
  else {
    double L[16];
    computeLorentzMatrix(uxin, uyin, uzin, L);
    double Ett, u0t, u1t, u2t;
    for (int p = 0; p < Np; p++)
    {
      u0(p) = p0*expDist(ext_rng);
      u1(p) = p0*expDist(ext_rng);
      u2(p) = p0*expDist(ext_rng);
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

//#define USE_BOOST
void SPECIE::callMaxwell(my_rng_generator& ext_rng, double Ta, double uxin, double uyin, double uzin) {

  my_normal_distribution myGaussian(0,sqrt(Ta));

#if defined(USE_BOOST)
  if(isQuiet){
    int Nshuffle = quietShuffle;
    std::srand ( unsigned ( std::time(0) ) );
    std::vector<int> myvector;

    // set some values:
    for (int i=0; i<Nshuffle; i++) myvector.push_back(i);
    int groups = Np/Nshuffle;
    int remain = Np%Nshuffle;

std::cout << "Nshuffle = " << Nshuffle << std::endl;

    boost::math::normal dist(0.0, sqrt(Ta));

    int dim_num = 6;
    int offset = 0;
    double randomU[dim_num];
    long long int seed=111111*(mygrid->myid+1);

    if ((uxin*uxin + uyin*uyin + uzin*uzin) < _VERY_SMALL_MOMENTUM*_VERY_SMALL_MOMENTUM) {
    /*  for (int p = 0; p < Np; p++){
        i8_sobol ( dim_num, &seed, randomU );
        u0(p) = uxin +quantile(dist, randomU[0+offset]);
        u1(p) = uyin +quantile(dist, randomU[1+offset]);
        u2(p) = uzin +quantile(dist, randomU[2+offset]);
      }
      */
      for(int nn=0; nn<groups;nn++){
        std::random_shuffle ( myvector.begin(), myvector.end() );
        for (std::vector<int>::iterator it=myvector.begin(); it!=myvector.end(); ++it){

          int p = nn*Nshuffle + *it;
          i8_sobol ( dim_num, &seed, randomU );
          u0(p) = uxin +quantile(dist, randomU[0+offset]);
          u1(p) = uyin +quantile(dist, randomU[1+offset]);
          u2(p) = uzin +quantile(dist, randomU[2+offset]);
        }
      }
      for (int p = groups*Nshuffle; p < Np; p++){
        i8_sobol ( dim_num, &seed, randomU );
        u0(p) = uxin +quantile(dist, randomU[0+offset]);
        u1(p) = uyin +quantile(dist, randomU[1+offset]);
        u2(p) = uzin +quantile(dist, randomU[2+offset]);
      }
    }
    else {
      std::cout<< "RELLATIVISTIC! uxin=" << uxin << std::endl;

      double L[16];
      computeLorentzMatrix(uxin, uyin, uzin, L);
      double Ett, u0t, u1t, u2t;
     /* for (int p = 0; p < Np; p++)
      {
        i8_sobol ( dim_num, &seed, randomU );
        u0(p) = quantile(dist, randomU[0+offset]);
        u1(p) = quantile(dist, randomU[1+offset]);
        u2(p) = quantile(dist, randomU[2+offset]);

        Ett = sqrt(1.0 + u0(p)*u0(p) + u1(p)*u1(p) + u2(p)*u2(p));

        u0t = L[1 * 4 + 0] * Ett + L[1 * 4 + 1] * u0(p) + L[1 * 4 + 2] * u1(p) + L[1 * 4 + 3] * u2(p);
        u1t = L[2 * 4 + 0] * Ett + L[2 * 4 + 1] * u0(p) + L[2 * 4 + 2] * u1(p) + L[2 * 4 + 3] * u2(p);
        u2t = L[3 * 4 + 0] * Ett + L[3 * 4 + 1] * u0(p) + L[3 * 4 + 2] * u1(p) + L[3 * 4 + 3] * u2(p);

        u0(p) = u0t;
        u1(p) = u1t;
        u2(p) = u2t;
      }
      */
      for(int nn=0; nn<groups;nn++){
        std::random_shuffle ( myvector.begin(), myvector.end() );
        for (std::vector<int>::iterator it=myvector.begin(); it!=myvector.end(); ++it){

          int p = nn*Nshuffle + *it;
          i8_sobol ( dim_num, &seed, randomU );
          u0(p) = quantile(dist, randomU[0+offset]);
          u1(p) = quantile(dist, randomU[1+offset]);
          u2(p) = quantile(dist, randomU[2+offset]);

          Ett = sqrt(1.0 + u0(p)*u0(p) + u1(p)*u1(p) + u2(p)*u2(p));

          u0t = L[1 * 4 + 0] * Ett + L[1 * 4 + 1] * u0(p) + L[1 * 4 + 2] * u1(p) + L[1 * 4 + 3] * u2(p);
          u1t = L[2 * 4 + 0] * Ett + L[2 * 4 + 1] * u0(p) + L[2 * 4 + 2] * u1(p) + L[2 * 4 + 3] * u2(p);
          u2t = L[3 * 4 + 0] * Ett + L[3 * 4 + 1] * u0(p) + L[3 * 4 + 2] * u1(p) + L[3 * 4 + 3] * u2(p);

          u0(p) = u0t;
          u1(p) = u1t;
          u2(p) = u2t;
        }
      }
      for (int p = groups*Nshuffle; p < Np; p++){
        i8_sobol ( dim_num, &seed, randomU );
        u0(p) = quantile(dist, randomU[0+offset]);
        u1(p) = quantile(dist, randomU[1+offset]);
        u2(p) = quantile(dist, randomU[2+offset]);

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
  else{
#endif

    if (uxin*uxin + uyin*uyin + uzin*uzin < _VERY_SMALL_MOMENTUM*_VERY_SMALL_MOMENTUM) {
      for (int p = 0; p < Np; p++)
      {
        u0(p) = uxin + myGaussian(ext_rng);
        u1(p) = uyin + myGaussian(ext_rng);
        u2(p) = uzin + myGaussian(ext_rng);
      }
    }
    else {
      double L[16];
      computeLorentzMatrix(uxin, uyin, uzin, L);
      double Ett, u0t, u1t, u2t;
      for (int p = 0; p < Np; p++)
      {
        u0(p) = myGaussian(ext_rng);
        u1(p) = myGaussian(ext_rng);
        u2(p) = myGaussian(ext_rng);

        Ett = sqrt(1.0 + u0(p)*u0(p) + u1(p)*u1(p) + u2(p)*u2(p));

        u0t = L[1 * 4 + 0] * Ett + L[1 * 4 + 1] * u0(p) + L[1 * 4 + 2] * u1(p) + L[1 * 4 + 3] * u2(p);
        u1t = L[2 * 4 + 0] * Ett + L[2 * 4 + 1] * u0(p) + L[2 * 4 + 2] * u1(p) + L[2 * 4 + 3] * u2(p);
        u2t = L[3 * 4 + 0] * Ett + L[3 * 4 + 1] * u0(p) + L[3 * 4 + 2] * u1(p) + L[3 * 4 + 3] * u2(p);

        u0(p) = u0t;
        u1(p) = u1t;
        u2(p) = u2t;
      }
    }
#if defined(USE_BOOST)
  }
#endif
}
void SPECIE::callJuttner(my_rng_generator& ext_rng, double Ta, double uxin, double uyin, double uzin) {
  //DA DEFINIRE
}
double densityFunctionMaxwell(double px, double alpha, double temp) {
  return exp(-(sqrt(alpha*alpha + px*px) - alpha) / temp);
}
void SPECIE::callSpecial(my_rng_generator& ext_rng, double Ta) {
  //double ptot, temp, cos_theta, sin_theta, segno, phi;
  double alpha;
  double auxPX, auxDF;
  my_uniform_real_distribution distDF(0, 1.0);
  my_uniform_real_distribution distPX(-50 * sqrt(Ta), 50 * sqrt(Ta));

  for (int p = 0; p < Np; p++)
  {
    double uperp2 = u1(p)*u1(p) - u2(p)*u2(p);
    alpha = sqrt(1 + uperp2);
    while (1) {
      auxDF = distDF(ext_rng);
      auxPX = distPX(ext_rng);
      if (densityFunctionMaxwell(auxPX, alpha, Ta) > auxDF)
        break;
    }
    u0(p) = auxPX;
  }
}

void SPECIE::add_momenta(my_rng_generator& ext_rng, double uxin, double uyin, double uzin, tempDistrib distribution)
{
  if (mygrid->withParticles == NO||isFrozen)
    return;

  if (!allocated) {
#ifndef NO_ALLOCATION
    std::cout << "Warning: species " << name << " is not allocated !" << std::endl;
#endif
    return;
  }

  if (!distribution.init) {
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
  if (mygrid->withParticles == NO||isFrozen)
    return;
  if (mygrid->isStretched()) {
    SPECIE::currentStretchedDepositionStandard(current);
    return;
  }

  double dt, gamma_i;
  int p;  // particle_int, component_int
  int i, j, i1, j1, k1, i2, j2, k2;
#ifdef OLDCURRENT
  int k;
#endif
  //int indexMaxQuadraticShape[]={1,4};
  int hii[3], wii[3];           // half integer index,   whole integer index
  double hiw[3][3], wiw[3][3];  // half integer weight,  whole integer weight
  double rr, rh, rr2, rh2;          // local coordinate to integer grid point and to half integer,     local coordinate squared
  double dvol, xx[3], vv[3];           // tensor_product,       absolute particle position

  dt = mygrid->dt;
  double* myCurrent = current->getDataPointer();
  int N_grid[3];
  current->writeN_grid(N_grid);
  int edge = mygrid->getEdge();
  if (!(mygrid->withCurrent == YES && (!isTestSpecies)))
  {
    for (p = 0; p < Np; p++)
    {
      double ux, uy, uz;
#ifdef _ACC_SINGLE_POINTER
      ux = pData[pIndex(3, p, Ncomp, Np)];
      uy = pData[pIndex(4, p, Ncomp, Np)];
      uz = pData[pIndex(5, p, Ncomp, Np)];
#else
      ux = val[3][p*Ncomp];
      uy = val[4][p*Ncomp];
      uz = val[5][p*Ncomp];
#endif
      gamma_i = 1. / sqrt(1 + ux*ux + uy*uy + uz*uz);

      for (int c = 0; c < mygrid->getDimensionality(); c++)
      {
#ifdef _ACC_SINGLE_POINTER
        vv[c] = gamma_i*pData[pIndex(c + 3, p, Ncomp, Np)];
        pData[pIndex(c, p, Ncomp, Np)] += dt*vv[c];
#else
        vv[c] = gamma_i*val[c+3][p];
        val[c][p] += dt*vv[c];
#endif
      }
    }
    return;
  }
  int Nx, Ny, Nz;
  Nx = N_grid[0];
  Ny = N_grid[1];
  Nz = N_grid[2];
  switch (mygrid->getDimensionality())
  {
  case 3:
    for (p = 0; p < Np; p++) {
      double ux, uy, uz;
#ifdef _ACC_SINGLE_POINTER
      ux = pData[pIndex(3, p, Ncomp, Np)];
      uy = pData[pIndex(4, p, Ncomp, Np)];
      uz = pData[pIndex(5, p, Ncomp, Np)];
#else
      ux = val[3][p*Ncomp];
      uy = val[4][p*Ncomp];
      uz = val[5][p*Ncomp];
#endif
      gamma_i = 1. / sqrt(1 + ux*ux + uy*uy + uz*uz);

      for (int c = 0; c < 3; c++)
      {
#ifdef _ACC_SINGLE_POINTER
        vv[c] = gamma_i*pData[pIndex(c + 3, p, Ncomp, Np)];
#else
        vv[c] = gamma_i*val[c + 3][p];
#endif
        hiw[c][1] = wiw[c][1] = 1;
        hiw[c][0] = wiw[c][0] = 0;
        hiw[c][2] = wiw[c][2] = 0;
        hii[c] = wii[c] = 0;
      }
      for (int c = 0; c < 3; c++)
      {
#ifdef _ACC_SINGLE_POINTER
        xx[c] = pData[pIndex(c, p, Ncomp, Np)] + 0.5*dt*vv[c];
        pData[pIndex(c, p, Ncomp, Np)] += dt*vv[c];
#else
        xx[c] = val[c][p] + 0.5*dt*vv[c];
        val[c][p] += dt*vv[c];
#endif
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

#ifdef _ACC_SINGLE_POINTER
      double myweight = pData[pIndex(6, p, Ncomp, Np)];
#else
      double myweight = val[6][p];
#endif
#ifdef OLDCURRENT

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

            double weight = pData[pIndex(6, p, Ncomp, Np)];
            double *JX, *JY, *JZ;
            JX = &myCurrent[my_indice(edge, 1, 1, 0, i2, j1, k1, N_grid[0], N_grid[1], N_grid[2], current->Ncomp)];
            JY = &myCurrent[my_indice(edge, 1, 1, 1, i1, j2, k1, N_grid[0], N_grid[1], N_grid[2], current->Ncomp)];
            JZ = &myCurrent[my_indice(edge, 1, 1, 2, i1, j1, k2, N_grid[0], N_grid[1], N_grid[2], current->Ncomp)];

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
#else
      {
        double *JXLLL, *JXCLL, *JXRLL, *JXLCL, *JXCCL, *JXRCL, *JXLRL, *JXCRL, *JXRRL;
        double *JYLLL, *JYCLL, *JYRLL, *JYLCL, *JYCCL, *JYRCL, *JYLRL, *JYCRL, *JYRRL;
        double *JZLLL, *JZCLL, *JZRLL, *JZLCL, *JZCCL, *JZRCL, *JZLRL, *JZCRL, *JZRRL;

        double *JXLLC, *JXCLC, *JXRLC, *JXLCC, *JXCCC, *JXRCC, *JXLRC, *JXCRC, *JXRRC;
        double *JYLLC, *JYCLC, *JYRLC, *JYLCC, *JYCCC, *JYRCC, *JYLRC, *JYCRC, *JYRRC;
        double *JZLLC, *JZCLC, *JZRLC, *JZLCC, *JZCCC, *JZRCC, *JZLRC, *JZCRC, *JZRRC;

        double *JXLLR, *JXCLR, *JXRLR, *JXLCR, *JXCCR, *JXRCR, *JXLRR, *JXCRR, *JXRRR;
        double *JYLLR, *JYCLR, *JYRLR, *JYLCR, *JYCCR, *JYRCR, *JYLRR, *JYCRR, *JYRRR;
        double *JZLLR, *JZCLR, *JZRLR, *JZLCR, *JZCCR, *JZRCR, *JZLRR, *JZCRR, *JZRRR;


        int iiL, iiC, iiR, hiL, hiC, hiR;
        iiL = wii[0] - 1;
        iiC = wii[0];
        iiR = wii[0] + 1;
        hiL = hii[0] - 1;
        hiC = hii[0];
        hiR = hii[0] + 1;

        int ijL, ijC, ijR, hjL, hjC, hjR;
        ijL = wii[1] - 1;
        ijC = wii[1];
        ijR = wii[1] + 1;
        hjL = hii[1] - 1;
        hjC = hii[1];
        hjR = hii[1] + 1;

        int ikL, ikC, ikR, hkL, hkC, hkR;
        ikL = wii[2] - 1;
        ikC = wii[2];
        ikR = wii[2] + 1;
        hkL = hii[2] - 1;
        hkC = hii[2];
        hkR = hii[2] + 1;


        {
          JXLLL = &myCurrent[my_indice(edge, 1, 1, 0, hiL, ijL, ikL, Nx, Ny, Nz, current->Ncomp)];
          JXCLL = &myCurrent[my_indice(edge, 1, 1, 0, hiC, ijL, ikL, Nx, Ny, Nz, current->Ncomp)];
          JXRLL = &myCurrent[my_indice(edge, 1, 1, 0, hiR, ijL, ikL, Nx, Ny, Nz, current->Ncomp)];
          JXLCL = &myCurrent[my_indice(edge, 1, 1, 0, hiL, ijC, ikL, Nx, Ny, Nz, current->Ncomp)];
          JXCCL = &myCurrent[my_indice(edge, 1, 1, 0, hiC, ijC, ikL, Nx, Ny, Nz, current->Ncomp)];
          JXRCL = &myCurrent[my_indice(edge, 1, 1, 0, hiR, ijC, ikL, Nx, Ny, Nz, current->Ncomp)];
          JXLRL = &myCurrent[my_indice(edge, 1, 1, 0, hiL, ijR, ikL, Nx, Ny, Nz, current->Ncomp)];
          JXCRL = &myCurrent[my_indice(edge, 1, 1, 0, hiC, ijR, ikL, Nx, Ny, Nz, current->Ncomp)];
          JXRRL = &myCurrent[my_indice(edge, 1, 1, 0, hiR, ijR, ikL, Nx, Ny, Nz, current->Ncomp)];

          JYLLL = &myCurrent[my_indice(edge, 1, 1, 1, iiL, hjL, ikL, Nx, Ny, Nz, current->Ncomp)];
          JYCLL = &myCurrent[my_indice(edge, 1, 1, 1, iiC, hjL, ikL, Nx, Ny, Nz, current->Ncomp)];
          JYRLL = &myCurrent[my_indice(edge, 1, 1, 1, iiR, hjL, ikL, Nx, Ny, Nz, current->Ncomp)];
          JYLCL = &myCurrent[my_indice(edge, 1, 1, 1, iiL, hjC, ikL, Nx, Ny, Nz, current->Ncomp)];
          JYCCL = &myCurrent[my_indice(edge, 1, 1, 1, iiC, hjC, ikL, Nx, Ny, Nz, current->Ncomp)];
          JYRCL = &myCurrent[my_indice(edge, 1, 1, 1, iiR, hjC, ikL, Nx, Ny, Nz, current->Ncomp)];
          JYLRL = &myCurrent[my_indice(edge, 1, 1, 1, iiL, hjR, ikL, Nx, Ny, Nz, current->Ncomp)];
          JYCRL = &myCurrent[my_indice(edge, 1, 1, 1, iiC, hjR, ikL, Nx, Ny, Nz, current->Ncomp)];
          JYRRL = &myCurrent[my_indice(edge, 1, 1, 1, iiR, hjR, ikL, Nx, Ny, Nz, current->Ncomp)];

          JZLLL = &myCurrent[my_indice(edge, 1, 1, 2, iiL, ijL, hkL, Nx, Ny, Nz, current->Ncomp)];
          JZCLL = &myCurrent[my_indice(edge, 1, 1, 2, iiC, ijL, hkL, Nx, Ny, Nz, current->Ncomp)];
          JZRLL = &myCurrent[my_indice(edge, 1, 1, 2, iiR, ijL, hkL, Nx, Ny, Nz, current->Ncomp)];
          JZLCL = &myCurrent[my_indice(edge, 1, 1, 2, iiL, ijC, hkL, Nx, Ny, Nz, current->Ncomp)];
          JZCCL = &myCurrent[my_indice(edge, 1, 1, 2, iiC, ijC, hkL, Nx, Ny, Nz, current->Ncomp)];
          JZRCL = &myCurrent[my_indice(edge, 1, 1, 2, iiR, ijC, hkL, Nx, Ny, Nz, current->Ncomp)];
          JZLRL = &myCurrent[my_indice(edge, 1, 1, 2, iiL, ijR, hkL, Nx, Ny, Nz, current->Ncomp)];
          JZCRL = &myCurrent[my_indice(edge, 1, 1, 2, iiC, ijR, hkL, Nx, Ny, Nz, current->Ncomp)];
          JZRRL = &myCurrent[my_indice(edge, 1, 1, 2, iiR, ijR, hkL, Nx, Ny, Nz, current->Ncomp)];
        }

        {
          JXLLC = &myCurrent[my_indice(edge, 1, 1, 0, hiL, ijL, ikC, Nx, Ny, Nz, current->Ncomp)];
          JXCLC = &myCurrent[my_indice(edge, 1, 1, 0, hiC, ijL, ikC, Nx, Ny, Nz, current->Ncomp)];
          JXRLC = &myCurrent[my_indice(edge, 1, 1, 0, hiR, ijL, ikC, Nx, Ny, Nz, current->Ncomp)];
          JXLCC = &myCurrent[my_indice(edge, 1, 1, 0, hiL, ijC, ikC, Nx, Ny, Nz, current->Ncomp)];
          JXCCC = &myCurrent[my_indice(edge, 1, 1, 0, hiC, ijC, ikC, Nx, Ny, Nz, current->Ncomp)];
          JXRCC = &myCurrent[my_indice(edge, 1, 1, 0, hiR, ijC, ikC, Nx, Ny, Nz, current->Ncomp)];
          JXLRC = &myCurrent[my_indice(edge, 1, 1, 0, hiL, ijR, ikC, Nx, Ny, Nz, current->Ncomp)];
          JXCRC = &myCurrent[my_indice(edge, 1, 1, 0, hiC, ijR, ikC, Nx, Ny, Nz, current->Ncomp)];
          JXRRC = &myCurrent[my_indice(edge, 1, 1, 0, hiR, ijR, ikC, Nx, Ny, Nz, current->Ncomp)];

          JYLLC = &myCurrent[my_indice(edge, 1, 1, 1, iiL, hjL, ikC, Nx, Ny, Nz, current->Ncomp)];
          JYCLC = &myCurrent[my_indice(edge, 1, 1, 1, iiC, hjL, ikC, Nx, Ny, Nz, current->Ncomp)];
          JYRLC = &myCurrent[my_indice(edge, 1, 1, 1, iiR, hjL, ikC, Nx, Ny, Nz, current->Ncomp)];
          JYLCC = &myCurrent[my_indice(edge, 1, 1, 1, iiL, hjC, ikC, Nx, Ny, Nz, current->Ncomp)];
          JYCCC = &myCurrent[my_indice(edge, 1, 1, 1, iiC, hjC, ikC, Nx, Ny, Nz, current->Ncomp)];
          JYRCC = &myCurrent[my_indice(edge, 1, 1, 1, iiR, hjC, ikC, Nx, Ny, Nz, current->Ncomp)];
          JYLRC = &myCurrent[my_indice(edge, 1, 1, 1, iiL, hjR, ikC, Nx, Ny, Nz, current->Ncomp)];
          JYCRC = &myCurrent[my_indice(edge, 1, 1, 1, iiC, hjR, ikC, Nx, Ny, Nz, current->Ncomp)];
          JYRRC = &myCurrent[my_indice(edge, 1, 1, 1, iiR, hjR, ikC, Nx, Ny, Nz, current->Ncomp)];

          JZLLC = &myCurrent[my_indice(edge, 1, 1, 2, iiL, ijL, hkC, Nx, Ny, Nz, current->Ncomp)];
          JZCLC = &myCurrent[my_indice(edge, 1, 1, 2, iiC, ijL, hkC, Nx, Ny, Nz, current->Ncomp)];
          JZRLC = &myCurrent[my_indice(edge, 1, 1, 2, iiR, ijL, hkC, Nx, Ny, Nz, current->Ncomp)];
          JZLCC = &myCurrent[my_indice(edge, 1, 1, 2, iiL, ijC, hkC, Nx, Ny, Nz, current->Ncomp)];
          JZCCC = &myCurrent[my_indice(edge, 1, 1, 2, iiC, ijC, hkC, Nx, Ny, Nz, current->Ncomp)];
          JZRCC = &myCurrent[my_indice(edge, 1, 1, 2, iiR, ijC, hkC, Nx, Ny, Nz, current->Ncomp)];
          JZLRC = &myCurrent[my_indice(edge, 1, 1, 2, iiL, ijR, hkC, Nx, Ny, Nz, current->Ncomp)];
          JZCRC = &myCurrent[my_indice(edge, 1, 1, 2, iiC, ijR, hkC, Nx, Ny, Nz, current->Ncomp)];
          JZRRC = &myCurrent[my_indice(edge, 1, 1, 2, iiR, ijR, hkC, Nx, Ny, Nz, current->Ncomp)];
        }

        {
          JXLLR = &myCurrent[my_indice(edge, 1, 1, 0, hiL, ijL, ikR, Nx, Ny, Nz, current->Ncomp)];
          JXCLR = &myCurrent[my_indice(edge, 1, 1, 0, hiC, ijL, ikR, Nx, Ny, Nz, current->Ncomp)];
          JXRLR = &myCurrent[my_indice(edge, 1, 1, 0, hiR, ijL, ikR, Nx, Ny, Nz, current->Ncomp)];
          JXLCR = &myCurrent[my_indice(edge, 1, 1, 0, hiL, ijC, ikR, Nx, Ny, Nz, current->Ncomp)];
          JXCCR = &myCurrent[my_indice(edge, 1, 1, 0, hiC, ijC, ikR, Nx, Ny, Nz, current->Ncomp)];
          JXRCR = &myCurrent[my_indice(edge, 1, 1, 0, hiR, ijC, ikR, Nx, Ny, Nz, current->Ncomp)];
          JXLRR = &myCurrent[my_indice(edge, 1, 1, 0, hiL, ijR, ikR, Nx, Ny, Nz, current->Ncomp)];
          JXCRR = &myCurrent[my_indice(edge, 1, 1, 0, hiC, ijR, ikR, Nx, Ny, Nz, current->Ncomp)];
          JXRRR = &myCurrent[my_indice(edge, 1, 1, 0, hiR, ijR, ikR, Nx, Ny, Nz, current->Ncomp)];

          JYLLR = &myCurrent[my_indice(edge, 1, 1, 1, iiL, hjL, ikR, Nx, Ny, Nz, current->Ncomp)];
          JYCLR = &myCurrent[my_indice(edge, 1, 1, 1, iiC, hjL, ikR, Nx, Ny, Nz, current->Ncomp)];
          JYRLR = &myCurrent[my_indice(edge, 1, 1, 1, iiR, hjL, ikR, Nx, Ny, Nz, current->Ncomp)];
          JYLCR = &myCurrent[my_indice(edge, 1, 1, 1, iiL, hjC, ikR, Nx, Ny, Nz, current->Ncomp)];
          JYCCR = &myCurrent[my_indice(edge, 1, 1, 1, iiC, hjC, ikR, Nx, Ny, Nz, current->Ncomp)];
          JYRCR = &myCurrent[my_indice(edge, 1, 1, 1, iiR, hjC, ikR, Nx, Ny, Nz, current->Ncomp)];
          JYLRR = &myCurrent[my_indice(edge, 1, 1, 1, iiL, hjR, ikR, Nx, Ny, Nz, current->Ncomp)];
          JYCRR = &myCurrent[my_indice(edge, 1, 1, 1, iiC, hjR, ikR, Nx, Ny, Nz, current->Ncomp)];
          JYRRR = &myCurrent[my_indice(edge, 1, 1, 1, iiR, hjR, ikR, Nx, Ny, Nz, current->Ncomp)];

          JZLLR = &myCurrent[my_indice(edge, 1, 1, 2, iiL, ijL, hkR, Nx, Ny, Nz, current->Ncomp)];
          JZCLR = &myCurrent[my_indice(edge, 1, 1, 2, iiC, ijL, hkR, Nx, Ny, Nz, current->Ncomp)];
          JZRLR = &myCurrent[my_indice(edge, 1, 1, 2, iiR, ijL, hkR, Nx, Ny, Nz, current->Ncomp)];
          JZLCR = &myCurrent[my_indice(edge, 1, 1, 2, iiL, ijC, hkR, Nx, Ny, Nz, current->Ncomp)];
          JZCCR = &myCurrent[my_indice(edge, 1, 1, 2, iiC, ijC, hkR, Nx, Ny, Nz, current->Ncomp)];
          JZRCR = &myCurrent[my_indice(edge, 1, 1, 2, iiR, ijC, hkR, Nx, Ny, Nz, current->Ncomp)];
          JZLRR = &myCurrent[my_indice(edge, 1, 1, 2, iiL, ijR, hkR, Nx, Ny, Nz, current->Ncomp)];
          JZCRR = &myCurrent[my_indice(edge, 1, 1, 2, iiC, ijR, hkR, Nx, Ny, Nz, current->Ncomp)];
          JZRRR = &myCurrent[my_indice(edge, 1, 1, 2, iiR, ijR, hkR, Nx, Ny, Nz, current->Ncomp)];
        }

        {

          double WSVx, WSVy, WSVz;
          WSVx = myweight*vv[0] * chargeSign;
          WSVy = myweight*vv[1] * chargeSign;
          WSVz = myweight*vv[2] * chargeSign;

          {
            *JXLLL += WSVx*hiw[0][0] * wiw[1][0] * wiw[2][0];
            *JXCLL += WSVx*hiw[0][1] * wiw[1][0] * wiw[2][0];
            *JXRLL += WSVx*hiw[0][2] * wiw[1][0] * wiw[2][0];
            *JXLCL += WSVx*hiw[0][0] * wiw[1][1] * wiw[2][0];
            *JXCCL += WSVx*hiw[0][1] * wiw[1][1] * wiw[2][0];
            *JXRCL += WSVx*hiw[0][2] * wiw[1][1] * wiw[2][0];
            *JXLRL += WSVx*hiw[0][0] * wiw[1][2] * wiw[2][0];
            *JXCRL += WSVx*hiw[0][1] * wiw[1][2] * wiw[2][0];
            *JXRRL += WSVx*hiw[0][2] * wiw[1][2] * wiw[2][0];

            *JXLLC += WSVx*hiw[0][0] * wiw[1][0] * wiw[2][1];
            *JXCLC += WSVx*hiw[0][1] * wiw[1][0] * wiw[2][1];
            *JXRLC += WSVx*hiw[0][2] * wiw[1][0] * wiw[2][1];
            *JXLCC += WSVx*hiw[0][0] * wiw[1][1] * wiw[2][1];
            *JXCCC += WSVx*hiw[0][1] * wiw[1][1] * wiw[2][1];
            *JXRCC += WSVx*hiw[0][2] * wiw[1][1] * wiw[2][1];
            *JXLRC += WSVx*hiw[0][0] * wiw[1][2] * wiw[2][1];
            *JXCRC += WSVx*hiw[0][1] * wiw[1][2] * wiw[2][1];
            *JXRRC += WSVx*hiw[0][2] * wiw[1][2] * wiw[2][1];

            *JXLLR += WSVx*hiw[0][0] * wiw[1][0] * wiw[2][2];
            *JXCLR += WSVx*hiw[0][1] * wiw[1][0] * wiw[2][2];
            *JXRLR += WSVx*hiw[0][2] * wiw[1][0] * wiw[2][2];
            *JXLCR += WSVx*hiw[0][0] * wiw[1][1] * wiw[2][2];
            *JXCCR += WSVx*hiw[0][1] * wiw[1][1] * wiw[2][2];
            *JXRCR += WSVx*hiw[0][2] * wiw[1][1] * wiw[2][2];
            *JXLRR += WSVx*hiw[0][0] * wiw[1][2] * wiw[2][2];
            *JXCRR += WSVx*hiw[0][1] * wiw[1][2] * wiw[2][2];
            *JXRRR += WSVx*hiw[0][2] * wiw[1][2] * wiw[2][2];
          }

          {
            *JYLLL += WSVy*wiw[0][0] * hiw[1][0] * wiw[2][0];
            *JYCLL += WSVy*wiw[0][1] * hiw[1][0] * wiw[2][0];
            *JYRLL += WSVy*wiw[0][2] * hiw[1][0] * wiw[2][0];
            *JYLCL += WSVy*wiw[0][0] * hiw[1][1] * wiw[2][0];
            *JYCCL += WSVy*wiw[0][1] * hiw[1][1] * wiw[2][0];
            *JYRCL += WSVy*wiw[0][2] * hiw[1][1] * wiw[2][0];
            *JYLRL += WSVy*wiw[0][0] * hiw[1][2] * wiw[2][0];
            *JYCRL += WSVy*wiw[0][1] * hiw[1][2] * wiw[2][0];
            *JYRRL += WSVy*wiw[0][2] * hiw[1][2] * wiw[2][0];

            *JYLLC += WSVy*wiw[0][0] * hiw[1][0] * wiw[2][1];
            *JYCLC += WSVy*wiw[0][1] * hiw[1][0] * wiw[2][1];
            *JYRLC += WSVy*wiw[0][2] * hiw[1][0] * wiw[2][1];
            *JYLCC += WSVy*wiw[0][0] * hiw[1][1] * wiw[2][1];
            *JYCCC += WSVy*wiw[0][1] * hiw[1][1] * wiw[2][1];
            *JYRCC += WSVy*wiw[0][2] * hiw[1][1] * wiw[2][1];
            *JYLRC += WSVy*wiw[0][0] * hiw[1][2] * wiw[2][1];
            *JYCRC += WSVy*wiw[0][1] * hiw[1][2] * wiw[2][1];
            *JYRRC += WSVy*wiw[0][2] * hiw[1][2] * wiw[2][1];

            *JYLLR += WSVy*wiw[0][0] * hiw[1][0] * wiw[2][2];
            *JYCLR += WSVy*wiw[0][1] * hiw[1][0] * wiw[2][2];
            *JYRLR += WSVy*wiw[0][2] * hiw[1][0] * wiw[2][2];
            *JYLCR += WSVy*wiw[0][0] * hiw[1][1] * wiw[2][2];
            *JYCCR += WSVy*wiw[0][1] * hiw[1][1] * wiw[2][2];
            *JYRCR += WSVy*wiw[0][2] * hiw[1][1] * wiw[2][2];
            *JYLRR += WSVy*wiw[0][0] * hiw[1][2] * wiw[2][2];
            *JYCRR += WSVy*wiw[0][1] * hiw[1][2] * wiw[2][2];
            *JYRRR += WSVy*wiw[0][2] * hiw[1][2] * wiw[2][2];
          }

          {
            *JZLLL += WSVz*wiw[0][0] * wiw[1][0] * hiw[2][0];
            *JZCLL += WSVz*wiw[0][1] * wiw[1][0] * hiw[2][0];
            *JZRLL += WSVz*wiw[0][2] * wiw[1][0] * hiw[2][0];
            *JZLCL += WSVz*wiw[0][0] * wiw[1][1] * hiw[2][0];
            *JZCCL += WSVz*wiw[0][1] * wiw[1][1] * hiw[2][0];
            *JZRCL += WSVz*wiw[0][2] * wiw[1][1] * hiw[2][0];
            *JZLRL += WSVz*wiw[0][0] * wiw[1][2] * hiw[2][0];
            *JZCRL += WSVz*wiw[0][1] * wiw[1][2] * hiw[2][0];
            *JZRRL += WSVz*wiw[0][2] * wiw[1][2] * hiw[2][0];

            *JZLLC += WSVz*wiw[0][0] * wiw[1][0] * hiw[2][1];
            *JZCLC += WSVz*wiw[0][1] * wiw[1][0] * hiw[2][1];
            *JZRLC += WSVz*wiw[0][2] * wiw[1][0] * hiw[2][1];
            *JZLCC += WSVz*wiw[0][0] * wiw[1][1] * hiw[2][1];
            *JZCCC += WSVz*wiw[0][1] * wiw[1][1] * hiw[2][1];
            *JZRCC += WSVz*wiw[0][2] * wiw[1][1] * hiw[2][1];
            *JZLRC += WSVz*wiw[0][0] * wiw[1][2] * hiw[2][1];
            *JZCRC += WSVz*wiw[0][1] * wiw[1][2] * hiw[2][1];
            *JZRRC += WSVz*wiw[0][2] * wiw[1][2] * hiw[2][1];

            *JZLLR += WSVz*wiw[0][0] * wiw[1][0] * hiw[2][2];
            *JZCLR += WSVz*wiw[0][1] * wiw[1][0] * hiw[2][2];
            *JZRLR += WSVz*wiw[0][2] * wiw[1][0] * hiw[2][2];
            *JZLCR += WSVz*wiw[0][0] * wiw[1][1] * hiw[2][2];
            *JZCCR += WSVz*wiw[0][1] * wiw[1][1] * hiw[2][2];
            *JZRCR += WSVz*wiw[0][2] * wiw[1][1] * hiw[2][2];
            *JZLRR += WSVz*wiw[0][0] * wiw[1][2] * hiw[2][2];
            *JZCRR += WSVz*wiw[0][1] * wiw[1][2] * hiw[2][2];
            *JZRRR += WSVz*wiw[0][2] * wiw[1][2] * hiw[2][2];
          }


        }
      }
    }
#endif
    break;

  case 2:
    for (p = 0; p < Np; p++)
    {
      double ux, uy, uz;
#ifdef _ACC_SINGLE_POINTER
      ux = pData[pIndex(3, p, Ncomp, Np)];
      uy = pData[pIndex(4, p, Ncomp, Np)];
      uz = pData[pIndex(5, p, Ncomp, Np)];
#else
      ux = val[3][p*Ncomp];
      uy = val[4][p*Ncomp];
      uz = val[5][p*Ncomp];
#endif


      gamma_i = 1. / sqrt(1 + ux*ux + uy*uy + uz*uz);
      for (int c = 0; c < 3; c++)
      {
#ifdef _ACC_SINGLE_POINTER
        vv[c] = gamma_i*pData[pIndex(c + 3, p, Ncomp, Np)];
#else
        vv[c] = gamma_i*val[c + 3][p];
#endif
        hiw[c][1] = wiw[c][1] = 1;
        hiw[c][0] = wiw[c][0] = 0;
        hiw[c][2] = wiw[c][2] = 0;
        hii[c] = wii[c] = 0;
      }
      for (int c = 0; c < 2; c++)
      {
#ifdef _ACC_SINGLE_POINTER
        xx[c] = pData[pIndex(c, p, Ncomp, Np)] + 0.5*dt*vv[c];
        pData[pIndex(c, p, Ncomp, Np)] += dt*vv[c];
#else
        xx[c] = val[c][p] + 0.5*dt*vv[c];
        val[c][p] += dt*vv[c];
#endif
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
#ifdef _ACC_SINGLE_POINTER
          double weight = pData[pIndex(6, p, Ncomp, Np)];
#else
          double weight = val[6][p];
#endif
          double *JX, *JY, *JZ;
          JX = &myCurrent[my_indice(edge, 1, 0, 0, i2, j1, k1, N_grid[0], N_grid[1], N_grid[2], current->Ncomp)];
          JY = &myCurrent[my_indice(edge, 1, 0, 1, i1, j2, k1, N_grid[0], N_grid[1], N_grid[2], current->Ncomp)];
          JZ = &myCurrent[my_indice(edge, 1, 0, 2, i1, j1, k2, N_grid[0], N_grid[1], N_grid[2], current->Ncomp)];

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
#ifdef _ACC_SINGLE_POINTER
      ux = pData[pIndex(3, p, Ncomp, Np)];
      uy = pData[pIndex(4, p, Ncomp, Np)];
      uz = pData[pIndex(5, p, Ncomp, Np)];
#else
      ux = val[3][p*Ncomp];
      uy = val[4][p*Ncomp];
      uz = val[5][p*Ncomp];
#endif

      gamma_i = 1. / sqrt(1 + ux*ux + uy*uy + uz*uz);
      for (int c = 0; c < 3; c++)
      {
#ifdef _ACC_SINGLE_POINTER
        vv[c] = gamma_i*pData[pIndex(c + 3, p, Ncomp, Np)];
#else
        vv[c] = gamma_i*val[c + 3][p];
#endif
        hiw[c][1] = wiw[c][1] = 1;
        hiw[c][0] = wiw[c][0] = 0;
        hiw[c][2] = wiw[c][2] = 0;
        hii[c] = wii[c] = 0;
      }
      for (int c = 0; c < 1; c++)
      {
#ifdef _ACC_SINGLE_POINTER
        xx[c] = pData[pIndex(c, p, Ncomp, Np)] + 0.5*dt*vv[c];
        pData[pIndex(c, p, Ncomp, Np)] += dt*vv[c];
#else
        xx[c] = val[c][p] + 0.5*dt*vv[c];
        val[c][p] += dt*vv[c];
#endif

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

#ifdef _ACC_SINGLE_POINTER
      double weight = pData[pIndex(6, p, Ncomp, Np)];
#else
      double weight = val[6][p];
#endif
        double *JX, *JY, *JZ;
        JX = &myCurrent[my_indice(edge, 0, 0, 0, i2, j1, k1, N_grid[0], N_grid[1], N_grid[2], current->Ncomp)];
        JY = &myCurrent[my_indice(edge, 0, 0, 1, i1, j2, k1, N_grid[0], N_grid[1], N_grid[2], current->Ncomp)];
        JZ = &myCurrent[my_indice(edge, 0, 0, 2, i1, j1, k2, N_grid[0], N_grid[1], N_grid[2], current->Ncomp)];

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

void SPECIE::debug_warning_particle_outside_boundaries(double x, double y, double z, int nump) {
  if (x < mygrid->rminloc[0]) {
    std::cout << "Particle " << nump << " of proc. " << mygrid->myid << " is at x = " << x << " (boundary is " << mygrid->rminloc[0] << ")" << std::endl;
    flush(std::cout);
  }

  if (x > mygrid->rmaxloc[0]) {
    std::cout << "Particle " << nump << " of proc. " << mygrid->myid << " is at x = " << x << " (boundary is " << mygrid->rmaxloc[0] << ")" << std::endl;
    flush(std::cout);
  }


  if (y < mygrid->rminloc[1]) {
    std::cout << "Particle " << nump << " of proc. " << mygrid->myid << " is at y = " << y << " (boundary is " << mygrid->rminloc[1] << ")" << std::endl;
    flush(std::cout);
  }

  if (y > mygrid->rmaxloc[1]) {
    std::cout << "Particle " << nump << " of proc. " << mygrid->myid << " is at y = " << y << " (boundary is " << mygrid->rmaxloc[1] << ")" << std::endl;
    flush(std::cout);
  }


  if (z < mygrid->rminloc[2]) {
    std::cout << "Particle " << nump << " of proc. " << mygrid->myid << " is at z = " << z << " (boundary is " << mygrid->rminloc[2] << ")" << std::endl;
    flush(std::cout);
  }

  if (z > mygrid->rmaxloc[2]) {
    std::cout << "Particle " << nump << " of proc. " << mygrid->myid << " is at z = " << z << " (boundary is " << mygrid->rmaxloc[2] << ")" << std::endl;
    flush(std::cout);
  }



}

void SPECIE::currentStretchedDepositionStandard(CURRENT *current)
{

  double dt, gamma_i;
  int p, c;  // particle_int, component_int
  int i, j, k, i1, j1, k1, i2, j2, k2;
  //int indexMaxQuadraticShape[]={1,4};
  int hii[3], wii[3];           // half integer index,   whole integer index
  double hiw[3][3], wiw[3][3];  // half integer weight,  whole integer weight
  double rr, rh, rr2, rh2;      // local coordinate to integer grid point and to half integer,     local coordinate squared
  double dvol, xx[3], vv[3];    // tensor_product,       absolute particle position
  double mydr[3], myweight;
  double mycsi[3];

  dt = mygrid->dt;
  double* myCurrent = current->getDataPointer();
  int N_grid[3];
  current->writeN_grid(N_grid);
  int edge = mygrid->getEdge();
  if (mygrid->withCurrent == YES && (!isTestSpecies))
  {
    for (p = 0; p < Np; p++)
    {
      //debug_warning_particle_outside_boundaries(r0(p), r1(p), r2(p), p);
      double ux, uy, uz;
#ifdef _ACC_SINGLE_POINTER
      ux = pData[pIndex(3, p, Ncomp, Np)];
      uy = pData[pIndex(4, p, Ncomp, Np)];
      uz = pData[pIndex(5, p, Ncomp, Np)];
#else
      ux = val[3][p*Ncomp];
      uy = val[4][p*Ncomp];
      uz = val[5][p*Ncomp];
#endif

      gamma_i = 1. / sqrt(1 + ux*ux + uy*uy + uz*uz);

      for (c = 0; c < 3; c++)
      {
#ifdef _ACC_SINGLE_POINTER
        vv[c] = gamma_i*pData[pIndex(c + 3, p, Ncomp, Np)];
#else
        vv[c] = gamma_i*val[c + 3][p];
#endif
        hiw[c][1] = wiw[c][1] = 1;
        hiw[c][0] = wiw[c][0] = 0;
        hiw[c][2] = wiw[c][2] = 0;
        hii[c] = wii[c] = 0;
      }
      for (c = 0; c < mygrid->getDimensionality(); c++)
      {
#ifdef _ACC_SINGLE_POINTER
        xx[c] = pData[pIndex(c, p, Ncomp, Np)] + 0.5*dt*vv[c];
        pData[pIndex(c, p, Ncomp, Np)] += dt*vv[c];
#else
        xx[c] = val[c][p] + 0.5*dt*vv[c];
        val[c][p] += dt*vv[c];
#endif

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
#ifdef _ACC_SINGLE_POINTER
              double weight = pData[pIndex(6, p, Ncomp, Np)];
#else
              double weight = val[6][p];
#endif
              double *JX, *JY, *JZ;
              JX = &myCurrent[my_indice(edge, 1, 1, 0, i2, j1, k1, N_grid[0], N_grid[1], N_grid[2], current->Ncomp)];
              JY = &myCurrent[my_indice(edge, 1, 1, 1, i1, j2, k1, N_grid[0], N_grid[1], N_grid[2], current->Ncomp)];
              JZ = &myCurrent[my_indice(edge, 1, 1, 2, i1, j1, k2, N_grid[0], N_grid[1], N_grid[2], current->Ncomp)];

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
#ifdef _ACC_SINGLE_POINTER
            double weight = pData[pIndex(6, p, Ncomp, Np)];
#else
            double weight = val[6][p];
#endif
            double *JX, *JY, *JZ;
            JX = &myCurrent[my_indice(edge, 1, 0, 0, i2, j1, k1, N_grid[0], N_grid[1], N_grid[2], current->Ncomp)];
            JY = &myCurrent[my_indice(edge, 1, 0, 1, i1, j2, k1, N_grid[0], N_grid[1], N_grid[2], current->Ncomp)];
            JZ = &myCurrent[my_indice(edge, 1, 0, 2, i1, j1, k2, N_grid[0], N_grid[1], N_grid[2], current->Ncomp)];

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
#ifdef _ACC_SINGLE_POINTER
          double weight = pData[pIndex(6, p, Ncomp, Np)];
#else
          double weight = val[6][p];
#endif
          double *JX, *JY, *JZ;
          JX = &myCurrent[my_indice(edge, 0, 0, 0, i2, j1, k1, N_grid[0], N_grid[1], N_grid[2], current->Ncomp)];
          JY = &myCurrent[my_indice(edge, 0, 0, 1, i1, j2, k1, N_grid[0], N_grid[1], N_grid[2], current->Ncomp)];
          JZ = &myCurrent[my_indice(edge, 0, 0, 2, i1, j1, k2, N_grid[0], N_grid[1], N_grid[2], current->Ncomp)];

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
#ifdef _ACC_SINGLE_POINTER
      ux = pData[pIndex(3, p, Ncomp, Np)];
      uy = pData[pIndex(4, p, Ncomp, Np)];
      uz = pData[pIndex(5, p, Ncomp, Np)];
#else
      ux = val[3][p*Ncomp];
      uy = val[4][p*Ncomp];
      uz = val[5][p*Ncomp];
#endif

      gamma_i = 1. / sqrt(1 + ux*ux + uy*uy + uz*uz);
      for (c = 0; c < mygrid->getDimensionality(); c++)
      {
#ifdef _ACC_SINGLE_POINTER
        vv[c] = gamma_i*pData[pIndex(c + 3, p, Ncomp, Np)];
        pData[pIndex(c, p, Ncomp, Np)] += dt*vv[c];
#else
        vv[c] = gamma_i*val[c+3][p];
        val[c][p] += dt*vv[c];
#endif

      }
    }
  }
}
void SPECIE::density_deposition_standard(CURRENT *current, bool withSign)
{
  if (mygrid->withParticles == NO) {
    return;
  }


  if (mygrid->isStretched()) {
    SPECIE::densityStretchedDepositionStandard(current, withSign);
    return;
  }

  int p, c;  // particle_int, component_int
  int i, j, k, i1, j1, k1;
  //int indexMaxQuadraticShape[]={1,4};
  int wii[3];           // whole integer index
  double wiw[3][3];     // whole integer weight
  double rr, rr2;       // local coordinate to integer grid point, local coordinate squared
  double dvol, xx[3];   // tensor_product, absolute particle position
  double _chargeSign;
  if(withSign)
    _chargeSign = chargeSign;
  else
    _chargeSign = 1;

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
              current->density(i1, j1, k1) += _chargeSign*w(p)*dvol;
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
            current->density(i1, j1, k1) += _chargeSign*w(p)*dvol;
        }
      }
      break;

    case 1:
      k1 = j1 = 0;
      for (i = 0; i < 3; i++)
      {
        i1 = i + wii[0] - 1;
        dvol = wiw[0][i],
          current->density(i1, j1, k1) += _chargeSign*w(p)*dvol;
      }
      break;
    }

  }
}
void SPECIE::densityStretchedDepositionStandard(CURRENT *current, bool withSign)
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

  double _chargeSign;
  if(withSign)
    _chargeSign = chargeSign;
  else
    _chargeSign = 1;

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
              current->density(i1, j1, k1) += _chargeSign*myweight*dvol;
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
            current->density(i1, j1, k1) += _chargeSign*myweight*dvol;
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
          current->density(i1, j1, k1) += _chargeSign*myweight*dvol;
      }
      break;
    }

  }
}

void SPECIE::setParticlesPerCellXYZ(int numX, int numY, int numZ) {
  numX = (numX <= 0) ? 1 : numX;
  numY = (numY <= 0) ? 1 : numY;
  numZ = (numZ <= 0) ? 1 : numZ;
  particlePerCellXYZ[0] = numX;
  particlePerCellXYZ[1] = numY;
  particlePerCellXYZ[2] = numZ;
}

void SPECIE::setName(std::string iname) {
  name = iname;
}

double SPECIE::getKineticEnergy() {
  if (mygrid->withParticles == NO) {
    return 0.0;
  }

  if (!allocated) {
    return 0.0;
  }
  double energy = 0.0;
  for (int p = 0; p < Np; p++) {
    energy += (sqrt(1.0 + (u0(p)*u0(p) + u1(p)*u1(p) + u2(p)*u2(p))) - 1.0)*w(p);
  }
  energy *= mygrid->dr[0] * mygrid->dr[1] * mygrid->dr[2] * mygrid->ref_den*M_PI / coupling*chargeSign;
  return energy;
}

//xmin, ymin, zmin, pxmin,pymin,pzmin,emin, xmax,ymax,zmax, pxmax,pymax,pzmax,emax
void SPECIE::computeKineticEnergyWExtrems() {

  if (mygrid->withParticles == NO) {
    return;
  }
  if (!allocated) {
    return;
  }

  if (energyExtremesFlag) {
    return;
  }

  const double VERY_BIG_NUM_POS = 1.0e30;
  const double VERY_BIG_NUM_NEG = -1.0e30;

  for (int i = 0; i < 7; i++) {
    minima[i] = VERY_BIG_NUM_POS;
    maxima[i] = VERY_BIG_NUM_NEG;
  }


  double energy = 0.0;
  totalEnergy = 0.0;
  totalMomentum[0] = 0;
  totalMomentum[1] = 0;
  totalMomentum[2] = 0;
  double gamma_minus_1 = 0;

  for (int p = 0; p < Np; p++) {
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
  for (int p = 0; p < Np; p++) {
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

void SPECIE::dump(std::ofstream &ff) {
  ff.write((char*)&Np, sizeof(Np));

  for (int i = 0; i < Np; i++) {
    for (int c = 0; c < Ncomp; c++) {
      ff.write((char*)&ru(c, i), sizeof(double));
    }
  }
}
void SPECIE::dumpBigBuffer(std::ofstream &ff) {
  ff.write((char*)&Np, sizeof(Np));
#ifdef _ACC_SINGLE_POINTER
  ff.write((char*)pData, sizeof(double)*Np*Ncomp);
#else
  for (int c = 0; c < Ncomp; c++) {
  ff.write((char*)val[c], sizeof(double)*Np);
  }
#endif
}
void SPECIE::debugDump(std::ofstream &ff) {
  ff << this->name << Np << std::endl;
}

void SPECIE::reloadDump(std::ifstream &ff) {
  ff.read((char*)&Np, sizeof(Np));
  SPECIE::reallocate_species();
#ifdef _ACC_SINGLE_POINTER
  ff.read((char*)pData, sizeof(double)*Np*Ncomp);
#else
  for (int c = 0; c < Ncomp; c++) {
  ff.read((char*)val[c], sizeof(double)*Np);
  }
#endif

}

void SPECIE::reloadBigBufferDump(std::ifstream &ff) {
  ff.read((char*)&Np, sizeof(Np));
  SPECIE::reallocate_species();
  ff.read((char*)&ru(0, 0), sizeof(double)*Np*Ncomp);
}

bool SPECIE::areEnergyExtremesAvailable() {
  return energyExtremesFlag;
}


uint64_t SPECIE::printParticleNumber() {
  if (mygrid->myid != mygrid->master_proc)
    return lastParticle;
  std::cout << name << " has " << lastParticle << " particles." << std::endl;
  return lastParticle;
}


