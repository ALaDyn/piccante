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

#ifndef __EMFIELD_H__
#define __EMFIELD_H__

#define _USE_MATH_DEFINES

#include <mpi.h>

#include "commons.h"
#include "structures.h"
#include "grid.h"
#include "current.h"

//#include <omp.h>
#ifdef _USE_FFTW_FILTER
#include <fftw3-mpi.h>
#endif

#if defined(_MSC_VER)
#include <ctime>
#include <cstdlib>
#include <io.h>
#include <process.h> /* for getpid() and the exec..() family */
#include <direct.h> /* for _getcwd() and _chdir() */

#define srandom srand
#define random rand

/* Values for the second argument to access. These may be OR'd together.  */
#define R_OK    4       /* Test for read permission.  */
#define W_OK    2       /* Test for write permission.  */
//#define   X_OK    1       /* execute permission - unsupported in windows*/
#define F_OK    0       /* Test for existence.  */

#define access _access
#define dup2 _dup2
#define execve _execve
#define ftruncate _chsize
#define unlink _unlink
#define fileno _fileno
#define getcwd _getcwd
#define chdir _chdir
#define isatty _isatty
#define lseek _lseek
/* read, write, and close are NOT being #defined here, because while there are file handle specific versions for Windows, they probably don't work for sockets. You need to look at your app and consider whether to call e.g. closesocket(). */

#define ssize_t int

#define STDIN_FILENO 0
#define STDOUT_FILENO 1
#define STDERR_FILENO 2
/* should be in some equivalent to <sys/types.h>
typedef __int8            int8_t;
typedef __int16           int16_t;
typedef __int32           int32_t;
typedef __int64           int64_t;
typedef unsigned __int8   uint8_t;
typedef unsigned __int16  uint16_t;
typedef unsigned __int32  uint32_t;
typedef unsigned __int64  uint64_t; */
#else
#include <sys/time.h>
#include <unistd.h>
#include <stdlib.h>
#endif

#include <stdint.h>



enum filterOptions {
  fltr_Ex = 1 << 0,
  fltr_Ey = 1 << 1,
  fltr_Ez = 1 << 2,
  fltr_Bx = 1 << 3,
  fltr_By = 1 << 4,
  fltr_Bz = 1 << 5
};

enum filterDir {
  dir_x = 1 << 0,
  dir_y = 1 << 1,
  dir_z = 1 << 2
};



class EM_FIELD {
public:
  double *getDataPointer();
  void writeN_grid(int *N_grid);
  double minima[6], maxima[8];  //14 utility values minima: Exmin Eymin, ..., Bzmin;     maxima: Exmax, Eymax, ..., Bzmax, Emax, Bmax
  double total_energy[7];  // Ex2, Ey2, Ez2, Bx2, By2, Bz2 E2+B2 (totalenergy)
  double total_momentum[3];  // pointing vector Sx Sy Sz


  EM_FIELD();
  ~EM_FIELD();
  EM_FIELD operator = (EM_FIELD &destro);

  void allocate(GRID *grid);
  void reallocate();
  void setAllValuesToZero();

  int getNcomp();
  integer_or_halfinteger getCompCoords(int c);
  bool amIAllocated();

  void difference(EM_FIELD *right);

  double getTotalCharge(CURRENT *current);
  double getChargeCorrection(double  totalCharge);
  void poissonSolver(CURRENT *current);
  void putNabla2ofB1inCurrentAux(CURRENT *current);
  double getErrorInPoissonEquation(CURRENT *current);

  void boundary_conditions();  // set on the ghost cells the boundary values  

  void new_halfadvance_B();
  void new_advance_E(CURRENT *current);

  static const int myWidth = 12;
  static const int myNarrowWidth = 6;
  void init_output_diag(std::ofstream &ff);
  void output_diag(int istep, std::ofstream &ff);
  void init_output_extrems(std::ofstream &ff);
  void output_extrems(int istep, std::ofstream &ff);

  void smooth_filter(int filter_points);

  void writeNewPulseInformation(laserPulse* pulse);
  void addPulse(laserPulse* pulse);
  void addFieldsFromFile(std::string name);

  void moveWindow();

  double getEBenergy(double* EEnergy, double* BEnergy);
  void computeEnergyAndExtremes();

  void openBoundariesE_1();
  void openBoundariesE_2();
  void openBoundariesB();

  void applyFilter(int flags, int dirflags);

  bool areEnergyExtremesAvailable();
  void dump(std::ofstream &ff);
  void debugDump(std::ofstream &ff);
  void reloadDump(std::ifstream &ff);

  void fftw_filter_Efield();

  //PUBLIC INLINE FUNCTIONS
  //  static inline int my_indice(int edge, int YGrid_factor, int ZGrid_factor, int c, int i, int j, int k, int Nx, int Ny, int Nz, int Nc){
  //    return (Nx*Ny*Nz*c + (i + edge) + YGrid_factor*Nx*(j + edge) + ZGrid_factor*Nx*Ny*(k + edge));
  //  }
  inline double & E0(int i, int j, int k) {
    return val[my_indice(mygrid->getEdge(), YGrid_factor, ZGrid_factor,
      0, i, j, k,
      N_grid[0], N_grid[1], N_grid[2], Ncomp)];
  }
  inline double & E1(int i, int j, int k) {
    return val[my_indice(mygrid->getEdge(), YGrid_factor, ZGrid_factor,
      1, i, j, k,
      N_grid[0], N_grid[1], N_grid[2], Ncomp)];
  }
  inline double & E2(int i, int j, int k) {
    return val[my_indice(mygrid->getEdge(), YGrid_factor, ZGrid_factor,
      2, i, j, k,
      N_grid[0], N_grid[1], N_grid[2], Ncomp)];
  }
  inline double & B0(int i, int j, int k) {
    return val[my_indice(mygrid->getEdge(), YGrid_factor, ZGrid_factor,
      3, i, j, k,
      N_grid[0], N_grid[1], N_grid[2], Ncomp)];
  }
  inline  double & B1(int i, int j, int k) {
    return val[my_indice(mygrid->getEdge(), YGrid_factor, ZGrid_factor,
      4, i, j, k,
      N_grid[0], N_grid[1], N_grid[2], Ncomp)];
  }
  inline double & B2(int i, int j, int k) {
    return val[my_indice(mygrid->getEdge(), YGrid_factor, ZGrid_factor,
      5, i, j, k,
      N_grid[0], N_grid[1], N_grid[2], Ncomp)];
  }
  inline double & VEB(int c, int i, int j, int k) {
    return val[my_indice(mygrid->getEdge(), YGrid_factor, ZGrid_factor,
      c, i, j, k,
      N_grid[0], N_grid[1], N_grid[2], Ncomp)];
  }

  //  inline double & E0(int i,int j,int k){
  //      int indice=(0+Ncomp*(i+acc.edge)+YGrid_factor*Ncomp*N_grid[0]*(j+acc.edge)+ZGrid_factor*Ncomp*N_grid[0]*N_grid[1]*(k+acc.edge));
  //      return val[indice];}
  //  inline double & E1(int i,int j,int k){
  //      int indice=(1+Ncomp*(i+acc.edge)+YGrid_factor*Ncomp*N_grid[0]*(j+acc.edge)+ZGrid_factor*Ncomp*N_grid[0]*N_grid[1]*(k+acc.edge));
  //      return val[indice];}
  //  inline double & E2(int i,int j,int k){
  //      int indice=(2+Ncomp*(i+acc.edge)+YGrid_factor*Ncomp*N_grid[0]*(j+acc.edge)+ZGrid_factor*Ncomp*N_grid[0]*N_grid[1]*(k+acc.edge));
  //      return val[indice];}
  //  inline double & B0(int i,int j,int k){
  //      int indice=(3+Ncomp*(i+acc.edge)+YGrid_factor*Ncomp*N_grid[0]*(j+acc.edge)+ZGrid_factor*Ncomp*N_grid[0]*N_grid[1]*(k+acc.edge));
  //      return val[indice];}
  //  inline  double & B1(int i,int j,int k){
  //      int indice=(4+Ncomp*(i+acc.edge)+YGrid_factor*Ncomp*N_grid[0]*(j+acc.edge)+ZGrid_factor*Ncomp*N_grid[0]*N_grid[1]*(k+acc.edge));
  //      return val[indice];}
  //  inline double & B2(int i,int j,int k){
  //      int indice=(5+Ncomp*(i+acc.edge)+YGrid_factor*Ncomp*N_grid[0]*(j+acc.edge)+ZGrid_factor*Ncomp*N_grid[0]*N_grid[1]*(k+acc.edge));
  //      return val[indice];}
  //  inline double & VEB(int c,int i,int j,int k){
  //      int indice=(c+Ncomp*(i+acc.edge)+YGrid_factor*Ncomp*N_grid[0]*(j+acc.edge)+ZGrid_factor*Ncomp*N_grid[0]*N_grid[1]*(k+acc.edge));
  //      return val[indice];}


private:
  int N_grid[3], Ncomp;
  uint64_t Ntot;
  int ZGrid_factor, YGrid_factor;
  double *val; //   THE BIG poiniter
  GRID *mygrid;         // pointer to the GIRD object 
  bool allocated;  //flag 1-0 allocaded-not alloc 

  bool EBEnergyExtremesFlag;

  void auxiliary_rotation(double xin, double yin, double &xp, double &yp, double xcenter, double theta);

  static double cos2_profile(double u);
  static double cos2_plateau_profile(double rise, double plateau, double x);
  static double cossin_profile(double u);

  int pbc_compute_alloc_size();
  void pbcExchangeAlongX(double* send_buffer, double* recv_buffer);
  void pbcExchangeAlongY(double* send_buffer, double* recv_buffer);
  void pbcExchangeAlongZ(double* send_buffer, double* recv_buffer);
  void pbc_EB();

  void gaussian_pulse(int dimensions, double xx, double yy, double zz,
    double tt, double lambda, double fwhm,
    double w0, double* field, pulsePolarization polarization);

  void laguerreGaussian_pulse(int dimensions, double xx, double yy, double zz,
    double tt, double lambda, double fwhm,
    double w0, double* field, pulsePolarization polarization,
                              int LG_l, int LG_m);

  void initialize_cos2_plane_wave_angle(double lambda0, double amplitude,
    double laser_pulse_initial_position,
    double t_FWHM, double xcenter, double angle, pulsePolarization polarization, double rise_time);

  void initialize_plane_wave_angle(double lambda0, double amplitude,
    double angle, pulsePolarization polarization);

  void initialize_gaussian_pulse_angle(double lambda0, double amplitude,
    double laser_pulse_initial_position, double t_FWHM,
    double waist, double focus_position,
    double xcenter, double angle, pulsePolarization polarization);

  void initialize_LG_pulse_angle(double lambda0, double amplitude,
    double laser_pulse_initial_position, double t_FWHM,
    double waist, double focus_position,
    double xcenter, double angle, pulsePolarization polarization, int LG_l, int LG_m);


  void filterDirSelect(int comp, int dirflags);
  void filterCompAlongX(int comp);
  void filterCompAlongY(int comp);
  void filterCompAlongZ(int comp);

  //PRIVATE INLINE FUNCTIONS


  bool checkIfFilterPossible();
};


#endif

