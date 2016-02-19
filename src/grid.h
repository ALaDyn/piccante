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

#ifndef __GRID_H__
#define __GRID_H__

#define _USE_MATH_DEFINES
#define _CRT_SECURE_NO_WARNINGS

#include <mpi.h>
#include "commons.h"
#include "structures.h"
#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <iomanip>
#include <stdint.h>
#include <ctime>
#if defined(_MSC_VER)
#include <cstdlib>
#else
#include <stdlib.h>
#endif

#if defined(USE_BOOST)
#include <boost/filesystem.hpp>
#endif

enum axisBoundaryConditions {
  _notAssigned,
  _PBC,
  _Open,
  _PML
};

enum boundaryConditions {
  xPBC = 1 << 0,
  yPBC = 1 << 1,
  zPBC = 1 << 2,
  xOpen = 1 << 3,
  yOpen = 1 << 4,
  zOpen = 1 << 5,
  xPML = 1 << 6,
  yPML = 1 << 7,
  zPML = 1 << 8
};


class GRID {
public:
  int NGridNodes[3], Nloc[3], istep;   //number of total grid points , number of local gridpoints
  int uniquePointsloc[3];
  int uniquePoints[3];
  double rmin[3], rmax[3], dr[3], dri[3], *cir[3], *chr[3];        // [in micron]
  double rminloc[3], rmaxloc[3], *cirloc[3], *chrloc[3];
  double time, dt;                   // [in micron]  

  bool shouldIMove;
  int imove_mw;
  double fmove_mw;
  MPI_Comm cart_comm;
  int  myid, rnproc[3], nproc;   //INPUT myid (COMM_WORLD), domain decomposition in the three dimensions, total number of porcessor
  int  rmyid[3], master_proc;  // INPUT my 3D int coordinate, master process

  double ref_den;   // reference density
  double den_factor;

  DUMP_CONTROL dumpControl;
  bool withParticles, withCurrent;
  int *rproc_imin[3], *rproc_imax[3]; // rproc_imax[ c ][ rid[c] ]
  int *rproc_NuniquePointsloc[3];   // rproc_NuniquePointsloc[ c ][ rid[c] ]
  int *proc_totUniquePoints;

  double *iStretchingDerivativeCorrection[3];
  double *hStretchingDerivativeCorrection[3];

  GRID(int dimensions);
  ~GRID();

  void setXrange(double min, double max);
  void setYrange(double min, double max);
  void setZrange(double min, double max);
  void setNCells(int xcells, int ycells, int zcells);
  void setNProcsAlongY(int nproc);
  void setNProcsAlongZ(int nproc);
  int getDimensionality();
  int getEdge();
  int getNexchange();
  int alloc_number(int *N_grid, int *N_loc);
  void setCourantFactor(double courant_factor);
  void setSimulationTime(double tot_time);
  void setMovingWindow(double start, double beta, int frequency_mw);
  void setStartMovingWindow(double start);
  void setBetaMovingWindow(double beta);
  void setFrequencyMovingWindow(int frequency_mw);
  void setMasterProc(int idMasterProc);
  int getTotalNumberOfTimesteps();
  void moveWindow();
  void printTStepEvery(int every);
  void initRNG(std::mt19937 &rng, uint32_t auxiliary_seed);
  void visualDiag();
  void mpi_grid_initialize(int *narg, char **args);
  void finalize();
  void computeTotUniquePoints();
  void setXandNxLeftStretchedGrid(double min, int N);
  void setYandNyLeftStretchedGrid(double min, int N);
  void setZandNzLeftStretchedGrid(double min, int N);
  void setXandNxRightStretchedGrid(double max, int N);
  void setYandNyRightStretchedGrid(double max, int N);
  void setZandNzRightStretchedGrid(double max, int N);
  void setNxLeftStretchedGrid(int N);
  void setNyLeftStretchedGrid(int N);
  void setNzLeftStretchedGrid(int N);
  void setNxUniformGrid(int N);
  void setNyUniformGrid(int N);
  void setNzUniformGrid(int N);
  //void setNxRightStretcheGrid(int N);
  //void setNyRightStretcheGrid(int N);
  //void setNzRightStretcheGrid(int N);
  double stretchingFunction(double csi, int c);
  double inverseStretchingFunction(double x, int c);
  double derivativeStretchingFunction(double csi, int c);
  double stretchGrid(double xi_x, int c);
  double unStretchGrid(double x, int c);
  void computeDerivativeCorrection();
  void enableStretchedGrid();
  void setBoundaries(int flags);
  bool isStretched();
  bool isStretchedAlong(int c);
  bool isLeftStretchedAlong(int c);
  bool isRightStretchedAlong(int c);
  axisBoundaryConditions getXBoundaryConditions();
  axisBoundaryConditions getYBoundaryConditions();
  axisBoundaryConditions getZBoundaryConditions();
  double getTotalTime();
  void dump(std::ofstream &ff);
  void debugDump(std::ofstream &ff);
  void reloadDump(std::ifstream &ff);
  double csimin[3], csimax[3];
  double csiminloc[3], csimaxloc[3];
  double getMarkMW();
  void setDumpPath(std::string _dumpDir);
  std::string composeDumpFileName(int dumpID);
  void enableRadiationFriction();
  void disableRadiationFriction();
  bool isRadiationFrictionEnabled();
  void setLambda0(double lambda0);
  double getLambda0();

private:
  int dimensions;
  static const int edge = 2;
  static const int Nexchange = 1;
  double courantFactor;
  bool isDimensionalitySet;

  int totalNumberOfTimesteps;

  double beta_mw, t_start_moving_mw, mark_mw;

  int frequency_mw_shifts;

  int  cyclic[3];   //cyclic conditions for MPI_CART

  double *rproc_rmin[3], *rproc_rmax[3]; //rminloc for each processor, rmaxloc for each processor in the 3D integer space
  int *rproc_Nloc[3]; //   rproc_Nloc[ c ][ rid[c] ]

  double totalTime;
  bool withMovingWindow;
  // =========== STRETCHED GRID ========
  bool flagLeftStretchedAlong[3], flagRightStretchedAlong[3];
  bool flagStretchedAlong[3], flagStretched;
  int NUniformGrid[3], NLeftStretcheGrid[3], NRightStretcheGrid[3];
  double leftAlphaStretch[3], rightAlphaStretch[3], rminUniformGrid[3], rmaxUniformGrid[3];
  double *rproc_csimin[3], *rproc_csimax[3]; //csiminloc for each processor, csimaxloc for each processor in the 3D integer space
  axisBoundaryConditions xBoundaryConditions, yBoundaryConditions, zBoundaryConditions;

  bool checkAssignBoundary(axisBoundaryConditions cond, axisBoundaryConditions* axisCond);

  time_t unix_time_start;
  std::string dumpPath;

  bool radiationFrictionFlag;
  double lambda0; //lunghezza fisica con cui sono normalizzate tutte le lunghezze

  bool isValidDimension();
  void checkDimensionality();
  void setGridDeltar();
  void setGridDeltarNormal();
  void setGridDeltarStretched();
  void printLogo();
  void printProcInformations();
  void checkProcNumber();
  void printGridProcessorInformation();
  void printParametersInformations();
  void allocateRProcQuantities();
  void computeRProcNloc();
  void computeRProcNuniquePointsLoc();
  void setRminRmax();
  void setRminRmaxStretched();
  void setCsiminCsimax();
  void setIminImax();
  void checkProcNumberInitialization();
  void setLocalExtrems();
  void initializeStretchParameters();
  void checkStretchedGridInitialization();
  void checkStretchedGridNpointsAlong(int c);
  void checkStretchedGridExtensionAlong(int c);
  void computeNStretchedPointsAlong(int c);
  void computeAlphaStretchAlong(int c);
  void computeAlphaStretchLeft(int c);
  void computeAlphaStretchRight(int c);
  void computeExtremsStretchedGrid(int c);
  void emergencyStop();

};

#endif

