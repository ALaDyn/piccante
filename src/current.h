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

#ifndef __CURRENT_H__
#define __CURRENT_H__

#define _USE_MATH_DEFINES

#include <mpi.h>
#include "commons.h"
#include "grid.h"
#include "structures.h"

#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <cstring>
#include <stdint.h>



class CURRENT {
public:
  double *getDataPointer();
  void writeN_grid(int *N_grid);
  int Ncomp; // N grid point including ghost cells,N_grid[0]*N_grid[1]*N_grid[2], comp number
  //double max_value[6],min_value[6];  //12 utility values

  CURRENT();
  ~CURRENT();
  void allocate(GRID *grid); //field allocation
  void reallocate();  //REALLOCATION only if load balancing is introduced
  void setAllValuesToZero();
  CURRENT operator = (CURRENT &destro);

  integer_or_halfinteger getJCoords(int c);
  integer_or_halfinteger getDensityCoords();

  void pbc();

  void eraseDensity();

  //PUBLIC INLINE FUNCTIONS

  inline double & Jx(int i, int j, int k) { return val[my_indice(mygrid->getEdge(), YGrid_factor, ZGrid_factor, 0, i, j*YGrid_factor, k*ZGrid_factor, N_grid[0], N_grid[1], N_grid[2], Ncomp)]; }
  inline double & Jy(int i, int j, int k) { return val[my_indice(mygrid->getEdge(), YGrid_factor, ZGrid_factor, 1, i, j*YGrid_factor, k*ZGrid_factor, N_grid[0], N_grid[1], N_grid[2], Ncomp)]; }
  inline double & Jz(int i, int j, int k) { return val[my_indice(mygrid->getEdge(), YGrid_factor, ZGrid_factor, 2, i, j*YGrid_factor, k*ZGrid_factor, N_grid[0], N_grid[1], N_grid[2], Ncomp)]; }
  inline double & density(int i, int j, int k) { return val[my_indice(mygrid->getEdge(), YGrid_factor, ZGrid_factor, 3, i, j*YGrid_factor, k*ZGrid_factor, N_grid[0], N_grid[1], N_grid[2], Ncomp)]; }
  inline double & aux(int i, int j, int k) { return val[my_indice(mygrid->getEdge(), YGrid_factor, ZGrid_factor, 0, i, j*YGrid_factor, k*ZGrid_factor, N_grid[0], N_grid[1], N_grid[2], Ncomp)]; }
  inline double & JJ(int c, int i, int j, int k) { return val[my_indice(mygrid->getEdge(), YGrid_factor, ZGrid_factor, c, i, j*YGrid_factor, k*ZGrid_factor, N_grid[0], N_grid[1], N_grid[2], Ncomp)]; }



private:
  int N_grid[3];
  uint64_t Ntot;
  int ZGrid_factor, YGrid_factor;
  double *val; //   THE BIG poiniter
  GRID *mygrid;         // pointer to the GIRD object 
  int allocated;  //flag 1-0 allocaded-not alloc


};

#endif

