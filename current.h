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

#ifndef __CURRENT_H__
#define __CURRENT_H__

#define _USE_MATH_DEFINES

#include <mpi.h>
#include <malloc.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>

#include "commons.h"
#include "grid.h"
#include "structures.h"

using namespace std;

class CURRENT{
public:
	int Ncomp; // N grid point including ghost cells,N_grid[0]*N_grid[1]*N_grid[2], comp number
	//double max_value[6],min_value[6];  //12 utility values

	CURRENT();
	void allocate(GRID *grid); //field allocation 
	void reallocate();	//REALLOCATION only if load balancing is introduced
	void setAllValuesToZero();
	CURRENT operator = (CURRENT &destro);

	integer_or_halfinteger getJCoords(int c);
	integer_or_halfinteger getDensityCoords();

	void pbc();
	void output_density(ofstream &ff);
	void output_2D_density(ofstream &ff);
	void output_2D_current(ofstream &ff);
	void output_2D_nth_comp(ofstream &ff, int c);
	void output_2D_nth_comp_human(ofstream &ff, int c);
	double integrate_field_nth_comp(int c);
	//void difference(CURRENT *right); 
	void set_const_field_nthcomp(double value, int c);
	void eraseDensity();

	//PUBLIC INLINE FUNCTIONS
	inline double & Jx(int i, int j, int k){ return val[my_indice(acc.edge, YGrid_factor, ZGrid_factor, 0, i, j*YGrid_factor, k*ZGrid_factor, N_grid[0], N_grid[1], N_grid[2], Ncomp)]; }
	inline double & Jy(int i, int j, int k){ return val[my_indice(acc.edge, YGrid_factor, ZGrid_factor, 1, i, j*YGrid_factor, k*ZGrid_factor, N_grid[0], N_grid[1], N_grid[2], Ncomp)]; }
	inline double & Jz(int i, int j, int k){ return val[my_indice(acc.edge, YGrid_factor, ZGrid_factor, 2, i, j*YGrid_factor, k*ZGrid_factor, N_grid[0], N_grid[1], N_grid[2], Ncomp)]; }
	inline double & density(int i, int j, int k){ return val[my_indice(acc.edge, YGrid_factor, ZGrid_factor, 3, i, j*YGrid_factor, k*ZGrid_factor, N_grid[0], N_grid[1], N_grid[2], Ncomp)]; }
	inline double & JJ(int c, int i, int j, int k){ return val[my_indice(acc.edge, YGrid_factor, ZGrid_factor, c, i, j*YGrid_factor, k*ZGrid_factor, N_grid[0], N_grid[1], N_grid[2], Ncomp)]; }



private:
	ACCESSO acc;    // object distinguishing 1-2-3 D
	double *val; //   THE BIG poiniter
	GRID *mygrid;         // pointer to the GIRD object 
	int allocated;  //flag 1-0 allocaded-not alloc
	int N_grid[3], Ntot;
	int ZGrid_factor, YGrid_factor;



	//PRIVATE INLINE FUNCTIONS
	inline int my_indice(int edge, int YGrid_factor, int ZGrid_factor, int c, int i, int j, int k, int Nx, int Ny, int Nz, int Nc){
		return (c + Nc*(i + edge) + YGrid_factor*Nc*Nx*(j + edge) + ZGrid_factor*Nc*Nx*Ny*(k + edge));
	}


};

#endif
