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

#include "commons.h"

#ifndef __ACCESS_H__
#define __ACCESS_H__


#if !defined DIMENSIONALITY
#error  It is required to define DIMENSIONALITY before including access
#endif


#if DIMENSIONALITY==3
ACCESSO::ACCESSO(){
    //edge = 2;
	dimensions = 3;
    //Nexchange = 1;
}
int ACCESSO::indice(int c, int i, int j, int k, int Nx, int Ny, int Nz, int Nc){
	return (c + Nc*(i + edge) + Nc*Nx*(j + edge) + Nc*Nx*Ny*(k + edge));
}
void ACCESSO::alloc_number(int *N_grid, int *N_loc){
	N_grid[0] = N_loc[0] + 2 * edge;
	N_grid[1] = N_loc[1] + 2 * edge;
	N_grid[2] = N_loc[2] + 2 * edge;
}

#elif DIMENSIONALITY==2
ACCESSO::ACCESSO(){
    //edge=2;
	dimensions=2;
    //Nexchange=1;
}
int ACCESSO::indice(int c, int i, int j, int k, int Nx, int Ny, int Nz, int Nc){
	return (c+Nc*(i+edge)+Nc*Nx*(j+edge));
}
void ACCESSO::alloc_number(int *N_grid, int *N_loc){
	N_grid[0]=N_loc[0]+2*edge;
	N_grid[1]=N_loc[1]+2*edge;
	N_grid[2]=1;
}
#elif DIMENSIONALITY==1
ACCESSO::ACCESSO(){
    //edge=2;
	dimensions=1;
    //Nexchange=1;
}
int ACCESSO::indice(int c, int i, int j, int k, int Nx, int Ny, int Nz, int Nc){
	return (c+Nc*(i+edge));}
void ACCESSO::alloc_number(int *N_grid, int *N_loc){
	N_grid[0]=N_loc[0]+2*edge;
	N_grid[1]=1;
	N_grid[2]=1;
}
#else
#error DIMENSIONALITY must be set to 1,2 or 3 !
#endif

#endif

