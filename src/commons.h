/*   Copyright 2014-2017 - Andrea Sgattoni, Luca Fedeli, Stefano Sinigardi   */

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

#ifndef __COMMONS_H__
#define __COMMONS_H__

//IMPORTANT! "preproc_defs.h" to be included as VERY FIRST
#include "preproc_defs.h"
#include <mpi.h>

#define _ACC_SINGLE_POINTER
#define _REORDER_MPI_CART_PROCESSES 1

#define FREQUENCY_STDOUT_STATUS 5
#if defined (USE_BOOST) && defined(NO_CXX11)
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/exponential_distribution.hpp>
#include <boost/random/normal_distribution.hpp>

typedef boost::random::mt19937 my_rng_generator;
typedef boost::random::minstd_rand  aux_rnd_generator;
typedef boost::random::uniform_real_distribution<double> my_uniform_real_distribution;
typedef boost::random::uniform_int_distribution<int> my_uniform_int_distribution;
typedef boost::random::exponential_distribution<double>  my_exponential_distribution;
typedef boost::random::normal_distribution<double>  my_normal_distribution;
#else
#include <random>
typedef std::mt19937 my_rng_generator;
typedef std::minstd_rand  aux_rnd_generator;
typedef std::uniform_real_distribution<double> my_uniform_real_distribution;
typedef std::uniform_int_distribution<int> my_uniform_int_distribution;
typedef std::exponential_distribution<double>  my_exponential_distribution;
typedef std::normal_distribution<double>  my_normal_distribution;
#endif

#include <string>
#if defined(_MSC_VER)
#include <cstdint>
#else
#include <stdint.h>
#endif


//*****VERSION*****
#define CURRENT_VERSION "version 1.5.3"

/***************  GENERAL CONSTANTS ***************/
#define NOT_DEFINED false // -2
#define NO false //-1
#define YES true

/*********** DERIVATIVES CONSTANTS ***********/
//boundaryConditionType:
/* #define PERIODIC_BOUNDARY_CONDITION 0 */
/* #define OPEN_BOUNDARY_CONDITION 1 */
/* #define REFLECTIVE_BOUNDARY_CONDITION 2 */

//*****PHYSICAL CONSTANTS*****
#define SPEED_OF_LIGHT   2.99792458e+10 // cm/s
#define CHARGE_ELECTRON 4.803242384e-10  // statCoulomb
#define MASS_ELECTRON  9.109534e-28    // g
#define MASS_NUCLEON   1.6726485e-24   // g
#define BOLTZMANN_CONSTANT 1.380662e-16     //  erg/K
#define CLASSICAL_ELECTRON_RADIUS 2.8179403267e-13//re in cm

#define FS2MICR 0.299792458  // 1fs=0.29979.. micron
#define MICR2FS 3.335640952  // 1micron = 3.3356.. fs


//*****CONVERSION CGS-GAUSS -> S.I.*****
#define THREE 2.99792458
#define THIRD 0.333564095

// from simulation units to c.g.s.  (nelle sim. i fields non sono giusti [k, (c)dt])
#define CONVERT_E2cgs (1e+4)    // da unita' di simulaz -> c.g.s.
#define CONVERT_B2cgs (1e+4)    // da unita' di simulaz -> c.g.s.

// from simulation units to S.I.
#define CONVERT_E2SI (THREE*1e+8)    
#define CONVERT_B2SI (1.)             
#define CONVERT_CHARGE2SI (THIRD*1e-9)  
#define CONVERT_ENERGY2SI (1e-7)

// from c.g.s. to simulation units
#define CONVERT_cgs2E (1e-4)  // da campi c.g.s. -> unita' di simulaz 
#define CONVERT_cgs2B (1e-4)  // da campi c.g.s. -> unita' di simulaz 

// from S.I. to simulation units
#define CONVERT_SI2E (THIRD*1e-8)
#define CONVERT_SI2B (1)      
#define CONVERT_SI2CHARGE (THREE*1e+9)

#define Dfloat float
#define MPI_Dfloat MPI_FLOAT
#define ONE_THIRD 0.333333333333333333

// ************** FIELD TYPE *************** (forse si può togliere)
//enum fieldType {DUMMY, ELECTRO_MAGNETIC_FIELD,DENSITY_CURRENT_FIELD };

// ************** PARTICLES TYPE ***************
enum particlesType { ELECTRON, POSITRON, ION };



//*****USEFUL FUNCTIONS*****
//inline int my_indice(int edge, int YGrid_factor, int ZGrid_factor, int c, int i, int j, int k, int Nx, int Ny, int Nz, int Nc){
//  return (Nx*Ny*Nz*c + (i + edge) + YGrid_factor*Nx*(j + edge) + ZGrid_factor*Nx*Ny*(k + edge));
//}
inline uint64_t my_indice(int edge, int YGrid_factor, int ZGrid_factor, int c, int i, int j, int k, int Nx, int Ny, int Nz, int Nc) {
  return (c + Nc*(i + edge) + YGrid_factor*Nc*Nx*(j + edge) + ZGrid_factor*Nc*Nx*Ny*(k + edge));
}
inline uint64_t pIndex(int c, int p, int Ncomp, int Npart) {
  return c + p*Ncomp;
}

#if defined (_GCC)
#define MAX(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

#define MIN(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })
#else
#define MAX(a,b) ((a > b) ? (a) : (b))
#define MIN(a,b) ((a < b) ? (a) : (b))
#endif

template <class T>const T& TMIN(const T& a, const T& b) {
  return (a < b) ? a : b;
}

template <class T>const T& TMAX(const T& a, const T& b) {
  return (a > b) ? a : b;
}

inline int factorial(int n)
{
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

inline void exitWithError(int error) {
  MPI_Finalize();
  exit(error);
}


//*****LOGOS******

const std::string common_logo =
"\n"
"         d8b                                    888                \n"
"         Y8P                                    888                \n"
"                                                888                \n"
"88888b.  888  .d8888b .d8888b  8888b.  88888b.  888888 .d88b.      \n"
"888 \"88b 888 d88P\"   d88P\"        \"88b 888 \"88b 888   d8P  Y8b\n"
"888  888 888 888     888      .d888888 888  888 888   88888888     \n"
"888 d88P 888 Y88b.   Y88b.    888  888 888  888 Y88b. Y8b.         \n"
"88888P\"  888  \"Y8888P \"Y8888P \"Y888888 888  888  \"Y888 \"Y8888\n"
"888                                                                \n"
"888                                                                \n"
"888                                                                \n";

//**********************************
//http://patorjk.com/software/taag
// character "bolger"
//************** PLASMA TYPES *******

typedef union{
    int iDat;
    double dDat;
} intdouble;

#endif

