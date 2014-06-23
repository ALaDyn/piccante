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

#ifndef __COMMONS_h__
#define __COMMONS_h__

#define _USE_MATH_DEFINES
//#define USE_HDF5
#define _REORDER_MPI_CART_PROCESSES 1

#include <string>
#if defined(_MSC_VER)
#include <cstdint>
#else
#include <stdint.h>
#endif

//*****VERSION*****
#define CURRENT_VERSION " version 1.2.4"

const std::string common_versionLine = "================  "CURRENT_VERSION"  ================";

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
#define UN_TERZO 0.333333333333333333

// ************** FIELD TYPE *************** (forse si può togliere)
//enum fieldType {DUMMY, ELECTRO_MAGNETIC_FIELD,DENSITY_CURRENT_FIELD };

// ************** PARTICLES TYPE ***************
enum particlesType{ ELECTRON, POSITRON, ION };



//*****USEFUL FUNCTIONS*****
template <class T>const T& MIN(const T& a, const T& b){
	return (a<b) ? a : b;
}

template <class T>const T& MAX(const T& a, const T& b){
	return (a>b) ? a : b;
}


//*****ACCESS PROTOTYPE*****
class ACCESSO{
public:
	ACCESSO();
    static const int edge = 2;
	int indice(int c, int i, int j, int k, int Nx, int Ny, int Nz, int Nc);
	void alloc_number(int *N_grid, int *N_loc);
	int dimensions;
    static const int Nexchange = 1;
};

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



#endif

