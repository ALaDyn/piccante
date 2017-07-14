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

#ifndef __UTILITIES_H__
#define __UTILITIES_H__

//IMPORTANT! "preproc_defs.h" to be included as VERY FIRST
#include "preproc_defs.h"

#include "grid.h"
#include "structures.h"
#include "em_field.h"
#include "particle_species.h"
#include "json/json.h"

#if defined(_MSC_VER)
#include <cstdint>
#else
#include <stdint.h>
#endif
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <vector>
#include <list>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include "commons.h"
#if defined(USE_BOOST)
#include <boost/filesystem.hpp>
#endif

class UTILITIES{

public:

    static void moveWindow(GRID* _mygrid, EM_FIELD* _myfield, std::vector<SPECIE*> _myspecies);

    static void considerRestartFromDump(GRID* grid, EM_FIELD* myfield, std::vector<SPECIE*> species);
    static void restartFromDump(GRID* _mygrid, EM_FIELD* _myfield, std::vector<SPECIE*> _myspecies);
    static void restartFromDump(int *dumpID, GRID* _mygrid, EM_FIELD* _myfield, std::vector<SPECIE*> _myspecies);
    static void considerDumpForRestart(GRID* grid, EM_FIELD* myfield, std::vector<SPECIE*> species);
    static void dumpFilesForRestart(GRID* _mygrid, EM_FIELD* _myfield, std::vector<SPECIE*> _myspecies);
    static void dumpFilesForRestart(int *dumpID, GRID* _mygrid, EM_FIELD* _myfield, std::vector<SPECIE*> _myspecies);
    static void dumpDebugFilesForRestart(int *dumpID, GRID* _mygrid, EM_FIELD* _myfield, std::vector<SPECIE*> _myspecies);

    static bool doesFileExist(const char *fileName);

    static void exitWithError(int error);
    static void splitCommGetRankNproc(MPI_Comm parentComm, MPI_Comm *childComm, int color, int *rank, int *NProcs);

    static void readAndAllocateSpheres(SPHERES &spheres, std::string filename, GRID &grid);
    static void fromCoordsToSpheresCoords(double &x, double min, double max);
    static bool isSphereInside(SPHERES& spheres, int index, GRID &grid);
    static void swapSpheres(SPHERES &spheres, int i, int j);
    static void selectSpheres(SPHERES &spheres, GRID &grid);

    static void readAndAllocateFFTplasma(FFTPLASMA &myfft, std::string filename, GRID &grid);

    static void allocateAccessibleKModes(GRIDmodes &gridModes, GRID &grid);
    static void writeGridModes(GRIDmodes &gridModes, GRID &grid);
    static void setKModesToBeInitialised(std::vector<KMODE> &myKModes, GRIDmodes &gridModes, double amplitude, double *centralK, double *sigmaK);
    static void newSetKModesToBeInitialised(LANGMUIRset &langmuirSet);
    static void exchangeKModesToBeInitialised(std::vector<KMODE> &myKModes, GRID &grid);
    static void writeKModesToBeInitialised(std::vector<KMODE> &myKModes, GRID &grid);
    static void writeLangmuirSetDataFile(LANGMUIRset &langmuirSet, GRID &grid);
    static void setLangmuirWaveSet(LANGMUIRset &langmuirSet, GRID &grid);

    static void OldMoveParticles(GRID* grid, SPECIE* specie, double amplitude,double lambda);
    static void moveParticles(GRID* grid, SPECIE* specie, std::vector<KMODE> myKModes, bool enableStandingWaves);
    static void setExternaField(EM_FIELD &exfield, GRID &mygrid, double time, LANGMUIRset &langmuirSet);

    static void printTotalNumberOfParticles(std::vector<SPECIE*> species, GRID &mygrid);
    static void launchPoissonSolver(EM_FIELD &myfield, std::vector<SPECIE*> species, GRID &grid, CURRENT &current);

};

#endif // UTILITIES_H

