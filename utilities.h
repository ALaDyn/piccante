#ifndef _UTILITIES_H_
#define _UTILITIES_H_

#include <iomanip>
#if defined(_MSC_VER)
#include <cstdint>
#else
#include <stdint.h>
#endif
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

#if defined(USE_LOCALBOOST)
#include "boost/filesystem.hpp"
#endif

#include "grid.h"
#include "structures.h"
#include "em_field.h"
#include "particle_species.h"

void moveWindow(GRID* _mygrid, EM_FIELD* _myfield, std::vector<SPECIE*> _myspecies);

void restartFromDump(int *dumpID, GRID* _mygrid, EM_FIELD* _myfield, std::vector<SPECIE*> _myspecies);
void dumpFilesForRestart(int *dumpID, GRID* _mygrid, EM_FIELD* _myfield, std::vector<SPECIE*> _myspecies);


#endif // UTILITIES_PLUS_H

