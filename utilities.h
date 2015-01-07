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

#ifndef __UTILITIES_H__
#define __UTILITIES_H__

#if defined(_MSC_VER)
#define _CRT_SECURE_NO_WARNINGS
#include <cstdint>
#else
#include <stdint.h>
#endif
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
#include "grid.h"
#include "structures.h"
#include "em_field.h"
#include "particle_species.h"
//#include "rapidjson/document.h"     // rapidjson's DOM-style API
//#include "rapidjson/filestream.h"   // wrapper of C stream for prettywriter as output
#include "json/json.h"

void moveWindow(GRID* _mygrid, EM_FIELD* _myfield, std::vector<SPECIE*> _myspecies);

void restartFromDump(int *dumpID, GRID* _mygrid, EM_FIELD* _myfield, std::vector<SPECIE*> _myspecies);
void dumpFilesForRestart(int *dumpID, GRID* _mygrid, EM_FIELD* _myfield, std::vector<SPECIE*> _myspecies);
void dumpDebugFilesForRestart(int *dumpID, GRID* _mygrid, EM_FIELD* _myfield, std::vector<SPECIE*> _myspecies);

void parseJsonInputFile(Json::Value &root, std::string nomeFile);
int getDimensionalityFromJson(Json::Value &document, int defaultDimensionality);
void setXrangeFromJson(Json::Value &parent, GRID *grid);
void setYrangeFromJson(Json::Value &parent,GRID *grid);
void setZrangeFromJson(Json::Value &parent,GRID *grid);

bool setIntFromJson(int *number, Json::Value &parent, const char* name);
bool setDoubleFromJson(double *number, Json::Value &parent,const char* name);
bool setBoolFromJson(bool *number, Json::Value &parent, const char* name);
bool setValueFromJson(Json::Value &child, Json::Value &parent, const char* name);
void setNCellsFromJson(Json::Value &parent, GRID *grid);
void setNprocsFromJson(Json::Value &document, GRID *grid);
void setSimulationTimeFromJson(Json::Value  &document,GRID *grid);
void setDumpControlFromJson(Json::Value  &parent, DUMP_CONTROL *myDumpControl);
void setStretchedGridFromJson(Json::Value &document, GRID *grid);
void setMovingWindowFromJson(Json::Value &document, GRID *grid);

bool setIntFromJson2(int *number, Json::Value &parent, const char* name);
bool setDoubleFromJson2(double *number, Json::Value &parent,const char* name);

bool setBoolFromJson2(bool *number, Json::Value &parent, const char* name);
int getDimensionalityFromJson2(Json::Value &parent, int defaultDimensionality);

#endif // UTILITIES_H

