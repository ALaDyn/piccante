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

#ifndef JSONPARSER_H
#define JSONPARSER_H

#if defined(_MSC_VER)
#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES
#endif

#include "grid.h"
#include "structures.h"
#include "em_field.h"
#include "particle_species.h"
#include "output_manager.h"
#include "utilities.h"
#include "json/json.h"
#if defined(_MSC_VER)
#include <cstdint>
#else
#include <stdint.h>
#endif
#include <cmath>
#include <iomanip>
#include <vector>
#include <list>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include "commons.h"
#if defined(USE_BOOST)
#include <boost/filesystem.hpp>
#endif

namespace jsonParser {
  extern bool isThisJsonMaster;
  extern int inputVersion;
  struct laserPulseBoolFlags {

    bool type, pol, waist, a, lambda, duration, initialPosition, focusPosition, rotation, riseTime;
    laserPulseBoolFlags() {
      type = pol = waist = a = lambda = duration = initialPosition = focusPosition = rotation = riseTime = false;
    }
  };

  struct outRequest {
    bool isFrom;
    bool isTo;
    bool isEvery;
    bool isAt;
    bool isDomainName;
    bool isSpecName;
    double from;
    double to;
    double every;
    double at;
    std::string domainName;
    std::string specName;
  };

  bool checkVersion(Json::Value &document, int &version);
  void lookForInputFile(int narg, char **args, std::string *inputFileName);

  std::string parseJsonInputFile(Json::Value &root, int narg, char **args);
  int getDimensionality(Json::Value &document, int defaultDimensionality);
  void setXrange(Json::Value &parent, GRID *grid);
  void setYrange(Json::Value &parent, GRID *grid);
  void setZrange(Json::Value &parent, GRID *grid);

  bool setInt(int *number, Json::Value &parent, const char* name);
  bool setDouble(double *number, Json::Value &parent, const char* name);
  bool setBool(bool *number, Json::Value &parent, const char* name);
  bool setString(std::string * number, Json::Value  &parent, const char* name);

  bool setValue(Json::Value &child, Json::Value &parent, const char* name);
  void setRadiationFriction(Json::Value &document, GRID *grid);
  bool getLambda0(Json::Value &document, double& lambda0);
  void setNCells(Json::Value &parent, GRID *grid);
  void setNprocs(Json::Value &document, GRID *grid);
  void setSimulationTime(Json::Value  &document, GRID *grid);
  void setMasterProc(Json::Value  &document, GRID *grid);
  void setCourantFactor(Json::Value  &document, GRID *grid);
  void setBoundaryConditions(Json::Value &parent, GRID *grid);
  void setDumpControl(Json::Value  &parent, GRID *mygrid);
  void setStretchedGrid(Json::Value &document, GRID *grid);
  void setMovingWindow(Json::Value &document, GRID *grid);
  void setLaserPulses(Json::Value &document, EM_FIELD *emfield);

  void setPoissonSolver(Json::Value &document, GRID *grid);

  void setPlasmas(Json::Value &document, std::map<std::string, PLASMA*> &map);
  bool checkSpecEssentials(Json::Value &child, std::map<std::string, PLASMA*> plasmas);
  bool addDistribution(std::string distName, Json::Value &child, my_rng_generator &ext_rng, SPECIE* spec);
  void setSpecies(Json::Value &document, std::vector<SPECIE*> &species, std::map<std::string, PLASMA*> plasmas, GRID* mygrid, my_rng_generator &ext_rng);

  bool setLaserType(laserPulse*, Json::Value&);
  bool setLaserPolarization(laserPulse*, Json::Value&);
  bool setLaserDurationFWHM(laserPulse*, Json::Value&);
  bool setLaserInitialPosition(laserPulse*, Json::Value&);
  bool setLaserAmplitude(laserPulse*, Json::Value&);
  bool setLaserWaist(laserPulse*, Json::Value&);
  bool setLaserFocusPosition(laserPulse*, Json::Value&);
  bool setLaserLambda(laserPulse*, Json::Value&);
  bool setLaserRotation(laserPulse*, Json::Value&);
  bool setLaserRiseTime(laserPulse*, Json::Value&);
  bool setLaserLG_l(laserPulse*, Json::Value&);
  bool setLaserLG_m(laserPulse*, Json::Value&);
  bool checkLaserBoolFlags(laserPulseBoolFlags, laserPulse*);
  int findPlasmaFunction(std::string);

  void setDomains(Json::Value &document, std::map<std::string, outDomain*>  &domains);

  bool isOutTypeOk(std::string type);
  bool isSpecNameOk(std::string specName, std::vector<SPECIE*> &species);

  bool addOutReq(std::string type, OUTPUT_MANAGER &manager, std::map<std::string, outDomain*>  &domains, jsonParser::outRequest outRequestVals);

  void setOutputRequests(Json::Value &document, OUTPUT_MANAGER &manager, std::map<std::string, outDomain*>  &domains, std::vector<SPECIE*> &species);
  void setOutputDirPath(Json::Value &document, OUTPUT_MANAGER &manager);
  void setOutputParameters(Json::Value &document, OUTPUT_MANAGER &manager);

}
#endif // JSONPARSER_H
