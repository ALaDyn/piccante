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

#ifndef JSONPARSER_H
#define JSONPARSER_H

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
#include <map>
#include "commons.h"
#if defined(USE_BOOST)
#include <boost/filesystem.hpp>
#endif
#include "grid.h"
#include "structures.h"
#include "em_field.h"
#include "particle_species.h"
#include "json/json.h"

namespace jsonParser{
    extern bool isThisJsonMaster;

    struct laserPulseBoolFlags{

        bool type, pol, waist, a, lambda, duration, initialPosition, focusPosition, rotation, riseTime;
        laserPulseBoolFlags(){
            type = pol = waist = a = lambda = duration = initialPosition = focusPosition = rotation = riseTime =false;
        }
    };

    bool checkVersion(Json::Value &document, int &version);
    void parseJsonInputFile(Json::Value &root, std::string nomeFile);
    int getDimensionalityFromJson(Json::Value &document, int defaultDimensionality);
    void setXrangeFromJson(Json::Value &parent, GRID *grid);
    void setYrangeFromJson(Json::Value &parent,GRID *grid);
    void setZrangeFromJson(Json::Value &parent,GRID *grid);

    bool setIntFromJson(int *number, Json::Value &parent, const char* name);
    bool setDoubleFromJson(double *number, Json::Value &parent,const char* name);
    bool setBoolFromJson(bool *number, Json::Value &parent, const char* name);
    bool setStringFromJson(std::string * number, Json::Value  &parent,const char* name);

    bool setValueFromJson(Json::Value &child, Json::Value &parent, const char* name);
    void setNCellsFromJson(Json::Value &parent, GRID *grid);
    void setNprocsFromJson(Json::Value &document, GRID *grid);
    void setSimulationTimeFromJson(Json::Value  &document,GRID *grid);
    void setBoundaryConditionsFromJson(Json::Value &parent,GRID *grid);
    void setDumpControlFromJson(Json::Value  &parent, DUMP_CONTROL *myDumpControl);
    void setStretchedGridFromJson(Json::Value &document, GRID *grid);
    void setMovingWindowFromJson(Json::Value &document, GRID *grid);
    void setLaserPulsesFromJson(Json::Value &document, EM_FIELD *emfield);
    void setPlasmasFromJson(Json::Value &document, std::map<std::string, PLASMA*> &map);
    bool checkSpecEssentials(Json::Value &child, std::map<std::string, PLASMA*> plasmas);
    bool addDistribution(std::string distName, Json::Value &child, gsl_rng* ext_rng, SPECIE* spec);
    void setSpeciesFromJson(Json::Value &document, std::vector<SPECIE*> &species, std::map<std::string, PLASMA*> plasmas, GRID* mygrid, gsl_rng* ext_rng);


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
    bool checkLaserBoolFlags(laserPulseBoolFlags, laserPulse*);
    int findPlasmaFunction(std::string);


}
#endif // JSONPARSER_H
