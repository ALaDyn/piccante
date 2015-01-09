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


#include "jsonparser.h"

using namespace jsonParser;

bool jsonParser::isThisJsonMaster=true;

bool jsonParser::checkVersion(Json::Value &document, int &version){
    setIntFromJson(&version, document, "VERSION");
    if(version==1)
        return true;
    else
        return false;
}

void jsonParser::parseJsonInputFile(Json::Value &root, std::string nomeFile){

    std::ifstream inputFile(nomeFile.c_str());
    std::stringstream buffer;
    buffer << inputFile.rdbuf();
    inputFile.close();
    Json::Reader reader;
    bool parsedSuccess = reader.parse(buffer.str().c_str(), root, false);

    int version = -1;
    if((!checkVersion(root, version))&&isThisJsonMaster){
        std::cout << "WARNING: version " << version <<
                     " may not be fully supported!" << std::endl;
    }

    if((!parsedSuccess)&&isThisJsonMaster){
        std::cout<<"Failed to parse JSON"<<std::endl
                <<reader.getFormatedErrorMessages()
               <<std::endl;
        exit(1);
    }
}


bool jsonParser::setIntFromJson(int *number, Json::Value &parent, const char* name){
    bool outFlag;
    if(outFlag = ( !parent[name].isNull() ) ){
        if( parent[name].isInt())
            *number = parent[name].asInt();
    }
    return outFlag;
}


bool jsonParser::setDoubleFromJson(double *number, Json::Value  &parent,const char* name){
    bool outFlag;
    if(outFlag = ( !parent[name].isNull() ) ){
        if( parent[name].isDouble())
            *number = parent[name].asDouble();
    }
    return outFlag;
}
bool jsonParser::setBoolFromJson(bool * number, Json::Value  &parent,const char* name){
    bool outFlag;
    if(outFlag = ( !parent[name].isNull() ) ){
        if( parent[name].isBool())
            *number = parent[name].asBool();
    }
    return outFlag;
}

bool jsonParser::setStringFromJson(std::string * number, Json::Value  &parent,const char* name){
    bool outFlag;
    if(outFlag = ( !parent[name].isNull() ) ){
        if( parent[name].isString())
            *number = parent[name].asString();
    }
    return outFlag;
}


bool jsonParser::setValueFromJson(Json::Value &child, Json::Value &parent, const char* name){
    bool outFlag;
    if(outFlag = ( !parent[name].isNull() ) ){
        child=parent[name];
    }
    return outFlag;
}


int jsonParser::getDimensionalityFromJson(Json::Value &document, int defaultDimensionality){
    int dim = defaultDimensionality;
    const char* name="dimensions";
    if((!setIntFromJson(&dim, document, name))&&isThisJsonMaster)
        std::cout << "dimensions not set in JSON input file!\n";
    return dim;
}
void jsonParser::setXrangeFromJson(Json::Value &parent,GRID *grid){
    double min=-1.0, max=1.0;
    const char* name="xRange";
    if(!parent[name].isNull()){
        if(parent[name].isArray()){
            min=parent[name][0].asDouble();
            max=parent[name][1].asDouble();
        }
    }
    grid->setXrange(min, max);
}
void jsonParser::setYrangeFromJson(Json::Value &parent,GRID *grid){
    double min=-1.0, max=1.0;
    const char* name="yRange";
    if(!parent[name].isNull()){
        if(parent[name].isArray()){
            min=parent[name][0].asDouble();
            max=parent[name][1].asDouble();
        }
    }
    grid->setYrange(min, max);
}
void jsonParser::setZrangeFromJson(Json::Value &parent,GRID *grid){
    double min=-1.0, max=1.0;
    const char* name="zRange";
    if(!parent[name].isNull()){
        if(parent[name].isArray()){
            min=parent[name][0].asDouble();
            max=parent[name][1].asDouble();
        }
    }
    grid->setZrange(min, max);
}

void jsonParser::setNCellsFromJson(Json::Value &parent,GRID *grid){
    int Nx, Ny, Nz;
    Nx=Ny=Nz=1;
    const char* name="NCells";
    if(!parent[name].isNull()){
        if(parent[name].isArray()){
            Nx=parent[name][0].asInt();
            Ny=parent[name][1].asInt();
            Nz=parent[name][2].asInt();
        }
    }
    grid->setNCells(Nx, Ny, Nz);

}
void jsonParser::setNprocsFromJson(Json::Value &document,GRID *grid){
    int nProcY=1, nProcZ=1;
    setIntFromJson(&nProcY, document,"nProcY");
    setIntFromJson(&nProcZ, document,"nProcZ");

    grid->setNProcsAlongY(nProcY);
    grid->setNProcsAlongZ(nProcZ);
}

void jsonParser::setSimulationTimeFromJson(Json::Value &document,GRID *grid){
    double simulationTime;
    setDoubleFromJson(&simulationTime,document,"simulationTime");
    grid->setSimulationTime(simulationTime);
}

void jsonParser::setBoundaryConditionsFromJson(Json::Value &parent,GRID *grid){
    std::string  name1="restart";
    std::string  xCondition, yCondition, zCondition;
    int xFlag, yFlag, zFlag;
    xFlag = xPBC;
    yFlag = yPBC;
    zFlag = zPBC;

    Json::Value boundaries;
    int boundariesFlag;
    if(setValueFromJson(boundaries,parent,name1.c_str())){
        if(boundaries.size()==3){
            xCondition = boundaries[0].asString();
            yCondition = boundaries[1].asString();
            zCondition = boundaries[2].asString();
            if(! xCondition.compare("periodic"))
                xFlag = xPBC;
            else if(! xCondition.compare("open"))
                xFlag = xOpen;
            else if(! xCondition.compare("pml"))
                xFlag = xPML;

            if(! yCondition.compare("periodic"))
                yFlag = yPBC;
            else if(! yCondition.compare("open"))
                yFlag = yOpen;
            else if(! yCondition.compare("pml"))
                yFlag = yPML;

            if(! zCondition.compare("periodic"))
                zFlag = zPBC;
            else if(! zCondition.compare("open"))
                zFlag = zOpen;
            else if(! zCondition.compare("pml"))
                zFlag = zPML;

        }
    }
    boundariesFlag = xFlag | yFlag | zFlag;
    grid->setBoundaries(boundariesFlag);
}

void jsonParser::setDumpControlFromJson(Json::Value &parent, DUMP_CONTROL *myDumpControl){
    myDumpControl->doRestart = false;
    myDumpControl->doDump = false;
    std::string  name1="restart";
    std::string  name2;
    Json::Value restartObject;
    if(setValueFromJson(restartObject,parent,name1.c_str())){

        name2 = "doRestart";
        setBoolFromJson(&myDumpControl->doRestart,restartObject,name2.c_str());

        name2 = "dumpEvery";
        setDoubleFromJson(&myDumpControl->dumpEvery,restartObject,name2.c_str());

        name2 = "doDump";
        setBoolFromJson(&myDumpControl->doDump,restartObject,name2.c_str());

        name2 = "restartFromDump";
        setIntFromJson(&myDumpControl->restartFromDump,restartObject,name2.c_str());

    }
}

void jsonParser::setStretchedGridFromJson(Json::Value &document,GRID *grid){
    std::string  name1="StretchedGrid";
    Json::Value  stretching;
    if(setValueFromJson(stretching, document, name1.c_str() ) ) {
        grid->enableStretchedGrid();
        std::string  name2;
        Json::Value stretching1D;

        name2="x";
        if(setValueFromJson(stretching1D, stretching,name2.c_str() ) ){
            std::string  name3="left";
            Json::Value stretchingLeft;
            if(setValueFromJson(stretchingLeft, stretching1D,name3.c_str() ) ){
                double limit;
                int NCells;
                setIntFromJson(&NCells,stretchingLeft,"NCells");
                setDoubleFromJson(&limit, stretchingLeft, "limit");
                grid->setXandNxLeftStretchedGrid(limit,NCells);
            }
            name3="right";
            Json::Value stretchingRight;
            if(setValueFromJson(stretchingRight, stretching1D,name3.c_str() ) ){
                double limit;
                int NCells;
                setIntFromJson(&NCells,stretchingRight,"NCells");
                setDoubleFromJson(&limit, stretchingRight, "limit");
                grid->setXandNxRightStretchedGrid(limit,NCells);
            }
        }

        name2="y";
        if(setValueFromJson(stretching1D, stretching,name2.c_str() ) ){
            std::string  name3="left";
            Json::Value stretchingLeft;
            if(setValueFromJson(stretchingLeft, stretching1D,name3.c_str() ) ){
                double limit;
                int NCells;
                setIntFromJson(&NCells,stretchingLeft,"NCells");
                setDoubleFromJson(&limit, stretchingLeft, "limit");
                grid->setYandNyLeftStretchedGrid(limit,NCells);
            }
            name3="right";
            Json::Value stretchingRight;
            if(setValueFromJson(stretchingRight, stretching1D,name3.c_str() ) ){
                double limit;
                int NCells;
                setIntFromJson(&NCells,stretchingRight,"NCells");
                setDoubleFromJson(&limit, stretchingRight, "limit");
                grid->setYandNyRightStretchedGrid(limit,NCells);
            }
        }

        name2="z";
        if(setValueFromJson(stretching1D, stretching,name2.c_str() ) ){
            std::string  name3="left";
            Json::Value stretchingLeft;
            if(setValueFromJson(stretchingLeft, stretching1D,name3.c_str() ) ){
                double limit;
                int NCells;
                setIntFromJson(&NCells,stretchingLeft,"NCells");
                setDoubleFromJson(&limit, stretchingLeft, "limit");
                grid->setZandNzLeftStretchedGrid(limit,NCells);
            }
            name3="right";
            Json::Value stretchingRight;
            if(setValueFromJson(stretchingRight, stretching1D,name3.c_str() ) ){
                double limit;
                int NCells;
                setIntFromJson(&NCells,stretchingRight,"NCells");
                setDoubleFromJson(&limit, stretchingRight, "limit");
                grid->setZandNzRightStretchedGrid(limit,NCells);
            }
        }
    }
}
void jsonParser::setMovingWindowFromJson(Json::Value  &document,GRID *grid){
    std::string  name1="MovingWindow";
    Json::Value movingWindow;
    if(setValueFromJson( movingWindow, document, name1.c_str() ) ) {
        std::string  name2;
        double start=0;
        name2= "start";
        setDoubleFromJson( &start, movingWindow, name2.c_str() );
        grid->setStartMovingWindow(start);

        name2= "beta";
        double beta;
        if(setDoubleFromJson( &beta, movingWindow, name2.c_str() ) ){
            grid->setBetaMovingWindow(beta);
        }
        name2= "frequency";
        int frequency;
        if(setIntFromJson( &frequency, movingWindow, name2.c_str() ) ){
            grid->setFrequencyMovingWindow(frequency);
        }
    }

}

bool jsonParser::setLaserType(laserPulse *pulse1, Json::Value  &mylaser){
    std::string name2="type";
    std::string type;
    bool flag=false;
    if(flag=setStringFromJson(&type,mylaser,name2.c_str())){
        if(type=="COS2_PLANE_WAVE"){
            pulse1->type = COS2_PLANE_WAVE;
        }
        else if (type=="GAUSSIAN"){
            pulse1->type = GAUSSIAN;
        }
        else if (type=="PLANE_WAVE"){
            pulse1->type = PLANE_WAVE;
        }
        else if (type=="COS2_PLATEAU_PLANE_WAVE"){
            pulse1->type = COS2_PLATEAU_PLANE_WAVE;
        }
        else
            flag=false;
    }
    return flag;
}

bool jsonParser::setLaserPolarization(laserPulse *pulse1, Json::Value  &mylaser){
    std::string name2="polarization";
    std::string polarization;
    bool flag=false;
    if(flag=setStringFromJson(&polarization,mylaser,name2.c_str())){
        if(polarization=="P"){
            pulse1->setPPolarization();
        }
        else if( polarization=="S"){
            pulse1->setSPolarization();
        }
        else if(polarization=="C"){
            pulse1->setCircularPolarization();
        }
        else
            flag=false;
    }
    return flag;
}

bool jsonParser::setLaserDurationFWHM(laserPulse *pulse1, Json::Value  &mylaser){
    std::string name2="durationFWHM";
    double durationFWHM;
    bool flag=false;
    if(flag=setDoubleFromJson(&durationFWHM,mylaser,name2.c_str())){
        pulse1->setDurationFWHM(durationFWHM);
    }
    return flag;
}

bool jsonParser::setLaserInitialPosition(laserPulse *pulse1, Json::Value  &mylaser){
    std::string name2="initialPosition";
    double initialPosition;
    bool flag=false;
    if(flag=setDoubleFromJson(&initialPosition,mylaser,name2.c_str())){
        pulse1->setPulseInitialPosition(initialPosition);
    }
    return flag;
}

bool jsonParser::setLaserAmplitude(laserPulse *pulse1, Json::Value  &mylaser){
    std::string name2="a";
    double amplitude;
    bool flag=false;
    if(flag=setDoubleFromJson(&amplitude,mylaser,name2.c_str())){
        pulse1->setNormalizedAmplitude(amplitude);
    }
    return flag;
}

bool jsonParser::setLaserWaist(laserPulse *pulse1, Json::Value  &mylaser){
    std::string name2="waist";
    double waist;
    bool flag=false;
    if(flag=setDoubleFromJson(&waist,mylaser,name2.c_str())){
        pulse1->setWaist(waist);
    }
    return flag;
}
bool jsonParser::setLaserFocusPosition(laserPulse *pulse1, Json::Value  &mylaser){
    std::string name2="focusPosition";
    double focusPosition;
    bool flag=false;
    if(flag=setDoubleFromJson(&focusPosition,mylaser,name2.c_str())){
        pulse1->setFocusPosition(focusPosition);
    }
    return flag;
}
bool jsonParser::setLaserLambda(laserPulse *pulse1, Json::Value  &mylaser){
    std::string name2="lambda";
    double lambda;
    bool flag=false;
    if(flag=setDoubleFromJson(&lambda,mylaser,name2.c_str())){
        pulse1->setLambda(lambda);
    }
    return flag;
}

bool jsonParser::setLaserRotation(laserPulse *pulse1, Json::Value  &mylaser){
    std::string name2="rotation";
    bool rotation;
    bool flag=false;
    if(flag=setBoolFromJson(&rotation,mylaser,name2.c_str())){

        double angle=0, center=0;
        name2="angle";
        setDoubleFromJson(&angle,mylaser,name2.c_str());

        name2="center";
        setDoubleFromJson(&center,mylaser,name2.c_str());

        pulse1->setRotationAngleAndCenter(2.0*M_PI*(angle / 360.0), center);
    }
    return flag;
}

bool jsonParser::setLaserRiseTime(laserPulse *pulse1, Json::Value  &mylaser){
    std::string name2="riseTime";
    double riseTime;
    bool flag=false;
    if(flag=setDoubleFromJson(&riseTime,mylaser,name2.c_str())){
        pulse1->setRiseTime(riseTime);
    }
    return flag;
}

bool jsonParser::checkLaserBoolFlags(laserPulseBoolFlags flags, laserPulse *pulse){
    if(!flags.type){
        return false;
    }
    switch (pulse->type){
    case GAUSSIAN:
        if(!(flags.a && flags.waist && flags.duration) ){
            return false;
        }
        break;
    case PLANE_WAVE:
        if(!(flags.a) ){
            return false;
        }
        break;
    case COS2_PLANE_WAVE:
        if(!(flags.a && flags.duration) ){
            return false;
        }
        break;
    case COS2_PLATEAU_PLANE_WAVE:
        if(!(flags.a && flags.duration) ){
            return false;
        }
        break;
    default:
        break;
    }
    return true;
}

void jsonParser::setLaserPulsesFromJson(Json::Value &document, EM_FIELD *emfield){
    std::string  name1="Laser";
    Json::Value lasers;

    if(setValueFromJson( lasers, document, name1.c_str() ) ) {
        if(lasers.isArray()){
            std::string  name2;

            for(unsigned int index=0; index<lasers.size(); index++){
                Json::Value myLaser = lasers[index];

                name2="enabled";
                bool enabled=false;
                setBoolFromJson(&enabled, myLaser, name2.c_str());
                if(enabled){
                    laserPulse *pulse1=new(laserPulse);

                    laserPulseBoolFlags flags;

                    flags.type = setLaserType(pulse1, myLaser);
                    flags.pol = setLaserPolarization(pulse1, myLaser);
                    flags.duration = setLaserDurationFWHM(pulse1, myLaser);
                    flags.initialPosition = setLaserInitialPosition(pulse1, myLaser);
                    flags.a = setLaserAmplitude(pulse1, myLaser);
                    flags.waist = setLaserWaist(pulse1, myLaser);
                    flags.focusPosition = setLaserFocusPosition(pulse1, myLaser);
                    flags.rotation = setLaserRotation(pulse1, myLaser);
                    flags.lambda = setLaserLambda(pulse1, myLaser);
                    flags.riseTime = setLaserRiseTime(pulse1, myLaser);

                    if(checkLaserBoolFlags(flags, pulse1)){
                        emfield->addPulse(pulse1);
                    }
                    else if(isThisJsonMaster){
                        std::cout << "warning: laser #" << index << " is incorrectly defined\n";
                    }

                    delete pulse1;
                }
            }
        }
    }
}

int jsonParser::findPlasmaFunction(std::string plasmaFunction){
    for (int i = 0; i < PLASMA::maxdF; i++){
        if (!plasmaFunction.compare(PLASMA::dFNames[i]))
            return i;
    }
    return -1;
}

void jsonParser::setPlasmasFromJson(Json::Value &document, std::map<std::string, PLASMA*> &map){
    std::string  name1="Plasma";
    Json::Value plasmas;

    if(setValueFromJson( plasmas, document, name1.c_str() ) && plasmas.isArray()){
        for(unsigned int index=0; index<plasmas.size(); index++){
            Json::Value myPlasma = plasmas[index];

            std::string plasmaName;
            setStringFromJson(&plasmaName,myPlasma,"name");

            if(map.find(plasmaName)!=map.end()&&isThisJsonMaster){
                std::cout << "Warning! Plasma " << plasmaName << " is defined multiple times!" << std::endl;
            }

            std::string plasmaFunction;
            bool isPlasmaFunctionSet=false;

            int plasmaFunctionIndex;
            if(setStringFromJson(&plasmaFunction, myPlasma, "densityFunction")){
                plasmaFunctionIndex = findPlasmaFunction(plasmaFunction);
                if(plasmaFunctionIndex >= 0)
                    isPlasmaFunctionSet=true;
            }
            if(isPlasmaFunctionSet){
                map[plasmaName] = new PLASMA();
                map[plasmaName]->density_function = PLASMA::dFPoint[plasmaFunctionIndex];

                double tdouble, t2double;
                Json::Value range;

                if(setValueFromJson(range, myPlasma, "XRangeBox") && (range.size() == 2))
                    map[plasmaName]->setXRangeBox(range[0].asDouble(), range[1].asDouble());

                if(setValueFromJson(range, myPlasma, "YRangeBox") && (range.size() == 2))
                    map[plasmaName]->setYRangeBox(range[0].asDouble(), range[1].asDouble());

                if(setValueFromJson(range, myPlasma, "ZRangeBox") && (range.size() == 2))
                    map[plasmaName]->setZRangeBox(range[0].asDouble(), range[1].asDouble());


                if(setDoubleFromJson(&tdouble, myPlasma, "DensityLambda")){
                    if(setDoubleFromJson(&t2double, myPlasma, "DensityCoefficient"))
                        map[plasmaName]->setDensityCoefficient(t2double,tdouble);
                }
                else{
                    if(setDoubleFromJson(&tdouble, myPlasma, "DensityCoefficient"))
                        map[plasmaName]->setDensityCoefficient(tdouble);
                }
                if(setDoubleFromJson(&tdouble, myPlasma, "RampLength"))
                    map[plasmaName]->setRampLength(tdouble);

                if(setDoubleFromJson(&tdouble, myPlasma, "RampMinDensity"))
                    map[plasmaName]->setRampMinDensity(tdouble);

                if(setDoubleFromJson(&tdouble, myPlasma, "LeftRampLength"))
                    map[plasmaName]->setLeftRampLength(tdouble);

                if(setDoubleFromJson(&tdouble, myPlasma, "RightRampLength"))
                    map[plasmaName]->setRightRampLength(tdouble);

                if(setDoubleFromJson(&tdouble, myPlasma, "ScaleLength"))
                    map[plasmaName]->setScaleLength(tdouble);

                if(setDoubleFromJson(&tdouble, myPlasma, "LeftScaleLength"))
                    map[plasmaName]->setLeftScaleLength(tdouble);

                if(setDoubleFromJson(&tdouble, myPlasma, "RightScaleLength"))
                    map[plasmaName]->setRightScaleLength(tdouble);

                if(setDoubleFromJson(&tdouble, myPlasma, "LeftRampMinDensity"))
                    map[plasmaName]->setLeftRampMinDensity(tdouble);

                if(setDoubleFromJson(&tdouble, myPlasma, "RightRampMinDensity"))
                    map[plasmaName]->setRightRampMinDensity(tdouble);

                if(PLASMA::isGrating(plasmaFunctionIndex)){
                    double g_depth = 0;
                    double g_period = 1.0;
                    double g_phase = 0.0;

                    double* additionalParams = new double[3];

                    setDoubleFromJson(&g_depth, myPlasma, "GratingDepth");
                    setDoubleFromJson(&g_period, myPlasma, "GratingPeriod");
                    setDoubleFromJson(&g_phase, myPlasma, "GratingPhase");

                    additionalParams[0] = g_depth;
                    additionalParams[1] = g_period;
                    additionalParams[2] = g_phase;

                    map[plasmaName]->setAdditionalParams(additionalParams);
                }
            }
            else{
                if(isThisJsonMaster)
                    std::cout << "Warning! Plasma " << plasmaName <<
                                 " has no valid density function. It will be ignored." << std::endl;
            }

        }
    }
}

bool jsonParser::checkSpecEssentials(Json::Value &child, std::map<std::string, PLASMA*> plasmas){
    bool isThereName=false;
    bool isThereType=false;
    bool isTherePlasma=false;
    bool isTherePPC=false;

    std::string dummy;

    isThereName=setStringFromJson(&dummy,child,"name");

    if(setStringFromJson(&dummy,child,"plasma")){
        if(plasmas.find(dummy)!=plasmas.end())
            isTherePlasma=true;
    }

    Json::Value ppc;

    if(setValueFromJson(ppc, child, "ParticlesPerCell")){
      if (ppc.isArray()&&(ppc.size()==3))
        if(ppc[0].asInt() >= 0 && ppc[1].asInt() >= 0 && ppc[2].asInt() >= 0)
            isTherePPC=true;
    }

    if(setStringFromJson(&dummy, child, "type")){
        if(dummy.compare("ELECTRON") ||
                dummy.compare("POSITRON")||
                dummy.compare("ION"))
            isThereType=true;
    }

    return isThereName&&isThereType&&isTherePlasma&&isTherePPC;
}

bool jsonParser::addDistribution(std::string distName, Json::Value &child, gsl_rng* ext_rng, SPECIE* spec){
  bool isThereMomentum = false;
  bool isDistOk = false;

  tempDistrib dist;

  Json::Value momentum;
  if(setValueFromJson(momentum, child, "distributionAddMomentum")&&momentum.isArray()&&(momentum.size()==3))
    isThereMomentum=true;

  Json::Value params;
  if(!distName.compare("Maxwell")){
    if(setValueFromJson(params, child, "distributionParams")&&params.isArray()&&(params.size()>=1)){
      dist.setMaxwell(params[0].asDouble());
      isDistOk=true;
    }
  }
  else if(!distName.compare("Juttner")){
    if(setValueFromJson(params, child, "distributionParams")&&params.isArray()&&(params.size()>=1)){
      dist.setJuttner(params[0].asDouble());
      isDistOk=true;
    }
  }
  else if(!distName.compare("Waterbag")){
    if(setValueFromJson(params, child, "distributionParams")&&params.isArray()&&(params.size()>=1)){
      dist.setWaterbag(params[0].asDouble());
      isDistOk=true;
    }
  }
  else if(!distName.compare("Waterbag3Temp")){
    if(setValueFromJson(params, child, "distributionParams")&&params.isArray()&&(params.size()>=3)){
      dist.setWaterbag3Temp(params[0].asDouble(),params[1].asDouble(),params[2].asDouble());
      isDistOk=true;
    }
  }
  else if(!distName.compare("UnifSphere")){
    if(setValueFromJson(params, child, "distributionParams")&&params.isArray()&&(params.size()>=1)){
      dist.setUnifSphere(params[0].asDouble());
      isDistOk=true;
    }
  }
  else if(!distName.compare("Supergaussian")){
    if(setValueFromJson(params, child, "distributionParams")&&params.isArray()&&(params.size()>=2)){
      dist.setSupergaussian(params[0].asDouble(),params[1].asDouble());
      isDistOk=true;
    }
  }
  else if(!distName.compare("Special")){
    if(setValueFromJson(params, child, "distributionParams")&&params.isArray()&&(params.size()>=1)){
      dist.setSpecial(params[0].asDouble());
      isDistOk=true;
    }
  }

  if(isDistOk){
    if(isThereMomentum){
      spec->add_momenta(ext_rng, momentum[0].asDouble(), momentum[1].asDouble(), momentum[2].asDouble(), dist);
    }
    else{
      spec->add_momenta(ext_rng, 0.0, 0.0, 0.0, dist);
    }
    return true;
  }

  return false;
}

void jsonParser::setSpeciesFromJson(Json::Value &document, std::vector<SPECIE*> &species,  std::map<std::string, PLASMA*> plasmas, GRID* myGrid, gsl_rng* ext_rng){
  std::string  name1="Species";
  Json::Value specList;
  if(setValueFromJson(specList, document, name1.c_str() ) && specList.isArray()){
    for(unsigned int index=0; index<specList.size(); index++){
      Json::Value mySpecies = specList[index];
      if(checkSpecEssentials(mySpecies, plasmas)){
        SPECIE *newSpec = new SPECIE(myGrid);

        std::string dummy;
        Json::Value ppc;

        setStringFromJson(&dummy, mySpecies, "plasma");

        newSpec->plasma = *plasmas[dummy];

        setStringFromJson(&dummy,mySpecies,"type");
        if(!dummy.compare("ION"))
           newSpec->type=ION;
        else if(!dummy.compare("POSITRON"))
          newSpec->type=POSITRON;
        else
          newSpec->type=ELECTRON;

        setStringFromJson(&dummy,mySpecies,"name");
        newSpec->setName(dummy);

        setValueFromJson(ppc, mySpecies, "ParticlesPerCell");
        newSpec->setParticlesPerCellXYZ(ppc[0].asInt(), ppc[1].asInt(), ppc[2].asInt());

        if(newSpec->type==ION){
          double A,Z;
          if(setDoubleFromJson(&A, mySpecies, "A"))
            newSpec->A = A;
          if(setDoubleFromJson(&Z, mySpecies, "Z"))
            newSpec->Z = Z;
        }

        int isMarker = 0;
        setIntFromJson(&isMarker, mySpecies, "isMarker");
        if(isMarker>0)
          newSpec->addMarker();

        newSpec->creation();

        if(setStringFromJson(&dummy, mySpecies, "distribution")){
          if ((!addDistribution(dummy, mySpecies, ext_rng, newSpec)) && isThisJsonMaster){
            std::cout << "Warning: temperature distribution for species "<< index
                      << " is incorrectly defined" << std::endl;
          }
        }

        species.push_back(newSpec);

      }
      else if(isThisJsonMaster){
        std::cout << "Warning: species " << index <<
                     " is not defined correctly."  << std::endl;
      }
    }
  }
}


