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
    setInt(&version, document, "VERSION");
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


bool jsonParser::setInt(int *number, Json::Value &parent, const char* name){
    bool outFlag;
    if(outFlag = ( !parent[name].isNull() ) ){
        if( parent[name].isInt())
            *number = parent[name].asInt();
    }
    return outFlag;
}


bool jsonParser::setDouble(double *number, Json::Value  &parent,const char* name){
    bool outFlag;
    if(outFlag = ( !parent[name].isNull() ) ){
        if( parent[name].isDouble())
            *number = parent[name].asDouble();
    }
    return outFlag;
}
bool jsonParser::setBool(bool * number, Json::Value  &parent,const char* name){
    bool outFlag;
    if(outFlag = ( !parent[name].isNull() ) ){
        if( parent[name].isBool())
            *number = parent[name].asBool();
    }
    return outFlag;
}

bool jsonParser::setString(std::string * number, Json::Value  &parent,const char* name){
    bool outFlag;
    if(outFlag = ( !parent[name].isNull() ) ){
        if( parent[name].isString())
            *number = parent[name].asString();
    }
    return outFlag;
}


bool jsonParser::setValue(Json::Value &child, Json::Value &parent, const char* name){
    bool outFlag;
    if(outFlag = ( !parent[name].isNull() ) ){
        child=parent[name];
    }
    return outFlag;
}


int jsonParser::getDimensionality(Json::Value &document, int defaultDimensionality){
    int dim = defaultDimensionality;
    const char* name="dimensions";
    if((!setInt(&dim, document, name))&&isThisJsonMaster)
        std::cout << "dimensions not set in JSON input file!\n";
    return dim;
}

bool jsonParser::getRadiationFriction(Json::Value &document){
  bool isFriction=false;
  setBool(&isFriction, document, "radiationFriction");
  return isFriction;
}

bool jsonParser::getLambda0(Json::Value &document, double& lambda0){
  return setDouble(&lambda0, document, "lambda0");
}

void jsonParser::setXrange(Json::Value &parent,GRID *grid){
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
void jsonParser::setYrange(Json::Value &parent,GRID *grid){
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
void jsonParser::setZrange(Json::Value &parent,GRID *grid){
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

void jsonParser::setNCells(Json::Value &parent,GRID *grid){
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
void jsonParser::setNprocs(Json::Value &document,GRID *grid){
    int nProcY=1, nProcZ=1;
    setInt(&nProcY, document,"nProcY");
    setInt(&nProcZ, document,"nProcZ");

    grid->setNProcsAlongY(nProcY);
    grid->setNProcsAlongZ(nProcZ);
}

void jsonParser::setSimulationTime(Json::Value &document,GRID *grid){
    double simulationTime;
    setDouble(&simulationTime,document,"simulationTime");
    grid->setSimulationTime(simulationTime);
}

void jsonParser::setBoundaryConditions(Json::Value &parent,GRID *grid){
    std::string  name1="boundaries";
    std::string  xCondition, yCondition, zCondition;
    int xFlag, yFlag, zFlag;
    xFlag = xPBC;
    yFlag = yPBC;
    zFlag = zPBC;

    Json::Value boundaries;
    int boundariesFlag;
    if(setValue(boundaries,parent,name1.c_str())){
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

void jsonParser::setDumpControl(Json::Value &parent, DUMP_CONTROL *myDumpControl){
    myDumpControl->doRestart = false;
    myDumpControl->doDump = false;
    std::string  name1="restart";
    std::string  name2;
    Json::Value restartObject;
    if(setValue(restartObject,parent,name1.c_str())){

        name2 = "doRestart";
        setBool(&myDumpControl->doRestart,restartObject,name2.c_str());

        name2 = "dumpEvery";
        setDouble(&myDumpControl->dumpEvery,restartObject,name2.c_str());

        name2 = "doDump";
        setBool(&myDumpControl->doDump,restartObject,name2.c_str());

        name2 = "restartFromDump";
        setInt(&myDumpControl->restartFromDump,restartObject,name2.c_str());

    }
}

void jsonParser::setStretchedGrid(Json::Value &document,GRID *grid){
    std::string  name1="StretchedGrid";
    Json::Value  stretching;
    if(setValue(stretching, document, name1.c_str() ) ) {
        grid->enableStretchedGrid();
        std::string  name2;
        Json::Value stretching1D;

        name2="x";
        if(setValue(stretching1D, stretching,name2.c_str() ) ){
            std::string  name3="left";
            Json::Value stretchingLeft;
            if(setValue(stretchingLeft, stretching1D,name3.c_str() ) ){
                double limit;
                int NCells;
                setInt(&NCells,stretchingLeft,"NCells");
                setDouble(&limit, stretchingLeft, "limit");
                grid->setXandNxLeftStretchedGrid(limit,NCells);
            }
            name3="right";
            Json::Value stretchingRight;
            if(setValue(stretchingRight, stretching1D,name3.c_str() ) ){
                double limit;
                int NCells;
                setInt(&NCells,stretchingRight,"NCells");
                setDouble(&limit, stretchingRight, "limit");
                grid->setXandNxRightStretchedGrid(limit,NCells);
            }
        }

        name2="y";
        if(setValue(stretching1D, stretching,name2.c_str() ) ){
            std::string  name3="left";
            Json::Value stretchingLeft;
            if(setValue(stretchingLeft, stretching1D,name3.c_str() ) ){
                double limit;
                int NCells;
                setInt(&NCells,stretchingLeft,"NCells");
                setDouble(&limit, stretchingLeft, "limit");
                grid->setYandNyLeftStretchedGrid(limit,NCells);
            }
            name3="right";
            Json::Value stretchingRight;
            if(setValue(stretchingRight, stretching1D,name3.c_str() ) ){
                double limit;
                int NCells;
                setInt(&NCells,stretchingRight,"NCells");
                setDouble(&limit, stretchingRight, "limit");
                grid->setYandNyRightStretchedGrid(limit,NCells);
            }
        }

        name2="z";
        if(setValue(stretching1D, stretching,name2.c_str() ) ){
            std::string  name3="left";
            Json::Value stretchingLeft;
            if(setValue(stretchingLeft, stretching1D,name3.c_str() ) ){
                double limit;
                int NCells;
                setInt(&NCells,stretchingLeft,"NCells");
                setDouble(&limit, stretchingLeft, "limit");
                grid->setZandNzLeftStretchedGrid(limit,NCells);
            }
            name3="right";
            Json::Value stretchingRight;
            if(setValue(stretchingRight, stretching1D,name3.c_str() ) ){
                double limit;
                int NCells;
                setInt(&NCells,stretchingRight,"NCells");
                setDouble(&limit, stretchingRight, "limit");
                grid->setZandNzRightStretchedGrid(limit,NCells);
            }
        }
    }
}
void jsonParser::setMovingWindow(Json::Value  &document,GRID *grid){
    std::string  name1="MovingWindow";
    Json::Value movingWindow;
    if(setValue( movingWindow, document, name1.c_str() ) ) {
        std::string  name2;
        double start=0;
        name2= "start";
        setDouble( &start, movingWindow, name2.c_str() );
        grid->setStartMovingWindow(start);

        name2= "beta";
        double beta;
        if(setDouble( &beta, movingWindow, name2.c_str() ) ){
            grid->setBetaMovingWindow(beta);
        }
        name2= "frequency";
        int frequency;
        if(setInt( &frequency, movingWindow, name2.c_str() ) ){
            grid->setFrequencyMovingWindow(frequency);
        }
    }

}

bool jsonParser::setLaserType(laserPulse *pulse1, Json::Value  &mylaser){
    std::string name2="type";
    std::string type;
    bool flag=false;
    if(flag=setString(&type,mylaser,name2.c_str())){
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
    if(flag=setString(&polarization,mylaser,name2.c_str())){
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
    if(flag=setDouble(&durationFWHM,mylaser,name2.c_str())){
        pulse1->setDurationFWHM(durationFWHM);
    }
    return flag;
}

bool jsonParser::setLaserInitialPosition(laserPulse *pulse1, Json::Value  &mylaser){
    std::string name2="initialPosition";
    double initialPosition;
    bool flag=false;
    if(flag=setDouble(&initialPosition,mylaser,name2.c_str())){
        pulse1->setPulseInitialPosition(initialPosition);
    }
    return flag;
}

bool jsonParser::setLaserAmplitude(laserPulse *pulse1, Json::Value  &mylaser){
    std::string name2="a";
    double amplitude;
    bool flag=false;
    if(flag=setDouble(&amplitude,mylaser,name2.c_str())){
        pulse1->setNormalizedAmplitude(amplitude);
    }
    return flag;
}

bool jsonParser::setLaserWaist(laserPulse *pulse1, Json::Value  &mylaser){
    std::string name2="waist";
    double waist;
    bool flag=false;
    if(flag=setDouble(&waist,mylaser,name2.c_str())){
        pulse1->setWaist(waist);
    }
    return flag;
}
bool jsonParser::setLaserFocusPosition(laserPulse *pulse1, Json::Value  &mylaser){
    std::string name2="focusPosition";
    double focusPosition;
    bool flag=false;
    if(flag=setDouble(&focusPosition,mylaser,name2.c_str())){
        pulse1->setFocusPosition(focusPosition);
    }
    return flag;
}
bool jsonParser::setLaserLambda(laserPulse *pulse1, Json::Value  &mylaser){
    std::string name2="lambda";
    double lambda;
    bool flag=false;
    if(flag=setDouble(&lambda,mylaser,name2.c_str())){
        pulse1->setLambda(lambda);
    }
    return flag;
}

bool jsonParser::setLaserRotation(laserPulse *pulse1, Json::Value  &mylaser){
    std::string name2="rotation";
    bool rotation;
    bool flag=false;
    if(flag=setBool(&rotation,mylaser,name2.c_str())){

        double angle=0, center=0;
        name2="angle";
        setDouble(&angle,mylaser,name2.c_str());

        name2="center";
        setDouble(&center,mylaser,name2.c_str());

        pulse1->setRotationAngleAndCenter(2.0*M_PI*(angle / 360.0), center);
    }
    return flag;
}

bool jsonParser::setLaserRiseTime(laserPulse *pulse1, Json::Value  &mylaser){
    std::string name2="riseTime";
    double riseTime;
    bool flag=false;
    if(flag=setDouble(&riseTime,mylaser,name2.c_str())){
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

void jsonParser::setLaserPulses(Json::Value &document, EM_FIELD *emfield){
    std::string  name1="Laser";
    Json::Value lasers;

    if(setValue( lasers, document, name1.c_str() ) ) {
        if(lasers.isArray()){
            std::string  name2;

            for(unsigned int index=0; index<lasers.size(); index++){
                Json::Value myLaser = lasers[index];

                name2="enabled";
                bool enabled=false;
                setBool(&enabled, myLaser, name2.c_str());
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

void jsonParser::setPlasmas(Json::Value &document, std::map<std::string, PLASMA*> &map){
    std::string  name1="Plasma";
    Json::Value plasmas;

    if(setValue( plasmas, document, name1.c_str() ) && plasmas.isArray()){
        for(unsigned int index=0; index<plasmas.size(); index++){
            Json::Value myPlasma = plasmas[index];

            std::string plasmaName;
            setString(&plasmaName,myPlasma,"name");

            if(map.find(plasmaName)!=map.end()&&isThisJsonMaster){
                std::cout << "Warning! Plasma " << plasmaName << " is defined multiple times!" << std::endl;
            }

            std::string plasmaFunction;
            bool isPlasmaFunctionSet=false;

            int plasmaFunctionIndex;
            if(setString(&plasmaFunction, myPlasma, "densityFunction")){
                plasmaFunctionIndex = findPlasmaFunction(plasmaFunction);
                if(plasmaFunctionIndex >= 0)
                    isPlasmaFunctionSet=true;
            }
            if(isPlasmaFunctionSet){
                map[plasmaName] = new PLASMA();
                map[plasmaName]->density_function = PLASMA::dFPoint[plasmaFunctionIndex];

                double tdouble, t2double;
                Json::Value range;

                if(setValue(range, myPlasma, "XRangeBox") && (range.size() == 2))
                    map[plasmaName]->setXRangeBox(range[0].asDouble(), range[1].asDouble());

                if(setValue(range, myPlasma, "YRangeBox") && (range.size() == 2))
                    map[plasmaName]->setYRangeBox(range[0].asDouble(), range[1].asDouble());

                if(setValue(range, myPlasma, "ZRangeBox") && (range.size() == 2))
                    map[plasmaName]->setZRangeBox(range[0].asDouble(), range[1].asDouble());


                if(setDouble(&tdouble, myPlasma, "DensityLambda")){
                    if(setDouble(&t2double, myPlasma, "DensityCoefficient"))
                        map[plasmaName]->setDensityCoefficient(t2double,tdouble);
                }
                else{
                    if(setDouble(&tdouble, myPlasma, "DensityCoefficient"))
                        map[plasmaName]->setDensityCoefficient(tdouble);
                }
                if(setDouble(&tdouble, myPlasma, "RampLength"))
                    map[plasmaName]->setRampLength(tdouble);

                if(setDouble(&tdouble, myPlasma, "RampMinDensity"))
                    map[plasmaName]->setRampMinDensity(tdouble);

                if(setDouble(&tdouble, myPlasma, "LeftRampLength"))
                    map[plasmaName]->setLeftRampLength(tdouble);

                if(setDouble(&tdouble, myPlasma, "RightRampLength"))
                    map[plasmaName]->setRightRampLength(tdouble);

                if(setDouble(&tdouble, myPlasma, "ScaleLength"))
                    map[plasmaName]->setScaleLength(tdouble);

                if(setDouble(&tdouble, myPlasma, "LeftScaleLength"))
                    map[plasmaName]->setLeftScaleLength(tdouble);

                if(setDouble(&tdouble, myPlasma, "RightScaleLength"))
                    map[plasmaName]->setRightScaleLength(tdouble);

                if(setDouble(&tdouble, myPlasma, "LeftRampMinDensity"))
                    map[plasmaName]->setLeftRampMinDensity(tdouble);

                if(setDouble(&tdouble, myPlasma, "RightRampMinDensity"))
                    map[plasmaName]->setRightRampMinDensity(tdouble);

                if(PLASMA::isGrating(plasmaFunctionIndex)){
                    double g_depth = 0;
                    double g_period = 1.0;
                    double g_phase = 0.0;

                    double* additionalParams = new double[3];

                    setDouble(&g_depth, myPlasma, "GratingDepth");
                    setDouble(&g_period, myPlasma, "GratingPeriod");
                    setDouble(&g_phase, myPlasma, "GratingPhase");

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

    isThereName=setString(&dummy,child,"name");

    if(setString(&dummy,child,"plasma")){
        if(plasmas.find(dummy)!=plasmas.end())
            isTherePlasma=true;
    }

    Json::Value ppc;

    if(setValue(ppc, child, "ParticlesPerCell")){
      if (ppc.isArray()&&(ppc.size()==3))
        if(ppc[0].asInt() >= 0 && ppc[1].asInt() >= 0 && ppc[2].asInt() >= 0)
            isTherePPC=true;
    }

    if(setString(&dummy, child, "type")){
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
  if(setValue(momentum, child, "distributionAddMomentum")&&momentum.isArray()&&(momentum.size()==3))
    isThereMomentum=true;

  Json::Value params;
  if(!distName.compare("Maxwell")){
    if(setValue(params, child, "distributionParams")&&params.isArray()&&(params.size()>=1)){
      dist.setMaxwell(params[0].asDouble());
      isDistOk=true;
    }
  }
  else if(!distName.compare("Juttner")){
    if(setValue(params, child, "distributionParams")&&params.isArray()&&(params.size()>=1)){
      dist.setJuttner(params[0].asDouble());
      isDistOk=true;
    }
  }
  else if(!distName.compare("Waterbag")){
    if(setValue(params, child, "distributionParams")&&params.isArray()&&(params.size()>=1)){
      dist.setWaterbag(params[0].asDouble());
      isDistOk=true;
    }
  }
  else if(!distName.compare("Waterbag3Temp")){
    if(setValue(params, child, "distributionParams")&&params.isArray()&&(params.size()>=3)){
      dist.setWaterbag3Temp(params[0].asDouble(),params[1].asDouble(),params[2].asDouble());
      isDistOk=true;
    }
  }
  else if(!distName.compare("UnifSphere")){
    if(setValue(params, child, "distributionParams")&&params.isArray()&&(params.size()>=1)){
      dist.setUnifSphere(params[0].asDouble());
      isDistOk=true;
    }
  }
  else if(!distName.compare("Supergaussian")){
    if(setValue(params, child, "distributionParams")&&params.isArray()&&(params.size()>=2)){
      dist.setSupergaussian(params[0].asDouble(),params[1].asDouble());
      isDistOk=true;
    }
  }
  else if(!distName.compare("Special")){
    if(setValue(params, child, "distributionParams")&&params.isArray()&&(params.size()>=1)){
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

void jsonParser::setSpecies(Json::Value &document, std::vector<SPECIE*> &species,  std::map<std::string, PLASMA*> plasmas, GRID* myGrid, gsl_rng* ext_rng){
  std::string  name1="Species";
  Json::Value specList;
  if(setValue(specList, document, name1.c_str() ) && specList.isArray()){
    for(unsigned int index=0; index<specList.size(); index++){
      Json::Value mySpecies = specList[index];
      if(checkSpecEssentials(mySpecies, plasmas)){
        SPECIE *newSpec = new SPECIE(myGrid);

        std::string dummy;
        Json::Value ppc;

        setString(&dummy, mySpecies, "plasma");

        newSpec->plasma = *plasmas[dummy];

        setString(&dummy,mySpecies,"type");
        if(!dummy.compare("ION"))
           newSpec->type=ION;
        else if(!dummy.compare("POSITRON"))
          newSpec->type=POSITRON;
        else
          newSpec->type=ELECTRON;

        setString(&dummy,mySpecies,"name");
        newSpec->setName(dummy);

        setValue(ppc, mySpecies, "ParticlesPerCell");
        newSpec->setParticlesPerCellXYZ(ppc[0].asInt(), ppc[1].asInt(), ppc[2].asInt());

        if(newSpec->type==ION){
          double A,Z;
          if(setDouble(&A, mySpecies, "A"))
            newSpec->A = A;
          if(setDouble(&Z, mySpecies, "Z"))
            newSpec->Z = Z;
        }

        int isMarker = 0;
        setInt(&isMarker, mySpecies, "isMarker");
        if(isMarker>0)
          newSpec->addMarker();

        newSpec->creation();

        if(setString(&dummy, mySpecies, "distribution")){
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

void jsonParser::setDomains(Json::Value &document, std::map<std::string, outDomain*> &domains){
  std::string  name1="Domains";
  Json::Value domainList;
  if(setValue(domainList, document, name1.c_str() ) && domainList.isArray()){
    for(unsigned int index=0; index<domainList.size(); index++){
      Json::Value domain = domainList[index];

      std::string name;
      if(setString(&name, domain, "name")&&(domains.find(name) == domains.end())){
        outDomain* newDomain = new outDomain();
        newDomain->setName(name);

        Json::Value freeDims;
        if(setValue(freeDims, domain, "freeDim")&&freeDims.isArray()&&(freeDims.size()==3)){
          newDomain->setFreeDimensions(freeDims[0].asBool(), freeDims[1].asBool(), freeDims[2].asBool());
        }

        Json::Value pointCoord;
        if(setValue(pointCoord, domain, "pointCoord")&&pointCoord.isArray()&&(pointCoord.size()==3)){
          newDomain->setPointCoordinate(pointCoord[0].asDouble(), pointCoord[1].asDouble(),  pointCoord[2].asDouble());

        }

        Json::Value xRange;
        if(setValue(xRange, domain, "xRange")&&xRange.isArray()&&(xRange.size()==2)){
          newDomain->setXRange(xRange[0].asDouble(), xRange[1].asDouble());
        }

        Json::Value yRange;
        if(setValue(yRange, domain, "yRange")&&yRange.isArray()&&(yRange.size()==2)){
          newDomain->setYRange(yRange[0].asDouble(), yRange[1].asDouble());
        }

        Json::Value zRange;
        if(setValue(zRange, domain, "zRange")&&zRange.isArray()&&(zRange.size()==2)){
          newDomain->setZRange(zRange[0].asDouble(), zRange[1].asDouble());
        }

        domains[name]=newDomain;
      }
      else if(isThisJsonMaster){
        std::cout << "Warning: output domain "<< index <<" is defined incorrectly. It will be ignored. " << std::endl;
      }
    }
  }
}

bool jsonParser::isOutTypeOk(std::string type){
  bool answ = (!type.compare("E")) || (!type.compare("B")) || (!type.compare("EB")) ||
      (!type.compare("Density")) || (!type.compare("Current")) || (!type.compare("PhSp")) ||( !type.compare("Diag"));
  return answ;
}

bool jsonParser::isSpecNameOk(std::string specName, std::vector<SPECIE*> &species){
  for (std::vector<SPECIE*>::const_iterator spec_iterator = species.begin(); spec_iterator != species.end(); spec_iterator++){
    if(!(*spec_iterator)->name.compare(specName))
      return true;
  }
  return false;
}

bool jsonParser::addOutReq(std::string type,OUTPUT_MANAGER &manager, std::map<std::string, outDomain*>  &domains, jsonParser::outRequest outRequestVals){
  if(outRequestVals.isDomainName){
    if(outRequestVals.isAt){
      if(!type.compare("E")){
         manager.addEFieldAt(domains[outRequestVals.domainName],outRequestVals.at);
      }
      else if(!type.compare("B")){
        manager.addBFieldAt(domains[outRequestVals.domainName],outRequestVals.at);
      }
      else if(!type.compare("EB")){
        manager.addEBFieldAt(domains[outRequestVals.domainName],outRequestVals.at);
      }
      else if(!type.compare("Density")){
        manager.addSpeciesDensityAt(domains[outRequestVals.domainName],outRequestVals.specName,outRequestVals.at);
      }
      else if(!type.compare("Current")){
        manager.addCurrentAt(domains[outRequestVals.domainName],outRequestVals.at);
      }
      else if(!type.compare("PhSp")){
        manager.addSpeciesPhaseSpaceAt(domains[outRequestVals.domainName],outRequestVals.specName,outRequestVals.at);
      }
      else if(!type.compare("Diag")){
        manager.addDiagAt(outRequestVals.at);
      }

    }

    if(outRequestVals.isEvery){
      double every = outRequestVals.every;
      double from;
      if(outRequestVals.isFrom)
        from=outRequestVals.from;
      else
        from = 0.0;

      if(outRequestVals.isTo){
        if(!type.compare("E")){
          manager.addEFieldFromTo(domains[outRequestVals.domainName], from, every, outRequestVals.to);
        }
        else if(!type.compare("B")){
          manager.addBFieldFromTo(domains[outRequestVals.domainName], from, every, outRequestVals.to);
        }
        else if(!type.compare("EB")){
          manager.addEBFieldFromTo(domains[outRequestVals.domainName], from, every, outRequestVals.to);
        }
        else if(!type.compare("Density")){
          manager.addSpeciesDensityFromTo(domains[outRequestVals.domainName], outRequestVals.specName, from, every, outRequestVals.to);
        }
        else if(!type.compare("Current")){
          manager.addCurrentFromTo(domains[outRequestVals.domainName], from, every, outRequestVals.to);
        }
        else if(!type.compare("PhSp")){
          manager.addSpeciesPhaseSpaceFromTo(domains[outRequestVals.domainName], outRequestVals.specName, from, every, outRequestVals.to);
        }
        else if(!type.compare("Diag")){
          manager.addDiagFromTo(from, every, outRequestVals.to);
        }

      }
      else{
        if(!type.compare("E")){
          manager.addEFieldFrom(domains[outRequestVals.domainName], from, every);
        }
        else if(!type.compare("B")){
          manager.addBFieldFrom(domains[outRequestVals.domainName], from, every);
        }
        else if(!type.compare("EB")){
          manager.addEBFieldFrom(domains[outRequestVals.domainName], from, every);
        }
        else if(!type.compare("Density")){
          manager.addSpeciesDensityFrom(domains[outRequestVals.domainName], outRequestVals.specName, from, every);
        }
        else if(!type.compare("Current")){
          manager.addCurrentFrom(domains[outRequestVals.domainName], from, every);
        }
        else if(!type.compare("PhSp")){
          manager.addSpeciesPhaseSpaceFrom(domains[outRequestVals.domainName], outRequestVals.specName, from, every);
        }
        else if(!type.compare("Diag")){
          manager.addDiagFrom(from, every);
        }

      }
    }
  }
  else{
    if(outRequestVals.isAt){
      if(!type.compare("E")){
         manager.addEFieldAt(outRequestVals.at);
      }
      else if(!type.compare("B")){
        manager.addBFieldAt(outRequestVals.at);
      }
      else if(!type.compare("EB")){
        manager.addEBFieldAt(outRequestVals.at);
      }
      else if(!type.compare("Density")){
        manager.addSpeciesDensityAt(outRequestVals.specName,outRequestVals.at);
      }
      else if(!type.compare("Current")){
        manager.addCurrentAt(outRequestVals.at);
      }
      else if(!type.compare("PhSp")){
        manager.addSpeciesPhaseSpaceAt(outRequestVals.specName,outRequestVals.at);
      }
      else if(!type.compare("Diag")){
        manager.addDiagAt(outRequestVals.at);
      }

    }
    if(outRequestVals.isEvery){
      double every = outRequestVals.every;
      double from;
      if(outRequestVals.isFrom)
        from=outRequestVals.from;
      else
        from = 0.0;

      if(outRequestVals.isTo){
        if(!type.compare("E")){
          manager.addEFieldFromTo(from, every, outRequestVals.to);
        }
        else if(!type.compare("B")){
          manager.addBFieldFromTo(from, every, outRequestVals.to);
        }
        else if(!type.compare("EB")){
          manager.addEBFieldFromTo(from, every, outRequestVals.to);
        }
        else if(!type.compare("Density")){
          manager.addSpeciesDensityFromTo(outRequestVals.specName, from, every, outRequestVals.to);
        }
        else if(!type.compare("Current")){
          manager.addCurrentFromTo(from, every, outRequestVals.to);
        }
        else if(!type.compare("PhSp")){
          manager.addSpeciesPhaseSpaceFromTo(outRequestVals.specName, from, every, outRequestVals.to);
        }
        else if(!type.compare("Diag")){
          manager.addDiagFromTo(from, every, outRequestVals.to);
        }

      }
      else{
        if(!type.compare("E")){
          manager.addEFieldFrom(from, every);
        }
        else if(!type.compare("B")){
          manager.addBFieldFrom(from, every);
        }
        else if(!type.compare("EB")){
          manager.addEBFieldFrom(from, every);
        }
        else if(!type.compare("Density")){
          manager.addSpeciesDensityFrom(outRequestVals.specName, from, every);
        }
        else if(!type.compare("Current")){
          manager.addCurrentFrom(from, every);
        }
        else if(!type.compare("PhSp")){
          manager.addSpeciesPhaseSpaceFrom(outRequestVals.specName, from, every);
        }
        else if(!type.compare("Diag")){
          manager.addDiagFrom(from, every);
        }


      }
    }
  }
  return true;
}



void jsonParser::setOutputRequests(Json::Value &document, OUTPUT_MANAGER &manager, std::map<std::string, outDomain*>  &domains, std::vector<SPECIE*> &species){
  std::string  name1="Output";
  Json::Value outputList;
  if(setValue(outputList, document, name1.c_str() ) && outputList.isArray()){
    for(unsigned int index=0; index<outputList.size(); index++){
      Json::Value outReq = outputList[index];

      std::string type;

      bool correct = false;
      if(setString(&type, outReq, "type")&&isOutTypeOk(type)){

        jsonParser::outRequest outRequestVals;

        outRequestVals.isFrom = setDouble(&outRequestVals.from, outReq, "from");
        outRequestVals.isTo = setDouble(&outRequestVals.to, outReq, "to");
        outRequestVals.isEvery = setDouble(&outRequestVals.every, outReq, "every");
        outRequestVals.isAt = setDouble(&outRequestVals.at, outReq, "at");

        outRequestVals.isSpecName = setString(&outRequestVals.specName, outReq, "spec");
        outRequestVals.isDomainName =setString(&outRequestVals.domainName, outReq, "in");

        bool isSpecRequired=false;
        if(!type.compare("density"))
          isSpecRequired=true;
        else if(!type.compare("PhSp"))
          isSpecRequired=true;

        bool isSpecValid = !isSpecRequired || (!outRequestVals.isSpecName) || isSpecNameOk(outRequestVals.specName, species);
        bool isDomainValid = (!outRequestVals.isDomainName) || (domains.find(outRequestVals.domainName)!=domains.end());
        bool isAtLeastOneTime=(outRequestVals.isAt) || (outRequestVals.isEvery);

        if(isAtLeastOneTime && isSpecValid && isDomainValid){
          correct=addOutReq(type, manager, domains, outRequestVals);
        }

      }

      if ((!correct) && isThisJsonMaster){
        std::cout << "Output request " << index << " has unrecognized format!" << std::endl;
      }

    }
  }
}


