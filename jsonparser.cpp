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


void parseJsonInputFile(Json::Value &root, std::string nomeFile){

  std::ifstream inputFile(nomeFile.c_str());
  std::stringstream buffer;
  buffer << inputFile.rdbuf();
  inputFile.close();
   Json::Reader reader;
   bool parsedSuccess = reader.parse(buffer.str().c_str(), root, false);

   if(not parsedSuccess){
     std::cout<<"Failed to parse JSON"<<std::endl
         <<reader.getFormatedErrorMessages()
         <<std::endl;
     exit(1);
   }
}


bool setIntFromJson(int *number, Json::Value &parent, const char* name){
  bool outFlag;
  if(outFlag = ( !parent[name].isNull() ) ){
    if( parent[name].isInt())
      *number = parent[name].asInt();
  }
  return outFlag;
}


bool setDoubleFromJson(double *number, Json::Value  &parent,const char* name){
  bool outFlag;
  if(outFlag = ( !parent[name].isNull() ) ){
    if( parent[name].isDouble())
      *number = parent[name].asDouble();
  }
  return outFlag;
}
bool setBoolFromJson(bool * number, Json::Value  &parent,const char* name){
  bool outFlag;
  if(outFlag = ( !parent[name].isNull() ) ){
    if( parent[name].isBool())
      *number = parent[name].asBool();
  }
  return outFlag;
}

bool setStringFromJson(std::string * number, Json::Value  &parent,const char* name){
  bool outFlag;
  if(outFlag = ( !parent[name].isNull() ) ){
    if( parent[name].isString())
      *number = parent[name].asString();
  }
  return outFlag;
}


bool setValueFromJson(Json::Value &child, Json::Value &parent, const char* name){
  bool outFlag;
  if(outFlag = ( !parent[name].isNull() ) ){
    child=parent[name];
  }
  return outFlag;
}


int getDimensionalityFromJson(Json::Value &document, int defaultDimensionality){
  int dim = defaultDimensionality;
  const char* name="dimensions";
  if(setIntFromJson(&dim, document, name))
  std::cout << "dimensions not set in JSON input file!\n";
return dim;
}
void setXrangeFromJson(Json::Value &parent,GRID *grid){
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
void setYrangeFromJson(Json::Value &parent,GRID *grid){
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
void setZrangeFromJson(Json::Value &parent,GRID *grid){
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

void setNCellsFromJson(Json::Value &parent,GRID *grid){
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
void setNprocsFromJson(Json::Value &document,GRID *grid){
  int nProcY=1, nProcZ=1;
  setIntFromJson(&nProcY, document,"nProcY");
  setIntFromJson(&nProcZ, document,"nProcZ");

  grid->setNProcsAlongY(nProcY);
  grid->setNProcsAlongZ(nProcZ);
}

void setSimulationTimeFromJson(Json::Value &document,GRID *grid){
  double simulationTime;
  setDoubleFromJson(&simulationTime,document,"simulationTime");
  grid->setSimulationTime(simulationTime);
}

void setDumpControlFromJson(Json::Value &parent, DUMP_CONTROL *myDumpControl){
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

void setStretchedGridFromJson(Json::Value &document,GRID *grid){
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
void setMovingWindowFromJson(Json::Value  &document,GRID *grid){
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

bool setLaserType(laserPulse *pulse1, Json::Value  &mylaser){
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

bool setLaserPolarization(laserPulse *pulse1, Json::Value  &mylaser){
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

bool setLaserDurationFWHM(laserPulse *pulse1, Json::Value  &mylaser){
  std::string name2="durationFWHM";
  double durationFWHM;
  bool flag=false;
  if(flag=setDoubleFromJson(&durationFWHM,mylaser,name2.c_str())){
  pulse1->setDurationFWHM(durationFWHM);
  }
  return flag;
}

bool setLaserInitialPosition(laserPulse *pulse1, Json::Value  &mylaser){
  std::string name2="initialPosition";
  double initialPosition;
  bool flag=false;
  if(flag=setDoubleFromJson(&initialPosition,mylaser,name2.c_str())){
  pulse1->setPulseInitialPosition(initialPosition);
  }
  return flag;
}

bool setLaserAmplitude(laserPulse *pulse1, Json::Value  &mylaser){
  std::string name2="a";
  double amplitude;
  bool flag=false;
  if(flag=setDoubleFromJson(&amplitude,mylaser,name2.c_str())){
  pulse1->setNormalizedAmplitude(amplitude);
  }
return flag;
}

bool setLaserWaist(laserPulse *pulse1, Json::Value  &mylaser){
  std::string name2="waist";
  double waist;
  bool flag=false;
  if(flag=setDoubleFromJson(&waist,mylaser,name2.c_str())){
  pulse1->setWaist(waist);
  }
return flag;
}
bool setLaserFocusPosition(laserPulse *pulse1, Json::Value  &mylaser){
  std::string name2="focusPosition";
  double focusPosition;
  bool flag=false;
  if(flag=setDoubleFromJson(&focusPosition,mylaser,name2.c_str())){
  pulse1->setFocusPosition(focusPosition);
  }
return flag;
}
bool setLaserLambda(laserPulse *pulse1, Json::Value  &mylaser){
  std::string name2="lambda";
  double lambda;
  bool flag=false;
  if(flag=setDoubleFromJson(&lambda,mylaser,name2.c_str())){
  pulse1->setLambda(lambda);
  }
return flag;
}

bool setLaserRotation(laserPulse *pulse1, Json::Value  &mylaser){
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

bool setLaserRiseTime(laserPulse *pulse1, Json::Value  &mylaser){
  std::string name2="riseTime";
  double riseTime;
  bool flag=false;
  if(flag=setDoubleFromJson(&riseTime,mylaser,name2.c_str())){
  pulse1->setRiseTime(riseTime);
  }
return flag;
}

struct laserPulseBoolFlags{

  bool type, pol, waist, a, lambda, duration, initialPosition, focusPosition, rotation, riseTime;
  laserPulseBoolFlags(){
    type = pol = waist = a = lambda = duration = initialPosition = focusPosition = rotation = riseTime =false;
  }
};

bool checkLaserBoolFlags(laserPulseBoolFlags flags, laserPulse *pulse){
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

void setLaserPulsesFromJson(Json::Value &document, EM_FIELD *emfield){
  std::string  name1="Laser";
  Json::Value lasers;
  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
  bool amIMasterProc = (myid==0);

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
          else if(amIMasterProc){
            std::cout << "warning: laser #" << index << " is incorrectly defined\n";
          }

          delete pulse1;
        }
      }
    }
  }
}

int findPlasmaFunction(std::string plasmaFunction){
    for (int i = 0; i < PLASMA::maxdF; i++){
        if (!plasmaFunction.compare(PLASMA::dFNames[i]))
            return i;
    }
    return -1;
}

void setPlasmasFromJson(Json::Value &document, std::map<std::string, PLASMA*> &map){
    std::string  name1="Plasma";
    Json::Value plasmas;
    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
    bool amIMasterProc = (myid==0);


    if(setValueFromJson( plasmas, document, name1.c_str() ) && plasmas.isArray()){
        for(unsigned int index=0; index<plasmas.size(); index++){
            Json::Value myPlasma = plasmas[index];

            std::string plasmaName;
            setStringFromJson(&plasmaName,myPlasma,"name");

            if(map.find(plasmaName)!=map.end()&&amIMasterProc){
                std::cout << "Warning! Plasma " << plasmaName << " is defined multiple times!" << std::endl;
            }



            std::string plasmaFunction;
            bool isPlasmaFunctionSet=false;

            int plasmaFunctionIndex;
            if(setStringFromJson(&plasmaFunction, myPlasma, "densityFunction")){
                plasmaFunctionIndex = findPlasmaFunction(plasmaFunction);
                std::cout << plasmaFunctionIndex << " " << plasmaFunction <<std::endl;
                if(plasmaFunctionIndex >= 0)
                   isPlasmaFunctionSet=true;
            }
            if(isPlasmaFunctionSet){
                map[plasmaName] = new PLASMA();
                map[plasmaName]->density_function = PLASMA::dFPoint[plasmaFunctionIndex];

                double tdouble, t2double;
                Json::Value range;

                if(setValueFromJson(range, myPlasma, "XRangeBox") && (range.size() == 2)){
                    map[plasmaName]->setXRangeBox(range[0].asDouble(), range[1].asDouble());
                }

                if(setValueFromJson(range, myPlasma, "YRangeBox") && (range.size() == 2)){
                    map[plasmaName]->setYRangeBox(range[0].asDouble(), range[1].asDouble());
                }

                if(setValueFromJson(range, myPlasma, "ZRangeBox") && (range.size() == 2)){
                    map[plasmaName]->setZRangeBox(range[0].asDouble(), range[1].asDouble());
                }

                if(setDoubleFromJson(&tdouble, myPlasma, "DensityLambda")){
                    if(setDoubleFromJson(&t2double, myPlasma, "DensityCoefficient")){
                        map[plasmaName]->setDensityCoefficient(t2double,tdouble);
                    }
                }
                else{

                    if(setDoubleFromJson(&tdouble, myPlasma, "DensityCoefficient")){
                        map[plasmaName]->setDensityCoefficient(tdouble);
                    }
                }

                if(setDoubleFromJson(&tdouble, myPlasma, "RampLength")){
                    map[plasmaName]->setRampLength(tdouble);
                }

                if(setDoubleFromJson(&tdouble, myPlasma, "RampMinDensity")){
                    map[plasmaName]->setRampMinDensity(tdouble);
                }

                if(setDoubleFromJson(&tdouble, myPlasma, "LeftRampLength")){
                    map[plasmaName]->setLeftRampLength(tdouble);
                }


                if(setDoubleFromJson(&tdouble, myPlasma, "RightRampLength")){
                    map[plasmaName]->setRightRampLength(tdouble);
                }


                if(setDoubleFromJson(&tdouble, myPlasma, "ScaleLength")){
                    map[plasmaName]->setScaleLength(tdouble);
                }

                if(setDoubleFromJson(&tdouble, myPlasma, "LeftScaleLength")){
                    map[plasmaName]->setLeftScaleLength(tdouble);
                }

                if(setDoubleFromJson(&tdouble, myPlasma, "RightScaleLength")){
                    map[plasmaName]->setRightScaleLength(tdouble);
                }

                if(setDoubleFromJson(&tdouble, myPlasma, "LeftRampMinDensity")){
                    map[plasmaName]->setLeftRampMinDensity(tdouble);
                }

                if(setDoubleFromJson(&tdouble, myPlasma, "RightRampMinDensity")){
                    map[plasmaName]->setRightRampMinDensity(tdouble);
                }


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
                if(amIMasterProc)
                     std::cout << "Warning! Plasma " << plasmaName <<
                                  " has no valid density function. It will be ignored." << std::endl;
            }

        }
    }
}

