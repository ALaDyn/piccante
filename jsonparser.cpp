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

void setLaserType(laserPulse *pulse1, Json::Value  &mylaser, bool amIMasterProc, int index){
  std::string name2="type";
  std::string type;
  if(setStringFromJson(&type,mylaser,name2.c_str())){
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
    else if(amIMasterProc){
          std::cout << "WARNING: badly defined laser type for pulse #" << index << std::endl;
    }
  }
  else{
    if(amIMasterProc)
      std::cout << "WARNING: badly defined laser type for pulse #" << index << std::endl;
  }
}

void setLaserPolarization(laserPulse *pulse1, Json::Value  &mylaser, bool amIMasterProc, int index){
  std::string name2="polarization";
  std::string polarization;
  if(setStringFromJson(&polarization,mylaser,name2.c_str())){
    if(polarization=="P"){
        pulse1->setPPolarization();
      }
    else if( polarization=="S"){
        pulse1->setSPolarization();
    }
    else if(polarization=="C"){
        pulse1->setCircularPolarization();
     }
    else if(amIMasterProc){
          std::cout << "WARNING: badly defined laser polarization for pulse #" << index << std::endl;
    }
  }
  else{
    if(amIMasterProc)
      std::cout << "laser polarization not defined for pulse #" << index << std::endl;
  }

}

void setLaserDurationFWHM(laserPulse *pulse1, Json::Value  &mylaser, bool amIMasterProc, int index){
  std::string name2="durationFWHM";
  double durationFWHM;
  if(setDoubleFromJson(&durationFWHM,mylaser,name2.c_str())){
  pulse1->setDurationFWHM(durationFWHM);
  }
  else{
    if(amIMasterProc)
      std::cout << "WARNING: badly defined laser duration for pulse #" << index << std::endl;
  }

}

void setLaserInitialPosition(laserPulse *pulse1, Json::Value  &mylaser, bool amIMasterProc, int index){
  std::string name2="initialPosition";
  double initialPosition;
  if(setDoubleFromJson(&initialPosition,mylaser,name2.c_str())){
  pulse1->setPulseInitialPosition(initialPosition);
  }
  else{
    if(amIMasterProc)
      std::cout << "WARNING: badly defined laser initialPosition for pulse #" << index << std::endl;
  }
}

void setLaserAmplitude(laserPulse *pulse1, Json::Value  &mylaser, bool amIMasterProc, int index){
  std::string name2="a";
  double amplitude;
  if(setDoubleFromJson(&amplitude,mylaser,name2.c_str())){
  pulse1->setNormalizedAmplitude(amplitude);
  }
  else{
    if(amIMasterProc)
      std::cout << "WARNING: badly defined laser NormalizedAmplitude for pulse #" << index << std::endl;
  }
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
        name2="enabled";
        bool enabled=false;
        Json::Value myLaser = lasers[index];

        setBoolFromJson(&enabled, myLaser, name2.c_str());
        if(enabled){
          laserPulse *pulse1=new(laserPulse);

          setLaserType(pulse1, myLaser, amIMasterProc, index);
          setLaserPolarization(pulse1, myLaser, amIMasterProc, index);
          setLaserDurationFWHM(pulse1, myLaser, amIMasterProc, index);
          setLaserInitialPosition(pulse1, myLaser, amIMasterProc, index);
          setLaserAmplitude(pulse1, myLaser, amIMasterProc, index);

          emfield->addPulse(pulse1);
        }


      }
    }
  }
}


