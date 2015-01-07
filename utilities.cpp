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


#include "utilities.h"

void moveWindow(GRID* _mygrid, EM_FIELD* _myfield, std::vector<SPECIE*> _myspecies){
  _mygrid->moveWindow();
  _myfield->moveWindow();
  for (std::vector<SPECIE*>::iterator spec_iterator = _myspecies.begin(); spec_iterator != _myspecies.end(); spec_iterator++){
    (*spec_iterator)->move_window();
  }
}

void restartFromDump(int *_dumpID, GRID* mygrid, EM_FIELD* myfield, std::vector<SPECIE*> species){
  int dumpID = _dumpID[0];
  std::ifstream dumpFile;
  MPI_Barrier(MPI_COMM_WORLD);

  if (mygrid->myid == mygrid->master_proc){
    time_t timer;
    std::time(&timer);  /* get current time; same as: timer = time(NULL)  */

    struct tm * now = localtime(&timer);

    printf("   restart from DUMP #%i ... %2.2i:%2.2i:%2.2i\n", (dumpID), now->tm_hour, now->tm_min, now->tm_sec);
    fflush(stdout);
  }
  dumpFile.open(mygrid->composeDumpFileName(dumpID).c_str());
  if (dumpFile.good()){
    mygrid->reloadDump(dumpFile);
    myfield->reloadDump(dumpFile);

    for (std::vector<SPECIE*>::iterator spec_iterator = species.begin(); spec_iterator != species.end(); spec_iterator++){
      (*spec_iterator)->reloadBigBufferDump(dumpFile);
    }
    dumpFile.close();
    dumpID++;
  }
  MPI_Barrier(MPI_COMM_WORLD);
  if (mygrid->myid == mygrid->master_proc){
    time_t timer;
    std::time(&timer);  /* get current time; same as: timer = time(NULL)  */

    struct tm * now = localtime(&timer);

    printf("  ... DONE %2.2i:%2.2i:%2.2i\n", now->tm_hour, now->tm_min, now->tm_sec);
    fflush(stdout);
  }
  _dumpID[0] = dumpID;
}

void dumpFilesForRestart(int *_dumpID, GRID* mygrid, EM_FIELD* myfield, std::vector<SPECIE*> species){
  int dumpID = _dumpID[0];
  std::ofstream dumpFile;
  dumpFile.open(mygrid->composeDumpFileName(dumpID).c_str());
  mygrid->dump(dumpFile);
  myfield->dump(dumpFile);
  for (std::vector<SPECIE*>::iterator spec_iterator = species.begin(); spec_iterator != species.end(); spec_iterator++){
    (*spec_iterator)->dumpBigBuffer(dumpFile);
  }
  dumpFile.close();
  dumpID++;
  _dumpID[0] = dumpID;
  MPI_Barrier(MPI_COMM_WORLD);
  if (mygrid->myid == mygrid->master_proc){
    printf("\t DUMP #%i done!\n", (dumpID - 1));
  }
}

void dumpDebugFilesForRestart(int *_dumpID, GRID* mygrid, EM_FIELD* myfield, std::vector<SPECIE*> species){
  int dumpID = _dumpID[0];
  std::ofstream dumpFile;
  dumpFile.open(mygrid->composeDumpFileName(dumpID).c_str());
  mygrid->debugDump(dumpFile);
  //myfield->debugDump(dumpFile);
  for (std::vector<SPECIE*>::iterator spec_iterator = species.begin(); spec_iterator != species.end(); spec_iterator++){
    (*spec_iterator)->debugDump(dumpFile);
  }
  dumpFile.close();
  dumpID++;
  _dumpID[0] = dumpID;
  MPI_Barrier(MPI_COMM_WORLD);
  if (mygrid->myid == mygrid->master_proc){
    printf("\t DUMP #%i done!\n", (dumpID - 1));
  }
}








void parseJsonInputFile2(Json::Value &root, std::string nomeFile){

  std::ifstream inputFile(nomeFile.c_str());
  std::stringstream buffer;
  buffer << inputFile.rdbuf();
  inputFile.close();
   Json::Reader reader;
   bool parsedSuccess = reader.parse(buffer.str().c_str(), root, false);

   if(not parsedSuccess){
     // Report failures and their locations
     // in the document.
     std::cout<<"Failed to parse JSON"<<std::endl
         <<reader.getFormatedErrorMessages()
         <<std::endl;
     exit(1);
   }
   const Json::Value array = root["array"];
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
  if(parent[name].isObject()){
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
  if(parent[name].isObject()){
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
  if(parent[name].isObject()){
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
  if(parent[name].isObject()){
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


int getDimensionalityFromJson2(Json::Value &parent, int defaultDimensionality){
  int dim = defaultDimensionality;
  const char* name="dimensions";
  if(( !parent[name].isNull() )){
    if( parent[name].type() == Json::intValue){
      dim = parent[name].asInt();
    }
  }
  return dim;
}
