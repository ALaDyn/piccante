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
  _mygrid->move_window();
  _myfield->move_window();
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

void parseJsonInputFile(rapidjson::Document &document, std::string nomeFile){
  //  FILE * pFile = fopen ("inputPiccante.json" , "r");
  //  rapidjson::FileStream is(pFile);
  //  rapidjson::Document document;
  //  document.ParseStream<0>(is);
  //  fclose(pFile);

  //  FILE* fp = fopen("big.json", "rb"); // non-Windows use "r"
  //  char readBuffer[65536];
  //  FileReadStream is(fp, readBuffer, sizeof(readBuffer));
  //  Document d;
  //  d.ParseStream(is);
  //  fclose(fp);
  //  std::ifstream inputFile("inputPiccante.json");
  //  std::stringstream buffer;
  //  buffer << inputFile.rdbuf();
  //  rapidjson::StringStream s(buffer.str().c_str());
  //  document.ParseStream(s);
  //  inputFile.close();

  //  int dim = DEFAULT_DIMENSIONALITY;
  //  if(document.HasMember("dimensions")){
  //    if(document["dimensions"].IsInt())
  //      dim = document["dimensions"].GetInt()
  //  }
  std::ifstream inputFile(nomeFile.c_str());
  std::stringstream buffer;
  buffer << inputFile.rdbuf();
  rapidjson::StringStream s(buffer.str().c_str());
  document.ParseStream(s);
  inputFile.close();
}

bool setIntFromJson(int *number, rapidjson::Value &document,const char* name){
  bool outFlag;
  if(outFlag=document.HasMember(name)){
    if(outFlag=document[name].IsInt())
      *number = document[name].GetInt();
  }
  return outFlag;
}
bool setDoubleFromJson(double *number, rapidjson::Value &document,const char* name){
  bool outFlag;
  if(outFlag=document.HasMember(name)){
      if(outFlag=document[name].IsNumber())
        *number = document[name].GetDouble();
    }
  return outFlag;
}
bool setBoolFromJson(bool * number, rapidjson::Value &document,const char* name){
  bool outFlag;
  if(outFlag=document.HasMember(name)){
    if(outFlag=document[name].IsBool())
      *number = document[name].GetBool();
  }
  return outFlag;
}

bool setValueFromJson(rapidjson::Value &jsonValue, rapidjson::Value &document, const char* name){
  bool outFlag;
  if(outFlag=document.HasMember(name)){
    jsonValue=document[name];
  }
  return outFlag;
}


int getDimensionalityFromJson(rapidjson::Document &document, int defaultDimensionality){
  int dim = defaultDimensionality;
  const char* name="dimensions";
  if(document.HasMember(name)){
    //if(document[name].IsInt())
      //dim = document[name].GetInt();
  }
return dim;
}
void setXrangeFromJson(rapidjson::Document &document,GRID *grid){
  double min=-1.0, max=1.0;
  const char* name="xRange";
  if(document.HasMember(name)){
    if(document[name].IsArray()){
      min=document[name][0].GetDouble();
      max=document[name][1].GetDouble();
    }
  }
  grid->setXrange(min, max);
}
void setYrangeFromJson(rapidjson::Document &document,GRID *grid){
  double min=-1.0, max=1.0;
  const char* name="yRange";
  if(document.HasMember(name)){
    if(document[name].IsArray()){
      min=document[name][0].GetDouble();
      max=document[name][1].GetDouble();
    }
  }
  grid->setYrange(min, max);
}
void setZrangeFromJson(rapidjson::Document &document,GRID *grid){
  double min=-1.0, max=1.0;
  const char* name="zRange";
  if(document.HasMember(name)){
    if(document[name].IsArray()){
      min=document[name][0].GetDouble();
      max=document[name][1].GetDouble();
    }
  }
  grid->setZrange(min, max);
}

void setNCellsFromJson(rapidjson::Document &document,GRID *grid){
  int Nx, Ny, Nz;
  Nx=Ny=Nz=1;
  const char* name="NCells";
  if(document.HasMember(name)){
    if(document[name].IsArray()){
      Nx=document[name][0].GetInt();
      Ny=document[name][1].GetInt();
      Nz=document[name][2].GetInt();
    }
  }
  grid->setNCells(Nx, Ny, Nz);

}
void setNprocsFromJson(rapidjson::Document &document,GRID *grid){
  int nProcY=1, nProcZ=1;
  setIntFromJson(&nProcY, document,"nProcY");
  setIntFromJson(&nProcZ, document,"nProcZ");

  grid->setNProcsAlongY(nProcY);
  grid->setNProcsAlongZ(nProcZ);
}

void setSimulationTimeFromJson(rapidjson::Document &document,GRID *grid){
  double simulationTime;
  setDoubleFromJson(&simulationTime,document,"simulationTime");
//  if(document.HasMember(name)){
//    if(document[name].IsNumber())
//      simulationTime = document[name].GetDouble();
//  }
  grid->setSimulationTime(simulationTime);
}

void setDumpControlFromJson(rapidjson::Document &document,DUMP_CONTROL *myDumpControl){
  myDumpControl->doRestart = false;
  myDumpControl->doDump = false;
  std::string  name1="restart";
  std::string  name2;
  if(document.HasMember(name1.c_str())){
    name2 = "doRestart";
    if(document[name1.c_str()].HasMember(name2.c_str())){
      if(document[name1.c_str()][name2.c_str()].IsBool()){
      myDumpControl->doRestart = document[name1.c_str()][name2.c_str()].GetBool();
      }
    }
    name2 = "dumpEvery";
    if(document[name1.c_str()].HasMember(name2.c_str())){
      if(document[name1.c_str()][name2.c_str()].IsNumber()){
      myDumpControl->dumpEvery = document[name1.c_str()][name2.c_str()].GetDouble();
      }
    }
    name2 = "doDump";
    if(document[name1.c_str()].HasMember(name2.c_str())){
      if(document[name1.c_str()][name2.c_str()].IsBool()){
      myDumpControl->doDump = document[name1.c_str()][name2.c_str()].GetBool();
      }
    }
    name2 = "restartFromDump";
    if(document[name1.c_str()].HasMember(name2.c_str())){
      if(document[name1.c_str()][name2.c_str()].IsInt()){
      myDumpControl->restartFromDump = document[name1.c_str()][name2.c_str()].GetInt();
      }
    }
  }

}

void setStretchedGridFromJson(rapidjson::Document &document,GRID *grid){
  std::string  name1="StretchedGrid";
  rapidjson::Value stretching;
  if(setValueFromJson(stretching, document, name1.c_str() ) ) {
    grid->enableStretchedGrid();
    std::string  name2;
    rapidjson::Value stretching1D;

    name2="x";
    if(setValueFromJson(stretching1D, stretching,name2.c_str() ) ){
      std::string  name3="left";
      rapidjson::Value stretchingLeft;
      if(setValueFromJson(stretchingLeft, stretching1D,name3.c_str() ) ){
        double limit;
        int NCells;
        setIntFromJson(&NCells,stretchingLeft,"NCells");
        setDoubleFromJson(&limit, stretchingLeft, "limit");
        grid->setXandNxLeftStretchedGrid(limit,NCells);
      }
      name3="right";
      rapidjson::Value stretchingRight;
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
      rapidjson::Value stretchingLeft;
      if(setValueFromJson(stretchingLeft, stretching1D,name3.c_str() ) ){
        double limit;
        int NCells;
        setIntFromJson(&NCells,stretchingLeft,"NCells");
        setDoubleFromJson(&limit, stretchingLeft, "limit");
        grid->setYandNyLeftStretchedGrid(limit,NCells);
      }
      name3="right";
      rapidjson::Value stretchingRight;
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
      rapidjson::Value stretchingLeft;
      if(setValueFromJson(stretchingLeft, stretching1D,name3.c_str() ) ){
        double limit;
        int NCells;
        setIntFromJson(&NCells,stretchingLeft,"NCells");
        setDoubleFromJson(&limit, stretchingLeft, "limit");
        grid->setZandNzLeftStretchedGrid(limit,NCells);
      }
      name3="right";
      rapidjson::Value stretchingRight;
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
void setMovingWindowFromJson(rapidjson::Document &document,GRID *grid){
  std::string  name1="MovingWindow";
  rapidjson::Value movingWindow;
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
