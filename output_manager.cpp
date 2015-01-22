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

#include "output_manager.h"

int is_big_endian()
{
  union {
    uint32_t i;
    char c[4];
  } bint = { 0x01020304 };

  return bint.c[0] == 1;
}

emProbe::emProbe(){
  coordinates[0] = coordinates[1] = coordinates[2] = 0;
  name = "NONAME";
}
bool emProbe::compareProbes(emProbe *rhs){
  return (coordinates[0] == rhs->coordinates[0] &&
      coordinates[1] == rhs->coordinates[1] &&
      coordinates[2] == rhs->coordinates[2] &&
      name.c_str() == rhs->name.c_str());

}
void emProbe::setPointCoordinate(double X, double Y, double Z){
  coordinates[0] = X;
  coordinates[1] = Y;
  coordinates[2] = Z;
}
void emProbe::setName(std::string iname){
  name = iname;
}


outDomain::outDomain(){
  coordinates[0] = coordinates[1] = coordinates[2] = 0;
  name = "";
  remainingCoord[0] = remainingCoord[1] = 1;
  remainingCoord[2] = 1;
  subselection = false;
  overrideFlag = false;
  followMovingWindowFlag = false;
  rmin[0] = rmin[1] = rmin[2] = -1e10;
  rmax[0] = rmax[1] = rmax[2] = +1e10;

}

void outDomain::followMovingWindow(){
  followMovingWindowFlag = true;
}

bool outDomain::compareDomains(outDomain *rhs){
  if (coordinates[0] == rhs->coordinates[0] &&
      coordinates[1] == rhs->coordinates[1] &&
      coordinates[2] == rhs->coordinates[2] &&
      remainingCoord[0] == rhs->remainingCoord[0] &&
      remainingCoord[1] == rhs->remainingCoord[1] &&
      remainingCoord[2] == rhs->remainingCoord[2] &&
      name.c_str() == rhs->name.c_str() &&
      subselection == rhs->subselection&&
      followMovingWindowFlag == rhs->followMovingWindowFlag&&
      overrideFlag == rhs->overrideFlag){
    if (!subselection)
      return true;
    else if (rmin[0] == rhs->rmin[0] &&
             rmin[1] == rhs->rmin[1] &&
             rmin[2] == rhs->rmin[2] &&
             rmax[0] == rhs->rmax[0] &&
             rmax[1] == rhs->rmax[1] &&
             rmax[2] == rhs->rmax[2]){
      return true;
    }
  }
  return false;
}
void outDomain::setFreeDimensions(bool flagX, bool flagY, bool flagZ){
  remainingCoord[0] = flagX;
  remainingCoord[1] = flagY;
  remainingCoord[2] = flagZ;
}
void outDomain::setPointCoordinate(double X, double Y, double Z){
  coordinates[0] = X;
  coordinates[1] = Y;
  coordinates[2] = Z;
}
void outDomain::setName(std::string iname){
  name = iname;
}
void outDomain::setXRange(double min, double max){
  subselection = true;
  rmin[0] = min;
  rmax[0] = max;
}

void outDomain::setYRange(double min, double max){
  subselection = true;
  rmin[1] = min;
  rmax[1] = max;
}

void outDomain::setZRange(double min, double max){
  subselection = true;
  rmin[2] = min;
  rmax[2] = max;
}

bool requestCompTime(const request &first, const request &second){
  if (first.itime != second.itime)
    return (first.itime < second.itime);
  else if (first.type != second.type)
    return (first.type < second.type);
  else if (first.target != second.target)
    return (first.target < second.target);
  else if (first.domain != second.domain)
    return (first.domain < second.domain);
  else
    return false;

}

bool requestCompUnique(const request &first, const request &second){
  return ((first.itime == second.itime) &&
          (first.type == second.type) &&
          (first.target == second.target) &&
          (first.domain == second.domain));
}
bool compOutput(const reqOutput &first, const reqOutput &second){
  if (first.p != second.p)
    return (first.p < second.p);
  else if (first.task != second.task)
    return (first.task < second.task);
  else
    return false;
}


OUTPUT_MANAGER::OUTPUT_MANAGER(GRID* _mygrid, EM_FIELD* _myfield, CURRENT* _mycurrent, std::vector<SPECIE*> _myspecies){
  mygrid = _mygrid;
  isThereGrid = (mygrid != NULL) ? true : false;

  myfield = _myfield;
  isThereField = (myfield != NULL) ? true : false;

  mycurrent = _mycurrent;
  isThereCurrent = (mycurrent != NULL) ? true : false;

  myspecies = _myspecies;
  isThereSpecList = (myspecies.size() > 0) ? true : false;

  amIInit = false;

  isThereDiag = false;
  isThereEMProbe =false;

  outDomain *domain1 = new outDomain;
  domain1->overrideFlag = true;
  myDomains.push_back(domain1);

  fieldGroupSize = FIELDS_GROUP_SIZE;
  multifileGroupSize = MACRO_CPUGROUP_FOR_MULTIFILE;
  particleGroupSize = PHASE_SPACE_GROUP_SIZE;
  particleBufferSize = NPARTICLE_BUFFER_SIZE;

  outputDir = "OUTPUT";
}

OUTPUT_MANAGER::~OUTPUT_MANAGER(){
  for (std::vector<outDomain*>::iterator it = myDomains.begin(); it != myDomains.end(); ++it)
    delete(*it);
}

void OUTPUT_MANAGER::setOutputPath(std::string dirName){
  outputDir = dirName;
}

void OUTPUT_MANAGER::createDiagFile(){
  if (checkGrid()){
    if (mygrid->myid != mygrid->master_proc)
      return;
  }
  else{
    if (mygrid->myid != 0)
      return;
  }

  std::ofstream of1;
  diagFileName = outputDir + "/diag.dat";
  of1.open(diagFileName.c_str(), std::ofstream::out | std::ofstream::trunc);
  of1 << " " << std::setw(diagNarrowWidth) << "#istep" << " " << std::setw(diagWidth) << "time";
  of1 << " " << std::setw(diagWidth) << "Etot";
  of1 << " " << std::setw(diagWidth) << "Ex2" << " " << std::setw(diagWidth) << "Ey2" << " " << std::setw(diagWidth) << "Ez2";
  of1 << " " << std::setw(diagWidth) << "Bx2" << " " << std::setw(diagWidth) << "By2" << " " << std::setw(diagWidth) << "Bz2";

  std::vector<SPECIE*>::const_iterator spec_iterator;

  for (spec_iterator = myspecies.begin(); spec_iterator != myspecies.end(); spec_iterator++){
    of1 << " " << std::setw(diagWidth) << (*spec_iterator)->name;
  }

  of1 << std::endl;
  of1.close();
}

void OUTPUT_MANAGER::createExtremaFiles(){
  //  if (checkGrid()){
  if (mygrid->myid != mygrid->master_proc)
    return;

  if (checkEMField()){
    extremaFieldFileName = outputDir + "/EXTREMES_EMfield.dat";
    std::ofstream of0;
    of0.open(extremaFieldFileName.c_str());
    myfield->init_output_extrems(of0);
    of0.close();
  }

  if (checkSpecies()){
    for (std::vector<SPECIE*>::iterator it = myspecies.begin(); it != myspecies.end(); it++){
      std::stringstream ss1;
      ss1 << outputDir << "/EXTREMES_" << (*it)->name << ".dat";
      extremaSpecFileNames.push_back(ss1.str());
      std::ofstream of1;
      of1.open(ss1.str().c_str());
      (*it)->init_output_extrems(of1);
      of1.close();
    }
  }
}

void OUTPUT_MANAGER::createEMProbeFiles(){
  int ii = 0;
  for (std::vector<emProbe*>::iterator it = myEMProbes.begin(); it != myEMProbes.end(); it++){
    std::stringstream ss0;
    ss0 << outputDir << "/EMProbe_" << (*it)->name << "_" << ii << ".txt";
    (*it)->fileName = ss0.str();
    ii++;
  }
  if (checkGrid()){
    if (mygrid->myid != mygrid->master_proc)
      return;
  }
  else{
    if (mygrid->myid != 0)
      return;
  }

  if (checkEMField()){
    for (std::vector<emProbe*>::iterator it = myEMProbes.begin(); it != myEMProbes.end(); it++){
      std::ofstream of0;
      of0.open((*it)->fileName.c_str());
      of0 << "# EM field probe at coordinates  ";
      of0 << (*it)->coordinates[0] << "  ";
      of0 << (*it)->coordinates[1] << "  ";
      of0 << (*it)->coordinates[2] << "  ";
      of0 << std::endl;
      of0.close();
    }
  }
}

void OUTPUT_MANAGER::initialize(std::string _outputDir){
#if defined (USE_BOOST)
  std::string _newoutputDir;
  std::stringstream ss;
  time_t timer;
  std::time(&timer);
  ss << _outputDir << "_" << (int)timer;
  _newoutputDir = ss.str();
  if (mygrid->myid == mygrid->master_proc){
    if (!boost::filesystem::exists(_outputDir)){
      boost::filesystem::create_directories(_outputDir);
    }
    else{
      boost::filesystem::rename(_outputDir, _newoutputDir);
      boost::filesystem::create_directories(_outputDir);
    }
  }
#endif
  outputDir = _outputDir;
  prepareOutputMap();

  if (!checkGrid()){
    return;
  }

  if (isThereDiag){
    createDiagFile();
    createExtremaFiles();
  }
  if (isThereEMProbe){
    createEMProbeFiles();
  }

  amIInit = true;
}
void OUTPUT_MANAGER::initialize(){
#if defined (USE_BOOST)
  std::string _newoutputDir;
  std::stringstream ss;
  time_t timer;
  std::time(&timer);
  ss << outputDir << "_" << (int)timer;
  _newoutputDir = ss.str();
  if (mygrid->myid == mygrid->master_proc){
    if (!boost::filesystem::exists(outputDir)){
      boost::filesystem::create_directories(outputDir);
    }
    else{
      boost::filesystem::rename(outputDir, _newoutputDir);
      boost::filesystem::create_directories(outputDir);
    }
  }
#endif

  prepareOutputMap();

  if (!checkGrid()){
    return;
  }

  if (isThereDiag){
    createDiagFile();
    createExtremaFiles();
  }
  if (isThereEMProbe){
    createEMProbeFiles();
  }

  amIInit = true;
}

void OUTPUT_MANAGER::close(){

}

bool OUTPUT_MANAGER::checkGrid(){
  if (!isThereGrid){
    if (mygrid->myid == mygrid->master_proc){
      std::cout << "WARNING! No valid GRID pointer provided. Output request will be ignored." << std::endl;
    }
  }
  return isThereGrid;
}

bool OUTPUT_MANAGER::checkEMField(){
  if (!isThereField){
    if (mygrid->myid == mygrid->master_proc){
      std::cout << "WARNING! No valid FIELD pointer provided. Output request will be ignored." << std::endl;
    }
  }
  return isThereField;
}

bool OUTPUT_MANAGER::checkCurrent(){
  if (!isThereCurrent){
    if (mygrid->myid  == mygrid->master_proc){
      std::cout << "WARNING! No valid CURRENT pointer provided. Output request will be ignored." << std::endl;
    }
  }
  return isThereCurrent;
}

bool OUTPUT_MANAGER::checkSpecies(){
  if (!isThereSpecList){
    if (mygrid->myid == mygrid->master_proc){
      std::cout << "WARNING! Empty SPECIES list provided. Output request will be ignored." << std::endl;
    }
  }
  return isThereSpecList;
}

int OUTPUT_MANAGER::getIntegerTime(double dtime){
  if (dtime == 0.0)
    return 0;

  int Nsteps = mygrid->getTotalNumberOfTimesteps();
  int step = (int)floor(dtime / mygrid->dt + 0.5);

  return (step <= Nsteps) ? (step) : (-1);
}

int OUTPUT_MANAGER::findSpecIndexInMyspeciesVector(std::string name){
  int pos = 0;
  for (std::vector<SPECIE*>::iterator it = myspecies.begin(); it != myspecies.end(); it++){
    if ((*it)->getName() == name){
      return pos;
    }
    pos++;
  }
  if (mygrid->myid == mygrid->master_proc){
    std::cout << "WARNING! Species" << name << " not in species list !" << std::endl;
  }
  return -1;
}
int OUTPUT_MANAGER::findProbeIndexInMyprobeVector(emProbe *newProbe){
  int pos = 0;
  for (std::vector<emProbe*>::iterator it = myEMProbes.begin(); it != myEMProbes.end(); it++){
    if ((*it)->compareProbes(newProbe)){
      return pos;
    }
    pos++;
  }
  return -1;

}
int OUTPUT_MANAGER::findDomainIndexInMydomainsVector(outDomain *newDomain){
  int pos = 0;
  for (std::vector<outDomain*>::iterator it = myDomains.begin(); it != myDomains.end(); it++){
    if ((*it)->compareDomains(newDomain)){
      return pos;
    }
    pos++;
  }
  return -1;

}

void OUTPUT_MANAGER::addRequestToList(std::list<request>& reqList, diagType type, int target, int domain, double startTime, double frequency, double endTime){
  std::list<request> tempList;

  for (double ttime = startTime; ttime <= endTime; ttime += frequency){
    request req;
    req.dtime = ttime;
    req.itime = getIntegerTime(ttime);
    req.type = type;
    req.target = target;
    req.domain = domain;
    tempList.push_back(req);
  }
  tempList.sort(requestCompTime);
  reqList.merge(tempList, requestCompTime);
}


void OUTPUT_MANAGER::addEBFieldFrom(double startTime, double frequency){
  addEBField(NULL, startTime, frequency, mygrid->getTotalTime(), WHICH_E_AND_B);
}

void OUTPUT_MANAGER::addEBFieldAt(double atTime){
  addEBField(NULL, atTime, 1.0, atTime, WHICH_E_AND_B);
}

void OUTPUT_MANAGER::addEBFieldFromTo(double startTime, double frequency, double endTime){
  addEBField(NULL, startTime, frequency, endTime, WHICH_E_AND_B);
}

void OUTPUT_MANAGER::addEBFieldFrom(outDomain* _domain, double startTime, double frequency){
  addEBField(_domain, startTime, frequency, mygrid->getTotalTime(), WHICH_E_AND_B);
}

void OUTPUT_MANAGER::addEBFieldAt(outDomain* _domain, double atTime){
  addEBField(_domain, atTime, 1.0, atTime, WHICH_E_AND_B);
}

void OUTPUT_MANAGER::addEBFieldFromTo(outDomain* _domain, double startTime, double frequency, double endTime){
  addEBField(_domain, startTime, frequency, endTime, WHICH_E_AND_B);
}

void OUTPUT_MANAGER::addEFieldFrom(double startTime, double frequency){
  addEBField(NULL, startTime, frequency, mygrid->getTotalTime(), WHICH_E_ONLY);
}

void OUTPUT_MANAGER::addEFieldAt(double atTime){
  addEBField(NULL, atTime, 1.0, atTime, WHICH_E_ONLY);
}

void OUTPUT_MANAGER::addEFieldFromTo(double startTime, double frequency, double endTime){
  addEBField(NULL, startTime, frequency, endTime, WHICH_E_ONLY);
}

void OUTPUT_MANAGER::addBFieldFrom(double startTime, double frequency){
  addEBField(NULL, startTime, frequency, mygrid->getTotalTime(), WHICH_B_ONLY);
}

void OUTPUT_MANAGER::addBFieldAt(double atTime){
  addEBField(NULL, atTime, 1.0, atTime, WHICH_B_ONLY);
}

void OUTPUT_MANAGER::addBFieldFromTo(double startTime, double frequency, double endTime){
  addEBField(NULL, startTime, frequency, endTime, WHICH_B_ONLY);
}
// SELECTION
void OUTPUT_MANAGER::addEFieldFrom(outDomain* _domain, double startTime, double frequency){
  addEBField(_domain, startTime, frequency, mygrid->getTotalTime(), WHICH_E_ONLY);
}

void OUTPUT_MANAGER::addEFieldAt(outDomain* _domain, double atTime){
  addEBField(_domain, atTime, 1.0, atTime, WHICH_E_ONLY);
}

void OUTPUT_MANAGER::addEFieldFromTo(outDomain* _domain, double startTime, double frequency, double endTime){
  addEBField(_domain, startTime, frequency, endTime, WHICH_E_ONLY);
}

void OUTPUT_MANAGER::addBFieldFrom(outDomain* _domain, double startTime, double frequency){
  addEBField(_domain, startTime, frequency, mygrid->getTotalTime(), WHICH_B_ONLY);
}

void OUTPUT_MANAGER::addBFieldAt(outDomain* _domain, double atTime){
  addEBField(_domain, atTime, 1.0, atTime, WHICH_B_ONLY);
}

void OUTPUT_MANAGER::addBFieldFromTo(outDomain* _domain, double startTime, double frequency, double endTime){
  addEBField(_domain, startTime, frequency, endTime, WHICH_B_ONLY);
}

void OUTPUT_MANAGER::addEBField(outDomain* _domain, double startTime, double frequency, double endTime, whichFieldOut whichOut){
  if (!(checkGrid() && checkEMField()))
    return;

  int domainID = 0;
  if (_domain != NULL){
    domainID = findDomainIndexInMydomainsVector(_domain);
    if (domainID < 0){
      myDomains.push_back(_domain);
      domainID = myDomains.size() - 1;
    }
  }

  if (whichOut == WHICH_E_AND_B){
    addRequestToList(requestList, OUT_E_FIELD, 0, domainID, startTime, frequency, endTime);
    addRequestToList(requestList, OUT_B_FIELD, 0, domainID, startTime, frequency, endTime);
  }
  else if (whichOut == WHICH_E_ONLY){
    addRequestToList(requestList, OUT_E_FIELD, 0, domainID, startTime, frequency, endTime);
  }
  else if (whichOut == WHICH_B_ONLY){
    addRequestToList(requestList, OUT_B_FIELD, 0, domainID, startTime, frequency, endTime);
  }

}

// EM Probe ///////////////////////////////////////
void OUTPUT_MANAGER::addEBFieldProbeFrom(emProbe* Probe, double startTime, double frequency){
  addEBFieldProbe(Probe, startTime, frequency, mygrid->getTotalTime());
}

void OUTPUT_MANAGER::addEBFieldProbeAt(emProbe* Probe, double atTime){
  addEBFieldProbe(Probe, atTime, 1.0, atTime);
}

void OUTPUT_MANAGER::addEBFieldProbeFromTo(emProbe* Probe, double startTime, double frequency, double endTime){
  addEBFieldProbe(Probe, startTime, frequency, endTime);
}
void OUTPUT_MANAGER::addEBFieldProbe(emProbe* Probe, double startTime, double frequency, double endTime){
  if (!(checkGrid() && checkEMField()))
    return;
  int domainID = findProbeIndexInMyprobeVector(Probe);
  if (domainID < 0){
    myEMProbes.push_back(Probe);
    isThereEMProbe = true;
    domainID = myEMProbes.size() - 1;
  }
  addRequestToList(requestList, OUT_EB_PROBE, 0, domainID, startTime, frequency, endTime);
}

// ++++++++++++++++++++++++++++     species density
void OUTPUT_MANAGER::addSpeciesDensityFrom(std::string name, double startTime, double frequency){
  addSpeciesDensity(NULL, name, startTime, frequency, mygrid->getTotalTime());
}

void OUTPUT_MANAGER::addSpeciesDensityAt(std::string name, double atTime){
  addSpeciesDensity(NULL, name, atTime, 1.0, atTime);
}

void OUTPUT_MANAGER::addSpeciesDensityFromTo(std::string name, double startTime, double frequency, double endTime){
  addSpeciesDensity(NULL, name, startTime, frequency, endTime);
}

void OUTPUT_MANAGER::addSpeciesDensityFrom(outDomain* domain_in, std::string name, double startTime, double frequency){
  addSpeciesDensity(domain_in, name, startTime, frequency, mygrid->getTotalTime());
}

void OUTPUT_MANAGER::addSpeciesDensityAt(outDomain* domain_in, std::string name, double atTime){
  addSpeciesDensity(domain_in, name, atTime, 1.0, atTime);
}

void OUTPUT_MANAGER::addSpeciesDensityFromTo(outDomain* domain_in, std::string name, double startTime, double frequency, double endTime){
  addSpeciesDensity(domain_in, name, startTime, frequency, endTime);
}

void OUTPUT_MANAGER::addSpeciesDensity(outDomain* domain_in, std::string name, double startTime, double frequency, double endTime){
  if (!(checkGrid() && checkSpecies() && checkCurrent()))
    return;
  int specNum = findSpecIndexInMyspeciesVector(name);
  if (specNum < 0)
    return;
  int domainID = 0;
  if (domain_in != NULL){
    domainID = findDomainIndexInMydomainsVector(domain_in);
    if (domainID < 0){
      myDomains.push_back(domain_in);
      domainID = myDomains.size() - 1;
    }
  }
  addRequestToList(requestList, OUT_SPEC_DENSITY, specNum, domainID, startTime, frequency, endTime);
}

// ++++++++++++++++++++++++++++     current
void OUTPUT_MANAGER::addCurrentFrom(double startTime, double frequency){
  addCurrent(NULL, startTime, frequency, mygrid->getTotalTime());
}

void OUTPUT_MANAGER::addCurrentAt(double atTime){
  addCurrent(NULL, atTime, 1.0, atTime);
}

void OUTPUT_MANAGER::addCurrentFromTo(double startTime, double frequency, double endTime){
  addCurrent(NULL, startTime, frequency, endTime);
}

void OUTPUT_MANAGER::addCurrentFrom(outDomain* domain_in, double startTime, double frequency){
  addCurrent(domain_in, startTime, frequency, mygrid->getTotalTime());
}
void OUTPUT_MANAGER::addCurrentAt(outDomain* domain_in, double atTime){
  addCurrent(domain_in, atTime, 1.0, atTime);
}

void OUTPUT_MANAGER::addCurrentFromTo(outDomain* domain_in, double startTime, double frequency, double endTime){
  addCurrent(domain_in, startTime, frequency, endTime);
}

void OUTPUT_MANAGER::addCurrent(outDomain* domain_in, double startTime, double frequency, double endTime){
  if (!(checkGrid() && checkCurrent()))
    return;
  int domainID = 0;

  if (domain_in != NULL){

    domainID = findDomainIndexInMydomainsVector(domain_in);
    if (domainID < 0){
      myDomains.push_back(domain_in);
      domainID = myDomains.size() - 1;
    }
  }
  addRequestToList(requestList, OUT_CURRENT, 0, domainID, startTime, frequency, endTime);
}

// ++++++++++++++++++++++++++++     species binary
void OUTPUT_MANAGER::addSpeciesPhaseSpaceFrom(std::string name, double startTime, double frequency){
  addSpeciesPhaseSpace(NULL, name, startTime, frequency, mygrid->getTotalTime());
}

void OUTPUT_MANAGER::addSpeciesPhaseSpaceAt(std::string name, double atTime){
  addSpeciesPhaseSpace(NULL, name, atTime, 1.0, atTime);
}

void OUTPUT_MANAGER::addSpeciesPhaseSpaceFromTo(std::string name, double startTime, double frequency, double endTime){
  addSpeciesPhaseSpace(NULL, name, startTime, frequency, endTime);
}

void OUTPUT_MANAGER::addSpeciesPhaseSpaceFrom(outDomain* domain_in, std::string name, double startTime, double frequency){
  addSpeciesPhaseSpace(domain_in, name, startTime, frequency, mygrid->getTotalTime());
}

void OUTPUT_MANAGER::addSpeciesPhaseSpaceAt(outDomain* domain_in, std::string name, double atTime){
  addSpeciesPhaseSpace(domain_in, name, atTime, 1.0, atTime);
}

void OUTPUT_MANAGER::addSpeciesPhaseSpaceFromTo(outDomain* domain_in, std::string name, double startTime, double frequency, double endTime){
  addSpeciesPhaseSpace(domain_in, name, startTime, frequency, endTime);
}

void OUTPUT_MANAGER::addSpeciesPhaseSpace(outDomain* domain_in, std::string name, double startTime, double frequency, double endTime){
  if (!(checkGrid() && checkSpecies()))
    return;
  int specNum = findSpecIndexInMyspeciesVector(name);
  if (specNum < 0)
    return;
  int domainID = 0;

  if (domain_in != NULL){
    domainID = findDomainIndexInMydomainsVector(domain_in);
    if (domainID < 0){
      myDomains.push_back(domain_in);
      domainID = myDomains.size() - 1;
    }
  }
  addRequestToList(requestList, OUT_SPEC_PHASE_SPACE, specNum, domainID, startTime, frequency, endTime);
}

// ++++++++++++++++++++++++++++     diag

void OUTPUT_MANAGER::addDiagFrom(double startTime, double frequency){
  addDiag(startTime, frequency, mygrid->getTotalTime());
}

void OUTPUT_MANAGER::addDiagAt(double atTime){
  addDiag(atTime, 1.0, atTime);
}

void OUTPUT_MANAGER::addDiagFromTo(double startTime, double frequency, double endTime){
  addDiag(startTime, frequency, endTime);
}

void OUTPUT_MANAGER::addDiag(double startTime, double frequency, double endTime){
  if (!(checkGrid() && checkCurrent() && checkEMField()))
    return;
  addRequestToList(requestList, OUT_DIAG, 0, 0, startTime, frequency, endTime);
  isThereDiag = true;
}

void OUTPUT_MANAGER::prepareOutputMap(){

  requestList.sort(requestCompTime);
  requestList.unique(requestCompUnique);

  std::map< int, std::vector<request> >::iterator itMap;

  for (std::list<request>::iterator itList = requestList.begin(); itList != requestList.end(); itList++){
    itMap = allOutputs.find(itList->itime);
    if (itMap != allOutputs.end()){
      itMap->second.push_back(*itList);
    }
    else{
      std::vector<request> newTimeVec;
      newTimeVec.push_back(*itList);
      allOutputs.insert(std::pair<int, std::vector<request> >(itList->itime, newTimeVec));
      timeList.push_back(itList->itime);
    }
  }

}


void OUTPUT_MANAGER::processOutputEntry(request req){
  switch (req.type){

    case OUT_E_FIELD:
      callEMFieldDomain(req);
      break;
    case OUT_B_FIELD:
      callEMFieldDomain(req);
      break;

    case OUT_EB_PROBE:
      callEMFieldProbe(req);
      break;

    case OUT_SPEC_DENSITY:
      callSpecDensity(req);
      break;

    case OUT_CURRENT:
      callCurrent(req);
      break;

    case OUT_SPEC_PHASE_SPACE:
      callSpecPhaseSpace(req);
      break;

    case OUT_DIAG:
      callDiag(req);
      break;

    default:
      break;
  }
}

void OUTPUT_MANAGER::callDiags(int istep){
  std::map< int, std::vector<request> >::iterator itMap;
  itMap = allOutputs.find(istep);

  if (itMap == allOutputs.end())
    return;

  std::vector<request> diagList = itMap->second;

  for (std::vector<request>::iterator it = diagList.begin(); it != diagList.end(); it++){
    processOutputEntry(*it);
  }
}

std::string OUTPUT_MANAGER::composeOutputName(std::string dir, std::string out, std::string opt, double time, std::string ext){
  std::stringstream ss;
  ss << dir << "/"
     << out << "_";
  if (opt != "")
    ss << opt << "_";
  ss << std::setw(OUTPUT_SIZE_TIME_DIAG) << std::setfill('0') << std::fixed << std::setprecision(3) << time;
  ss << ext;
  return ss.str();
}

std::string OUTPUT_MANAGER::composeOutputName(std::string dir, std::string out, std::string opt1, std::string opt2, int domain, double time, std::string ext){
  std::stringstream ss;
  ss << dir << "/"
     << out << "_";
  if (opt1 != "")
    ss << opt1 << "_";
  if (opt2 != "")
    ss << opt2 << "_";
  if (domain != 0)
    ss << domain << "_";
  ss << std::setw(OUTPUT_SIZE_TIME_DIAG) << std::setfill('0') << std::fixed << std::setprecision(3) << time;
  ss << ext;
  return ss.str();
}



#if defined(USE_HDF5)

void OUTPUT_MANAGER::writeEMFieldBinaryHDF5(std::string fileName, request req){
  int dimensionality = mygrid->accesso.dimensions;
  int Ncomp = myfield->getNcomp();

  MPI_Info info = MPI_INFO_NULL;
  char nomi[6][3] = { "Ex", "Ey", "Ez", "Bx", "By", "Bz" };
  char* nomefile = new char[fileName.size() + 4];
  nomefile[fileName.size()] = 0;
  sprintf(nomefile, "%s.h5", fileName.c_str());

  hid_t       file_id, dset_id[6];         /* file and dataset identifiers */
  hid_t       filespace[6], memspace[6];      /* file and memory dataspace identifiers */
  hsize_t     dimsf[3];                 /* dataset dimensions */
  hsize_t	count[3];	          /* hyperslab selection parameters */
  hsize_t	offset[3];
  hid_t	plist_id;                 /* property list identifier */
  int         i;
  herr_t	status;
  /*
       * Set up file access property list with parallel I/O access
       */
  plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, info);

  /*
       * Create a new file collectively and release property list identifier.
       */
  file_id = H5Fcreate(nomefile, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
  H5Pclose(plist_id);
  /*
           * Create the dataspace for the dataset.
           */
  if (dimensionality == 1){
    dimsf[0] = mygrid->uniquePoints[0];
  }
  else if (dimensionality == 2){
    dimsf[0] = mygrid->uniquePoints[1];
    dimsf[1] = mygrid->uniquePoints[0];
  }
  else{
    dimsf[0] = mygrid->uniquePoints[2];
    dimsf[1] = mygrid->uniquePoints[1];
    dimsf[2] = mygrid->uniquePoints[0];

  }
  for (int i = 0; i < 6; i++)
    filespace[i] = H5Screate_simple(dimensionality, dimsf, NULL);

  /*
       * Create the dataset with default properties and close filespace.
       */
  for (int i = 0; i < 6; i++){
    dset_id[i] = H5Dcreate1(file_id, nomi[i], H5T_NATIVE_FLOAT, filespace[i],
                            H5P_DEFAULT);
    H5Sclose(filespace[i]);
  }
  /*
          * Each process defines dataset in memory and writes it to the hyperslab
          * in the file.
          */
  if (dimensionality == 1){
    offset[0] = mygrid->rproc_imin[0][mygrid->rmyid[0]];
    count[0] = mygrid->uniquePointsloc[0];
  }
  else if (dimensionality == 2){
    offset[1] = mygrid->rproc_imin[0][mygrid->rmyid[0]];
    offset[0] = mygrid->rproc_imin[1][mygrid->rmyid[1]];
    count[1] = mygrid->uniquePointsloc[0];
    count[0] = mygrid->uniquePointsloc[1];
  }
  else{
    offset[2] = mygrid->rproc_imin[0][mygrid->rmyid[0]];
    offset[1] = mygrid->rproc_imin[1][mygrid->rmyid[1]];
    offset[0] = mygrid->rproc_imin[2][mygrid->rmyid[2]];
    count[2] = mygrid->uniquePointsloc[0];
    count[1] = mygrid->uniquePointsloc[1];
    count[0] = mygrid->uniquePointsloc[2];
  }
  /*
       * Select hyperslab in the file.
       */
  for (int i = 0; i < 6; i++){
    memspace[i] = H5Screate_simple(dimensionality, count, NULL);
    filespace[i] = H5Dget_space(dset_id[i]);
    H5Sselect_hyperslab(filespace[i], H5S_SELECT_SET, offset, NULL, count, NULL);
  }
  //+++++++++++ Start CPU Field Values  +++++++++++++++++++++
  {
    float *Ex, *Ey, *Ez, *Bx, *By, *Bz;
    int Nx, Ny, Nz;
    Nx = mygrid->uniquePointsloc[0];
    Ny = mygrid->uniquePointsloc[1];
    Nz = mygrid->uniquePointsloc[2];
    int size = Nx*Ny*Nz;
    Ex = new float[size];
    Ey = new float[size];
    Ez = new float[size];
    Bx = new float[size];
    By = new float[size];
    Bz = new float[size];
    for (int k = 0; k < Nz; k++)
      for (int j = 0; j < Ny; j++)
        for (int i = 0; i < Nx; i++){
          Ex[i + j*Nx + k*Ny*Nx] = (float)myfield->VEB(0, i, j, k);
          Ey[i + j*Nx + k*Ny*Nx] = (float)myfield->VEB(1, i, j, k);
          Ez[i + j*Nx + k*Ny*Nx] = (float)myfield->VEB(2, i, j, k);
          Bx[i + j*Nx + k*Ny*Nx] = (float)myfield->VEB(3, i, j, k);
          By[i + j*Nx + k*Ny*Nx] = (float)myfield->VEB(4, i, j, k);
          Bz[i + j*Nx + k*Ny*Nx] = (float)myfield->VEB(5, i, j, k);
        }
    /*
          * Create property list for collective dataset write.
          */
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    status = H5Dwrite(dset_id[0], H5T_NATIVE_FLOAT, memspace[0], filespace[0],
        plist_id, Ex);
    status = H5Dwrite(dset_id[1], H5T_NATIVE_FLOAT, memspace[1], filespace[1],
        plist_id, Ey);
    status = H5Dwrite(dset_id[2], H5T_NATIVE_FLOAT, memspace[2], filespace[2],
        plist_id, Ez);
    status = H5Dwrite(dset_id[3], H5T_NATIVE_FLOAT, memspace[3], filespace[3],
        plist_id, Bx);
    status = H5Dwrite(dset_id[4], H5T_NATIVE_FLOAT, memspace[4], filespace[4],
        plist_id, By);
    status = H5Dwrite(dset_id[5], H5T_NATIVE_FLOAT, memspace[5], filespace[5],
        plist_id, Bz);

    delete[]Ex;
    delete[]Ey;
    delete[]Ez;
    delete[]Bx;
    delete[]By;
    delete[]Bz;
  }

  for (int i = 0; i < 6; i++){
    H5Dclose(dset_id[i]);
    H5Sclose(filespace[i]);
    H5Sclose(memspace[i]);
  }
  H5Pclose(plist_id);
  H5Fclose(file_id);
  delete[] nomefile;
  //////////////////////////// END of collective binary file write
}
#endif

void OUTPUT_MANAGER::prepareIntegerBigHeader(int *itodo, int uniqueN[3], int slice_rNproc[], int Ncomp){
  itodo[0] = is_big_endian();
  itodo[1] = uniqueN[0];
  itodo[2] = uniqueN[1];
  itodo[3] = uniqueN[2];
  itodo[4] = slice_rNproc[0];
  itodo[5] = slice_rNproc[1];
  itodo[6] = slice_rNproc[2];
  itodo[7] = Ncomp;
}

void OUTPUT_MANAGER::prepareFloatCoordinatesHeader(float *fcir[3], int uniqueN[3], int imin[]){
  for (int c = 0; c < 3; c++){
    for (int m = 0; m < uniqueN[c]; m++){
      fcir[c][m] = (float)mygrid->cir[c][m + imin[c]];
    }
  }
}
void OUTPUT_MANAGER::writeBigHeader(MPI_File thefile, int uniqueN[3], int imin[3], int slice_rNproc[3], int Ncomp){
  MPI_Status status;
  int itodo[8];
  float *fcir[3];
  for (int c = 0; c < 3; c++){
    fcir[c] = new float[uniqueN[c]];
  }
  prepareIntegerBigHeader(itodo, uniqueN, slice_rNproc, Ncomp);
  prepareFloatCoordinatesHeader(fcir, uniqueN, imin);
#ifndef DEBUG_NO_MPI_FILE_WRITE
  MPI_File_write(thefile, itodo, 8, MPI_INT, &status);
#endif
  for (int c = 0; c < 3; c++){
#ifndef DEBUG_NO_MPI_FILE_WRITE
    MPI_File_write(thefile, fcir[c], uniqueN[c], MPI_FLOAT, &status);
#endif
  }

  for (int c = 0; c < 3; c++){
    delete[] fcir[c];
  }

}
void OUTPUT_MANAGER::writeBigHeaderSingleFile(std::string fileName, int uniqueN[3], int imin[3], int slice_rNproc[3], int Ncomp){
  int itodo[8];
  float *fcir[3];
  std::ofstream thefile;
  thefile.open(fileName.c_str(), std::ios::app);
  for (int c = 0; c < 3; c++){
    fcir[c] = new float[uniqueN[c]];
  }
  prepareIntegerBigHeader(itodo, uniqueN, slice_rNproc, Ncomp);
  prepareFloatCoordinatesHeader(fcir, uniqueN, imin);

  thefile.write((char*)itodo, 8 * sizeof(int));
  for (int c = 0; c < 3; c++){
    thefile.write((char*)fcir[c], uniqueN[c] * sizeof(float));
  }

  for (int c = 0; c < 3; c++){
    delete[] fcir[c];
  }
  thefile.close();
}

void OUTPUT_MANAGER::prepareIntegerSmallHeader(int *itodo, int uniqueLocN[], int imin[3], int remains[3]){
  for (int c = 0; c<3; c++){
    if (remains[c]){
      if (mygrid->rproc_imin[c][mygrid->rmyid[c]]>imin[c])
        itodo[c] = mygrid->rproc_imin[c][mygrid->rmyid[c]] - imin[c];
      else
        itodo[c] = 0;
      itodo[c + 3] = uniqueLocN[c];
    }
    else{
      itodo[c] = 0;
      itodo[c + 3] = 1;
    }
  }
}

void OUTPUT_MANAGER::writeSmallHeader(MPI_File thefile, int uniqueLocN[], int imin[], int remains[]){
  MPI_Status status;
  int itodo[6];
  prepareIntegerSmallHeader(itodo, uniqueLocN, imin, remains);
#ifndef DEBUG_NO_MPI_FILE_WRITE
#if defined (USE_MPI_FILE_WRITE_ALL)
  MPI_File_write_all(thefile, itodo, 6, MPI_INT, &status);
#else
  MPI_File_write(thefile, itodo, 6, MPI_INT, &status);
#endif
#endif
}

void OUTPUT_MANAGER::writeSmallHeaderSingleFile(std::string fileName, int uniqueLocN[], int imin[], int remains[]){
  int itodo[6];
  std::ofstream thefile;
  thefile.open(fileName.c_str(), std::ios::app);
  prepareIntegerSmallHeader(itodo, uniqueLocN, imin, remains);
  thefile.write((char*)itodo, 6 * sizeof(int));
  thefile.close();
}

void OUTPUT_MANAGER::prepareFloatField(float *todo, int NN[3], int origin[3], request req){
  int offset = 0, Ncomp = 3;
  if (req.type == OUT_E_FIELD)
    offset = 0;
  if (req.type == OUT_B_FIELD)
    offset = 3;

  int ii, jj, kk;
  int Nx, Ny, Nz;
  Nx = NN[0];
  Ny = NN[1];
  Nz = NN[2];

  if ((req.type == OUT_E_FIELD) || (req.type == OUT_B_FIELD)){
    Ncomp = 3;
    for (int k = 0; k < Nz; k++){
      kk = k + origin[2];
      for (int j = 0; j < Ny; j++){
        jj = j + origin[1];
        for (int i = 0; i < Nx; i++){
          ii = i + origin[0];
          for (int c = 0; c < Ncomp; c++)
            todo[c + i*Ncomp + j*Nx*Ncomp + k*Ny*Nx*Ncomp] = (float)myfield->VEB(c + offset, ii, jj, kk);
        }
      }
    }
  }
  else if (req.type == OUT_SPEC_DENSITY){
    Ncomp = 1;
    for (int k = 0; k < Nz; k++){
      kk = k + origin[2];
      for (int j = 0; j < Ny; j++){
        jj = j + origin[1];
        for (int i = 0; i < Nx; i++){
          ii = i + origin[0];
          todo[i + j*Nx + k*Ny*Nx] = (float)mycurrent->density(ii, jj, kk);
        }
      }
    }
  }
  else if (req.type == OUT_CURRENT){
    Ncomp = 3;
    for (int k = 0; k < Nz; k++){
      kk = k + origin[2];
      for (int j = 0; j < Ny; j++){
        jj = j + origin[1];
        for (int i = 0; i < Nx; i++){
          ii = i + origin[0];
          for (int c = 0; c < Ncomp; c++)
            todo[c + i*Ncomp + j*Nx + k*Ny*Nx] = (float)mycurrent->JJ(c, ii, jj, kk);
        }
      }
    }
  }

}

void OUTPUT_MANAGER::findDispForSetView(MPI_Offset *disp, int myOutputID, int *totUniquePoints, int big_header, int small_header, int Ncomp){
  MPI_Offset mydisp = big_header;
  for (int rank = 0; rank < myOutputID; rank++)
    mydisp += small_header + totUniquePoints[rank] * sizeof(float)*Ncomp;
  if (mydisp < 0){
    std::cout << "a problem occurred when trying to mpi_file_set_view in writeEMFieldBinary" << std::endl;
    std::cout << "myrank=" << mygrid->myid << " disp=" << mydisp << std::endl;
    exit(33);
  }
  disp[0] = mydisp;
}

void OUTPUT_MANAGER::setLocalOutputOffset(int *origin, int locimin[], int ri[], int remains[]){
  for (int c = 0; c < 3; c++){
    if (remains[c]){
      origin[c] = locimin[c];
    }
    else{
      origin[c] = ri[c];
    }
  }
}

void OUTPUT_MANAGER::writeCPUFieldValues(MPI_File thefile, int uniqueLocN[], int locimin[], int remains[3], request req){
  MPI_Status status;
  int Ncomp = 3;
  if ((req.type == OUT_E_FIELD) || (req.type == OUT_B_FIELD))
    Ncomp = 3;
  else if (req.type == OUT_SPEC_DENSITY)
    Ncomp = 1;
  else if (req.type == OUT_CURRENT)
    Ncomp = 3;

  int origin[3];
  int size = Ncomp*uniqueLocN[0] * uniqueLocN[1] * uniqueLocN[2];
  float *todo;
  todo = new float[size];
  int ri[3], globalri[3];
  nearestInt(myDomains[req.domain]->coordinates, ri, globalri);

  setLocalOutputOffset(origin, locimin, ri, remains);
  prepareFloatField(todo, uniqueLocN, origin, req);
#ifndef DEBUG_NO_MPI_FILE_WRITE
#if defined (USE_MPI_FILE_WRITE_ALL)
  MPI_File_write_all(thefile, todo, size, MPI_FLOAT, &status);
#else
  MPI_File_write(thefile, todo, size, MPI_FLOAT, &status);
#endif
#endif
  delete[]todo;
}

void OUTPUT_MANAGER::prepareCPUFieldValues(float *buffer, int uniqueLocN[], int imin[], int locimin[], int remains[3], request req){
  int Ncomp = 3;
  if ((req.type == OUT_E_FIELD) || (req.type == OUT_B_FIELD))
    Ncomp = 3;
  else if (req.type == OUT_SPEC_DENSITY)
    Ncomp = 1;
  else if (req.type == OUT_CURRENT)
    Ncomp = 3;

  int origin[3];
  int ri[3], globalri[3];
  nearestInt(myDomains[req.domain]->coordinates, ri, globalri);

  setLocalOutputOffset(origin, locimin, ri, remains);
  prepareIntegerSmallHeader((int*)buffer, uniqueLocN, imin, remains);
  prepareFloatField(buffer + 6, uniqueLocN, origin, req);
}

void OUTPUT_MANAGER::writeCPUParticlesValues(MPI_File thefile, double rmin[3], double rmax[3], SPECIE* spec){
  MPI_Status status;

  float *buf;
  int dimensione = 100000;
  buf = new float[dimensione*spec->Ncomp];
  int counter = 0;
  double rr[3];
  for (int p = 0; p < spec->Np; p++){
    rr[0] = spec->r0(p);
    rr[1] = spec->r1(p);
    rr[2] = spec->r2(p);
    if (rmax[0] >= rr[0] && rmin[0] < rr[0]){
      if (mygrid->getDimensionality() < 2 || (rmax[1] >= rr[1] && rmin[1] < rr[1])){
        if (mygrid->getDimensionality() < 3 || (rmax[2] >= rr[2] && rmin[2] < rr[2])){
          for (int c = 0; c < spec->Ncomp; c++){
            buf[c + counter*spec->Ncomp] = (float)spec->ru(c, p);
          }
          counter++;
        }
      }
    }
    if (counter == dimensione){
#ifndef DEBUG_NO_MPI_FILE_WRITE
      MPI_File_write(thefile, buf, counter*spec->Ncomp, MPI_FLOAT, &status);
#endif
      counter = 0;
    }
  }
  if (counter > 0){
#ifndef DEBUG_NO_MPI_FILE_WRITE
    MPI_File_write(thefile, buf, counter*spec->Ncomp, MPI_FLOAT, &status);
#endif
  }
  delete[]buf;

}

void OUTPUT_MANAGER::writeCPUParticlesValuesSingleFile(std::string  fileName, double rmin[3], double rmax[3], SPECIE* spec){
  std::ofstream thefile;
  thefile.open(fileName.c_str(), std::ios::app);
  float *buf;
  int dimensione = 100000;
  buf = new float[dimensione*spec->Ncomp];
  int counter = 0;
  double rr[3];
  for (int p = 0; p < spec->Np; p++){
    rr[0] = spec->r0(p);
    rr[1] = spec->r1(p);
    rr[2] = spec->r2(p);
    if (rmax[0] >= rr[0] && rmin[0] < rr[0]){
      if (mygrid->getDimensionality() < 2 || (rmax[1] >= rr[1] && rmin[1] < rr[1])){
        if (mygrid->getDimensionality() < 3 || (rmax[2] >= rr[2] && rmin[2] < rr[2])){
          for (int c = 0; c < spec->Ncomp; c++){
            buf[c + counter*spec->Ncomp] = (float)spec->ru(c, p);
          }
          counter++;
        }
      }
    }
    if (counter == dimensione){
      thefile.write((char*)buf, counter*spec->Ncomp*sizeof(float));
      counter = 0;
    }
  }
  if (counter > 0){
    thefile.write((char*)buf, counter*spec->Ncomp*sizeof(float));
  }
  delete[]buf;
  thefile.close();
}

void OUTPUT_MANAGER::writeCPUFieldValuesSingleFile(std::string  fileName, int uniqueLocN[], int locimin[], int remains[], request req){
  std::ofstream thefile;
  thefile.open(fileName.c_str(), std::ios::app);
  int Ncomp = 3;
  if ((req.type == OUT_E_FIELD) || (req.type == OUT_B_FIELD))
    Ncomp = 3;
  else if (req.type == OUT_SPEC_DENSITY)
    Ncomp = 1;
  else if (req.type == OUT_CURRENT)
    Ncomp = 3;

  int origin[3];
  int size = Ncomp*uniqueLocN[0] * uniqueLocN[1] * uniqueLocN[2];
  float *todo;
  todo = new float[size];
  int ri[3], globalri[3];
  nearestInt(myDomains[req.domain]->coordinates, ri, globalri);

  setLocalOutputOffset(origin, locimin, ri, remains);
  prepareFloatField(todo, uniqueLocN, origin, req);
  thefile.write((char*)todo, size*sizeof(float));
  delete[]todo;
  thefile.close();
}

void OUTPUT_MANAGER::writeGridFieldSubDomain(std::string fileName, request req){
  int Ncomp = 3;
  if ((req.type == OUT_E_FIELD) || (req.type == OUT_B_FIELD))
    Ncomp = 3;
  else if (req.type == OUT_SPEC_DENSITY)
    Ncomp = 1;
  else if (req.type == OUT_CURRENT)
    Ncomp = 3;

  int isInMyHyperplane = false, shouldIWrite = false;
  int uniqueN[3], uniqueLocN[3], slice_rNproc[3];
  int remains[3];
  int imin[3], imax[3], locimin[3], locimax[3];

  setAndCheckRemains(remains, myDomains[req.domain]->remainingCoord);
  findGlobalIntegerBoundaries(myDomains[req.domain]->rmin, myDomains[req.domain]->rmax, imin, imax);
  findLocalIntegerBoundaries(myDomains[req.domain]->rmin, myDomains[req.domain]->rmax, locimin, locimax);
  findNumberOfProcsWithinSubdomain(slice_rNproc, imin, imax, remains);

  findGlobalSubdomainUniquePointsNumber(uniqueN, imin, imax, remains);
  findLocalSubdomainUniquePointsNumber(uniqueLocN, locimin, locimax, remains);

  isInMyHyperplane = isThePointInMyDomain(myDomains[req.domain]->coordinates);

  MPI_Comm outputCommunicator;


  if(shouldICreateHyperplane(remains)){
    MPI_Comm sliceCommunicator;
    int mySliceID, sliceNProc;
    MPI_Cart_sub(mygrid->cart_comm, remains, &sliceCommunicator);
    MPI_Comm_rank(sliceCommunicator, &mySliceID);
    MPI_Comm_size(sliceCommunicator, &sliceNProc);
    MPI_Allreduce(MPI_IN_PLACE, &isInMyHyperplane, 1, MPI_INT, MPI_LOR, sliceCommunicator);

    if (isInMyHyperplane)
      shouldIWrite = amIInTheSubDomain(req);

    MPI_Comm_free(&sliceCommunicator);
    MPI_Comm_split(mygrid->cart_comm, shouldIWrite, 0, &outputCommunicator);
  }
  else{
    shouldIWrite = amIInTheSubDomain(req);
    MPI_Comm_split(mygrid->cart_comm, shouldIWrite, 0, &outputCommunicator);
  }
  int myOutputID, outputNProc;
  MPI_Comm_rank(outputCommunicator, &myOutputID);
  MPI_Comm_size(outputCommunicator, &outputNProc);


#ifdef FIELDS_USE_MPI_FILE_OUTPUT
  int *totUniquePoints = new int[outputNProc];
  totUniquePoints[myOutputID] = uniqueLocN[0] * uniqueLocN[1] * uniqueLocN[2];
  MPI_Allgather(MPI_IN_PLACE, 1, MPI_INT, totUniquePoints, 1, MPI_INT, outputCommunicator);

  MPI_Offset disp = 0;
  const int smallHeaderSize = (3 + 3) * sizeof(int);
  const int bigHeaderSize = (1 + 3 + 3 + 1)*sizeof(int) + (uniqueN[0] + uniqueN[1] + uniqueN[2])*sizeof(float);

  MPI_File thefile;
  char* nomefile = new char[fileName.length() + 1];
  strcpy(nomefile, fileName.c_str());

  if (shouldIWrite){
    MPI_File_open(outputCommunicator, nomefile, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &thefile);

    if (myOutputID == 0){
      MPI_File_set_view(thefile, 0, MPI_FLOAT, MPI_FLOAT, (char *) "native", MPI_INFO_NULL);
      writeBigHeader(thefile, uniqueN, imin, slice_rNproc, Ncomp);
    }
    else{
      findDispForSetView(&disp, myOutputID, totUniquePoints, bigHeaderSize, smallHeaderSize, Ncomp);
      MPI_File_set_view(thefile, disp, MPI_FLOAT, MPI_FLOAT, (char *) "native", MPI_INFO_NULL);
    }

    writeSmallHeader(thefile, uniqueLocN, imin, remains);
    writeCPUFieldValues(thefile, uniqueLocN, locimin, remains, req);

    MPI_File_close(&thefile);
  }
  delete[] nomefile;
#elif defined(FIELDS_USE_INDIVIDUAL_FILE_OUTPUT)
  //MPI_File thefile;
  char* nomefile = new char[fileName.length() + 1];
  strcpy(nomefile, fileName.c_str());
  if (shouldIWrite){
    std::stringstream myFileName;
    myFileName << fileName << "." << std::setfill('0') << std::setw(5) << myOutputID;
    if (myOutputID == 0){
      writeBigHeaderSingleFile(myFileName.str(), uniqueN, imin, slice_rNproc, Ncomp);
    }

    writeSmallHeaderSingleFile(myFileName.str(), uniqueLocN, imin, remains);
    writeCPUFieldValuesSingleFile(myFileName.str(), uniqueLocN, locimin, remains, req);

  }
#elif defined(FIELDS_USE_SEPARATE_FILES_MACROGROUPS)


  if(shouldIWrite){
    int fileCommunicatorID = myOutputID/multifileGroupSize;
    std::stringstream myFileName;
    myFileName << fileName << "." << std::setfill('0') << std::setw(5) << fileCommunicatorID;
    char *nomefile = new char[myFileName.str().size() + 1];
    strcpy(nomefile, myFileName.str().c_str());

    MPI_Comm FileCommunicator;
    MPI_Comm_split(outputCommunicator, fileCommunicatorID, 0, &FileCommunicator);
    int myFileId, fileNproc;
    MPI_Comm_size(FileCommunicator, &fileNproc);
    MPI_Comm_rank(FileCommunicator, &myFileId);

    int *totUniquePoints = new int[fileNproc];
    totUniquePoints[myFileId] = uniqueLocN[0] * uniqueLocN[1] * uniqueLocN[2];
    MPI_Allgather(MPI_IN_PLACE, 1, MPI_INT, totUniquePoints, 1, MPI_INT, FileCommunicator);

    MPI_Offset disp = 0;
    const int smallHeaderSize = (3 + 3) * sizeof(int);
    const int bigHeaderSize = (1 + 3 + 3 + 1)*sizeof(int) + (uniqueN[0] + uniqueN[1] + uniqueN[2])*sizeof(float);

    MPI_File thefile;

    MPI_Comm groupCommunicator;
    MPI_Comm_split(FileCommunicator, (myFileId/fieldGroupSize), 0, &groupCommunicator);
    int myGroupId, groupNproc;
    MPI_Comm_size(groupCommunicator, &groupNproc);
    MPI_Comm_rank(groupCommunicator, &myGroupId);

    MPI_Comm MPIFileCommunicator;
    MPI_Comm_split(FileCommunicator, (myFileId%fieldGroupSize), 0, &MPIFileCommunicator);
    int myMPIFileId, MPIFileNproc;
    MPI_Comm_size(MPIFileCommunicator, &MPIFileNproc);
    MPI_Comm_rank(MPIFileCommunicator, &myMPIFileId);

    int *groupBufferSize = new int[groupNproc];
    int maxBufferSize, myBufferSize = (totUniquePoints[myFileId] * Ncomp + 6);
    groupBufferSize[myGroupId] = myBufferSize;
    MPI_Allgather(MPI_IN_PLACE, 1, MPI_INT, groupBufferSize, 1, MPI_INT, groupCommunicator);
    MPI_Allreduce(&myBufferSize, &maxBufferSize, 1, MPI_INT, MPI_MAX, groupCommunicator);

    float *databuf=new float[maxBufferSize];

    int tag = 11;
    if(myGroupId !=0){
      prepareCPUFieldValues(databuf, uniqueLocN, imin, locimin, remains, req);
      MPI_Send( databuf, maxBufferSize, MPI_FLOAT, 0, tag,groupCommunicator);
    }
    else{
      MPI_Status status;
      MPI_File_open(MPIFileCommunicator, nomefile, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &thefile);

      //BEGIN OF NEW WRITE ALL CHANGE
//      int maxGroupNproc;
//      MPI_Allreduce(&groupNproc, &maxGroupNproc, 1, MPI_INT, MPI_MAX, MPIFileCommunicator);
//      int rest = maxGroupNproc - groupNproc;
//      //END
      if (myOutputID == 0){
        MPI_File_set_view(thefile, 0, MPI_FLOAT, MPI_FLOAT, (char *) "native", MPI_INFO_NULL);
        writeBigHeader(thefile, uniqueN, imin, slice_rNproc, Ncomp);
      }
      else{
        findDispForSetView(&disp, myFileId, totUniquePoints, bigHeaderSize, smallHeaderSize, Ncomp);
        if(fileCommunicatorID!=0)
          disp-=bigHeaderSize;
        MPI_File_set_view(thefile, disp, MPI_DOUBLE, MPI_DOUBLE, (char *) "native", MPI_INFO_NULL);
      }

      //la mia roba: la preparo e la scrivo
      prepareCPUFieldValues(databuf, uniqueLocN, imin, locimin, remains, req);
      MPI_File_write(thefile, databuf, groupBufferSize[0], MPI_FLOAT, &status);

      for (int procID = 1; procID < (groupNproc); procID++){
        MPI_Recv(databuf, maxBufferSize, MPI_FLOAT, procID, tag, groupCommunicator, &status);
        MPI_File_write(thefile, databuf, groupBufferSize[procID], MPI_FLOAT, &status);
      }
      //BEGIN OF NEW WRITE ALL CHANGE
//      for(int i=0; i< rest; i++){
//        MPI_File_write_all(thefile, databuf, 0, MPI_FLOAT, &status);
//      }
      //END
      MPI_File_close(&thefile);
    }
    delete[] nomefile;
    delete[] databuf;
    delete[] totUniquePoints;
    MPI_Comm_free(&groupCommunicator);
    MPI_Comm_free(&MPIFileCommunicator);
    MPI_Comm_free(&FileCommunicator);

  }


#elif defined(FIELDS_USE_MULTI_FILE)
  int *totUniquePoints = new int[outputNProc];
  totUniquePoints[myOutputID] = uniqueLocN[0] * uniqueLocN[1] * uniqueLocN[2];
  MPI_Allgather(MPI_IN_PLACE, 1, MPI_INT, totUniquePoints, 1, MPI_INT, outputCommunicator);

  if (shouldIWrite){
    MPI_Comm groupCommunicator;
    MPI_Comm_split(outputCommunicator, (myOutputID/fieldGroupSize), 0, &groupCommunicator);

    int myGroupId, groupNproc;
    MPI_Comm_size(groupCommunicator, &groupNproc);
    MPI_Comm_rank(groupCommunicator, &myGroupId);

    int *groupBufferSize = new int[groupNproc];
    int maxBufferSize, myBufferSize = (totUniquePoints[myOutputID] * Ncomp + 6);
    groupBufferSize[myGroupId] = myBufferSize;
    MPI_Allgather(MPI_IN_PLACE, 1, MPI_INT, groupBufferSize, 1, MPI_INT, groupCommunicator);
    MPI_Allreduce(&myBufferSize, &maxBufferSize, 1, MPI_INT, MPI_MAX, groupCommunicator);

    float **databuf=new float*[2];
    for(int i=0;i<2;i++){
      databuf[i]=new float[maxBufferSize];
    }

    int tag[]={11,12};
    if(myGroupId !=0){
      MPI_Status status;
      MPI_Request request[2];
      prepareCPUFieldValues(databuf[0], uniqueLocN, imin, locimin, remains, req);
      MPI_Isend(databuf[0], maxBufferSize, MPI_FLOAT, 0, tag[myGroupId%2],
          groupCommunicator, &request[myGroupId%2]);
      MPI_Wait(&request[myGroupId%2], &status);
    }
    else{
      std::stringstream myFileName;
      myFileName << fileName << "." << std::setfill('0') << std::setw(5) << myOutputID;
      writeBigHeaderSingleFile(myFileName.str(), uniqueN, imin, slice_rNproc, Ncomp);

      std::ofstream thefile;
      thefile.open(myFileName.str().c_str(), std::ios::app);

      prepareCPUFieldValues(databuf[0], uniqueLocN, imin, locimin, remains, req);
      thefile.write((char*)databuf[0], groupBufferSize[0]*sizeof(float));

      if(groupNproc==2){
        MPI_Status status;
        MPI_Request request[2];
        int procID=1;
        MPI_Irecv(databuf[procID%2], maxBufferSize, MPI_FLOAT, procID, tag[procID%2], groupCommunicator, &request[procID%2]);
        MPI_Wait(&request[procID%2], &status);
        thefile.write((char*)databuf[procID%2], groupBufferSize[procID%2]*sizeof(float));
      }

      if(groupNproc>2){
        int procID;
        MPI_Status status;
        MPI_Request request[2];
        procID = 1;
        MPI_Irecv(databuf[1], maxBufferSize, MPI_FLOAT, procID, tag[procID%2], groupCommunicator, &request[procID%2]);
#ifndef _ASYNCHRONOUS
        MPI_Wait(&request[procID%2], &status);
#endif
        procID++;
        MPI_Irecv(databuf[0], maxBufferSize, MPI_FLOAT, procID, tag[procID%2], groupCommunicator, &request[procID%2]);
#ifndef _ASYNCHRONOUS
        MPI_Wait(&request[procID%2], &status);
#endif

        for (procID = 1; procID < (groupNproc-2); procID++){
          MPI_Wait(&request[procID%2], &status);
          thefile.write((char*)databuf[procID%2], groupBufferSize[procID%2]*sizeof(float));
          MPI_Irecv(databuf[procID%2], maxBufferSize, MPI_FLOAT, procID+2, tag[procID%2], groupCommunicator, &request[procID%2]);
#ifndef _ASYNCHRONOUS
          MPI_Wait(&request[procID%2], &status);
#endif

        }

        MPI_Wait(&request[procID%2], &status);
        thefile.write((char*)databuf[procID%2], groupBufferSize[procID%2]*sizeof(float));
        procID++;
        MPI_Wait(&request[(procID)%2], &status);
        thefile.write((char*)databuf[procID%2], groupBufferSize[procID%2]*sizeof(float));
      }
      thefile.close();
    }
    delete[] databuf[0];
    delete[] databuf[1];

    MPI_Comm_free(&groupCommunicator);
  }
  delete[] totUniquePoints;
#else
  int *totUniquePoints = new int[outputNProc];
  totUniquePoints[myOutputID] = uniqueLocN[0] * uniqueLocN[1] * uniqueLocN[2];
  MPI_Allgather(MPI_IN_PLACE, 1, MPI_INT, totUniquePoints, 1, MPI_INT, outputCommunicator);

  MPI_Offset disp = 0;
  const int smallHeaderSize = (3 + 3) * sizeof(int);
  const int bigHeaderSize = (1 + 3 + 3 + 1)*sizeof(int) + (uniqueN[0] + uniqueN[1] + uniqueN[2])*sizeof(float);

  MPI_File thefile;
  char* nomefile = new char[fileName.length() + 1];
  strcpy(nomefile, fileName.c_str());
  if (shouldIWrite){
    MPI_Comm groupCommunicator;
    MPI_Comm MPIFileCommunicator;
    MPI_Comm_split(outputCommunicator, (myOutputID/fieldGroupSize), 0, &groupCommunicator);
    MPI_Comm_split(outputCommunicator, (myOutputID%fieldGroupSize), 0, &MPIFileCommunicator);

    int myGroupId, groupNproc;
    MPI_Comm_size(groupCommunicator, &groupNproc);
    MPI_Comm_rank(groupCommunicator, &myGroupId);
    int myMPIFileId, MPIFileNproc;
    MPI_Comm_size(MPIFileCommunicator, &MPIFileNproc);
    MPI_Comm_rank(MPIFileCommunicator, &myMPIFileId);

    int *groupBufferSize = new int[groupNproc];
    int maxBufferSize, myBufferSize = (totUniquePoints[myOutputID] * Ncomp + 6);
    groupBufferSize[myGroupId] = myBufferSize;
    MPI_Allgather(MPI_IN_PLACE, 1, MPI_INT, groupBufferSize, 1, MPI_INT, groupCommunicator);
    MPI_Allreduce(&myBufferSize, &maxBufferSize, 1, MPI_INT, MPI_MAX, groupCommunicator);

    float **databuf=new float*[2];
    for(int i=0;i<2;i++){
      databuf[i]=new float[maxBufferSize];
    }

    int tag[]={11,12};
    if(myGroupId !=0){
      MPI_Status status;
      MPI_Request request[2];
      prepareCPUFieldValues(databuf[0], uniqueLocN, imin, locimin, remains, req);
      MPI_Isend(databuf[0], maxBufferSize, MPI_FLOAT, 0, tag[myGroupId%2],
          groupCommunicator, &request[myGroupId%2]);
      MPI_Wait(&request[myGroupId%2], &status);
    }
    else{
      MPI_Status status;
      MPI_File_open(MPIFileCommunicator, nomefile, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &thefile);

      if (myOutputID == 0){
        MPI_File_set_view(thefile, 0, MPI_FLOAT, MPI_FLOAT, (char *) "native", MPI_INFO_NULL);
        writeBigHeader(thefile, uniqueN, imin, slice_rNproc, Ncomp);
      }
      else{
        findDispForSetView(&disp, myOutputID, totUniquePoints, bigHeaderSize, smallHeaderSize, Ncomp);
        MPI_File_set_view(thefile, disp, MPI_DOUBLE, MPI_DOUBLE, (char *) "native", MPI_INFO_NULL);
      }

      //la mia roba: la preparo e la scrivo
      prepareCPUFieldValues(databuf[0], uniqueLocN, imin, locimin, remains, req);
      MPI_File_write(thefile, databuf[0], groupBufferSize[0], MPI_FLOAT, &status);

      if(groupNproc==2){
        MPI_Status status;
        MPI_Request request[2];
        int procID=1;
        MPI_Irecv(databuf[procID%2], maxBufferSize, MPI_FLOAT, procID, tag[procID%2], groupCommunicator, &request[procID%2]);
        MPI_Wait(&request[procID%2], &status);
        MPI_File_write(thefile, databuf[procID%2], groupBufferSize[procID], MPI_FLOAT, &status);
      }
      if(groupNproc>2){
        int procID;
        MPI_Status status;
        MPI_Request request[2];
        procID = 1;
        MPI_Irecv(databuf[1], maxBufferSize, MPI_FLOAT, procID, tag[procID%2], groupCommunicator, &request[procID%2]);
#ifndef _ASYNCHRONOUS
        MPI_Wait(&request[procID%2], &status);
#endif
        procID++;
        MPI_Irecv(databuf[0], maxBufferSize, MPI_FLOAT, procID, tag[procID%2], groupCommunicator, &request[procID%2]);
#ifndef _ASYNCHRONOUS
        MPI_Wait(&request[procID%2], &status);
#endif

        for (procID = 1; procID < (groupNproc-2); procID++){
          MPI_Wait(&request[procID%2], &status);
          MPI_File_write(thefile, databuf[procID%2], groupBufferSize[procID], MPI_FLOAT, &status);
          MPI_Irecv(databuf[procID%2], maxBufferSize, MPI_FLOAT, procID+2, tag[procID%2], groupCommunicator, &request[procID%2]);
#ifndef _ASYNCHRONOUS
          MPI_Wait(&request[procID%2], &status);
#endif
        }

        MPI_Wait(&request[procID%2], &status);
        MPI_File_write(thefile, databuf[procID%2], groupBufferSize[procID], MPI_FLOAT, &status);
        procID++;
        MPI_Wait(&request[(procID)%2], &status);
        MPI_File_write(thefile, databuf[(procID)%2], groupBufferSize[procID], MPI_FLOAT, &status);
      }
      MPI_File_close(&thefile);
    }
    delete[] nomefile;
    delete[] databuf[0];
    delete[] databuf[1];

    MPI_Comm_free(&groupCommunicator);
    MPI_Comm_free(&MPIFileCommunicator);
  }
  delete[] totUniquePoints;
#endif

  MPI_Comm_free(&outputCommunicator);

}


void OUTPUT_MANAGER::callEMFieldDomain(request req){
  std::stringstream groupname;
  groupname << "Group" << fieldGroupSize;
  if (req.type == OUT_E_FIELD){
        std::string nameBin = composeOutputName(outputDir, "E_FIELD", myDomains[req.domain]->name, "", req.domain, req.dtime, ".bin");
    //std::string nameBin = composeOutputName(outputDir, "E_FIELD", myDomains[req.domain]->name, groupname.str(), req.domain, req.dtime, ".bin");
    writeGridFieldSubDomain(nameBin, req);
  }
  else if (req.type == OUT_B_FIELD){
        std::string nameBin = composeOutputName(outputDir, "B_FIELD", myDomains[req.domain]->name, "", req.domain, req.dtime, ".bin");
    //std::string nameBin = composeOutputName(outputDir, "B_FIELD", myDomains[req.domain]->name, groupname.str(), req.domain, req.dtime, ".bin");
    writeGridFieldSubDomain(nameBin, req);
  }

}

void OUTPUT_MANAGER::interpolateEBFieldsToPosition(double position[3], double E[3], double B[3]){
  int hii[3], wii[3];
  double hiw[3][3], wiw[3][3];
  double rr, rh, rr2, rh2;
  int i1, j1, k1, i2, j2, k2;
  double dvol;
  double mycsi[3];

  for (int c = 0; c < 3; c++){
    hiw[c][1] = wiw[c][1] = 1;
    hii[c] = wii[c] = 0;
  }
  if (mygrid->isStretched()){
    for (int c = 0; c < mygrid->getDimensionality(); c++){
      mycsi[c] = mygrid->unStretchGrid(position[c], c);
      rr = mygrid->dri[c] * (mycsi[c] - mygrid->csiminloc[c]);

      rh = rr - 0.5;
      wii[c] = (int)floor(rr + 0.5); //whole integer int
      hii[c] = (int)floor(rr);     //half integer int
      rr -= wii[c];
      rh -= hii[c];
      rr2 = rr*rr;
      rh2 = rh*rh;

      wiw[c][1] = 0.75 - rr2;
      wiw[c][2] = 0.5*(0.25 + rr2 + rr);
      wiw[c][0] = 1. - wiw[c][1] - wiw[c][2];

      hiw[c][1] = 0.75 - rh2;
      hiw[c][2] = 0.5*(0.25 + rh2 + rh);
      hiw[c][0] = 1. - hiw[c][1] - hiw[c][2];
    }
  }
  else{
    for (int c = 0; c < mygrid->getDimensionality(); c++){
      rr = mygrid->dri[c] * (position[c] - mygrid->rminloc[c]);
      rh = rr - 0.5;
      wii[c] = (int)floor(rr + 0.5); //whole integer int
      hii[c] = (int)floor(rr);     //half integer int
      rr -= wii[c];
      rh -= hii[c];
      rr2 = rr*rr;
      rh2 = rh*rh;

      wiw[c][1] = 0.75 - rr2;
      wiw[c][2] = 0.5*(0.25 + rr2 + rr);
      wiw[c][0] = 1. - wiw[c][1] - wiw[c][2];

      hiw[c][1] = 0.75 - rh2;
      hiw[c][2] = 0.5*(0.25 + rh2 + rh);
      hiw[c][0] = 1. - hiw[c][1] - hiw[c][2];
    }
  }

  E[0] = E[1] = E[2] = B[0] = B[1] = B[2] = 0;

  switch (mygrid->getDimensionality())
  {
    case 3:
      for (int k = 0; k < 3; k++)
      {
        k1 = k + wii[2] - 1;
        k2 = k + hii[2] - 1;
        for (int j = 0; j < 3; j++)
        {
          j1 = j + wii[1] - 1;
          j2 = j + hii[1] - 1;
          for (int i = 0; i < 3; i++)
          {
            i1 = i + wii[0] - 1;
            i2 = i + hii[0] - 1;
            dvol = hiw[0][i] * wiw[1][j] * wiw[2][k],
                E[0] += myfield->E0(i2, j1, k1)*dvol;  //Ex
            dvol = wiw[0][i] * hiw[1][j] * wiw[2][k],
                E[1] += myfield->E1(i1, j2, k1)*dvol;  //Ey
            dvol = wiw[0][i] * wiw[1][j] * hiw[2][k],
                E[2] += myfield->E2(i1, j1, k2)*dvol;  //Ez

            dvol = wiw[0][i] * hiw[1][j] * hiw[2][k],
                B[0] += myfield->B0(i1, j2, k2)*dvol;  //Bx
            dvol = hiw[0][i] * wiw[1][j] * hiw[2][k],
                B[1] += myfield->B1(i2, j1, k2)*dvol;  //By
            dvol = hiw[0][i] * hiw[1][j] * wiw[2][k],
                B[2] += myfield->B2(i2, j2, k1)*dvol;  //Bz
          }
        }
      }
      break;

    case 2:
      k1 = k2 = 0;
      for (int j = 0; j < 3; j++)
      {
        j1 = j + wii[1] - 1;
        j2 = j + hii[1] - 1;
        for (int i = 0; i < 3; i++)
        {
          i1 = i + wii[0] - 1;
          i2 = i + hii[0] - 1;
          dvol = hiw[0][i] * wiw[1][j],
              E[0] += myfield->E0(i2, j1, k1)*dvol;  //Ex
          dvol = wiw[0][i] * hiw[1][j],
              E[1] += myfield->E1(i1, j2, k1)*dvol;  //Ey
          dvol = wiw[0][i] * wiw[1][j],
              E[2] += myfield->E2(i1, j1, k2)*dvol;  //Ez

          dvol = wiw[0][i] * hiw[1][j],
              B[0] += myfield->B0(i1, j2, k2)*dvol;  //Bx
          dvol = hiw[0][i] * wiw[1][j],
              B[1] += myfield->B1(i2, j1, k2)*dvol;  //By
          dvol = hiw[0][i] * hiw[1][j],
              B[2] += myfield->B2(i2, j2, k1)*dvol;  //Bz
        }
      }
      break;

    case 1:
      k1 = k2 = j1 = j2 = 0;
      for (int i = 0; i < 3; i++)
      {
        i1 = i + wii[0] - 1;
        i2 = i + hii[0] - 1;
        dvol = hiw[0][i],
            E[0] += myfield->E0(i2, j1, k1)*dvol;  //Ex
        dvol = wiw[0][i],
            E[1] += myfield->E1(i1, j2, k1)*dvol;  //Ey
        dvol = wiw[0][i],
            E[2] += myfield->E2(i1, j1, k2)*dvol;  //Ez

        dvol = wiw[0][i],
            B[0] += myfield->B0(i1, j2, k2)*dvol;  //Bx
        dvol = hiw[0][i],
            B[1] += myfield->B1(i2, j1, k2)*dvol;  //By
        dvol = hiw[0][i],
            B[2] += myfield->B2(i2, j2, k1)*dvol;  //Bz
      }
      break;
  }

}

void OUTPUT_MANAGER::callEMFieldProbe(request req){
  double rr[3], EE[3], BB[3];
  rr[0] = myEMProbes[req.domain]->coordinates[0];
  rr[1] = myEMProbes[req.domain]->coordinates[1];
  rr[2] = myEMProbes[req.domain]->coordinates[2];

  if (rr[0] >= mygrid->rminloc[0] && rr[0] < mygrid->rmaxloc[0]){
    if (mygrid->getDimensionality() < 2 || (rr[1] >= mygrid->rminloc[1] && rr[1] < mygrid->rmaxloc[1])){
      if (mygrid->getDimensionality() < 3 || (rr[2] >= mygrid->rminloc[2] && rr[2] < mygrid->rmaxloc[2])){
        interpolateEBFieldsToPosition(rr, EE, BB);
        std::ofstream of0;
        of0.open(myEMProbes[req.domain]->fileName.c_str(), std::ios::app);
        of0 << " " << std::setw(diagNarrowWidth) << req.itime << " " << std::setw(diagWidth) << req.dtime;
        of0 << " " << std::setw(diagWidth) << EE[0] << " " << std::setw(diagWidth) << EE[1] << " " << std::setw(diagWidth) << EE[2];
        of0 << " " << std::setw(diagWidth) << BB[0] << " " << std::setw(diagWidth) << BB[1] << " " << std::setw(diagWidth) << BB[2];
        of0 << "\n";
        of0.close();
      }
    }
  }
}

void OUTPUT_MANAGER::callSpecDensity(request req){
  mycurrent->eraseDensity();
  myspecies[req.target]->density_deposition_standard(mycurrent);
  mycurrent->pbc();

  std::string nameBin = composeOutputName(outputDir, "DENS", myspecies[req.target]->name, myDomains[req.domain]->name, req.domain, req.dtime, ".bin");
  writeGridFieldSubDomain(nameBin, req);
}

void OUTPUT_MANAGER::writeCurrent(std::string fileName, request req){
  int Ncomp = 3;//myfield->getNcomp();
  int *totUniquePoints;
  int shouldIWrite = false;
  int uniqueN[3], slice_rNproc[3];
  double rr[3] = { myDomains[req.domain]->coordinates[0], myDomains[req.domain]->coordinates[1], myDomains[req.domain]->coordinates[2] };
  int ri[3];
  int remains[3] = { myDomains[req.domain]->remainingCoord[0], myDomains[req.domain]->remainingCoord[1], myDomains[req.domain]->remainingCoord[2] };

  if (myDomains[req.domain]->overrideFlag){
    rr[0] = 0.5*(mygrid->rmin[0] + mygrid->rmax[0]);
    rr[1] = 0.5*(mygrid->rmin[1] + mygrid->rmax[1]);
    rr[2] = 0.5*(mygrid->rmin[2] + mygrid->rmax[2]);
    remains[0] = 1;
    remains[1] = 1;
    remains[2] = 1;
  }


  for (int c = 0; c < 3; c++){
    if (remains[c]){
      uniqueN[c] = mygrid->uniquePoints[c];
      slice_rNproc[c] = mygrid->rnproc[c];
    }
    else{
      uniqueN[c] = 1;
      slice_rNproc[c] = 1;
    }
  }

  shouldIWrite = isThePointInMyDomain(rr);

  MPI_Comm sliceCommunicator;
  int mySliceID, sliceNProc;

  int dimension;
  MPI_Cart_sub(mygrid->cart_comm, remains, &sliceCommunicator);
  MPI_Comm_rank(sliceCommunicator, &mySliceID);
  MPI_Comm_size(sliceCommunicator, &sliceNProc);
  MPI_Cartdim_get(sliceCommunicator, &dimension);
  MPI_Allreduce(MPI_IN_PLACE, &shouldIWrite, 1, MPI_INT, MPI_LOR, sliceCommunicator);

  totUniquePoints = new int[sliceNProc];
  for (int rank = 0; rank < sliceNProc; rank++){
    int rid[3], idbookmark = 0;
    MPI_Cart_coords(sliceCommunicator, rank, dimension, rid);
    totUniquePoints[rank] = 1;
    for (int c = 0; c < 3; c++){
      if (remains[c]){
        totUniquePoints[rank] *= mygrid->rproc_NuniquePointsloc[c][rid[idbookmark]];
        idbookmark++;
      }
    }
  }

  MPI_Offset disp = 0;
  int small_header = 6 * sizeof(int);
  int big_header = (1 + 3 + 3 + 1)*sizeof(int)
      + (uniqueN[0] + uniqueN[1] + uniqueN[2])*sizeof(float);

  MPI_File thefile;
  MPI_Status status;

  char* nomefile = new char[fileName.size() + 1];

  nomefile[fileName.size()] = 0;
  sprintf(nomefile, "%s", fileName.c_str());

  if (shouldIWrite){
    MPI_File_open(sliceCommunicator, nomefile, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &thefile);
    int globalri[3];
    nearestInt(rr, ri, globalri);
    //+++++++++++ FILE HEADER  +++++++++++++++++++++
    if (mySliceID == 0){
      MPI_File_set_view(thefile, 0, MPI_FLOAT, MPI_FLOAT, (char *) "native", MPI_INFO_NULL);
      int itodo[8];
      itodo[0] = is_big_endian();
      itodo[1] = uniqueN[0];
      itodo[2] = uniqueN[1];
      itodo[3] = uniqueN[2];
      itodo[4] = slice_rNproc[0];
      itodo[5] = slice_rNproc[1];
      itodo[6] = slice_rNproc[2];
      itodo[7] = Ncomp;
#ifndef DEBUG_NO_MPI_FILE_WRITE
      MPI_File_write(thefile, itodo, 8, MPI_INT, &status);
#endif

      float *fcir[3];
      for (int c = 0; c < 3; c++){
        fcir[c] = new float[uniqueN[c]];
        for (int m = 0; m < uniqueN[c]; m++){
          fcir[c][m] = (float)mygrid->cir[c][m];
        }
      }
      for (int c = 0; c < 3; c++){
        if (remains[c]){
#ifndef DEBUG_NO_MPI_FILE_WRITE
          MPI_File_write(thefile, fcir[c], uniqueN[c], MPI_FLOAT, &status);
#endif
        }
        else{
#ifndef DEBUG_NO_MPI_FILE_WRITE
          MPI_File_write(thefile, &fcir[c][globalri[c]], 1, MPI_FLOAT, &status);
#endif
        }
      }
      for (int c = 0; c < 3; c++){
        delete[] fcir[c];
      }

    }
    //*********** END HEADER *****************

    disp = big_header;
    for (int rank = 0; rank < mySliceID; rank++)
      disp += small_header + totUniquePoints[rank] * sizeof(float)*Ncomp;
    if (disp < 0){
      std::cout << "a problem occurred when trying to mpi_file_set_view in writeEMFieldBinary" << std::endl;
      std::cout << "myrank=" << mygrid->myid << " disp=" << disp << std::endl;
      exit(33);
    }
    if (mySliceID != 0){
      MPI_File_set_view(thefile, disp, MPI_FLOAT, MPI_FLOAT, (char *) "native", MPI_INFO_NULL);
    }

    //+++++++++++ Start CPU HEADER  +++++++++++++++++++++
    {
      int itodo[6];
      for (int c = 0; c < 3; c++){
        if (remains[c]){
          itodo[c] = mygrid->rproc_imin[c][mygrid->rmyid[c]];
          itodo[c + 3] = mygrid->uniquePointsloc[c];
        }
        else{
          itodo[c] = 0;
          itodo[c + 3] = 1;
        }
      }
#ifndef DEBUG_NO_MPI_FILE_WRITE
      MPI_File_write(thefile, itodo, 6, MPI_INT, &status);
#endif
    }
    //+++++++++++ Start CPU Field Values  +++++++++++++++++++++
    {
      float *todo;
      int NN[3], Nx, Ny, Nz, origin[3];
      for (int c = 0; c < 3; c++){
        if (remains[c]){
          NN[c] = mygrid->uniquePointsloc[c];
          origin[c] = 0;
        }
        else{
          NN[c] = 1;
          origin[c] = ri[c];
        }
      }
      Nx = NN[0];
      Ny = NN[1];
      Nz = NN[2];
      int size = Ncomp*NN[0] * NN[1] * NN[2];
      todo = new float[size];
      int ii, jj, kk;
      for (int k = 0; k < Nz; k++){
        kk = k + origin[2];
        for (int j = 0; j < Ny; j++){
          jj = j + origin[1];
          for (int i = 0; i < Nx; i++){
            ii = i + origin[0];
            for (int c = 0; c < Ncomp; c++)
              todo[c + i*Ncomp + j*Nx*Ncomp + k*Ny*Nx*Ncomp] = (float)mycurrent->JJ(c, i, j, k);
          }
        }
      }
#ifndef DEBUG_NO_MPI_FILE_WRITE
      MPI_File_write(thefile, todo, size, MPI_FLOAT, &status);
#endif
      delete[]todo;
    }
    MPI_File_close(&thefile);
    delete[] nomefile;
    delete[] totUniquePoints;
  }
  MPI_Comm_free(&sliceCommunicator);
}

void  OUTPUT_MANAGER::callCurrent(request req){

  std::string nameBin = composeOutputName(outputDir, "J", "", myDomains[req.domain]->name, req.domain, req.dtime, ".bin");
  writeGridFieldSubDomain(nameBin, req);
  //writeCurrent(nameBin, req);


}

void OUTPUT_MANAGER::writeCPUParticlesValues(MPI_File thefile, SPECIE* spec, bool flagMarker){

  int outputNComp, NCompFloat;
  if (flagMarker){
    NCompFloat = spec->Ncomp - 1;
    outputNComp = spec->Ncomp + 1;

  }
  else{
    NCompFloat = spec->Ncomp;
    outputNComp = spec->Ncomp;
  }
  MPI_Status status;
  int dimensione, passaggi, resto;

  float *buf;
  dimensione = 100000;
  buf = new float[dimensione*outputNComp];
  passaggi = spec->Np / dimensione;
  resto = spec->Np % dimensione;

  for (int i = 0; i < passaggi; i++){
    for (int p = 0; p < dimensione; p++){
      int c;
      for (c = 0; c < NCompFloat; c++){
        buf[c + p*outputNComp] = (float)spec->ru(c, p + dimensione*i);
      }
      if (flagMarker){
        buf[c + p*outputNComp] = *((float*)(&(spec->marker(p + dimensione*i))));
        buf[c + 1 + p*outputNComp] = *((float*)(&(spec->marker(p + dimensione*i))) + 1);
      }
    }
#ifndef DEBUG_NO_MPI_FILE_WRITE
    MPI_File_write(thefile, buf, dimensione*outputNComp, MPI_FLOAT, &status);
#endif
  }
  for (int p = 0; p < resto; p++){
    int c;
    for (c = 0; c < NCompFloat; c++){
      buf[c + p*outputNComp] = (float)spec->ru(c, p + dimensione*passaggi);
    }
    if (flagMarker){
      buf[c + p*outputNComp] = *((float*)(&(spec->marker(p + dimensione*passaggi))));
      buf[c + 1 + p*outputNComp] = *((float*)(&(spec->marker(p + dimensione*passaggi))) + 1);
    }
#ifndef DEBUG_NO_MPI_FILE_WRITE
    MPI_File_write(thefile, buf, resto*outputNComp, MPI_FLOAT, &status);
#endif
    delete[]buf;
  }
}

void OUTPUT_MANAGER::writeCPUParticlesValues(MPI_File thefile, SPECIE* spec){

  MPI_Status status;
  int dimensione, passaggi, resto;
  int Ncomp = spec->Ncomp;
  float *buf;
  dimensione = particleBufferSize;
  buf = new float[dimensione*Ncomp];
  passaggi = spec->Np / dimensione;
  resto = spec->Np % dimensione;

  for (int i = 0; i < passaggi; i++){
    for (int p = 0; p < dimensione; p++){
      int c;
      for (c = 0; c < Ncomp; c++){
        buf[c + p*Ncomp] = (float)spec->ru(c, p + dimensione*i);
      }
    }
#ifndef DEBUG_NO_MPI_FILE_WRITE
    MPI_File_write(thefile, buf, dimensione*Ncomp, MPI_FLOAT, &status);
#endif
  }
  for (int p = 0; p < resto; p++){
    int c;
    for (c = 0; c < Ncomp; c++){
      buf[c + p*Ncomp] = (float)spec->ru(c, p + dimensione*passaggi);
    }
  }
#ifndef DEBUG_NO_MPI_FILE_WRITE
  MPI_File_write(thefile, buf, resto*Ncomp, MPI_FLOAT, &status);
#endif
  delete[]buf;
}

void OUTPUT_MANAGER::writeAllCPUParticlesValues(MPI_File thefile, SPECIE* spec, int maxNfloatLoc){
  MPI_Status status;
  int dimensione, passaggi, maxPassaggi, resto;
  int Ncomp = spec->Ncomp;
  float *buf;
  dimensione = particleBufferSize;
  buf = new float[dimensione*Ncomp];
  passaggi = spec->Np / dimensione;
  maxPassaggi = (maxNfloatLoc / Ncomp) / dimensione;
  resto = spec->Np % dimensione;

  for (int i = 0; i < passaggi; i++){
    for (int p = 0; p < dimensione; p++){
      int c;
      for (c = 0; c < Ncomp; c++){
        buf[c + p*Ncomp] = (float)spec->ru(c, p + dimensione*i);
      }
    }
#ifndef DEBUG_NO_MPI_FILE_WRITE
    MPI_File_write_all(thefile, buf, dimensione*Ncomp, MPI_FLOAT, &status);
#endif
  }
  for (int p = 0; p < resto; p++){
    int c;
    for (c = 0; c < Ncomp; c++){
      buf[c + p*Ncomp] = (float)spec->ru(c, p + dimensione*passaggi);
    }
  }
#ifndef DEBUG_NO_MPI_FILE_WRITE
  MPI_File_write_all(thefile, buf, resto*Ncomp, MPI_FLOAT, &status);

  for (int i = 0; i < (maxPassaggi - passaggi); i++){
    MPI_File_write_all(thefile, buf, 0, MPI_FLOAT, &status);
  }
  #endif
  delete[]buf;
}

void OUTPUT_MANAGER::writeCPUParticlesValuesSingleFile(std::string  fileName, SPECIE* spec){
  std::ofstream thefile;
  thefile.open(fileName.c_str(), std::ios::app);
  int dimensione, passaggi, resto;
  int Ncomp = spec->Ncomp;
  float *buf;
  dimensione = 100000;
  buf = new float[dimensione*Ncomp];
  passaggi = spec->Np / dimensione;
  resto = spec->Np % dimensione;

  for (int i = 0; i < passaggi; i++){
    for (int p = 0; p < dimensione; p++){
      int c;
      for (c = 0; c < Ncomp; c++){
        buf[c + p*Ncomp] = (float)spec->ru(c, p + dimensione*i);
      }
    }
    thefile.write((char*)buf, dimensione*Ncomp*sizeof(float));
  }
  for (int p = 0; p < resto; p++){
    int c;
    for (c = 0; c < Ncomp; c++){
      buf[c + p*Ncomp] = (float)spec->ru(c, p + dimensione*passaggi);
    }
  }
  thefile.write((char*)buf, resto*Ncomp*sizeof(float));
  delete[]buf;
  thefile.close();
}

int OUTPUT_MANAGER::packageSize(int bufsize, int* groupProcNumData, int procID, int packageNumber){
  if (packageNumber < groupProcNumData[procID]/bufsize){
    return bufsize;
  }
  else if (packageNumber == (groupProcNumData[procID]/bufsize)){
    return groupProcNumData[procID]%bufsize;
  }
  else{
    return 0;
  }
}

int OUTPUT_MANAGER::numPackages(int bufsize, int* groupProcNumData, int procID){
  return groupProcNumData[procID]/bufsize + 1;
}

void OUTPUT_MANAGER::fillRequestList(int bufsize, int* groupProcNumData, int groupNproc, std::vector<reqOutput> &reqList){
  for(int i = 1; i < groupNproc; i++){
    for(int p = 0; p < numPackages(bufsize, groupProcNumData,i);p++){
      reqOutput newReq;
      newReq.task = i;
      newReq.p = p;
      newReq.packageSize = packageSize(bufsize, groupProcNumData, i, p);
      reqList.push_back(newReq);
    }
  }
  std::sort(reqList.begin(),reqList.end(),compOutput);
}

void OUTPUT_MANAGER::writeAllSeparateFilesParticlesValues(std::string fileName, SPECIE* spec){

  MPI_File thefile;
  std::stringstream myFileName;
  int fileCommunicatorID = mygrid->myid/multifileGroupSize;
  myFileName << fileName << "." << std::setfill('0') << std::setw(5) << fileCommunicatorID;
  char *nomefile = new char[myFileName.str().size() + 1];
  nomefile[myFileName.str().size()] = 0;
  sprintf(nomefile, "%s", myFileName.str().c_str());

  MPI_Comm FileCommunicator;
  MPI_Comm_split(MPI_COMM_WORLD, fileCommunicatorID, 0, &FileCommunicator);
  int fileMyid, fileNproc;
  MPI_Comm_size(FileCommunicator, &fileNproc);
  MPI_Comm_rank(FileCommunicator, &fileMyid);

  int* NfloatLoc = new int[fileNproc];
  NfloatLoc[fileMyid] = spec->Np*spec->Ncomp;
  MPI_Allgather(MPI_IN_PLACE, 1, MPI_INT, NfloatLoc, 1, MPI_INT, FileCommunicator);

  int maxNfloatLoc = 0;
  for (int pp = 0; pp < fileNproc; pp++) {
    maxNfloatLoc = MAX(maxNfloatLoc, NfloatLoc[pp]);
  }

  MPI_Offset disp = 0;
  for (int pp = 0; pp < fileMyid; pp++)
    disp += (MPI_Offset)(NfloatLoc[pp] * sizeof(float));

  MPI_File_open(FileCommunicator, nomefile,
                MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &thefile);
  MPI_File_set_view(thefile, disp, MPI_FLOAT, MPI_FLOAT, (char *) "native", MPI_INFO_NULL);
  writeAllCPUParticlesValues(thefile, spec, maxNfloatLoc);
  MPI_File_close(&thefile);
  MPI_Comm_free(&FileCommunicator);
}

void OUTPUT_MANAGER::writeCPUParticlesValuesWritingGroups(std::string  fileName, SPECIE* spec){
  const int groupsize = particleGroupSize;
  const int bufsize = particleBufferSize*spec->Ncomp;

  MPI_File thefile;
  MPI_Status status;
  char *nomefile = new char[fileName.size() + 1];
  nomefile[fileName.size()] = 0;
  sprintf(nomefile, "%s", fileName.c_str());

  MPI_Comm groupCommunicator;
  MPI_Comm MPIFileCommunicator;
  MPI_Comm_split(MPI_COMM_WORLD, mygrid->myid / groupsize, 0, &groupCommunicator);
  MPI_Comm_split(MPI_COMM_WORLD, mygrid->myid%groupsize, 0, &MPIFileCommunicator);
  int groupMyid, groupNproc;
  MPI_Comm_size(groupCommunicator, &groupNproc);
  MPI_Comm_rank(groupCommunicator, &groupMyid);
  int mpiFileMyid, mpiFileNproc;
  MPI_Comm_size(MPIFileCommunicator, &mpiFileNproc);
  MPI_Comm_rank(MPIFileCommunicator, &mpiFileMyid);

  int Ncomp = spec->Ncomp;

  int totalFloatNumber = spec->Np*Ncomp;
  int* groupProcNumData = new int[groupNproc];
  groupProcNumData[groupMyid] = totalFloatNumber;
  MPI_Allgather(MPI_IN_PLACE, 1, MPI_INT, groupProcNumData, 1, MPI_INT, groupCommunicator);

  MPI_Offset disp = 0;
  long long int groupNumData = 0;
  if (groupMyid == 0){
    for (int i = 0; i < groupNproc; i++){
      groupNumData += groupProcNumData[i];
    }
    long long int *allGroupNumData = new long long int[mpiFileNproc];
    allGroupNumData[mpiFileMyid] = groupNumData;
    MPI_Allgather(MPI_IN_PLACE, 1, MPI_LONG_LONG_INT, allGroupNumData, 1, MPI_LONG_LONG_INT, MPIFileCommunicator);
    for (int i = 0; i < mpiFileMyid; i++){
      disp += allGroupNumData[i] * sizeof(float);
    }
    delete[] allGroupNumData;
  }

  float* data = new float[bufsize];

  if (groupMyid != 0){
    int numPackages = spec->Np/particleBufferSize;
    int resto = spec->Np%particleBufferSize;

    for (int i = 0; i < numPackages; i++){
      for (int p = 0; p < particleBufferSize; p++){
        int c;
        for (c = 0; c < Ncomp; c++){
          data[c + p*Ncomp] = (float)spec->ru(c, p + particleBufferSize*i);
        }
      }
      MPI_Send(data, bufsize, MPI_FLOAT, 0, i, groupCommunicator);
    }
    for (int p = 0; p < resto; p++){
      int c;
      for (c = 0; c < Ncomp; c++){
        data[c + p*Ncomp] = (float)spec->ru(c, p + particleBufferSize*numPackages);
      }
    }
    MPI_Send(data, resto*Ncomp, MPI_FLOAT, 0, numPackages, groupCommunicator);
  }
  else{
    MPI_File_open(MPIFileCommunicator, nomefile, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &thefile);
    MPI_File_set_view(thefile, disp, MPI_FLOAT, MPI_FLOAT, (char *) "native", MPI_INFO_NULL);

    int numPackages = spec->Np / particleBufferSize;
    int resto = spec->Np%particleBufferSize;
    for (int i = 0; i < numPackages; i++){
      for (int p = 0; p < particleBufferSize; p++){
        int c;
        for (c = 0; c < Ncomp; c++){
          data[c + p*Ncomp] = (float)spec->ru(c, p + particleBufferSize*i);
        }
      }
#ifndef DEBUG_NO_MPI_FILE_WRITE
      MPI_File_write(thefile, data, bufsize, MPI_FLOAT, &status);
#endif
    }
    for (int p = 0; p < resto; p++){
      int c;
      for (c = 0; c < Ncomp; c++){
        data[c + p*Ncomp] = (float)spec->ru(c, p + particleBufferSize*numPackages);
      }
    }
#ifndef DEBUG_NO_MPI_FILE_WRITE   
    MPI_File_write(thefile, data, resto*Ncomp, MPI_FLOAT, &status);
#endif
    MPI_Status status;

    std::vector<reqOutput> reqList;
    fillRequestList(bufsize, groupProcNumData, groupNproc, reqList);

    for (int i = 0; i < reqList.size(); i++){
      MPI_Recv(data, bufsize, MPI_FLOAT, reqList[i].task,  reqList[i].p, groupCommunicator, &status);
#ifndef DEBUG_NO_MPI_FILE_WRITE
      MPI_File_write(thefile, data, reqList[i].packageSize, MPI_FLOAT, &status);
#endif
    }


    MPI_File_close(&thefile);
  }
  MPI_Comm_free(&groupCommunicator);
  MPI_Comm_free(&MPIFileCommunicator);
  delete[] data;
  delete[] groupProcNumData;
}

void OUTPUT_MANAGER::writeCPUParticlesValuesFewFilesWritingGroups(std::string  fileName, SPECIE* spec, int NParticleToWrite, MPI_Comm outputCommunicator){
  const int groupsize = particleGroupSize;
  const int bufsize = particleBufferSize*spec->Ncomp;

  int myOutId;
  MPI_Comm_rank(outputCommunicator, &myOutId);
  int fileCommunicatorID = myOutId/multifileGroupSize;
  MPI_Comm FileCommunicator;
  MPI_Comm_split(outputCommunicator, fileCommunicatorID, 0, &FileCommunicator);
  int fileMyid, fileNproc;
  MPI_Comm_size(FileCommunicator, &fileNproc);
  MPI_Comm_rank(FileCommunicator, &fileMyid);

  MPI_Comm groupCommunicator;
  MPI_Comm_split(FileCommunicator, fileMyid / groupsize, 0, &groupCommunicator);
  int groupMyid, groupNproc;
  MPI_Comm_size(groupCommunicator, &groupNproc);
  MPI_Comm_rank(groupCommunicator, &groupMyid);

  MPI_Comm MPIFileCommunicator;
  MPI_Comm_split(FileCommunicator, fileMyid % groupsize, 0, &MPIFileCommunicator);
  int mpiFileMyid, mpiFileNproc;
  MPI_Comm_size(MPIFileCommunicator, &mpiFileNproc);
  MPI_Comm_rank(MPIFileCommunicator, &mpiFileMyid);

  MPI_File thefile;
  MPI_Status status;
  std::stringstream myFileName;
  myFileName << fileName << "." << std::setfill('0') << std::setw(5) << fileCommunicatorID;
  char *nomefile = new char[myFileName.str().size() + 1];
  nomefile[myFileName.str().size()] = 0;
  sprintf(nomefile, "%s", myFileName.str().c_str());

  const int Ncomp = spec->Ncomp;

  int totalFloatNumber = NParticleToWrite*Ncomp;
  int* groupProcNumData = new int[groupNproc];
  groupProcNumData[groupMyid] = totalFloatNumber;
  MPI_Allgather(MPI_IN_PLACE, 1, MPI_INT, groupProcNumData, 1, MPI_INT, groupCommunicator);
  long long int groupNumData = 0;
  for (int i = 0; i < groupNproc; i++){
    groupNumData += groupProcNumData[i];
  }

  MPI_Offset disp = 0;
  if (groupMyid == 0){
    long long int *allGroupNumData = new long long int[mpiFileNproc];
    allGroupNumData[mpiFileMyid] = groupNumData;
    MPI_Allgather(MPI_IN_PLACE, 1, MPI_LONG_LONG_INT, allGroupNumData, 1, MPI_LONG_LONG_INT, MPIFileCommunicator);
    for (int i = 0; i < mpiFileMyid; i++){
      disp += allGroupNumData[i] * sizeof(float);
    }
    delete[] allGroupNumData;

    MPI_File_open(MPIFileCommunicator, nomefile, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &thefile);
    MPI_File_set_view(thefile, disp, MPI_FLOAT, MPI_FLOAT, (char *) "native", MPI_INFO_NULL);

  }

  float* data = new float[bufsize];

  if (groupMyid != 0){
    int numPackages = NParticleToWrite/particleBufferSize;
    int resto = NParticleToWrite%particleBufferSize;

    for (int i = 0; i < numPackages; i++){
      for (int p = 0; p < particleBufferSize; p++){
        int c;
        for (c = 0; c < Ncomp; c++){
          data[c + p*Ncomp] = (float)spec->ru(c, p + particleBufferSize*i);
        }
      }
      MPI_Send(data, bufsize, MPI_FLOAT, 0, i, groupCommunicator);
    }
    for (int p = 0; p < resto; p++){
      int c;
      for (c = 0; c < Ncomp; c++){
        data[c + p*Ncomp] = (float)spec->ru(c, p + particleBufferSize*numPackages);
      }
    }
    MPI_Send(data, resto*Ncomp, MPI_FLOAT, 0, numPackages, groupCommunicator);
  }
  else{

    int numPackages = NParticleToWrite / particleBufferSize;
    int resto = NParticleToWrite%particleBufferSize;
    for (int i = 0; i < numPackages; i++){
      for (int p = 0; p < particleBufferSize; p++){
        int c;
        for (c = 0; c < Ncomp; c++){
          data[c + p*Ncomp] = (float)spec->ru(c, p + particleBufferSize*i);
        }
      }
#ifndef DEBUG_NO_MPI_FILE_WRITE
      MPI_File_write(thefile, data, bufsize, MPI_FLOAT, &status);
#endif
    }
    for (int p = 0; p < resto; p++){
      int c;
      for (c = 0; c < Ncomp; c++){
        data[c + p*Ncomp] = (float)spec->ru(c, p + particleBufferSize*numPackages);
      }
    }
#ifndef DEBUG_NO_MPI_FILE_WRITE
    MPI_File_write(thefile, data, resto*Ncomp, MPI_FLOAT, &status);
#endif
    MPI_Status status;

    std::vector<reqOutput> reqList;
    fillRequestList(bufsize, groupProcNumData, groupNproc, reqList);

    for (int i = 0; i < reqList.size(); i++){
      MPI_Recv(data, bufsize, MPI_FLOAT, reqList[i].task,  reqList[i].p, groupCommunicator, &status);
#ifndef DEBUG_NO_MPI_FILE_WRITE
      MPI_File_write(thefile, data, reqList[i].packageSize, MPI_FLOAT, &status);
#endif
    }


    MPI_File_close(&thefile);
  }
  MPI_Comm_free(&groupCommunicator);
  MPI_Comm_free(&MPIFileCommunicator);
  MPI_Comm_free(&FileCommunicator);
  delete[] data;
  delete[] groupProcNumData;
}

void OUTPUT_MANAGER::writeSpecPhaseSpace(std::string fileName, request req){

  SPECIE* spec = myspecies[req.target];
  //  int outputNComp, NCompFloat;
  //  bool flagMarker = spec->amIWithMarker();
  //  if (flagMarker){
  //    NCompFloat = spec->Ncomp - 1;
  //    outputNComp = spec->Ncomp + 1;
  //    NfloatLoc[mygrid->myid] = spec->Np*outputNComp;
  //  }
  //  else{
  //    NCompFloat = spec->Ncomp;
  //    outputNComp = spec->Ncomp;
  //    NfloatLoc[mygrid->myid] = spec->Np*outputNComp;
  //  }

  char *nomefile = new char[fileName.size() + 1];
  nomefile[fileName.size()] = 0;
  sprintf(nomefile, "%s", fileName.c_str());

  {
#if defined(PHASE_SPACE_USE_MPI_FILE_WRITE_ALL)
    int* NfloatLoc = new int[mygrid->nproc];
    int maxNfloatLoc = 0;
    NfloatLoc[mygrid->myid] = spec->Np*spec->Ncomp;
    MPI_Allgather(MPI_IN_PLACE, 1, MPI_INT, NfloatLoc, 1, MPI_INT, MPI_COMM_WORLD);

    for (int pp = 0; pp < mygrid->nproc; pp++) {
      maxNfloatLoc = MAX(maxNfloatLoc, NfloatLoc[pp]);
    }
    MPI_Offset disp = 0;
    for (int pp = 0; pp < mygrid->myid; pp++)
      disp += (MPI_Offset)(NfloatLoc[pp] * sizeof(float));
MPI_File thefile;
    MPI_File_open(MPI_COMM_WORLD, nomefile,
                  MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &thefile);
    MPI_File_set_view(thefile, disp, MPI_FLOAT, MPI_FLOAT, (char *) "native", MPI_INFO_NULL);
    writeAllCPUParticlesValues(thefile, spec, maxNfloatLoc);
    MPI_File_close(&thefile);
delete[] NfloatLoc;
#elif defined(PHASE_SPACE_USE_SEPARATE_FILES_MPI_FILE_WRITE_ALL)

    writeAllSeparateFilesParticlesValues(nomefile,spec);

#elif defined(PHASE_SPACE_USE_OUTPUT_WRITING_GROUPS)

    writeCPUParticlesValuesWritingGroups(nomefile,spec);

#elif defined(PHASE_SPACE_USE_HYBRID_OUTPUT)

    writeCPUParticlesValuesFewFilesWritingGroups(nomefile,spec, spec->Np, MPI_COMM_WORLD);

#elif defined(PHASE_SPACE_USE_MULTIFILE_OUTPUT)

    std::stringstream myFileName;
    myFileName << fileName << "." << std::setfill('0') << std::setw(5) << mygrid->myid;
    writeCPUParticlesValuesSingleFile(myFileName.str(),  spec);

#else
    int* NfloatLoc = new int[mygrid->nproc];
    int maxNfloatLoc = 0;
    NfloatLoc[mygrid->myid] = spec->Np*spec->Ncomp;
    MPI_Allgather(MPI_IN_PLACE, 1, MPI_INT, NfloatLoc, 1, MPI_INT, MPI_COMM_WORLD);

    for (int pp = 0; pp < mygrid->nproc; pp++) {
      maxNfloatLoc = MAX(maxNfloatLoc, NfloatLoc[pp]);
    }
    MPI_Offset disp = 0;
    for (int pp = 0; pp < mygrid->myid; pp++)
      disp += (MPI_Offset)(NfloatLoc[pp] * sizeof(float));
    MPI_File thefile;
    MPI_File_open(MPI_COMM_WORLD, nomefile,
                  MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &thefile);
    MPI_File_set_view(thefile, disp, MPI_FLOAT, MPI_FLOAT, (char *) "native", MPI_INFO_NULL);
    writeCPUParticlesValues(thefile, spec);
    MPI_File_close(&thefile);
delete[] NfloatLoc;
#endif


  }

  delete[] nomefile;

}
int OUTPUT_MANAGER::findNumberOfParticlesInSubdomain(request req){
  double rmin[3], rmax[3];
  for (int c = 0; c < 3; c++){
    rmin[c] = myDomains[req.domain]->rmin[c];
    rmax[c] = myDomains[req.domain]->rmax[c];
  }
  int counter = 0;
  double rr[3];
  SPECIE* spec = myspecies[req.target];
  for (int p = 0; p < spec->Np; p++){
    rr[0] = spec->r0(p);
    rr[1] = spec->r1(p);
    rr[2] = spec->r2(p);

    if (rmax[0] >= rr[0] && rmin[0] < rr[0]){
      if (mygrid->getDimensionality() < 2 || (rmax[1] >= rr[1] && rmin[1] < rr[1])){
        if (mygrid->getDimensionality() < 3 || (rmax[2] >= rr[2] && rmin[2] < rr[2])){
          counter++;
        }
      }
    }
  }
  return counter;
}

int OUTPUT_MANAGER::findNumberOfParticlesInSubdomainAndReorder(request req){
  double rmin[3], rmax[3];
  for (int c = 0; c < 3; c++){
    rmin[c] = myDomains[req.domain]->rmin[c];
    rmax[c] = myDomains[req.domain]->rmax[c];
  }
  int counter = 0;
  double rr[3], buffer;
  SPECIE* spec = myspecies[req.target];
  for (int p = 0; p < spec->Np; p++){
    rr[0] = spec->r0(p);
    rr[1] = spec->r1(p);
    rr[2] = spec->r2(p);

    if (rmax[0] >= rr[0] && rmin[0] < rr[0]){
      if (mygrid->getDimensionality() < 2 || (rmax[1] >= rr[1] && rmin[1] < rr[1])){
        if (mygrid->getDimensionality() < 3 || (rmax[2] >= rr[2] && rmin[2] < rr[2])){
          for(int c=0; c<spec->Ncomp; c++){
            buffer = spec->ru(c,counter);
            spec->ru(c,counter) = spec->ru(c,p);
            spec->ru(c,p) = buffer;
          }
          counter++;
        }
      }
    }
  }
  return counter;
}

void OUTPUT_MANAGER::writeSpecPhaseSpaceSubDomain(std::string fileName, request req){
  double rmin[3], rmax[3];
  for (int c = 0; c < 3; c++){
    rmin[c] = myDomains[req.domain]->rmin[c];
    rmax[c] = myDomains[req.domain]->rmax[c];
  }

  SPECIE* spec = myspecies[req.target];
  int shouldIWrite = false;
  shouldIWrite = amIInTheSubDomain(req);

  MPI_Comm outputCommunicator;
  MPI_Comm_split(MPI_COMM_WORLD, shouldIWrite, 0, &outputCommunicator);
  int myOutputID, outputNProc;
  MPI_Comm_rank(outputCommunicator, &myOutputID);
  MPI_Comm_size(outputCommunicator, &outputNProc);

  int outputNPart = findNumberOfParticlesInSubdomainAndReorder(req);
  char *nomefile = new char[fileName.size() + 1];
  nomefile[fileName.size()] = 0;
  sprintf(nomefile, "%s", fileName.c_str());

  if (shouldIWrite){
#if defined(PHASE_SPACE_USE_HYBRID_OUTPUT)
        writeCPUParticlesValuesFewFilesWritingGroups(nomefile,spec, outputNPart, outputCommunicator);

#elif defined(PHASE_SPACE_USE_MULTIFILE_OUTPUT)
    std::stringstream myFileName;
    myFileName << fileName << "." << std::setfill('0') << std::setw(5) << myOutputID;
    writeCPUParticlesValuesSingleFile(myFileName.str(), rmin, rmax, spec);
delete[] NfloatLoc;
#else

    int* NfloatLoc = new int[outputNProc];
    NfloatLoc[myOutputID] = outputNPart*spec->Ncomp;
    MPI_Offset disp = 0;
    for (int pp = 0; pp < myOutputID; pp++)
      disp += (MPI_Offset)(NfloatLoc[pp] * sizeof(float));
    MPI_File thefile;
    MPI_Status status;
    MPI_Allgather(MPI_IN_PLACE, 1, MPI_INT, NfloatLoc, 1, MPI_INT, outputCommunicator);

    MPI_File_open(outputCommunicator, nomefile, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &thefile);
    MPI_File_set_view(thefile, disp, MPI_FLOAT, MPI_FLOAT, (char *) "native", MPI_INFO_NULL);

    writeCPUParticlesValues(thefile, rmin, rmax, spec);
    MPI_File_close(&thefile);
    delete[] NfloatLoc;
#endif
  }
  MPI_Comm_free(&outputCommunicator);

  delete[] nomefile;
}

void OUTPUT_MANAGER::callSpecPhaseSpace(request req){
  std::string name = myspecies[req.target]->name;
  std::string nameBin;
  std::string outputName;
  if (myspecies[req.target]->amIWithMarker()){
    outputName = "PHASESPACE_WM";
  }
  else{
    outputName = "PHASESPACE";
  }

  if (req.domain == 0){
    nameBin = composeOutputName(outputDir, outputName, name, req.dtime, ".bin");
    writeSpecPhaseSpace(nameBin, req);
  }
  else{
    nameBin = composeOutputName(outputDir, outputName, name, myDomains[req.domain]->name, req.domain, req.dtime, ".bin");
    writeSpecPhaseSpaceSubDomain(nameBin, req);
  }
}


void OUTPUT_MANAGER::callDiag(request req){
  std::vector<SPECIE*>::const_iterator spec_iterator;
  double * ekinSpecies;
  size_t specie = 0;

  double tw = req.dtime;

  ekinSpecies = new double[myspecies.size()];

  double EE[3], BE[3];

  myfield->computeEnergyAndExtremes();
  double etotFields = myfield->total_energy[6];

  EE[0] = myfield->total_energy[0];
  EE[1] = myfield->total_energy[1];
  EE[2] = myfield->total_energy[2];
  BE[0] = myfield->total_energy[3];
  BE[1] = myfield->total_energy[4];
  BE[2] = myfield->total_energy[5];

  specie = 0;
  double etotKin = 0.0;
  for (spec_iterator = myspecies.begin(); spec_iterator != myspecies.end(); spec_iterator++){
    (*spec_iterator)->computeKineticEnergyWExtrems();
    ekinSpecies[specie++] = (*spec_iterator)->totalEnergy;
    etotKin += (*spec_iterator)->totalEnergy;
  }
  double etot = etotKin + etotFields;

  if (mygrid->myid == mygrid->master_proc){
    std::ofstream outStat;
    outStat.open(diagFileName.c_str(), std::ios::app);

    outStat << " " << std::setw(diagNarrowWidth) << req.itime << " " << std::setw(diagWidth) << tw;
    outStat << " " << std::setw(diagWidth) << etot;
    outStat << " " << std::setw(diagWidth) << EE[0] << " " << std::setw(diagWidth) << EE[1] << " " << std::setw(diagWidth) << EE[2];
    outStat << " " << std::setw(diagWidth) << BE[0] << " " << std::setw(diagWidth) << BE[1] << " " << std::setw(diagWidth) << BE[2];

    for (specie = 0; specie < myspecies.size(); specie++){
      outStat << " " << std::setw(diagWidth) << ekinSpecies[specie];
    }
    outStat << std::endl;

    outStat.close();

    std::ofstream ofField;
    ofField.open(extremaFieldFileName.c_str(), std::ios_base::app);
    myfield->output_extrems(req.itime, ofField);
    ofField.close();

  }

  for (spec_iterator = myspecies.begin(); spec_iterator != myspecies.end(); spec_iterator++){
    std::stringstream ss1;
    ss1 << outputDir << "/EXTREMES_" << (*spec_iterator)->name << ".dat";
    std::ofstream ofSpec;
    if (mygrid->myid == mygrid->master_proc)
      ofSpec.open(ss1.str().c_str(), std::ios_base::app);
    (*spec_iterator)->output_extrems(req.itime, ofSpec);
    ofSpec.close();
  }

  for (spec_iterator = myspecies.begin(); spec_iterator != myspecies.end(); spec_iterator++){
    std::string outNameSpec = composeOutputName(outputDir, "SPECTRUM", (*spec_iterator)->name, req.dtime, ".dat");
    std::ofstream ofSpec;
    if (mygrid->myid == mygrid->master_proc)
      ofSpec.open(outNameSpec.c_str());
    (*spec_iterator)->outputSpectrum(ofSpec);
    ofSpec.close();
  }

  delete[] ekinSpecies;
}

bool OUTPUT_MANAGER::isThePointInMyDomain(double rr[3]){
  if (rr[0] >= mygrid->rminloc[0] && rr[0] < mygrid->rmaxloc[0]){
    if (mygrid->getDimensionality() < 2 || (rr[1] >= mygrid->rminloc[1] && rr[1] < mygrid->rmaxloc[1])){
      if (mygrid->getDimensionality() < 3 || (rr[2] >= mygrid->rminloc[2] && rr[2] < mygrid->rmaxloc[2])){
        return true;
      }
    }
  }
  return false;
}
bool OUTPUT_MANAGER::shouldICreateHyperplane(int remains[3]){
  bool shouldICreateAnHyperplane=false;

  for(int c=0; c<mygrid->getDimensionality(); c++)
    shouldICreateAnHyperplane |= (remains[c]==0);
  return shouldICreateAnHyperplane;
}

bool OUTPUT_MANAGER::amIInTheSubDomain(request req){
  double rmin[3], rmax[3];
  for (int c = 0; c < 3; c++){
    rmin[c] = myDomains[req.domain]->rmin[c];
    rmax[c] = myDomains[req.domain]->rmax[c];
  }
  if (rmax[0] >= mygrid->rminloc[0] && rmin[0] < mygrid->rmaxloc[0]){
    if (mygrid->getDimensionality() < 2 || (rmax[1] >= mygrid->rminloc[1] && rmin[1] < mygrid->rmaxloc[1])){
      if (mygrid->getDimensionality() < 3 || (rmax[2] >= mygrid->rminloc[2] && rmin[2] < mygrid->rmaxloc[2])){
        return true;
      }
    }
  }
  return false;
}


void OUTPUT_MANAGER::nearestInt(double rr[], int *ri, int *globalri){
  int c;
  if (mygrid->isStretched()){
    for (c = 0; c < mygrid->getDimensionality(); c++){
      double mycsi = mygrid->unStretchGrid(rr[c], c);
      double xx = mygrid->dri[c] * (mycsi - mygrid->csiminloc[c]);
      ri[c] = (int)floor(xx + 0.5); //whole integer int
      mycsi = mygrid->unStretchGrid(rr[c], c);
      xx = mygrid->dri[c] * (mycsi - mygrid->csimin[c]);
      globalri[c] = (int)floor(xx + 0.5);
    }
    for (; c < 3; c++){
      ri[c] = globalri[c] = 0;
    }
  }
  else{
    for (c = 0; c < mygrid->getDimensionality(); c++){
      double xx = mygrid->dri[c] * (rr[c] - mygrid->rminloc[c]);
      ri[c] = (int)floor(xx + 0.5); //whole integer int
      xx = mygrid->dri[c] * (rr[c] - mygrid->rmin[c]);
      globalri[c] = (int)floor(xx + 0.5);
    }
    for (; c < 3; c++){
      ri[c] = globalri[c] = 0;
    }
  }
}
int OUTPUT_MANAGER::findLeftNeightbourPoint(double val, double* coords, int numcoords){
  if (numcoords <= 1)
    return 0;
  if (val <= coords[0])
    return 0;
  for (int i = 1; i < numcoords; i++){
    if (val < coords[i])
      return (i - 1);
  }
  return 0;
}

int OUTPUT_MANAGER::findRightNeightbourPoint(double val, double* coords, int numcoords){
  int indexMax = 0;
  if (numcoords <= 1)
    return 0;
  if (val >= coords[numcoords - 1])
    return (numcoords - 1);
  for (int i = (numcoords - 1); i >= 0; i--){
    if (val > coords[i])
      return (i + 1);
  }
  return 0;
}
void OUTPUT_MANAGER::setAndCheckRemains(int *remains, bool remainingCoord[3]){
  for (int c = 0; c < 3; c++){
    if (c < mygrid->getDimensionality())
      remains[c] = remainingCoord[c];
    else
      remains[c] = 0;
  }
}


void OUTPUT_MANAGER::findLocalIntegerBoundaries(double rmin[3], double rmax[3], int *imin, int *imax){
  for (int c = 0; c < 3; c++){
    imin[c] = findLeftNeightbourPoint(rmin[c], mygrid->cirloc[c], mygrid->Nloc[c]);
    imax[c] = findRightNeightbourPoint(rmax[c], mygrid->cirloc[c], mygrid->Nloc[c]);
  }
}
void OUTPUT_MANAGER::findGlobalIntegerBoundaries(double rmin[3], double rmax[3], int *imin, int *imax){
  for (int c = 0; c < 3; c++){
    imin[c] = findLeftNeightbourPoint(rmin[c], mygrid->cir[c], mygrid->NGridNodes[c]);
    imax[c] = findRightNeightbourPoint(rmax[c], mygrid->cir[c], mygrid->NGridNodes[c]);
  }
}
void OUTPUT_MANAGER::findNumberOfProcsWithinSubdomain(int *Nproc, int imin[3], int imax[3], int remains[3]){
  int mymin, mymax;
  for (int c = 0; c < 3; c++){
    if (remains[c]){
      if (imax[c] >= (mygrid->NGridNodes[c] - 1))
        mymax = mygrid->rnproc[c] - 1;
      for (int proc = 0; proc < mygrid->rnproc[c]; proc++){
        if (mygrid->rproc_imin[c][proc] <= imin[c] && mygrid->rproc_imax[c][proc] > imin[c])
          mymin = proc;
        if (mygrid->rproc_imin[c][proc] <= imax[c] && mygrid->rproc_imax[c][proc] > imax[c])
          mymax = proc;
      }

      Nproc[c] = mymax - mymin + 1;
    }
    else
      Nproc[c] = 1;
  }
}

void OUTPUT_MANAGER::findGlobalSubdomainUniquePointsNumber(int *uniqueN, int imin[3], int imax[3], int remains[3]){
  for (int c = 0; c < 3; c++){
    if (remains[c]){
      if (imin[c] < 0){
        printf("ERROR: subdomain minimum is wrong (component=%i, imin[c]=%i\n", c, imin[c]);
        exit(118);
      }

      if (imax[c] < (mygrid->uniquePoints[c]))
        uniqueN[c] = imax[c] - imin[c] + 1;
      else
        uniqueN[c] = mygrid->uniquePoints[c] - imin[c];
    }
    else{
      uniqueN[c] = 1;
    }
  }
}

void OUTPUT_MANAGER::findLocalSubdomainUniquePointsNumber(int *uniqueLocN, int locimin[3], int locimax[3], int remains[3]){
  for (int c = 0; c < 3; c++){
    if (remains[c]){
      if (locimin[c] < 0){
        printf("ERROR: subdomain local minimum is wrong (component=%i, locimin[c]=%i\n", c, locimin[c]);
        exit(118);
      }

      if (locimax[c] < (mygrid->uniquePointsloc[c]))
        uniqueLocN[c] = locimax[c] - locimin[c] + 1;
      else
        uniqueLocN[c] = mygrid->uniquePointsloc[c] - locimin[c];
    }
    else{
      uniqueLocN[c] = 1;
    }
  }

}

void OUTPUT_MANAGER::autoVisualDiag(){
  if (mygrid->myid  == mygrid->master_proc){
    std::cout << "*******OUTPUT MANAGER DEBUG***********" << std::endl;
    std::map< int, std::vector<request> >::iterator itMap;
    for (std::vector<int>::iterator itD = timeList.begin(); itD != timeList.end(); itD++) {
      std::map< int, std::vector<request> >::iterator itMapD;
      itMapD = allOutputs.find(*itD);
      if (itMapD != allOutputs.end()){
        std::cout << *itD << " ";
        for (std::vector<request>::iterator itR = itMapD->second.begin(); itR != itMapD->second.end(); itR++)
          std::cout << "(" << itR->type << "," << itR->target << ")";
        std::cout << std::endl;
      }

    }
    std::cout << "**************************************" << std::endl;
  }
}


int OUTPUT_MANAGER::getFieldGroupSize(){
  return fieldGroupSize;
}

int OUTPUT_MANAGER::getParticleGroupSize(){
  return particleGroupSize;
}

int OUTPUT_MANAGER::getParticleBufferSize(){
  return particleBufferSize;
}

void OUTPUT_MANAGER::setFieldGroupSize(int gsize){
  fieldGroupSize = gsize;
}

void OUTPUT_MANAGER::setParticleGroupSize(int gsize){
  particleGroupSize = gsize;
}

void OUTPUT_MANAGER::setParticleBufferSize(int bsize){
  particleBufferSize = bsize;
}
