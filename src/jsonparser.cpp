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


#include "jsonparser.h"
#include "jsonnames.h"
using namespace jsonParser;

bool jsonParser::isThisJsonMaster = true;
int jsonParser::inputVersion = 1;
bool jsonParser::checkVersion(Json::Value &document, int &version) {
  setInt(&version, document, "VERSION");
  if (version == 1)
    return true;
  else
    return false;
}

void jsonParser::lookForInputFile(int narg, char **args, std::string *inputFileName) {

  for (int i = 0; i < narg; i++) {
    if (std::string(args[i]) == _INPUT_LINE_INPUT_FILE_FLAG) {
      *inputFileName = std::string(args[i + 1]);
      break;
    }
  }

  if (doesFileExist(inputFileName->c_str())) {
    return;
  }
  else {
    *inputFileName = _DEFAULT_INPUT_FILE_NAME;
    if (doesFileExist(inputFileName->c_str())) {
      return;
    }
    else {
      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      if (rank == 0) {
        std::cout << "ERROR: no input file found!\n";
        std::cout.flush();
      }
      exitWithError(12);
    }

  }

}

std::string jsonParser::parseJsonInputFile(Json::Value &root, int narg, char **args) {
  std::string nomeFile;
  lookForInputFile(narg, args, &nomeFile);

  std::ifstream inputFile(nomeFile.c_str());
  std::stringstream buffer;
  buffer << inputFile.rdbuf();
  inputFile.close();
  Json::Reader reader;
  bool parsedSuccess = reader.parse(buffer.str().c_str(), root, false);

  if ((!parsedSuccess) && isThisJsonMaster) {
    std::cout << "Failed to parse JSON" << std::endl
      << reader.getFormattedErrorMessages()
      << std::endl;
    exit(1);
  }

  int masterProc = 0;
  setInt(&masterProc, root, _JSON_INT_MASTERPROC);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == masterProc)
    isThisJsonMaster = true;
  else
    isThisJsonMaster = false;

  int version = 1;
  if (!setInt(&version, root, _JSON_INT_VERSION)) {
    std::cout << "WARNING: version undefined, version = " << version <<
      " will be used as defautl!" << std::endl;
  }
  inputVersion = version;

  return nomeFile;
}


bool jsonParser::setInt(int *number, Json::Value &parent, const char* name) {
  bool outFlag= (!parent[name].isNull());
  if (outFlag) {
    if (parent[name].isInt())
      *number = parent[name].asInt();
  }
  return outFlag;
}


bool jsonParser::setDouble(double *number, Json::Value  &parent, const char* name) {
  bool outFlag= (!parent[name].isNull());
  if (outFlag) {
    if (parent[name].isDouble())
      *number = parent[name].asDouble();
  }
  return outFlag;
}
bool jsonParser::setBool(bool * number, Json::Value  &parent, const char* name) {
  bool outFlag= (!parent[name].isNull());
  if (outFlag) {
    if (parent[name].isBool())
      *number = parent[name].asBool();
  }
  return outFlag;
}

bool jsonParser::setString(std::string * number, Json::Value  &parent, const char* name) {
  bool outFlag= (!parent[name].isNull());
  if (outFlag) {
    if (parent[name].isString())
      *number = parent[name].asString();
  }
  return outFlag;
}


bool jsonParser::setValue(Json::Value &child, Json::Value &parent, const char* name) {
  bool outFlag= (!parent[name].isNull());
  if (outFlag) {
    child = parent[name];
  }
  return outFlag;
}


int jsonParser::getDimensionality(Json::Value &document, int defaultDimensionality) {
  int dim = defaultDimensionality;
  const char* name = INT_DIMENSIONS;
  if ((!setInt(&dim, document, name)) && isThisJsonMaster)
    std::cout << "dimensions not set in JSON input file!\n";
  return dim;
}

void jsonParser::setRadiationFriction(Json::Value &document, GRID *grid) {
  bool isFriction = false, isLambda = false;
  double lambda0;
  isLambda = setDouble(&lambda0, document, _JSON_DOUBLE_FRICTION_LENGTH_CM);
  setBool(&isFriction, document, _JSON_BOOL_IS_FRICTION);

  if (isFriction&&isLambda) {
    grid->enableRadiationFriction();
    grid->setLambda0(lambda0);
  }
}

void jsonParser::setXrange(Json::Value &parent, GRID *grid) {
  double min = -1.0, max = 1.0;
  const char* name = _JSON_DOUBLEARRAY_X_RANGE;
  if (!parent[name].isNull()) {
    if (parent[name].isArray()) {
      min = parent[name][0].asDouble();
      max = parent[name][1].asDouble();
    }
  }
  grid->setXrange(min, max);
}
void jsonParser::setYrange(Json::Value &parent, GRID *grid) {
  double min = -1.0, max = 1.0;
  const char* name = _JSON_DOUBLEARRAY_Y_RANGE;
  if (!parent[name].isNull()) {
    if (parent[name].isArray()) {
      min = parent[name][0].asDouble();
      max = parent[name][1].asDouble();
    }
  }
  grid->setYrange(min, max);
}
void jsonParser::setZrange(Json::Value &parent, GRID *grid) {
  double min = -1.0, max = 1.0;
  const char* name = _JSON_DOUBLEARRAY_Z_RANGE;
  if (!parent[name].isNull()) {
    if (parent[name].isArray()) {
      min = parent[name][0].asDouble();
      max = parent[name][1].asDouble();
    }
  }
  grid->setZrange(min, max);
}

void jsonParser::setNCells(Json::Value &parent, GRID *grid) {
  int Nx, Ny, Nz;
  Nx = Ny = Nz = 1;
  const char* name = _JSON_INT_N_CELLS_;
  if (!parent[name].isNull()) {
    if (parent[name].isArray()) {
      Nx = parent[name][0].asInt();
      Ny = parent[name][1].asInt();
      Nz = parent[name][2].asInt();
    }
  }
  grid->setNCells(Nx, Ny, Nz);

}
void jsonParser::setNprocs(Json::Value &document, GRID *grid) {
  int nProcY = 1, nProcZ = 1;
  setInt(&nProcY, document, _JSON_INT_N_PROC_Y_);
  setInt(&nProcZ, document, _JSON_INT_N_PROC_Z_);

  grid->setNProcsAlongY(nProcY);
  grid->setNProcsAlongZ(nProcZ);
}

void jsonParser::setSimulationTime(Json::Value &document, GRID *grid) {
  double simulationTime;
  setDouble(&simulationTime, document, _JSON_DOUBLE_SIMULATION_TIME_);
  grid->setSimulationTime(simulationTime);
}
void jsonParser::setMasterProc(Json::Value  &document, GRID *grid) {
  int masterProc = 0;
  setInt(&masterProc, document, _JSON_INT_MASTERPROC);
  grid->setMasterProc(masterProc);
}

void jsonParser::setCourantFactor(Json::Value  &document, GRID *grid) {
  double courantFactor = 0.98;
  setDouble(&courantFactor, document, _JSON_DOUBLE_COURANT_FACTOR);
  grid->setCourantFactor(courantFactor);
}
void jsonParser::setBoundaryConditions(Json::Value &parent, GRID *grid) {
  std::string  name1 = _JSON_OBJ_BOUNDARIES_;
  std::string  xCondition, yCondition, zCondition;
  int xFlag, yFlag, zFlag;
  xFlag = xPBC;
  yFlag = yPBC;
  zFlag = zPBC;

  Json::Value boundaries;
  int boundariesFlag;
  if (setValue(boundaries, parent, name1.c_str())) {
    if (boundaries.size() == 3) {
      xCondition = boundaries[0].asString();
      yCondition = boundaries[1].asString();
      zCondition = boundaries[2].asString();
      if (!xCondition.compare(_JSON_TAG_PERIODIC_BC_))
        xFlag = xPBC;
      else if (!xCondition.compare(_JSON_TAG_OPEN_BC_))
        xFlag = xOpen;
      else if (!xCondition.compare(_JSON_TAG_PML_))
        xFlag = xPML;

      if (!yCondition.compare(_JSON_TAG_PERIODIC_BC_))
        yFlag = yPBC;
      else if (!yCondition.compare(_JSON_TAG_OPEN_BC_))
        yFlag = yOpen;
      else if (!yCondition.compare(_JSON_TAG_PML_))
        yFlag = yPML;

      if (!zCondition.compare(_JSON_TAG_PERIODIC_BC_))
        zFlag = zPBC;
      else if (!zCondition.compare(_JSON_TAG_OPEN_BC_))
        zFlag = zOpen;
      else if (!zCondition.compare(_JSON_TAG_PML_))
        zFlag = zPML;

    }
  }
  boundariesFlag = xFlag | yFlag | zFlag;
  grid->setBoundaries(boundariesFlag);
}

void jsonParser::setDumpControl(Json::Value &parent, GRID *mygrid) {
  mygrid->dumpControl.doRestart = false;
  mygrid->dumpControl.doDump = false;

  Json::Value restartObject;
  if (setValue(restartObject, parent, _JSON_OBJ_RESTART_)) {

    std::string folderName;
    if (setString(&folderName, restartObject, _JSON_STRING_DUMP_FOLDER_NAME)) {
      mygrid->setDumpPath(folderName);
    }
    else {
      mygrid->setDumpPath(_STRING_DUMP_DEFAULT_PATH);
    }

    setBool(&mygrid->dumpControl.doRestart, restartObject, _JSON_BOOL_RESTART_);

    setDouble(&mygrid->dumpControl.dumpEvery, restartObject, _JSON_DOUBLE_DUMPEVERY_);

    setBool(&mygrid->dumpControl.doDump, restartObject, _JSON_BOOL_DODUMP_);

    setInt(&mygrid->dumpControl.restartFromDump, restartObject, _JSON_INT_RESTART_FROM_DUMP_);

  }
}


void jsonParser::setStretchedGrid(Json::Value &document, GRID *grid) {
  Json::Value  stretching;
  bool isEnabled;
  if (setValue(stretching, document, _JSON_OBJ_STRETCHED_GRID_)) {
    if (inputVersion == 1)
      isEnabled = true;

    setBool(&isEnabled, stretching, _JSON_BOOL_ENABLED);
    if (isEnabled) {

      grid->enableStretchedGrid();
    }
    std::string  name2;
    Json::Value stretching1D;

    name2 = _JSON_OBJ_X_STTETCHING_;
    if (setValue(stretching1D, stretching, name2.c_str())) {
      std::string  name3 = _JSON_OBJ_LEFT_STRETCHING_;
      Json::Value stretchingLeft;
      if (setValue(stretchingLeft, stretching1D, name3.c_str())) {
        double limit;
        int NCells;
        setInt(&NCells, stretchingLeft, _JSON_INT_NCELL_STRETCHING_);
        setDouble(&limit, stretchingLeft, _JSON_DOUBLE_LIMIT_STRETCHING);
        grid->setXandNxLeftStretchedGrid(limit, NCells);
      }
      name3 = _JSON_OBJ_RIGHT_STRETCHING_;
      Json::Value stretchingRight;
      if (setValue(stretchingRight, stretching1D, name3.c_str())) {
        double limit;
        int NCells;
        setInt(&NCells, stretchingRight, _JSON_INT_NCELL_STRETCHING_);
        setDouble(&limit, stretchingRight, _JSON_DOUBLE_LIMIT_STRETCHING);
        grid->setXandNxRightStretchedGrid(limit, NCells);
      }
    }

    name2 = _JSON_OBJ_Y_STTETCHING_;
    if (setValue(stretching1D, stretching, name2.c_str())) {
      std::string  name3 = _JSON_OBJ_LEFT_STRETCHING_;
      Json::Value stretchingLeft;
      if (setValue(stretchingLeft, stretching1D, name3.c_str())) {
        double limit;
        int NCells;
        setInt(&NCells, stretchingLeft, _JSON_INT_NCELL_STRETCHING_);
        setDouble(&limit, stretchingLeft, _JSON_DOUBLE_LIMIT_STRETCHING);
        grid->setYandNyLeftStretchedGrid(limit, NCells);
      }
      name3 = _JSON_OBJ_RIGHT_STRETCHING_;
      Json::Value stretchingRight;
      if (setValue(stretchingRight, stretching1D, name3.c_str())) {
        double limit;
        int NCells;
        setInt(&NCells, stretchingRight, _JSON_INT_NCELL_STRETCHING_);
        setDouble(&limit, stretchingRight, _JSON_DOUBLE_LIMIT_STRETCHING);
        grid->setYandNyRightStretchedGrid(limit, NCells);
      }
    }

    name2 = _JSON_OBJ_Z_STTETCHING_;
    if (setValue(stretching1D, stretching, name2.c_str())) {
      std::string  name3 = _JSON_OBJ_LEFT_STRETCHING_;
      Json::Value stretchingLeft;
      if (setValue(stretchingLeft, stretching1D, name3.c_str())) {
        double limit;
        int NCells;
        setInt(&NCells, stretchingLeft, _JSON_INT_NCELL_STRETCHING_);
        setDouble(&limit, stretchingLeft, _JSON_DOUBLE_LIMIT_STRETCHING);
        grid->setZandNzLeftStretchedGrid(limit, NCells);
      }
      name3 = _JSON_OBJ_RIGHT_STRETCHING_;
      Json::Value stretchingRight;
      if (setValue(stretchingRight, stretching1D, name3.c_str())) {
        double limit;
        int NCells;
        setInt(&NCells, stretchingRight, _JSON_INT_NCELL_STRETCHING_);
        setDouble(&limit, stretchingRight, _JSON_DOUBLE_LIMIT_STRETCHING);
        grid->setZandNzRightStretchedGrid(limit, NCells);
      }
    }
  }
}
void jsonParser::setMovingWindow(Json::Value  &document, GRID *grid) {
  Json::Value movingWindow;
  bool isEnambled;
  if (setValue(movingWindow, document, _JSON_OBJ_MOVING_WINDOW_)) {
    if (inputVersion == 1)
      isEnambled = true;

    setBool(&isEnambled, movingWindow, _JSON_BOOL_ENABLED);
    if (isEnambled) {
      double start = 0;
      setDouble(&start, movingWindow, _JSON_DOUBLE_START_MW_);
      grid->setStartMovingWindow(start);

      double beta;
      if (setDouble(&beta, movingWindow, _JSON_DOUBLE_BETA_MW_)) {
        grid->setBetaMovingWindow(beta);
      }
      int frequency;
      if (setInt(&frequency, movingWindow, _JSON_DOUBLE_FREQUENCY_MW_)) {
        grid->setFrequencyMovingWindow(frequency);
      }
    }
  }

}

bool jsonParser::setLaserType(laserPulse *pulse1, Json::Value  &mylaser) {
  std::string name2 = _JSON_STRING_TYPE_;
  std::string type;
  bool flag  = setString(&type, mylaser, name2.c_str());
  if (flag) {
    if (type == _LASERTYPEVALUE_COS_PLANE_WAVE_) {
      pulse1->type = COS2_PLANE_WAVE;
    }
    else if (type == _LASERTYPEVALUE_GAUSSIAN_) {
      pulse1->type = GAUSSIAN;
    }
    else if (type == _LASERTYPEVALUE_PLANE_WAVE_) {
      pulse1->type = PLANE_WAVE;
    }
    else if (type == _LASERTYPEVALUE_COS2_PLATEAU_PLANE_WAVE_) {
      pulse1->type = COS2_PLATEAU_PLANE_WAVE;
    }
    else if (type == _LASERTYPEVALUE_LG_) {
      pulse1->type = LAGUERRE_GAUSSIAN;
    }
    else
      flag = false;
  }
  return flag;
}

bool jsonParser::setLaserPolarization(laserPulse *pulse1, Json::Value  &mylaser) {
  std::string name2 = _JSON_STRING_POLARIZATION_;
  std::string polarization;
  bool flag = setString(&polarization, mylaser, name2.c_str());
  if (flag ) {
    if (polarization == _LASERPOLARIZATIONVALUE_P_) {
      pulse1->setPPolarization();
    }
    else if (polarization == _LASERPOLARIZATIONVALUE_S_) {
      pulse1->setSPolarization();
    }
    else if (polarization == _LASERPOLARIZATIONVALUE_C_) {
      pulse1->setCircularPolarization();
    }
    else
      flag = false;
  }
  return flag;
}

bool jsonParser::setLaserDurationFWHM(laserPulse *pulse1, Json::Value  &mylaser) {
  std::string name2 = _JSON_DOUBLE_LASER_DURATION_FWHM_;
  double durationFWHM;
  bool flag = setDouble(&durationFWHM, mylaser, name2.c_str());
  if (flag ) {
    pulse1->setDurationFWHM(durationFWHM);
  }
  return flag;
}

bool jsonParser::setLaserInitialPosition(laserPulse *pulse1, Json::Value  &mylaser) {
  std::string name2 = _JSON_DOUBLE_LASER_INITIAL_POSITION_;
  double initialPosition;
  bool flag  = setDouble(&initialPosition, mylaser, name2.c_str());
  if (flag) {
    pulse1->setPulseInitialPosition(initialPosition);
  }
  return flag;
}

bool jsonParser::setLaserAmplitude(laserPulse *pulse1, Json::Value  &mylaser) {
  std::string name2 = _JSON_DOUBLE_LASER_A_;
  double amplitude;
  bool flag = setDouble(&amplitude, mylaser, name2.c_str());
  if (flag ) {
    pulse1->setNormalizedAmplitude(amplitude);
  }
  return flag;
}

bool jsonParser::setLaserWaist(laserPulse *pulse1, Json::Value  &mylaser) {
  std::string name2 = _JSON_DOUBLE_LASER_WAIST_;
  double waist;
  bool flag = setDouble(&waist, mylaser, name2.c_str());
  if (flag ) {
    pulse1->setWaist(waist);
  }
  return flag;
}
bool jsonParser::setLaserFocusPosition(laserPulse *pulse1, Json::Value  &mylaser) {
  std::string name2 = _JSON_DOUBLE_LASER_FOCUS_POSITION_;
  double focusPosition;
  bool flag = setDouble(&focusPosition, mylaser, name2.c_str());
  if (flag ) {
    pulse1->setFocusPosition(focusPosition);
  }
  return flag;
}
bool jsonParser::setLaserLambda(laserPulse *pulse1, Json::Value  &mylaser) {
  std::string name2 = _JSON_DOUBLE_LASER_LAMBDA_;
  double lambda;
  bool flag = setDouble(&lambda, mylaser, name2.c_str());
  if (flag ) {
    pulse1->setLambda(lambda);
  }
  return flag;
}

bool jsonParser::setLaserRotation(laserPulse *pulse1, Json::Value  &mylaser) {
  std::string name2 = _JSON_BOOL_LASER_ROTATION_;
  bool rotation;
  bool flag = false;
  flag = setBool(&rotation, mylaser, name2.c_str());
  if (rotation) {

    double angle = 0, center = 0;
    name2 = _JSON_DOUBLE_ROTATION_ANGLE_;
    setDouble(&angle, mylaser, name2.c_str());

    name2 = _JSON_DOUBLE_ROTATION_CENTER_;
    setDouble(&center, mylaser, name2.c_str());

    pulse1->setRotationAngleAndCenter(2.0*M_PI*(angle / 360.0), center);
  }
  return flag;
}

bool jsonParser::setLaserRiseTime(laserPulse *pulse1, Json::Value  &mylaser) {
  std::string name2 = _JSON_DOUBLE_LASER_RISE_TIME_;
  double riseTime;
  bool flag = setDouble(&riseTime, mylaser, name2.c_str());
  if (flag) {
    pulse1->setRiseTime(riseTime);
  }
  return flag;
}

bool jsonParser::setLaserLG_l(laserPulse *pulse1, Json::Value  &mylaser) {
  std::string name2 = _JSON_INT_LG_L_;
  int LG_l;
  bool flag = setInt(&LG_l, mylaser, name2.c_str());
  if (flag) {
    pulse1->setLaguerreGaussian_l(LG_l);
  }
  return flag;
}

bool jsonParser::setLaserLG_m(laserPulse *pulse1, Json::Value  &mylaser) {
  std::string name2 = _JSON_INT_LG_M_;
  int LG_m;
  bool flag = setInt(&LG_m, mylaser, name2.c_str());
  if (flag) {
    pulse1->setLaguerreGaussian_m(LG_m);
  }
  return flag;
}

bool jsonParser::checkLaserBoolFlags(laserPulseBoolFlags flags, laserPulse *pulse) {
      if (!flags.type) {
    return false;
  }
  switch (pulse->type) {
  case GAUSSIAN:
    if (!(flags.a && flags.waist && flags.duration)) {
      return false;
    }
    break;
  case PLANE_WAVE:
    if (!(flags.a)) {
      return false;
    }
    break;
  case COS2_PLANE_WAVE:
    if (!(flags.a && flags.duration)) {
      return false;
    }
    break;
  case COS2_PLATEAU_PLANE_WAVE:
    if (!(flags.a && flags.duration)) {
      return false;
    }
    break;
  default:
    break;
  }  
  return true;
}

void jsonParser::setLaserPulses(Json::Value &document, EM_FIELD *emfield) {
  std::string  name1 = _JSON_OBJ_ARRAY_LASER_;
  Json::Value lasers;



  if (setValue(lasers, document, name1.c_str())) {
    if (lasers.isArray()) {
      std::string  name2;

      for (unsigned int index = 0; index < lasers.size(); index++) {

        Json::Value myLaser = lasers[index];

        name2 = _JSON_BOOL_ENABLED;
        bool enabled = false;
        setBool(&enabled, myLaser, name2.c_str());
        if (enabled) {
          laserPulse *pulse1 = new(laserPulse);

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

          setLaserLG_l(pulse1, myLaser);
          setLaserLG_m(pulse1, myLaser);


          if (checkLaserBoolFlags(flags, pulse1)) {
            emfield->addPulse(pulse1);
          }
          else if (isThisJsonMaster) {
            std::cout << "warning: laser #" << index << " is incorrectly defined\n";
          }

          delete pulse1;
        }
      }
    }
  }
}

void jsonParser::setPoissonSolver(Json::Value &document, GRID *grid){

  bool withPoissonSolver = false;
  setBool(&withPoissonSolver, document, _JSON_BOOL_WITH_POISSON);
  if (withPoissonSolver)
    grid->setWithPoisson();

  bool neutraliseDensity = false;
  setBool(&neutraliseDensity, document, _JSON_BOOL_NEUTRALISE_DENSITY_FOR_POISSON);
  if (neutraliseDensity)
    grid->setAutoNeutraliseDensity();

}


int jsonParser::findPlasmaFunction(std::string plasmaFunction) {
  for (int i = 0; i < PLASMA::maxdF; i++) {
    //std::cout << "-" << plasmaFunction << "-  PLASMA::dFNames[" << i << "] = " <<  PLASMA::dFNames[i]<< std::endl;
    if (!plasmaFunction.compare(PLASMA::dFNames[i])) {
      //std::cout << "i = " << i<< "  YES!\n";
      return i;
    }
  }
  return -1;
}

void jsonParser::setPlasmas(Json::Value &document, std::map<std::string, PLASMA*> &map) {
  std::string  name1 = _JSON_OBJARRAY_PLASMA_;
  Json::Value plasmas;

  if (setValue(plasmas, document, name1.c_str()) && plasmas.isArray()) {
    for (unsigned int index = 0; index < plasmas.size(); index++) {
      Json::Value myPlasma = plasmas[index];

      std::string plasmaName;
      setString(&plasmaName, myPlasma, _JSON_STRING_PLASMA_NAME);

      if (map.find(plasmaName) != map.end() && isThisJsonMaster) {
        std::cout << "Warning! Plasma " << plasmaName << " is defined multiple times!" << std::endl;
      }

      std::string plasmaFunction;
      bool isPlasmaFunctionSet = false;

      int plasmaFunctionIndex;
      if (setString(&plasmaFunction, myPlasma, _JSON_STRING_PLASMA_DENSITYFUNCTION)) {
        plasmaFunctionIndex = findPlasmaFunction(plasmaFunction);
        if (plasmaFunctionIndex >= 0)
          isPlasmaFunctionSet = true;
      }
      if (isPlasmaFunctionSet) {
        map[plasmaName] = new PLASMA();
        map[plasmaName]->density_function = PLASMA::dFPoint[plasmaFunctionIndex];

        double tdouble, t2double;
        Json::Value range;

        if (setValue(range, myPlasma, _JSON_DOUBLEARRAY_PLASMA_X_RANGEBOX) && (range.size() == 2))
          map[plasmaName]->setXRangeBox(range[0].asDouble(), range[1].asDouble());

        if (setValue(range, myPlasma, _JSON_DOUBLEARRAY_PLASMA_Y_RANGEBOX) && (range.size() == 2))
          map[plasmaName]->setYRangeBox(range[0].asDouble(), range[1].asDouble());

        if (setValue(range, myPlasma, _JSON_DOUBLEARRAY_PLASMA_Z_RANGEBOX) && (range.size() == 2))
          map[plasmaName]->setZRangeBox(range[0].asDouble(), range[1].asDouble());


        if (setDouble(&tdouble, myPlasma, _JSON_DOUBLE_DENSITY_LAMBDA)) {
          if (setDouble(&t2double, myPlasma, _JSON_DOUBLE_PLASMA_DENSITY))
            map[plasmaName]->setDensityCoefficient(t2double, tdouble);
        }
        else {
          if (setDouble(&tdouble, myPlasma, _JSON_DOUBLE_PLASMA_DENSITY))
            map[plasmaName]->setDensityCoefficient(tdouble);
        }
        if (setDouble(&tdouble, myPlasma, _JSON_DOUBLE_PLASMA_RAMP_LENGTH))
          map[plasmaName]->setRampLength(tdouble);

        if (setDouble(&tdouble, myPlasma, _JSON_DOUBLE_PLASMA_RAMP_MIN_DENSITY))
          map[plasmaName]->setRampMinDensity(tdouble);

        if (setDouble(&tdouble, myPlasma, _JSON_DOUBLE_PLASMA_LEFT_RAMP_LENGTH))
          map[plasmaName]->setLeftRampLength(tdouble);

        if (setDouble(&tdouble, myPlasma, _JSON_DOUBLE_PLASMA_RIGHT_RAMP_LENGTH))
          map[plasmaName]->setRightRampLength(tdouble);

        if (setDouble(&tdouble, myPlasma, _JSON_DOUBLE_PLASMA_SCALE_LENGTH))
          map[plasmaName]->setScaleLength(tdouble);

        if (setDouble(&tdouble, myPlasma, _JSON_DOUBLE_PLASMA_LEFT_SCALE_LENGTH))
          map[plasmaName]->setLeftScaleLength(tdouble);

        if (setDouble(&tdouble, myPlasma, _JSON_DOUBLE_PLASMA_RIGHT_SCALE_LENGTH))
          map[plasmaName]->setRightScaleLength(tdouble);

        if (setDouble(&tdouble, myPlasma, _JSON_DOUBLE_PLASMA_LEFT_RAMP_MIN_DENSITY))
          map[plasmaName]->setLeftRampMinDensity(tdouble);

        if (setDouble(&tdouble, myPlasma, _JSON_DOUBLE_PLASMA_RIGHT_RAMP_MIN_DENSITY))
          map[plasmaName]->setRightRampMinDensity(tdouble);

        if (PLASMA::isGrating(plasmaFunctionIndex)) {
          double g_depth = 0;
          double g_period = 1.0;
          double g_phase = 0.0;

          double* additionalParams = new double[3];

          setDouble(&g_depth, myPlasma, _JSON_DOUBLE_PLASMA_GRATING_DEPTH);
          setDouble(&g_period, myPlasma, _JSON_DOUBLE_PLASMA_GRATING_PERIOD);
          setDouble(&g_phase, myPlasma, _JSON_DOUBLE_PLASMA_GRATING_PHASE);

          additionalParams[0] = g_depth;
          additionalParams[1] = g_period;
          additionalParams[2] = g_phase;

          map[plasmaName]->setAdditionalParams(additionalParams);
        }

        if (PLASMA::isPillar2D(plasmaFunctionIndex)) {
          double p2d_dx = 1.0;
          double p2d_dy = 1.0;
          double p2d_r = 0.1;

          double* additionalParams = new double[3];


          setDouble(&p2d_dx, myPlasma, _JSON_DOUBLE_PLASMA_PILLARS2D_DX);
          setDouble(&p2d_dy, myPlasma, _JSON_DOUBLE_PLASMA_PILLARS2D_DY);
          setDouble(&p2d_r,  myPlasma, _JSON_DOUBLE_PLASMA_PILLARS2D_R);

          additionalParams[0] = p2d_dx;
          additionalParams[1] = p2d_dy;
          additionalParams[2] = p2d_r;

          map[plasmaName]->setAdditionalParams(additionalParams);
        }
        if (PLASMA::isNanotubes2D(plasmaFunctionIndex)) {
          double n2d_width = 1.0;
          double n2d_dist = 1.0;
          double n2d_depth = 0.1;

          double* additionalParams = new double[3];


          setDouble(&n2d_width, myPlasma, _JSON_DOUBLE_PLASMA_NANOTUBES2D_WIDTH);
          setDouble(&n2d_dist,  myPlasma, _JSON_DOUBLE_PLASMA_NANOTUBES2D_DIST);
          setDouble(&n2d_depth, myPlasma, _JSON_DOUBLE_PLASMA_NANOTUBES2D_DEPTH);

          additionalParams[0] = n2d_width;
          additionalParams[1] = n2d_dist;
          additionalParams[2] = n2d_depth;

          map[plasmaName]->setAdditionalParams(additionalParams);
        }

        if (PLASMA::isFoils2D(plasmaFunctionIndex)) {
          double f2d_width = 1.0;
          double f2d_dist = 1.0;

          double* additionalParams = new double[2];


          setDouble(&f2d_width, myPlasma, _JSON_DOUBLE_PLASMA_FOILS2D_WIDTH);
          setDouble(&f2d_dist, myPlasma, _JSON_DOUBLE_PLASMA_FOILS2D_DIST);

          additionalParams[0] = f2d_width;
          additionalParams[1] = f2d_dist;

          map[plasmaName]->setAdditionalParams(additionalParams);
        }
        if (PLASMA::isUser1(plasmaFunctionIndex) || PLASMA::isUser2(plasmaFunctionIndex)) {
          double width = 0;
          double depth = 0;
          double position = 0;
          double temperature = 0;

          double* additionalParams = new double[4];

          setDouble(&width, myPlasma, _JSON_DOUBLE_PLASMA_USER1_WIDTH);
          setDouble(&depth, myPlasma, _JSON_DOUBLE_PLASMA_USER1_DEPTH);
          setDouble(&position, myPlasma, _JSON_DOUBLE_PLASMA_USER1_POSITION);
          setDouble(&temperature, myPlasma, _JSON_DOUBLE_PLASMA_USER1_TEMPERATURE);

          additionalParams[0] = width;
          additionalParams[1] = depth;
          additionalParams[2] = position;
          additionalParams[3] = temperature;

          map[plasmaName]->setAdditionalParams(additionalParams);
        }
      }
      else {
        if (isThisJsonMaster)
          std::cout << "Warning! Plasma " << plasmaName <<
          " has no valid density function. It will be ignored." << std::endl;
      }

    }
  }
}

bool jsonParser::checkSpecEssentials(Json::Value &child, std::map<std::string, PLASMA*> plasmas) {
  bool isEnambled = false;
  bool isThereName = false;
  bool isThereType = false;
  bool isTherePlasma = false;
  bool isTherePPC = false;

  std::string dummy;

  if (inputVersion == 1)
    isEnambled = true;

  setBool(&isEnambled, child, _JSON_BOOL_ENABLED);

  isThereName = setString(&dummy, child, _JSON_STRING_SPECIE_NAME);

  if (setString(&dummy, child, _JSON_STRING_SPECIE_PLASMANAME)) {
    if (plasmas.find(dummy) != plasmas.end())
      isTherePlasma = true;
  }

  Json::Value ppc;

  if (setValue(ppc, child, _JSON_INTARRAY3_PARTICLES_PER_CELL)) {
    if (ppc.isArray() && (ppc.size() == 3))
      if (ppc[0].asInt() >= 0 && ppc[1].asInt() >= 0 && ppc[2].asInt() >= 0)
        isTherePPC = true;
  }

  if (setString(&dummy, child, _JSON_STRING_SPECIES_TYPE)) {
    if (dummy.compare(SPECIES_TYPEVALUE_ELECTRON) ||
      dummy.compare(SPECIES_TYPEVALUE_POSITRON) ||
      dummy.compare(SPECIES_TYPEVALUE_ION))
      isThereType = true;
  }

  return isEnambled&&isThereName&&isThereType&&isTherePlasma&&isTherePPC;
}

bool jsonParser::addDistribution(std::string distName, Json::Value &child, my_rng_generator& ext_rng, SPECIE* spec) {
  bool isThereMomentum = false;
  bool isDistOk = false;

  tempDistrib dist;

  Json::Value momentum;
  if (setValue(momentum, child, _JSON_DOUBLEARRAY_ADD_MOMENTUM) && momentum.isArray() && (momentum.size() == 3))
    isThereMomentum = true;

  Json::Value params;
  if (!distName.compare(DISTRIBUTIONVALUE_MAXWELL)) {
    if (setValue(params, child, _JSON_DOUBLEARRAY_DISTRIBUTION_PARAMS) && params.isArray() && (params.size() >= 1)) {
      dist.setMaxwell(params[0].asDouble());
      isDistOk = true;
    }
  }
  else if (!distName.compare(DISTRIBUTIONVALUE_JUTTNER)) {
    if (setValue(params, child, _JSON_DOUBLEARRAY_DISTRIBUTION_PARAMS) && params.isArray() && (params.size() >= 1)) {
      dist.setJuttner(params[0].asDouble());
      isDistOk = true;
    }
  }
  else if (!distName.compare(DISTRIBUTIONVALUE_WATERBAG)) {
    if (setValue(params, child, _JSON_DOUBLEARRAY_DISTRIBUTION_PARAMS) && params.isArray() && (params.size() >= 1)) {
      dist.setWaterbag(params[0].asDouble());
      isDistOk = true;
    }
  }
  else if (!distName.compare(DISTRIBUTIONVALUE_WATERBAG3TEMP)) {
    if (setValue(params, child, _JSON_DOUBLEARRAY_DISTRIBUTION_PARAMS) && params.isArray() && (params.size() >= 3)) {
      dist.setWaterbag3Temp(params[0].asDouble(), params[1].asDouble(), params[2].asDouble());
      isDistOk = true;
    }
  }
  else if (!distName.compare(DISTRIBUTIONVALUE_UNFORM_SPHERE)) {
    if (setValue(params, child, _JSON_DOUBLEARRAY_DISTRIBUTION_PARAMS) && params.isArray() && (params.size() >= 1)) {
      dist.setUnifSphere(params[0].asDouble());
      isDistOk = true;
    }
  }
  else if (!distName.compare(DISTRIBUTIONVALUE_SUPERGAUSSIAN)) {
    if (setValue(params, child, _JSON_DOUBLEARRAY_DISTRIBUTION_PARAMS) && params.isArray() && (params.size() >= 2)) {
      dist.setSupergaussian(params[0].asDouble(), params[1].asDouble());
      isDistOk = true;
    }
  }
  else if (!distName.compare(DISTRIBUTIONVALUE_SPECIAL)) {
    if (setValue(params, child, _JSON_DOUBLEARRAY_DISTRIBUTION_PARAMS) && params.isArray() && (params.size() >= 1)) {
      dist.setSpecial(params[0].asDouble());
      isDistOk = true;
    }
  }

  if (isDistOk) {
    if (isThereMomentum) {
      spec->add_momenta(ext_rng, momentum[0].asDouble(), momentum[1].asDouble(), momentum[2].asDouble(), dist);
    }
    else {
      spec->add_momenta(ext_rng, 0.0, 0.0, 0.0, dist);
    }
    return true;
  }

  return false;
}

void jsonParser::setSpecies(Json::Value &document, std::vector<SPECIE*> &species, std::map<std::string, PLASMA*> plasmas, GRID* myGrid, my_rng_generator& ext_rng) {
  std::cout.flush();
  Json::Value specList;
  if (setValue(specList, document, _JSON_OBJARRAY_SPECIES) && specList.isArray()) {
    for (unsigned int index = 0; index < specList.size(); index++) {
      Json::Value mySpecies = specList[index];
      if (checkSpecEssentials(mySpecies, plasmas)) {
        SPECIE *newSpec = new SPECIE(myGrid);

        std::string dummy;
        Json::Value ppc;
        setString(&dummy, mySpecies, _JSON_STRING_SPECIE_PLASMANAME);

        newSpec->plasma = *plasmas[dummy];

        setString(&dummy, mySpecies, _JSON_STRING_SPECIES_TYPE);

        if (!dummy.compare(SPECIES_TYPEVALUE_ION))
          newSpec->type = ION;
        else if (!dummy.compare(SPECIES_TYPEVALUE_POSITRON))
          newSpec->type = POSITRON;
        else
          newSpec->type = ELECTRON;

        setString(&dummy, mySpecies, _JSON_STRING_SPECIE_NAME);
        newSpec->setName(dummy);

        setValue(ppc, mySpecies, _JSON_INTARRAY3_PARTICLES_PER_CELL);
        newSpec->setParticlesPerCellXYZ(ppc[0].asInt(), ppc[1].asInt(), ppc[2].asInt());

        if (newSpec->type == ION) {
          double A, Z;
          if (setDouble(&A, mySpecies, _JSON_DOUBLE_IONS_A))
            newSpec->A = A;
          if (setDouble(&Z, mySpecies, _JSON_DOUBLE_IONS_Z))
            newSpec->Z = Z;
        }

        bool isMarker = false;
        setBool(&isMarker, mySpecies, _JSON_BOOL_IS_MARKER);
        if (isMarker)
          newSpec->addMarker();

        bool isTest = false;
        setBool(&isTest, mySpecies, _JSON_BOOL_IS_TEST);
        if (isTest)
          newSpec->setTestSpecies();

        bool isFrozen = false;
        setBool(&isFrozen, mySpecies, _JSON_BOOL_IS_FROZEN);
        if (isFrozen)
          newSpec->setFrozenSpecies();

        bool isQuiet = false;
        setBool(&isQuiet, mySpecies, _JSON_BOOL_IS_QUIET);
        if (isQuiet)
          newSpec->setQuietStart();

        int quiteShuffle = 1;
        setInt(&quiteShuffle, mySpecies, _JSON_INT_QUIET_SHUFFLE);
        newSpec->setQuietShuffle(quiteShuffle);

        newSpec->creation();

        if (setString(&dummy, mySpecies, _JSON_STRING_DISTRIBUTION)) {
          if ((!addDistribution(dummy, mySpecies, ext_rng, newSpec)) && isThisJsonMaster) {
            std::cout << "Warning: temperature distribution for species " << index
              << " is incorrectly defined" << std::endl;
          }
        }

        species.push_back(newSpec);

      }
      else if (isThisJsonMaster) {
        std::cout << "Warning: species " << index <<
          " is not defined correctly." << std::endl;
      }
    }
  }
}

void jsonParser::setDomains(Json::Value &document, std::map<std::string, outDomain*> &domains) {
  std::string  name1 = _JSON_OBJECTARRAY_DOMAINS;
  Json::Value domainList;
  if (setValue(domainList, document, name1.c_str()) && domainList.isArray()) {
    for (unsigned int index = 0; index < domainList.size(); index++) {
      Json::Value domain = domainList[index];

      std::string name;
      if (setString(&name, domain, _JSON_STRING_DOMAIN_NAME) && (domains.find(name) == domains.end())) {
        outDomain* newDomain = new outDomain();
        newDomain->setName(name);

        Json::Value freeDims;
        if (setValue(freeDims, domain, _JSON_INTARRAY3_FREE_DIM) && freeDims.isArray() && (freeDims.size() == 3)) {
          newDomain->setFreeDimensions(freeDims[0].asBool(), freeDims[1].asBool(), freeDims[2].asBool());
        }

        Json::Value pointCoord;
        if (setValue(pointCoord, domain, _JSON_DOUBLEARRAY3_POINT_COORDINATE) && pointCoord.isArray() && (pointCoord.size() == 3)) {
          newDomain->setPointCoordinate(pointCoord[0].asDouble(), pointCoord[1].asDouble(), pointCoord[2].asDouble());

        }

        Json::Value xRange;
        if (setValue(xRange, domain, _JSON_INTARRAY2_DOMAIN_X_RANGE) && xRange.isArray() && (xRange.size() == 2)) {
          newDomain->setXRange(xRange[0].asDouble(), xRange[1].asDouble());
        }

        Json::Value yRange;
        if (setValue(yRange, domain, _JSON_INTARRAY2_DOMAIN_Y_RANGE) && yRange.isArray() && (yRange.size() == 2)) {
          newDomain->setYRange(yRange[0].asDouble(), yRange[1].asDouble());
        }

        Json::Value zRange;
        if (setValue(zRange, domain, _JSON_INTARRAY2_DOMAIN_Z_RANGE) && zRange.isArray() && (zRange.size() == 2)) {
          newDomain->setZRange(zRange[0].asDouble(), zRange[1].asDouble());
        }

        domains[name] = newDomain;
      }
      else if (isThisJsonMaster) {
        std::cout << "Warning: output domain " << index << " is defined incorrectly. It will be ignored. " << std::endl;
      }
    }
  }
}

bool jsonParser::isOutTypeOk(std::string type) {
  bool answ = (!type.compare(OUTPUTTYPEVALUE_E)) || (!type.compare(OUTPUTTYPEVALUE_B)) || (!type.compare(OUTPUTTYPEVALUE_EB)) ||
    (!type.compare(OUTPUTTYPEVALUE_DENSITY)) || (!type.compare(OUTPUTTYPEVALUE_CURRENT)) || (!type.compare(OUTPUTTYPEVALUE_PHASESPACE)) || (!type.compare(OUTPUTTYPEVALUE_DIAG));
  return answ;
}

bool jsonParser::isSpecNameOk(std::string specName, std::vector<SPECIE*> &species) {
  for (std::vector<SPECIE*>::const_iterator spec_iterator = species.begin(); spec_iterator != species.end(); spec_iterator++) {
    if (!(*spec_iterator)->name.compare(specName))
      return true;
  }
  return false;
}

bool jsonParser::addOutReq(std::string type, OUTPUT_MANAGER &manager, std::map<std::string, outDomain*>  &domains, jsonParser::outRequest outRequestVals) {
  if (outRequestVals.isDomainName) {
    if (outRequestVals.isAt) {
      if (!type.compare(OUTPUTTYPEVALUE_E)) {
        manager.addEFieldAt(domains[outRequestVals.domainName], outRequestVals.at);
      }
      else if (!type.compare(OUTPUTTYPEVALUE_B)) {
        manager.addBFieldAt(domains[outRequestVals.domainName], outRequestVals.at);
      }
      else if (!type.compare(OUTPUTTYPEVALUE_EB)) {
        manager.addEBFieldAt(domains[outRequestVals.domainName], outRequestVals.at);
      }
      else if (!type.compare(OUTPUTTYPEVALUE_DENSITY)) {
        manager.addSpeciesDensityAt(domains[outRequestVals.domainName], outRequestVals.specName, outRequestVals.at);
      }
      else if (!type.compare(OUTPUTTYPEVALUE_CURRENT)) {
        manager.addCurrentAt(domains[outRequestVals.domainName], outRequestVals.at);
      }
      else if (!type.compare(OUTPUTTYPEVALUE_PHASESPACE)) {
        manager.addSpeciesPhaseSpaceAt(domains[outRequestVals.domainName], outRequestVals.specName, outRequestVals.at);
      }
      else if (!type.compare(OUTPUTTYPEVALUE_DIAG)) {
        manager.addDiagAt(outRequestVals.at);
      }

    }

    if (outRequestVals.isEvery) {
      double every = outRequestVals.every;
      double from;
      if (outRequestVals.isFrom)
        from = outRequestVals.from;
      else
        from = 0.0;

      if (outRequestVals.isTo) {
        if (!type.compare(OUTPUTTYPEVALUE_E)) {
          manager.addEFieldFromTo(domains[outRequestVals.domainName], from, every, outRequestVals.to);
        }
        else if (!type.compare(OUTPUTTYPEVALUE_B)) {
          manager.addBFieldFromTo(domains[outRequestVals.domainName], from, every, outRequestVals.to);
        }
        else if (!type.compare(OUTPUTTYPEVALUE_EB)) {
          manager.addEBFieldFromTo(domains[outRequestVals.domainName], from, every, outRequestVals.to);
        }
        else if (!type.compare(OUTPUTTYPEVALUE_DENSITY)) {
          manager.addSpeciesDensityFromTo(domains[outRequestVals.domainName], outRequestVals.specName, from, every, outRequestVals.to);
        }
        else if (!type.compare(OUTPUTTYPEVALUE_CURRENT)) {
          manager.addCurrentFromTo(domains[outRequestVals.domainName], from, every, outRequestVals.to);
        }
        else if (!type.compare(OUTPUTTYPEVALUE_PHASESPACE)) {
          manager.addSpeciesPhaseSpaceFromTo(domains[outRequestVals.domainName], outRequestVals.specName, from, every, outRequestVals.to);
        }
        else if (!type.compare(OUTPUTTYPEVALUE_DIAG)) {
          manager.addDiagFromTo(from, every, outRequestVals.to);
        }

      }
      else {
        if (!type.compare(OUTPUTTYPEVALUE_E)) {
          manager.addEFieldFrom(domains[outRequestVals.domainName], from, every);
        }
        else if (!type.compare(OUTPUTTYPEVALUE_B)) {
          manager.addBFieldFrom(domains[outRequestVals.domainName], from, every);
        }
        else if (!type.compare(OUTPUTTYPEVALUE_EB)) {
          manager.addEBFieldFrom(domains[outRequestVals.domainName], from, every);
        }
        else if (!type.compare(OUTPUTTYPEVALUE_DENSITY)) {
          manager.addSpeciesDensityFrom(domains[outRequestVals.domainName], outRequestVals.specName, from, every);
        }
        else if (!type.compare(OUTPUTTYPEVALUE_CURRENT)) {
          manager.addCurrentFrom(domains[outRequestVals.domainName], from, every);
        }
        else if (!type.compare(OUTPUTTYPEVALUE_PHASESPACE)) {
          manager.addSpeciesPhaseSpaceFrom(domains[outRequestVals.domainName], outRequestVals.specName, from, every);
        }
        else if (!type.compare(OUTPUTTYPEVALUE_DIAG)) {
          manager.addDiagFrom(from, every);
        }

      }
    }
  }
  else {
    if (outRequestVals.isAt) {
      if (!type.compare(OUTPUTTYPEVALUE_E)) {
        manager.addEFieldAt(outRequestVals.at);
      }
      else if (!type.compare(OUTPUTTYPEVALUE_B)) {
        manager.addBFieldAt(outRequestVals.at);
      }
      else if (!type.compare(OUTPUTTYPEVALUE_EB)) {
        manager.addEBFieldAt(outRequestVals.at);
      }
      else if (!type.compare(OUTPUTTYPEVALUE_DENSITY)) {
        manager.addSpeciesDensityAt(outRequestVals.specName, outRequestVals.at);
      }
      else if (!type.compare(OUTPUTTYPEVALUE_CURRENT)) {
        manager.addCurrentAt(outRequestVals.at);
      }
      else if (!type.compare(OUTPUTTYPEVALUE_PHASESPACE)) {
        manager.addSpeciesPhaseSpaceAt(outRequestVals.specName, outRequestVals.at);
      }
      else if (!type.compare(OUTPUTTYPEVALUE_DIAG)) {
        manager.addDiagAt(outRequestVals.at);
      }

    }
    if (outRequestVals.isEvery) {
      double every = outRequestVals.every;
      double from;
      if (outRequestVals.isFrom)
        from = outRequestVals.from;
      else
        from = 0.0;

      if (outRequestVals.isTo) {
        if (!type.compare(OUTPUTTYPEVALUE_E)) {
          manager.addEFieldFromTo(from, every, outRequestVals.to);
        }
        else if (!type.compare(OUTPUTTYPEVALUE_B)) {
          manager.addBFieldFromTo(from, every, outRequestVals.to);
        }
        else if (!type.compare(OUTPUTTYPEVALUE_EB)) {
          manager.addEBFieldFromTo(from, every, outRequestVals.to);
        }
        else if (!type.compare(OUTPUTTYPEVALUE_DENSITY)) {
          manager.addSpeciesDensityFromTo(outRequestVals.specName, from, every, outRequestVals.to);
        }
        else if (!type.compare(OUTPUTTYPEVALUE_CURRENT)) {
          manager.addCurrentFromTo(from, every, outRequestVals.to);
        }
        else if (!type.compare(OUTPUTTYPEVALUE_PHASESPACE)) {
          manager.addSpeciesPhaseSpaceFromTo(outRequestVals.specName, from, every, outRequestVals.to);
        }
        else if (!type.compare(OUTPUTTYPEVALUE_DIAG)) {
          manager.addDiagFromTo(from, every, outRequestVals.to);
        }

      }
      else {
        if (!type.compare(OUTPUTTYPEVALUE_E)) {
          manager.addEFieldFrom(from, every);
        }
        else if (!type.compare(OUTPUTTYPEVALUE_B)) {
          manager.addBFieldFrom(from, every);
        }
        else if (!type.compare(OUTPUTTYPEVALUE_EB)) {
          manager.addEBFieldFrom(from, every);
        }
        else if (!type.compare(OUTPUTTYPEVALUE_DENSITY)) {
          manager.addSpeciesDensityFrom(outRequestVals.specName, from, every);
        }
        else if (!type.compare(OUTPUTTYPEVALUE_CURRENT)) {
          manager.addCurrentFrom(from, every);
        }
        else if (!type.compare(OUTPUTTYPEVALUE_PHASESPACE)) {
          manager.addSpeciesPhaseSpaceFrom(outRequestVals.specName, from, every);
        }
        else if (!type.compare(OUTPUTTYPEVALUE_DIAG)) {
          manager.addDiagFrom(from, every);
        }


      }
    }
  }
  return true;
}



void jsonParser::setOutputRequests(Json::Value &document, OUTPUT_MANAGER &manager, std::map<std::string, outDomain*>  &domains, std::vector<SPECIE*> &species) {
  std::string  name1 = _JSON_OBJARRAY_OUTPUT;
  Json::Value outputList;
  if (setValue(outputList, document, name1.c_str()) && outputList.isArray()) {
    for (unsigned int index = 0; index < outputList.size(); index++) {
      Json::Value outReq = outputList[index];

      std::string type;

      bool correct = false;
      if (setString(&type, outReq, _JSON_STRING_OUTPUT_TYPE) && isOutTypeOk(type)) {

        jsonParser::outRequest outRequestVals;

        outRequestVals.isFrom = setDouble(&outRequestVals.from, outReq, _JSON_DOUBLE_OUTPUT_FROM);
        outRequestVals.isTo = setDouble(&outRequestVals.to, outReq, _JSON_DOUBLE_OUTPUT_TO);
        outRequestVals.isEvery = setDouble(&outRequestVals.every, outReq, _JSON_DOUBLE_OUTPUT_EVERY);
        outRequestVals.isAt = setDouble(&outRequestVals.at, outReq, _JSON_DOUBLE_OUTPUT_AT);

        outRequestVals.isSpecName = setString(&outRequestVals.specName, outReq, _JSON_STRING_OUTPUT_SPECIES_NAME);
        outRequestVals.isDomainName = setString(&outRequestVals.domainName, outReq, _JSON_STRING_IN_OUTPUT_DOMAIN_NAME);

        bool isSpecRequired = false;
        if (!type.compare(OUTPUTTYPEVALUE_DENSITY))
          isSpecRequired = true;
        else if (!type.compare(OUTPUTTYPEVALUE_PHASESPACE))
          isSpecRequired = true;

        bool isSpecValid = !isSpecRequired || (!outRequestVals.isSpecName) || isSpecNameOk(outRequestVals.specName, species);
        bool isDomainValid = (!outRequestVals.isDomainName) || (domains.find(outRequestVals.domainName) != domains.end());
        bool isAtLeastOneTime = (outRequestVals.isAt) || (outRequestVals.isEvery);

        if (isAtLeastOneTime && isSpecValid && isDomainValid) {
          correct = addOutReq(type, manager, domains, outRequestVals);
        }

      }

      if ((!correct) && isThisJsonMaster) {
        std::cout << "Output request " << index << " has unrecognized format!" << std::endl;
      }

    }
  }
}

void jsonParser::setOutputDirPath(Json::Value &document, OUTPUT_MANAGER &manager) {
  std::string  folderName;
  if (setString(&folderName, document, _JSON_STRING_OUTPUT_FOLDER_NAME))
    manager.setOutputPath(folderName);
  else {
    manager.setOutputPath(_STRING_OUTPUT_DIR_DEFAULT_PATH);
  }
}

void jsonParser::setOutputParameters(Json::Value &document, OUTPUT_MANAGER &manager) {
  int bufferInt;
  if (setInt(&bufferInt, document, _JSON_INT_OUTPUT_MULTIFILE_GROUP_SIZE))
    manager.setMultifileGroupSize(bufferInt);

  if (setInt(&bufferInt, document, _JSON_INT_OUTPUT_FIELD_GROUP_SIZE))
    manager.setFieldGroupSize(bufferInt);

  if (setInt(&bufferInt, document, _JSON_INT_OUTPUT_PARTICLES_GROUP_SIZE))
    manager.setParticleGroupSize(bufferInt);

  if (setInt(&bufferInt, document, _JSON_INT_OUTPUT_PARTICLES_BUFFER_SIZE))
    manager.setParticleBufferSize(bufferInt);

}



