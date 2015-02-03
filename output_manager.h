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

#ifndef __OUTPUT_MANAGER_H__
#define __OUTPUT_MANAGER_H__

#define _USE_MATH_DEFINES
#define _CRT_SECURE_NO_WARNINGS

//#define PHASE_SPACE_USE_MPI_FILE_WRITE_ALL
//#define PHASE_SPACE_USE_SEPARATE_FILES_MPI_FILE_WRITE_ALL
//#define PHASE_SPACE_USE_OUTPUT_WRITING_GROUPS
#define PHASE_SPACE_USE_HYBRID_OUTPUT
//#define PHASE_SPACE_USE_MULTIFILE_OUTPUT
#define PHASE_SPACE_GROUP_SIZE 128 //1024 was the best case for 4096 MPI_TASKS on 1024 BlueGeneQ cores
#define NPARTICLE_BUFFER_SIZE 1000000

#define FIELDS_USE_SEPARATE_FILES_MACROGROUPS

//#define FIELDS_TEST_SEPARATE_FILES_MACROGROUPS
//#define FIELDS_USE_MPI_FILE_OUTPUT
//#define FIELDS_USE_MPI_FILE_WRITE_ALL
//#define FIELDS_USE_OUTPUT_WRITING_GROUPS
//#define FIELDS_USE_INDIVIDUAL_FILE_OUTPUT
//#define FIELDS_USE_MULTI_FILE
#define FIELDS_GROUP_SIZE 64
#define MACRO_CPUGROUP_FOR_MULTIFILE 1024

#include <mpi.h>
#include <iomanip>
#if defined(_MSC_VER)
#include <cstdint>
#include <cstdlib>
#else
#include <stdint.h>
#include <stdlib.h>
#endif
#include <vector>
#include <list>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include "commons.h"
#if defined(USE_BOOST)
#include <boost/filesystem.hpp>
#endif
#if defined(USE_HDF5)
#include <hdf5.h>
#endif
#include "grid.h"
#include "structures.h"
#include "current.h"
#include "em_field.h"
#include "particle_species.h"
#include "utilities.h"



#define OUTPUT_SIZE_TIME_DIAG 7

#define SPEC_DIAG_COMP 14
#define FIELD_DIAG_COMP 14

//#define DEBUG_NO_MPI_FILE_WRITE

enum diagType{
  OUT_E_FIELD,
  OUT_B_FIELD,
  OUT_SPEC_PHASE_SPACE,
  OUT_DIAG,
  OUT_CURRENT,
  OUT_EB_PROBE,
  OUT_SPEC_DENSITY
};

enum whichFieldOut{
  WHICH_E_ONLY,
  WHICH_B_ONLY,
  WHICH_E_AND_B
};


struct request{
  double dtime;
  int itime;
  diagType type;
  int target;
  int domain;
};

struct emProbe{
  double coordinates[3];
  std::string name;
  std::string fileName;
  emProbe();
  bool compareProbes(emProbe* rhs);
  void setPointCoordinate(double X, double Y, double Z);
  void setName(std::string iname);
};

struct reqOutput{
  int task;
  int p;
  int packageSize;
} ;

struct outDomain{
  double coordinates[3];
  bool remainingCoord[3], subselection;
  double rmin[3], rmax[3];
  bool overrideFlag;
  bool followMovingWindowFlag;
  std::string name;
  outDomain();
  bool compareDomains(outDomain* rhs);
  void setFreeDimensions(bool flagX, bool flagY, bool flagZ);
  void setPointCoordinate(double X, double Y, double Z);
  void setName(std::string iname);
  void setXRange(double min, double max);
  void setYRange(double min, double max);
  void setZRange(double min, double max);
  void followMovingWindow();
};

bool requestCompTime(const request &first, const request &second);
bool requestCompUnique(const request &first, const request &second);
bool compOutput(const reqOutput &first, const reqOutput &second);


class OUTPUT_MANAGER
{

public:
  OUTPUT_MANAGER(GRID* _mygrid, EM_FIELD* _myfield, CURRENT* _mycurrent, std::vector<SPECIE*> _myspecies);
  ~OUTPUT_MANAGER();

  void initialize(std::string _outputDir);
  void initialize();
  void close();

  void setOutputPath(std::string dirName);

  void addEBFieldFrom(double startTime, double frequency);
  void addEBFieldAt(double atTime);
  void addEBFieldFromTo(double startTime, double frequency, double endTime);

  void addEBFieldFrom(outDomain* _domain, double startTime, double frequency);
  void addEBFieldAt(outDomain* _domain, double atTime);
  void addEBFieldFromTo(outDomain* _domain, double startTime, double frequency, double endTime);

  void addEFieldFrom(double startTime, double frequency);
  void addEFieldAt(double atTime);
  void addEFieldFromTo(double startTime, double frequency, double endTime);

  void addBFieldFrom(double startTime, double frequency);
  void addBFieldAt(double atTime);
  void addBFieldFromTo(double startTime, double frequency, double endTime);

  void addEFieldFrom(outDomain* _domain, double startTime, double frequency);
  void addEFieldAt(outDomain* _domain, double atTime);
  void addEFieldFromTo(outDomain* _domain, double startTime, double frequency, double endTime);

  void addBFieldFrom(outDomain* _domain, double startTime, double frequency);
  void addBFieldAt(outDomain* _domain, double atTime);
  void addBFieldFromTo(outDomain* _domain, double startTime, double frequency, double endTime);

  void addEBFieldProbeFrom(emProbe* Probe, double startTime, double frequency);
  void addEBFieldProbeAt(emProbe* Probe, double atTime);
  void addEBFieldProbeFromTo(emProbe* Probe, double startTime, double frequency, double endTime);

  void addSpeciesDensityFrom(std::string name, double startTime, double frequency);
  void addSpeciesDensityAt(std::string name, double atTime);
  void addSpeciesDensityFromTo(std::string name, double startTime, double frequency, double endTime);

  void addSpeciesDensityFrom(outDomain* Plane, std::string name, double startTime, double frequency);
  void addSpeciesDensityAt(outDomain* Plane, std::string name, double atTime);
  void addSpeciesDensityFromTo(outDomain* Plane, std::string name, double startTime, double frequency, double endTime);

  void addCurrentFrom(double startTime, double frequency);
  void addCurrentAt(double atTime);
  void addCurrentFromTo(double startTime, double frequency, double endTime);

  void addCurrentFrom(outDomain* Plane, double startTime, double frequency);
  void addCurrentAt(outDomain* Plane, double atTime);
  void addCurrentFromTo(outDomain* Plane, double startTime, double frequency, double endTime);

  void addSpeciesPhaseSpaceFrom(outDomain* domain_in, std::string name, double startTime, double frequency);
  void addSpeciesPhaseSpaceAt(outDomain* domain_in, std::string name, double atTime);
  void addSpeciesPhaseSpaceFromTo(outDomain* domain_in, std::string name, double startTime, double frequency, double endTime);

  void addSpeciesPhaseSpaceFrom(std::string name, double startTime, double frequency);
  void addSpeciesPhaseSpaceAt(std::string name, double atTime);
  void addSpeciesPhaseSpaceFromTo(std::string name, double startTime, double frequency, double endTime);

  void addDiagFrom(double startTime, double frequency);
  void addDiagAt(double atTime);
  void addDiagFromTo(double startTime, double frequency, double endTime);


  void autoVisualDiag();

  void callDiags(int istep);

  int getFieldGroupSize();
  int getParticleGroupSize();
  int getParticleBufferSize();
  void setFieldGroupSize(int gsize);
  void setParticleGroupSize(int gsize);
  void setParticleBufferSize(int bsize);

private:
  GRID* mygrid;
  EM_FIELD* myfield;
  CURRENT* mycurrent;
  std::vector<SPECIE*> myspecies;
  std::vector<emProbe*> myEMProbes;
  std::vector<outDomain*> myDomains;

  bool isThereGrid;
  bool isThereCurrent;
  bool isThereField;
  bool isThereSpecList;

  std::string outputDir;
  bool amIInit;

  bool isThereDiag;
  bool isThereEMProbe;

  int fieldGroupSize;
  int multifileGroupSize;
  int particleGroupSize;
  int particleBufferSize;

  bool checkGrid();
  bool checkEMField();
  bool checkCurrent();
  bool checkSpecies();
  bool isThePointInMyDomain(double rr[3]);
  bool shouldICreateHyperplane(int remains[]);
  bool amIInTheSubDomain(request req);
  void nearestInt(double rr[3], int *ri, int *globalri);
  void setAndCheckRemains(int *remains, bool remainingCoord[]);
  void findLocalIntegerBoundaries(double rmin[3], double rmax[3], int *imin, int *imax);
  void findGlobalIntegerBoundaries(double rmin[3], double rmax[3], int *imin, int *imax);
  void findNumberOfProcsWithinSubdomain(int *Nproc, int imin[3], int imax[3], int remains[3]);
  void findGlobalSubdomainUniquePointsNumber(int *uniqueN, int imin[3], int imax[3], int remains[3]);
  void findLocalSubdomainUniquePointsNumber(int *uniqueLocN, int locimin[3], int locimax[3], int remains[3]);
  int findSpecIndexInMyspeciesVector(std::string name);
  int findProbeIndexInMyprobeVector(emProbe *newProbe);
  int findDomainIndexInMydomainsVector(outDomain *newDomain);
  std::string diagFileName;
  std::string extremaFieldFileName;
  std::vector<std::string> extremaSpecFileNames;

  std::list<request> requestList;
  std::vector<int> timeList;
  std::map< int, std::vector<request> > allOutputs;

  static const int diagWideWidth = 16;
  static const int diagNarrowWidth = 6;
  static const int diagWidth = 10;

  int getIntegerTime(double dtime);

  void addEBField(outDomain* _domain, double startTime, double frequency, double endTime, whichFieldOut whichOut);
  void addEBFieldProbe(emProbe* Probe, double startTime, double frequency, double endTime);
  void addSpeciesDensity(outDomain* domain_in, std::string name, double startTime, double frequency, double endTime);
  void addCurrent(outDomain* domain_in, double startTime, double frequency, double endTime);
  void addSpeciesPhaseSpace(outDomain* domain_in, std::string name, double startTime, double frequency, double endTime);
  void addDiag(double startTime, double frequency, double endTime);

  void addRequestToList(std::list<request>& timeList, diagType type, int target, int domain, double startTime, double frequency, double endTime);

  void prepareOutputMap();

  void createDiagFile();
  void createExtremaFiles();
  void createEMProbeFiles();

  void processOutputEntry(request req);

  std::string composeOutputName(std::string dir, std::string out, std::string opt, double time, std::string ext);
  std::string composeOutputName(std::string dir, std::string out, std::string opt1, std::string opt2, int domain, double time, std::string ext);
  void appendIDtoFileName(char *nomefile, std::string fileName, int ID);
  void writeEMFieldBinaryHDF5(std::string fileName, request req);

  void callEMFieldProbe(request req);
  void interpolateEBFieldsToPosition(double pos[3], double E[3], double B[3]);

  void writeGridFieldSubDomain(std::string fileName, request req);
  void callEMFieldDomain(request req);

  void writeSpecDensity(std::string fileName, request req);
  void writeSpecDensitySubDomain(std::string fileName, request req);
  void callSpecDensity(request req);

  void writeCurrent(std::string fileName, request req);
  void callCurrent(request req);

  void writeCPUParticlesValues(MPI_File thefile, double rmin[3], double rmax[3], SPECIE* spec);
  void writeCPUParticlesValuesSingleFile(std::string  fileName, double rmin[3], double rmax[3], SPECIE* spec);
  void writeCPUParticlesValues(MPI_File thefile, SPECIE* spec, bool flagMarker);
  void writeCPUParticlesValuesSingleFile(std::string  fileName, SPECIE* spec, bool flagMarker);
  void writeCPUParticlesValues(MPI_File thefile, SPECIE* spec);
  void writeAllCPUParticlesValues(MPI_File thefile, SPECIE* spec, int maxNfloatLoc);
  void writeCPUParticlesValuesSingleFile(std::string  fileName, SPECIE* spec);
  int packageSize(int bufsize, int* groupProcNumData, int procID, int packageNumber);
  int numPackages(int bufsize, int* groupProcNumData, int procID);
  void fillRequestList(int bufsize, int* groupProcNumData, int groupNproc, std::vector<reqOutput> &reqList);
  void writeAllSeparateFilesParticlesValues(std::string fileName, SPECIE* spec);
  void writeCPUParticlesValuesWritingGroups(std::string fileName, SPECIE* spec);
  void writeCPUParticlesValuesFewFilesWritingGroups(std::string fileName, SPECIE* spec, int NParticleToWrite, MPI_Comm outputCommunicator);

  void writeSpecPhaseSpace(std::string fileName, request req);
  void writeSpecPhaseSpaceSubDomain(std::string fileName, request req);
  void callSpecPhaseSpace(request req);


  void callDiag(request req);

  int findLeftNeightbourPoint(double val, double* coords, int numcoords);
  int findRightNeightbourPoint(double val, double* coords, int numcoords);

  void prepareIntegerBigHeader(int *itodo, int uniqueN[3], int slice_rNproc[3], int Ncomp);
  void prepareFloatCoordinatesHeader(float *fcir[3], int uniqueN[3], int imin[3]);
  void writeBigHeader(MPI_File thefile, int uniqueN[3], int imin[3], int slice_rNproc[3], int Ncomp);
  void prepareIntegerSmallHeader(int *itodo, int uniqueLocN[3], int imin[3], int remains[3]);
  void writeSmallHeader(MPI_File thefile, int uniqueLocN[3], int imin[3], int remains[3]);
  void prepareFloatField(float *todo, int uniqueLocN[3], int origin[3], request req);

  void prepareCPUFieldValues(float *buffer, int uniqueLocN[], int imin[], int locimin[], int remains[3], request req);

  void findDispForSetView(MPI_Offset *disp, int myOutputID, int *totUniquePoints, int big_header, int small_header, int Ncomp);
  void findDispForSetView(MPI_Offset *disp, int myOutputID, int *bufferSize);
  void setLocalOutputOffset(int *origin, int locimin[3], int ri[3], int remains[3]);
  void writeCPUFieldValues(MPI_File thefile, int uniqueLocN[3], int locimin[3], int remains[3], request req);
  int findNumberOfParticlesInSubdomain(request req);
  int findNumberOfParticlesInSubdomainAndReorder(request req);

  void writeBigHeaderSingleFile(std::string  fileName, int uniqueN[3], int imin[3], int slice_rNproc[3], int Ncomp);
  void writeSmallHeaderSingleFile(std::string  fileName, int uniqueLocN[3], int imin[3], int remains[3]);
  void writeCPUFieldValuesSingleFile(std::string  fileName, int uniqueLocN[3], int locimin[3], int remains[3], request req);

};

#endif

