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

#ifndef _OUTPUT_MANAGER_H_
#define _OUTPUT_MANAGER_H_

#define _USE_MATH_DEFINES
#define _CRT_SECURE_NO_WARNINGS

#include "stdlib.h"
#include <mpi.h>
#include <iomanip>
#if defined(_MSC_VER)
#include <cstdint>
#else
#include <stdint.h>
#endif
#include <vector>
#include <list>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
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

using namespace std;

#define OUTPUT_SIZE_TIME_DIAG 7

#define SPEC_DIAG_COMP 14
#define FIELD_DIAG_COMP 14

enum diagType{
    OUT_E_FIELD,
    OUT_B_FIELD,
    OUT_SPEC_PHASE_SPACE,
	OUT_DIAG,
    OUT_CURRENT,
    OUT_EB_PROBE,
    OUT_SPEC_DENSITY
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

struct outDomain{
    double coordinates[3];
    bool remainingCoord[3], subselection;
    double rmin[3], rmax[3];
    std::string name;
    outDomain();
    bool compareDomains(outDomain* rhs);
    void setFreeDimensions(bool flagX, bool flagY, bool flagZ);
    void setPointCoordinate(double X, double Y, double Z);
    void setName(std::string iname);
    void setXRange(double min, double max);
    void setYRange(double min, double max);
    void setZRange(double min, double max);


};

bool requestCompTime(const request &first, const request &second);
bool requestCompUnique(const request &first, const request &second);

class OUTPUT_MANAGER
{

public:
	OUTPUT_MANAGER(GRID* _mygrid, EM_FIELD* _myfield, CURRENT* _mycurrent, std::vector<SPECIE*> _myspecies);
	~OUTPUT_MANAGER();

	void initialize(std::string _outputDir);
	void close();

    void addEBFieldFrom(double startTime, double frequency);
    void addEBFieldAt(double atTime);
    void addEBFieldFromTo(double startTime, double frequency, double endTime);

    void addEBFieldFrom(outDomain* _domain,double startTime, double frequency);
    void addEBFieldAt(outDomain* _domain,double atTime);
    void addEBFieldFromTo(outDomain* _domain,double startTime, double frequency, double endTime);

    void addEFieldFrom(double startTime, double frequency);
    void addEFieldAt(double atTime);
    void addEFieldFromTo(double startTime, double frequency, double endTime);

    void addBFieldFrom(double startTime, double frequency);
    void addBFieldAt(double atTime);
    void addBFieldFromTo(double startTime, double frequency, double endTime);

    void addEFieldFrom(outDomain* _domain,double startTime, double frequency);
    void addEFieldAt(outDomain* _domain,double atTime);
    void addEFieldFromTo(outDomain* _domain,double startTime, double frequency, double endTime);

    void addBFieldFrom(outDomain* _domain,double startTime, double frequency);
    void addBFieldAt(outDomain* _domain,double atTime);
    void addBFieldFromTo(outDomain* Plane,double startTime, double frequency, double endTime);

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

    void addSpeciesPhaseSpaceFrom(std::string name, double startTime, double frequency);
    void addSpeciesPhaseSpaceAt(std::string name, double atTime);
    void addSpeciesPhaseSpaceFromTo(std::string name, double startTime, double frequency, double endTime);

	void addDiagFrom(double startTime, double frequency);
	void addDiagAt(double atTime);
	void addDiagFromTo(double startTime, double frequency, double endTime);


	void autoVisualDiag();

	void callDiags(int istep);


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

	bool checkGrid();
	bool checkEMField();
	bool checkCurrent();
	bool checkSpecies();
    bool isInMyDomain(double *rr);
    bool amIInTheSubDomain(request req);
    void nearestInt(double *rr, int *ri, int *globalri);
    void findIntLocalBoundaries(double *rmin, double *rmax, int *imin, int *imax);
    void findIntGlobalBoundaries(double *rmin, double *rmax, int *imin, int *imax);
    void findNumberOfProc(int *Nproc, int *imin, int *imax);
    int findSpecName(std::string name);
    int returnDomainIfProbeIsInList(emProbe *newProbe);
    int returnDomainIDIfDomainIsInList(outDomain *newDomain);
    std::string diagFileName;
	std::string extremaFieldFileName;
	std::vector<std::string> extremaSpecFileNames;

	std::list<request> requestList;
	std::vector<int> timeList;
	std::map< int, std::vector<request> > allOutputs;

	static const int diagWideWidth = 16;
	static const int diagNarrowWidth = 6;
	static const int diagWidth = 10;

	int getIntTime(double dtime);

    void addRequestToList(std::list<request>& timeList, diagType type, int target,  int domain, double startTime, double frequency, double endTime);

	void prepareOutputMap();

	void createDiagFile();
	void createExtremaFiles();
    void createEMProbeFiles();

	void processOutputEntry(request req);

    std::string composeOutputName(std::string dir, std::string out, std::string opt, double time, std::string ext);
    std::string composeOutputName(std::string dir, std::string out, std::string opt1, std::string opt2, int domain, double time, std::string ext);
    void writeEMFieldMap(std::ofstream &output, request req);
	void writeEMFieldBinary(std::string fileName, request req);
    void writeNewEMFieldBinary(std::string fileName, request req, int comp);
    void writeEMFieldBinaryHDF5(std::string fileName, request req);
    void callEMFieldOld(request req);

    void callEMFieldProbe(request req);
    void writeEBFieldDomain(std::string fileName,  request req);
    void writeEBFieldSubDomain(std::string fileName,  request req);
    void callEMFieldDomain(request req);

    void writeSpecDensityMap(std::ofstream &output, request req);
    void writeSpecDensityOld(std::string fileName, request req);
    void writeSpecDensity(std::string fileName, request req);
    void writeSpecDensitySubDomain(std::string fileName,  request req);
    void callSpecDensity(request req);

	void writeCurrentMap(std::ofstream &output, request req);
    void writeCurrentOld(std::string fileName, request req);
    void writeCurrentNew(std::string fileName, request req);
    void callCurrent(request req);

    void writeSpecPhaseSpace(std::string fileName, request req);
    void callSpecPhaseSpace(request req);


	void callDiag(request req);
    void interpEB(double pos[3], double E[3], double B[3]);
};

#endif // DIAG_MANAGER_PLUS_H

