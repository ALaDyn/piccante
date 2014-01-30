#ifndef _DIAG_MANAGER_H_
#define _DIAG_MANAGER_H_

#define _USE_MATH_DEFINES
#define _CRT_SECURE_NO_WARNINGS

#include <mpi.h>
#include <iomanip>
#if defined(_MSC_VER)
#include <cstdint>
#else
#include <stdint.h>
#endif
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include "current.h"
#include "em_field.h"
#include "particle_species.h"
#include "commons.h"
#include "grid.h"
#include "structures.h"

using namespace std;


enum diagType{
	DIAG_EM_FIELD_BINARY, DIAG_EM_FIELD, DIAG_CURRENT,
	DIAG_SPEC_DENSITY, DIAG_SPEC_PHASE_SPACE, DIAG_SPEC_PHASE_SPACE_BINARY,
	DIAG_STAT, DIAG_STAT_SINGLE, DIAG_SPEC_STAT,
	DIAG_SPEC_EXTREMS, DIAG_EM_FIELD_STAT, DIAG_SPEC_DENSITY_BINARY,
	DIAG_CURRENT_BINARY, DIAG_EMJPROBE, DIAG_1D_EM, DIAG_1D_SPEC
};

struct diagEntry{
	diagType type;
	void* diagPointer;
	int start;
	int frequency;
	bool multiple;
	int tstepout;
	int lastTimeStep;
	void* auxPointer;

	std::string myOutputFile[2];
	double position[3];
	std::string auxName;
};

#define MAX_PROC_NUM_DIAG 4
#define SIZE_TIME_DIAG 7
#define SPEC_DIAG_COMP 14
#define FIELD_DIAG_COMP 14
#define FIELD_DIAG_MIN 6

//RISOLVERE CASI IN CUI (es. con atTime) DEVO APRIRE UN FILE GIA' APERTO...
//IMPLEMENTARE CHIUSURA DEI FILES BENE

class DIAG_MANAGER{
public:
	bool amIInit;
	DIAG_MANAGER();
	~DIAG_MANAGER();

	void addEMFieldBinary(int start, int frequency, EM_FIELD* field_pointer);
	void addEMFieldBinaryAtTime(int tstepout, EM_FIELD* field_pointer);

	void addEMField(int start, int frequency, EM_FIELD* field_pointer);
	void addEMFieldAtTime(int tstepout, EM_FIELD* field_pointer);

	void addEMFieldStat(int start, int frequency, EM_FIELD* field_pointer);
	void addEMFieldStatAtTime(int tstepout, EM_FIELD* field_pointer);

	void addCurrent(int start, int frequency, CURRENT* current_pointer);
	void addCurrentAtTime(int tstepout, CURRENT* current_pointer);

	void addSpecDensity(int start, int frequency, SPECIE* spec, CURRENT* aux);
	void addSpecDensityAtTime(int tstepout, SPECIE* spec, CURRENT* aux);

	void addSpecPhaseSpace(int start, int frequency, SPECIE* spec);
	void addSpecPhaseSpaceAtTime(int tstepout, SPECIE* spec);

	void addSpecPhaseSpaceBinary(int start, int frequency, SPECIE* spec);
	void addSpecPhaseSpaceBinaryAtTime(int tstepout, SPECIE* spec);

	void addSpecStat(int start, int frequency, SPECIE* spec);
	void addSpecStatAtTime(int tstepout, SPECIE* spec);

	void addStat(int start, int frequency, EM_FIELD* aux, std::vector<SPECIE*>* spec_vec, bool simple);
	void addStatAtTime(int tstepout, EM_FIELD* aux, std::vector<SPECIE*>* spec_vec, bool simple);

	void addStatSingle(int start, int frequency, EM_FIELD* aux, std::vector<SPECIE*>* spec_vec);
	void addStatSingleAtTime(int tstepout, EM_FIELD* aux, std::vector<SPECIE*>* spec_vec);

	void addSpecDensityBinary(int start, int frequency, SPECIE* spec, CURRENT* aux);
	void addSpecDensityBinaryAtTime(int tstepout, SPECIE* spec, CURRENT* aux);

	void addCurrentBinary(int start, int frequency, CURRENT* current_pointer);
	void addCurrentBinaryAtTime(int tstepout, CURRENT* current_pointer);

	void addEMJProbe(int start, int frequency, double x, double y, double z, EM_FIELD* field_pointer, CURRENT* current_pointer);
	void addEMJProbeAtTime(int tstepout, double x, double y, double z, EM_FIELD* field_pointer, CURRENT* current_pointer);

	void add1DEMField(int start, int frequency, EM_FIELD* field_pointer);
	void add1DEMFieldAtTime(int tstepout, EM_FIELD* field_pointer);

	void add1DSpecDen(int start, int frequency, SPECIE* spec, CURRENT* curaux);
	void add1DSpecDenAtTime(int tstepout, SPECIE* spec, CURRENT* curaux);

	void initialize(std::string ioutdir, GRID* igrid);
	void callDiags(int tstep);
	void close();
private:
	std::vector<diagEntry> diagList;
	std::ofstream outStat;
	std::ofstream outExtr;
	std::ofstream outStatSingle;
	std::ofstream outDiagIndex;
	std::string outDir;
	std::ofstream tempOut;
	bool isThereDiag;
	bool isThereDiagSingle;
	bool simpleExtrems;
	std::string outFileString;
	std::string outExtrString;
	GRID* mygrid;
	std::vector<SPECIE*> vsp;
	static const int diagWideWidth = 16;
	static const int diagNarrowWidth = 6;
	static const int diagWidth = 10;
	ACCESSO acc;

	void callEMFieldBinary(int istep, EM_FIELD* fieldp);
	void callEMField(int istep, EM_FIELD* fieldp);
	void callCurrent(int istep, CURRENT* cur);
	void callSpecDen(int istep, SPECIE* spec, CURRENT* curaux);

	void callSpecPhaseSpace(int istep, SPECIE* spec);

	void callSpecPhaseSpaceBinary(int istep, SPECIE* spec);

	void callStat(int istep, EM_FIELD* fieldp, std::vector<SPECIE*>* specp);
	void callStatSingle(int istep, EM_FIELD* fieldp, std::vector<SPECIE*>* specp);
	void callEMDiag(int istep, EM_FIELD* fieldp, std::string myOutputFile);
	void callEMExtrems(int istep, EM_FIELD* fieldp, std::string myOutputFile);
	void callSpecDiag(int istep, SPECIE* spec, std::string myOutputFile);
	void callSpecExtrems(int istep, SPECIE* spec, std::string myOutputFile);
	void callSpecStat(int istep, SPECIE* spec, std::string fileDiag, std::string fileExtrem);
	void callSpecDensityBinary(int istep, SPECIE* spec, CURRENT* curaux);

	void callCurrentBinary(int tstep, CURRENT* cur);

	void interpEB(double pos[3], EM_FIELD* ebfield, double E[3], double B[3]);
	void interpJ(double pos[3], CURRENT* curp, double J[3]);

	void callEMJProbe(int tstep, double position[3], EM_FIELD* fieldp,
		CURRENT* curp, std::string fileName);
};


#endif /* DIAG_MANAGER_H_ */
