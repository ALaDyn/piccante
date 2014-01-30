#include"diag_manager.h"


int is_big_endian(void)
{
	union {
		uint32_t i;
		char c[4];
	} bint = { 0x01020304 };

	return bint.c[0] == 1;
}

DIAG_MANAGER::DIAG_MANAGER(){
	isThereDiag = false;
	isThereDiagSingle = false;
	amIInit = false;
}

DIAG_MANAGER::~DIAG_MANAGER(){
}

void DIAG_MANAGER::addEMFieldBinary(int start, int frequency, EM_FIELD* field_pointer){
	diagEntry entry;
	entry.type = DIAG_EM_FIELD_BINARY;
	entry.start = start;
	entry.frequency = frequency;
	entry.diagPointer = field_pointer;
	entry.multiple = true;
	diagList.push_back(entry);
}

void DIAG_MANAGER::addEMFieldBinaryAtTime(int tstepout, EM_FIELD* field_pointer){
	diagEntry entry;
	entry.type = DIAG_EM_FIELD_BINARY;
	entry.tstepout = tstepout;
	entry.diagPointer = field_pointer;
	entry.multiple = false;
	diagList.push_back(entry);
}

void DIAG_MANAGER::addEMField(int start, int frequency, EM_FIELD* field_pointer){
	diagEntry entry;
	entry.type = DIAG_EM_FIELD;
	entry.start = start;
	entry.frequency = frequency;
	entry.diagPointer = field_pointer;
	entry.multiple = true;
	diagList.push_back(entry);
}

void DIAG_MANAGER::addEMFieldAtTime(int tstepout, EM_FIELD* field_pointer){
	diagEntry entry;
	entry.type = DIAG_EM_FIELD;
	entry.tstepout = tstepout;
	entry.diagPointer = field_pointer;
	entry.multiple = false;
	diagList.push_back(entry);
}

void DIAG_MANAGER::addEMFieldStat(int start, int frequency, EM_FIELD* field_pointer){
	diagEntry entry;
	entry.type = DIAG_EM_FIELD_STAT;
	entry.start = start;
	entry.frequency = frequency;
	entry.diagPointer = field_pointer;
	entry.multiple = true;
	//entry.myOutputFile=NULL;
	diagList.push_back(entry);
}

void DIAG_MANAGER::addEMFieldStatAtTime(int tstepout, EM_FIELD* field_pointer){
	diagEntry entry;
	entry.type = DIAG_EM_FIELD_STAT;
	entry.tstepout = tstepout;
	entry.diagPointer = field_pointer;
	entry.multiple = false;
	//entry.myOutputFile=NULL;
	diagList.push_back(entry);
}

void DIAG_MANAGER::addCurrent(int start, int frequency, CURRENT* current_pointer){
	diagEntry entry;
	entry.type = DIAG_CURRENT;
	entry.start = start;
	entry.frequency = frequency;
	entry.diagPointer = current_pointer;
	entry.multiple = true;
	diagList.push_back(entry);
}

void DIAG_MANAGER::addCurrentAtTime(int tstepout, CURRENT* current_pointer){
	diagEntry entry;
	entry.type = DIAG_CURRENT;
	entry.tstepout = tstepout;
	entry.diagPointer = current_pointer;
	entry.multiple = false;
	diagList.push_back(entry);
}

void DIAG_MANAGER::addSpecDensity(int start, int frequency, SPECIE* spec, CURRENT* aux){
	diagEntry entry;
	entry.type = DIAG_SPEC_DENSITY;
	entry.start = start;
	entry.frequency = frequency;
	entry.diagPointer = spec;
	entry.auxPointer = aux;
	entry.multiple = true;
	diagList.push_back(entry);
}

void DIAG_MANAGER::addSpecDensityAtTime(int tstepout, SPECIE* spec, CURRENT* aux){
	diagEntry entry;
	entry.type = DIAG_SPEC_DENSITY;
	entry.tstepout = tstepout;
	entry.diagPointer = spec;
	entry.auxPointer = aux;
	entry.multiple = false;
	diagList.push_back(entry);
}

void DIAG_MANAGER::addSpecPhaseSpace(int start, int frequency, SPECIE* spec){
	diagEntry entry;
	entry.type = DIAG_SPEC_PHASE_SPACE;
	entry.start = start;
	entry.frequency = frequency;
	entry.diagPointer = spec;
	entry.multiple = true;
	diagList.push_back(entry);
}

void DIAG_MANAGER::addSpecPhaseSpaceAtTime(int tstepout, SPECIE* spec){
	diagEntry entry;
	entry.type = DIAG_SPEC_PHASE_SPACE;
	entry.tstepout = tstepout;
	entry.diagPointer = spec;
	entry.multiple = false;
	diagList.push_back(entry);
}

void DIAG_MANAGER::addSpecPhaseSpaceBinary(int start, int frequency, SPECIE* spec){
	diagEntry entry;
	entry.type = DIAG_SPEC_PHASE_SPACE_BINARY;
	entry.start = start;
	entry.frequency = frequency;
	entry.diagPointer = spec;
	entry.multiple = true;
	diagList.push_back(entry);
}

void DIAG_MANAGER::addSpecPhaseSpaceBinaryAtTime(int tstepout, SPECIE* spec){
	diagEntry entry;
	entry.type = DIAG_SPEC_PHASE_SPACE_BINARY;
	entry.tstepout = tstepout;
	entry.diagPointer = spec;
	entry.multiple = false;
	diagList.push_back(entry);
}

void DIAG_MANAGER::addSpecStat(int start, int frequency, SPECIE* spec){
	diagEntry entry;
	entry.type = DIAG_SPEC_STAT;
	entry.start = start;
	entry.frequency = frequency;
	entry.diagPointer = spec;
	entry.multiple = true;
	//	entry.myOutputFile=NULL;
	diagList.push_back(entry);
}

void DIAG_MANAGER::addSpecStatAtTime(int tstepout, SPECIE* spec){
	diagEntry entry;
	entry.type = DIAG_SPEC_STAT;
	entry.tstepout = tstepout;
	entry.diagPointer = spec;
	entry.multiple = false;
	//	entry.myOutputFile=NULL;
	diagList.push_back(entry);
}

void DIAG_MANAGER::addStat(int start, int frequency, EM_FIELD* aux, std::vector<SPECIE*>* spec_vec, bool simple){
	isThereDiag = true;
	diagEntry entry;
	entry.type = DIAG_STAT;
	entry.start = start;
	entry.frequency = frequency;
	entry.diagPointer = spec_vec;
	entry.auxPointer = aux;
	entry.multiple = true;
	simpleExtrems = simple;
	diagList.push_back(entry);
	vsp = *spec_vec;
}

void DIAG_MANAGER::addStatAtTime(int tstepout, EM_FIELD* aux, std::vector<SPECIE*>* spec_vec, bool simple){
	isThereDiag = true;
	diagEntry entry;
	entry.type = DIAG_STAT;
	entry.tstepout = tstepout;
	entry.diagPointer = spec_vec;
	entry.auxPointer = aux;
	entry.multiple = false;
	simpleExtrems = simple;
	diagList.push_back(entry);
}

void DIAG_MANAGER::addStatSingle(int start, int frequency, EM_FIELD* aux, std::vector<SPECIE*>* spec_vec){
	isThereDiagSingle = true;
	diagEntry entry;
	entry.type = DIAG_STAT_SINGLE;
	entry.start = start;
	entry.frequency = frequency;
	entry.diagPointer = spec_vec;
	entry.auxPointer = aux;
	entry.multiple = true;
	diagList.push_back(entry);
	vsp = *spec_vec;
}

void DIAG_MANAGER::addStatSingleAtTime(int tstepout, EM_FIELD* aux, std::vector<SPECIE*>* spec_vec){
	isThereDiagSingle = true;
	diagEntry entry;
	entry.type = DIAG_STAT_SINGLE;
	entry.tstepout = tstepout;
	entry.diagPointer = spec_vec;
	entry.auxPointer = aux;
	entry.multiple = false;
	diagList.push_back(entry);
}

void DIAG_MANAGER::addSpecDensityBinary(int start, int frequency, SPECIE* spec, CURRENT* aux){
	diagEntry entry;
	entry.type = DIAG_SPEC_DENSITY_BINARY;
	entry.start = start;
	entry.frequency = frequency;
	entry.diagPointer = spec;
	entry.auxPointer = aux;
	entry.multiple = true;
	diagList.push_back(entry);
}

void DIAG_MANAGER::addSpecDensityBinaryAtTime(int tstepout, SPECIE* spec, CURRENT* aux){
	diagEntry entry;
	entry.type = DIAG_SPEC_DENSITY_BINARY;
	entry.tstepout = tstepout;
	entry.diagPointer = spec;
	entry.auxPointer = aux;
	entry.multiple = false;
	diagList.push_back(entry);
}


void DIAG_MANAGER::addCurrentBinary(int start, int frequency, CURRENT* current_pointer){
	diagEntry entry;
	entry.type = DIAG_CURRENT_BINARY;
	entry.start = start;
	entry.frequency = frequency;
	entry.diagPointer = current_pointer;
	entry.multiple = true;
	diagList.push_back(entry);
}

void DIAG_MANAGER::addCurrentBinaryAtTime(int tstepout, CURRENT* current_pointer){
	diagEntry entry;
	entry.type = DIAG_CURRENT_BINARY;
	entry.tstepout = tstepout;
	entry.diagPointer = current_pointer;
	entry.multiple = false;
	diagList.push_back(entry);
}

void DIAG_MANAGER::addEMJProbe(int start, int frequency, double x, double y, double z, EM_FIELD* field_pointer, CURRENT* current_pointer){
	diagEntry entry;
	entry.type = DIAG_EMJPROBE;
	entry.diagPointer = field_pointer;
	entry.auxPointer = current_pointer;
	entry.position[0] = x;
	entry.position[1] = y;
	entry.position[2] = z;
	entry.start = start;
	entry.frequency = frequency;
	entry.multiple = true;
	diagList.push_back(entry);
}

void DIAG_MANAGER::addEMJProbeAtTime(int tstepout, double x, double y, double z, EM_FIELD* field_pointer, CURRENT* current_pointer){
	diagEntry entry;
	entry.type = DIAG_EMJPROBE;
	entry.diagPointer = field_pointer;
	entry.auxPointer = current_pointer;
	entry.position[0] = x;
	entry.position[1] = y;
	entry.position[2] = z;
	entry.tstepout = tstepout;
	entry.multiple = false;
	diagList.push_back(entry);
}

void DIAG_MANAGER::add1DSpecDen(int start, int frequency, SPECIE* spec, CURRENT* aux){
	diagEntry entry;
	entry.type = DIAG_1D_SPEC;
	entry.start = start;
	entry.frequency = frequency;
	entry.diagPointer = spec;
	entry.auxPointer = aux;
	entry.multiple = true;
	diagList.push_back(entry);
}

void DIAG_MANAGER::add1DSpecDenAtTime(int tstepout, SPECIE* spec, CURRENT* aux){
	diagEntry entry;
	entry.type = DIAG_1D_SPEC;
	entry.tstepout = tstepout;
	entry.diagPointer = spec;
	entry.auxPointer = aux;
	entry.multiple = false;
	diagList.push_back(entry);
}

void DIAG_MANAGER::add1DEMField(int start, int frequency, EM_FIELD* field_pointer){
	diagEntry entry;
	entry.type = DIAG_1D_EM;
	entry.start = start;
	entry.frequency = frequency;
	entry.diagPointer = field_pointer;
	entry.multiple = true;
	diagList.push_back(entry);
}

void DIAG_MANAGER::add1DEMFieldAtTime(int tstepout, EM_FIELD* field_pointer){
	diagEntry entry;
	entry.type = DIAG_1D_EM;
	entry.tstepout = tstepout;
	entry.diagPointer = field_pointer;
	entry.multiple = false;
	diagList.push_back(entry);
}


void DIAG_MANAGER::initialize(std::string ioutdir, GRID* igrid){
	mygrid = igrid;
	outDir = ioutdir;

	std::stringstream ssIndex;
	ssIndex << outDir << "/SIMULATION.index";

	if (mygrid->myid == mygrid->master_proc){
		outDiagIndex.open(ssIndex.str().c_str(), std::ofstream::out | std::ofstream::trunc);
		outDiagIndex << "!===========================================================================" << std::endl;
		outDiagIndex << "!Lines with '!' are comments. Lines with '+' are additional infos for reader" << std::endl;
		outDiagIndex << "!Only collective output is listed in this file." << std::endl;
		outDiagIndex << "!===========================================================================" << std::endl;
		outDiagIndex << "!" << std::endl;
		outDiagIndex << "+DIAG_INDEX_START" << std::endl;
		outDiagIndex.flush();
	}


	if (isThereDiag){
		if (mygrid->myid == mygrid->master_proc){
			std::stringstream ss;
			ss << outDir << "/diag.txt";
			outFileString = ss.str();
			outStat.open(outFileString.c_str(), std::ofstream::out | std::ofstream::trunc);

			outDiagIndex << "MAIN_DIAG" << "\t" << "-1" << "\t" << outFileString << std::endl;

			outStat << " " << setw(diagNarrowWidth) << "#istep" << " " << setw(diagWidth) << "time";
			outStat << " " << setw(diagWidth) << "Etot";
			outStat << " " << setw(diagWidth) << "Ex2" << " " << setw(diagWidth) << "Ey2" << " " << setw(diagWidth) << "Ez2";
			outStat << " " << setw(diagWidth) << "Bx2" << " " << setw(diagWidth) << "By2" << " " << setw(diagWidth) << "Bz2";

			std::vector<SPECIE*>::const_iterator spec_iterator;
			outDiagIndex << "+SPECS" << "\t";
			outDiagIndex << vsp.size();

			for (spec_iterator = vsp.begin(); spec_iterator != vsp.end(); spec_iterator++){
				outDiagIndex << "\t" << (*spec_iterator)->name;
				outStat << " " << setw(diagWidth) << (*spec_iterator)->name;
			}

			outDiagIndex << std::endl;

			outStat << std::endl;

			std::stringstream ssB;
			ssB << outDir << "/extrema.txt";
			outExtrString = ssB.str();
			outExtr.open(outExtrString.c_str(), std::ofstream::out | std::ofstream::trunc);

			outDiagIndex << "EXTREMA" << "\t" << "-1" << "\t" << outExtrString << std::endl;

			if (simpleExtrems){
				outExtr << " " << setw(diagNarrowWidth) << "#istep" << " " << setw(diagWidth) << "time";
				outExtr << " " << setw(diagWidth) << "Emax" << " " << setw(diagWidth) << "Bmax";

				for (spec_iterator = vsp.begin(); spec_iterator != vsp.end(); spec_iterator++){
					outExtr << " " << setw(diagNarrowWidth) << (*spec_iterator)->name;
					outExtr << setw(diagNarrowWidth) << " Kmax";
				}
			}
			else{
				outExtr << "#STEP\t\tt\t\tExmin\t\tExmax\t\tEymin\t\tEymax\t\tEzmin\t\tEzmax\t\tBxmin\t\tBxmax\t\tBymin\t\tBymax\t\tBzmin\t\tBzmax\t\tEtmax\t\tBtmax";

				for (spec_iterator = vsp.begin(); spec_iterator != vsp.end(); spec_iterator++){
					outExtr << "\t\t" << (*spec_iterator)->name << "\t\txmin\t\txmax\t\tymin\t\tymax\t\tzmin\t\tzmax\t\tpxmin\t\tpxmax\t\tpymin\t\tpymax\t\tpzmin\t\tpzmax\t\temin\t\temax";
				}
			}
			outExtr << std::endl;
		}
	}

	if (isThereDiagSingle){
		std::stringstream ss;
		ss << outDir << "/p" << setfill('0') << setw(MAX_PROC_NUM_DIAG) << mygrid->myid << "_diag.txt";
		outFileString = ss.str();
		outStatSingle.open(outFileString.c_str(), std::ofstream::out | std::ofstream::trunc);

		outStatSingle << "#STEP    t     Etot    Ex     Ey     Ez     Bx    By     Bz";

		std::vector<SPECIE*>::const_iterator spec_iterator;
		for (spec_iterator = vsp.begin(); spec_iterator != vsp.end(); spec_iterator++){
			outStatSingle << "     " << (*spec_iterator)->name;
		}

		outStatSingle << std::endl;
	}

	for (size_t i = 0; i < diagList.size(); i++){
		diagList[i].lastTimeStep = -1;

		if (diagList[i].type == DIAG_EM_FIELD_STAT){
			if (mygrid->myid == mygrid->master_proc){
				std::stringstream myname0;
				myname0 << outDir << "/EMfield_diag.txt";

				outDiagIndex << "EM_DIAG" << "\t" << "-1" << "\t" << myname0.str() << std::endl;

				diagList[i].myOutputFile[0] = myname0.str();
				std::ofstream of1;
				of1.open(diagList[i].myOutputFile[0].c_str());

				((EM_FIELD*)diagList[i].diagPointer)->init_output_diag(of1);

				of1.close();

				std::ofstream of2;

				std::stringstream myname1;
				myname1 << outDir << "/EMfield_extrems.txt";

				outDiagIndex << "EM_EXTRM" << "\t" << "-1" << "\t" << myname1.str() << std::endl;

				diagList[i].myOutputFile[1] = myname1.str();
				of2.open(diagList[i].myOutputFile[1].c_str());

				((EM_FIELD*)diagList[i].diagPointer)->init_output_extrems(of2);
				of2.close();
			}
		}

		if (diagList[i].type == DIAG_SPEC_STAT){
			if (mygrid->myid == mygrid->master_proc){
				std::stringstream myname0;
				myname0 << outDir << "/SP_" << ((SPECIE*)diagList[i].diagPointer)->name << "_diag.txt";

				outDiagIndex << "SP_DIAG" << "\t" << ((SPECIE*)diagList[i].diagPointer)->name << "\t" << "-1" << "\t" << myname0.str() << std::endl;

				std::ofstream of1;
				std::ofstream of2;

				diagList[i].myOutputFile[0] = myname0.str();

				of1.open(diagList[i].myOutputFile[0].c_str());
				((SPECIE*)diagList[i].diagPointer)->init_output_diag(of1);
				of1.close();


				std::stringstream myname1;
				myname1 << outDir << "/SP_" << ((SPECIE*)diagList[i].diagPointer)->name << "_extrems.txt";

				outDiagIndex << "SP_EXTRM" << "\t" << ((SPECIE*)diagList[i].diagPointer)->name << "\t" << "-1" << "\t" << myname1.str() << std::endl;

				diagList[i].myOutputFile[1] = myname1.str();
				of2.open(diagList[i].myOutputFile[1].c_str());

				((SPECIE*)diagList[i].diagPointer)->init_output_extrems(of2);

				of2.close();
			}
		}

		if (diagList[i].type == DIAG_EMJPROBE){

			std::stringstream myname;
			myname << outDir << "/EMJ_"
				<< std::fixed << std::setprecision(5) << diagList[i].position[0] << "-"
				<< std::fixed << std::setprecision(5) << diagList[i].position[1] << "-"
				<< std::fixed << std::setprecision(5) << diagList[i].position[2]
				<< ".prb";

			diagList[i].auxName = myname.str();

			if (mygrid->myid == mygrid->master_proc){
				std::ofstream outf;
				outf.open(myname.str().c_str(), std::ios_base::trunc);
				outDiagIndex << "EMJ" << "\t" << "-1" << "\t" << myname.str() << std::endl;
				outf << "# EMJ probe at position ( "
					<< diagList[i].position[0] << " , "
					<< diagList[i].position[1] << " , "
					<< diagList[i].position[2] << " )" << std::endl
					<< "#E,B calculated with particle shapefunction" << std::endl
					<< "#J calculated at nearest integer point" << std::endl
					<< std::endl
					<< "#" << setw(diagNarrowWidth) << "tstep"
					<< " " << setw(diagWidth) << "time"
					<< " " << setw(diagWidth) << "Ex"
					<< " " << setw(diagWidth) << "Ey"
					<< " " << setw(diagWidth) << "Ez"
					<< " " << setw(diagWidth) << "Bx"
					<< " " << setw(diagWidth) << "By"
					<< " " << setw(diagWidth) << "Bz"
					<< " " << setw(diagWidth) << "Jx"
					<< " " << setw(diagWidth) << "Jy"
					<< " " << setw(diagWidth) << "Jz"
					<< std::endl;
				outf.close();
			}
		}
	}


	amIInit = true;
}

void DIAG_MANAGER::callDiags(int tstep){

	for (size_t i = 0; i < diagList.size(); i++){
		if (((!diagList[i].multiple) && (tstep == diagList[i].tstepout)) ||
			(diagList[i].multiple && (tstep >= diagList[i].start) && ((tstep - diagList[i].start) % diagList[i].frequency == 0)))

		{
			if (diagList[i].lastTimeStep != tstep){
				diagList[i].lastTimeStep = tstep;
				switch (diagList[i].type){
				case DIAG_EM_FIELD:
					callEMField(tstep, (EM_FIELD*)diagList[i].diagPointer);
					break;

				case DIAG_EM_FIELD_BINARY:
					callEMFieldBinary(tstep, (EM_FIELD*)diagList[i].diagPointer);
					break;

				case DIAG_CURRENT:
					callCurrent(tstep, (CURRENT*)diagList[i].diagPointer);
					break;

				case DIAG_SPEC_DENSITY:
					callSpecDen(tstep, (SPECIE*)diagList[i].diagPointer, (CURRENT*)diagList[i].auxPointer);
					break;

				case DIAG_SPEC_PHASE_SPACE:
					callSpecPhaseSpace(tstep, (SPECIE*)diagList[i].diagPointer);
					break;

				case DIAG_SPEC_PHASE_SPACE_BINARY:
					callSpecPhaseSpaceBinary(tstep, (SPECIE*)diagList[i].diagPointer);
					break;

				case	DIAG_STAT:
					callStat(tstep, (EM_FIELD*)diagList[i].auxPointer, (std::vector<SPECIE*>*) diagList[i].diagPointer);
					break;

				case	DIAG_STAT_SINGLE:
					callStatSingle(tstep, (EM_FIELD*)diagList[i].auxPointer, (std::vector<SPECIE*>*) diagList[i].diagPointer);
					break;

				case	DIAG_EM_FIELD_STAT:
					callEMDiag(tstep, (EM_FIELD*)diagList[i].diagPointer, diagList[i].myOutputFile[0]);
					callEMExtrems(tstep, (EM_FIELD*)diagList[i].diagPointer, diagList[i].myOutputFile[1]);
					break;

				case	DIAG_SPEC_STAT:
					//callSpecDiag(tstep, (SPECIE*) diagList[i].diagPointer,  diagList[i].myOutputFile[0]);
					//callSpecExtrems(tstep, (SPECIE*) diagList[i].diagPointer,  diagList[i].myOutputFile[1]);
					callSpecStat(tstep, (SPECIE*)diagList[i].diagPointer, diagList[i].myOutputFile[0], diagList[i].myOutputFile[1]);
					break;

				case DIAG_SPEC_DENSITY_BINARY:
					callSpecDensityBinary(tstep, (SPECIE*)diagList[i].diagPointer, (CURRENT*)diagList[i].auxPointer);
					break;

				case DIAG_CURRENT_BINARY:
					callCurrentBinary(tstep, (CURRENT*)diagList[i].diagPointer);

					break;

				case DIAG_EMJPROBE:
					callEMJProbe(tstep, diagList[i].position,
						(EM_FIELD*)diagList[i].diagPointer,
						(CURRENT*)diagList[i].auxPointer,
						diagList[i].auxName);

					break;

				default:
					break;
				}
			}
		}
	}
}

void DIAG_MANAGER::close(){
	outDiagIndex << "+DIAG_INDEX_END" << std::endl;
	outStat.close();
	outStatSingle.close();
	outDiagIndex.close();
	outExtr.close();
}


void DIAG_MANAGER::callEMFieldBinary(int istep, EM_FIELD* fieldp)
{
	int Ncomp = fieldp->getNcomp();

	if (mygrid->myid == mygrid->master_proc){
		std::stringstream namedat;
		namedat << outDir << "/" << "EMfield_";
		namedat << std::setfill('0') << std::setw(SIZE_TIME_DIAG) << std::fixed << std::setprecision(3) << mygrid->time << ".map";

		outDiagIndex << "EM_MAP" << "\t" << istep << "\t" << std::fixed << std::setprecision(3)
			<< mygrid->time << "\t" << namedat.str() << std::endl;

		tempOut.open(namedat.str().c_str());
		int uniqueN[3];
		uniqueN[0] = mygrid->uniquePoints[0];
		uniqueN[1] = mygrid->uniquePoints[1];
		uniqueN[2] = mygrid->uniquePoints[2];

		for (int c = 0; c < 3; c++)
			tempOut << uniqueN[c] << "\t";
		tempOut << endl;
		for (int c = 0; c < 3; c++)
			tempOut << mygrid->rnproc[c] << "\t";
		tempOut << endl;

		tempOut << "Ncomp: " << Ncomp << std::endl;
		for (int c = 0; c < Ncomp; c++){
			integer_or_halfinteger crd = fieldp->getCompCoords(c);
			tempOut << c << ":\t" << (int)crd.x << "\t"
				<< (int)crd.y << "\t"
				<< (int)crd.z << std::endl;
		}

		for (int c = 0; c < 3; c++){
			for (int i = 0; i < uniqueN[c]; i++)
				tempOut << mygrid->cir[c][i] << "  ";
			tempOut << endl;
		}
		for (int c = 0; c < 3; c++){
			for (int i = 0; i < uniqueN[c]; i++)
				tempOut << mygrid->chr[c][i] << "  ";
			tempOut << endl;
		}


		tempOut.close();


	}

	//////////////////////////// END ofserial data file write

	std::stringstream  namebin;

	namebin << outDir << "/" << "EMfield_";
	namebin << std::setfill('0') << std::setw(SIZE_TIME_DIAG) << std::fixed << std::setprecision(3) << mygrid->time << ".bin";


	if (mygrid->myid == mygrid->master_proc)
		outDiagIndex << "EM_BIN" << "\t" << istep << "\t" << std::fixed << std::setprecision(3)
		<< mygrid->time << "\t" << namebin.str() << std::endl;

	int uniqueN[3];
	uniqueN[0] = mygrid->uniquePoints[0];
	uniqueN[1] = mygrid->uniquePoints[1];
	uniqueN[2] = mygrid->uniquePoints[2];
	long disp = 0;
	int small_header = 6 * sizeof(int);
	int big_header = 3 * sizeof(int)+3 * sizeof(int)
		+2 * (uniqueN[0] + uniqueN[1] + uniqueN[2])*sizeof(double)
		+sizeof(int)+fieldp->getNcomp() * 3 * sizeof(char);

	MPI_File thefile;
	MPI_Status status;
	char *nomefile = new char[namebin.str().size() + 1];
	nomefile[namebin.str().size()] = 0;
	sprintf(nomefile, "%s", namebin.str().c_str());
	MPI_File_open(MPI_COMM_WORLD, nomefile,
		MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &thefile);
	//+++++++++++ FILE HEADER  +++++++++++++++++++++
	if (mygrid->myid == 0){
		MPI_File_write(thefile, uniqueN, 3, MPI_INT, &status);
		MPI_File_write(thefile, mygrid->rnproc, 3, MPI_INT, &status);
		MPI_File_write(thefile, &Ncomp, 1, MPI_INT, &status);
		for (int c = 0; c < Ncomp; c++){
			integer_or_halfinteger crd = fieldp->getCompCoords(c);
			char tp[3] = { crd.x, crd.y, crd.z };
			MPI_File_write(thefile, tp, 3, MPI_CHAR, &status);
		}
		for (int c = 0; c < 3; c++)
			MPI_File_write(thefile, mygrid->cir[c], uniqueN[c], MPI_DOUBLE, &status);
		for (int c = 0; c < 3; c++)
			MPI_File_write(thefile, mygrid->chr[c], uniqueN[c], MPI_DOUBLE, &status);


	}

	//*********** END HEADER *****************


	disp = big_header;
	for (int rank = 0; rank < mygrid->myid; rank++)
		disp += small_header + mygrid->proc_totUniquePoints[rank] * sizeof(float)*Ncomp;
	MPI_File_set_view(thefile, disp, MPI_FLOAT, MPI_FLOAT, (char *) "native", MPI_INFO_NULL);

	//+++++++++++ Start CPU HEADER  +++++++++++++++++++++
	{
		int todo[6];
		todo[0] = mygrid->rproc_imin[0][mygrid->rmyid[0]];
		todo[1] = mygrid->rproc_imin[1][mygrid->rmyid[1]];
		todo[2] = mygrid->rproc_imin[2][mygrid->rmyid[2]];
		todo[3] = mygrid->uniquePointsloc[0];
		todo[4] = mygrid->uniquePointsloc[1];
		todo[5] = mygrid->uniquePointsloc[2];
		MPI_File_write(thefile, todo, 6, MPI_INT, &status);
	}
	//+++++++++++ Start CPU Field Values  +++++++++++++++++++++
	{
		float *todo;
		int Nx, Ny, Nz, Ncomp;
		Nx = mygrid->uniquePointsloc[0];
		Ny = mygrid->uniquePointsloc[1];
		Nz = mygrid->uniquePointsloc[2];
		Ncomp = fieldp->getNcomp();
		int size = Ncomp*Nx*Ny*Nz;
		todo = new float[size];
		for (int k = 0; k < Nz; k++)
		for (int j = 0; j < Ny; j++)
		for (int i = 0; i < Nx; i++)
		for (int c = 0; c < Ncomp; c++)
			todo[c + i*Ncomp + j*Nx*Ncomp + k*Ny*Nx*Ncomp] = (float)fieldp->VEB(c, i, j, k);
		MPI_File_write(thefile, todo, size, MPI_FLOAT, &status);
		delete[]todo;
	}

	MPI_File_close(&thefile);
	//////////////////////////// END of collective binary file write


}
void DIAG_MANAGER::callEMField(int istep, EM_FIELD* fieldp)
{

	std::stringstream ss1;
	ss1 << outDir << "/p" << std::setfill('0') << std::setw(MAX_PROC_NUM_DIAG) << mygrid->myid << "_E_";
	ss1 << std::setfill('0') << std::setw(SIZE_TIME_DIAG) << std::fixed << std::setprecision(3) << mygrid->time;

	std::stringstream ss2;
	ss2 << outDir << "/p" << std::setfill('0') << std::setw(MAX_PROC_NUM_DIAG) << mygrid->myid << "_B_";
	ss2 << std::setfill('0') << std::setw(SIZE_TIME_DIAG) << std::fixed << std::setprecision(3) << mygrid->time;

	tempOut.open(ss1.str().c_str());
	fieldp->output_2D_E(tempOut);
	tempOut.close();

	tempOut.open(ss2.str().c_str());
	fieldp->output_2D_B(tempOut);
	tempOut.close();
}

void DIAG_MANAGER::callCurrent(int istep, CURRENT* cur){
	double tw = mygrid->dt * istep;
	std::stringstream ss1;
	ss1 << outDir << "/p" << std::setfill('0') << std::setw(MAX_PROC_NUM_DIAG) << mygrid->myid << "_J_";
	ss1 << std::setfill('0') << std::setw(SIZE_TIME_DIAG) << std::fixed << std::setprecision(3) << tw;

	tempOut.open(ss1.str().c_str());
	cur->output_2D_current(tempOut);
	tempOut.close();

}
void DIAG_MANAGER::callSpecDen(int istep, SPECIE* spec, CURRENT* curaux){
	double tw = mygrid->dt * istep;
	std::stringstream ss1;
	ss1 << outDir << "/p" << std::setfill('0') << std::setw(MAX_PROC_NUM_DIAG) << mygrid->myid << "_" << spec->name << "_den_";
	ss1 << std::setfill('0') << std::setw(SIZE_TIME_DIAG) << std::fixed << std::setprecision(3) << tw;

	tempOut.open(ss1.str().c_str());
	curaux->setAllValuesToZero();
	spec->density_deposition_standard(curaux);
	curaux->pbc();
	curaux->output_2D_density(tempOut);
	tempOut.close();
}

void DIAG_MANAGER::callSpecPhaseSpace(int istep, SPECIE* spec){
	double tw = mygrid->dt * istep;
	std::stringstream ss1;
	ss1 << outDir << "/p" << std::setfill('0') << std::setw(MAX_PROC_NUM_DIAG) << mygrid->myid << "_" << spec->name << "_pspace_";
	ss1 << std::setfill('0') << std::setw(SIZE_TIME_DIAG) << std::fixed << std::setprecision(3) << tw;

	tempOut.open(ss1.str().c_str());
	spec->output(tempOut);
	tempOut.close();

}

void DIAG_MANAGER::callSpecPhaseSpaceBinary(int istep, SPECIE* spec){
	double tw = mygrid->time;
	int* NfloatLoc = new int[mygrid->nproc];

	std::stringstream ss1, namebin, namedat;

	ss1 << outDir << "/" << "SP_" << spec->name << "_pspace_";
	ss1 << std::setfill('0') << std::setw(SIZE_TIME_DIAG) << std::fixed << std::setprecision(3) << tw;
	namebin << ss1.str() << ".bin";
	namedat << ss1.str() << ".dat";

	if (mygrid->myid == mygrid->master_proc){
		outDiagIndex << "SP_PSDAT" << "\t" << spec->name << "\t" << istep << "\t" << std::fixed << std::setprecision(3)
			<< mygrid->time << "\t" << namedat.str() << std::endl; outDiagIndex.flush();

		int iparam[20];
		float fparam[20];
		for (int i = 0; i < 20; i++){
			iparam[i] = 0;
			fparam[i] = 0;
		}
		iparam[0] = mygrid->nproc;
		iparam[17] = spec->Ncomp;
		iparam[18] = iparam[19] = is_big_endian() + 1;
		tempOut.open(namedat.str().c_str());
		tempOut << "Integer Parameters" << endl;
		for (int i = 0; i < 20; i++){
			tempOut << iparam[i] << "\t";
		}
		tempOut << endl << "Real Parameters" << endl;
		for (int i = 0; i < 20; i++){
			tempOut << fparam[i] << "\t";
		}
		tempOut.close();

		outDiagIndex << "SP_PSBIN" << "\t" << spec->name << "\t" << istep << "\t" << std::fixed << std::setprecision(3)
			<< mygrid->time << "\t" << namebin.str() << std::endl; outDiagIndex.flush();
	}




	NfloatLoc[mygrid->myid] = spec->Np*spec->Ncomp;
	MPI_Allgather(MPI_IN_PLACE, 1, MPI_INT, NfloatLoc, 1, MPI_INT, MPI_COMM_WORLD);

	long disp = 0;
	for (int pp = 0; pp < mygrid->myid; pp++)
		disp += (long)(NfloatLoc[pp] * sizeof(float));
	MPI_File thefile;
	char *nomefile = new char[namebin.str().size() + 1];
	nomefile[namebin.str().size()] = 0;
	sprintf(nomefile, "%s", namebin.str().c_str());
	MPI_File_open(MPI_COMM_WORLD, nomefile,
		MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &thefile);

	MPI_File_set_view(thefile, disp, MPI_FLOAT, MPI_FLOAT, (char *) "native", MPI_INFO_NULL);

	float *buf;
	MPI_Status status;
	int dimensione, passaggi, resto;
	dimensione = 100000;
	buf = new float[dimensione*spec->Ncomp];
	passaggi = spec->Np / dimensione;
	resto = spec->Np % dimensione;
	for (int i = 0; i < passaggi; i++){
		for (int p = 0; p < dimensione; p++){
			for (int c = 0; c < spec->Ncomp; c++){
				buf[c + p*spec->Ncomp] = spec->ru(c, p + dimensione*i);
			}
		}
		MPI_File_write(thefile, buf, dimensione*spec->Ncomp, MPI_FLOAT, &status);
	}
	for (int p = 0; p < resto; p++){
		for (int c = 0; c < spec->Ncomp; c++){
			buf[c + p*spec->Ncomp] = spec->ru(c, p + dimensione*passaggi);
		}
	}
	MPI_File_write(thefile, buf, resto*spec->Ncomp, MPI_FLOAT, &status);
	MPI_File_close(&thefile);
	delete[]buf;
}



void DIAG_MANAGER::callStat(int istep, EM_FIELD* fieldp, std::vector<SPECIE*>* specp){
	std::vector<SPECIE*> vsp = *specp;
	std::vector<SPECIE*>::const_iterator spec_iterator;
	double * ekinSpecies;
	size_t specie = 0;
	double tw = mygrid->time;

	ekinSpecies = new double[vsp.size()];

	double t_spec_extrema[SPEC_DIAG_COMP];
	double field_extrema[FIELD_DIAG_COMP];
	//double* spec_extrema = new double[SPEC_DIAG_COMP*vsp.size()];	

	double EE[3], BE[3];
	double etotFields = fieldp->computeEnergyAndExtremes(EE, BE, field_extrema);
	double etotKin = 0;
	double tKin;

	specie = 0;
	for (spec_iterator = vsp.begin(); spec_iterator != vsp.end(); spec_iterator++){
		tKin = (*spec_iterator)->getKineticEnergyWExtrems(t_spec_extrema);
		ekinSpecies[specie++] = tKin;
		etotKin += tKin;
	}
	double etot = etotKin + etotFields;
	MPI_Allreduce(MPI_IN_PLACE, &etot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, EE, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, BE, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, ekinSpecies, vsp.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	if (mygrid->myid == mygrid->master_proc){
		outStat << " " << setw(diagNarrowWidth) << istep << " " << setw(diagWidth) << tw;
		outStat << " " << setw(diagWidth) << etot;
		outStat << " " << setw(diagWidth) << EE[0] << " " << setw(diagWidth) << EE[1] << " " << setw(diagWidth) << EE[2];
		outStat << " " << setw(diagWidth) << BE[0] << " " << setw(diagWidth) << BE[1] << " " << setw(diagWidth) << BE[2];

		for (specie = 0; specie < vsp.size(); specie++){
			outStat << " " << setw(diagWidth) << ekinSpecies[specie];
		}
		outStat << std::endl;

		outExtr << " " << setw(diagNarrowWidth) << istep << " " << setw(diagWidth) << tw;
		if (simpleExtrems){
			outExtr << " " << setw(diagWidth) << fieldp->maxima[6] << " " << setw(diagWidth) << fieldp->maxima[7];

			for (spec_iterator = vsp.begin(); spec_iterator != vsp.end(); spec_iterator++){
				outExtr << " " << setw(diagNarrowWidth * 2) << (*spec_iterator)->maxima[6];
			}
		}
		else{
			for (int c = 0; c < 6; c++){
				outExtr << " " << setw(diagWidth) << fieldp->minima[c] << " " << setw(diagWidth) << fieldp->maxima[c];
			}
			outExtr << " " << setw(diagWidth) << fieldp->maxima[6] << " " << setw(diagWidth) << fieldp->maxima[7];

			for (spec_iterator = vsp.begin(); spec_iterator != vsp.end(); spec_iterator++){
				for (int c = 0; c < 7; c++){
					outExtr << " " << setw(diagWideWidth) << (*spec_iterator)->minima[c] << " " << setw(diagWidth) << (*spec_iterator)->maxima[c];
				}
			}
		}
		outExtr << std::endl;
	}
}
void DIAG_MANAGER::callStatSingle(int istep, EM_FIELD* fieldp, std::vector<SPECIE*>* specp){
	std::vector<SPECIE*> vsp = *specp;
	std::vector<SPECIE*>::const_iterator spec_iterator;
	double * ekinSpecies;
	size_t specie = 0;
	double tw = mygrid->dt * istep;

	ekinSpecies = new double[vsp.size()];
	outStatSingle << istep << "\t" << tw;

	double EE[3], BE[3];
	double etotFields = fieldp->getEBenergy(EE, BE);
	double etotKin = 0;
	double tKin;

	specie = 0;
	for (spec_iterator = vsp.begin(); spec_iterator != vsp.end(); spec_iterator++){
		tKin = (*spec_iterator)->getKineticEnergy();
		ekinSpecies[specie++] = tKin;
		etotKin += tKin;
	}
	outStatSingle << "\t" << etotKin + etotFields;
	outStatSingle << "\t" << EE[0] << "\t" << EE[1] << "\t" << EE[2];
	outStatSingle << "\t" << BE[0] << "\t" << BE[1] << "\t" << BE[2];

	for (specie = 0; specie < vsp.size(); specie++){
		outStatSingle << "\t" << ekinSpecies[specie];
	}

	outStatSingle << std::endl;
}
void DIAG_MANAGER::callEMDiag(int istep, EM_FIELD* fieldp, std::string outputName)
{
	std::ofstream myOutputFile;
	if (mygrid->myid == mygrid->master_proc){
		myOutputFile.open(outputName.c_str(), std::ios_base::app);
	}
	fieldp->output_diag(istep, myOutputFile);
	myOutputFile.close();
}
void DIAG_MANAGER::callEMExtrems(int istep, EM_FIELD* fieldp, std::string outputName)
{
	std::ofstream myOutputFile;
	if (mygrid->myid == mygrid->master_proc){
		myOutputFile.open(outputName.c_str(), std::ios_base::app);

	}
	fieldp->output_extrems(istep, myOutputFile);
	myOutputFile.close();
}
void DIAG_MANAGER::callSpecDiag(int istep, SPECIE* spec, std::string outputName)
{
	std::ofstream myOutputFile;
	if (mygrid->myid == mygrid->master_proc){
		myOutputFile.open(outputName.c_str(), std::ios_base::app);
	}
	spec->output_diag(istep, myOutputFile);
	myOutputFile.close();
}
void DIAG_MANAGER::callSpecExtrems(int istep, SPECIE* spec, std::string outputName)
{
	std::ofstream myOutputFile;
	if (mygrid->myid == mygrid->master_proc){
		myOutputFile.open(outputName.c_str(), std::ios_base::app);

	}
	spec->output_extrems(istep, myOutputFile);
	myOutputFile.close();
}
void DIAG_MANAGER::callSpecStat(int istep, SPECIE* spec, std::string fileDiagName, std::string fileExtremName)
{
	std::ofstream fileDiag;
	std::ofstream fileExtrem;

	if (mygrid->myid == mygrid->master_proc){
		fileDiag.open(fileDiagName.c_str(), std::ios::app);
		fileExtrem.open(fileExtremName.c_str(), std::ios::app);
	}


	std::stringstream myname;
	myname << outDir << "/SP_" << spec->name << "_";
	myname << std::setfill('0') << std::setw(SIZE_TIME_DIAG) << std::fixed << std::setprecision(3) << mygrid->time << "_spectrum.txt";

	if (mygrid->myid == mygrid->master_proc)
		tempOut.open(myname.str().c_str());

	spec->output_stat(istep, fileDiag, fileExtrem, tempOut);

	if (mygrid->myid == mygrid->master_proc)
		tempOut.close();

}

void DIAG_MANAGER::callSpecDensityBinary(int istep, SPECIE* spec, CURRENT* curaux){

	curaux->setAllValuesToZero();
	spec->density_deposition_standard(curaux);
	curaux->pbc();

	int uniqueN[3];
	uniqueN[0] = mygrid->uniquePoints[0];
	uniqueN[1] = mygrid->uniquePoints[1];
	uniqueN[2] = mygrid->uniquePoints[2];

	//****************SERIAL DATA FILE WRITE
	if (mygrid->myid == mygrid->master_proc){
		std::stringstream namedat;
		namedat << outDir << "/SP_" << spec->name << "_den_";
		namedat << std::setfill('0') << std::setw(SIZE_TIME_DIAG)
			<< std::fixed << std::setprecision(3)
			<< mygrid->time << ".map";

		outDiagIndex << "SP_DNMAP" << "\t" << spec->name << "\t" << istep << "\t" << std::fixed << std::setprecision(3)
			<< mygrid->time << "\t" << namedat.str() << std::endl; outDiagIndex.flush();

		tempOut.open(namedat.str().c_str());

		for (int c = 0; c < 3; c++)
			tempOut << uniqueN[c] << "\t";
		tempOut << std::endl;
		for (int c = 0; c < 3; c++)
			tempOut << mygrid->rnproc[c] << "\t";
		tempOut << std::endl;

		tempOut << "Ncomp: " << 1 << std::endl;
		integer_or_halfinteger crd = curaux->getDensityCoords();
		tempOut << 0 << ":\t" << (int)crd.x << "\t"
			<< (int)crd.y << "\t"
			<< (int)crd.z << std::endl;
		for (int c = 0; c < 3; c++){
			for (int i = 0; i < uniqueN[c]; i++)
				tempOut << mygrid->cir[c][i] << "  ";
			tempOut << std::endl;
		}

		for (int c = 0; c < 3; c++){
			for (int i = 0; i < uniqueN[c]; i++)
				tempOut << mygrid->chr[c][i] << "  ";
			tempOut << std::endl;
		}

		tempOut.close();
	}
	//****************END OF SERIAL DATA FILE WRITE

	std::stringstream namebin;
	namebin << outDir << "/SP_" << spec->name << "_den_";
	namebin << std::setfill('0') << std::setw(SIZE_TIME_DIAG)
		<< std::fixed << std::setprecision(3)
		<< mygrid->time << ".bin";

	if (mygrid->myid == mygrid->master_proc)
		outDiagIndex << "SP_DNBIN" << "\t" << spec->name << "\t" << istep << "\t" << std::fixed << std::setprecision(3)
		<< mygrid->time << "\t" << namebin.str() << std::endl; outDiagIndex.flush();


	long disp = 0;

	int small_header = 6 * sizeof(int);
	int big_header = 3 * sizeof(int)+3 * sizeof(int)
		+2 * (uniqueN[0] + uniqueN[1] + uniqueN[2])*sizeof(double)
		+sizeof(int)+3 * sizeof(char);

	MPI_File thefile;
	MPI_Status status;

	char *nomefile = new char[namebin.str().size() + 1];
	nomefile[namebin.str().size()] = 0;
	sprintf(nomefile, "%s", namebin.str().c_str());

	MPI_File_open(MPI_COMM_WORLD, nomefile, MPI_MODE_CREATE | MPI_MODE_WRONLY,
		MPI_INFO_NULL, &thefile);

	//****************FILE HEADER

	if (mygrid->myid == mygrid->master_proc){
		MPI_File_write(thefile, uniqueN, 3, MPI_INT, &status);
		MPI_File_write(thefile, mygrid->rnproc, 3, MPI_INT, &status);

		int tNcomp = 1;
		MPI_File_write(thefile, &tNcomp, 1, MPI_INT, &status);
		integer_or_halfinteger crd = curaux->getDensityCoords();
		char tp[3] = { crd.x, crd.y, crd.z };
		MPI_File_write(thefile, tp, 3, MPI_CHAR, &status);

		for (int c = 0; c < 3; c++)
			MPI_File_write(thefile, mygrid->cir[c], uniqueN[c],
			MPI_DOUBLE, &status);

		for (int c = 0; c < 3; c++)
			MPI_File_write(thefile, mygrid->chr[c], uniqueN[c],
			MPI_DOUBLE, &status);

	}




	//****************END OF FILE HEADER

	disp = big_header;

	for (int rank = 0; rank < mygrid->myid; rank++)
		disp += small_header +
		mygrid->proc_totUniquePoints[rank] * sizeof(float);

	MPI_File_set_view(thefile, disp, MPI_FLOAT, MPI_FLOAT,
		(char*) "native", MPI_INFO_NULL);

	//****************CPU HEADER
	{
		int todo[6];
		todo[0] = mygrid->rproc_imin[0][mygrid->rmyid[0]];
		todo[1] = mygrid->rproc_imin[1][mygrid->rmyid[1]];
		todo[2] = mygrid->rproc_imin[2][mygrid->rmyid[2]];
		todo[3] = mygrid->uniquePointsloc[0];
		todo[4] = mygrid->uniquePointsloc[1];
		todo[5] = mygrid->uniquePointsloc[2];

		MPI_File_write(thefile, todo, 6, MPI_INT, &status);
	}
	//****************END OF CPU HEADER

	//****************CPU DENSITY VALUES
	{
		float *todo;
		int Nx, Ny, Nz;

		Nx = mygrid->uniquePointsloc[0];
		Ny = mygrid->uniquePointsloc[1];
		Nz = mygrid->uniquePointsloc[2];

		int size = Nx*Ny*Nz;

		todo = new float[size];

		for (int k = 0; k < Nz; k++)
		for (int j = 0; j < Ny; j++)
		for (int i = 0; i < Nx; i++)
			todo[i + j*Nx + k*Nx*Ny] =
			(float)curaux->density(i, j, k);

		MPI_File_write(thefile, todo, size, MPI_FLOAT, &status);

		delete[]todo;
	}
	//****************END OF CPU DENSITY VALUES

	MPI_File_close(&thefile);
}

void DIAG_MANAGER::callCurrentBinary(int tstep, CURRENT* cur){

	int Ncomp = cur->Ncomp;

	int uniqueN[3];
	uniqueN[0] = mygrid->uniquePoints[0];
	uniqueN[1] = mygrid->uniquePoints[1];
	uniqueN[2] = mygrid->uniquePoints[2];

	//****************SERIAL DATA FILE WRITE
	if (mygrid->myid == mygrid->master_proc){
		std::stringstream namedat;
		namedat << outDir << "/" << "J_";
		namedat << std::setfill('0') << std::setw(SIZE_TIME_DIAG)
			<< std::fixed << std::setprecision(3)
			<< mygrid->time << ".map";

		outDiagIndex << "J_MAP" << "\t" << tstep << "\t" << std::fixed << std::setprecision(3)
			<< mygrid->time << "\t" << namedat.str() << std::endl; outDiagIndex.flush();

		tempOut.open(namedat.str().c_str());

		for (int c = 0; c < 3; c++)
			tempOut << uniqueN[c] << "\t";
		tempOut << std::endl;
		for (int c = 0; c < 3; c++)
			tempOut << mygrid->rnproc[c] << "\t";
		tempOut << std::endl;

		for (int c = 0; c < 3; c++){
			for (int i = 0; i < uniqueN[c]; i++)
				tempOut << mygrid->cir[c][i] << "  ";
			tempOut << std::endl;
		}

		for (int c = 0; c < 3; c++){
			for (int i = 0; i < uniqueN[c]; i++)
				tempOut << mygrid->chr[c][i] << "  ";
			tempOut << std::endl;
		}

		tempOut.close();
	}
	//****************END OF SERIAL DATA FILE WRITE

	std::stringstream namebin;
	namebin << outDir << "/" << "J_";
	namebin << std::setfill('0') << std::setw(SIZE_TIME_DIAG)
		<< std::fixed << std::setprecision(3)
		<< mygrid->time << ".bin";

	if (mygrid->myid == mygrid->master_proc)
		outDiagIndex << "J_BIN" << "\t" << tstep << "\t" << std::fixed << std::setprecision(3)
		<< mygrid->time << "\t" << namebin.str() << std::endl; outDiagIndex.flush();

	long disp = 0;

	int small_header = 6 * sizeof(int);
	int big_header = 3 * sizeof(int)+3 * sizeof(int)+
		2 * (uniqueN[0] + uniqueN[1] + uniqueN[2])*sizeof(double)+
		sizeof(int)+3 * 3 * sizeof(char);

	MPI_File thefile;
	MPI_Status status;

	char *nomefile = new char[namebin.str().size() + 1];
	nomefile[namebin.str().size()] = 0;
	sprintf(nomefile, "%s", namebin.str().c_str());

	MPI_File_open(MPI_COMM_WORLD, nomefile, MPI_MODE_CREATE | MPI_MODE_WRONLY,
		MPI_INFO_NULL, &thefile);

	//****************FILE HEADER

	if (mygrid->myid == mygrid->master_proc){
		MPI_File_write(thefile, uniqueN, 3, MPI_INT, &status);
		MPI_File_write(thefile, mygrid->rnproc, 3, MPI_INT, &status);

		int tNcomp = 3;
		MPI_File_write(thefile, &tNcomp, 1, MPI_INT, &status);
		for (int c = 0; c < 3; c++){
			integer_or_halfinteger crd = cur->getJCoords(c);
			char tp[3] = { crd.x, crd.y, crd.z };
			MPI_File_write(thefile, tp, 3, MPI_CHAR, &status);
		}


		for (int c = 0; c < 3; c++)
			MPI_File_write(thefile, mygrid->cir[c], uniqueN[c],
			MPI_DOUBLE, &status);
		for (int c = 0; c < 3; c++)
			MPI_File_write(thefile, mygrid->chr[c], uniqueN[c],
			MPI_DOUBLE, &status);
	}

	//****************END OF FILE HEADER

	disp = big_header;

	for (int rank = 0; rank < mygrid->myid; rank++)
		disp += small_header +
		mygrid->proc_totUniquePoints[rank] * sizeof(float)*Ncomp;

	MPI_File_set_view(thefile, disp, MPI_FLOAT, MPI_FLOAT,
		(char*) "native", MPI_INFO_NULL);

	//****************CPU HEADER
	{
		int todo[6];
		todo[0] = mygrid->rproc_imin[0][mygrid->rmyid[0]];
		todo[1] = mygrid->rproc_imin[1][mygrid->rmyid[1]];
		todo[2] = mygrid->rproc_imin[2][mygrid->rmyid[2]];
		todo[3] = mygrid->uniquePointsloc[0];
		todo[4] = mygrid->uniquePointsloc[1];
		todo[5] = mygrid->uniquePointsloc[2];

		MPI_File_write(thefile, todo, 6, MPI_INT, &status);
	}
	//****************END OF CPU HEADER


	//****************CPU CURRENT VALUES
	{
		float *todo;
		int Nx, Ny, Nz;

		int Nc = Ncomp;

		Nx = mygrid->uniquePointsloc[0];
		Ny = mygrid->uniquePointsloc[1];
		Nz = mygrid->uniquePointsloc[2];

		int size = Nx*Ny*Nz*Nc;

		todo = new float[size];

		for (int k = 0; k < Nz; k++)
		for (int j = 0; j < Ny; j++)
		for (int i = 0; i < Nx; i++)
		for (int c = 0; c < Nc; c++)
			todo[c + i*Nc + j*Nx*Nc + k*Nx*Ny*Nc] =
			(float)cur->JJ(c, i, j, k);

		MPI_File_write(thefile, todo, size, MPI_FLOAT, &status);

		delete[]todo;
	}
	//****************END OF CPU CURRENT VALUES

	MPI_File_close(&thefile);

}

void DIAG_MANAGER::interpEB(double pos[3], EM_FIELD* ebfield, double E[3], double B[3]){
	int hii[3], wii[3];
	double hiw[3][3], wiw[3][3];
	double rr, rh, rr2, rh2;
	int i1, j1, k1, i2, j2, k2;
	double dvol;

	for (int c = 0; c < 3; c++){
		hiw[c][1] = wiw[c][1] = 1;
		hii[c] = wii[c] = 0;
	}
	for (int c = 0; c < acc.dimensions; c++)
	{
		rr = mygrid->dri[c] * (pos[c] - mygrid->rminloc[c]);
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
	E[0] = E[1] = E[2] = B[0] = B[1] = B[2] = 0;

	switch (acc.dimensions)
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
						E[0] += ebfield->E0(i2, j1, k1)*dvol;  //Ex
					dvol = wiw[0][i] * hiw[1][j] * wiw[2][k],
						E[1] += ebfield->E1(i1, j2, k1)*dvol;  //Ey
					dvol = wiw[0][i] * wiw[1][j] * hiw[2][k],
						E[2] += ebfield->E2(i1, j1, k2)*dvol;  //Ez

					dvol = wiw[0][i] * hiw[1][j] * hiw[2][k],
						B[0] += ebfield->B0(i1, j2, k2)*dvol;  //Bx
					dvol = hiw[0][i] * wiw[1][j] * hiw[2][k],
						B[1] += ebfield->B1(i2, j1, k2)*dvol;  //By
					dvol = hiw[0][i] * hiw[1][j] * wiw[2][k],
						B[2] += ebfield->B2(i2, j2, k1)*dvol;  //Bz
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
					E[0] += ebfield->E0(i2, j1, k1)*dvol;  //Ex
				dvol = wiw[0][i] * hiw[1][j],
					E[1] += ebfield->E1(i1, j2, k1)*dvol;  //Ey
				dvol = wiw[0][i] * wiw[1][j],
					E[2] += ebfield->E2(i1, j1, k2)*dvol;  //Ez

				dvol = wiw[0][i] * hiw[1][j],
					B[0] += ebfield->B0(i1, j2, k2)*dvol;  //Bx
				dvol = hiw[0][i] * wiw[1][j],
					B[1] += ebfield->B1(i2, j1, k2)*dvol;  //By
				dvol = hiw[0][i] * hiw[1][j],
					B[2] += ebfield->B2(i2, j2, k1)*dvol;  //Bz
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
				E[0] += ebfield->E0(i2, j1, k1)*dvol;  //Ex
			dvol = wiw[0][i],
				E[1] += ebfield->E1(i1, j2, k1)*dvol;  //Ey
			dvol = wiw[0][i],
				E[2] += ebfield->E2(i1, j1, k2)*dvol;  //Ez

			dvol = wiw[0][i],
				B[0] += ebfield->B0(i1, j2, k2)*dvol;  //Bx
			dvol = hiw[0][i],
				B[1] += ebfield->B1(i2, j1, k2)*dvol;  //By
			dvol = hiw[0][i],
				B[2] += ebfield->B2(i2, j2, k1)*dvol;  //Bz
		}
		break;
	}

}

void DIAG_MANAGER::interpJ(double pos[3], CURRENT* curp, double J[3]){
	int wii[3];
	double rr;

	wii[0] = wii[1] = wii[2] = 0;

	for (int c = 0; c < acc.dimensions; c++)
	{
		rr = mygrid->dri[c] * (pos[c] - mygrid->rminloc[c]);
		wii[c] = (int)floor(rr + 0.5); //whole integer int
	}

	J[0] = curp->Jx(wii[0], wii[1], wii[2]);
	J[1] = curp->Jy(wii[0], wii[1], wii[2]);
	J[2] = curp->Jz(wii[0], wii[1], wii[2]);
}

void DIAG_MANAGER::callEMJProbe(int tstep, double position[3], EM_FIELD* fieldp,
	CURRENT* curp, std::string fileName){

	if (position[0] >= mygrid->rminloc[0] && position[0] < mygrid->rmaxloc[0] &&
		position[1] >= mygrid->rminloc[1] && position[1] < mygrid->rmaxloc[1] &&
		position[2] >= mygrid->rminloc[2] && position[2] < mygrid->rmaxloc[2]){


		double E[3];
		double B[3];
		double J[3];
		interpEB(position, fieldp, E, B);
		interpJ(position, curp, J);

		std::ofstream outf;
		outf.open(fileName.c_str(), std::ios_base::app);

		outf << " " << setw(diagNarrowWidth) << tstep
			<< std::fixed << std::setprecision(7)
			<< " " << setw(diagWidth) << mygrid->time
			<< " " << setw(diagWidth) << E[0]
			<< " " << setw(diagWidth) << E[1]
			<< " " << setw(diagWidth) << E[2]
			<< " " << setw(diagWidth) << B[0]
			<< " " << setw(diagWidth) << B[1]
			<< " " << setw(diagWidth) << B[2]
			<< " " << setw(diagWidth) << J[0]
			<< " " << setw(diagWidth) << J[1]
			<< " " << setw(diagWidth) << J[2]
			<< std::endl;

		outf.close();

	}



}

