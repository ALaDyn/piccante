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

bool requestCompTime(const request &first, const request &second){
	if (first.itime != second.itime)
		return (first.itime < second.itime);
	else if (first.type != second.type)
		return (first.type < second.type);
	else if (first.target != second.target)
		return (first.target < second.target);
	else
		return false;

}

bool requestCompUnique(const request &first, const request &second){
	return ((first.itime == second.itime) &&
		(first.type == second.type) &&
		(first.target == second.target));
}


OUTPUT_MANAGER::OUTPUT_MANAGER(GRID* _mygrid, EM_FIELD* _myfield, CURRENT* _mycurrent, std::vector<SPECIE*> _myspecies)
{
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
}

OUTPUT_MANAGER::~OUTPUT_MANAGER(){

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

	std::stringstream ss;
	ss << outputDir << "/diag.dat";

	diagFileName = ss.str();

	std::ofstream of1;

	of1.open(diagFileName.c_str(), std::ofstream::out | std::ofstream::trunc);

	of1 << " " << setw(diagNarrowWidth) << "#istep" << " " << setw(diagWidth) << "time";
	of1 << " " << setw(diagWidth) << "Etot";
	of1 << " " << setw(diagWidth) << "Ex2" << " " << setw(diagWidth) << "Ey2" << " " << setw(diagWidth) << "Ez2";
	of1 << " " << setw(diagWidth) << "Bx2" << " " << setw(diagWidth) << "By2" << " " << setw(diagWidth) << "Bz2";

	std::vector<SPECIE*>::const_iterator spec_iterator;

	for (spec_iterator = myspecies.begin(); spec_iterator != myspecies.end(); spec_iterator++){
		of1 << " " << setw(diagWidth) << (*spec_iterator)->name;
	}

	of1 << std::endl;

	of1.close();
}

void OUTPUT_MANAGER::createExtremaFiles(){

	if (checkGrid()){
		if (mygrid->myid != mygrid->master_proc)
			return;
	}
	else{
		if (mygrid->myid != 0)
			return;
	}

	if (checkEMField()){
		std::stringstream ss0;
		ss0 << outputDir << "/EXTREMES_EMfield.dat";
		extremaFieldFileName = ss0.str();

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

void OUTPUT_MANAGER::initialize(std::string _outputDir){

    if ( !boost::filesystem::exists(_outputDir) ){
        boost::filesystem::create_directories(_outputDir);
    }
	outputDir = _outputDir;
	prepareOutputMap();

	if (isThereDiag){
		createDiagFile();
		createExtremaFiles();
	}

	amIInit = true;
}

void OUTPUT_MANAGER::close(){

}

bool OUTPUT_MANAGER::checkGrid(){
	if (!isThereGrid){
		int myid;
		MPI_Comm_rank(MPI_COMM_WORLD, &myid);
		if (myid == 0){
			std::cout << "WARNING! No valid GRID pointer provided. Output request will be ignored." << std::endl;
		}

	}
	return isThereGrid;
}

bool OUTPUT_MANAGER::checkEMField(){
	if (!isThereField){
		int myid;
		MPI_Comm_rank(MPI_COMM_WORLD, &myid);
		if (myid == 0){
			std::cout << "WARNING! No valid FIELD pointer provided. Output request will be ignored." << std::endl;
		}
	}
	return isThereField;
}

bool OUTPUT_MANAGER::checkCurrent(){
	if (!isThereCurrent){
		int myid;
		MPI_Comm_rank(MPI_COMM_WORLD, &myid);
		if (myid == 0){
			std::cout << "WARNING! No valid CURRENT pointer provided. Output request will be ignored." << std::endl;
		}
	}
	return isThereCurrent;
}

bool OUTPUT_MANAGER::checkSpecies(){
	if (!isThereSpecList){
		int myid;
		MPI_Comm_rank(MPI_COMM_WORLD, &myid);
		if (myid == 0){
			std::cout << "WARNING! Empty SPECIES list provided. Output request will be ignored." << std::endl;
		}
	}
	return isThereSpecList;
}

int OUTPUT_MANAGER::getIntTime(double dtime){
	if (dtime == 0.0)
		return 0;

	int Nsteps = mygrid->getTotalNumberOfTimesteps();
	int step = (int)floor(dtime / mygrid->dt + 0.5);

	return (step <= Nsteps) ? (step) : (-1);
}

int OUTPUT_MANAGER::findSpecName(std::string name){
	int pos = 0;
	for (std::vector<SPECIE*>::iterator it = myspecies.begin(); it != myspecies.end(); it++){
		if ((*it)->getName() == name){
			return pos;
		}
		pos++;
	}
	int myid;
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	if (myid == 0){
		std::cout << "WARNING! Species" << name << " not in species list !" << std::endl;
	}
	return -1;
}

void OUTPUT_MANAGER::addRequestToList(std::list<request>& reqList, diagType type, int target, double startTime, double frequency, double endTime){
	std::list<request> tempList;

	for (double ttime = startTime; ttime <= endTime; ttime += frequency){
		request req;
		req.dtime = ttime;
		req.itime = getIntTime(ttime);
		req.type = type;
		req.target = target;
		tempList.push_back(req);
	}
	tempList.sort(requestCompTime);
	reqList.merge(tempList, requestCompTime);
}

void OUTPUT_MANAGER::addEMFieldBinaryFrom(double startTime, double frequency){
	if (!(checkGrid() && checkEMField()))
		return;
	double endSimTime = mygrid->dt * mygrid->getTotalNumberOfTimesteps();
	addRequestToList(requestList, OUT_EM_FIELD_BINARY, 0, startTime, frequency, endSimTime);

}

void OUTPUT_MANAGER::addEMFieldBinaryAt(double atTime){
	if (!(checkGrid() && checkEMField()))
		return;
	addRequestToList(requestList, OUT_EM_FIELD_BINARY, 0, atTime, 1.0, atTime);
}

void OUTPUT_MANAGER::addEMFieldBinaryFromTo(double startTime, double frequency, double endTime){
	if (!(checkGrid() && checkEMField()))
		return;
	addRequestToList(requestList, OUT_EM_FIELD_BINARY, 0, startTime, frequency, endTime);
}

void OUTPUT_MANAGER::addSpecDensityBinaryFrom(std::string name, double startTime, double frequency){
	if (!(checkGrid() && checkSpecies() && checkCurrent()))
		return;
	double endSimTime = mygrid->dt * mygrid->getTotalNumberOfTimesteps();
	int specNum = findSpecName(name);
	if (specNum < 0)
		return;
	addRequestToList(requestList, OUT_SPEC_DENSITY_BINARY, specNum, startTime, frequency, endSimTime);

}

void OUTPUT_MANAGER::addSpecDensityBinaryAt(std::string name, double atTime){
	if (!(checkGrid() && checkSpecies() && checkCurrent()))
		return;
	int specNum = findSpecName(name);
	if (specNum < 0)
		return;
	addRequestToList(requestList, OUT_SPEC_DENSITY_BINARY, specNum, atTime, 1.0, atTime);
}

void OUTPUT_MANAGER::addSpecDensityBinaryFromTo(std::string name, double startTime, double frequency, double endTime){
	if (!(checkGrid() && checkSpecies() && checkCurrent()))
		return;
	int specNum = findSpecName(name);
	if (specNum < 0)
		return;
	addRequestToList(requestList, OUT_SPEC_DENSITY_BINARY, specNum, startTime, frequency, endTime);
}

void OUTPUT_MANAGER::addCurrentBinaryFrom(double startTime, double frequency){
	if (!(checkGrid() && checkCurrent()))
		return;
	double endSimTime = mygrid->dt * mygrid->getTotalNumberOfTimesteps();
	addRequestToList(requestList, OUT_CURRENT_BINARY, 0, startTime, frequency, endSimTime);

}

void OUTPUT_MANAGER::addCurrentBinaryAt(double atTime){
	if (!(checkGrid() && checkCurrent()))
		return;
	addRequestToList(requestList, OUT_CURRENT_BINARY, 0, atTime, 1.0, atTime);
}

void OUTPUT_MANAGER::addCurrentBinaryFromTo(double startTime, double frequency, double endTime){
	if (!(checkGrid() && checkCurrent()))
		return;
	addRequestToList(requestList, OUT_CURRENT_BINARY, 0, startTime, frequency, endTime);
}

void OUTPUT_MANAGER::addSpecPhaseSpaceBinaryFrom(std::string name, double startTime, double frequency){
	if (!(checkGrid() && checkSpecies()))
		return;
	double endSimTime = mygrid->dt * mygrid->getTotalNumberOfTimesteps();
	int specNum = findSpecName(name);
	if (specNum < 0)
		return;
	addRequestToList(requestList, OUT_SPEC_PHASE_SPACE_BINARY, specNum, startTime, frequency, endSimTime);

}

void OUTPUT_MANAGER::addSpecPhaseSpaceBinaryAt(std::string name, double atTime){
	if (!(checkGrid() && checkSpecies()))
		return;
	int specNum = findSpecName(name);
	if (specNum < 0)
		return;
	addRequestToList(requestList, OUT_SPEC_PHASE_SPACE_BINARY, specNum, atTime, 1.0, atTime);
}

void OUTPUT_MANAGER::addSpecPhaseSpaceBinaryFromTo(std::string name, double startTime, double frequency, double endTime){
	if (!(checkGrid() && checkSpecies()))
		return;
	int specNum = findSpecName(name);
	if (specNum < 0)
		return;
	addRequestToList(requestList, OUT_SPEC_PHASE_SPACE_BINARY, specNum, startTime, frequency, endTime);
}

void OUTPUT_MANAGER::addDiagFrom(double startTime, double frequency){
	if (!(checkGrid() && checkCurrent() && checkSpecies() && checkEMField()))
		return;
	double endSimTime = mygrid->dt * mygrid->getTotalNumberOfTimesteps();
	addRequestToList(requestList, OUT_DIAG, 0, startTime, frequency, endSimTime);
	isThereDiag = true;
}

void OUTPUT_MANAGER::addDiagAt(double atTime){
	if (!(checkGrid() && checkCurrent() && checkSpecies() && checkEMField()))
		return;
	addRequestToList(requestList, OUT_DIAG, 0, atTime, 1.0, atTime);
	isThereDiag = true;
}

void OUTPUT_MANAGER::addDiagFromTo(double startTime, double frequency, double endTime){
	if (!(checkGrid() && checkCurrent() && checkSpecies() && checkEMField()))
		return;
	addRequestToList(requestList, OUT_DIAG, 0, startTime, frequency, endTime);
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

void OUTPUT_MANAGER::autoVisualDiag(){
	int myid;
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	if (myid == 0){
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

void OUTPUT_MANAGER::processOutputEntry(request req){
	switch (req.type){
	case OUT_EM_FIELD_BINARY:
		callEMFieldBinary(req);
		break;

	case OUT_SPEC_DENSITY_BINARY:
		callSpecDensityBinary(req);
		break;

	case OUT_CURRENT_BINARY:
		callCurrentBinary(req);
		break;

	case OUT_SPEC_PHASE_SPACE_BINARY:
		callSpecPhaseSpaceBinary(req);
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

void OUTPUT_MANAGER::writeEMFieldMap(std::ofstream &output, request req){
	int uniqueN[3];
	uniqueN[0] = mygrid->uniquePoints[0];
	uniqueN[1] = mygrid->uniquePoints[1];
	uniqueN[2] = mygrid->uniquePoints[2];

	int Ncomp = myfield->getNcomp();

	for (int c = 0; c < 3; c++)
		output << uniqueN[c] << "\t";

	output << std::endl;

	for (int c = 0; c < 3; c++)
		output << mygrid->rnproc[c] << "\t";

	output << std::endl;

	output << "Ncomp: " << Ncomp << std::endl;

	for (int c = 0; c < Ncomp; c++){
		integer_or_halfinteger crd = myfield->getCompCoords(c);
		output << c << ":\t" << (int)crd.x << "\t"
			<< (int)crd.y << "\t"
			<< (int)crd.z << std::endl;
	}

	for (int c = 0; c < 3; c++){
		for (int i = 0; i < uniqueN[c]; i++)
			output << mygrid->cir[c][i] << "  ";
		output << std::endl;
	}
	for (int c = 0; c < 3; c++){
		for (int i = 0; i < uniqueN[c]; i++)
			output << mygrid->chr[c][i] << "  ";
		output << std::endl;
	}
}

void OUTPUT_MANAGER::writeEMFieldBinary(std::string fileName, request req){
	int Ncomp = myfield->getNcomp();
	int uniqueN[3];
	uniqueN[0] = mygrid->uniquePoints[0];
	uniqueN[1] = mygrid->uniquePoints[1];
	uniqueN[2] = mygrid->uniquePoints[2];
    MPI_Offset disp = 0;
	int small_header = 6 * sizeof(int);
	int big_header = 3 * sizeof(int)+3 * sizeof(int)
		+2 * (uniqueN[0] + uniqueN[1] + uniqueN[2])*sizeof(double)
		+sizeof(int)+myfield->getNcomp() * 3 * sizeof(char);

	MPI_File thefile;
	MPI_Status status;

	char* nomefile = new char[fileName.size() + 1];

	nomefile[fileName.size()] = 0;
	sprintf(nomefile, "%s", fileName.c_str());

    int openTest;
    openTest=MPI_File_open(MPI_COMM_WORLD, nomefile, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &thefile);
    if(openTest=0){
        std::cout<<"a problem occurred when trying to mpi_file_open in writeEMFieldBinary"<<std::endl;
        std::cout<<"myrank="<<mygrid->myid<<" status="<<openTest<<std::endl;
        exit(33);
    }
	//+++++++++++ FILE HEADER  +++++++++++++++++++++
    if (mygrid->myid == mygrid->master_proc){
        MPI_File_set_view(thefile, 0, MPI_FLOAT, MPI_FLOAT, (char *) "native", MPI_INFO_NULL);
		MPI_File_write(thefile, uniqueN, 3, MPI_INT, &status);
		MPI_File_write(thefile, mygrid->rnproc, 3, MPI_INT, &status);
		MPI_File_write(thefile, &Ncomp, 1, MPI_INT, &status);
		for (int c = 0; c < Ncomp; c++){
			integer_or_halfinteger crd = myfield->getCompCoords(c);
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
    if(disp<0)
    {
        std::cout<<"a problem occurred when trying to mpi_file_set_view in writeEMFieldBinary"<<std::endl;
        std::cout<<"myrank="<<mygrid->myid<<" disp="<<disp<<std::endl;
        exit(33);
    }
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
		int Nx, Ny, Nz;
		Nx = mygrid->uniquePointsloc[0];
		Ny = mygrid->uniquePointsloc[1];
		Nz = mygrid->uniquePointsloc[2];
		int size = Ncomp*Nx*Ny*Nz;
		todo = new float[size];
		for (int k = 0; k < Nz; k++)
			for (int j = 0; j < Ny; j++)
				for (int i = 0; i < Nx; i++)
					for (int c = 0; c < Ncomp; c++)
						todo[c + i*Ncomp + j*Nx*Ncomp + k*Ny*Nx*Ncomp] = (float)myfield->VEB(c, i, j, k);
		MPI_File_write(thefile, todo, size, MPI_FLOAT, &status);
		delete[]todo;
	}

	MPI_File_close(&thefile);
	//////////////////////////// END of collective binary file write
}

void OUTPUT_MANAGER::callEMFieldBinary(request req){

	if (mygrid->myid == mygrid->master_proc){
		std::string nameMap = composeOutputName(outputDir, "EMfield", "", req.dtime, ".map");
		std::ofstream of1;
		of1.open(nameMap.c_str());
		writeEMFieldMap(of1, req);
		of1.close();
	}

	std::string nameBin = composeOutputName(outputDir, "EMfield", "", req.dtime, ".bin");

	writeEMFieldBinary(nameBin, req);
}

void OUTPUT_MANAGER::writeSpecDensityMap(std::ofstream &output, request req){
	int uniqueN[3];
	uniqueN[0] = mygrid->uniquePoints[0];
	uniqueN[1] = mygrid->uniquePoints[1];
	uniqueN[2] = mygrid->uniquePoints[2];

	for (int c = 0; c < 3; c++)
		output << uniqueN[c] << "\t";
	output << std::endl;
	for (int c = 0; c < 3; c++)
		output << mygrid->rnproc[c] << "\t";
	output << std::endl;

	output << "Ncomp: " << 1 << std::endl;
	integer_or_halfinteger crd = mycurrent->getDensityCoords();
	output << 0 << ":\t" << (int)crd.x << "\t"
		<< (int)crd.y << "\t"
		<< (int)crd.z << std::endl;
	for (int c = 0; c < 3; c++){
		for (int i = 0; i < uniqueN[c]; i++)
			output << mygrid->cir[c][i] << "  ";
		output << std::endl;
	}

	for (int c = 0; c < 3; c++){
		for (int i = 0; i < uniqueN[c]; i++)
			output << mygrid->chr[c][i] << "  ";
		output << std::endl;
	}

}

void OUTPUT_MANAGER::writeSpecDensityBinary(std::string fileName, request req){
	int uniqueN[3];
	uniqueN[0] = mygrid->uniquePoints[0];
	uniqueN[1] = mygrid->uniquePoints[1];
	uniqueN[2] = mygrid->uniquePoints[2];

    MPI_Offset disp = 0;

	int small_header = 6 * sizeof(int);
	int big_header = 3 * sizeof(int)+3 * sizeof(int)
		+2 * (uniqueN[0] + uniqueN[1] + uniqueN[2])*sizeof(double)
		+sizeof(int)+3 * sizeof(char);

	MPI_File thefile;
	MPI_Status status;

	char* nomefile = new char[fileName.size() + 1];

	nomefile[fileName.size()] = 0;
	sprintf(nomefile, "%s", fileName.c_str());

	MPI_File_open(MPI_COMM_WORLD, nomefile, MPI_MODE_CREATE | MPI_MODE_WRONLY,
		MPI_INFO_NULL, &thefile);

	//****************FILE HEADER

	if (mygrid->myid == mygrid->master_proc){
        MPI_File_set_view(thefile, 0, MPI_FLOAT, MPI_FLOAT, (char *) "native", MPI_INFO_NULL);
        MPI_File_write(thefile, uniqueN, 3, MPI_INT, &status);
		MPI_File_write(thefile, mygrid->rnproc, 3, MPI_INT, &status);

		int tNcomp = 1;
		MPI_File_write(thefile, &tNcomp, 1, MPI_INT, &status);
		integer_or_halfinteger crd = mycurrent->getDensityCoords();
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

    if(disp<0)
    {
        std::cout<<"a problem occurred when trying to mpi_file_set_view in writeSpecDensityBinary"<<std::endl;
        std::cout<<"myrank="<<mygrid->myid<<" disp="<<disp<<std::endl;
        exit(33);
    }
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
					(float)mycurrent->density(i, j, k);

		MPI_File_write(thefile, todo, size, MPI_FLOAT, &status);

		delete[]todo;
	}
	//****************END OF CPU DENSITY VALUES

	MPI_File_close(&thefile);
	delete[] nomefile;

}

void OUTPUT_MANAGER::callSpecDensityBinary(request req){
	mycurrent->eraseDensity();
	myspecies[req.target]->density_deposition_standard(mycurrent);
	mycurrent->pbc();

	std::string name = myspecies[req.target]->name;

	if (mygrid->myid == mygrid->master_proc){
		std::string nameMap = composeOutputName(outputDir, "DENS", name, req.dtime, ".map");
		std::ofstream of1;
		of1.open(nameMap.c_str());
		writeSpecDensityMap(of1, req);
		of1.close();
	}

	std::string nameBin = composeOutputName(outputDir, "DENS", name, req.dtime, ".bin");

	writeSpecDensityBinary(nameBin, req);

}

void  OUTPUT_MANAGER::writeCurrentMap(std::ofstream &output, request req){
	int Ncomp = mycurrent->Ncomp;

	int uniqueN[3];
	uniqueN[0] = mygrid->uniquePoints[0];
	uniqueN[1] = mygrid->uniquePoints[1];
	uniqueN[2] = mygrid->uniquePoints[2];

	for (int c = 0; c < 3; c++)
		output << uniqueN[c] << "\t";
	output << std::endl;
	for (int c = 0; c < 3; c++)
		output << mygrid->rnproc[c] << "\t";
	output << std::endl;

	for (int c = 0; c < 3; c++){
		for (int i = 0; i < uniqueN[c]; i++)
			output << mygrid->cir[c][i] << "  ";
		output << std::endl;
	}

	for (int c = 0; c < 3; c++){
		for (int i = 0; i < uniqueN[c]; i++)
			output << mygrid->chr[c][i] << "  ";
		output << std::endl;
	}
}

void  OUTPUT_MANAGER::writeCurrentBinary(std::string fileName, request req){

	int Ncomp = 3;

	int uniqueN[3];
	uniqueN[0] = mygrid->uniquePoints[0];
	uniqueN[1] = mygrid->uniquePoints[1];
	uniqueN[2] = mygrid->uniquePoints[2];


    MPI_Offset disp = 0;

	int small_header = 6 * sizeof(int);
	int big_header = 3 * sizeof(int)+3 * sizeof(int)+
		2 * (uniqueN[0] + uniqueN[1] + uniqueN[2])*sizeof(double)+
		sizeof(int)+3 * 3 * sizeof(char);

	MPI_File thefile;
	MPI_Status status;

	char *nomefile = new char[fileName.size() + 1];
	nomefile[fileName.size()] = 0;
	sprintf(nomefile, "%s", fileName.c_str());

	MPI_File_open(MPI_COMM_WORLD, nomefile, MPI_MODE_CREATE | MPI_MODE_WRONLY,
		MPI_INFO_NULL, &thefile);

	//****************FILE HEADER

	if (mygrid->myid == mygrid->master_proc){
        MPI_File_set_view(thefile, 0, MPI_FLOAT, MPI_FLOAT, (char *) "native", MPI_INFO_NULL);
        MPI_File_write(thefile, uniqueN, 3, MPI_INT, &status);
		MPI_File_write(thefile, mygrid->rnproc, 3, MPI_INT, &status);

		int tNcomp = 3;
		MPI_File_write(thefile, &tNcomp, 1, MPI_INT, &status);
		for (int c = 0; c < 3; c++){
			integer_or_halfinteger crd = mycurrent->getJCoords(c);
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
        disp += small_header + mygrid->proc_totUniquePoints[rank] * sizeof(float)*Ncomp;
    if(disp<0)
    {
        std::cout<<"a problem occurred when trying to mpi_file_set_view in writeCurrentBinary"<<std::endl;
        std::cout<<"myrank="<<mygrid->myid<<" disp="<<disp<<std::endl;
        exit(33);
    }

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
						(float)mycurrent->JJ(c, i, j, k);

		MPI_File_write(thefile, todo, size, MPI_FLOAT, &status);

		delete[]todo;
	}
	//****************END OF CPU CURRENT VALUES

	MPI_File_close(&thefile);

	delete[] nomefile;
}

void  OUTPUT_MANAGER::callCurrentBinary(request req){
	if (mygrid->myid == mygrid->master_proc){
		std::string nameMap = composeOutputName(outputDir, "J", "", req.dtime, ".map");
		std::ofstream of1;
		of1.open(nameMap.c_str());
		writeCurrentMap(of1, req);
		of1.close();
	}
	std::string nameBin = composeOutputName(outputDir, "J", "", req.dtime, ".bin");

	writeCurrentBinary(nameBin, req);

}

void OUTPUT_MANAGER::writeSpecPhaseSpaceBinary(std::string fileName, request req){

	SPECIE* spec = myspecies[req.target];

	int* NfloatLoc = new int[mygrid->nproc];
	NfloatLoc[mygrid->myid] = spec->Np*spec->Ncomp;

	MPI_Allgather(MPI_IN_PLACE, 1, MPI_INT, NfloatLoc, 1, MPI_INT, MPI_COMM_WORLD);

    MPI_Offset disp = 0;
	for (int pp = 0; pp < mygrid->myid; pp++)
        disp += (MPI_Offset)(NfloatLoc[pp] * sizeof(float));
	MPI_File thefile;

	char *nomefile = new char[fileName.size() + 1];
	nomefile[fileName.size()] = 0;
	sprintf(nomefile, "%s", fileName.c_str());

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
				buf[c + p*spec->Ncomp] = (float) spec->ru(c, p + dimensione*i);
			}
		}
		MPI_File_write(thefile, buf, dimensione*spec->Ncomp, MPI_FLOAT, &status);
	}
	for (int p = 0; p < resto; p++){
		for (int c = 0; c < spec->Ncomp; c++){
			buf[c + p*spec->Ncomp] = (float) spec->ru(c, p + dimensione*passaggi);
		}
	}
	MPI_File_write(thefile, buf, resto*spec->Ncomp, MPI_FLOAT, &status);
	MPI_File_close(&thefile);
	delete[]buf;
	delete[] NfloatLoc;

}
void OUTPUT_MANAGER::callSpecPhaseSpaceBinary(request req){
	std::string name = myspecies[req.target]->name;

	std::string nameBin = composeOutputName(outputDir, "PHASESPACE", name, req.dtime, ".bin");

	writeSpecPhaseSpaceBinary(nameBin, req);
}


void OUTPUT_MANAGER::callDiag(request req){
	std::vector<SPECIE*>::const_iterator spec_iterator;
	double * ekinSpecies;
	size_t specie = 0;

	double tw = req.dtime;

	ekinSpecies = new double[myspecies.size()];

	//double t_spec_extrema[SPEC_DIAG_COMP];
	// double field_extrema[FIELD_DIAG_COMP];

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
		ekinSpecies[specie++] = (*spec_iterator)->total_energy;
		etotKin += (*spec_iterator)->total_energy;
	}
	double etot = etotKin + etotFields;

	if (mygrid->myid == mygrid->master_proc){
		std::ofstream outStat;
		outStat.open(diagFileName.c_str(), std::ios::app);

		outStat << " " << setw(diagNarrowWidth) << req.itime << " " << setw(diagWidth) << tw;
		outStat << " " << setw(diagWidth) << etot;
		outStat << " " << setw(diagWidth) << EE[0] << " " << setw(diagWidth) << EE[1] << " " << setw(diagWidth) << EE[2];
		outStat << " " << setw(diagWidth) << BE[0] << " " << setw(diagWidth) << BE[1] << " " << setw(diagWidth) << BE[2];

		for (specie = 0; specie < myspecies.size(); specie++){
			outStat << " " << setw(diagWidth) << ekinSpecies[specie];
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
