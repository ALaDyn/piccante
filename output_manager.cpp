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

bool OUTPUT_MANAGER::isInMyDomain(double rr[3]){
    if (rr[0] >= mygrid->rminloc[0] && rr[0] < mygrid->rmaxloc[0]){
        if (mygrid->accesso.dimensions < 2 || (rr[1] >= mygrid->rminloc[1] && rr[1] < mygrid->rmaxloc[1])){
            if (mygrid->accesso.dimensions < 3 || (rr[2] >= mygrid->rminloc[2] && rr[2] < mygrid->rmaxloc[2])){
                return true;
            }
        }
    }
    return false;
}
bool OUTPUT_MANAGER::amIInTheSubDomain(request req){
    double rmin[3], rmax[3];
    for (int c = 0; c < 3; c++){
        rmin[c] = myDomains[req.domain]->rmin[c];
        rmax[c] = myDomains[req.domain]->rmax[c];
    }
    if (rmax[0] >= mygrid->rminloc[0] && rmin[0] < mygrid->rmaxloc[0]){
        if (mygrid->accesso.dimensions < 2 || (rmax[1] >= mygrid->rminloc[1] && rmin[1] < mygrid->rmaxloc[1])){
            if (mygrid->accesso.dimensions < 3 || (rmax[2] >= mygrid->rminloc[2] && rmin[2] < mygrid->rmaxloc[2])){
                return true;
            }
        }
    }
    return false;
}


void OUTPUT_MANAGER::nearestInt(double rr[], int *ri, int *globalri){
    int c;
    if (mygrid->isStretched()){
        for (c = 0; c < mygrid->accesso.dimensions; c++){
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
        for (c = 0; c < mygrid->accesso.dimensions; c++){
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
int OUTPUT_MANAGER::findIndexMin(double val, double* coords, int numcoords){
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

int OUTPUT_MANAGER::findIndexMax(double val, double* coords, int numcoords){
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
        if (c < mygrid->accesso.dimensions)
            remains[c] = remainingCoord[c];
        else
            remains[c] = 0;
    }
}


void OUTPUT_MANAGER::findIntLocalBoundaries(double rmin[3], double rmax[3], int *imin, int *imax){
    for (int c = 0; c < 3; c++){
        imin[c] = findIndexMin(rmin[c], mygrid->cirloc[c], mygrid->Nloc[c]);
        imax[c] = findIndexMax(rmax[c], mygrid->cirloc[c], mygrid->Nloc[c]);
    }
}
void OUTPUT_MANAGER::findIntGlobalBoundaries(double rmin[3], double rmax[3], int *imin, int *imax){
    for (int c = 0; c < 3; c++){
        imin[c] = findIndexMin(rmin[c], mygrid->cir[c], mygrid->NGridNodes[c]);
        imax[c] = findIndexMax(rmax[c], mygrid->cir[c], mygrid->NGridNodes[c]);
    }
}
void OUTPUT_MANAGER::findNumberOfProc(int *Nproc, int imin[3], int imax[3], int remains[3]){
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

void OUTPUT_MANAGER::findGlobalSubdomainUniquePoints(int *uniqueN, int imin[3], int imax[3], int remains[3]){
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

void OUTPUT_MANAGER::findLocalSubdomainUniquePoints(int *uniqueLocN, int locimin[3], int locimax[3], int remains[3]){
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

    outDomain *domain1 = new outDomain;
    domain1->overrideFlag = true;
    myDomains.push_back(domain1);

}

OUTPUT_MANAGER::~OUTPUT_MANAGER(){
    for (std::vector<outDomain*>::iterator it = myDomains.begin(); it != myDomains.end(); ++it)
        delete(*it);


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
    _newoutputDir=ss.str();
    if (mygrid->myid == mygrid->master_proc){
        if ( !boost::filesystem::exists(_outputDir) ){
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
int OUTPUT_MANAGER::returnDomainIfProbeIsInList(emProbe *newProbe){
    int pos = 0;
    for (std::vector<emProbe*>::iterator it = myEMProbes.begin(); it != myEMProbes.end(); it++){
        if ((*it)->compareProbes(newProbe)){
            return pos;
        }
        pos++;
    }
    return -1;

}
int OUTPUT_MANAGER::returnDomainIDIfDomainIsInList(outDomain *newDomain){
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
        req.itime = getIntTime(ttime);
        req.type = type;
        req.target = target;
        req.domain = domain;
        tempList.push_back(req);
    }
    tempList.sort(requestCompTime);
    reqList.merge(tempList, requestCompTime);
}

void OUTPUT_MANAGER::addEBField(outDomain* _domain, double startTime, double frequency, double endTime, whichFieldOut whichOut){
    if (!(checkGrid() && checkEMField()))
        return;

    int domainID = 0;
    if (_domain != NULL){
        domainID = returnDomainIDIfDomainIsInList(_domain);
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


// EM Probe ///////////////////////////////////////

void OUTPUT_MANAGER::addEBFieldProbe(emProbe* Probe, double startTime, double frequency, double endTime){
    if (!(checkGrid() && checkEMField()))
        return;
    int domainID = returnDomainIfProbeIsInList(Probe);
    if (domainID < 0){
        myEMProbes.push_back(Probe);
        isThereEMProbe = true;
        domainID = myEMProbes.size() - 1;
    }
    addRequestToList(requestList, OUT_EB_PROBE, 0, domainID, startTime, frequency, endTime);
}

void OUTPUT_MANAGER::addEBFieldProbeFrom(emProbe* Probe, double startTime, double frequency){
    addEBFieldProbe(Probe, startTime, frequency, mygrid->getTotalTime());
}

void OUTPUT_MANAGER::addEBFieldProbeAt(emProbe* Probe, double atTime){
    addEBFieldProbe(Probe, atTime, 1.0, atTime);
}

void OUTPUT_MANAGER::addEBFieldProbeFromTo(emProbe* Probe, double startTime, double frequency, double endTime){
    addEBFieldProbe(Probe, startTime, frequency, endTime);
}

// ++++++++++++++++++++++++++++     species density

void OUTPUT_MANAGER::addSpeciesDensity(outDomain* domain_in, std::string name, double startTime, double frequency, double endTime){
    if (!(checkGrid() && checkSpecies() && checkCurrent()))
        return;
    int specNum = findSpecName(name);
    if (specNum < 0)
        return;
    int domainID = 0;
    if (domain_in != NULL){
        domainID = returnDomainIDIfDomainIsInList(domain_in);
        if (domainID < 0){
            myDomains.push_back(domain_in);
            domainID = myDomains.size() - 1;
        }
    }
    addRequestToList(requestList, OUT_SPEC_DENSITY, specNum, domainID, startTime, frequency, endTime);
}

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

// ++++++++++++++++++++++++++++     current

void OUTPUT_MANAGER::addCurrent(outDomain* domain_in, double startTime, double frequency, double endTime){
    if (!(checkGrid() && checkCurrent()))
        return;
    int domainID = 0;

    if (domain_in != NULL){

        domainID = returnDomainIDIfDomainIsInList(domain_in);
        if (domainID < 0){
            myDomains.push_back(domain_in);
            domainID = myDomains.size() - 1;
        }
    }
    addRequestToList(requestList, OUT_CURRENT, 0, domainID, startTime, frequency, endTime);
}

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

// ++++++++++++++++++++++++++++     species binary

void OUTPUT_MANAGER::addSpeciesPhaseSpace(outDomain* domain_in, std::string name, double startTime, double frequency, double endTime){
    if (!(checkGrid() && checkSpecies()))
        return;
    int specNum = findSpecName(name);
    if (specNum < 0)
        return;
    int domainID = 0;

    if (domain_in != NULL){
        domainID = returnDomainIDIfDomainIsInList(domain_in);
        if (domainID < 0){
            myDomains.push_back(domain_in);
            domainID = myDomains.size() - 1;
        }
    }
    addRequestToList(requestList, OUT_SPEC_PHASE_SPACE, specNum, domainID, startTime, frequency, endTime);
}

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

// ++++++++++++++++++++++++++++     diag

void OUTPUT_MANAGER::addDiag(double startTime, double frequency, double endTime){
    if (!(checkGrid() && checkCurrent() && checkEMField()))
        return;
    addRequestToList(requestList, OUT_DIAG, 0, 0, startTime, frequency, endTime);
    isThereDiag = true;
}

void OUTPUT_MANAGER::addDiagFrom(double startTime, double frequency){
    addDiag(startTime, frequency, mygrid->getTotalTime());
}

void OUTPUT_MANAGER::addDiagAt(double atTime){
    addDiag(atTime, 1.0, atTime);
}

void OUTPUT_MANAGER::addDiagFromTo(double startTime, double frequency, double endTime){
    addDiag(startTime, frequency, endTime);
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


void OUTPUT_MANAGER::writeEMFieldBinary(std::string fileName, request req){
    int Ncomp = myfield->getNcomp();
    int uniqueN[3];
    uniqueN[0] = mygrid->uniquePoints[0];
    uniqueN[1] = mygrid->uniquePoints[1];
    uniqueN[2] = mygrid->uniquePoints[2];
    MPI_Offset disp = 0;
    int small_header = 6 * sizeof(int);
    int big_header = 3 * sizeof(int) + 3 * sizeof(int)
        + 2 * (uniqueN[0] + uniqueN[1] + uniqueN[2])*sizeof(double)
        + sizeof(int) + myfield->getNcomp() * 3 * sizeof(int);

    MPI_File thefile;
    MPI_Status status;

    char* nomefile = new char[fileName.size() + 1];

    nomefile[fileName.size()] = 0;
    sprintf(nomefile, "%s", fileName.c_str());

    MPI_File_open(MPI_COMM_WORLD, nomefile, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &thefile);

    //+++++++++++ FILE HEADER  +++++++++++++++++++++
    // We are not using the mygrid->master_proc to write them since it would require a double MPI_File_set_view
    // in case it is not the #0. Doing a double MPI_File_set_view from the same task to the same file is not guaranteed to work.
    if (mygrid->myid == 0){
        MPI_File_set_view(thefile, 0, MPI_FLOAT, MPI_FLOAT, (char *) "native", MPI_INFO_NULL);
        MPI_File_write(thefile, uniqueN, 3, MPI_INT, &status);
        MPI_File_write(thefile, mygrid->rnproc, 3, MPI_INT, &status);
        MPI_File_write(thefile, &Ncomp, 1, MPI_INT, &status);
        for (int c = 0; c < Ncomp; c++){
            integer_or_halfinteger crd = myfield->getCompCoords(c);
            int tp[3] = { (int)crd.x, (int)crd.y, (int)crd.z };
            MPI_File_write(thefile, tp, 3, MPI_INT, &status);
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
    if (disp < 0)
    {
        std::cout << "a problem occurred when trying to mpi_file_set_view in writeEMFieldBinary" << std::endl;
        std::cout << "myrank=" << mygrid->myid << " disp=" << disp << std::endl;
        exit(33);
    }
    if (mygrid->myid != 0){
        MPI_File_set_view(thefile, disp, MPI_FLOAT, MPI_FLOAT, (char *) "native", MPI_INFO_NULL);
    }


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
    delete[] nomefile;
    //////////////////////////// END of collective binary file write
}



#if defined(USE_HDF5)

void OUTPUT_MANAGER::writeEMFieldBinaryHDF5(std::string fileName, request req){
    int dimensionality=mygrid->accesso.dimensions;
    int Ncomp = myfield->getNcomp();

    MPI_Info info  = MPI_INFO_NULL;
    char nomi[6][3]={"Ex", "Ey", "Ez", "Bx", "By", "Bz"};
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
    if(dimensionality==1){
        dimsf[0] = mygrid->uniquePoints[0];
    }
    else if (dimensionality==2){
        dimsf[0] = mygrid->uniquePoints[1];
        dimsf[1] = mygrid->uniquePoints[0];
    }
    else{
        dimsf[0] = mygrid->uniquePoints[2];
        dimsf[1] = mygrid->uniquePoints[1];
        dimsf[2] = mygrid->uniquePoints[0];

    }
    for(int i=0;i<6;i++)
        filespace[i] = H5Screate_simple(dimensionality, dimsf, NULL);

    /*
         * Create the dataset with default properties and close filespace.
         */
    for(int i=0;i<6;i++){
        dset_id[i] = H5Dcreate1(file_id, nomi[i], H5T_NATIVE_FLOAT, filespace[i],
            H5P_DEFAULT);
        H5Sclose(filespace[i]);
    }
    /*
            * Each process defines dataset in memory and writes it to the hyperslab
            * in the file.
            */
    if(dimensionality==1){
        offset[0] = mygrid->rproc_imin[0][mygrid->rmyid[0]];
        count[0] = mygrid->uniquePointsloc[0];
    }
    else if (dimensionality==2){
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
    for(int i=0;i<6;i++){
        memspace[i] = H5Screate_simple(dimensionality, count, NULL);
        filespace[i] = H5Dget_space(dset_id[i]);
        H5Sselect_hyperslab(filespace[i], H5S_SELECT_SET, offset, NULL, count, NULL);
    }
    //+++++++++++ Start CPU Field Values  +++++++++++++++++++++
    {
        float *Ex,*Ey,*Ez,*Bx,*By,*Bz;
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

    for(int i=0;i<6;i++){
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


void OUTPUT_MANAGER::writeEBFieldDomain(std::string fileName, request req){
    int Ncomp = 3;//myfield->getNcomp();
    int offset = 0;
    if (req.type == OUT_E_FIELD)
        offset = 0;
    if (req.type == OUT_B_FIELD)
        offset = 3;
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


    shouldIWrite = isInMyDomain(rr);

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
            MPI_File_write(thefile, itodo, 8, MPI_INT, &status);

            float *fcir[3];
            for (int c = 0; c < 3; c++){
                fcir[c] = new float[uniqueN[c]];
                for (int m = 0; m < uniqueN[c]; m++){
                    fcir[c][m] = (float)mygrid->cir[c][m];
                }
            }
            for (int c = 0; c < 3; c++){
                if (remains[c])
                    MPI_File_write(thefile, fcir[c], uniqueN[c], MPI_FLOAT, &status);
                else
                    MPI_File_write(thefile, &fcir[c][globalri[c]], 1, MPI_FLOAT, &status);
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
            MPI_File_write(thefile, itodo, 6, MPI_INT, &status);
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
                            todo[c + i*Ncomp + j*Nx*Ncomp + k*Ny*Nx*Ncomp] = (float)myfield->VEB(c + offset, ii, jj, kk);
                    }
                }
            }
            MPI_File_write(thefile, todo, size, MPI_FLOAT, &status);
            delete[]todo;
        }
        MPI_File_close(&thefile);
        delete[] nomefile;
        delete[] totUniquePoints;
    }
    MPI_Comm_free(&sliceCommunicator);
}

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

    MPI_File_write(thefile, itodo, 8, MPI_INT, &status);
    for (int c = 0; c < 3; c++){
        MPI_File_write(thefile, fcir[c], uniqueN[c], MPI_FLOAT, &status);
    }

    for (int c = 0; c < 3; c++){
        delete[] fcir[c];
    }

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
    MPI_File_write(thefile, itodo, 6, MPI_INT, &status);
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

    int origin[3];
    int size = Ncomp*uniqueLocN[0] * uniqueLocN[1] * uniqueLocN[2];
    float *todo;
    todo = new float[size];
    int ri[3], globalri[3];
    nearestInt(myDomains[req.domain]->coordinates, ri, globalri);

    setLocalOutputOffset(origin, locimin, ri, remains);
    prepareFloatField(todo, uniqueLocN, origin, req);
    MPI_File_write(thefile, todo, size, MPI_FLOAT, &status);
    delete[]todo;
}

void OUTPUT_MANAGER::writeGridFieldSubDomain(std::string fileName, request req){
    int Ncomp = 3;
    if ((req.type == OUT_E_FIELD) || (req.type == OUT_B_FIELD))
        Ncomp = 3;
    else if (req.type == OUT_SPEC_DENSITY)
        Ncomp = 1;

    int *totUniquePoints;
    int isInMyHyperplane = false, shouldIWrite = false;
    int uniqueN[3], uniqueLocN[3], slice_rNproc[3];
    int remains[3];
    int imin[3], imax[3], locimin[3], locimax[3];

    setAndCheckRemains(remains, myDomains[req.domain]->remainingCoord);
    findIntGlobalBoundaries(myDomains[req.domain]->rmin, myDomains[req.domain]->rmax, imin, imax);
    findIntLocalBoundaries(myDomains[req.domain]->rmin, myDomains[req.domain]->rmax, locimin, locimax);
    findNumberOfProc(slice_rNproc, imin, imax, remains);

    findGlobalSubdomainUniquePoints(uniqueN, imin, imax, remains);
    findLocalSubdomainUniquePoints(uniqueLocN, locimin, locimax, remains);

    isInMyHyperplane = isInMyDomain(myDomains[req.domain]->coordinates);

    MPI_Comm sliceCommunicator;
    int mySliceID, sliceNProc;
    MPI_Cart_sub(mygrid->cart_comm, remains, &sliceCommunicator);
    MPI_Comm_rank(sliceCommunicator, &mySliceID);
    MPI_Comm_size(sliceCommunicator, &sliceNProc);
    MPI_Allreduce(MPI_IN_PLACE, &isInMyHyperplane, 1, MPI_INT, MPI_LOR, sliceCommunicator);
    if (isInMyHyperplane)
        shouldIWrite = amIInTheSubDomain(req);

    MPI_Comm outputCommunicator;
    MPI_Comm_split(mygrid->cart_comm, shouldIWrite, 0, &outputCommunicator);
    int myOutputID, outputNProc;
    MPI_Comm_rank(outputCommunicator, &myOutputID);
    MPI_Comm_size(outputCommunicator, &outputNProc);

    totUniquePoints = new int[outputNProc];
    totUniquePoints[myOutputID] = uniqueLocN[0] * uniqueLocN[1] * uniqueLocN[2];
    MPI_Allgather(MPI_IN_PLACE, 1, MPI_INT, totUniquePoints, 1, MPI_INT, outputCommunicator);

    MPI_Offset disp = 0;
    const int small_header = 6 * sizeof(int);
    const int big_header = (1 + 3 + 3 + 1)*sizeof(int) + (uniqueN[0] + uniqueN[1] + uniqueN[2])*sizeof(float);

    MPI_File thefile;
    //MPI_Status status;
    char* nomefile = new char[fileName.size() + 1];
    nomefile[fileName.size()] = 0;
    sprintf(nomefile, "%s", fileName.c_str());

    if (shouldIWrite){
        MPI_File_open(outputCommunicator, nomefile, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &thefile);

        //+++++++++++ FILE HEADER  +++++++++++++++++++++
        if (myOutputID == 0){
            MPI_File_set_view(thefile, 0, MPI_FLOAT, MPI_FLOAT, (char *) "native", MPI_INFO_NULL);
            writeBigHeader(thefile, uniqueN, imin, slice_rNproc, Ncomp);
        }
        //*********** END HEADER *****************


        //*********** SET VIEW *****************
        if (myOutputID != 0){
            findDispForSetView(&disp, myOutputID, totUniquePoints, big_header, small_header, Ncomp);
            MPI_File_set_view(thefile, disp, MPI_FLOAT, MPI_FLOAT, (char *) "native", MPI_INFO_NULL);
        }

        //+++++++++++ Write CPU HEADER  +++++++++++++++++++++
        {
            writeSmallHeader(thefile, uniqueLocN, imin, remains);
        }
        //+++++++++++ Write CPU Field Values  +++++++++++++++++++++
        {
            writeCPUFieldValues(thefile, uniqueLocN, locimin, remains, req);

        }
        MPI_File_close(&thefile);
        delete[] nomefile;
        delete[] totUniquePoints;
    }
    MPI_Comm_free(&sliceCommunicator);
    MPI_Comm_free(&outputCommunicator);
}


void OUTPUT_MANAGER::callEMFieldDomain(request req){

    if (req.type == OUT_E_FIELD){
        std::string nameBin = composeOutputName(outputDir, "E_FIELD", myDomains[req.domain]->name, "", req.domain, req.dtime, ".bin");
        if (!myDomains[req.domain]->subselection)
            writeEBFieldDomain(nameBin, req);
        else
            writeGridFieldSubDomain(nameBin, req);
    }
    else if (req.type == OUT_B_FIELD){
        std::string nameBin = composeOutputName(outputDir, "B_FIELD", myDomains[req.domain]->name, "", req.domain, req.dtime, ".bin");
        if (!myDomains[req.domain]->subselection)
            writeEBFieldDomain(nameBin, req);
        else
            writeGridFieldSubDomain(nameBin, req);
    }

}

void OUTPUT_MANAGER::interpEB(double pos[3], double E[3], double B[3]){
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
        for (int c = 0; c < mygrid->accesso.dimensions; c++){
            mycsi[c] = mygrid->unStretchGrid(pos[c], c);
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
        for (int c = 0; c < mygrid->accesso.dimensions; c++){
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
    }

    E[0] = E[1] = E[2] = B[0] = B[1] = B[2] = 0;

    switch (mygrid->accesso.dimensions)
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
        if (mygrid->accesso.dimensions < 2 || (rr[1] >= mygrid->rminloc[1] && rr[1] < mygrid->rmaxloc[1])){
            if (mygrid->accesso.dimensions < 3 || (rr[2] >= mygrid->rminloc[2] && rr[2] < mygrid->rmaxloc[2])){
                interpEB(rr, EE, BB);
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




void OUTPUT_MANAGER::writeSpecDensity(std::string fileName, request req){
    int Ncomp = 1;//myfield->getNcomp();
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

    shouldIWrite = isInMyDomain(rr);

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
            MPI_File_write(thefile, itodo, 8, MPI_INT, &status);

            float *fcir[3];
            for (int c = 0; c < 3; c++){
                fcir[c] = new float[uniqueN[c]];
                for (int m = 0; m < uniqueN[c]; m++){
                    fcir[c][m] = (float)mygrid->cir[c][m];
                }
            }
            for (int c = 0; c < 3; c++){
                //if(remains[c])
                MPI_File_write(thefile, fcir[c], uniqueN[c], MPI_FLOAT, &status);
                //else
                //  MPI_File_write(thefile, &fcir[c][globalri[c]], 1, MPI_FLOAT, &status);
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
            MPI_File_write(thefile, itodo, 6, MPI_INT, &status);
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
                            todo[c + i*Ncomp + j*Nx*Ncomp + k*Ny*Nx*Ncomp] = mycurrent->density(i, j, k);
                    }
                }
            }
            MPI_File_write(thefile, todo, size, MPI_FLOAT, &status);
            delete[]todo;
        }
        MPI_File_close(&thefile);
        delete[] nomefile;
        delete[] totUniquePoints;
    }
    MPI_Comm_free(&sliceCommunicator);
}

void OUTPUT_MANAGER::writeSpecDensitySubDomain(std::string fileName, request req){
    int Ncomp = 1;//myfield->getNcomp();
    int *totUniquePoints;
    int shouldIWrite = false;
    int isInMyHyperplane = false;
    int uniqueN[3], uniqueLocN[3], slice_rNproc[3];
    double rr[3] = { myDomains[req.domain]->coordinates[0], myDomains[req.domain]->coordinates[1], myDomains[req.domain]->coordinates[2] };
    int ri[3];
    int remains[3] = { myDomains[req.domain]->remainingCoord[0], myDomains[req.domain]->remainingCoord[1], myDomains[req.domain]->remainingCoord[2] };
    int imin[3], imax[3], locimin[3], locimax[3];

    findIntGlobalBoundaries(myDomains[req.domain]->rmin, myDomains[req.domain]->rmax, imin, imax);
    findIntLocalBoundaries(myDomains[req.domain]->rmin, myDomains[req.domain]->rmax, locimin, locimax);
    findNumberOfProc(slice_rNproc, imin, imax, remains);

    for (int c = 0; c < 3; c++){
        if (remains[c]){
            if (imax[c] < (mygrid->uniquePoints[c]))
                uniqueN[c] = imax[c] - imin[c] + 1;
            else
                uniqueN[c] = mygrid->uniquePoints[c] - imin[c];
        }
        else{
            uniqueN[c] = 1;
        }
    }
    for (int c = 0; c < 3; c++){
        if (remains[c]){
            if (locimax[c] < (mygrid->uniquePointsloc[c]))
                uniqueLocN[c] = locimax[c] - locimin[c] + 1;
            else
                uniqueLocN[c] = mygrid->uniquePointsloc[c] - locimin[c];
        }
        else{
            uniqueLocN[c] = 1;
        }
    }

    isInMyHyperplane = isInMyDomain(rr);

    MPI_Comm sliceCommunicator;
    int mySliceID, sliceNProc;

    int dimension;
    MPI_Cart_sub(mygrid->cart_comm, remains, &sliceCommunicator);
    MPI_Comm_rank(sliceCommunicator, &mySliceID);
    MPI_Comm_size(sliceCommunicator, &sliceNProc);
    MPI_Cartdim_get(sliceCommunicator, &dimension);
    MPI_Allreduce(MPI_IN_PLACE, &isInMyHyperplane, 1, MPI_INT, MPI_LOR, sliceCommunicator);
    if (isInMyHyperplane)
        shouldIWrite = amIInTheSubDomain(req);

    MPI_Comm outputCommunicator;
    MPI_Comm_split(mygrid->cart_comm, shouldIWrite, 0, &outputCommunicator);
    int myOutputID, outputNProc;
    MPI_Comm_rank(outputCommunicator, &myOutputID);
    MPI_Comm_size(outputCommunicator, &outputNProc);

    totUniquePoints = new int[outputNProc];
    totUniquePoints[myOutputID] = uniqueLocN[0] * uniqueLocN[1] * uniqueLocN[2];


    MPI_Allgather(MPI_IN_PLACE, 1, MPI_INT, totUniquePoints, 1, MPI_INT, outputCommunicator);


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
        MPI_File_open(outputCommunicator, nomefile, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &thefile);
        int globalri[3];
        nearestInt(rr, ri, globalri);
        //+++++++++++ FILE HEADER  +++++++++++++++++++++
        if (myOutputID == 0){
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
            MPI_File_write(thefile, itodo, 8, MPI_INT, &status);

            float *fcir[3];
            for (int c = 0; c < 3; c++){
                fcir[c] = new float[uniqueN[c]];
                for (int m = 0; m < uniqueN[c]; m++){
                    fcir[c][m] = (float)mygrid->cir[c][m + imin[c]];
                }
            }
            for (int c = 0; c < 3; c++){
                MPI_File_write(thefile, fcir[c], uniqueN[c], MPI_FLOAT, &status);
            }
            for (int c = 0; c < 3; c++){
                delete[] fcir[c];
            }

        }
        //*********** END HEADER *****************

        disp = big_header;
        for (int rank = 0; rank < myOutputID; rank++)
            disp += small_header + totUniquePoints[rank] * sizeof(float)*Ncomp;
        if (disp < 0){
            std::cout << "a problem occurred when trying to mpi_file_set_view in writeEMFieldBinary" << std::endl;
            std::cout << "myrank=" << mygrid->myid << " disp=" << disp << std::endl;
            exit(33);
        }
        if (myOutputID != 0){
            MPI_File_set_view(thefile, disp, MPI_FLOAT, MPI_FLOAT, (char *) "native", MPI_INFO_NULL);
        }

        //+++++++++++ Start CPU HEADER  +++++++++++++++++++++
        {
            int itodo[6];
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
            MPI_File_write(thefile, itodo, 6, MPI_INT, &status);
        }
        //+++++++++++ Start CPU Field Values  +++++++++++++++++++++
            {
                float *todo;
                int NN[3], Nx, Ny, Nz, origin[3];
                for (int c = 0; c < 3; c++){
                    if (remains[c]){
                        NN[c] = uniqueLocN[c];
                        origin[c] = locimin[c];
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
                                todo[c + i*Ncomp + j*Nx*Ncomp + k*Ny*Nx*Ncomp] = (float)mycurrent->density(ii, jj, kk);
                        }
                    }
                }
                MPI_File_write(thefile, todo, size, MPI_FLOAT, &status);
                delete[]todo;
            }
        MPI_File_close(&thefile);
        delete[] nomefile;
        delete[] totUniquePoints;
    }
    MPI_Comm_free(&sliceCommunicator);
    MPI_Comm_free(&outputCommunicator);
}

void OUTPUT_MANAGER::callSpecDensity(request req){
    mycurrent->eraseDensity();
    myspecies[req.target]->density_deposition_standard(mycurrent);
    mycurrent->pbc();

    //    if (mygrid->myid == mygrid->master_proc){
    //        std::string nameMap = composeOutputName(outputDir, "DENS", name, req.dtime, ".map");
    //        std::ofstream of1;
    //        of1.open(nameMap.c_str());
    //        writeSpecDensityMap(of1, req);
    //        of1.close();
    //    }

    std::string nameBin = composeOutputName(outputDir, "DENS", myspecies[req.target]->name, myDomains[req.domain]->name, req.domain, req.dtime, ".bin");

    if (!myDomains[req.domain]->subselection)
        writeSpecDensity(nameBin, req);
    else
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

    shouldIWrite = isInMyDomain(rr);

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
            MPI_File_write(thefile, itodo, 8, MPI_INT, &status);

            float *fcir[3];
            for (int c = 0; c < 3; c++){
                fcir[c] = new float[uniqueN[c]];
                for (int m = 0; m < uniqueN[c]; m++){
                    fcir[c][m] = (float)mygrid->cir[c][m];
                }
            }
            for (int c = 0; c < 3; c++){
                if (remains[c])
                    MPI_File_write(thefile, fcir[c], uniqueN[c], MPI_FLOAT, &status);
                else
                    MPI_File_write(thefile, &fcir[c][globalri[c]], 1, MPI_FLOAT, &status);
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
            MPI_File_write(thefile, itodo, 6, MPI_INT, &status);
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
            MPI_File_write(thefile, todo, size, MPI_FLOAT, &status);
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

    writeCurrent(nameBin, req);

}

void OUTPUT_MANAGER::writeSpecPhaseSpace(std::string fileName, request req){

    SPECIE* spec = myspecies[req.target];
    int* NfloatLoc = new int[mygrid->nproc];
    int outputNComp, NCompFloat;
    bool flagMarker = spec->amIWithMarker();
    if (flagMarker){
        NCompFloat = spec->Ncomp - 1;
        outputNComp = spec->Ncomp + 1;
        NfloatLoc[mygrid->myid] = spec->Np*outputNComp;
    }
    else{
        NCompFloat = spec->Ncomp;
        outputNComp = spec->Ncomp;
        NfloatLoc[mygrid->myid] = spec->Np*outputNComp;
    }
    MPI_Allgather(MPI_IN_PLACE, 1, MPI_INT, NfloatLoc, 1, MPI_INT, MPI_COMM_WORLD);

    MPI_Offset disp = 0;
    for (int pp = 0; pp < mygrid->myid; pp++)
        disp += (MPI_Offset)(NfloatLoc[pp] * sizeof(float));
    MPI_File thefile;
    //printf("writePHACESPACE NCompFloat = %i     outputNComp = %i\n", NCompFloat, outputNComp);
    //printf("NfloatLoc[%i] = %i \n", mygrid->myid,NfloatLoc[mygrid->myid] );
    //printf("flagMarker=%i        disp=%i\n", flagMarker, disp);

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
    buf = new float[dimensione*outputNComp];
    passaggi = spec->Np / dimensione;
    resto = spec->Np % dimensione;
    //printf("writePHACESPACE NCompFloat = %i     outputNComp = %i\n", NCompFloat, outputNComp);

    for (int i = 0; i < passaggi; i++){
        for (int p = 0; p < dimensione; p++){
            int c;
            for (c = 0; c < NCompFloat; c++){
                buf[c + p*outputNComp] = (float)spec->ru(c, p + dimensione*i);
            }
            if (flagMarker){
                buf[c + p*outputNComp] = *((float*)(&(spec->marker(p + dimensione*i))));
                buf[c + 1 + p*outputNComp] = *((float*)(&(spec->marker(p + dimensione*i))) + 1);
                //*((double*)(&(buf[c + p*outputNComp]))) = spec->pettorale(p+ dimensione*i;
            }
        }
        MPI_File_write(thefile, buf, dimensione*outputNComp, MPI_FLOAT, &status);
    }
    for (int p = 0; p < resto; p++){
        int c;
        for (c = 0; c < NCompFloat; c++){
            buf[c + p*outputNComp] = (float)spec->ru(c, p + dimensione*passaggi);
        }
        if (flagMarker){
            buf[c + p*outputNComp] = *((float*)(&(spec->marker(p + dimensione*passaggi))));
            buf[c + 1 + p*outputNComp] = *((float*)(&(spec->marker(p + dimensione*passaggi))) + 1);
            //*((double*)(&(buf[c + p*outputNComp]))) = spec->pettorale(p+ dimensione*i;
        }
    }
    MPI_File_write(thefile, buf, resto*outputNComp, MPI_FLOAT, &status);
    MPI_File_close(&thefile);
    delete[]buf;
    delete[] NfloatLoc;
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
            if (mygrid->accesso.dimensions < 2 || (rmax[1] >= rr[1] && rmin[1] < rr[1])){
                if (mygrid->accesso.dimensions < 3 || (rmax[2] >= rr[2] && rmin[2] < rr[2])){
                    counter++;
                }
            }
        }
    }
    return counter;
}

void OUTPUT_MANAGER::writeSpecPhaseSpaceSubDomain(std::string fileName, request req){

    SPECIE* spec = myspecies[req.target];
    int shouldIWrite = false;
    shouldIWrite = amIInTheSubDomain(req);

    MPI_Comm outputCommunicator;
    MPI_Comm_split(mygrid->cart_comm, shouldIWrite, 0, &outputCommunicator);
    int myOutputID, outputNProc;
    MPI_Comm_rank(outputCommunicator, &myOutputID);
    MPI_Comm_size(outputCommunicator, &outputNProc);

    int outputNPart = findNumberOfParticlesInSubdomain(req);
    int* NfloatLoc = new int[outputNProc];
    NfloatLoc[myOutputID] = outputNPart*spec->Ncomp;

    MPI_Allgather(MPI_IN_PLACE, 1, MPI_INT, NfloatLoc, 1, MPI_INT, outputCommunicator);

    MPI_Offset disp = 0;
    for (int pp = 0; pp < myOutputID; pp++)
        disp += (MPI_Offset)(NfloatLoc[pp] * sizeof(float));
    MPI_File thefile;
    MPI_Status status;

    char *nomefile = new char[fileName.size() + 1];
    nomefile[fileName.size()] = 0;
    sprintf(nomefile, "%s", fileName.c_str());

    if (shouldIWrite){
        MPI_File_open(outputCommunicator, nomefile, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &thefile);
        MPI_File_set_view(thefile, disp, MPI_FLOAT, MPI_FLOAT, (char *) "native", MPI_INFO_NULL);

        float *buf;
        int dimensione = 100000;
        buf = new float[dimensione*spec->Ncomp];

        double rr[3];
        int counter = 0;
        double rmin[3], rmax[3];
        for (int c = 0; c < 3; c++){
            rmin[c] = myDomains[req.domain]->rmin[c];
            rmax[c] = myDomains[req.domain]->rmax[c];
        }

        for (int p = 0; p < spec->Np; p++){
            rr[0] = spec->r0(p);
            rr[1] = spec->r1(p);
            rr[2] = spec->r2(p);
            if (rmax[0] >= rr[0] && rmin[0] < rr[0]){
                if (mygrid->accesso.dimensions < 2 || (rmax[1] >= rr[1] && rmin[1] < rr[1])){
                    if (mygrid->accesso.dimensions < 3 || (rmax[2] >= rr[2] && rmin[2] < rr[2])){
                        for (int c = 0; c < spec->Ncomp; c++){
                            buf[c + counter*spec->Ncomp] = (float)spec->ru(c, p);
                        }
                        counter++;
                    }
                }
            }
            if (counter == dimensione){
                MPI_File_write(thefile, buf, counter*spec->Ncomp, MPI_FLOAT, &status);
                counter = 0;
            }
        }
        if (counter > 0){
            MPI_File_write(thefile, buf, counter*spec->Ncomp, MPI_FLOAT, &status);
        }
        MPI_File_close(&thefile);
        delete[]buf;
    }
    MPI_Comm_free(&outputCommunicator);
    delete[] NfloatLoc;
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

