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


#include "utilities.h"


void UTILITIES::moveWindow(GRID* _mygrid, EM_FIELD* _myfield, std::vector<SPECIE*> _myspecies) {
  _mygrid->moveWindow();
  _myfield->moveWindow();
  for (std::vector<SPECIE*>::iterator spec_iterator = _myspecies.begin(); spec_iterator != _myspecies.end(); spec_iterator++) {
    (*spec_iterator)->move_window();
  }
}

void UTILITIES::considerRestartFromDump(GRID* grid, EM_FIELD* myfield, std::vector<SPECIE*> species) {
  if (grid->dumpControl.doRestart) {
    grid->dumpControl.currentDumpID = grid->dumpControl.restartFromDump;
    UTILITIES::restartFromDump(grid, myfield, species);
  }
}

void UTILITIES::restartFromDump(GRID* mygrid, EM_FIELD* myfield, std::vector<SPECIE*> species) {
  int dumpID = mygrid->dumpControl.currentDumpID;
  std::ifstream dumpFile;
  MPI_Barrier(MPI_COMM_WORLD);

  if (mygrid->myid == mygrid->master_proc) {
    time_t timer;
    std::time(&timer);  /* get current time; same as: timer = time(NULL)  */

    struct tm * now = localtime(&timer);

    printf("   restart from DUMP #%i ... %2.2i:%2.2i:%2.2i\n", (dumpID), now->tm_hour, now->tm_min, now->tm_sec);
    fflush(stdout);
  }
  dumpFile.open(mygrid->composeDumpFileName(dumpID).c_str());
  if (dumpFile.good()) {
    mygrid->reloadDump(dumpFile);
    myfield->reloadDump(dumpFile);

    for (std::vector<SPECIE*>::iterator spec_iterator = species.begin(); spec_iterator != species.end(); spec_iterator++) {
      (*spec_iterator)->reloadBigBufferDump(dumpFile);
    }
    dumpFile.close();
    dumpID++;
  }
  MPI_Barrier(MPI_COMM_WORLD);
  if (mygrid->myid == mygrid->master_proc) {
    time_t timer;
    std::time(&timer);  /* get current time; same as: timer = time(NULL)  */

    struct tm * now = localtime(&timer);

    printf("  ... DONE %2.2i:%2.2i:%2.2i\n", now->tm_hour, now->tm_min, now->tm_sec);
    fflush(stdout);
  }
  mygrid->dumpControl.currentDumpID = dumpID;
}

void UTILITIES::restartFromDump(int *_dumpID, GRID* mygrid, EM_FIELD* myfield, std::vector<SPECIE*> species) {
  int dumpID = mygrid->dumpControl.currentDumpID;
  std::ifstream dumpFile;
  MPI_Barrier(MPI_COMM_WORLD);

  if (mygrid->myid == mygrid->master_proc) {
    time_t timer;
    std::time(&timer);  /* get current time; same as: timer = time(NULL)  */

    struct tm * now = localtime(&timer);

    printf("   restart from DUMP #%i ... %2.2i:%2.2i:%2.2i\n", (dumpID), now->tm_hour, now->tm_min, now->tm_sec);
    fflush(stdout);
  }
  dumpFile.open(mygrid->composeDumpFileName(dumpID).c_str());
  if (dumpFile.good()) {
    mygrid->reloadDump(dumpFile);
    myfield->reloadDump(dumpFile);

    for (std::vector<SPECIE*>::iterator spec_iterator = species.begin(); spec_iterator != species.end(); spec_iterator++) {
      (*spec_iterator)->reloadBigBufferDump(dumpFile);
    }
    dumpFile.close();
    dumpID++;
  }
  MPI_Barrier(MPI_COMM_WORLD);
  if (mygrid->myid == mygrid->master_proc) {
    time_t timer;
    std::time(&timer);  /* get current time; same as: timer = time(NULL)  */

    struct tm * now = localtime(&timer);

    printf("  ... DONE %2.2i:%2.2i:%2.2i\n", now->tm_hour, now->tm_min, now->tm_sec);
    fflush(stdout);
  }
  _dumpID[0] = dumpID;
}

void UTILITIES::considerDumpForRestart(GRID* grid, EM_FIELD* myfield, std::vector<SPECIE*> species){
  if (grid->dumpControl.doDump) {
    if (grid->istep != 0 && !(grid->istep % ((int)(grid->dumpControl.dumpEvery / grid->dt)))) {
      UTILITIES::dumpFilesForRestart(grid, myfield, species);
    }
  }
}


void UTILITIES::dumpFilesForRestart(GRID* mygrid, EM_FIELD* myfield, std::vector<SPECIE*> species) {
  int dumpID = mygrid->dumpControl.currentDumpID;
  std::ofstream dumpFile;
  dumpFile.open(mygrid->composeDumpFileName(dumpID).c_str());
  mygrid->dump(dumpFile);
  myfield->dump(dumpFile);
  for (std::vector<SPECIE*>::iterator spec_iterator = species.begin(); spec_iterator != species.end(); spec_iterator++) {
    (*spec_iterator)->dumpBigBuffer(dumpFile);
  }
  dumpFile.close();
  dumpID++;
  mygrid->dumpControl.currentDumpID = dumpID;
  MPI_Barrier(MPI_COMM_WORLD);
  if (mygrid->myid == mygrid->master_proc) {
    printf("\t DUMP #%i done!\n", (dumpID - 1));
  }
}
void UTILITIES::dumpFilesForRestart(int *_dumpID, GRID* mygrid, EM_FIELD* myfield, std::vector<SPECIE*> species) {
  int dumpID = _dumpID[0];
  std::ofstream dumpFile;
  dumpFile.open(mygrid->composeDumpFileName(dumpID).c_str());
  mygrid->dump(dumpFile);
  myfield->dump(dumpFile);
  for (std::vector<SPECIE*>::iterator spec_iterator = species.begin(); spec_iterator != species.end(); spec_iterator++) {
    (*spec_iterator)->dumpBigBuffer(dumpFile);
  }
  dumpFile.close();
  dumpID++;
  _dumpID[0] = dumpID;
  MPI_Barrier(MPI_COMM_WORLD);
  if (mygrid->myid == mygrid->master_proc) {
    printf("\t DUMP #%i done!\n", (dumpID - 1));
  }
}



void UTILITIES::dumpDebugFilesForRestart(int *_dumpID, GRID* mygrid, EM_FIELD* myfield, std::vector<SPECIE*> species) {
  int dumpID = _dumpID[0];
  std::ofstream dumpFile;
  dumpFile.open(mygrid->composeDumpFileName(dumpID).c_str());
  mygrid->debugDump(dumpFile);
  //myfield->debugDump(dumpFile);
  for (std::vector<SPECIE*>::iterator spec_iterator = species.begin(); spec_iterator != species.end(); spec_iterator++) {
    (*spec_iterator)->debugDump(dumpFile);
  }
  dumpFile.close();
  dumpID++;
  _dumpID[0] = dumpID;
  MPI_Barrier(MPI_COMM_WORLD);
  if (mygrid->myid == mygrid->master_proc) {
    printf("\t DUMP #%i done!\n", (dumpID - 1));
  }
}

bool UTILITIES::doesFileExist(const char *fileName)
{
  std::ifstream infile(fileName);
  return infile.good();
}

void UTILITIES::exitWithError(int error) {
  MPI_Finalize();
  exit(error);
}

void UTILITIES::splitCommGetRankNproc(MPI_Comm parentComm, MPI_Comm *childComm, int color, int *rank, int *NProcs) {
  MPI_Comm_split(parentComm, color, 0, childComm);
  MPI_Comm_size(*childComm, NProcs);
  MPI_Comm_rank(*childComm, rank);
}

void UTILITIES::readAndAllocateSpheres(SPHERES &spheres, std::string filename, GRID &grid) {

  if (grid.myid == grid.master_proc) {
    std::ifstream myfile;
    myfile.open(filename.c_str());
    if (!myfile.good()) {
      std::cout << "     spheres file: \"" << filename << "\" does not exists" << std::endl;
    }
    int Nsph = 0;

    myfile.read((char*)&Nsph, sizeof(int));

    spheres.NSpheres = Nsph;
    spheres.coords = new float[spheres.NSpheres * 4];

    myfile.read((char*)&(spheres.fillingFactor), sizeof(float));
    myfile.read((char*)spheres.rmin, sizeof(float) * 3);
    myfile.read((char*)spheres.rmax, sizeof(float) * 3);
    myfile.read((char*)spheres.coords, sizeof(float)*spheres.NSpheres * 4);
    myfile.close();
  }


  MPI_Bcast(&(spheres.NSpheres), 1, MPI_INT, grid.master_proc, MPI_COMM_WORLD);

  if (grid.myid != grid.master_proc)
    spheres.coords = new float[spheres.NSpheres * 4];

  MPI_Bcast(&(spheres.fillingFactor), 1, MPI_FLOAT, grid.master_proc, MPI_COMM_WORLD);
  MPI_Bcast(spheres.rmin, 3, MPI_FLOAT, grid.master_proc, MPI_COMM_WORLD);
  MPI_Bcast(spheres.rmax, 3, MPI_FLOAT, grid.master_proc, MPI_COMM_WORLD);
  MPI_Bcast(spheres.coords, 4 * spheres.NSpheres, MPI_FLOAT, grid.master_proc, MPI_COMM_WORLD);

}

void UTILITIES::fromCoordsToSpheresCoords(double &x, double min, double max) {
  double box = max - min;
  if (x < min) {
    x += box * ((int)((max - x) / box));
  }
  if (x > max) {
    x -= box* ((int)((x - min) / box));
  }
}

bool UTILITIES::isSphereInside(SPHERES& spheres, int index, GRID &grid) {
  float y = spheres.coords[index * 4 + 1];
  float z = spheres.coords[index * 4 + 2];
  float r = spheres.coords[index * 4 + 3];

  double yl = grid.rminloc[1] - r;
  double zl = grid.rminloc[2] - r;

  double yr = grid.rmaxloc[1] + r;
  double zr = grid.rmaxloc[2] + r;

  UTILITIES::fromCoordsToSpheresCoords(yl, spheres.rmin[1], spheres.rmax[1]);
  UTILITIES::fromCoordsToSpheresCoords(zl, spheres.rmin[2], spheres.rmax[2]);
  UTILITIES::fromCoordsToSpheresCoords(yr, spheres.rmin[1], spheres.rmax[1]);
  UTILITIES::fromCoordsToSpheresCoords(zr, spheres.rmin[2], spheres.rmax[2]);

  bool chkX, chkY, chkZ;
  chkX = chkY = chkZ = true;

  if (grid.getDimensionality() >= 2) {
    if ((grid.rmaxloc[1] - grid.rminloc[1]) >= (spheres.rmax[1] - spheres.rmin[1]))
      chkY = true;
    else {
      if (yl <= yr)
        chkY = ((y >= yl) && (y <= yr));
      else
        chkY = ((y >= yr) || (y <= yl));
    }
  }

  if (grid.getDimensionality() == 3) {
    if ((grid.rmaxloc[2] - grid.rminloc[2]) >= (spheres.rmax[2] - spheres.rmin[2]))
      chkZ = true;
    else {

      if (zl <= zr)
        chkZ = ((z >= zl) && (z <= zr));
      else
        chkZ = ((z >= zr) || (z <= zl));
    }
  }

  return chkX && chkY && chkZ;

}

void UTILITIES::swapSpheres(SPHERES &spheres, int i, int j) {
  float dummyCoords[4];
  dummyCoords[0] = spheres.coords[i * 4];
  dummyCoords[1] = spheres.coords[i * 4 + 1];
  dummyCoords[2] = spheres.coords[i * 4 + 2];
  dummyCoords[3] = spheres.coords[i * 4 + 3];

  spheres.coords[i * 4] = spheres.coords[j * 4];
  spheres.coords[i * 4 + 1] = spheres.coords[j * 4 + 1];
  spheres.coords[i * 4 + 2] = spheres.coords[j * 4 + 2];
  spheres.coords[i * 4 + 3] = spheres.coords[j * 4 + 3];

  spheres.coords[j * 4] = dummyCoords[0];
  spheres.coords[j * 4 + 1] = dummyCoords[1];
  spheres.coords[j * 4 + 2] = dummyCoords[2];
  spheres.coords[j * 4 + 3] = dummyCoords[3];
}

void UTILITIES::selectSpheres(SPHERES &spheres, GRID &grid) {
  int counter = 0;
  for (int i = 0; i < spheres.NSpheres; i++) {
    if (UTILITIES::isSphereInside(spheres, i, grid)) {
      UTILITIES::swapSpheres(spheres, i, counter);
      counter++;
    }
  }
  spheres.NSpheres = counter;
  spheres.coords = (float*)realloc(spheres.coords, counter * 4 * sizeof(float));
}

void UTILITIES::readAndAllocateFFTplasma(FFTPLASMA &myfft, std::string filename, GRID &grid) {

  if (grid.myid == grid.master_proc) {
    std::ifstream myfile;
    myfile.open(filename.c_str());
    if (!myfile.good()) {
      std::cout << "     FFTplasma file: \"" << filename << "\" does not exists" << std::endl;
    }


    int intbuf;
    double doubuf;
    myfile.read((char*)&intbuf, sizeof(int));
    myfile.read((char*)&doubuf, sizeof(double));
    myfft.numcomp=intbuf;
    myfft.shift=doubuf;

    myfft.kx = new double[myfft.numcomp];
    myfft.ky = new double[myfft.numcomp];
    myfft.phi = new double[myfft.numcomp];
    myfft.cc = new double[myfft.numcomp];

    myfile.read((char*)(myfft.kx), sizeof(double)*myfft.numcomp);
    myfile.read((char*)(myfft.ky), sizeof(double)*myfft.numcomp);
    myfile.read((char*)(myfft.phi), sizeof(double)*myfft.numcomp);
    myfile.read((char*)(myfft.cc), sizeof(double)*myfft.numcomp);

    myfile.close();
  }


  MPI_Bcast(&(myfft.numcomp), 1, MPI_INT, grid.master_proc, MPI_COMM_WORLD);
  MPI_Bcast(&(myfft.shift), 1, MPI_DOUBLE, grid.master_proc, MPI_COMM_WORLD);

  if (grid.myid != grid.master_proc){
      myfft.kx = new double[myfft.numcomp];
      myfft.ky = new double[myfft.numcomp];
      myfft.phi = new double[myfft.numcomp];
      myfft.cc = new double[myfft.numcomp];

  }

    MPI_Bcast(myfft.kx, myfft.numcomp, MPI_DOUBLE, grid.master_proc, MPI_COMM_WORLD);
    MPI_Bcast(myfft.ky, myfft.numcomp, MPI_DOUBLE, grid.master_proc, MPI_COMM_WORLD);
    MPI_Bcast(myfft.phi, myfft.numcomp, MPI_DOUBLE, grid.master_proc, MPI_COMM_WORLD);
    MPI_Bcast(myfft.cc, myfft.numcomp, MPI_DOUBLE, grid.master_proc, MPI_COMM_WORLD);


}

void UTILITIES::allocateAccessibleKModes(GRIDmodes &gridModes, GRID &grid){

  double Lr[3];

  {
    int c;
    for(c=0; c<grid.getDimensionality();c++){
      Lr[c] = grid.rmax[c] - grid.rmin[c];
    }
    for(;c<3;c++){
      Lr[c] = 0;
    }
  }
  for(int c=0; c<3;c++){
    gridModes.Nk[c] = grid.uniquePoints[c];
    if(Lr[c]>0){
      gridModes.dk[c] = 2*M_PI/Lr[c];
    }
    else{
      gridModes.dk[c] = 0;
    }

    gridModes.kmin[c] = -(gridModes.Nk[c]/2)*gridModes.dk[c];
    gridModes.kmax[c] = gridModes.kmin[c] + (gridModes.Nk[c]-1)*gridModes.dk[c];
    gridModes.kvalues[c] = (double*)malloc(sizeof(double)*gridModes.Nk[c]);
    for(int i=0; i<gridModes.Nk[c]; i++){
      gridModes.kvalues[c][i] = gridModes.kmin[c] + i*gridModes.dk[c];
    }
  }
}

void UTILITIES::writeGridModes(GRIDmodes &gridModes, GRID &grid){
  std::stringstream message;
  char letter[]={'x','y', 'z'};
  for(int c=0; c<3; c++){
    message << " # k"<< letter[c] <<"x modes:" << std::endl;
    message << " Nk"<< letter[c] <<": " << gridModes.Nk[c] << std::endl;
    message << " dk"<< letter[c] <<": " << gridModes.dk[c] << std::endl;
    message << " k"<< letter[c] <<" range: [ " << gridModes.kmin[c] << " : " << gridModes.kmax[c] << " ] " << std::endl;
    for(int i=0; i<gridModes.Nk[c]; i++){
      message << "| " << gridModes.kvalues[c][i] << " |";
    }
    message << "\n ######################" << std::endl;
  }
  std::string fileName("gridModes.txt");
  grid.printInfoFile(fileName,message.str());
}

void UTILITIES::setKModesToBeInitialised(std::vector<KMODE> &myKModes, GRIDmodes &gridModes, double amplitude, double* centralK, double* sigmaK){
  my_rng_generator mt_rng;
  my_uniform_real_distribution randPhase(0, 2.0*M_PI);

  double kz, ky, kx;
  for(int k=0; k<gridModes.Nk[2]; k++){
    kz = gridModes.kvalues[2][k];
    for(int j=0; j<gridModes.Nk[1]; j++){
      ky = gridModes.kvalues[1][j];
        for(int i=0; i<gridModes.Nk[0]; i++){
          kx = gridModes.kvalues[0][i];
        double amp = 1;
        amp *= exp( -(kx-centralK[0])*(kx-centralK[0])/(sigmaK[0]*sigmaK[0]) );
        amp *= exp( -(ky-centralK[1])*(ky-centralK[1])/(sigmaK[1]*sigmaK[1]) );
        amp *= exp( -(kz-centralK[2])*(kz-centralK[2])/(sigmaK[2]*sigmaK[2]) );
        if(amp>1e-2){
          amp *= amplitude;
          KMODE newMode;
          newMode.amplitude=amp;
          newMode.k[0] = kx;
          newMode.k[1] = ky;
          newMode.k[2] = kz;
          newMode.phase = randPhase(mt_rng);
          myKModes.push_back(newMode);
        }

      }
    }
  }

}

void UTILITIES::newSetKModesToBeInitialised(LANGMUIRset &langmuirSet){
  my_rng_generator mt_rng;
  mt_rng.seed(54890);
  my_uniform_real_distribution randPhase(0, 2.0*M_PI);

  double kz, ky, kx;
  for(int k=0; k<langmuirSet.gridModes.Nk[2]; k++){
    kz = langmuirSet.gridModes.kvalues[2][k];
    for(int j=0; j<langmuirSet.gridModes.Nk[1]; j++){
      ky = langmuirSet.gridModes.kvalues[1][j];
        for(int i=0; i<langmuirSet.gridModes.Nk[0]; i++){
          kx = langmuirSet.gridModes.kvalues[0][i];
        double amp = 1;
        amp *= exp( -(kx-langmuirSet.centralK[0])*(kx-langmuirSet.centralK[0])/(langmuirSet.sigmaK[0]*langmuirSet.sigmaK[0]) );
        amp *= exp( -(ky-langmuirSet.centralK[1])*(ky-langmuirSet.centralK[1])/(langmuirSet.sigmaK[1]*langmuirSet.sigmaK[1]) );
        amp *= exp( -(kz-langmuirSet.centralK[2])*(kz-langmuirSet.centralK[2])/(langmuirSet.sigmaK[2]*langmuirSet.sigmaK[2]) );
        if(amp>1e-2){
          amp *= langmuirSet.amplitude;
          KMODE newMode;
          newMode.amplitude=amp;
          newMode.k[0] = kx;
          newMode.k[1] = ky;
          newMode.k[2] = kz;
          newMode.phase = randPhase(mt_rng);
          langmuirSet.myKModes.push_back(newMode);
        }

      }
    }
  }

}

void UTILITIES::exchangeKModesToBeInitialised(std::vector<KMODE> &myKModes, GRID &grid){
  int Nmodes=myKModes.size();
  double *buffer=new double[Nmodes*5];
  for(int i=0; i<Nmodes; i++){
    buffer[i*5+0] = myKModes[i].k[0];
    buffer[i*5+1] = myKModes[i].k[1];
    buffer[i*5+2] = myKModes[i].k[2];
    buffer[i*5+3] = myKModes[i].amplitude;
    buffer[i*5+4] = myKModes[i].phase;
  }
  MPI_Bcast(buffer, Nmodes*5, MPI_DOUBLE, grid.master_proc, MPI_COMM_WORLD);
  myKModes.clear();
  for(int i=0; i<Nmodes; i++){
    KMODE newMode;
    newMode.k[0]      = buffer[i*5+0];
    newMode.k[1]      = buffer[i*5+1];
    newMode.k[2]      = buffer[i*5+2];
    newMode.amplitude = buffer[i*5+3];
    newMode.phase     = buffer[i*5+4];
    myKModes.push_back(newMode);
  }
  delete[] buffer;

}

void UTILITIES::writeKModesToBeInitialised(std::vector<KMODE> &myKModes, GRID &grid){
  int Nmodes=myKModes.size();
  std::stringstream message;
  message << "#KModes are " << myKModes.size() << std::endl;
  message << "#" << std::setw(10) << " kx " << std::setw(11) << " ky " << std::setw(11) << " kz " << std::setw(11) << " amp " << std::setw(11) << " phase" << std::endl;
  for(int i=0; i<Nmodes; i++){
    message << std::setw(10) << myKModes[i].k[0] << " ";
    message << std::setw(10) << myKModes[i].k[1] << " ";
    message << std::setw(10) << myKModes[i].k[2] << " ";
    message << std::setw(10) << myKModes[i].amplitude << " ";
    message << std::setw(10) << myKModes[i].phase << std::endl;
  }
  std::string fileName("kmodes.txt");
  grid.printInfoFile(fileName,message.str());
}

void UTILITIES::writeLangmuirSetDataFile(LANGMUIRset &langmuirSet, GRID &grid){
  std::stringstream message;
  message << " langmuirSet Validity=" << langmuirSet.checkLangmuirSetValidity << std::endl;
  message << " amplitude=" << langmuirSet.amplitude << std::endl;
  message << " refTemperature=" << langmuirSet.refTemp << std::endl;
  message << " refDensity=" << langmuirSet.refDens << std::endl;
  message << " centralKx=" << langmuirSet.centralK[0] << std::endl;
  message << " sigmaKx=" << langmuirSet.sigmaK[0] << std::endl;
  message << " endTime=" << langmuirSet.endTime << std::endl;

  grid.printInfoFile("langmuirSetDataFile.txt",message.str());
}

void UTILITIES::setLangmuirWaveSet(LANGMUIRset &langmuirSet, GRID &grid){
  if(langmuirSet.checkLangmuirSetValidity){
  UTILITIES::allocateAccessibleKModes(langmuirSet.gridModes, grid);
  UTILITIES::newSetKModesToBeInitialised(langmuirSet);
  UTILITIES::exchangeKModesToBeInitialised(langmuirSet.myKModes, grid);

  }
  else{
    grid.printMessage("NO Langmuir Wave Set is possible");
  }
  UTILITIES::writeLangmuirSetDataFile(langmuirSet, grid);
  UTILITIES::writeGridModes(langmuirSet.gridModes, grid);
  UTILITIES::writeKModesToBeInitialised(langmuirSet.myKModes, grid);

}

void UTILITIES::OldMoveParticles(GRID* grid, SPECIE* specie, double amplitude,double lambda){
  int Npart = specie->Np;
  double kdx = 2*M_PI/lambda;
  double density=specie->plasma.params.density_coefficient;
  double deltaX = amplitude;
  double deltaV = amplitude*2*M_PI*sqrt(density);
  double oldX;

  if(false){
    std::cout<< "sposto le particelle che sono" << Npart << std::endl;
    std::cout<< "density = " << density << std::endl;
    std::cout<< " ===================== "<< std::endl;
  }
  for(int n=0;n<Npart;n++){
    oldX = specie->r0(n);
    specie->r0(n) += deltaX*cos(kdx*oldX);
    specie->u0(n) += deltaV*sin(kdx*oldX+2*M_PI*sqrt(density)*grid->dt*0.5);
  }
}

void UTILITIES::moveParticles(GRID* grid, SPECIE* specie, std::vector<KMODE> myKModes){
  int Npart=specie->Np;
  double density=specie->plasma.params.density_coefficient;
  double kL, phi, dr;
  double kx, dx, dVx, oldx;
  double ky, dy, dVy, oldy;
  double kz, dz, dVz, oldz;
  for(int n=0;n<Npart;n++){
    oldx = specie->r0(n);
    oldy = specie->r1(n);
    oldz = specie->r2(n);
    for(size_t m=0; m < myKModes.size(); m++){
      kx  = myKModes[m].k[0];
      ky  = myKModes[m].k[1];
      kz  = myKModes[m].k[2];
      kL  = sqrt(kx*kx + ky*ky + kz*kz);

      if(fabs(kL)>1e-2){
        dr  = myKModes[m].amplitude/kL;
        phi = myKModes[m].phase;

        dx  = dr*kx/kL;
        dy  = dr*ky/kL;
        dz  = dr*kz/kL;
        dVx  = dx*2*M_PI*sqrt(density);
        dVy  = dy*2*M_PI*sqrt(density);
        dVz  = dz*2*M_PI*sqrt(density);

        phi += kx*oldx + ky*oldy + kz*oldz;
        specie->r0(n) += dx*cos(phi);
        specie->r1(n) += dy*cos(phi);
        specie->r2(n) += dz*cos(phi);

        phi += 2*M_PI*sqrt(density)*grid->dt*0.5;
        //specie->u0(n) += dVx*sin(phi);
        //specie->u1(n) += dVy*sin(phi);
        //specie->u2(n) += dVz*sin(phi);
      }
    }
  }
}

void UTILITIES::setExternaField(EM_FIELD &exfield, GRID &mygrid, double time, LANGMUIRset &langmuirSet){
  exfield.setAllValuesToZero();

  langmuirSet.keepForcing =  (time <= langmuirSet.endTime);

  if(!langmuirSet.keepForcing){
    return;
  }
  double rampa=langmuirSet.growthRate;
  double amplitude =  1.0 - exp(-(time)/rampa);
  double Temperature=langmuirSet.refTemp;
  for (int k = 0; k < mygrid.Nloc[2]; k++){
    double zz = mygrid.cirloc[2][k];
    for (int j = 0; j < mygrid.Nloc[1]; j++){
      double yy = mygrid.cirloc[1][j];
      for (int i = 0; i < mygrid.Nloc[0]; i++){
        double xx = mygrid.cirloc[0][i];
        for(size_t m=0; m < langmuirSet.myKModes.size(); m++){
          double kL, phi, dE;
          double kx, dEx;
          double ky, dEy;
          double kz, dEz;
          double wL;
          kx  = langmuirSet.myKModes[m].k[0];
          ky  = langmuirSet.myKModes[m].k[1];
          kz  = langmuirSet.myKModes[m].k[2];
          kL  = sqrt(kx*kx + ky*ky + kz*kz);
          wL  = sqrt(4*M_PI*M_PI+3*kL*kL*Temperature);
          if(fabs(kL)>1e-2){
//            dE = amplitude*4*M_PI*langmuirSet.myKModes[m].amplitude/(kL);
            dE = amplitude*langmuirSet.myKModes[m].amplitude/(kL);

            phi = langmuirSet.myKModes[m].phase;
            phi += kx*xx + ky*yy + kz*zz;
            phi -= wL*time;

            dEx  = dE*kx/kL*cos(phi);
            dEy  = dE*ky/kL*cos(phi);
            dEz  = dE*kz/kL*cos(phi);

            exfield.E0(i,j,k) +=dEx;
            if(mygrid.getDimensionality()>=2)
              exfield.E1(i,j,k) +=dEy;
            if(mygrid.getDimensionality()>=3)
              exfield.E2(i,j,k) +=dEz;
          }
        }
      }
    }
  }
}

void UTILITIES::printTotalNumberOfParticles(std::vector<SPECIE*> species, GRID &mygrid){
  uint64_t totPartNum = 0;
  for (std::vector<SPECIE*>::const_iterator spec_iterator = species.begin(); spec_iterator != species.end(); spec_iterator++) {
    totPartNum += (*spec_iterator)->printParticleNumber();
  }

  std::stringstream messaggio;
  messaggio << "Total particle number: " << totPartNum << std::endl;
  mygrid.printMessage(messaggio.str());
}

void UTILITIES::launchPoissonSolver(EM_FIELD &myfield, std::vector<SPECIE*> species, GRID &grid, CURRENT &current){
  if(grid.isWithPoisson()){
    std::stringstream messaggio;
    messaggio << " evaluating density..." << std::endl;
    grid.printMessage(messaggio.str());

    bool withSign = true;
    current.setAllValuesToZero();
    for (std::vector<SPECIE*>::const_iterator spec_iterator = species.begin(); spec_iterator != species.end(); spec_iterator++) {
      (*spec_iterator)->density_deposition_standard(&current, withSign);
    }
    current.pbc();

    messaggio.str(std::string());
    messaggio << "   done... now into Poisson solver" << std::endl;
    grid.printMessage(messaggio.str());

    myfield.poissonSolver(&current);
  }
  current.setAllValuesToZero();
}
