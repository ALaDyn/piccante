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

void UTILITIES::restartFromDump(int *_dumpID, GRID* mygrid, EM_FIELD* myfield, std::vector<SPECIE*> species) {
  int dumpID = _dumpID[0];
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
  float x = spheres.coords[index * 4];
  float y = spheres.coords[index * 4 + 1];
  float z = spheres.coords[index * 4 + 2];
  float r = spheres.coords[index * 4 + 3];

  double xl = grid.rminloc[0] - r;
  double yl = grid.rminloc[1] - r;
  double zl = grid.rminloc[2] - r;

  double xr = grid.rmaxloc[0] + r;
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
