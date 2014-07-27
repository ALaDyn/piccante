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


#include "utilities.h"

void moveWindow(GRID* _mygrid, EM_FIELD* _myfield, std::vector<SPECIE*> _myspecies){
    _mygrid->move_window();
    _myfield->move_window();
    for (std::vector<SPECIE*>::iterator spec_iterator = _myspecies.begin(); spec_iterator != _myspecies.end(); spec_iterator++){
        (*spec_iterator)->move_window();
    }
}

void restartFromDump(int *_dumpID, GRID* mygrid, EM_FIELD* myfield, std::vector<SPECIE*> species){
    int dumpID=_dumpID[0];
    std::ifstream dumpFile;
    dumpFile.open( mygrid->composeDumpFileName(dumpID).c_str() );
    if( dumpFile.good()){
        mygrid->reloadDump(dumpFile);
        myfield->reloadDump(dumpFile);
        for (std::vector<SPECIE*>::iterator spec_iterator = species.begin(); spec_iterator != species.end(); spec_iterator++){
            (*spec_iterator)->reloadDump(dumpFile);
        }
        dumpFile.close();
        dumpID++;
    }
    _dumpID[0]=dumpID;
}

void dumpFilesForRestart(int *_dumpID, GRID* mygrid, EM_FIELD* myfield, std::vector<SPECIE*> species){
    int dumpID=_dumpID[0];
    std::ofstream dumpFile;
    dumpFile.open( mygrid->composeDumpFileName(dumpID).c_str() );
    mygrid->dump(dumpFile);
    myfield->dump(dumpFile);
    for (std::vector<SPECIE*>::iterator spec_iterator = species.begin(); spec_iterator != species.end(); spec_iterator++){
        (*spec_iterator)->dump(dumpFile);
    }
    dumpFile.close();
    dumpID++;
    _dumpID[0]=dumpID;
    MPI_Barrier(MPI_COMM_WORLD);
    if(mygrid->myid==mygrid->master_proc){
        printf("\t DUMP #%i done!\n", (dumpID-1));
    }
}

