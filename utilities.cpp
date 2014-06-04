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
}

