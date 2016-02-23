COMPILER = mpicxx
EXE = piccante.exe
MAINFILE = main-piccante.cpp

FILES = grid.cpp \
        structures.cpp \
        current.cpp \
        em_field.cpp \
        particle_species.cpp \
        output_manager.cpp \
        utilities.cpp \
        jsoncpp.cpp \
        jsonparser.cpp

SRC_FOLDER = src
OBJ_FOLDER = obj
SRC = $(addprefix $(SRC_FOLDER)/, $(FILES))
MAIN = $(addprefix $(SRC_FOLDER)/, $(MAINFILE))
OBJ = $(addprefix $(OBJ_FOLDER)/, $(addsuffix .o, $(basename $(FILES))))
OBJ += $(addprefix $(OBJ_FOLDER)/, $(addsuffix .o, $(basename $(MAINFILE))))
#LIB = -lgsl -lgslcblas
LIB =
OPT = -O3 -std=c++11

#GSL_LIB = $(SRC_FOLDER)
#GSL_INC = $(SRC_FOLDER)
BOOST_LIB = $(SRC_FOLDER)
BOOST_INC = $(SRC_FOLDER)
HDF5_INC = $(SRC_FOLDER)
HDF5_LIB = $(SRC_FOLDER)

all: $(EXE)

boost: OPT += -DUSE_BOOST
boost: LIB += -lboost_filesystem -lboost_system
boost: all

cnaf: all

cnaf-intel: boost
cnaf-intel: COMPILER = mpiicpc
cnaf-intel: OPT += -axSSE4.2,AVX -ipo
cnaf-intel: OPT_REPORT += -vec-report -opt-report 3 
cnaf-intel: all

brew: boost
brew: BOOST_LIB = /usr/local/Cellar/boost/1.60.0_1/lib
brew: BOOST_INC = /usr/local/Cellar/boost/1.60.0_1/include
brew: all

hdf5: boost
hdf5: OPT += -DUSE_HDF5
hdf5: HDF5_INC = /usr/lib/hdf5-1.8.12/hdf5/include
hdf5: HDF5_LIB = /usr/lib/hdf5-1.8.12/hdf5/lib
hdf5: LIB += -lhdf5
hdf5: all

warn: OPT += -Wall -Winline -Wextra
warn: all

debug: OPT = -O0 -g 
debug: all

profiling: OPT += -g -pg
profiling: all

scalasca: COMPILER = scalasca -instrument mpicxx
scalasca: OPT += -g
scalasca: all 

sse2-vec: OPT += -ftree-vectorize -msse2 -ftree-vectorizer-verbose=5
sse2-vec: all

galileo: boost
galileo: OPT += -xCORE-AVX2 -ipo
#galileo: GSL_INC = /cineca/prod/libraries/gsl/1.16/intel--cs-xe-2015--binary/include
#galileo: GSL_LIB = /cineca/prod/libraries/gsl/1.16/intel--cs-xe-2015--binary/lib
galileo: BOOST_INC = /cineca/prod/libraries/boost/1.57.0/intel--cs-xe-2015--binary/include
galileo: BOOST_LIB = /cineca/prod/libraries/boost/1.57.0/intel--cs-xe-2015--binary/lib
galileo: all

fermi: boost
fermi: COMPILER = mpixlcxx
fermi: OPT = -qlanglvl=extended0x -qstrict -O5 -qipa=partition=large -qarch=qp -qtune=qp -qmaxmem=-1
#fermi: OPT = -qipa=level=2 -qipa=partition=large -O5 -qstrict -qinline -qhot -qlibmpi -qarch=qp -qtune=qp -qmaxmem=-1
#fermi: OPT = -O5 -qstrict -qinline -qhot -qlibmpi -qarch=qp -qtune=qp -qmaxmem=-1
#fermi: GSL_LIB = /cineca/prod/libraries/gsl/1.15/bgq-xl--1.0/lib/
#fermi: GSL_INC = /cineca/prod/libraries/gsl/1.15/bgq-xl--1.0/include/
fermi: BOOST_LIB = /cineca/prod/libraries/boost/1.51.0/bgq-xl--1.0/lib/
fermi: BOOST_INC = /cineca/prod/libraries/boost/1.51.0/bgq-xl--1.0/include/
fermi: all

juqueen: fermi
juqueen: LIB = -lgsl -lgslcblas -lboost_system-1_47 -lboost_filesystem-1_47
#juqueen: GSL_LIB = /bgsys/local/gsl/1.15_O3g/lib
#juqueen: GSL_INC = /bgsys/local/gsl/1.15_O3g/include
juqueen: BOOST_LIB = /bgsys/local/boost/1.47.0/lib
juqueen: BOOST_INC = /bgsys/local/boost/1.47.0
juqueen: all


fermi-scalasca: fermi
fermi-scalasca: COMPILER = scalasca -instrument mpixlcxx
fermi-scalasca: EXE = piccante.scalasca
fermi-scalasca: all

fermi-debug: fermi
fermi-debug: OPT = -g -qfullpath -qcheck -qflttrap -qinitauto=FF -qkeepparm
fermi-debug: EXE = piccante.debug
fermi-debug: all

fermi-debug-ipa: fermi
fermi-debug-ipa: OPT = -qipa=partition=large -qarch=qp -qtune=qp -qmaxmem=-1 -g -qfullpath -qcheck -qflttrap -qinitauto=FF
fermi-debug-ipa: EXE = piccante.debug
fermi-debug-ipa: all


$(EXE): $(OBJ)
	$(COMPILER) $(OPT)  -L$(BOOST_LIB) -L$(HDF5_LIB) -o $(EXE) $(OBJ) $(LIB)

$(OBJ_FOLDER)/%.o: $(SRC_FOLDER)/%.cpp
	$(COMPILER) $(OPT)  -I$(BOOST_INC) -I$(HDF5_INC) -c -o $@ $<

clean:
	rm -f $(OBJ) *~

cleanall:
	rm -f $(OBJ) $(EXE) $(EXE).debug *~

