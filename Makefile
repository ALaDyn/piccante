COMPILER = mpicxx
EXE = piccante.exe
MAINFILE = main-piccante.cpp
MAINFILE_DEV = main-devel.cpp

FILES = grid.cpp \
        sobol.cpp \
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
OBJDEV = $(addprefix $(OBJ_FOLDER)/, $(addsuffix .o, $(basename $(MAINFILE_DEV))))
LIB =
RPATH = 
OPT = -O3 -std=c++11

BOOST_LIB = $(SRC_FOLDER)
BOOST_INC = $(SRC_FOLDER)
HDF5_INC = $(SRC_FOLDER)
HDF5_LIB = $(SRC_FOLDER)

ifneq (,$(findstring devel,$(config)))
	MAINFILE = $(MAINFILE_DEV)
endif

all: $(EXE)

boost: OPT += -DUSE_BOOST
boost: LIB += -lboost_filesystem -lboost_system
boost: all

nocpp11: boost
nocpp11: OPT = -O3 -DNO_CXX11 -lboost_random
nocpp11: all

cnaf: boost
cnaf: BOOST_LIB = /shared/software/BOOST/boost_1_56_0/lib
cnaf: BOOST_INC = /shared/software/BOOST/boost_1_56_0/include
cnaf: all

cnaf-intel: boost
cnaf-intel: COMPILER = mpiicpc
cnaf-intel: OPT += -axSSE4.2,AVX -ipo
cnaf-intel: OPT_REPORT += -vec-report -opt-report 3 
cnaf-intel: BOOST_LIB = /shared/software/BOOST/boost_1_56_0/lib
cnaf-intel: BOOST_INC = /shared/software/BOOST/boost_1_56_0/include
cnaf-intel: all

cnaf-phi: COMPILER = mpiicpc
cnaf-phi: OPT = -O3 -std=c++11 -mmic 
cnaf-phi: RPATH = -Wl,-rpath=/shared/software/compilers/intel/compilers_and_libraries_2016.0.109/linux/compiler/lib/intel64_lin_mic 
cnaf-phi: all

brew: 
	@echo ' '
	@echo '"brew" is deprecated: use "make boost" instead'
	@echo 'if it fails, check Makefile and adapt the rule "brewold"'
	@echo ' '

brewold: boost
brewold: BOOST_LIB = /usr/local/Cellar/boost/1.60.0_2/lib
brewold: BOOST_INC = /usr/local/Cellar/boost/1.60.0_2/include
brewold: all

hdf5: boost
hdf5: OPT += -DUSE_HDF5
hdf5: HDF5_INC = /usr/lib/hdf5-1.8.12/hdf5/include
hdf5: HDF5_LIB = /usr/lib/hdf5-1.8.12/hdf5/lib
hdf5: LIB += -lhdf5
hdf5: all

warn: OPT += -Wall -Winline -Wextra
warn: all

debug: OPT = -O0 -g -std=c++11 -Wall
debug: all

profiling: OPT += -g -pg
profiling: all

scalasca: COMPILER = scalasca -instrument mpicxx
scalasca: OPT += -g
scalasca: all 

sse2-vec: OPT += -ftree-vectorize -msse2 -ftree-vectorizer-verbose=5
sse2-vec: all

marconi: boost
marconi: OPT += -xCORE-AVX2 -ipo
marconi: BOOST_LIB = /cineca/prod/opt/libraries/boost/1.61.0/intelmpi--5.1--binary/lib
marconi: BOOST_INC = /cineca/prod/opt/libraries/boost/1.61.0/intelmpi--5.1--binary/include
marconi: all

galileo: boost
galileo: OPT += -xCORE-AVX2 -ipo
galileo: BOOST_INC = /cineca/prod/libraries/boost/1.57.0/intel--cs-xe-2015--binary/include
galileo: BOOST_LIB = /cineca/prod/libraries/boost/1.57.0/intel--cs-xe-2015--binary/lib
galileo: all

galileo-phi: COMPILER = mpiicpc
galileo-phi: OPT = -O3 -std=c++11 -mmic
galileo-phi: all

fermi: COMPILER = mpixlcxx
fermi: LIB = -lboost_filesystem -lboost_system -lboost_random
fermi: BOOST_LIB = /cineca/prod/libraries/boost/1.51.0/bgq-xl--1.0/lib/
fermi: BOOST_INC = /cineca/prod/libraries/boost/1.51.0/bgq-xl--1.0/include/
fermi: all

fermi-perf: fermi
fermi-perf: OPT = -DUSE_BOOST -DNO_CXX11 -qstrict -O5 -qipa=partition=large -qarch=qp -qtune=qp -qmaxmem=-1
#fermi-perf: OPT = -DUSE_BOOST -DNO_CXX11 -qipa=level=2 -qipa=partition=large -O5 -qstrict -qinline -qhot -qlibmpi -qarch=qp -qtune=qp -qmaxmem=-1
#fermi-perf: OPT = -DUSE_BOOST -DNO_CXX11 -O5 -qstrict -qinline -qhot -qlibmpi -qarch=qp -qtune=qp -qmaxmem=-1
fermi-perf: all

fermi-scalasca: fermi
fermi-scalasca: COMPILER = scalasca -instrument mpixlcxx
fermi-scalasca: OPT = -DUSE_BOOST -DNO_CXX11
fermi-scalasca: EXE = piccante.scalasca
fermi-scalasca: all

fermi-debug: fermi
fermi-debug: OPT = -DUSE_BOOST -DNO_CXX11 -g -qfullpath -qcheck -qflttrap -qinitauto=FF -qkeepparm
fermi-debug: EXE = piccante.debug
fermi-debug: all

fermi-debug-ipa: fermi
fermi-debug-ipa: OPT = -DUSE_BOOST -DNO_CXX11 -qipa=partition=large -qarch=qp -qtune=qp -qmaxmem=-1 -g -qfullpath -qcheck -qflttrap -qinitauto=FF
fermi-debug-ipa: EXE = piccante.debug
fermi-debug-ipa: all

juqueen: nocpp11
juqueen: COMPILER = mpixlcxx
juqueen: OPT = -DUSE_BOOST -DNO_CXX11 -qstrict -O5 -qipa=partition=large -qarch=qp -qtune=qp -qmaxmem=-1
juqueen: LIB = -lboost_system-1_47 -lboost_filesystem-1_47 -lboost_random-1_47
juqueen: BOOST_LIB = /bgsys/local/boost/1.47.0/lib
juqueen: BOOST_INC = /bgsys/local/boost/1.47.0
juqueen: all


$(EXE): $(OBJ)
	$(COMPILER) $(OPT)  -L$(BOOST_LIB) -L$(HDF5_LIB) $(RPATH) -o $(EXE) $(OBJ) $(LIB)

$(OBJ_FOLDER)/%.o: $(SRC_FOLDER)/%.cpp $(SRC_FOLDER)/%.h $(SRC_FOLDER)/preproc_defs.h
	$(COMPILER) $(OPT)  -I$(BOOST_INC) -I$(HDF5_INC) -c -o $@ $<

clean:
	rm -f $(OBJ) $(OBJDEV) *~

cleanall:
	rm -f $(OBJ) $(OBJDEV) $(EXE) $(EXE).debug *~

help: 
	@echo 'Usage: make config=OPTIONS'
	@echo '	    OPTIONS is a string composed of one or more of:'
	@echo '	        devel      : to compile using main-devel.cpp instead of default'
	@echo ' examples:'
	@echo '     make config=devel'

