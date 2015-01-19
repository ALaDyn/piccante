COMPILER = mpic++
EXE = piccante

OBJ = main-1.o grid.o structures.o current.o em_field.o particle_species.o output_manager.o utilities.o jsoncpp.o jsonparser.o

OPT = -O3 

LFLAGS = -Wall

LIB = -lgsl -lgslcblas 

all : $(EXE)

boost : OPT = -O3 -DUSE_BOOST
boost : LIB = -lgsl -lgslcblas  -lboost_filesystem-mt -lboost_system-mt
boost : $(EXE)

hdf5 : OPT = -O3  -DUSE_BOOST -I/usr/lib/hdf5-1.8.12/hdf5/include -DUSE_HDF5
hdf5 : LIB = -lgsl -lgslcblas -lboost_filesystemt -lboost_system -lhdf5  -L/usr/lib/hdf5-1.8.12/hdf5/lib
hdf5 : $(EXE)

warn : OPT = -O3 -Wall -Winline -Wextra
warn : $(EXE)

debug : OPT = -O0 -g 
debug : $(EXE)

perf : OPT = -O3 -g -pg
perf : $(EXE)

scal : COMPILER = scalasca -instrument mpic++
scal : OPT = -O3 -g
scal : $(EXE) 

vec : OPT = -O3 -ftree-vectorize -msse2 -ftree-vectorizer-verbose=5
vec : $(EXE)

$(EXE) : $(OBJ)
				$(COMPILER) -o $(EXE) $(OPT) $(OBJ)  $(LIB) 

main-1.o: main-1.cpp 
				$(COMPILER) $(OPT) -c main-1.cpp

grid.o : grid.cpp 
				$(COMPILER) $(OPT) -c grid.cpp

structures.o : structures.cpp 
					$(COMPILER) $(OPT) -c structures.cpp 

current.o : current.cpp
					$(COMPILER) $(OPT) -c current.cpp

em_field.o: em_field.cpp
						$(COMPILER) $(OPT) -c em_field.cpp

particle_species.o: particle_species.cpp 
						$(COMPILER) $(OPT) -c particle_species.cpp 

output_manager.o:  output_manager.cpp 
		$(COMPILER) $(OPT) -c output_manager.cpp

utilities.o: utilities.cpp
		$(COMPILER) $(OPT) -c utilities.cpp
jsoncpp.o: jsoncpp.cpp
		$(COMPILER) $(OPT) -c jsoncpp.cpp
jsonparser.o: jsonparser.cpp
		$(COMPILER) $(OPT) -c jsonparser.cpp

clean :
				rm -f $(OBJ)


