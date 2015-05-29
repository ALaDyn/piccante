COMPILER = mpic++
EXE = piccante
MAIN = main-piccante.cpp

OBJ = main-piccante.o grid.o structures.o current.o em_field.o particle_species.o output_manager.o utilities.o jsoncpp.o jsonparser.o

OPT = -O3 
OPT_REPORT =

LFLAGS = -Wall

LIB = -lgsl -lgslcblas 

all : $(EXE)

specific : $(EXE)
specific : MAIN = main-1.cpp

intelcnaf : COMPILER = mpiicpc
intelcnaf : OPT += -DUSE_BOOST -axSSE4.2,AVX -ipo
intelcnaf : OPT_REPORT += -vec-report -opt-report 3 
intelcnaf : LIB += -lboost_filesystem -lboost_system
intelcnaf : $(EXE)

boost : OPT = -O3 -DUSE_BOOST
boost : LIB = -lgsl -lgslcblas  -lboost_filesystem -lboost_system
boost : $(EXE)

hdf5 : OPT = -O3  -DUSE_BOOST -I/usr/lib/hdf5-1.8.12/hdf5/include -DUSE_HDF5
hdf5 : LIB = -lgsl -lgslcblas -lboost_filesystem -lboost_system -lhdf5  -L/usr/lib/hdf5-1.8.12/hdf5/lib
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
				$(COMPILER) $(OPT) -o $(EXE) $(OBJ)  $(LIB) 

main-piccante.o: $(MAIN) 
				$(COMPILER) $(OPT) -c $(MAIN) $(OPT_REPORT)

grid.o : grid.cpp 
				$(COMPILER) $(OPT) -c grid.cpp $(OPT_REPORT)

structures.o : structures.cpp 
					$(COMPILER) $(OPT) -c structures.cpp $(OPT_REPORT) 

current.o : current.cpp
					$(COMPILER) $(OPT) -c current.cpp $(OPT_REPORT)

em_field.o: em_field.cpp
						$(COMPILER) $(OPT) -c em_field.cpp $(OPT_REPORT)

particle_species.o: particle_species.cpp 
						$(COMPILER) $(OPT) -c particle_species.cpp  $(OPT_REPORT)

output_manager.o:  output_manager.cpp 
		$(COMPILER) $(OPT) -c output_manager.cpp $(OPT_REPORT)

utilities.o: utilities.cpp
		$(COMPILER) $(OPT) -c utilities.cpp $(OPT_REPORT)

jsoncpp.o: jsoncpp.cpp
		$(COMPILER) $(OPT) -c jsoncpp.cpp 

jsonparser.o: jsonparser.cpp
		$(COMPILER) $(OPT) -c jsonparser.cpp 

clean :
				rm -f $(OBJ)


