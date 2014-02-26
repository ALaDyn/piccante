COMPILER = mpic++
EXE = piccante

OBJ = main-1.o grid.o structures.o current.o em_field.o particle_species.o output_manager.o 

OPT = -O3  

LFLAGS = -Wall

LIB = -lgsl -lgslcblas -lboost_filesystem -lboost_system

all : $(EXE)

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

clean :
				rm $(OBJ)


