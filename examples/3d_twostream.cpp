#include<stdio.h>
#include<iostream>
#include<fstream>
#include<sstream>
#include<malloc.h>
#include<math.h>
#include<iomanip>
#include<string.h>
#include<ctime>       /* time */
#include"mpi.h"
#include <gsl/gsl_rng.h> // gnu scientific linux per generatore di numeri casuali
#include <gsl/gsl_randist.h>
#include <cstdarg> //Per chiamare funzioni con numero variabile di argomenti
#include <vector>

using namespace std;

#define DIMENSIONALITY 3
#include "access.h"
#include "commons.h"
#include "grid.h"
#include "structures.h"
#include "current.h"
#include "em_field.h"
#include "particle_species.h"
//#include "diag_manager.h"
#include "output_manager.h"

#define NPROC_ALONG_Y 16
#define NPROC_ALONG_Z 8


#define DIRECTORY_OUTPUT "TEST"
#define RANDOM_NUMBER_GENERATOR_SEED 5489

double left_grating(double, double, double, PLASMAparams, int, int);

int main(int narg, char **args)
{
    GRID grid;
    EM_FIELD myfield;
    CURRENT current;
    std::vector<SPECIE*> species;
    vector<SPECIE*>::const_iterator spec_iterator;
    int istep;
    gsl_rng* rng = gsl_rng_alloc(gsl_rng_ranlxd1);

    //*******************************************INIZIO DEFINIZIONE GRIGLIA*******************************************************

    grid.setXrange(-2.0,2.0);
    grid.setYrange(-2.0,2.0);
    grid.setZrange(-0.5 ,+0.5);

    grid.setNCells(400,400,100);
    grid.setNProcsAlongY(NPROC_ALONG_Y);
    grid.setNProcsAlongZ(NPROC_ALONG_Z);

    //grid.enableStretchedGrid();
    grid.setXandNxLeftStretchedGrid(-20.0,250);
    grid.setXandNxRightStretchedGrid(20.0,250);
    grid.setYandNyLeftStretchedGrid(-20.0,250);
    grid.setYandNyRightStretchedGrid(20.0,250);

    grid.setBoundaries(xPBC|yPBC|zPBC);
    grid.mpi_grid_initialize(&narg, args);
    grid.setCourantFactor(0.98);

    grid.setSimulationTime(100.05);

    grid.with_particles = YES;//NO;
    grid.with_current   = YES;//YES;
    //double start, beta_mw;	int frequency_of_shifts;
    //grid.setMovingWindow(start=0, beta_mw=0.0, frequency_of_shifts=10);

    grid.setMasterProc(0);

    grid.finalize();

    srand (time(NULL));
    grid.initRNG(rng,RANDOM_NUMBER_GENERATOR_SEED);

    grid.visualDiag();
    if(grid.myid==grid.master_proc){
      std::ofstream of1;
        of1.open("merda_dx.txt");
        for(int i=0; i<grid.Nloc[0]; i++){
            of1<<i<<"\t"<<grid.dr[0]/grid.iStretchingDerivativeCorrection[0][i]<<"\n";
        }
        of1.close();

        of1.open("merda_dy.txt");
        for(int i=0; i<grid.Nloc[1]; i++){
            of1<<i<<"\t"<<grid.dr[1]/grid.iStretchingDerivativeCorrection[1][i]<<"\n";
        }
        of1.close();

        of1.open("merda_x.txt");
        for(int i=0; i<grid.NGridNodes[0]; i++){
            of1<<i<<"\t"<<grid.cir[0][i]<<"\n";
        }
        of1.close();

        of1.open("merda_y.txt");
        for(int i=0; i<grid.NGridNodes[1]; i++){
            of1<<i<<"\t"<<grid.cir[1][i]<<"\n";
        }
        of1.close();
    }
    //********************************************FINE DEFINIZIONE GRIGLIA********************************************************

    //*******************************************INIZIO DEFINIZIONE CAMPI*********************************************************
    myfield.allocate(&grid);
    myfield.setAllValuesToZero();    

    laserPulse pulse1;
    pulse1.type = GAUSSIAN;                        //Opzioni : GAUSSIAN, PLANE_WAVE, COS2_PLANE_WAVE
    pulse1.polarization = CIRCULAR_POLARIZATION;
    pulse1.t_FWHM = 5.0;
    pulse1.waist = 3.0;
    pulse1.focus_position = 0.0;
    pulse1.laser_pulse_initial_position = -15;
    pulse1.lambda0 = 1.0;
    pulse1.normalized_amplitude = 10.0;
    pulse1.rotation = false;
    pulse1.angle = 2.0*M_PI*(45.0/360.0);
    pulse1.rotation_center_along_x = 0.0;

    laserPulse pulse2;
    pulse2 = pulse1;
    pulse2.angle = 2.0*M_PI*(15.0/360.0);
    pulse2.laser_pulse_initial_position = -10.0;

    //myfield.addPulse(&pulse1);
    //myfield.addPulse(&pulse2);

    myfield.boundary_conditions();
    //myfield.smooth_filter(10);   

    current.allocate(&grid);
    current.setAllValuesToZero();
    //*******************************************FINE DEFINIZIONE CAMPI***********************************************************

    //*******************************************INIZIO DEFINIZIONE SPECIE*********************************************************
    PLASMA plasma1;
    plasma1.density_function = box;      //Opzioni: box, left_linear_ramp, left_soft_ramp, left_grating
    plasma1.setMinBox(-10.0 , -10.0, grid.rmin[2]);   //double (* distrib_function)(double x, double y, double z, PLASMAparams plist, int Z, int A)
    plasma1.setMaxBox(10.0, 10.0, grid.rmax[2]);    //PLASMAparams: rminbox[3], rmaxbox[3], ramp_length, density_coefficient,
    plasma1.setRampLength(0.2);                       //ramp_min_density,void *additional_params
    plasma1.setDensityCoefficient(1.0);         // Per grating double g_depth = paramlist[0];double g_lambda = paramlist[1];
    plasma1.setRampMinDensity(0.001);                 //double g_phase = paramlist[2];

  

    // double gratingParams[3] = {0.2, sqrt(2), 0.0};
    //plasma1.setAdditionalParams(gratingParams);

    PLASMA plasma2;
    plasma2.density_function = box;
    plasma2.setMinBox(0.3 , -60.0, grid.rmin[2]);
    plasma2.setMaxBox(1.0, 60.0, grid.rmax[2]);
    plasma2.setDensityCoefficient(100);

    PLASMA plasma3;
    plasma3.density_function = box;
    plasma3.setMinBox(1.0 , -60.0, grid.rmin[2]);
    plasma3.setMaxBox(4.0, 60.0, grid.rmax[2]);
    plasma3.setDensityCoefficient(100);

    PLASMA plasma4;
    plasma4.density_function = left_linear_ramp;      //Opzioni: box, left_linear_ramp, left_soft_ramp, left_grating                                                                                        
    plasma4.setMinBox(0.0 , -60.0, grid.rmin[2]);   //double (* distrib_function)(double x, double y, double z, PLASMAparams plist, int Z, int A)                                                           
    plasma4.setMaxBox(0.2, 60.0, grid.rmax[2]);    //PLASMAparams: rminbox[3], rmaxbox[3], ramp_length, density_coefficient,                                                                                
    plasma4.setRampLength(0.2);                       //ramp_min_density,void *additional_params                                                                                                            
    plasma4.setDensityCoefficient(1.0);         // Per grating double g_depth = paramlist[0];double g_lambda = paramlist[1];                                                                              
    plasma4.setRampMinDensity(0.0);

    SPECIE  electrons1(&grid);
    electrons1.plasma = plasma1;
    electrons1.setParticlesPerCellXYZ(3, 3, 3);       //Se < 1 il nPPC viene sostituito con 1
    electrons1.setName("ELE1");
    electrons1.type=ELECTRON;
    electrons1.creation();                            //electrons.isTestSpecies=true disabilita deposizione corrente.
    species.push_back(&electrons1);


    SPECIE electrons2(&grid);
    electrons2.plasma=plasma1;
    electrons2.setParticlesPerCellXYZ(3, 3, 3);
    electrons2.setName("ELE2");
    electrons2.type=ELECTRON;
    electrons2.creation();                            //electrons.isTestSpecies=true disabilita deposizione corrente.
    species.push_back(&electrons2);


    SPECIE electrons3(&grid);
    electrons3.plasma=plasma3;
    electrons3.setParticlesPerCellXYZ(1, 1, 1);
    electrons3.setName("ELE3");
    electrons3.type=ELECTRON;
    //electrons3.creation();                            //electrons.isTestSpecies=true disabilita deposizione corrente.                                                                                        
    //species.push_back(&electrons3);

    SPECIE electrons4(&grid);
    electrons4.plasma=plasma4;
    electrons4.setParticlesPerCellXYZ(5, 5, 1);
    electrons4.setName("ELE4");
    electrons4.type=ELECTRON;
    //electrons4.creation();                            //electrons.isTestSpecies=true disabilita deposizione corrente.                                                                                                                                                                                                                                                                                                   
    //species.push_back(&electrons4);




    SPECIE ions1(&grid);
    ions1.plasma = plasma1;
    ions1.setParticlesPerCellXYZ(10, 10, 1);
    ions1.setName("ION1");
    ions1.type=ION;
    ions1.Z= 6.0;
    ions1.A =12.0;
    //ions1.creation();
    //species.push_back(&ions1);


    SPECIE  ions2(&grid);
    ions2.plasma = plasma2;
    ions2.setParticlesPerCellXYZ(2, 2, 1);       //Se < 1 il nPPC viene sostituito con 1
    ions2.setName("ION2");
    ions2.type=ION;
    ions2.Z= 6.0;
    ions2.A =12.0;
    //    ions2.creation();                            //electrons.isTestSpecies=true disabilita deposizione corrente.
    //species.push_back(&ions2);


    SPECIE ions3(&grid);
    ions3.plasma = plasma3;
    ions3.setParticlesPerCellXYZ(1, 1, 1);
    ions3.setName("ION3");
    ions3.type=ION;
    ions3.Z= 6.0;
    ions3.A =12.0;
    //ions3.creation();
    //species.push_back(&ions3);

    SPECIE ions4(&grid);
    ions4.plasma = plasma4;
    ions4.setParticlesPerCellXYZ(4, 4, 1);
    ions4.setName("ION4");
    ions4.type=ION;
    ions4.Z= 1.0;
    ions4.A =1.0;
    //ions4.creation();
    //species.push_back(&ions4);

    // WATERBAG        : [P0] -P0< px < + P0, -P0< py < + P0, -P0< pz < + P0 uniforme
    // WATERBAG_3TEMP  : [P0_X,P0_Y,P0_Z] -P0_X< px < + P0_X, -P0_Y< py < + P0_Y, -P0_Z< pz < + P0_Z uniforme
    // UNIF_SPHERE     : [P0] -P0 < p < +P0 uniforme
    // SUPERGAUSSIAN   : [P0, ALPHA] f(p) = C*exp(-abs(p/P0)^(ALPHA))
    // MAXWELL         : [Ta] Maxwell alla Macchi
    // JUTTNER         : [a] f(p) = C*exp(-a*gamma(p)) [DA RISCRIVERE]
    //si può ancora usare electrons.add_momenta(1.0, 0.0, 0);

    tempDistrib distribution;
    distribution.setMaxwell(1.0e-5);

    electrons1.add_momenta(rng,0.0, 0.0, -1.0, distribution);
    electrons2.add_momenta(rng,0.0, 0.0, 1.0, distribution);
    //     ions1.add_momenta(rng,0.0, 0.0, 0.0, distribution);
//    ions2.add_momenta(rng,0.0, 0.0, 0.0, distribution);
//    //*******************************************FINE DEFINIZIONE CAMPI***********************************************************

    //*******************************************INIZIO DEFINIZIONE DIAGNOSTICHE**************************************************

   OUTPUT_MANAGER manager(&grid, &myfield, &current, species);

    manager.addEMFieldBinaryFrom(5.0,5.0);

    manager.addSpecDensityBinaryFrom("ELE1", 5.0,5.0);
    manager.addSpecDensityBinaryFrom("ELE2", 5.0,5.0);
    //manager.addSpecDensityBinaryFrom("ELE3", 0.0,5.0);
    //manager.addSpecDensityBinaryFrom("ELE4", 0.0,5.0);
 
    //manager.addSpecDensityBinaryFrom("ION1", 0.0,5.0);
    //manager.addSpecDensityBinaryFrom("ION2", 0.0,5.0);
    //manager.addSpecDensityBinaryFrom("ION3", 0.0,5.0);
    //manager.addSpecDensityBinaryFrom("ION4", 0.0,5.0);

    manager.addCurrentBinaryFrom(5.0,5.0);


    manager.addSpecPhaseSpaceBinaryFrom("ELE1", 5.0,5.0);
    manager.addSpecPhaseSpaceBinaryFrom("ELE2", 5.0,5.0);
    //manager.addSpecPhaseSpaceBinaryFrom("ION1", 0.0,5.0);  
    //manager.addSpecPhaseSpaceBinaryFrom("ION4", 0.0,5.0);

    manager.addDiagFrom(0.0,2.0);
	
    manager.initialize(DIRECTORY_OUTPUT);

    manager.autoVisualDiag();



    //*******************************************FINE DEFINIZIONE DIAGNOSTICHE**************************************************

    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ CICLO PRINCIPALE (NON MODIFICARE) @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    if(grid.myid==grid.master_proc){
        printf("----- START temporal cicle -----\n");
        fflush(stdout);
    }

    int Nstep=grid.getTotalNumberOfTimesteps();
    for(istep=0;istep<=Nstep;istep++)
    {

        grid.istep = istep;

        grid.printTStepEvery(1);


        manager.callDiags(istep);  /// deve tornare all'inizo del ciclo


        myfield.openBoundariesE_1();
        myfield.new_halfadvance_B();
        myfield.boundary_conditions();

	current.setAllValuesToZero();

        for(spec_iterator=species.begin(); spec_iterator!=species.end(); spec_iterator++){
            (*spec_iterator)->current_deposition_standard(&current);
        }

        current.pbc();


        for(spec_iterator=species.begin(); spec_iterator!=species.end(); spec_iterator++){
            (*spec_iterator)->position_parallel_pbc();
        }

        myfield.openBoundariesB();
        myfield.new_advance_E(&current);

        myfield.boundary_conditions();

        myfield.openBoundariesE_2();
        myfield.new_halfadvance_B();

        myfield.boundary_conditions();

        for(spec_iterator=species.begin(); spec_iterator!=species.end(); spec_iterator++){
            (*spec_iterator)->momenta_advance(&myfield);
        }


        grid.time+=grid.dt;

	myfield.move_window();
	for(spec_iterator=species.begin(); spec_iterator!=species.end(); spec_iterator++){
	  (*spec_iterator)->move_window();
	}


    }

    manager.close();
    MPI_Finalize();

}
