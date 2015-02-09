/* Copyright 2014 - Andrea Sgattoni, Luca Fedeli, Stefano Sinigardi */

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

#ifndef JSONNAMES_H
#define JSONNAMES_H

#define _INPUT_LINE_INPUT_FILE_FLAG "-inputfile"
#define _DEFAULT_INPUT_FILE_NAME "inputPiccante.json"
#define _JSON_INT_VERSION "VERSION"
#define INT_DIMENSIONS "dimensions"
#define _JSON_BOOL_ENABLED "enabled"
#define _JSON_INT_MASTERPROC "masterProc"
#define _JSON_DOUBLE_COURANT_FACTOR "courantFactor"

#define _JSON_DOUBLE_FRICTION_LENGTH_CM "lambda0"
#define _JSON_BOOL_IS_FRICTION "radiationFriction"
#define DOUBLEARRAY_X_RANGE "xRange"
#define DOUBLEARRAY_Y_RANGE "yRange"
#define DOUBLEARRAY_Z_RANGE "zRange"
#define _INT_N_CELLS_ "NCells"
#define _INT_N_PROC_Y_ "nProcY"
#define _INT_N_PROC_Z_ "nProcZ"
#define _DOUBLE_SIMULATION_TIME_ "simulationTime"

#define _OBJ_BOUNDARIES_ "boundaries"
#define _TAG_PERIODIC_BC_ "periodic"
#define _TAG_OPEN_BC_ "open"
#define _TAG_PML_ "pml"

#define _OBJ_RESTART_ "restart"
#define _JSON_STRING_DUMP_FOLDER_NAME "DumpPath"
#define _BOOL_RESTART_ "doRestart"
#define _BOOL_DODUMP_ "doDump"
#define _DOUBLE_DUMPEVERY_ "dumpEvery"
#define _INT_RESTART_FROM_DUMP_  "restartFromDump"

#define _OBJ_STRETCHED_GRID_ "StretchedGrid"
#define _OBJ_X_STTETCHING_ "x"
#define _OBJ_Y_STTETCHING_ "y"
#define _OBJ_Z_STTETCHING_ "z"
#define _OBJ_LEFT_STRETCHING_ "left"
#define _OBJ_RIGHT_STRETCHING_ "right"
#define _INT_NCELL_STRETCHING_ "NCells"
#define _DOUBLE_LIMIT_STRETCHING "limit"

#define _OBJ_MOVING_WINDOW_ "MovingWindow"
#define _DOUBLE_START_MW_ "start"
#define _DOUBLE_BETA_MW_ "beta"
#define _DOUBLE_FREQUENCY_MW_ "frequency"

#define _OBJ_ARRAY_LASER_ "Laser"

#define _STRING_TYPE_ "type"
#define _LASERTYPEVALUE_COS_PLANE_WAVE_ "COS2_PLANE_WAVE"
#define _LASERTYPEVALUE_GAUSSIAN_ "GAUSSIAN"
#define _LASERTYPEVALUE_PLANE_WAVE_ "PLANE_WAVE"
#define _LASERTYPEVALUE_COS2_PLATEAU_PLANE_WAVE_ "COS2_PLATEAU_PLANE_WAVE"

#define _STRING_POLARIZATION_ "polarization"
#define _LASERPOLARIZATIONVALUE_P_ "P"
#define _LASERPOLARIZATIONVALUE_S_ "S"
#define _LASERPOLARIZATIONVALUE_C_ "C"

#define _DOUBLE_LASER_DURATION_FWHM_ "durationFWHM"
#define _DOUBLE_LASER_INITIAL_POSITION_ "initialPosition"
#define _DOUBLE_LASER_A_ "a"
#define _DOUBLE_LASER_WAIST_ "waist"
#define _DOUBLE_LASER_FOCUS_POSITION_ "focusPosition"
#define _DOUBLE_LASER_LAMBDA_ "lambda"
#define _BOOL_LASER_ROTATION_ "rotation"
#define _DOUBLE_ROTATION_ANGLE_ "angle"
#define _DOUBLE_ROTATION_CENTER_ "center"
#define _DOUBLE_LASER_RISE_TIME_ "riseTime"

#define _JSON_OBJARRAY_PLASMA_ "Plasma"
#define _JSON_STRING_PLASMA_NAME "name"
#define _JSON_STRING_PLASMA_DENSITYFUNCTION "densityFunction"
#define _JSON_DOUBLEARRAY_PLASMA_X_RANGEBOX "XRangeBox"
#define _JSON_DOUBLEARRAY_PLASMA_Y_RANGEBOX "YRangeBox"
#define _JSON_DOUBLEARRAY_PLASMA_Z_RANGEBOX "ZRangeBox"
#define _JSON_DOUBLE_DENSITY_LAMBDA "DensityLambda"
#define _JSON_DOUBLE_PLASMA_DENSITY "DensityCoefficient"
#define _JSON_DOUBLE_PLASMA_RAMP_LENGTH "RampLength"
#define _JSON_DOUBLE_PLASMA_RAMP_MIN_DENSITY "RampMinDensity"
#define _JSON_DOUBLE_PLASMA_LEFT_RAMP_LENGTH "LeftRampLength"
#define _JSON_DOUBLE_PLASMA_RIGHT_RAMP_LENGTH "RightRampLength"
#define _JSON_DOUBLE_PLASMA_SCALE_LENGTH "ScaleLength"
#define _JSON_DOUBLE_PLASMA_LEFT_SCALE_LENGTH "LeftScaleLength"
#define _JSON_DOUBLE_PLASMA_RIGHT_SCALE_LENGTH "RightScaleLength"
#define _JSON_DOUBLE_PLASMA_LEFT_RAMP_MIN_DENSITY "LeftRampMinDensity"
#define _JSON_DOUBLE_PLASMA_RIGHT_RAMP_MIN_DENSITY "RightRampMinDensity"
#define _JSON_DOUBLE_PLASMA_GRATING_DEPTH "GratingDepth"
#define _JSON_DOUBLE_PLASMA_GRATING_PERIOD "GratingPeriod"
#define _JSON_DOUBLE_PLASMA_GRATING_PHASE "GratingPhase"

#define _JSON_OBJARRAY_SPECIES "Species"
#define _JSON_STRING_SPECIE_NAME "name"
#define _JSON_STRING_SPECIE_PLASMANAME "plasma"
#define _JSON_INTARRAY3_PARTICLES_PER_CELL "ParticlesPerCell"
#define _JSON_DOUBLE_IONS_A "A"
#define _JSON_DOUBLE_IONS_Z "Z"
#define _JSON_BOOL_IS_MARKER "isMarker"
#define _JSON_BOOL_IS_TEST "isTest"
#define _JSON_STRING_DISTRIBUTION "distribution"

#define _JSON_STRING_SPEIES_TYPE "type"
#define SPECIES_TYPEVALUE_ELECTRON "ELECTRON"
#define SPECIES_TYPEVALUE_POSITRON "POSITRON"
#define SPECIES_TYPEVALUE_ION "ION"

#define DISTRIBUTIONVALUE_MAXWELL "Maxwell"
#define DISTRIBUTIONVALUE_JUTTNER "Juttner"
#define DISTRIBUTIONVALUE_WATERBAG "Waterbag"
#define DISTRIBUTIONVALUE_WATERBAG3TEMP "Waterbag3Temp"
#define DISTRIBUTIONVALUE_UNFORM_SPHERE "UnifSphere"
#define DISTRIBUTIONVALUE_SUPERGAUSSIAN "Supergaussian"
#define DISTRIBUTIONVALUE_SPECIAL "Special"
#define _JSON_DOUBLEARRAY_ADD_MOMENTUM "distributionAddMomentum"
#define _JSON_DOUBLEARRAY_DISTRIBUTION_PARAMS "distributionParams"

#define _JSON_OBJECTARRAY_DOMAINS "Domains"
#define _JSON_STRING_DOMAIN_NAME "name"
#define _JSON_INTARRAY3_FREE_DIM "freeDim"
#define _JSON_DOUBLEARRAY3_POINT_COORDINATE "pointCoord"
#define _JSON_INTARRAY2_DOMAIN_X_RANGE "xRange"
#define _JSON_INTARRAY2_DOMAIN_Y_RANGE "yRange"
#define _JSON_INTARRAY2_DOMAIN_Z_RANGE "zRange"

#define _JSON_STRING_OUTPUT_FOLDER_NAME "OutputPath"
#define _JSON_OBJARRAY_OUTPUT "Output"
#define _JSON_STRING_OUTPUT_TYPE "type"
#define OUTPUTTYPEVALUE_E "E"
#define OUTPUTTYPEVALUE_EB "EB"
#define OUTPUTTYPEVALUE_B "B"
#define OUTPUTTYPEVALUE_DENSITY "Density"
#define OUTPUTTYPEVALUE_CURRENT "Current"
#define OUTPUTTYPEVALUE_PHASESPACE "PhSp"
#define OUTPUTTYPEVALUE_DIAG "Diag"




#define _JSON_DOUBLE_OUTPUT_FROM "from"
#define _JSON_DOUBLE_OUTPUT_TO "to"
#define _JSON_DOUBLE_OUTPUT_AT "at"
#define _JSON_DOUBLE_OUTPUT_EVERY "every"
#define _JSON_STRING_IN_OUTPUT_DOMAIN_NAME "in"
#define _JSON_STRING_OUTPUT_SPECIES_NAME "spec"





#endif // JSONNAMES_H
