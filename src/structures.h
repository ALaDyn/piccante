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

#ifndef __STRUCTURES_H__
#define __STRUCTURES_H__

#define _USE_MATH_DEFINES
#define NUMBER_OF_PLASMA_FUNCTIONS 26

#include <math.h>
#include "commons.h"
#include <cstdio>
#include <iostream>

struct SPHERES {
  float *coords;
  int NSpheres;
  float rmin[3], rmax[3];
  float fillingFactor;
};

struct PLASMAparams {
  double rminbox[3];
  double rmaxbox[3];
  double left_ramp_length;
  double right_ramp_length;
  double left_scale_length;
  double right_scale_length;
  double density_coefficient;
  double left_ramp_min_density;
  double right_ramp_min_density;
  SPHERES *spheres;
  void *additional_params;
  PLASMAparams operator=(const PLASMAparams& p1);

};
#define NBIN_SPECTRUM 1000;
struct SPECIEspectrum {
  double Kmax;
  double Dk;
  int Nbin;
  double *values;
};

struct DUMP_CONTROL {
  bool doRestart, doDump;
  int restartFromDump;
  double dumpEvery;
};

typedef double(*distrib_function)(double x, double y, double z, PLASMAparams plist, double Z, double A);

class PLASMA
{

public:
  PLASMAparams params;
  distrib_function density_function;

  PLASMA();
  PLASMA(const PLASMA& other);
  PLASMA operator=(const PLASMA& p1);
  void setRampLength(double rlength);
  void setLeftRampLength(double rlength);
  void setRightRampLength(double rlength);
  void setScaleLength(double slength);
  void setLeftScaleLength(double slength);
  void setRightScaleLength(double slength);
  void setDensityCoefficient(double dcoeff);
  void setDensityCoefficient(double dcoeff, double lambda);
  void setRampMinDensity(double minden);
  void setLeftRampMinDensity(double minden);
  void setRightRampMinDensity(double minden);
  void setAdditionalParams(void* addpar);
  void setMinBox(double xmin, double ymin, double zmin);
  void setMaxBox(double xmax, double ymax, double zmax);
  void setXRangeBox(double xmin, double xmax);
  void setYRangeBox(double ymin, double ymax);
  void setZRangeBox(double zmin, double zmax);

  ~PLASMA();

  static const int maxdF = NUMBER_OF_PLASMA_FUNCTIONS;
  static const std::string dFNames[];
  static const distrib_function dFPoint[];

  static bool isGrating(int dfIndex);
  static bool isPillar2D(int dfIndex);
  static bool isNanotubes2D(int dfIndex);
  static bool isFoils2D(int dfIndex);
  static bool isUser1(int dfIndex);
  static bool isUser2(int dfIndex);

};


//Pre-defined density functions
double box(double x, double y, double z, PLASMAparams plist, double Z, double A);

double left_linear_ramp(double x, double y, double z, PLASMAparams plist, double Z, double A);
double right_linear_ramp(double x, double y, double z, PLASMAparams plist, double Z, double A);
double left_right_linear_ramp(double x, double y, double z, PLASMAparams plist, double Z, double A);
double left_fixed_exp_ramp(double x, double y, double z, PLASMAparams plist, double Z, double A);
double right_fixed_exp_ramp(double x, double y, double z, PLASMAparams plist, double Z, double A);
double left_right_fixed_exp_ramp(double x, double y, double z, PLASMAparams plist, double Z, double A);
double left_free_exp_ramp(double x, double y, double z, PLASMAparams plist, double Z, double A);
double right_free_exp_ramp(double x, double y, double z, PLASMAparams plist, double Z, double A);
double left_right_free_exp_ramp(double x, double y, double z, PLASMAparams plist, double Z, double A);

double left_soft_ramp(double x, double y, double z, PLASMAparams plist, double Z, double A);

double* rough_box_prepareAdditionalParams(my_rng_generator* rng, double roughness, double shift);
double rough_box(double x, double y, double z, PLASMAparams plist, double Z, double A);
double rough_box_edgeCalc(double x0, double y0, double x, double y, int order, double* FFT_params, double roughness);

double box_minus_box(double x, double y, double z, PLASMAparams plist, double Z, double A);

double left_grating(double x, double y, double z, PLASMAparams plist, double Z, double A);
double left_blazed_grating(double x, double y, double z, PLASMAparams plist, double Z, double A);

double square_func(double x);
double left_square_grating(double x, double y, double z, PLASMAparams plist, double Z, double A);

double guide(double x, double y, double z, PLASMAparams plist, double Z, double A);

double modGrat(double x, double y, double z, PLASMAparams plist, double Z, double A);

double spoofGrat(double x, double y, double z, PLASMAparams plist, double Z, double A);
double spheres(double x, double y, double z, PLASMAparams plist, double Z, double A);

double pillars2D (double x, double y, double z, PLASMAparams plist, double Z, double A);
double nanotubes2D (double x, double y, double z, PLASMAparams plist, double Z, double A);
double foils2D (double x, double y, double z, PLASMAparams plist, double Z, double A);

double demo_2D_resonator (double x, double y, double z, PLASMAparams plist, double Z, double A);
double user1(double x, double y, double z, PLASMAparams plist, double Z, double A);
double user2(double x, double y, double z, PLASMAparams plist, double Z, double A);
//************** LASER PULSE TYPES *******
enum laserPulseType { DEFAULT_PULSE, GAUSSIAN, PLANE_WAVE, COS2_PLANE_WAVE, COS2_PLATEAU_PLANE_WAVE, LAGUERRE_GAUSSIAN};
enum pulsePolarization { P_POLARIZATION, S_POLARIZATION, CIRCULAR_POLARIZATION };

class laserPulse
{
public:
  laserPulseType type;
  pulsePolarization polarization;
  double t_FWHM;
  double waist;
  double focus_position;
  double laser_pulse_initial_position;
  double normalized_amplitude;
  double lambda0;
  bool rotation;
  double angle;
  double rotation_center_along_x;
  double rise_time;
  int LG_m;
  int LG_l;


  laserPulse();
  ~laserPulse();
  laserPulse(const laserPulse& other);
  laserPulse operator=(const laserPulse& p1);
  void setFocusPosition(double _focus_position);
  void setPulseInitialPosition(double _laser_pulse_initial_position);
  void setLambda(double _lambda0);
  void setWaist(double _waist);
  void setDurationFWHM(double _t_FWHM);
  void setNormalizedAmplitude(double _normalized_amplitude);
  void setRiseTime(double _rise_time);
  void setRotationAngleAndCenter(double _angle, double _rotation_center_along_x);
  void setGaussianPulse(double _waist, double _t_FWHM, double _normalized_amplitude);
  void setPlaneWave(double _normalized_amplitude);
  void setCos2PlaneWave(double _t_FWHM, double _normalized_amplitude);
  void setCos2PlateauPlaneWave(double _t_FWHM, double _rise_time, double _normalized_amplitude);
  void setGaussianPulse();
  void setPlaneWave();
  void setCos2PlaneWave();
  void setCos2PlateauPlaneWave();
  void setLaguerreGaussian_m(int m);
  void setLaguerreGaussian_l(int l);

  void setPPolarization();
  void setSPolarization();
  void setCircularPolarization();

};

//************** FIELD COORDINATES *******

#define _INTG_CRD 0
#define _HALF_CRD 1;
#define _NULL_CRD -1;

struct integer_or_halfinteger {
  char x;
  char y;
  char z;
};

//************** PARTICLES DISTRIBUTION FUNCTION *******
enum tempDistribType { WATERBAG, WATERBAG_3TEMP, UNIF_SPHERE, SUPERGAUSSIAN, MAXWELL, JUTTNER, SPECIAL };

class tempDistrib {
public:

  tempDistribType type;

  tempDistrib();
  void setWaterbag(double _p0);
  void setWaterbag3Temp(double _p0_x, double _p0_y, double _p0_z);
  void setUnifSphere(double _p0);
  void setSupergaussian(double _p0, double _alpha);
  void setMaxwell(double _temp);
  void setJuttner(double _a);
  void setSpecial(double _a);

  bool isInit();

  friend class SPECIE;
private:
  double p0;
  double p0_x;
  double p0_y;
  double p0_z;
  double alpha;
  double temp;
  double a;
  bool init;

};
//**************** SPHERES *****************/

#endif

