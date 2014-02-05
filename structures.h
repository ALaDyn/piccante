#ifndef __STRUCTURES_H__
#define __STRUCTURES_H__

#define _USE_MATH_DEFINES

#include <cmath>
#include "commons.h"
#if defined(_MSC_VER)
#include "gsl/gsl_rng.h" // gnu scientific linux per generatore di numeri casuali
#include "gsl/gsl_randist.h"
#else
#include <gsl/gsl_rng.h> // gnu scientific linux per generatore di numeri casuali
#include <gsl/gsl_randist.h>
#endif

struct PLASMAparams{
	double rminbox[3];
	double rmaxbox[3];
	double ramp_length;
	double density_coefficient;
	double ramp_min_density;
	void *additional_params;
};
#define NBIN_SPECTRUM 1000;
struct SPECIEspectrum{
	double Kmax;
	double Dk;
	int Nbin;
	double *values;
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
	void setDensityCoefficient(double dcoeff);
	void setRampMinDensity(double minden);
	void setAdditionalParams(void* addpar);
	void setMinBox(double xmin, double ymin, double zmin);
	void setMaxBox(double xmax, double ymax, double zmax);
	void setXRangeBox(double xmin, double xmax);
	void setYRangeBox(double ymin, double ymax);
	void setZRangeBox(double zmin, double zmax);
	~PLASMA();
};

//Pre-defined density functions
double box(double x, double y, double z, PLASMAparams plist, double Z, double A);

double left_linear_ramp(double x, double y, double z, PLASMAparams plist, double Z, double A);

double left_soft_ramp(double x, double y, double z, PLASMAparams plist, double Z, double A);

double* rough_box_prepareAdditionalParams(gsl_rng* rng, double roughness, double shift);
double rough_box(double x, double y, double z, PLASMAparams plist, double Z, double A);
double rough_box_edgeCalc(double x0, double y0, double x, double y, int order, double* FFT_params, double roughness);

double box_minus_box(double x, double y, double z, PLASMAparams plist, double Z, double A);

double left_grating(double x, double y, double z, PLASMAparams plist, double Z, double A);

double square_func(double x);
double left_square_grating(double x, double y, double z, PLASMAparams plist, double Z, double A);

double saw_func(double x);
double left_saw_grating(double x, double y, double z, PLASMAparams plist, double Z, double A);


//************** LASER PULSE TYPES *******
enum laserPulseType{ DEFAULT_PULSE, GAUSSIAN, PLANE_WAVE, COS2_PLANE_WAVE };
enum pulsePolarization{ P_POLARIZATION, S_POLARIZATION, CIRCULAR_POLARIZATION };

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

	laserPulse();
	~laserPulse();
	laserPulse(const laserPulse& other);
	laserPulse operator=(const laserPulse& p1);
};

//************** FIELD COORDINATES *******

#define _INTG_CRD 0
#define _HALF_CRD 1;
#define _NULL_CRD -1;

struct integer_or_halfinteger{
	char x;
	char y;
	char z;
};

//************** PARTICLES DISTRIBUTION FUNCTION *******
enum tempDistribType{ WATERBAG, WATERBAG_3TEMP, UNIF_SPHERE, SUPERGAUSSIAN, MAXWELL, JUTTNER };

class tempDistrib{
public:

	tempDistribType type;


	tempDistrib();
	void setWaterbag(double _p0);
	void setWaterbag3Temp(double _p0_x, double _p0_y, double _p0_z);
	void setUnifSphere(double _p0);
	void setSupergaussian(double _p0, double _alpha);
	void setMaxwell(double _temp);
	void setJuttner(double _a);

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

#endif
