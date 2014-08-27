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


#include "structures.h"


/***************************************************************
Se volessi cambiare gli estremi e la definizione della griglia....
definire ghost_region come numero di punti totali in comune
in ACCESSO "N_grid[0]=N_loc[0]+4" cambierebbe per diventare
N_grid[0]=N_loc[0]+ghost_region
****************************************************************/
/***************************************************************
DA FARE:
- diag manager non ha bisogno di farsi passare istep, visto che ha già un puntatore a grid
- i flag devono essere bool !
***************************************************************/


//*************************PLASMA******************************
PLASMA::PLASMA(){
    params.rminbox[0] = params.rminbox[1] = params.rminbox[2] = 0.0;
    params.rmaxbox[0] = params.rmaxbox[1] = params.rmaxbox[2] = 0.0;
    params.left_ramp_length = 0.0;
    params.right_ramp_length = 0.0;
    params.left_scale_length = 1.0;
    params.right_scale_length = 1.0;
    params.density_coefficient = 0.0;
    params.left_ramp_min_density = 0.0;
    params.right_ramp_min_density = 0.0;
    params.additional_params = NULL;
    density_function = NULL;
}

PLASMA::PLASMA(const PLASMA& other)
{
    params.left_ramp_length = other.params.left_ramp_length;
    params.right_ramp_length = other.params.right_ramp_length;
    params.left_scale_length = other.params.left_scale_length;
    params.right_scale_length = other.params.right_scale_length;
    params.density_coefficient = other.params.density_coefficient;
    params.left_ramp_min_density = other.params.left_ramp_min_density;
    params.right_ramp_min_density = other.params.right_ramp_min_density;
    params.additional_params = other.params.additional_params;
    density_function = other.density_function;
    params.rminbox[0] = other.params.rminbox[0];
    params.rminbox[1] = other.params.rminbox[1];
    params.rminbox[2] = other.params.rminbox[2];
    params.rmaxbox[0] = other.params.rmaxbox[0];
    params.rmaxbox[1] = other.params.rmaxbox[1];
    params.rmaxbox[2] = other.params.rmaxbox[2];
}

PLASMA PLASMA::operator=(const PLASMA& p1){
    params.left_ramp_length = p1.params.left_ramp_length;
    params.right_ramp_length = p1.params.right_ramp_length;
    params.left_scale_length = p1.params.left_scale_length;
    params.right_scale_length = p1.params.right_scale_length;
    params.density_coefficient = p1.params.density_coefficient;
    params.left_ramp_min_density = p1.params.left_ramp_min_density;
    params.right_ramp_min_density = p1.params.right_ramp_min_density;
    params.additional_params = p1.params.additional_params;
    density_function = p1.density_function;
    params.rminbox[0] = p1.params.rminbox[0];
    params.rminbox[1] = p1.params.rminbox[1];
    params.rminbox[2] = p1.params.rminbox[2];
    params.rmaxbox[0] = p1.params.rmaxbox[0];
    params.rmaxbox[1] = p1.params.rmaxbox[1];
    params.rmaxbox[2] = p1.params.rmaxbox[2];
    return *this;
}

void PLASMA::setRampLength(double rlength){
    params.left_ramp_length = rlength;
}

void PLASMA::setLeftRampLength(double rlength){
    params.left_ramp_length = rlength;
}

void PLASMA::setRightRampLength(double rlength){
    params.right_ramp_length = rlength;
}

void PLASMA::setScaleLength(double slength){
    params.left_scale_length = slength;
}

void PLASMA::setLeftScaleLength(double slength){
    params.left_scale_length = slength;
}

void PLASMA::setRightScaleLength(double slength){
    params.right_scale_length = slength;
}

void PLASMA::setDensityCoefficient(double dcoeff){
    params.density_coefficient = dcoeff;
}

void PLASMA::setDensityCoefficient(double dcoeff, double lambda){
    params.density_coefficient = dcoeff / (lambda*lambda);
}

void PLASMA::setRampMinDensity(double minden){
    params.left_ramp_min_density = minden;
}

void PLASMA::setLeftRampMinDensity(double minden){
    params.left_ramp_min_density = minden;
}

void PLASMA::setRightRampMinDensity(double minden){
    params.right_ramp_min_density = minden;
}

void PLASMA::setAdditionalParams(void* addpar){
    params.additional_params = addpar;
}

void PLASMA::setMinBox(double xmin, double ymin, double zmin){
    params.rminbox[0] = xmin;
    params.rminbox[1] = ymin;
    params.rminbox[2] = zmin;
}

void PLASMA::setMaxBox(double xmax, double ymax, double zmax){
    params.rmaxbox[0] = xmax;
    params.rmaxbox[1] = ymax;
    params.rmaxbox[2] = zmax;
}

void PLASMA::setXRangeBox(double xmin, double xmax){
    params.rminbox[0] = xmin;
    params.rmaxbox[0] = xmax;
}

void PLASMA::setYRangeBox(double ymin, double ymax){
    params.rminbox[1] = ymin;
    params.rmaxbox[1] = ymax;
}

void PLASMA::setZRangeBox(double zmin, double zmax){
    params.rminbox[2] = zmin;
    params.rmaxbox[2] = zmax;
}


PLASMA::~PLASMA(){
}


double box(double x, double y, double z, PLASMAparams plist, double Z, double A){
    if ((plist.rminbox[0] <= x) && (x <= plist.rmaxbox[0]) &&
        (plist.rminbox[1] <= y) && (y <= plist.rmaxbox[1]) &&
        (plist.rminbox[2] <= z) && (z <= plist.rmaxbox[2])){
        return plist.density_coefficient;
    }
    else{
        return -1;
    }
}

double left_linear_ramp(double x, double y, double z, PLASMAparams plist, double Z, double A){
    if ((plist.rminbox[0] <= x) && (x <= plist.rmaxbox[0]) &&
        (plist.rminbox[1] <= y) && (y <= plist.rmaxbox[1]) &&
        (plist.rminbox[2] <= z) && (z <= plist.rmaxbox[2])){
        if ((x - plist.rminbox[0]) <= plist.left_ramp_length){
            return (plist.density_coefficient - plist.left_ramp_min_density)*(x - plist.rminbox[0]) / plist.left_ramp_length + plist.left_ramp_min_density;
        }
        else{
            return plist.density_coefficient;
        }
    }
    else{
        return -1;
    }
}

double right_linear_ramp(double x, double y, double z, PLASMAparams plist, double Z, double A){
    if ((plist.rminbox[0] <= x) && (x <= plist.rmaxbox[0]) &&
        (plist.rminbox[1] <= y) && (y <= plist.rmaxbox[1]) &&
        (plist.rminbox[2] <= z) && (z <= plist.rmaxbox[2])){
        if ((x) <= (plist.rmaxbox[0] - plist.right_ramp_length)){
            return plist.density_coefficient;
        }
        else{
            return (plist.density_coefficient - plist.right_ramp_min_density)*(plist.rmaxbox[0] - x) / plist.right_ramp_length + plist.right_ramp_min_density;
        }
    }
    else{
        return -1;
    }
}

double left_right_linear_ramp(double x, double y, double z, PLASMAparams plist, double Z, double A){
    if ((plist.rminbox[0] <= x) && (x <= plist.rmaxbox[0]) &&
        (plist.rminbox[1] <= y) && (y <= plist.rmaxbox[1]) &&
        (plist.rminbox[2] <= z) && (z <= plist.rmaxbox[2])){
        if ((x - plist.rminbox[0]) <= plist.left_ramp_length){
            return (plist.density_coefficient - plist.left_ramp_min_density)*(x - plist.rminbox[0]) / plist.left_ramp_length + plist.left_ramp_min_density;
        }
        else if ((x) <= (plist.rmaxbox[0] - plist.right_ramp_length)){
            return plist.density_coefficient;
        }
        else{
            return (plist.density_coefficient - plist.right_ramp_min_density)*(plist.rmaxbox[0] - x) / plist.right_ramp_length + plist.right_ramp_min_density;
        }
    }
    else{
        return -1;
    }
}

double left_fixed_exp_ramp(double x, double y, double z, PLASMAparams plist, double Z, double A){
    if ((plist.rminbox[0] <= x) && (x <= plist.rmaxbox[0]) &&
        (plist.rminbox[1] <= y) && (y <= plist.rmaxbox[1]) &&
        (plist.rminbox[2] <= z) && (z <= plist.rmaxbox[2])){
        if ((x - plist.rminbox[0]) <= plist.left_ramp_length){
            double xx = (x - plist.rminbox[0] - plist.left_ramp_length);
            double densDiff = (plist.density_coefficient - plist.left_ramp_min_density);
            double alpha = densDiff / (1 - exp(-plist.left_ramp_length / plist.left_scale_length));
            double kk = plist.density_coefficient - alpha;
            return (alpha*exp(xx / plist.left_scale_length) + kk);
        }
        else{
            return plist.density_coefficient;
        }
    }
    else{
        return -1;
    }
}

double right_fixed_exp_ramp(double x, double y, double z, PLASMAparams plist, double Z, double A){
    if ((plist.rminbox[0] <= x) && (x <= plist.rmaxbox[0]) &&
        (plist.rminbox[1] <= y) && (y <= plist.rmaxbox[1]) &&
        (plist.rminbox[2] <= z) && (z <= plist.rmaxbox[2])){
        if ((x) <= (plist.rmaxbox[0] - plist.right_ramp_length)){
            return plist.density_coefficient;
        }
        else{
            double xx = (plist.rmaxbox[0] - x - plist.right_ramp_length);
            double densDiff = (plist.density_coefficient - plist.right_ramp_min_density);
            double alpha = densDiff / (1 - exp(-plist.right_ramp_length / plist.right_scale_length));
            double kk = plist.density_coefficient - alpha;
            return (alpha*exp(xx / plist.right_scale_length) + kk);
        }
    }
    else{
        return -1;
    }
}

double left_right_fixed_exp_ramp(double x, double y, double z, PLASMAparams plist, double Z, double A){
    if ((plist.rminbox[0] <= x) && (x <= plist.rmaxbox[0]) &&
        (plist.rminbox[1] <= y) && (y <= plist.rmaxbox[1]) &&
        (plist.rminbox[2] <= z) && (z <= plist.rmaxbox[2])){
        if ((x - plist.rminbox[0]) <= plist.left_ramp_length){
            double xx = (x - plist.rminbox[0] - plist.left_ramp_length);
            double densDiff = (plist.density_coefficient - plist.left_ramp_min_density);
            double alpha = densDiff / (1 - exp(-plist.left_ramp_length / plist.left_scale_length));
            double kk = plist.density_coefficient - alpha;
            return (alpha*exp(xx / plist.left_scale_length) + kk);
        }
        else if ((x) <= (plist.rmaxbox[0] - plist.right_ramp_length)){
            return plist.density_coefficient;
        }
        else{
            double xx = (plist.rmaxbox[0] - x - plist.right_ramp_length);
            double densDiff = (plist.density_coefficient - plist.right_ramp_min_density);
            double alpha = densDiff / (1 - exp(-plist.right_ramp_length / plist.right_scale_length));
            double kk = plist.density_coefficient - alpha;
            return (alpha*exp(xx / plist.right_scale_length) + kk);
        }
    }
    else{
        return -1;
    }
}

double left_free_exp_ramp(double x, double y, double z, PLASMAparams plist, double Z, double A){
    if ((plist.rminbox[0] <= x) && (x <= plist.rmaxbox[0]) &&
        (plist.rminbox[1] <= y) && (y <= plist.rmaxbox[1]) &&
        (plist.rminbox[2] <= z) && (z <= plist.rmaxbox[2])){
        if ((x - plist.rminbox[0]) <= plist.left_ramp_length){
            double xx = (x - plist.rminbox[0] - plist.left_ramp_length);
            return (plist.density_coefficient*exp(xx / plist.left_scale_length));
        }
        else{
            return plist.density_coefficient;
        }
    }
    else{
        return -1;
    }
}

double right_free_exp_ramp(double x, double y, double z, PLASMAparams plist, double Z, double A){
    if ((plist.rminbox[0] <= x) && (x <= plist.rmaxbox[0]) &&
        (plist.rminbox[1] <= y) && (y <= plist.rmaxbox[1]) &&
        (plist.rminbox[2] <= z) && (z <= plist.rmaxbox[2])){
        if ((x) <= (plist.rmaxbox[0] - plist.right_ramp_length)){
            return plist.density_coefficient;
        }
        else{
            double xx = (plist.rmaxbox[0] - x - plist.right_ramp_length);
            return (plist.density_coefficient*exp(xx / plist.left_scale_length));
        }
    }
    else{
        return -1;
    }
}

double left_right_free_exp_ramp(double x, double y, double z, PLASMAparams plist, double Z, double A){
    if ((plist.rminbox[0] <= x) && (x <= plist.rmaxbox[0]) &&
        (plist.rminbox[1] <= y) && (y <= plist.rmaxbox[1]) &&
        (plist.rminbox[2] <= z) && (z <= plist.rmaxbox[2])){
        if ((x - plist.rminbox[0]) <= plist.left_ramp_length){
            double xx = (x - plist.rminbox[0] - plist.left_ramp_length);
            return (plist.density_coefficient*exp(xx / plist.left_scale_length));
        }
        else if ((x) <= (plist.rmaxbox[0] - plist.right_ramp_length)){
            return plist.density_coefficient;
        }
        else{
            double xx = (plist.rmaxbox[0] - x - plist.right_ramp_length);
            return (plist.density_coefficient*exp(+xx / plist.right_scale_length));
        }
    }
    else{
        return -1;
    }
}


double left_soft_ramp(double x, double y, double z, PLASMAparams plist, double Z, double A){
    double lng;
    if ((plist.rminbox[0] <= x) && (x <= plist.rmaxbox[0]) &&
        (plist.rminbox[1] <= y) && (y <= plist.rmaxbox[1]) &&
        (plist.rminbox[2] <= z) && (z <= plist.rmaxbox[2])){
        if ((x - plist.rminbox[0]) <= plist.left_ramp_length){
            lng = ((x - plist.rminbox[0]) / plist.left_ramp_length)*0.5*M_PI;
            return (plist.density_coefficient - plist.left_ramp_min_density)*sin(lng)*sin(lng) + plist.left_ramp_min_density;
        }
        else{
            return plist.density_coefficient;
        }
    }
    else{
        return -1;
    }
}

double left_grating(double x, double y, double z, PLASMAparams plist, double Z, double A){
    double g_y0 = (plist.rmaxbox[1] - plist.rminbox[1])*0.5;
    double* paramlist = (double*)plist.additional_params;
    double g_depth = paramlist[0] * 0.5;
    double g_lambda = paramlist[1];
    double g_phase = paramlist[2];

    double phase = 2.0 * M_PI * ((y - g_y0) + g_phase) / g_lambda;
    double xminbound = plist.rminbox[0] + g_depth*(1.0 - cos(phase));

    if ((xminbound <= x) && (x <= plist.rmaxbox[0]) &&
        (plist.rminbox[1] <= y) && (y <= plist.rmaxbox[1]) &&
        (plist.rminbox[2] <= z) && (z <= plist.rmaxbox[2])){
        if ((x - xminbound) <= plist.left_ramp_length){
            return (plist.density_coefficient - plist.left_ramp_min_density)*(x - xminbound) / plist.left_ramp_length + plist.left_ramp_min_density;
        }
        else{
            return plist.density_coefficient;
        }
    }
    else{
        return -1;
    }
}

double* rough_box_prepareAdditionalParams(gsl_rng* rng, double roughness, double shift){
    const int order = 20;
    double* res = new double[order * 2 + 2];
    for (int i = 0; i < 2 * order; i++){
        res[i] = gsl_ran_flat(rng, 0.0, 2.0*M_PI);
    }
    res[2 * order] = roughness;
    res[2 * order + 1] = shift;
    return res;
}

double rough_box_edgeCalc(double x0, double y0, double x, double y, int order, double* FFT_params, double roughness, double shift){
    x += shift;
    double red_factor;
    for (int i = 0; i < order; i++){
        red_factor = (order - i)*1.0 / (order*1.0);
        x += roughness*red_factor*x0*sin(FFT_params[i] + (y / y0)*2.0*M_PI*(i + 1));
    }
    return x;
}

double rough_box(double x, double y, double z, PLASMAparams plist, double Z, double A){
    const int order = 20;
    double* params = (double*)plist.additional_params;
    double roughness = params[2 * order];
    double shift = params[2 * order + 1];

    double xlimit_left;
    double xlimit_right;

    double xsize = plist.rmaxbox[0] - plist.rminbox[0];
    double ysize = plist.rmaxbox[1] - plist.rminbox[1];

    xlimit_left = rough_box_edgeCalc(xsize,
        ysize,
        plist.rminbox[0], y,
        order,
        params,
        roughness,
        shift);

    xlimit_right = rough_box_edgeCalc(xsize,
        ysize,
        plist.rmaxbox[0], y,
        order,
        params + order,
        roughness,
        -shift);

    if ((xlimit_left <= x) && (x <= xlimit_right) &&
        (plist.rminbox[1] <= y) && (y <= plist.rmaxbox[1]) &&
        (plist.rminbox[2] <= z) && (z <= plist.rmaxbox[2])){
        return plist.density_coefficient;
    }
    else{
        return -1;
    }

}

double box_minus_box(double x, double y, double z, PLASMAparams plist, double Z, double A){
    double* mbox_extrems = (double*)plist.additional_params;
    double mbox_xmin = mbox_extrems[0];
    double mbox_xmax = mbox_extrems[1];
    double mbox_ymin = mbox_extrems[2];
    double mbox_ymax = mbox_extrems[3];
    double mbox_zmin = mbox_extrems[4];
    double mbox_zmax = mbox_extrems[5];

    if ((plist.rminbox[0] <= x) && (x <= plist.rmaxbox[0]) &&
        (plist.rminbox[1] <= y) && (y <= plist.rmaxbox[1]) &&
        (plist.rminbox[2] <= z) && (z <= plist.rmaxbox[2]) &&
        (!((mbox_xmin <= x) && (x <= mbox_xmax) &&
        (mbox_ymin <= y) && (y <= mbox_ymax) &&
        (mbox_zmin <= z) && (z <= mbox_zmax)))
        ){
        return plist.density_coefficient;
    }
    else{
        return -1;
    }
}

double square_func(double x){
    if (cos(x) > 0){
        return 1.0;
    }
    else if (cos(x) < 0){
        return -1;
    }
    else{
        return 0;
    }
}

double left_square_grating(double x, double y, double z, PLASMAparams plist, double Z, double A){
    double g_y0 = (plist.rmaxbox[1] - plist.rminbox[1])*0.5;
    double* paramlist = (double*)plist.additional_params;
    double g_depth = paramlist[0] * 0.5;
    double g_lambda = paramlist[1];
    double g_phase = paramlist[2];

    double phase = 2.0 * M_PI * ((y - g_y0) + g_phase) / g_lambda;
    double xminbound = plist.rminbox[0] + g_depth*(1.0 - square_func(phase));

    if ((xminbound <= x) && (x <= plist.rmaxbox[0]) &&
        (plist.rminbox[1] <= y) && (y <= plist.rmaxbox[1]) &&
        (plist.rminbox[2] <= z) && (z <= plist.rmaxbox[2])){
        if ((x - xminbound) <= plist.left_ramp_length){
            return (plist.density_coefficient - plist.left_ramp_min_density)*(x - xminbound) / plist.left_ramp_length + plist.left_ramp_min_density;
        }
        else{
            return plist.density_coefficient;
        }
    }
    else{
        return -1;
    }
}



//*************************END_PLASMA*****************************
//*************************LASER_PULSE***************************

laserPulse::laserPulse(){
    type = DEFAULT_PULSE;
    polarization = P_POLARIZATION;
    t_FWHM = 0.0;
    waist = 0.0;
    focus_position = 0.0;
    laser_pulse_initial_position = 0.0;
    normalized_amplitude = 0.0;
    lambda0 = 1.0;
    rotation = false;
    angle = 0.0;
    rotation_center_along_x = 0.0;
    rise_time = 0.0;
}

laserPulse::~laserPulse(){}

laserPulse::laserPulse(const laserPulse& other)
    :type(other.type), t_FWHM(other.t_FWHM), waist(other.waist), focus_position(other.focus_position),
    laser_pulse_initial_position(other.laser_pulse_initial_position), normalized_amplitude(other.normalized_amplitude),
    lambda0(other.lambda0), rotation(other.rotation), angle(other.angle),
    rotation_center_along_x(other.rotation_center_along_x), rise_time(other.rise_time){}


laserPulse laserPulse::operator=(const laserPulse& p1){
    type = p1.type;
    polarization = p1.polarization;
    t_FWHM = p1.t_FWHM;
    waist = p1.waist;
    focus_position = p1.focus_position;
    laser_pulse_initial_position = p1.laser_pulse_initial_position;
    normalized_amplitude = p1.normalized_amplitude;
    lambda0 = p1.lambda0;
    rotation = p1.rotation;
    angle = p1.angle;
    rotation_center_along_x = p1.rotation_center_along_x;
    rise_time = p1.rise_time;
    return *this;
}

void laserPulse::setFocusPosition(double _focus_position){
    focus_position = _focus_position;
}

void laserPulse::setPulseInitialPosition(double _laser_pulse_initial_position){
    laser_pulse_initial_position = _laser_pulse_initial_position;
}

void laserPulse::setLambda(double _lambda0){
    lambda0 = _lambda0;
}
void laserPulse::setWaist(double _waist){
    waist = _waist;
}

void laserPulse::setDurationFWHM(double _t_FWHM){
    t_FWHM = _t_FWHM;
}

void laserPulse::setNormalizedAmplitude(double _normalized_amplitude){
    normalized_amplitude = _normalized_amplitude;
}

void laserPulse::setRiseTime(double _rise_time){
    rise_time = _rise_time;
}

void laserPulse::setRotationAngleAndCenter(double _angle, double _rotation_center_along_x){
    angle = _angle;
    rotation_center_along_x = _rotation_center_along_x;
    rotation = true;
}

void laserPulse::setGaussianPulse(double _waist, double _t_FWHM, double _normalized_amplitude){
    type = GAUSSIAN;
    waist = _waist;
    t_FWHM = _t_FWHM;
    normalized_amplitude = _normalized_amplitude;
}

void laserPulse::setPlaneWave(double _normalized_amplitude){
    type = PLANE_WAVE;
    normalized_amplitude = _normalized_amplitude;
}

void laserPulse::setCos2PlaneWave(double _t_FWHM, double _normalized_amplitude){
    type = COS2_PLANE_WAVE;
    t_FWHM = _t_FWHM;
    normalized_amplitude = _normalized_amplitude;
}

void laserPulse::setCos2PlateauPlaneWave(double _t_FWHM, double _rise_time, double _normalized_amplitude){
    type = COS2_PLATEAU_PLANE_WAVE;
    t_FWHM = _t_FWHM;
    normalized_amplitude = _normalized_amplitude;
    rise_time = _rise_time;
}

void laserPulse::setGaussianPulse(){
    type = GAUSSIAN;
}

void laserPulse::setPlaneWave(){
    type = PLANE_WAVE;
}

void laserPulse::setCos2PlaneWave(){
    type = COS2_PLANE_WAVE;
}

void laserPulse::setCos2PlateauPlaneWave(){
    type = COS2_PLATEAU_PLANE_WAVE;
}

void laserPulse::setPPolarization(){
    polarization = P_POLARIZATION;
}

void laserPulse::setSPolarization(){
    polarization = S_POLARIZATION;
}

void laserPulse::setCircularPolarization(){
    polarization = CIRCULAR_POLARIZATION;
}


//************** DISTRIBUTION_FUNCTION **********

tempDistrib::tempDistrib(){
    init = false;
}

bool tempDistrib::isInit(){
    return init;
}

void tempDistrib::setWaterbag(double _p0){
    type = WATERBAG;
    p0 = _p0;
    init = true;
}

void tempDistrib::setWaterbag3Temp(double _p0_x, double _p0_y, double _p0_z){
    type = WATERBAG_3TEMP;
    p0_x = _p0_x;
    p0_y = _p0_y;
    p0_z = _p0_z;
    init = true;
}

void tempDistrib::setUnifSphere(double _p0){
    type = UNIF_SPHERE;
    p0 = _p0;
    init = true;
}

void tempDistrib::setSupergaussian(double _p0, double _alpha){
    type = SUPERGAUSSIAN;
    p0 = _p0;
    alpha = _alpha;
    init = true;
}

void tempDistrib::setMaxwell(double _temp){
    type = MAXWELL;
    temp = _temp;
    init = true;
}

void tempDistrib::setJuttner(double _a){
    type = JUTTNER;
    a = _a;
    init = true;
}

void tempDistrib::setSpecial(double _a){
    type = SPECIAL;
    a = _a;
    init = true;
}

//************** END DISTRIBUTION_FUNCTION ******

