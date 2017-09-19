/*   Copyright 2014-2017 - Andrea Sgattoni, Luca Fedeli, Stefano Sinigardi   */

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


#include "structures.h"
using namespace MUTILS;


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

//*************************PLASMAparams*****************************
PLASMAparams PLASMAparams::operator=(const PLASMAparams& p1) {
  left_ramp_length = p1.left_ramp_length;
  right_ramp_length = p1.right_ramp_length;
  left_scale_length = p1.left_scale_length;
  right_scale_length = p1.right_scale_length;
  density_coefficient = p1.density_coefficient;
  left_ramp_min_density = p1.left_ramp_min_density;
  right_ramp_min_density = p1.right_ramp_min_density;
  additional_params = p1.additional_params;
  rminbox[0] = p1.rminbox[0];
  rminbox[1] = p1.rminbox[1];
  rminbox[2] = p1.rminbox[2];
  rmaxbox[0] = p1.rmaxbox[0];
  rmaxbox[1] = p1.rmaxbox[1];
  rmaxbox[2] = p1.rmaxbox[2];
  spheres = p1.spheres;
  FFTplasma = p1.FFTplasma;
  amIRNDwire = p1.amIRNDwire;
  return *this;
}

//*************************PLASMA******************************

const std::string PLASMA::dFNames[] = {
  "box",
  "left_linear_ramp",
  "right_linear_ramp",
  "left_right_linear_ramp",
  "left_fixed_exp_ramp",
  "right_fixed_exp_ramp",
  "left_right_fixed_exp_ramp",
  "left_free_exp_ramp",
  "right_free_exp_ramp",
  "left_right_free_exp_ramp",
  "left soft ramp",
  "rough_box",
  "box_minus_box",
  "left_grating",
  "left_square_grating",
  "guide",
  "modGrat",
  "spoofGrat",
  "spheres",
  "left_blazed_grating",
  "pillars2D",
  "nanotubes2D",
  "foils2D",
  "res2D",
  "user1",
  "user2",
  "fftplasma",
  "cylinder",
  "left_grating_exp_ramp",
  "rand_wires",
  "pillars3D"
};
const distrib_function PLASMA::dFPoint[] = {
  box,
  left_linear_ramp,
  right_linear_ramp,
  left_right_linear_ramp,
  left_fixed_exp_ramp,
  right_fixed_exp_ramp,
  left_right_fixed_exp_ramp,
  left_free_exp_ramp,
  right_free_exp_ramp,
  left_right_free_exp_ramp,
  left_soft_ramp,
  rough_box,
  box_minus_box,
  left_grating,
  left_square_grating,
  guide,
  modGrat,
  spoofGrat,
  spheres,
  left_blazed_grating,
  pillars2D,
  nanotubes2D,
  foils2D,
  demo_2D_resonator,
  user1,
  user2,
  fftplasma,
  cylinder,
  left_grating_exp_ramp,
  rand_wires,
  pillars3D
};

bool PLASMA::isRndWir(int dfIndex){
    return (dfIndex == 29);

}

bool PLASMA::isGrating(int dfIndex) {
  if (dfIndex == 13 || dfIndex == 14 || dfIndex == 19 || dfIndex == 28)
    return true;
  else
    return false;
}

bool PLASMA::isBlazedGrating(int dfIndex) {
  if (dfIndex == 19)
    return true;
  else
    return false;
}

bool PLASMA::isPillar2D(int dfIndex) {
  if (dfIndex == 20)
    return true;
  else
    return false;
}


bool PLASMA::isNanotubes2D(int dfIndex) {
  if (dfIndex == 21)
    return true;
  else
    return false;
}

bool PLASMA::isFoils2D(int dfIndex) {
  if (dfIndex == 22)
    return true;
  else
    return false;
}

bool PLASMA::isr1r1(int dfIndex) {
  if (dfIndex == 24)
    return true;
  else
    return false;
}

bool PLASMA::isUser1(int dfIndex) {
  if (dfIndex == 24)
    return true;
  else
    return false;
}


bool PLASMA::isUser2(int dfIndex) {
  if (dfIndex == 25)
    return true;
  else
    return false;
}

bool PLASMA::isPillar3D(int dfIndex){
    if (dfIndex == 30)
      return true;
    else
      return false;
}



PLASMA::PLASMA() {
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
  params.amIRNDwire = false;
  density_function = NULL;
}

PLASMA::PLASMA(const PLASMA& other)
{
  density_function = other.density_function;
  params = other.params;

}

PLASMA PLASMA::operator=(const PLASMA& p1) {
  density_function = p1.density_function;
  params = p1.params; 
  return *this;
}

void PLASMA::setRampLength(double rlength) {
  params.left_ramp_length = rlength;
}

void PLASMA::setLeftRampLength(double rlength) {
  params.left_ramp_length = rlength;
}

void PLASMA::setRightRampLength(double rlength) {
  params.right_ramp_length = rlength;
}


void PLASMA::setScaleLength(double slength) {
  params.left_scale_length = slength;
}

void PLASMA::setLeftScaleLength(double slength) {
  params.left_scale_length = slength;
}

void PLASMA::setRightScaleLength(double slength) {
  params.right_scale_length = slength;
}

void PLASMA::setDensityCoefficient(double dcoeff) {
  params.density_coefficient = dcoeff;
}

void PLASMA::setDensityCoefficient(double dcoeff, double lambda) {
  params.density_coefficient = dcoeff / (lambda*lambda);
}

void PLASMA::setRampMinDensity(double minden) {
  params.left_ramp_min_density = minden;
}

void PLASMA::setLeftRampMinDensity(double minden) {
  params.left_ramp_min_density = minden;
}

void PLASMA::setRightRampMinDensity(double minden) {
  params.right_ramp_min_density = minden;
}

void PLASMA::setAdditionalParams(void* addpar) {
  params.additional_params = addpar;
}

void PLASMA::setMinBox(double xmin, double ymin, double zmin) {
  params.rminbox[0] = xmin;
  params.rminbox[1] = ymin;
  params.rminbox[2] = zmin;
}

void PLASMA::setMaxBox(double xmax, double ymax, double zmax) {
  params.rmaxbox[0] = xmax;
  params.rmaxbox[1] = ymax;
  params.rmaxbox[2] = zmax;
}

void PLASMA::setXRangeBox(double xmin, double xmax) {
  params.rminbox[0] = xmin;
  params.rmaxbox[0] = xmax;
}

void PLASMA::setYRangeBox(double ymin, double ymax) {
  params.rminbox[1] = ymin;
  params.rmaxbox[1] = ymax;
}

void PLASMA::setZRangeBox(double zmin, double zmax) {
  params.rminbox[2] = zmin;
  params.rmaxbox[2] = zmax;
}

void PLASMA::setRNDwire(bool val){
    params.amIRNDwire = val;
}

bool PLASMA::getRNDwire(){
    return params.amIRNDwire;
}



PLASMA::~PLASMA() {    
}



double guide(double x, double y, double z, PLASMAparams plist, double Z, double A) {
  double g_x0 = 0.0;
  double g_x1 = 10.0;
  double g_x2 = 20.0;

  double g_depth = 0.250 * 0.5;
  double g_lambda = 2.0;

  double phase = 2.0*M_PI*(x - g_x0) / g_lambda;
  double yminbound1 = -10.0;
  double ymaxbound1 = -9.0 + g_depth*(1.0 - cos(phase));
  double ymaxbound2 = +10.0;
  double yminbound2 = +9.0 + g_depth*(1.0 - cos(phase));

  if ((plist.rminbox[0] <= x) && (x <= plist.rmaxbox[0]) &&
    (plist.rminbox[1] <= y) && (y <= plist.rmaxbox[1]) &&
    (plist.rminbox[2] <= z) && (z <= plist.rmaxbox[2])) {
    if ((x > g_x0) && (x < g_x1)) {
      if (y< ymaxbound1 && y > yminbound1)
        return plist.density_coefficient;
      else if (y < ymaxbound2 && y > yminbound2)
        return plist.density_coefficient;
    }
    if ((x > g_x1) && (x < g_x2)) {
      double uplimitmax = 10.0 + (g_x1 - x) / (g_x2 - g_x1)*9.0;
      double uplimitmin = 9.0 + (g_x1 - x) / (g_x2 - g_x1)*9.0;
      double downlimitmin = -10.0 + (x - g_x1) / (g_x2 - g_x1)*9.0;
      double downlimitmax = -9.0 + (x - g_x1) / (g_x2 - g_x1)*9.0;
      if (y < uplimitmax && y > uplimitmin)
        return plist.density_coefficient;
      if (y > downlimitmin && y < downlimitmax)
        return plist.density_coefficient;
    }
  }
  return -1.0;
}

double modGrat(double x, double y, double z, PLASMAparams plist, double Z, double A) {
  double g_y0 = (plist.rmaxbox[1] - plist.rminbox[1])*0.5;
  double g_depth_1 = 0.125;
  double g_lambda_1 = 1.0;
  double g_depth_2 = 0.1;
  double g_lambda_2 = 0.1;


  double phase1 = 2.0 * M_PI * ((y - g_y0)) / g_lambda_1;
  double phase2 = 2.0 * M_PI * ((y - g_y0)) / g_lambda_2;

  double xminbound = plist.rminbox[0] + g_depth_1*(1.0 - cos(phase1)) + g_depth_2*(1.0 - cos(phase2));

  if ((xminbound <= x) && (x <= plist.rmaxbox[0]) &&
    (plist.rminbox[1] <= y) && (y <= plist.rmaxbox[1]) &&
    (plist.rminbox[2] <= z) && (z <= plist.rmaxbox[2])) {
    return plist.density_coefficient;
  }
  else {
    return -1;
  }

}

double spoofGrat(double x, double y, double z, PLASMAparams plist, double Z, double A) {
  double g_y0 = (plist.rmaxbox[1] - plist.rminbox[1])*0.5;
  double g_depth = 0.5;
  double g_a = 0.125;
  double g_d = 0.25;

  bool flag = false;
  double rem = fmod(fabs(y - g_y0), (g_d*0.5));
  if (rem <= g_a / 2)
    flag = true;

  if ((plist.rminbox[0] <= x) && (x <= plist.rmaxbox[0]) &&
    (plist.rminbox[1] <= y) && (y <= plist.rmaxbox[1]) &&
    (plist.rminbox[2] <= z) && (z <= plist.rmaxbox[2])) {
    if (x - plist.rminbox[0] > g_depth)
      return plist.density_coefficient;
    else {
      if (!flag)
        return plist.density_coefficient;
      else
        return -1;
    }
  }
  else {
    return -1;
  }

}


double box(double x, double y, double z, PLASMAparams plist, double Z, double A) {
  if ((plist.rminbox[0] <= x) && (x <= plist.rmaxbox[0]) &&
    (plist.rminbox[1] <= y) && (y <= plist.rmaxbox[1]) &&
    (plist.rminbox[2] <= z) && (z <= plist.rmaxbox[2])) {
    return plist.density_coefficient;
  }
  else {
    return -1;
  }
}

double left_linear_ramp(double x, double y, double z, PLASMAparams plist, double Z, double A) {
  if ((plist.rminbox[0] <= x) && (x <= plist.rmaxbox[0]) &&
    (plist.rminbox[1] <= y) && (y <= plist.rmaxbox[1]) &&
    (plist.rminbox[2] <= z) && (z <= plist.rmaxbox[2])) {
    if ((x - plist.rminbox[0]) <= plist.left_ramp_length) {
      return (plist.density_coefficient - plist.left_ramp_min_density)*(x - plist.rminbox[0]) / plist.left_ramp_length + plist.left_ramp_min_density;
    }
    else {
      return plist.density_coefficient;
    }
  }
  else {
    return -1;
  }
}

double right_linear_ramp(double x, double y, double z, PLASMAparams plist, double Z, double A) {
  if ((plist.rminbox[0] <= x) && (x <= plist.rmaxbox[0]) &&
    (plist.rminbox[1] <= y) && (y <= plist.rmaxbox[1]) &&
    (plist.rminbox[2] <= z) && (z <= plist.rmaxbox[2])) {
    if ((x) <= (plist.rmaxbox[0] - plist.right_ramp_length)) {
      return plist.density_coefficient;
    }
    else {
      return (plist.density_coefficient - plist.right_ramp_min_density)*(plist.rmaxbox[0] - x) / plist.right_ramp_length + plist.right_ramp_min_density;
    }
  }
  else {
    return -1;
  }
}

double left_right_linear_ramp(double x, double y, double z, PLASMAparams plist, double Z, double A) {
  if ((plist.rminbox[0] <= x) && (x <= plist.rmaxbox[0]) &&
    (plist.rminbox[1] <= y) && (y <= plist.rmaxbox[1]) &&
    (plist.rminbox[2] <= z) && (z <= plist.rmaxbox[2])) {
    if ((x - plist.rminbox[0]) <= plist.left_ramp_length) {
      return (plist.density_coefficient - plist.left_ramp_min_density)*(x - plist.rminbox[0]) / plist.left_ramp_length + plist.left_ramp_min_density;
    }
    else if ((x) <= (plist.rmaxbox[0] - plist.right_ramp_length)) {
      return plist.density_coefficient;
    }
    else {
      return (plist.density_coefficient - plist.right_ramp_min_density)*(plist.rmaxbox[0] - x) / plist.right_ramp_length + plist.right_ramp_min_density;
    }
  }
  else {
    return -1;
  }
}

double left_fixed_exp_ramp(double x, double y, double z, PLASMAparams plist, double Z, double A) {
  if ((plist.rminbox[0] <= x) && (x <= plist.rmaxbox[0]) &&
    (plist.rminbox[1] <= y) && (y <= plist.rmaxbox[1]) &&
    (plist.rminbox[2] <= z) && (z <= plist.rmaxbox[2])) {
    if ((x - plist.rminbox[0]) <= plist.left_ramp_length) {
      double xx = (x - plist.rminbox[0] - plist.left_ramp_length);
      double densDiff = (plist.density_coefficient - plist.left_ramp_min_density);
      double alpha = densDiff / (1 - exp(-plist.left_ramp_length / plist.left_scale_length));
      double kk = plist.density_coefficient - alpha;
      return (alpha*exp(xx / plist.left_scale_length) + kk);
    }
    else {
      return plist.density_coefficient;
    }
  }
  else {
    return -1;
  }
}

double right_fixed_exp_ramp(double x, double y, double z, PLASMAparams plist, double Z, double A) {
  if ((plist.rminbox[0] <= x) && (x <= plist.rmaxbox[0]) &&
    (plist.rminbox[1] <= y) && (y <= plist.rmaxbox[1]) &&
    (plist.rminbox[2] <= z) && (z <= plist.rmaxbox[2])) {
    if ((x) <= (plist.rmaxbox[0] - plist.right_ramp_length)) {
      return plist.density_coefficient;
    }
    else {
      double xx = (plist.rmaxbox[0] - x - plist.right_ramp_length);
      double densDiff = (plist.density_coefficient - plist.right_ramp_min_density);
      double alpha = densDiff / (1 - exp(-plist.right_ramp_length / plist.right_scale_length));
      double kk = plist.density_coefficient - alpha;
      return (alpha*exp(xx / plist.right_scale_length) + kk);
    }
  }
  else {
    return -1;
  }
}


double cylinder(double x, double y, double z, PLASMAparams plist, double Z, double A) {
  if ((plist.rminbox[0] <= x) && (x <= plist.rmaxbox[0]) &&
    (plist.rminbox[1] <= y) && (y <= plist.rmaxbox[1]) &&
    (plist.rminbox[2] <= z) && (z <= plist.rmaxbox[2])) {
    double cy = 0.5*(plist.rminbox[1]+plist.rmaxbox[1]);
    double cz = 0.5*(plist.rminbox[2]+plist.rmaxbox[2]);
    double exty = plist.rmaxbox[1]-plist.rminbox[1];
     double extz = plist.rmaxbox[2]-plist.rminbox[2];
    double radius = 0.5*MIN(exty,extz);
    double rr = (y - cy)*(y - cy) + (z -cz)*(z-cz);
    if (rr <= radius*radius) {
      return plist.density_coefficient;
    }
    else {
      return -1;
    }
  }
  else {
    return -1;
  }
}

double left_right_fixed_exp_ramp(double x, double y, double z, PLASMAparams plist, double Z, double A) {
  if ((plist.rminbox[0] <= x) && (x <= plist.rmaxbox[0]) &&
    (plist.rminbox[1] <= y) && (y <= plist.rmaxbox[1]) &&
    (plist.rminbox[2] <= z) && (z <= plist.rmaxbox[2])) {
    if ((x - plist.rminbox[0]) <= plist.left_ramp_length) {
      double xx = (x - plist.rminbox[0] - plist.left_ramp_length);
      double densDiff = (plist.density_coefficient - plist.left_ramp_min_density);
      double alpha = densDiff / (1 - exp(-plist.left_ramp_length / plist.left_scale_length));
      double kk = plist.density_coefficient - alpha;
      return (alpha*exp(xx / plist.left_scale_length) + kk);
    }
    else if ((x) <= (plist.rmaxbox[0] - plist.right_ramp_length)) {
      return plist.density_coefficient;
    }
    else {
      double xx = (plist.rmaxbox[0] - x - plist.right_ramp_length);
      double densDiff = (plist.density_coefficient - plist.right_ramp_min_density);
      double alpha = densDiff / (1 - exp(-plist.right_ramp_length / plist.right_scale_length));
      double kk = plist.density_coefficient - alpha;
      return (alpha*exp(xx / plist.right_scale_length) + kk);
    }
  }
  else {
    return -1;
  }
}

double left_free_exp_ramp(double x, double y, double z, PLASMAparams plist, double Z, double A) {
  if ((plist.rminbox[0] <= x) && (x <= plist.rmaxbox[0]) &&
    (plist.rminbox[1] <= y) && (y <= plist.rmaxbox[1]) &&
    (plist.rminbox[2] <= z) && (z <= plist.rmaxbox[2])) {
    if ((x - plist.rminbox[0]) <= plist.left_ramp_length) {
      double xx = (x - plist.rminbox[0] - plist.left_ramp_length);
      return (plist.density_coefficient*exp(xx / plist.left_scale_length));
    }
    else {
      return plist.density_coefficient;
    }
  }
  else {
    return -1;
  }
}

double right_free_exp_ramp(double x, double y, double z, PLASMAparams plist, double Z, double A) {
  if ((plist.rminbox[0] <= x) && (x <= plist.rmaxbox[0]) &&
    (plist.rminbox[1] <= y) && (y <= plist.rmaxbox[1]) &&
    (plist.rminbox[2] <= z) && (z <= plist.rmaxbox[2])) {
    if ((x) <= (plist.rmaxbox[0] - plist.right_ramp_length)) {
      return plist.density_coefficient;
    }
    else {
      double xx = (plist.rmaxbox[0] - x - plist.right_ramp_length);
      return (plist.density_coefficient*exp(xx / plist.right_scale_length));
    }
  }
  else {
    return -1;
  }
}

double left_right_free_exp_ramp(double x, double y, double z, PLASMAparams plist, double Z, double A) {
  if ((plist.rminbox[0] <= x) && (x <= plist.rmaxbox[0]) &&
    (plist.rminbox[1] <= y) && (y <= plist.rmaxbox[1]) &&
    (plist.rminbox[2] <= z) && (z <= plist.rmaxbox[2])) {
    if ((x - plist.rminbox[0]) <= plist.left_ramp_length) {
      double xx = (x - plist.rminbox[0] - plist.left_ramp_length);
      return (plist.density_coefficient*exp(xx / plist.left_scale_length));
    }
    else if ((x) <= (plist.rmaxbox[0] - plist.right_ramp_length)) {
      return plist.density_coefficient;
    }
    else {
      double xx = (plist.rmaxbox[0] - x - plist.right_ramp_length);
      return (plist.density_coefficient*exp(+xx / plist.right_scale_length));
    }
  }
  else {
    return -1;
  }
}


double left_soft_ramp(double x, double y, double z, PLASMAparams plist, double Z, double A) {
  double lng;
  if ((plist.rminbox[0] <= x) && (x <= plist.rmaxbox[0]) &&
    (plist.rminbox[1] <= y) && (y <= plist.rmaxbox[1]) &&
    (plist.rminbox[2] <= z) && (z <= plist.rmaxbox[2])) {
    if ((x - plist.rminbox[0]) <= plist.left_ramp_length) {
      lng = ((x - plist.rminbox[0]) / plist.left_ramp_length)*0.5*M_PI;
      return (plist.density_coefficient - plist.left_ramp_min_density)*sin(lng)*sin(lng) + plist.left_ramp_min_density;
    }
    else {
      return plist.density_coefficient;
    }
  }
  else {
    return -1;
  }
}

double left_grating(double x, double y, double z, PLASMAparams plist, double Z, double A) {
  double g_y0 = (plist.rmaxbox[1] - plist.rminbox[1])*0.5;
  double* paramlist = (double*)plist.additional_params;
  double g_depth = paramlist[0] * 0.5;
  double g_lambda = paramlist[1];
  double g_phase = paramlist[2];

  double phase = 2.0 * M_PI * ((y - g_y0) + g_phase) / g_lambda;
  double xminbound = plist.rminbox[0] + g_depth*(1.0 - cos(phase));

  if ((xminbound <= x) && (x <= plist.rmaxbox[0]) &&
    (plist.rminbox[1] <= y) && (y <= plist.rmaxbox[1]) &&
    (plist.rminbox[2] <= z) && (z <= plist.rmaxbox[2])) {
    if ((x - xminbound) <= plist.left_ramp_length) {
      return (plist.density_coefficient - plist.left_ramp_min_density)*(x - xminbound) / plist.left_ramp_length + plist.left_ramp_min_density;
    }
    else {
      return plist.density_coefficient;
    }
  }
  else {
    return -1;
  }
}

double left_grating_exp_ramp(double x, double y, double z, PLASMAparams plist, double Z, double A) {
  double g_y0 = (plist.rmaxbox[1] - plist.rminbox[1])*0.5;
  double* paramlist = (double*)plist.additional_params;
  double g_depth = paramlist[0] * 0.5;
  double g_lambda = paramlist[1];
  double g_phase = paramlist[2];

  double phase = 2.0 * M_PI * ((y - g_y0) + g_phase) / g_lambda;
  double xminbound = plist.rminbox[0] + g_depth*(1.0 - cos(phase));

  if ((xminbound <= x) && (x <= plist.rmaxbox[0]) &&
    (plist.rminbox[1] <= y) && (y <= plist.rmaxbox[1]) &&
    (plist.rminbox[2] <= z) && (z <= plist.rmaxbox[2])) {
    if ((x - xminbound) <= plist.left_ramp_length) {
      double xx = (x - xminbound - plist.left_ramp_length);
      return (plist.density_coefficient*exp(xx / plist.left_scale_length));
    }
    else {
      return plist.density_coefficient;
    }
  }
  else {
    return -1;
  }
}

double left_blazed_grating(double x, double y, double z, PLASMAparams plist, double Z, double A) {
  double g_y0 = (plist.rmaxbox[1] - plist.rminbox[1])*0.5;
  double* paramlist = (double*)plist.additional_params;
  //double g_depth = paramlist[0];
  double g_lambda = paramlist[1];
  double g_alpha = paramlist[0]/ 180.0 * M_PI;
  double g_phase = paramlist[2];
  double mySign = 1;
  if(g_alpha<0){
    mySign=-1;
    g_alpha= -g_alpha;
  }
  g_y0 += g_phase/(2*M_PI)*g_lambda;
  double g_ah = g_lambda*cos(g_alpha)*cos(g_alpha);
  double g_ch = g_lambda -g_ah;
  double g_depth = g_lambda*cos(g_alpha)*sin(g_alpha);
  double yy = (mySign*y - g_y0) / g_lambda;
  yy = (yy - floor(yy))*g_lambda;

  double xminbound;

  if (yy < g_ah) {
    xminbound = plist.rminbox[0] + g_depth*(yy / g_ah);
  }
  else {
    xminbound = plist.rminbox[0] + g_depth*((g_lambda - yy) / g_ch);
  }

  if ((xminbound <= x) && (x <= plist.rmaxbox[0]) &&
    (plist.rminbox[1] <= y) && (y <= plist.rmaxbox[1]) &&
    (plist.rminbox[2] <= z) && (z <= plist.rmaxbox[2])) {
    return plist.density_coefficient;
  }
  else {
    return -1;
  }
}


double* rough_box_prepareAdditionalParams(my_rng_generator& rng, double roughness, double shift) {
  const int order = 20;
  double* res = new double[order * 2 + 2];
  my_uniform_real_distribution dist(0.0, 2.0*M_PI);
  for (int i = 0; i < 2 * order; i++) {
    res[i] = dist(rng);
  }
  res[2 * order] = roughness;
  res[2 * order + 1] = shift;
  return res;
}

double rough_box_edgeCalc(double x0, double y0, double x, double y, int order, double* FFT_params, double roughness, double shift) {
  x += shift;
  double red_factor;
  for (int i = 0; i < order; i++) {
    red_factor = (order - i)*1.0 / (order*1.0);
    x += roughness*red_factor*x0*sin(FFT_params[i] + (y / y0)*2.0*M_PI*(i + 1));
  }
  return x;
}

double rough_box(double x, double y, double z, PLASMAparams plist, double Z, double A) {
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
    (plist.rminbox[2] <= z) && (z <= plist.rmaxbox[2])) {
    return plist.density_coefficient;
  }
  else {
    return -1;
  }

}

double box_minus_box(double x, double y, double z, PLASMAparams plist, double Z, double A) {
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
    ) {
    return plist.density_coefficient;
  }
  else {
    return -1;
  }
}

double square_func(double x) {
  if (cos(x) > 0) {
    return 1.0;
  }
  else if (cos(x) < 0) {
    return -1;
  }
  else {
    return 0;
  }
}

double left_square_grating(double x, double y, double z, PLASMAparams plist, double Z, double A) {
  double g_y0 = (plist.rmaxbox[1] - plist.rminbox[1])*0.5;
  double* paramlist = (double*)plist.additional_params;
  double g_depth = paramlist[0] * 0.5;
  double g_lambda = paramlist[1];
  double g_phase = paramlist[2];

  double phase = 2.0 * M_PI * ((y - g_y0) + g_phase) / g_lambda;
  double xminbound = plist.rminbox[0] + g_depth*(1.0 - square_func(phase));

  if ((xminbound <= x) && (x <= plist.rmaxbox[0]) &&
    (plist.rminbox[1] <= y) && (y <= plist.rmaxbox[1]) &&
    (plist.rminbox[2] <= z) && (z <= plist.rmaxbox[2])) {
    if ((x - xminbound) <= plist.left_ramp_length) {
      return (plist.density_coefficient - plist.left_ramp_min_density)*(x - xminbound) / plist.left_ramp_length + plist.left_ramp_min_density;
    }
    else {
      return plist.density_coefficient;
    }
  }
  else {
    return -1;
  }
}

void setCoordWithinBoundaries(double &x, double min, double max) {
  double box = max - min;
  if (x < min) {
    x += box * ((int)((max - x) / box));
  }
  if (x > max) {
    x -= box* ((int)((x - min) / box));
  }
}

double spheres(double x, double y, double z, PLASMAparams plist, double Z, double A) {
  SPHERES *myspheres = plist.spheres;
  int Nspheres = myspheres->NSpheres;
  double spDensity = plist.density_coefficient / myspheres->fillingFactor;
  double value = 0, xsp, ysp, zsp, radius;
  //double extrems[2];
  //std::cout << " Nspheres = " << Nspheres << std::endl;

  //std::cout << " density_coefficient = " << plist.density_coefficient << std::endl;
  //std::cout << "  fillingFactor = " << plist.spheres->fillingFactor << std::endl;
  //std::cout << " spheres density = " << spDensity << std::endl;
  if ((plist.rminbox[0] <= x) && (x <= plist.rmaxbox[0])) {
    if ((plist.rminbox[1] <= y) && (y <= plist.rmaxbox[1])) {
      if ((plist.rminbox[2] <= z) && (z <= plist.rmaxbox[2])) {
        setCoordWithinBoundaries(y, myspheres->rmin[1], myspheres->rmax[1]);
        setCoordWithinBoundaries(z, myspheres->rmin[2], myspheres->rmax[2]);
        for (int i = 0; i < Nspheres; i++) {
          radius = myspheres->coords[i * 4 + 3];
          xsp = plist.rmaxbox[0] - myspheres->coords[i * 4 + 0] - radius;
          ysp = myspheres->coords[i * 4 + 1];
          zsp = myspheres->coords[i * 4 + 2];
          double distance = (xsp - x)*(xsp - x) + (ysp - y)*(ysp - y) + (zsp - z)*(zsp - z);
          if (distance <= (radius*radius)){
            value = spDensity;
            break;
          }

        }
      }
    }
  }
  else {
    value = -1;
  }
  return value;
}


double pillars2D (double x, double y, double z, PLASMAparams plist, double Z, double A){
    double *paramlist = (double*)plist.additional_params;
    double dx = paramlist[0];
    double dy = paramlist[1];
    double r = paramlist[2];
    double rho = -1;

    int Nx;
    int Ny;

    Nx = (plist.rmaxbox[0]-plist.rminbox[0]-2*r+dx)/dx;
    Ny = (plist.rmaxbox[1]-plist.rminbox[1]-2*r+dy)/dy;

    if ( (plist.rminbox[0] <= x) && (x <= plist.rmaxbox[0]) &&
         (plist.rminbox[1] <= y) && (y <= plist.rmaxbox[1]) &&
         (plist.rminbox[2] <= z) && (z <= plist.rmaxbox[2]) ){

        for(int i=0; i<Nx; i++){
            for(int j=0; j<Ny; j++){
                if( fabs( pow(x-plist.rminbox[0]-r-i*dx,2) + pow(y-plist.rminbox[1]-r-j*dy,2)  ) <= r*r ){
                    rho = plist.density_coefficient;
                }
            }
        }
    }
    return rho;
}


double nanotubes2D (double x, double y, double z, PLASMAparams plist, double Z, double A){
    double *paramlist = (double*)plist.additional_params;
    double width = paramlist[0];
    double dist = paramlist[1];
    double depth = paramlist[2];
    double rho = -1;
    if( (plist.rminbox[0] <= x) && (x <= plist.rmaxbox[0]) &&
        (plist.rminbox[1] <= y) && (y <= plist.rmaxbox[1]) &&
        (plist.rminbox[2] <= z) && (z <= plist.rmaxbox[2]) ){

        int r = (plist.rmaxbox[1] - plist.rminbox[1])/(width+dist);
        int i = 1;
            while ( i <= (r+1) ){
                if( (x < plist.rminbox[0]+depth) && fabs( (y-plist.rminbox[1]) - ((2*i-1)*width*0.5+(i-1)*dist) ) <= width*0.5 ){
                    rho = plist.density_coefficient;
                }
                i+=1;
            }
            if(x >= (plist.rminbox[0]+depth)){
                rho = plist.density_coefficient;
            }
    }
    return rho;
}

double foils2D (double x, double y, double z, PLASMAparams plist, double Z, double A){
    double *paramlist = (double*)plist.additional_params;
    double width = paramlist[0];
    double dist = paramlist[1];
    double rho = -1;
    if( (plist.rminbox[0] <= x) && (x <= plist.rmaxbox[0]) &&
        (plist.rminbox[1] <= y) && (y <= plist.rmaxbox[1])  &&
        (plist.rminbox[2] <= z) && (z <= plist.rmaxbox[2]) ){

        int r = (plist.rmaxbox[0] - plist.rminbox[0])/(width+dist);
        int i = 1;
            while ( i <= (r+1) ){
                if( fabs( (x-plist.rminbox[0]) - ((2*i-1)*width*0.5+(i-1)*dist) ) <= width*0.5 ){
                    rho = plist.density_coefficient;
                }
                i+=1;
            }
    }
    return rho;
}

double demo_2D_resonator (double x, double y, double z, PLASMAparams plist, double Z, double A){
    double rho = -1;
    if( (plist.rminbox[0] <= x) && (x <= plist.rmaxbox[0]) &&
        (plist.rminbox[1] <= y) && (y <= plist.rmaxbox[1])  &&
        (plist.rminbox[2] <= z) && (z <= plist.rmaxbox[2]) ){

        if(x> 0.0 && x < 0.6 && y < 1 && y > -1){
            rho = -1;
        }
        else if( ( (x-3.25)*(x-3.25)+y*y < 3*3 ) && x < 3.25){
            rho = -1;
        }
        else if (x> 3.25 && x < 6.25 && y < 3 && y > -3 ){
            rho = -1;
        }
        else if (x > 6.25 && ( (x-6.25)*(x-6.25)+y*y < 3*3) ){
            rho = -1;
        }
        else{
            rho = plist.density_coefficient;
        }
    }
    return rho;
}

double user1(double x, double y, double z, PLASMAparams plist, double Z, double A){

  double rho = -1;
  double *paramlist = (double*)plist.additional_params;
  double width = paramlist[0];
  double depth = paramlist[1];
  double position = paramlist[2];
  double temperature = paramlist[3];


  if( (plist.rminbox[0] <= x) && (x <= plist.rmaxbox[0]) &&
      (plist.rminbox[1] <= y) && (y <= plist.rmaxbox[1])  &&
      (plist.rminbox[2] <= z) && (z <= plist.rmaxbox[2]) ){
    if(Z==-1){
      double XX = (x-position);
      double myexp = depth*exp(-(XX*XX)/(width*width));
      rho = plist.density_coefficient*(1-myexp);
    }
    else{
      double XX = (x-position);
      double myexp = depth*exp(-(XX*XX)/(width*width));
      double ne = plist.density_coefficient*(1-myexp);
      double dne = plist.density_coefficient*2*XX/(width*width)*myexp;
      double d2ne = plist.density_coefficient*2/(width*width)*( 1 - 2*XX*XX/(width*width) )*myexp;

      rho = ne + temperature/((2 * M_PI)*(2 * M_PI))*( (dne*dne)/(ne*ne) - d2ne/ne );
    }
  }
  return rho;
}

double user2(double x, double y, double z, PLASMAparams plist, double Z, double A){

  double rho = -1;
  double *paramlist = (double*)plist.additional_params;
  double width = paramlist[0];
  double depth = paramlist[1];
  double position = paramlist[2];
  double temperature = paramlist[3];


  if( (plist.rminbox[0] <= x) && (x <= plist.rmaxbox[0]) &&
      (plist.rminbox[1] <= y) && (y <= plist.rmaxbox[1])  &&
      (plist.rminbox[2] <= z) && (z <= plist.rmaxbox[2]) ){
    if(Z==-1){
      double XX = (x-position);
      double YY = y;
      double myexp = depth*exp(-(XX*XX+YY*YY)/(width*width));
      rho = plist.density_coefficient*(1-myexp);
    }
    else{
      double XX = (x-position);
      double YY = y;
      double myexp = depth*exp(-(XX*XX+YY*YY)/(width*width));
      double ne = plist.density_coefficient*(1-myexp);
      double dne = plist.density_coefficient*2*XX/(width*width)*myexp;
      dne += plist.density_coefficient*2*YY/(width*width)*myexp;
      double d2ne = plist.density_coefficient*2/(width*width)*( 1 - 2*XX*XX/(width*width) )*myexp;
      d2ne += plist.density_coefficient*2/(width*width)*( 1 - 2*YY*YY/(width*width) )*myexp;

      rho = ne + temperature/((2 * M_PI)*(2 * M_PI))*( (dne*dne)/(ne*ne) - d2ne/ne );
    }
  }
  return rho;
}

double fftplasma(double x, double y, double z, PLASMAparams plist, double Z, double A) {
  FFTPLASMA *myfft = plist.FFTplasma;
  int numcomp = myfft->numcomp;


  double val;

  if ((plist.rminbox[0] <= x) && (x <= plist.rmaxbox[0])) {
    if ((plist.rminbox[1] <= y) && (y <= plist.rmaxbox[1])) {
      if ((plist.rminbox[2] <= z) && (z <= plist.rmaxbox[2])) {
         val = myfft->shift;

        for (int i = 0; i < numcomp; i++) {
           val += myfft->cc[i] * sin( 2.0*M_PI*(myfft->kx[i]*(x-plist.rminbox[0]) + myfft->ky[i]*(y-plist.rminbox[1]) ) + myfft->phi[i]);

        }
      }
    }
  }
  else {
    val = -1;
  }
  return val;
}

double rand_wires(double x, double y, double z, PLASMAparams plist, double Z, double A){
    ALLWIRS* wirs = (ALLWIRS*)plist.additional_params;

    double val = -1;

    if ((plist.rminbox[0] <= x) && (x <= plist.rmaxbox[0])) {
      if ((plist.rminbox[1] <= y) && (y <= plist.rmaxbox[1])) {
        if ((plist.rminbox[2] <= z) && (z <= plist.rmaxbox[2])) {

            for(int i = 0; i < wirs->num; i++){

                double dist;

                double vx = wirs->x2[i] -  wirs->x1[i];
                double vy = wirs->y2[i] -  wirs->y1[i];
                double vz = wirs->z2[i] -  wirs->z1[i];

                double wx = x - wirs->x1[i];
                double wy = y - wirs->y1[i];
                double wz = z - wirs->z1[i];

                double c1 = wx*vx + wy*vy + wz*vz;
                if (c1 <= 0)
                    dist = sqrt(wx*wx + wy*wy + wz*wz);
                else{
                    double c2 = vx*vx + vy*vy + vz*vz;
                    if(c2 <= c1){
                        dist = sqrt((x - wirs->x2[i])*(x - wirs->x2[i]) + (y - wirs->y2[i])*(y - wirs->y2[i]) + (z - wirs->z2[i])*(z - wirs->z2[i]) );
                    }
                    else{
                        double b = c1/c2;
                        double bx = wirs->x1[i] + b*vx;
                        double by = wirs->y1[i] + b*vy;
                        double bz = wirs->z1[i] + b*vz;

                        dist = sqrt((x - bx)*(x - bx) + (y - by)*(y - by) + (z - bz)*(z - bz) );
                    }
                }

                if(dist <= wirs->radius){
                    val = plist.density_coefficient;
                    break;
                }
            }

          }
        }
      }
    else {
      val = -1;
    }
    return val;


}

void PLASMA::trimWirs(double llimits[3], double rlimits[3]){
    if(!params.amIRNDwire)
        return;

    ALLWIRS* wirs = (ALLWIRS*)params.additional_params;

    if(!wirs->init)
        return;

    int num = wirs->num;
    double radius = wirs->radius;

    double* x1 = wirs->x1;
    double* x2 = wirs->x2;
    double* y1 = wirs->y1;
    double* y2 = wirs->y2;
    double* z1 = wirs->z1;
    double* z2 = wirs->z2;

    //SELECT ONLY RELEVANT WIRES

    int del = 0;

    for(int i = 0; i < (num-del); i++){

        const double incfact = 1.01;

        double xleft = llimits[0] - radius*incfact;
        double xright =rlimits[0] + radius*incfact;

        double yleft = llimits[1] - radius*incfact;
        double yright =rlimits[1] + radius*incfact;

        double zleft = llimits[2] - radius*incfact;
        double zright =rlimits[2] + radius*incfact;

        bool isIn = false;


        double xa = x1[i];
        double xb = x2[i];
        double ya = y1[i];
        double yb = y2[i];
        double za = z1[i];
        double zb = z2[i];

        if( INRANGE(xa,xleft,xright)&&INRANGE(ya,yleft,yright)&&INRANGE(za,zleft,zright)&&INRANGE(xb,xleft,xright)&&INRANGE(yb,yleft,yright)&&INRANGE(zb,zleft,zright)){
                    isIn = true;
                }
                else{
                    VECT its(0,0,0);

                    if(SEG_PLANE_INTERSECT(VECT(xleft,0,0),VECT(1,0,0),VECT(xa,ya,za),VECT(xb,yb,zb),its)){
                        if(INRANGE(its.y,yleft,yright)&&INRANGE(its.z,zleft,zright))
                            isIn = true;
                    }

                    if(!isIn && SEG_PLANE_INTERSECT(VECT(xright,0,0),VECT(-1,0,0),VECT(xa,ya,za),VECT(xb,yb,zb),its)){
                        if(INRANGE(its.y,yleft,yright)&&INRANGE(its.z,zleft,zright))
                            isIn = true;
                    }

                    if(!isIn && SEG_PLANE_INTERSECT(VECT(0,yleft,0),VECT(0,1,0),VECT(xa,ya,za),VECT(xb,yb,zb),its)){
                        if(INRANGE(its.x,xleft,xright)&&INRANGE(its.z,zleft,zright))
                            isIn = true;
                    }

                    if(!isIn && SEG_PLANE_INTERSECT(VECT(0,yright,0),VECT(0,-1,0),VECT(xa,ya,za),VECT(xb,yb,zb),its)){
                        if(INRANGE(its.x,xleft,xright)&&INRANGE(its.z,zleft,zright))
                            isIn = true;
                    }

                    if(!isIn && SEG_PLANE_INTERSECT(VECT(0,0,zleft),VECT(0,0,1),VECT(xa,ya,za),VECT(xb,yb,zb),its)){
                        if(INRANGE(its.x,xleft,xright)&&INRANGE(its.y,yleft,yright))
                            isIn = true;
                    }

                    if(!isIn && SEG_PLANE_INTERSECT(VECT(0,0,zright),VECT(0,0,-1),VECT(xa,ya,za),VECT(xb,yb,zb),its)){
                        if(INRANGE(its.x,xleft,xright)&&INRANGE(its.y,yleft,yright))
                            isIn = true;
                    }

                }



        if(!isIn){
            double x1t = x1[num-del-1];
            double x2t = x2[num-del-1];
            double y1t = y1[num-del-1];
            double y2t = y2[num-del-1];
            double z1t = z1[num-del-1];
            double z2t = z2[num-del-1];

            x1[num-del-1] = x1[i];
            x2[num-del-1] = x2[i];
            y1[num-del-1] = y1[i];
            y2[num-del-1] = y2[i];
            z1[num-del-1] = z1[i];
            z2[num-del-1] = z2[i];

            x1[i] = x1t;
            x2[i] = x2t;
            y1[i] = y1t;
            y2[i] = y2t;
            z1[i] = z1t;
            z2[i] = z2t;

            i--;
            del++;
        }


    }

    double* x1r = new double[num-del];
    double* x2r = new double[num-del];
    double* y1r = new double[num-del];
    double* y2r = new double[num-del];
    double* z1r = new double[num-del];
    double* z2r = new double[num-del];

    std::memcpy(x1r, x1, sizeof(double)*(num-del));
    std::memcpy(x2r, x2, sizeof(double)*(num-del));
    std::memcpy(y1r, y1, sizeof(double)*(num-del));
    std::memcpy(y2r, y2, sizeof(double)*(num-del));
    std::memcpy(z1r, z1, sizeof(double)*(num-del));
    std::memcpy(z2r, z2, sizeof(double)*(num-del));

    delete[] x1;
    delete[] x2;
    delete[] y1;
    delete[] y2;
    delete[] z1;
    delete[] z2;

    wirs->x1 = x1r;
    wirs->x2 = x2r;
    wirs->y1 = y1r;
    wirs->y2 = y2r;
    wirs->z1 = z1r;
    wirs->z2 = z2r;

    wirs->num = num-del;

    //***SELECT ONLY RELEVANT WIRES***

}


double pillars3D (double x, double y, double z, PLASMAparams plist, double Z, double A){
    double *paramlist = (double*)plist.additional_params;
    double dy = paramlist[0];
    double dz = paramlist[1];
    double r = paramlist[2];
    double h = paramlist[3];
    double rhoDefault = -1;

    double middleY = 0.5*(plist.rminbox[1]+plist.rmaxbox[1]);
    double middleZ = 0.5*(plist.rminbox[2]+plist.rmaxbox[2]);

    if ( (plist.rminbox[0] <= x) && (x <= plist.rmaxbox[0]) &&
         (plist.rminbox[1] <= y) && (y <= plist.rmaxbox[1]) &&
         (plist.rminbox[2] <= z) && (z <= plist.rmaxbox[2]) ){

        if((plist.rmaxbox[0]-x) > h)
            return rhoDefault;

        double yy = y - middleY;
        double zz = z - middleZ;

        yy = yy - (dy)*round(yy/(dz));
        zz = zz - (dz)*round(zz/(dy));

        double rad2 = (yy - 0.5*dy)*(yy - 0.5*dy) + (zz - 0.5*dz)*(zz - 0.5*dz);

        if(rad2 <= r)
            return plist.density_coefficient;
        else
            return rhoDefault;

    }
    return rhoDefault;
}


//*************************END_PLASMA*****************************
//*************************LASER_PULSE***************************

laserPulse::laserPulse() {
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
  LG_l = 0;
  LG_m = 0;
}

laserPulse::~laserPulse() {}

laserPulse::laserPulse(const laserPulse& other)
  :type(other.type), t_FWHM(other.t_FWHM), waist(other.waist), focus_position(other.focus_position),
  laser_pulse_initial_position(other.laser_pulse_initial_position), normalized_amplitude(other.normalized_amplitude),
  lambda0(other.lambda0), rotation(other.rotation), angle(other.angle),
  rotation_center_along_x(other.rotation_center_along_x), rise_time(other.rise_time), LG_l(other.LG_l), LG_m(other.LG_m) {}


laserPulse laserPulse::operator=(const laserPulse& p1) {
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
  LG_l = p1.LG_l;
  LG_m = p1.LG_m;
  return *this;
}

void laserPulse::setFocusPosition(double _focus_position) {
  focus_position = _focus_position;
}

void laserPulse::setPulseInitialPosition(double _laser_pulse_initial_position) {
  laser_pulse_initial_position = _laser_pulse_initial_position;
}

void laserPulse::setLambda(double _lambda0) {
  lambda0 = _lambda0;
}
void laserPulse::setWaist(double _waist) {
  waist = _waist;
}

void laserPulse::setDurationFWHM(double _t_FWHM) {
  t_FWHM = _t_FWHM;
}

void laserPulse::setNormalizedAmplitude(double _normalized_amplitude) {
  normalized_amplitude = _normalized_amplitude;
}

void laserPulse::setRiseTime(double _rise_time) {
  rise_time = _rise_time;
}

void laserPulse::setRotationAngleAndCenter(double _angle, double _rotation_center_along_x) {
  angle = _angle;
  rotation_center_along_x = _rotation_center_along_x;
  rotation = true;
}

void laserPulse::setGaussianPulse(double _waist, double _t_FWHM, double _normalized_amplitude) {
  type = GAUSSIAN;
  waist = _waist;
  t_FWHM = _t_FWHM;
  normalized_amplitude = _normalized_amplitude;
}

void laserPulse::setPlaneWave(double _normalized_amplitude) {
  type = PLANE_WAVE;
  normalized_amplitude = _normalized_amplitude;
}

void laserPulse::setCos2PlaneWave(double _t_FWHM, double _normalized_amplitude) {
  type = COS2_PLANE_WAVE;
  t_FWHM = _t_FWHM;
  normalized_amplitude = _normalized_amplitude;
}

void laserPulse::setCos2PlateauPlaneWave(double _t_FWHM, double _rise_time, double _normalized_amplitude) {
  type = COS2_PLATEAU_PLANE_WAVE;
  t_FWHM = _t_FWHM;
  normalized_amplitude = _normalized_amplitude;
  rise_time = _rise_time;
}

void laserPulse::setGaussianPulse() {
  type = GAUSSIAN;
}

void laserPulse::setPlaneWave() {
  type = PLANE_WAVE;
}

void laserPulse::setCos2PlaneWave() {
  type = COS2_PLANE_WAVE;
}

void laserPulse::setCos2PlateauPlaneWave() {
  type = COS2_PLATEAU_PLANE_WAVE;
}

void laserPulse::setPPolarization() {
  polarization = P_POLARIZATION;
}

void laserPulse::setSPolarization() {
  polarization = S_POLARIZATION;
}

void laserPulse::setCircularPolarization() {
  polarization = CIRCULAR_POLARIZATION;
}

void laserPulse::setLaguerreGaussian_l(int l){
  LG_l = l;
}

void laserPulse::setLaguerreGaussian_m(int m){
  LG_m = m;
}

void laserPulse::setConstFieldComponent(int comp){
  component = comp;
}

//************** DISTRIBUTION_FUNCTION **********

tempDistrib::tempDistrib() {
  init = false;
  temp = 0;
}
tempDistrib tempDistrib::operator = (tempDistrib &destro){
  type = destro.type;
  p0 = destro.p0;
  p0_x = destro.p0_x;
  p0_y = destro.p0_y;
  p0_z = destro.p0_z;
  alpha = destro.alpha;
  temp = destro.temp;
  a = destro.a;
  init = destro.init;
  return *this;
}

bool tempDistrib::isInit() {
  return init;
}
 double tempDistrib::getTemperature(){
   return temp;
 }

void tempDistrib::setWaterbag(double _p0) {
  type = WATERBAG;
  p0 = _p0;
  init = true;
}

void tempDistrib::setWaterbag3Temp(double _p0_x, double _p0_y, double _p0_z) {
  type = WATERBAG_3TEMP;
  p0_x = _p0_x;
  p0_y = _p0_y;
  p0_z = _p0_z;
  init = true;
}

void tempDistrib::setUnifSphere(double _p0) {
  type = UNIF_SPHERE;
  p0 = _p0;
  init = true;
}

void tempDistrib::setSupergaussian(double _p0, double _alpha) {
  type = SUPERGAUSSIAN;
  p0 = _p0;
  alpha = _alpha;
  init = true;
}

void tempDistrib::setMaxwell(double _temp) {
  type = MAXWELL;
  temp = _temp;
  init = true;
}

void tempDistrib::setJuttner(double _a) {
  type = JUTTNER;
  a = _a;
  init = true;
}

void tempDistrib::setSpecial(double _a) {
  type = SPECIAL;
  a = _a;
  init = true;
}




//************** END DISTRIBUTION_FUNCTION ******

MUTILS::VECT::VECT(double x, double y, double z){
        this->x = x;
        this->y = y;
        this->z = z;
}
VECT MUTILS::VECT::operator+(const VECT& rhs){
          return VECT(x+rhs.x, y+rhs.y,z+rhs.z);
}
VECT MUTILS::VECT::operator-(const VECT& rhs){
        return VECT(x-rhs.x, y-rhs.y,z-rhs.z);
    }
double MUTILS::VECT::operator*(const VECT& rhs){
        return x*rhs.x + y*rhs.y + z*rhs.z;
    }
VECT MUTILS::VECT::operator*(const double& coeff){
        return VECT(coeff*x, coeff*y,coeff*z);
    }



bool MUTILS::SEG_PLANE_INTERSECT(VECT v0, VECT n, VECT xa, VECT xb, VECT& inters){
    double den = (n*(xb-xa));
    if(den == 0)
        return false;
    double s = (n*(v0 - xa))/den;
    if (INRANGE(s,0.0,1.0)){
        inters = xa + (xb-xa)*s;
        return true;
    }
    return false;
}


