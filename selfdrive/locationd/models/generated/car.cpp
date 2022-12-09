#include "car.h"

namespace {
#define DIM 9
#define EDIM 9
#define MEDIM 9
typedef void (*Hfun)(double *, double *, double *);

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}
const static double MAHA_THRESH_25 = 3.8414588206941227;
const static double MAHA_THRESH_24 = 5.991464547107981;
const static double MAHA_THRESH_30 = 3.8414588206941227;
const static double MAHA_THRESH_26 = 3.8414588206941227;
const static double MAHA_THRESH_27 = 3.8414588206941227;
const static double MAHA_THRESH_29 = 3.8414588206941227;
const static double MAHA_THRESH_28 = 3.8414588206941227;
const static double MAHA_THRESH_31 = 3.8414588206941227;

/******************************************************************************
 *                       Code generated with sympy 1.9                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_4330747659891609478) {
   out_4330747659891609478[0] = delta_x[0] + nom_x[0];
   out_4330747659891609478[1] = delta_x[1] + nom_x[1];
   out_4330747659891609478[2] = delta_x[2] + nom_x[2];
   out_4330747659891609478[3] = delta_x[3] + nom_x[3];
   out_4330747659891609478[4] = delta_x[4] + nom_x[4];
   out_4330747659891609478[5] = delta_x[5] + nom_x[5];
   out_4330747659891609478[6] = delta_x[6] + nom_x[6];
   out_4330747659891609478[7] = delta_x[7] + nom_x[7];
   out_4330747659891609478[8] = delta_x[8] + nom_x[8];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_3461638426050684139) {
   out_3461638426050684139[0] = -nom_x[0] + true_x[0];
   out_3461638426050684139[1] = -nom_x[1] + true_x[1];
   out_3461638426050684139[2] = -nom_x[2] + true_x[2];
   out_3461638426050684139[3] = -nom_x[3] + true_x[3];
   out_3461638426050684139[4] = -nom_x[4] + true_x[4];
   out_3461638426050684139[5] = -nom_x[5] + true_x[5];
   out_3461638426050684139[6] = -nom_x[6] + true_x[6];
   out_3461638426050684139[7] = -nom_x[7] + true_x[7];
   out_3461638426050684139[8] = -nom_x[8] + true_x[8];
}
void H_mod_fun(double *state, double *out_560952803384681099) {
   out_560952803384681099[0] = 1.0;
   out_560952803384681099[1] = 0;
   out_560952803384681099[2] = 0;
   out_560952803384681099[3] = 0;
   out_560952803384681099[4] = 0;
   out_560952803384681099[5] = 0;
   out_560952803384681099[6] = 0;
   out_560952803384681099[7] = 0;
   out_560952803384681099[8] = 0;
   out_560952803384681099[9] = 0;
   out_560952803384681099[10] = 1.0;
   out_560952803384681099[11] = 0;
   out_560952803384681099[12] = 0;
   out_560952803384681099[13] = 0;
   out_560952803384681099[14] = 0;
   out_560952803384681099[15] = 0;
   out_560952803384681099[16] = 0;
   out_560952803384681099[17] = 0;
   out_560952803384681099[18] = 0;
   out_560952803384681099[19] = 0;
   out_560952803384681099[20] = 1.0;
   out_560952803384681099[21] = 0;
   out_560952803384681099[22] = 0;
   out_560952803384681099[23] = 0;
   out_560952803384681099[24] = 0;
   out_560952803384681099[25] = 0;
   out_560952803384681099[26] = 0;
   out_560952803384681099[27] = 0;
   out_560952803384681099[28] = 0;
   out_560952803384681099[29] = 0;
   out_560952803384681099[30] = 1.0;
   out_560952803384681099[31] = 0;
   out_560952803384681099[32] = 0;
   out_560952803384681099[33] = 0;
   out_560952803384681099[34] = 0;
   out_560952803384681099[35] = 0;
   out_560952803384681099[36] = 0;
   out_560952803384681099[37] = 0;
   out_560952803384681099[38] = 0;
   out_560952803384681099[39] = 0;
   out_560952803384681099[40] = 1.0;
   out_560952803384681099[41] = 0;
   out_560952803384681099[42] = 0;
   out_560952803384681099[43] = 0;
   out_560952803384681099[44] = 0;
   out_560952803384681099[45] = 0;
   out_560952803384681099[46] = 0;
   out_560952803384681099[47] = 0;
   out_560952803384681099[48] = 0;
   out_560952803384681099[49] = 0;
   out_560952803384681099[50] = 1.0;
   out_560952803384681099[51] = 0;
   out_560952803384681099[52] = 0;
   out_560952803384681099[53] = 0;
   out_560952803384681099[54] = 0;
   out_560952803384681099[55] = 0;
   out_560952803384681099[56] = 0;
   out_560952803384681099[57] = 0;
   out_560952803384681099[58] = 0;
   out_560952803384681099[59] = 0;
   out_560952803384681099[60] = 1.0;
   out_560952803384681099[61] = 0;
   out_560952803384681099[62] = 0;
   out_560952803384681099[63] = 0;
   out_560952803384681099[64] = 0;
   out_560952803384681099[65] = 0;
   out_560952803384681099[66] = 0;
   out_560952803384681099[67] = 0;
   out_560952803384681099[68] = 0;
   out_560952803384681099[69] = 0;
   out_560952803384681099[70] = 1.0;
   out_560952803384681099[71] = 0;
   out_560952803384681099[72] = 0;
   out_560952803384681099[73] = 0;
   out_560952803384681099[74] = 0;
   out_560952803384681099[75] = 0;
   out_560952803384681099[76] = 0;
   out_560952803384681099[77] = 0;
   out_560952803384681099[78] = 0;
   out_560952803384681099[79] = 0;
   out_560952803384681099[80] = 1.0;
}
void f_fun(double *state, double dt, double *out_4520015740428216527) {
   out_4520015740428216527[0] = state[0];
   out_4520015740428216527[1] = state[1];
   out_4520015740428216527[2] = state[2];
   out_4520015740428216527[3] = state[3];
   out_4520015740428216527[4] = state[4];
   out_4520015740428216527[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] - 9.8000000000000007*state[8] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_4520015740428216527[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_4520015740428216527[7] = state[7];
   out_4520015740428216527[8] = state[8];
}
void F_fun(double *state, double dt, double *out_3495565834374361759) {
   out_3495565834374361759[0] = 1;
   out_3495565834374361759[1] = 0;
   out_3495565834374361759[2] = 0;
   out_3495565834374361759[3] = 0;
   out_3495565834374361759[4] = 0;
   out_3495565834374361759[5] = 0;
   out_3495565834374361759[6] = 0;
   out_3495565834374361759[7] = 0;
   out_3495565834374361759[8] = 0;
   out_3495565834374361759[9] = 0;
   out_3495565834374361759[10] = 1;
   out_3495565834374361759[11] = 0;
   out_3495565834374361759[12] = 0;
   out_3495565834374361759[13] = 0;
   out_3495565834374361759[14] = 0;
   out_3495565834374361759[15] = 0;
   out_3495565834374361759[16] = 0;
   out_3495565834374361759[17] = 0;
   out_3495565834374361759[18] = 0;
   out_3495565834374361759[19] = 0;
   out_3495565834374361759[20] = 1;
   out_3495565834374361759[21] = 0;
   out_3495565834374361759[22] = 0;
   out_3495565834374361759[23] = 0;
   out_3495565834374361759[24] = 0;
   out_3495565834374361759[25] = 0;
   out_3495565834374361759[26] = 0;
   out_3495565834374361759[27] = 0;
   out_3495565834374361759[28] = 0;
   out_3495565834374361759[29] = 0;
   out_3495565834374361759[30] = 1;
   out_3495565834374361759[31] = 0;
   out_3495565834374361759[32] = 0;
   out_3495565834374361759[33] = 0;
   out_3495565834374361759[34] = 0;
   out_3495565834374361759[35] = 0;
   out_3495565834374361759[36] = 0;
   out_3495565834374361759[37] = 0;
   out_3495565834374361759[38] = 0;
   out_3495565834374361759[39] = 0;
   out_3495565834374361759[40] = 1;
   out_3495565834374361759[41] = 0;
   out_3495565834374361759[42] = 0;
   out_3495565834374361759[43] = 0;
   out_3495565834374361759[44] = 0;
   out_3495565834374361759[45] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_3495565834374361759[46] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_3495565834374361759[47] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_3495565834374361759[48] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_3495565834374361759[49] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_3495565834374361759[50] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_3495565834374361759[51] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_3495565834374361759[52] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_3495565834374361759[53] = -9.8000000000000007*dt;
   out_3495565834374361759[54] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_3495565834374361759[55] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_3495565834374361759[56] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_3495565834374361759[57] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_3495565834374361759[58] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_3495565834374361759[59] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_3495565834374361759[60] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_3495565834374361759[61] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_3495565834374361759[62] = 0;
   out_3495565834374361759[63] = 0;
   out_3495565834374361759[64] = 0;
   out_3495565834374361759[65] = 0;
   out_3495565834374361759[66] = 0;
   out_3495565834374361759[67] = 0;
   out_3495565834374361759[68] = 0;
   out_3495565834374361759[69] = 0;
   out_3495565834374361759[70] = 1;
   out_3495565834374361759[71] = 0;
   out_3495565834374361759[72] = 0;
   out_3495565834374361759[73] = 0;
   out_3495565834374361759[74] = 0;
   out_3495565834374361759[75] = 0;
   out_3495565834374361759[76] = 0;
   out_3495565834374361759[77] = 0;
   out_3495565834374361759[78] = 0;
   out_3495565834374361759[79] = 0;
   out_3495565834374361759[80] = 1;
}
void h_25(double *state, double *unused, double *out_2384350295709118751) {
   out_2384350295709118751[0] = state[6];
}
void H_25(double *state, double *unused, double *out_1152488062000513114) {
   out_1152488062000513114[0] = 0;
   out_1152488062000513114[1] = 0;
   out_1152488062000513114[2] = 0;
   out_1152488062000513114[3] = 0;
   out_1152488062000513114[4] = 0;
   out_1152488062000513114[5] = 0;
   out_1152488062000513114[6] = 1;
   out_1152488062000513114[7] = 0;
   out_1152488062000513114[8] = 0;
}
void h_24(double *state, double *unused, double *out_3173356781225767583) {
   out_3173356781225767583[0] = state[4];
   out_3173356781225767583[1] = state[5];
}
void H_24(double *state, double *unused, double *out_1020161537004986452) {
   out_1020161537004986452[0] = 0;
   out_1020161537004986452[1] = 0;
   out_1020161537004986452[2] = 0;
   out_1020161537004986452[3] = 0;
   out_1020161537004986452[4] = 1;
   out_1020161537004986452[5] = 0;
   out_1020161537004986452[6] = 0;
   out_1020161537004986452[7] = 0;
   out_1020161537004986452[8] = 0;
   out_1020161537004986452[9] = 0;
   out_1020161537004986452[10] = 0;
   out_1020161537004986452[11] = 0;
   out_1020161537004986452[12] = 0;
   out_1020161537004986452[13] = 0;
   out_1020161537004986452[14] = 1;
   out_1020161537004986452[15] = 0;
   out_1020161537004986452[16] = 0;
   out_1020161537004986452[17] = 0;
}
void h_30(double *state, double *unused, double *out_3786314753859575719) {
   out_3786314753859575719[0] = state[4];
}
void H_30(double *state, double *unused, double *out_3670821020507761741) {
   out_3670821020507761741[0] = 0;
   out_3670821020507761741[1] = 0;
   out_3670821020507761741[2] = 0;
   out_3670821020507761741[3] = 0;
   out_3670821020507761741[4] = 1;
   out_3670821020507761741[5] = 0;
   out_3670821020507761741[6] = 0;
   out_3670821020507761741[7] = 0;
   out_3670821020507761741[8] = 0;
}
void h_26(double *state, double *unused, double *out_8787988021520421364) {
   out_8787988021520421364[0] = state[7];
}
void H_26(double *state, double *unused, double *out_2589015256873543110) {
   out_2589015256873543110[0] = 0;
   out_2589015256873543110[1] = 0;
   out_2589015256873543110[2] = 0;
   out_2589015256873543110[3] = 0;
   out_2589015256873543110[4] = 0;
   out_2589015256873543110[5] = 0;
   out_2589015256873543110[6] = 0;
   out_2589015256873543110[7] = 1;
   out_2589015256873543110[8] = 0;
}
void h_27(double *state, double *unused, double *out_5005018669045744043) {
   out_5005018669045744043[0] = state[3];
}
void H_27(double *state, double *unused, double *out_1496057708707336830) {
   out_1496057708707336830[0] = 0;
   out_1496057708707336830[1] = 0;
   out_1496057708707336830[2] = 0;
   out_1496057708707336830[3] = 1;
   out_1496057708707336830[4] = 0;
   out_1496057708707336830[5] = 0;
   out_1496057708707336830[6] = 0;
   out_1496057708707336830[7] = 0;
   out_1496057708707336830[8] = 0;
}
void h_29(double *state, double *unused, double *out_6254772500205749244) {
   out_6254772500205749244[0] = state[1];
}
void H_29(double *state, double *unused, double *out_4181052364822153925) {
   out_4181052364822153925[0] = 0;
   out_4181052364822153925[1] = 1;
   out_4181052364822153925[2] = 0;
   out_4181052364822153925[3] = 0;
   out_4181052364822153925[4] = 0;
   out_4181052364822153925[5] = 0;
   out_4181052364822153925[6] = 0;
   out_4181052364822153925[7] = 0;
   out_4181052364822153925[8] = 0;
}
void h_28(double *state, double *unused, double *out_7042683887249958427) {
   out_7042683887249958427[0] = state[0];
}
void H_28(double *state, double *unused, double *out_901346652247376649) {
   out_901346652247376649[0] = 1;
   out_901346652247376649[1] = 0;
   out_901346652247376649[2] = 0;
   out_901346652247376649[3] = 0;
   out_901346652247376649[4] = 0;
   out_901346652247376649[5] = 0;
   out_901346652247376649[6] = 0;
   out_901346652247376649[7] = 0;
   out_901346652247376649[8] = 0;
}
void h_31(double *state, double *unused, double *out_7698958344398103958) {
   out_7698958344398103958[0] = state[8];
}
void H_31(double *state, double *unused, double *out_3215223359106894586) {
   out_3215223359106894586[0] = 0;
   out_3215223359106894586[1] = 0;
   out_3215223359106894586[2] = 0;
   out_3215223359106894586[3] = 0;
   out_3215223359106894586[4] = 0;
   out_3215223359106894586[5] = 0;
   out_3215223359106894586[6] = 0;
   out_3215223359106894586[7] = 0;
   out_3215223359106894586[8] = 1;
}
#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;

  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);

  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H);

  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();


    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;

  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);

  // update cov
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}




}
extern "C" {

void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
}
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<2, 3, 0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
}
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
}
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
}
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
}
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
}
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
}
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_31, H_31, NULL, in_z, in_R, in_ea, MAHA_THRESH_31);
}
void car_err_fun(double *nom_x, double *delta_x, double *out_4330747659891609478) {
  err_fun(nom_x, delta_x, out_4330747659891609478);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_3461638426050684139) {
  inv_err_fun(nom_x, true_x, out_3461638426050684139);
}
void car_H_mod_fun(double *state, double *out_560952803384681099) {
  H_mod_fun(state, out_560952803384681099);
}
void car_f_fun(double *state, double dt, double *out_4520015740428216527) {
  f_fun(state,  dt, out_4520015740428216527);
}
void car_F_fun(double *state, double dt, double *out_3495565834374361759) {
  F_fun(state,  dt, out_3495565834374361759);
}
void car_h_25(double *state, double *unused, double *out_2384350295709118751) {
  h_25(state, unused, out_2384350295709118751);
}
void car_H_25(double *state, double *unused, double *out_1152488062000513114) {
  H_25(state, unused, out_1152488062000513114);
}
void car_h_24(double *state, double *unused, double *out_3173356781225767583) {
  h_24(state, unused, out_3173356781225767583);
}
void car_H_24(double *state, double *unused, double *out_1020161537004986452) {
  H_24(state, unused, out_1020161537004986452);
}
void car_h_30(double *state, double *unused, double *out_3786314753859575719) {
  h_30(state, unused, out_3786314753859575719);
}
void car_H_30(double *state, double *unused, double *out_3670821020507761741) {
  H_30(state, unused, out_3670821020507761741);
}
void car_h_26(double *state, double *unused, double *out_8787988021520421364) {
  h_26(state, unused, out_8787988021520421364);
}
void car_H_26(double *state, double *unused, double *out_2589015256873543110) {
  H_26(state, unused, out_2589015256873543110);
}
void car_h_27(double *state, double *unused, double *out_5005018669045744043) {
  h_27(state, unused, out_5005018669045744043);
}
void car_H_27(double *state, double *unused, double *out_1496057708707336830) {
  H_27(state, unused, out_1496057708707336830);
}
void car_h_29(double *state, double *unused, double *out_6254772500205749244) {
  h_29(state, unused, out_6254772500205749244);
}
void car_H_29(double *state, double *unused, double *out_4181052364822153925) {
  H_29(state, unused, out_4181052364822153925);
}
void car_h_28(double *state, double *unused, double *out_7042683887249958427) {
  h_28(state, unused, out_7042683887249958427);
}
void car_H_28(double *state, double *unused, double *out_901346652247376649) {
  H_28(state, unused, out_901346652247376649);
}
void car_h_31(double *state, double *unused, double *out_7698958344398103958) {
  h_31(state, unused, out_7698958344398103958);
}
void car_H_31(double *state, double *unused, double *out_3215223359106894586) {
  H_31(state, unused, out_3215223359106894586);
}
void car_predict(double *in_x, double *in_P, double *in_Q, double dt) {
  predict(in_x, in_P, in_Q, dt);
}
void car_set_mass(double x) {
  set_mass(x);
}
void car_set_rotational_inertia(double x) {
  set_rotational_inertia(x);
}
void car_set_center_to_front(double x) {
  set_center_to_front(x);
}
void car_set_center_to_rear(double x) {
  set_center_to_rear(x);
}
void car_set_stiffness_front(double x) {
  set_stiffness_front(x);
}
void car_set_stiffness_rear(double x) {
  set_stiffness_rear(x);
}
}

const EKF car = {
  .name = "car",
  .kinds = { 25, 24, 30, 26, 27, 29, 28, 31 },
  .feature_kinds = {  },
  .f_fun = car_f_fun,
  .F_fun = car_F_fun,
  .err_fun = car_err_fun,
  .inv_err_fun = car_inv_err_fun,
  .H_mod_fun = car_H_mod_fun,
  .predict = car_predict,
  .hs = {
    { 25, car_h_25 },
    { 24, car_h_24 },
    { 30, car_h_30 },
    { 26, car_h_26 },
    { 27, car_h_27 },
    { 29, car_h_29 },
    { 28, car_h_28 },
    { 31, car_h_31 },
  },
  .Hs = {
    { 25, car_H_25 },
    { 24, car_H_24 },
    { 30, car_H_30 },
    { 26, car_H_26 },
    { 27, car_H_27 },
    { 29, car_H_29 },
    { 28, car_H_28 },
    { 31, car_H_31 },
  },
  .updates = {
    { 25, car_update_25 },
    { 24, car_update_24 },
    { 30, car_update_30 },
    { 26, car_update_26 },
    { 27, car_update_27 },
    { 29, car_update_29 },
    { 28, car_update_28 },
    { 31, car_update_31 },
  },
  .Hes = {
  },
  .sets = {
    { "mass", car_set_mass },
    { "rotational_inertia", car_set_rotational_inertia },
    { "center_to_front", car_set_center_to_front },
    { "center_to_rear", car_set_center_to_rear },
    { "stiffness_front", car_set_stiffness_front },
    { "stiffness_rear", car_set_stiffness_rear },
  },
  .extra_routines = {
  },
};

ekf_init(car);
