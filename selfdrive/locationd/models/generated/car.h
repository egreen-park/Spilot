#pragma once
#include "rednose/helpers/common_ekf.h"
extern "C" {
void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_err_fun(double *nom_x, double *delta_x, double *out_4330747659891609478);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_3461638426050684139);
void car_H_mod_fun(double *state, double *out_560952803384681099);
void car_f_fun(double *state, double dt, double *out_4520015740428216527);
void car_F_fun(double *state, double dt, double *out_3495565834374361759);
void car_h_25(double *state, double *unused, double *out_2384350295709118751);
void car_H_25(double *state, double *unused, double *out_1152488062000513114);
void car_h_24(double *state, double *unused, double *out_3173356781225767583);
void car_H_24(double *state, double *unused, double *out_1020161537004986452);
void car_h_30(double *state, double *unused, double *out_3786314753859575719);
void car_H_30(double *state, double *unused, double *out_3670821020507761741);
void car_h_26(double *state, double *unused, double *out_8787988021520421364);
void car_H_26(double *state, double *unused, double *out_2589015256873543110);
void car_h_27(double *state, double *unused, double *out_5005018669045744043);
void car_H_27(double *state, double *unused, double *out_1496057708707336830);
void car_h_29(double *state, double *unused, double *out_6254772500205749244);
void car_H_29(double *state, double *unused, double *out_4181052364822153925);
void car_h_28(double *state, double *unused, double *out_7042683887249958427);
void car_H_28(double *state, double *unused, double *out_901346652247376649);
void car_h_31(double *state, double *unused, double *out_7698958344398103958);
void car_H_31(double *state, double *unused, double *out_3215223359106894586);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}