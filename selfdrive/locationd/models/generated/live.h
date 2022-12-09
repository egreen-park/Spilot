#pragma once
#include "rednose/helpers/common_ekf.h"
extern "C" {
void live_update_4(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_9(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_10(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_12(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_32(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_13(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_14(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_33(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_H(double *in_vec, double *out_1500569955103598890);
void live_err_fun(double *nom_x, double *delta_x, double *out_6370410215467224596);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_4404289754720028568);
void live_H_mod_fun(double *state, double *out_65844624761511818);
void live_f_fun(double *state, double dt, double *out_4639761141299267529);
void live_F_fun(double *state, double dt, double *out_2080679607502804155);
void live_h_4(double *state, double *unused, double *out_6128971225621953021);
void live_H_4(double *state, double *unused, double *out_1892393878140881810);
void live_h_9(double *state, double *unused, double *out_771804716078600861);
void live_H_9(double *state, double *unused, double *out_996467674139197532);
void live_h_10(double *state, double *unused, double *out_8834859403002773647);
void live_H_10(double *state, double *unused, double *out_6626926868353906761);
void live_h_12(double *state, double *unused, double *out_4194655737861473773);
void live_H_12(double *state, double *unused, double *out_1271294853093288143);
void live_h_31(double *state, double *unused, double *out_7667153106358443581);
void live_H_31(double *state, double *unused, double *out_1474268179231725566);
void live_h_32(double *state, double *unused, double *out_5319208611029596904);
void live_H_32(double *state, double *unused, double *out_7275572551439338124);
void live_h_13(double *state, double *unused, double *out_6041845049562932408);
void live_H_13(double *state, double *unused, double *out_595429171797323065);
void live_h_14(double *state, double *unused, double *out_771804716078600861);
void live_H_14(double *state, double *unused, double *out_996467674139197532);
void live_h_33(double *state, double *unused, double *out_5377665374765292145);
void live_H_33(double *state, double *unused, double *out_4624825183870583170);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}