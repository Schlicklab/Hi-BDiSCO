// **************************************************************************************//
//											 //
//		 		Brownian Dynamics Simulation Algorithm			 //
// 		     Copyright Zilong Li, Tamar Schlick and New York University 	 //
//					    April 2020					 //
//                                                                                       //
// **************************************************************************************//


#ifndef FUNC
#define FUNC


#include <vector>
#include <algorithm>
#include <functional>
#include <list>
#include <iterator>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <stdio.h>
#include <cmath>
#include <cctype>

#include "mt.h"
#include "utilities.h"
#include "constants.h"

using namespace std;

double twist_angle_phi(int w);

double Energy_Stretching(double h, std::vector<double> coord1, std::vector<double> coord2, double lo);

double Energy_Bending(double g, double beta, double beta_o);

double Energy_Twisting(double s, double alpha, double gamma, double phi_o, double lo);

void Force_Stretching(double h, double* coord1, double* coord2, double lo, double* force_projection1, double* force_projection2, double& Energy);

void Force_Bending(double g, double beta, double beta_o, double* coord1, double* coord2, double* coord3, double* force_projection1, double* force_projection2, double* force_projection3, double& Energy);

void Bending_force_projection(double g, double beta, double beta_b, double length, double* a_f, double* a_b, double* a, double* force_projection1, double* force_projection2, double& Energy);

void Twisting_force_projection(double s, double alpha, double beta, double gamma, double phi_o, double length, double alpha_b, double beta_b, double gamma_b, double phi_o_b, double gamma_n, double* b, double* c, double* force_projection1, double* force_projection2, double& Energy);

double Energy_Electrostatics(double q1, double q2, double epslon, double kappa, std::vector<double> coord1, std::vector<double> coord2);

void Force_Electrostatics(double q1, double q2, double epslon, double kappa, double* coord1, double* coord2, double* force_projection1, double* force_projection2i, double& Energy);

double Energy_Exclude_Volume(double k_ev, double sigma, std::vector<double> coord1, std::vector<double> coord2);

void Force_Exclude_Volume(double k_ev, double sigma, double* coord1, double* coord2, double* force_projection1, double* force_projection2i, double& Energy);

void Force_Ele_Vdw(double q1, double q2, double epslon, double kappa, double k_ev, double sigma, double* coord1, double* coord2, double* force_projection1, double* force_projection2, double& Energy);

void torque_due_to_force(double* force, double* coord_f, double* coord_c, double* a, double* b, double* c, double* torque);

void torque_due_to_force_relative(double* force, double* comp, double* a, double* b, double* c, double* torque);

void update_phi_o(int n, int* type, double* phi_o);

void update_Euler_Angle(int n_c, int nc3, int n, int n3, int* type, double* r, double* a, double* b, double* c, double* alpha, double* beta, double* gamma, double* length, double* a_dna, double* b_dna, double* c_dna, double* alpha_p, double* beta_p, double* gamma_p);

void first_coord(int t, double* r, double* a, double* b, double* c, double* r_f);

void second_coord(int t, double* r, double* a, double* b, double* c, double* r_s);

void Diffusion_Tensor(int n, int n3, double* r, double a1, double a2, double* rad, std::vector<std::vector<double>> &D, std::vector<std::vector<double>> &Chol);

int Diffusion_Tensor_check(int n, int n3, double* r, double a1, double a2, double* rad, std::vector<std::vector<double>> &D, std::vector<std::vector<double>> &Chol);

void mechanical_force_and_torque(int n_c, int nc3, int n, int n3, int* type, double* r, double* a, double* b, double* c, double* alpha, double* beta, double* gamma, double* length, double* a_dna, double* b_dna, double* c_dna, double* alpha_p, double* beta_p, double* gamma_p, double h, double g, double s, double* phi_o, double* force, double* torque, double& Energy, int rank, int comm_size);

void Electrostatic_and_Excluded_volume_force(int n, int n3, int n_c, int nc3, int* type, double* r, double* a, double* b, double* c, double debyell, double debye, double q_l, double k_e, double k_ex, double k_h1, double sigma_DNA_DNA, double sigma_DNA_Core, double sigma_Core_Core, int Nq, int Nq3, std::vector<double> core_pos, std::vector<double> core_q, double* force, double* torque, double& Energy, int rank, int comm_size);

void rotate(int n, int n3, double* a, double* b, double* c, double* a_n, double* b_n, double* c_n, double* d_theta, double dt);

void build_tail(int n, int* type, int n_t, std::vector<double> tail_pos, std::vector<double> tail_q, std::vector<int> tail_grp, std::vector<int> tail_fix, std::vector<double> tail_rad, std::vector<double> tail_bond_c, std::vector<double> tail_bond_v, std::vector<double> tail_angle_c, std::vector<double> tail_angle_v, double* r, double* a, double* b, double* c, int n_tail, int* nc_t_flag, double* r_t, double* h_t, double* g_t, double* lo_t, double* beta_o_t, double* t_q, double* t_rad, int* t_grp, int* t_fix);

void tail_local(int n, int* type, int n_t, double* r, double* a, double* b, double* c, int* nc_t_flag, int n_tail, double* r_t, double* r_t_local);

void tail_grow(int n, int* type, int n_t, double* r, double* a, double* b, double* c, int* nc_t_flag, int n_tail, double* r_t, double* r_t_local);

void update_tail_beta(int n_tail, double* r_t, double* beta_t, int* t_grp);

void tail_force(int n, int n_t, int* type, int n_tail, int n_tail3, std::vector<double> tail_pos, std::vector<int> tail_fix, int* nc_t_flag, double* r_t, double* beta_t, double* h_t, double* g_t, double* lo_t, double* beta_o_t, double h, double k_e, double debye, double k_ex, double* t_q, double* t_rad, int* t_grp, int* t_fix, double* r, double* a, double* b, double* c, double q_l, int Nq, int Nq3, std::vector<double> core_pos, std::vector<double> core_q, double* t_force, double* force, double* torque, double sigma_Tail_Tail, double sigma_Tail_Linker, double sigma_Tail_Core, double& Energy, int rank, int comm_size);

void build_LH(int n, int* type, int n_LH, std::vector<double> LH_n_pos, std::vector<double> LH_c_pos, std::vector<double> LH_radius, double* r, double* a, double* b, double* c, int n_lh_n, int n_lh_c, int* nc_lh_flag, double* r_lh, double* lh_rad);

void update_LH_beta(int n_LH, int n_lh_c, double* r_lh, double* beta_lh);

void Linker_Histone_force(int n, int n_tail, int n_lh_n, int n_lh_g, int n_lh_c, int* t_grp, int* t_fix, int* type, int n_LH, int n_LH3, std::vector<double> LH_g_pos, std::vector<int> LH_conn, int* nc_lh_flag, double* beta_lh, double* r_lh, std::vector<double> LH_q, double k_e, double debye, double k_ex, std::vector<double> LH_vdw_hh, std::vector<double> LH_vdw_hc, std::vector<double> LH_vdw_hl, std::vector<double> LH_vdw_ht, std::vector<double> LH_kstr, std::vector<double> LH_kben, std::vector<double> LH_streq, std::vector<double> LH_betaeq, double* r, double* a, double* b, double* c, double q_l, int Nq, int Nq3, std::vector<double> core_pos, std::vector<double> core_q, double* r_t, double* t_q, double* t_force, double* LH_force, double* force, double* torque, double& Energy, int rank, int comm_size);

#endif
