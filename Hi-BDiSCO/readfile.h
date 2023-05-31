// **************************************************************************************//
//											 //
//		 		Brownian Dynamics Simulation Algorithm			 //
// 		     Copyright Zilong Li, Tamar Schlick and New York University 	 //
//					    April 2020					 //
//                                                                                       //
// **************************************************************************************//


#ifndef READFILE
#define READFILE

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

#include "constants.h"

using namespace std;

std::vector<double> Read_const(string filename);

void Read_initial_struc(int &n, int &n_c, int &n3, int &nc3, std::vector<double> &tmp_r, std::vector<double> &tmp_a, std::vector<double> &tmp_b, std::vector<double> &tmp_c, std::vector<double> &tmp_rad, std::vector<double> &tmp_Er, std::vector<int> &tmp_type, string filename);

void Read_core(int &Nq, int &Nq3, std::vector<double> &core_pos, std::vector<double> &core_q, string filename);

void Read_extra(int &ex_n, std::vector<int> &ex_m, string filename);

void Read_tail(int &n_t, std::vector<double> &tail_pos, std::vector<double> &tail_q, std::vector<int> &tail_grp, std::vector<int> &tail_fix, std::vector<double> &tail_rad, std::vector<double> &tail_bond_c, std::vector<double> &tail_bond_v, std::vector<double> &tail_angle_c, std::vector<double> &tail_angle_v, string filename);

void Read_LH(int &n_lh_g, int &n_lh_n, int &n_lh_c, std::vector<double> &LH_g_pos, std::vector<double> &LH_n_pos, std::vector<double> &LH_c_pos, std::vector<int> &LH_conn, std::vector<double> &LH_q, std::vector<double> &LH_vdw_hh, std::vector<double> &LH_vdw_hc, std::vector<double> &LH_vdw_hl, std::vector<double> &LH_vdw_ht, std::vector<double> &LH_kstr, std::vector<double> &LH_kben, std::vector<double> &LH_streq, std::vector<double> &LH_betaeq, std::vector<double> &LH_radius, string filename);

void write_xyz_append(int n, int n3, int n_c, int Nq, int Nq3, int* type, double* r, double* a, double* b, double* c, std::vector<double> core_pos, string filename);

void write_restart(int n3, double* r, double* a, double* b, double* c);

void Read_restart_ini(int n3, double* r, double* a, double* b, double* c);

void Read_restart_tail(int n_tail3, double* r_t);

void Read_restart_LH(int n_LH3, double* r_lh);

void Read_LH_core(int n_c, int &nc_lh, int* nc_lh_flag);

#endif
