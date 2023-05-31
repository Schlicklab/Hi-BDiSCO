// **************************************************************************************//
//                                                                                       //
//                              Brownian Dynamics Simulation Algorithm                   //
//                   Copyright Zilong Li, Tamar Schlick and New York University          //
//                                          April 2020                                   //
//                                                                                       //
// **************************************************************************************//


#ifndef FUNC_CUDA
#define FUNC_CUDA


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
#include <cublas_v2.h>
#include <cusolverDn.h>
#include <cuda_runtime_api.h> 

#include "mt.h"
#include "utilities.h"
#include "constants.h"

using namespace std;

extern "C++" void cuda_application_init_D_Chol(int n3);

extern "C++" void cuda_application_init_data(int n_c, int nc3, int n, int n3, int* type, double* r, double* a, double* b, double* c, double* alpha, double* beta, double* gamma, double* length, double* a_dna, double* b_dna, double* c_dna, double* alpha_p, double* beta_p, double* gamma_p, double h, double g, double s, double* phi_o, double debyell, double debye, double q_l, double k_e, double k_ex, double k_h1, double sigma_DNA_DNA, double sigma_DNA_Core, double sigma_Core_Core, double sigma_Tail_Tail, double sigma_Tail_Linker, double sigma_Tail_Core, int Nq, int Nq3, double* core_pos, double* core_q, double* force, double* torque, double* Energy, double* r_all, double* rad_all, int ex_n, int* ex_force_m);

extern "C++" void main_cuda(int n_c, int nc3, int ex_n, int step, int number_of_steps, double time_step, double del, int frequency_RP, int frequency_of_sampling, double h, double g, double s, double debyell, double debye, double q_l, double k_e, double k_ex, double k_h1, double sigma_DNA_DNA, double sigma_DNA_Core, double sigma_Core_Core, double sigma_Tail_Tail, double sigma_Tail_Linker, double sigma_Tail_Core, int Nq, int Nq3, int n, int n3, double a1, double a2, double s2dt, double* rr, double* p, double* Energy, double* h_r, double* h_a, double* h_b, double* h_c, double* h_rad_all);

extern "C++" void free_all();

#endif

