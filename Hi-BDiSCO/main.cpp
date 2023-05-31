// **************************************************************************************//
//                                                                                       //
//                              Brownian Dynamics Simulation Algorithm                   //
//                   Copyright Zilong Li, Tamar Schlick and New York University          //
//                                          April 2020                                   //
//                                                                                       //
// **************************************************************************************//


#include "constants.h"
#include "readfile.h"
#include "func.h"
#include "func_cuda.h"
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <iomanip>
using namespace std;

int main(int argc, char *argv[]){

	std::vector<double> input_consts;

        input_consts = Read_const("setup.txt");

        double T = input_consts[0];
        int number_of_steps = (int)input_consts[1];
        int frequency_of_sampling = (int)input_consts[2];
        int frequency_RP = (int)input_consts[3];
        double sigma_DNA_DNA = input_consts[4];
        double sigma_DNA_Core = input_consts[5];
        double sigma_Core_Core = input_consts[6];
        double sigma_Tail_Tail = input_consts[7];
        double sigma_Tail_Linker = input_consts[8];
        double sigma_Tail_Core = input_consts[9];
        int Mg_flag = (int)input_consts[10];
        double Cs = input_consts[11];
        int restart = (int)input_consts[12];
        double unbind = input_consts[13];
        double Pbind = input_consts[14];

	//random seed

        time_t now = time(0);

        init_genrand((unsigned long) now);

	//Other constants/parameters


        double kb = 1.380649e-5; // (nm^2*kg)/(s^2*K)
        double kbt = kb*300.0; // (nm^2*kg)/s^2

        double q_l = lo * (-5.8824) *
                (1.901912136835033e-8 * pow((Cs * 1000), 3.0) + -8.211102025728157e-6 * pow((Cs * 1000), 2.0)
                + 7.554037628581672e-3 * (Cs * 1000) + 3.524292543853884e-1);
        double per;
        double debye = 0.736*sqrt((Cs/0.05)*(298/T));
        double debyell;

        if (Mg_flag==0){
                per = 50;
                debyell = debye;
        }else{
                per = 30;
                debyell = 2.5;
        }
        double k_e = 0.4151*kb*300;
        double k_ex = 0.001*kbt;
        double k_h1 = 10.0*kbt;
        double a1 = kbt/(6*PI*eta);
        double a2 = kbt/(8*PI*eta);
        double h = 100*kbt/(lo*lo);
        double hd2 = 0.5*h;
        double g = per*kbt/lo;
        double gd2 = 0.5*g;
        double s = 72.429*kbt/lo;
        double sd2 = 0.5*s;

        double time_step = 1e-12;
        double del = time_step/kbt;
        double s2dt = sqrt(2*time_step);


	//Initial structure

        int n, n_c, n3, nc3, ex_n;
        std::vector<double> tmp_r, tmp_a, tmp_b, tmp_c, tmp_rad, tmp_Er;
        std::vector<int> tmp_type, ex_m;

        Read_initial_struc(n, n_c, n3, nc3, tmp_r, tmp_a, tmp_b, tmp_c, tmp_rad, tmp_Er, tmp_type, "data_mod");

	Read_extra(ex_n, ex_m, "ex_force.txt");


	//Place holder

        double r[n3], a[n3], b[n3], c[n3], Er[n3];
        int type[n], ex_force_m[ex_n*2];
        double rad[n];
        double alpha[n], beta[n], gamma[n], phi_o[n];
        double a_dna[nc3], b_dna[nc3], c_dna[nc3];
        double alpha_p[n_c], beta_p[n_c], gamma_p[n_c];
        double length[n];

        double rr[n3], r_var[n3], d_theta[n3];
        double force[n3], torque[n3], force_test[n3], torque_test[n3];
        double force_global[n3], torque_global[n3];
        double E[1];
        double Energy;
        double Energy_global;

        //Place holder for intermediate configurations

        double r_n[n3], a_n[n3], b_n[n3], c_n[n3], r_n_test[n3];
        double alpha_n[n], beta_n[n], gamma_n[n];
        double a_dna_n[nc3], b_dna_n[nc3], c_dna_n[nc3];
        double alpha_p_n[n_c], beta_p_n[n_c], gamma_p_n[n_c];
        double length_n[n];
        double force_n[n3], torque_n[n3];
        double force_n_global[n3], torque_n_global[n3], force_global_tmp[n3];

	//Counting
        int i, j, k;
        int i1, i2, i3, j1, j2, j3, k1, k2, k3;

        double tmp_dummy;

        if (restart==0){
                for (i=0; i<n3; i++){
                        r[i] = tmp_r[i];
                        a[i] = tmp_a[i];
                        b[i] = tmp_b[i];
                        c[i] = tmp_c[i];
                        Er[i] = tmp_Er[i];
                }
        }else{
                Read_restart_ini(n3, r, a, b, c);
        }

        for (i=0; i<n; i++){
                rad[i] = tmp_rad[i];
                type[i] = tmp_type[i];
        }

	for (i=0;i<ex_n*2;i++){
		ex_force_m[i]=ex_m[i];
	}

        update_phi_o(n, type, phi_o);

	//Update director frame and Euler Angles

        update_Euler_Angle(n_c, nc3, n, n3, type, r, a, b, c, alpha, beta, gamma, length, a_dna, b_dna, c_dna, alpha_p, beta_p, gamma_p);


        for (i=0; i<n3; i++){
                a_n[i] = a[i];
                b_n[i] = b[i];
                c_n[i] = c[i];
        }


        //Place holder for charge beads, tail beads, linker histone beads.

        //Nucleosome charge beads

        int Nq, Nq3;
        std::vector<double> core_pos, core_q;

        Read_core(Nq, Nq3, core_pos, core_q, "core_data.reg.150mM");

        double core_pos_d[Nq3], core_q_d[Nq];


	// Prepare parameter for random rotation

        for (i=0; i<n; i++){
                if (type[i]!=0){
                        r_var[i*3] = 2.0*time_step*kbt/(8*PI*eta*125.0);
                        r_var[i*3+1] = 2.0*time_step*kbt/(8*PI*eta*125.0);
                        r_var[i*3+2] = 2.0*time_step*kbt/(8*PI*eta*125.0);
                }else{
                        r_var[i*3] = 2.0*time_step*kbt/(4*PI*eta*r_h*r_h*lo);
                        r_var[i*3+1] = 0.0;
                        r_var[i*3+2] = 0.0;
                }
        }


	double r_all[n3],rad_all[n];

	for (i=0;i<n;i++){
		r_all[i*3] = r[i*3];
                r_all[i*3+1] = r[i*3+1];
                r_all[i*3+2] = r[i*3+2];
                rad_all[i] = rad[i];
	}

	int n_D, n_D3;

        n_D = n;
        n_D3 = n_D*3;

	double p[n_D3],rd[n_D3];

	for (i=0;i<Nq3;i++){
                core_pos_d[i] = core_pos[i];
        }
	for (i=0;i<Nq;i++){
                core_q_d[i] = core_q[i];
        }

	cuda_application_init_D_Chol(n_D3);

	cuda_application_init_data(n_c, nc3, n, n3, type, r, a, b, c, alpha, beta, gamma, length, a_dna, b_dna, c_dna, alpha_p, beta_p, gamma_p, h, g, s, phi_o, debyell, debye, q_l, k_e, k_ex, k_h1, sigma_DNA_DNA, sigma_DNA_Core, sigma_Core_Core, sigma_Tail_Tail, sigma_Tail_Linker, sigma_Tail_Core, Nq, Nq3, core_pos_d, core_q_d, force, torque, E, r_all, rad_all, ex_n, ex_force_m);

	for (int step = 0; step < number_of_steps; step++){
	
		for (i=0; i<n; i++){
                	rad_all[i] = rad[i];
                }

                for (i = 0; i < n3; i++){
                        rr[i] = rand_normal(0,r_var[i]);
                }


                for (i = 0; i < n3; i++){
                        p[i] = rand_normal(0,1);
                }

		main_cuda(n_c, nc3, ex_n, step, number_of_steps, time_step, del, frequency_RP, frequency_of_sampling, h, g, s, debyell, debye, q_l, k_e, k_ex, k_h1, sigma_DNA_DNA, sigma_DNA_Core, sigma_Core_Core, sigma_Tail_Tail, sigma_Tail_Linker, sigma_Tail_Core, Nq, Nq3, n, n3, a1, a2, s2dt, rr, p, E, r, a, b, c, rad_all);	
		ofstream energyfile;
                energyfile.open("Energy.txt", std::ios_base::app);
                energyfile << E[0] << endl;
                energyfile.close();


                if (step%frequency_of_sampling == 0){
                        write_xyz_append(n, n3, n_c, Nq, Nq3, type, r, a, b, c, core_pos, "out.xyz");
                }

		if (step == number_of_steps-1){
			write_restart(n3, r, a, b, c);
		}
	}

	free_all();
	time_t end = time(0);

        cout << "Wall time: " << end-now << " seconds." << endl;

}
