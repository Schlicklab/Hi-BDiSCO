// **************************************************************************************//
//											 //
//		 		Brownian Dynamics Simulation Algorithm			 //
// 		     Copyright Zilong Li, Tamar Schlick and New York University 	 //
//					    April 2020					 //
//                                                                                       //
// **************************************************************************************//

#include "func.h"

using namespace std;

double twist_angle_phi(int w){

        double rise = 0.34;
        double B_DNA_twist = 34.95;
        double number_of_base_pairs = 0.0;
        double number_of_turns = 0.0;
        double DNA_rotation = 0.0;
        double rotation_per_base_pair = 0.0;
        double DNA_twist = 0.0;
        double whole_linker_twist = 0.0;
        double twist_per_segment = 0.0;
        double twist_per_segment_rad = 0.0;

        double twist_angle;

        if (w <= 11){
                switch (w){
                case 0:{
                           cout << "Wrong DNA linker length!" << endl;
                           exit(0);
                }
                case 1:{
                           cout << "Wrong DNA linker length!" << endl;
                           exit(0);
                }
                case 2:{
                           twist_angle = 0.9007;
                           break;
                }
                case 3:{
                           twist_angle = -0.6754;
                           break;
                }
                case 4:{
                           twist_angle = -0.3519;
                           break;
                }
                case 5:{
                           twist_angle = -0.1463;
                           break;
                }
                case 6:{
                           twist_angle = 0.0215;
                           break;
                }
                case 7:{
                           twist_angle = 0.1178;
                           break;
                }
                case 8:{
                           twist_angle = 0.2025;
                           break;
                }

                case 9:{
                           twist_angle = 0.2726;
                           break;
                }

                case 10:{
                            twist_angle = -0.2415;
                            break;
                }

                case 11:{
                            twist_angle = -0.15;
                            break;
                }


                }
            }
            else{
                number_of_base_pairs = (w + 1) * 3 / rise;
                number_of_turns = number_of_base_pairs / 10.3;
                DNA_rotation = round(number_of_turns);
                rotation_per_base_pair = DNA_rotation * 360 / number_of_base_pairs;
                DNA_twist = rotation_per_base_pair - B_DNA_twist;
                whole_linker_twist = number_of_base_pairs*DNA_twist;
                twist_per_segment = whole_linker_twist / (w + 1);
                twist_per_segment_rad = PI*twist_per_segment / 180.00;
                twist_angle = twist_per_segment_rad;
            }

        return twist_angle;

}

double Energy_Stretching(double h, std::vector<double> coord1, std::vector<double> coord2, double lo){

        double energy = 0.0;

        double distance = 0.0;

        distance = sqrt((coord1[0]-coord2[0])*(coord1[0]-coord2[0]) + (coord1[1]-coord2[1])*(coord1[1]-coord2[1]) + (coord1[2]-coord2[2])*(coord1[2]-coord2[2]));

        energy = h*(distance-lo)*(distance-lo)/2;

        return energy;
}

double Energy_Bending(double g, double beta, double beta_o){

        double energy = 0.0;

        energy = g*(beta-beta_o)*(beta-beta_o)/2;

        return energy;

}

double Energy_Twisting(double s, double alpha, double gamma, double phi_o, double lo){

        double energy = 0.0;

        energy = s*((alpha+gamma-phi_o)*(alpha+gamma-phi_o))/2;

        return energy;

}

void Force_Stretching(double h, double* coord1, double* coord2, double lo, double* force_projection1, double* force_projection2, double& Energy){

        double force = 0.0;

        double distance = 0.0;

        distance = sqrt((coord1[0]-coord2[0])*(coord1[0]-coord2[0]) + (coord1[1]-coord2[1])*(coord1[1]-coord2[1]) + (coord1[2]-coord2[2])*(coord1[2]-coord2[2]));

	Energy = Energy + h*(distance-lo)*(distance-lo)/2;

        force = -h*(distance-lo);

	for (int i = 0; i < 3; i++){
                force_projection1[i] = force*(coord1[i] - coord2[i])/distance;
		force_projection2[i] = -force*(coord1[i] - coord2[i])/distance;
        }

}

void Force_Bending(double g, double beta, double beta_o, double* coord1, double* coord2, double* coord3, double* force_projection1, double* force_projection2, double* force_projection3, double& Energy){

        double force = 0.0;

	double distance1, distance2;
        double norm_ri, norm_rk;

	std::vector<double> ji, jk, kj, ri, rk;

	Energy = Energy + g*(beta-beta_o)*(beta-beta_o)/2;

	force = -g*(beta-beta_o);

	distance1 = sqrt((coord1[0]-coord2[0])*(coord1[0]-coord2[0]) + (coord1[1]-coord2[1])*(coord1[1]-coord2[1]) + (coord1[2]-coord2[2])*(coord1[2]-coord2[2]));
        distance2 = sqrt((coord2[0]-coord3[0])*(coord2[0]-coord3[0]) + (coord2[1]-coord3[1])*(coord2[1]-coord3[1]) + (coord2[2]-coord3[2])*(coord2[2]-coord3[2]));


        ji = {coord2[0]-coord1[0], coord2[1]-coord1[1], coord2[2]-coord1[2]};
        jk = {coord2[0]-coord3[0], coord2[1]-coord3[1], coord2[2]-coord3[2]};
        kj = {coord3[0]-coord2[0], coord3[1]-coord2[1], coord3[2]-coord2[2]};

        ri = cross_product(ji, jk);
//	printf("host %f %f %f %f %f %f %f %f %f\n", ji[0], ji[1], ji[2],jk[0],jk[1],jk[2], ri[0],ri[1],ri[2]);
        ri = cross_product(ji, ri);

        norm_ri = norm(ri);

        for (int i = 0; i < 3; i++){
                force_projection1[i] = (force/distance1)*ri[i]/norm_ri;
		force_projection2[i] = -(force/distance1)*ri[i]/norm_ri;
        }

        rk = cross_product(ji, jk);
        rk = cross_product(kj, rk);

        norm_rk = norm(rk);

        for (int j = 0; j < 3; j++){
                force_projection2[j] = force_projection2[j] - (force/distance2)*rk[j]/norm_rk;
		force_projection3[j] = (force/distance2)*rk[j]/norm_rk;
        }

//	printf("host %f %f %f %f %f \n", distance1, distance2, ji[0], ri[0], norm_ri);
}

void Bending_force_projection(double g, double beta, double beta_b, double length, double* a_f, double* a_b, double* a, double* force_projection1, double* force_projection2, double& Energy){

	double Ai[3], Bi[3];
	double c1, c2, g1, g2;

	Energy = Energy + g*(beta-beta_b)*(beta-beta_b)/2;

	if (beta >= 1e-10){
		g1 = beta/(sin(beta)*length);
	}else{
		g1 = 1.0/length;
	}
	c1 = cos(beta);
	
	Ai[0] = g1*(a_f[0]-c1*a[0]);
	Ai[1] = g1*(a_f[1]-c1*a[1]);
	Ai[2] = g1*(a_f[2]-c1*a[2]);	

	if (beta_b >= 1e-10){
                g2 = beta_b/(sin(beta_b)*length);
        }else{
                g2 = 1.0/length;
        }
        c2 = cos(beta_b);

        Bi[0] = g2*(a_b[0]-c2*a[0]);
        Bi[1] = g2*(a_b[1]-c2*a[1]);
        Bi[2] = g2*(a_b[2]-c2*a[2]);

	for (int i = 0; i < 3; i++){
                force_projection1[i] = -g*(Ai[i]+Bi[i]);
                force_projection2[i] = g*(Ai[i]+Bi[i]);
        }

}

void Twisting_force_projection(double s, double alpha, double beta, double gamma, double phi_o, double length, double alpha_b, double beta_b, double gamma_b, double phi_o_b, double gamma_n, double* b, double* c, double* force_projection1, double* force_projection2, double& Energy){

	double Chi[3], Chim1[3], Zhi[3], Zhim1[3];
	double g1, g2, c1, c2, s1, s2;

	Energy = Energy + s*((alpha+gamma-phi_o)*(alpha+gamma-phi_o))/2;

	g1 = (alpha+gamma-phi_o)*tan(0.5*beta)/length;
	c1 = cos(alpha);
	s1 = sin(alpha);

	g2 = (alpha_b+gamma_b-phi_o_b)*tan(0.5*beta_b)/length;
	c2 = cos(gamma_n);
	s2 = sin(gamma_n);

	Chi[0] = g1*(c1*c[0]-s1*b[0]);
        Chi[1] = g1*(c1*c[1]-s1*b[1]);
        Chi[2] = g1*(c1*c[2]-s1*b[2]);
        Zhi[0] = g2*(c2*c[0]+s2*b[0]);
        Zhi[1] = g2*(c2*c[1]+s2*b[1]);
        Zhi[2] = g2*(c2*c[2]+s2*b[2]);

	for (int i = 0; i < 3; i++){
                force_projection1[i] = s*(Chi[i]+Zhi[i]);
                force_projection2[i] = -s*(Chi[i]+Zhi[i]);
        }
}

double Energy_Electrostatics(double q1, double q2, double epslon, double kappa, std::vector<double> coord1, std::vector<double> coord2){

        double energy = 0.0;
        double Rcut = 7.0;
        double distance = 0.0;

        distance = sqrt((coord1[0]-coord2[0])*(coord1[0]-coord2[0]) + (coord1[1]-coord2[1])*(coord1[1]-coord2[1]) + (coord1[2]-coord2[2])*(coord1[2]-coord2[2]));

        if (distance < Rcut){
                energy = (q1*q2/(4*PI*epslon*distance))*exp(-kappa*distance);
        }

        return energy;

}

void Force_Electrostatics(double q1, double q2, double epslon, double kappa, double* coord1, double* coord2, double* force_projection1, double* force_projection2, double& Energy){

        double force = 0.0;
        double Rcut = 7.0;
        double distance = 0.0;

        distance = sqrt((coord1[0]-coord2[0])*(coord1[0]-coord2[0]) + (coord1[1]-coord2[1])*(coord1[1]-coord2[1]) + (coord1[2]-coord2[2])*(coord1[2]-coord2[2]));

        if (distance < Rcut){
                force = ((q1*q2*(kappa*distance+1))/(4*PI*epslon*distance*distance))*exp(-kappa*distance);
		Energy = Energy + (q1*q2/(4*PI*epslon*distance))*exp(-kappa*distance);
        }

	for (int i = 0; i < 3; i++){
                force_projection1[i] = force*(coord1[i] - coord2[i])/distance;
                force_projection2[i] = -force*(coord1[i] - coord2[i])/distance;
		if (force_projection1[i]>0.9){ force_projection1[i]=0.9; }
		if (force_projection2[i]>0.9){ force_projection2[i]=0.9; }
		if (force_projection1[i]<-0.9){ force_projection1[i]=-0.9; }
                if (force_projection2[i]<-0.9){ force_projection2[i]=-0.9; }
        }

}

double Energy_Exclude_Volume(double k_ev, double sigma, std::vector<double> coord1, std::vector<double> coord2){

        double energy = 0.0;
        double vdw_cut = 4.0;
        double distance = 0.0;

        distance = sqrt((coord1[0]-coord2[0])*(coord1[0]-coord2[0]) + (coord1[1]-coord2[1])*(coord1[1]-coord2[1]) + (coord1[2]-coord2[2])*(coord1[2]-coord2[2]));

        if (distance < vdw_cut){
                energy = k_ev*(pow(sigma/distance,12) - pow(sigma/distance,6));
        }

        return energy;

}

void Force_Exclude_Volume(double k_ev, double sigma, double* coord1, double* coord2, double* force_projection1, double* force_projection2, double& Energy){

        double force = 0.0;
        double vdw_cut = 4.0;
        double distance = 0.0;

        distance = sqrt((coord1[0]-coord2[0])*(coord1[0]-coord2[0]) + (coord1[1]-coord2[1])*(coord1[1]-coord2[1]) + (coord1[2]-coord2[2])*(coord1[2]-coord2[2]));

        if (distance < vdw_cut){
                force = -k_ev*((6*pow(sigma,6))/pow(distance,7) - (12*pow(sigma,12))/pow(distance,13));
		Energy = Energy + k_ev*(pow(sigma/distance,12) - pow(sigma/distance,6));
        }

	for (int i = 0; i < 3; i++){
                force_projection1[i] = force*(coord1[i] - coord2[i])/distance;
                force_projection2[i] = -force*(coord1[i] - coord2[i])/distance;
		if (force_projection1[i]>0.9){ force_projection1[i]=0.9; }
                if (force_projection2[i]>0.9){ force_projection2[i]=0.9; }
                if (force_projection1[i]<-0.9){ force_projection1[i]=-0.9; }
                if (force_projection2[i]<-0.9){ force_projection2[i]=-0.9; }
        }

}

void Force_Ele_Vdw(double q1, double q2, double epslon, double kappa, double k_ev, double sigma, double* coord1, double* coord2, double* force_projection1, double* force_projection2, double& Energy){

        double force = 0.0;
	double force2 = 0.0;
	double Rcut = 7.0;
        double vdw_cut = 4.0;
        double distance = 0.0;

        distance = sqrt((coord1[0]-coord2[0])*(coord1[0]-coord2[0]) + (coord1[1]-coord2[1])*(coord1[1]-coord2[1]) + (coord1[2]-coord2[2])*(coord1[2]-coord2[2]));

	if (distance < Rcut){
                force = ((q1*q2*(kappa*distance+1))/(4*PI*epslon*distance*distance))*exp(-kappa*distance);
		Energy = Energy + (q1*q2/(4*PI*epslon*distance))*exp(-kappa*distance);
        }

        if (distance < vdw_cut){
                force2 = -k_ev*((6*pow(sigma,6))/pow(distance,7) - (12*pow(sigma,12))/pow(distance,13));
		Energy = Energy + k_ev*(pow(sigma/distance,12) - pow(sigma/distance,6));
        }

        for (int i = 0; i < 3; i++){
                force_projection1[i] = (force+force2)*(coord1[i] - coord2[i])/distance;
                force_projection2[i] = -(force+force2)*(coord1[i] - coord2[i])/distance;
		if (force_projection1[i]>0.9){ force_projection1[i]=0.9; }
                if (force_projection2[i]>0.9){ force_projection2[i]=0.9; }
                if (force_projection1[i]<-0.9){ force_projection1[i]=-0.9; }
                if (force_projection2[i]<-0.9){ force_projection2[i]=-0.9; }
        }

}


void torque_due_to_force(double* force, double* coord_f, double* coord_c, double* a, double* b, double* c, double* torque){

	//Given exact coordinates of the bead with force applied (coord_f) and the exact coordinate of the center bead (coord_c)

        double fa, fb, fc;
	double comp[3];

	comp[0] = (coord_f[0]-coord_c[0])*a[0] + (coord_f[1]-coord_c[1])*a[1] + (coord_f[2]-coord_c[2])*a[2];
	comp[1] = (coord_f[0]-coord_c[0])*b[0] + (coord_f[1]-coord_c[1])*b[1] + (coord_f[2]-coord_c[2])*b[2];
	comp[2] = (coord_f[0]-coord_c[0])*c[0] + (coord_f[1]-coord_c[1])*c[1] + (coord_f[2]-coord_c[2])*c[2];


        fa = a[0]*force[0] + a[1]*force[1] + a[2]*force[2];
        fb = b[0]*force[0] + b[1]*force[1] + b[2]*force[2];
        fc = c[0]*force[0] + c[1]*force[1] + c[2]*force[2];

        torque[0] = fc*comp[1] - fb*comp[2];
        torque[1] = fa*comp[2] - fc*comp[0];
        torque[2] = fb*comp[0] - fa*comp[1];

}

void torque_due_to_force_relative(double* force, double* comp, double* a, double* b, double* c, double* torque){

	//Given relative position (comp) of the bead with force applied

        double fa, fb, fc;

        fa = a[0]*force[0] + a[1]*force[1] + a[2]*force[2];
        fb = b[0]*force[0] + b[1]*force[1] + b[2]*force[2];
        fc = c[0]*force[0] + c[1]*force[1] + c[2]*force[2];

        torque[0] = fc*comp[1] - fb*comp[2];
        torque[1] = fa*comp[2] - fc*comp[0];
        torque[2] = fb*comp[0] - fa*comp[1];

}

void update_phi_o(int n, int* type, double* phi_o){

        int k, l;

        for (int i = 0; i < n; i++){
                if (type[i] == 1){
                        k = 1;
                        while (i+k < n and type[i+k] == 0){
                                k = k+1;
                        }

                        for (l = 0; l < k; l++){
                                phi_o[i+l] = twist_angle_phi(k-1);
                        }


                }
        }

}

void update_Euler_Angle(int n_c, int nc3, int n, int n3, int* type, double* r, double* a, double* b, double* c, double* alpha, double* beta, double* gamma, double* length, double* a_dna, double* b_dna, double* c_dna, double* alpha_p, double* beta_p, double* gamma_p){

	double r_forw[3];
	double mi;
	double da[3], a_old[3];
	double a_m[3], b_m[3];
	double Ac, apg, f1, f2, ada, bda, si, co;
	double sa, ca, sb, cb, sg, cg;	
	double R21, R22, R23, R31, R32, R33;
	int i, i1,i2,i3, if1,if2,if3, nm1, ic,ic1,ic2,ic3;
	int count;
	
	si = sin(theta);
	co = cos(theta);
	
	for (i = 0; i < n-1; i++){
		i1 = 3*i;
		i2 = i1+1;
		i3 = i2+1;
		if1 = i1+3;
		if2 = if1+1;
		if3 = if2+1;

		if (type[i] == 0){
			a_old[0] = a[i1];
			a_old[1] = a[i2];
			a_old[2] = a[i3];
			if (type[i+1] == 0){
				r_forw[0] = r[if1] - r[i1];
				r_forw[1] = r[if2] - r[i2];
				r_forw[2] = r[if3] - r[i3];
			}else{
				b_m[0] = -si*a[if1] + co*b[if1];
				b_m[1] = -si*a[if2] + co*b[if2];
				b_m[2] = -si*a[if3] + co*b[if3];
				r_forw[0] = r[if1] - ro*b_m[0] + d1*c[if1] - r[i1];
				r_forw[1] = r[if2] - ro*b_m[1] + d1*c[if2] - r[i2];
				r_forw[2] = r[if3] - ro*b_m[2] + d1*c[if3] - r[i3];
			}
			length[i] = sqrt(r_forw[0]*r_forw[0] + r_forw[1]*r_forw[1] + r_forw[2]*r_forw[2]);
			mi = 1.0/length[i];
			a[i1] = mi*r_forw[0];
			a[i2] = mi*r_forw[1];
			a[i3] = mi*r_forw[2];
			da[0] = a[i1] - a_old[0];
			da[1] = a[i2] - a_old[1];
			da[2] = a[i3] - a_old[2];
			
			bda = b[i1]*da[0] +b[i2]*da[1] + b[i3]*da[2];
			b[i1] = b[i1] - bda*a_old[0];
			b[i2] = b[i2] - bda*a_old[1];
			b[i3] = b[i3] - bda*a_old[2];

			bda = b[i1]*a[i1] + b[i2]*a[i2] + b[i3]*a[i3];
			b[i1] = b[i1] - bda*a[i1];
                        b[i2] = b[i2] - bda*a[i2];
                        b[i3] = b[i3] - bda*a[i3];

			mi = 1.0/sqrt(b[i1]*b[i1] + b[i2]*b[i2] + b[i3]*b[i3]);
			b[i1] = mi*b[i1];
			b[i2] = mi*b[i2];
			b[i3] = mi*b[i3];

			c[i1] = a[i2]*b[i3] - a[i3]*b[i2];
			c[i2] = a[i3]*b[i1] - a[i1]*b[i3];
			c[i3] = a[i1]*b[i2] - a[i2]*b[i1];
		}else{
			r_forw[0] = r[if1] - (r[i1]-ro*b[i1]-d1*c[i1]);
			r_forw[1] = r[if2] - (r[i2]-ro*b[i2]-d1*c[i2]);
			r_forw[2] = r[if3] - (r[i3]-ro*b[i3]-d1*c[i3]);
			length[i] = sqrt(r_forw[0]*r_forw[0] + r_forw[1]*r_forw[1] + r_forw[2]*r_forw[2]);
		}
	}

	for (i = 0; i < n-1; i++){
		i1 = 3*i;
                i2 = i1+1;
                i3 = i2+1;
                if1 = i1+3;
                if2 = if1+1;
                if3 = if2+1;

		if (type[i]==0){
			if(type[i+1]==0){
				ada = a[i1]*a[if1]+a[i2]*a[if2]+a[i3]*a[if3];
				if (ada > 1.0) ada = 1.0;
				if (ada < -1.0) ada = -1.0;
				beta[i] = acos(ada);
				sb = sin(beta[i]);
				if (beta[i] > 1e-10){
					f1 = (a[if1]*b[i1]+a[if2]*b[i2]+a[if3]*b[i3])/sb;
				}else{
					f1 = (b[if1]*b[i1]+b[if2]*b[i2]+b[if3]*b[i3]);
				}
				if (f1 > 1.0) f1 = 1.0;
				if (f1 < -1.0) f1 = -1.0;
				Ac = acos(f1);
				f2 = a[if1]*c[i1] + a[if2]*c[i2] + a[if3]*c[i3];
				if (f2 >= 0){
					alpha[i] = Ac;
				}else{
					alpha[i] = -Ac;
				}

				f1 = (b[i1]*b[if1]+b[i2]*b[if2]+b[i3]*b[if3]+c[i1]*c[if1]+c[i2]*c[if2]+c[i3]*c[if3])/(1.0 + ada);
				if (f1 > 1.0) f1 = 1.0;
                                if (f1 < -1.0) f1 = -1.0;
				apg = acos(f1);
				f2 = (c[i1]*b[if1]+c[i2]*b[if2]+c[i3]*b[if3]-(b[i1]*c[if1]+b[i2]*c[if2]+b[i3]*c[if3]))/(1.0 + ada);
				if (f2 >= 0.0){
                                        gamma[i] = apg - alpha[i];
                                }else{
                                        gamma[i] = -apg - alpha[i];
                                }
			}else{
				a_m[0] = co*a[if1] + si*b[if1];
				a_m[1] = co*a[if2] + si*b[if2];
				a_m[2] = co*a[if3] + si*b[if3];
				b_m[0] = -si*a[if1] + co*b[if1];
				b_m[1] = -si*a[if2] + co*b[if2];
				b_m[2] = -si*a[if3] + co*b[if3];

				ada = a[i1]*a_m[0]+a[i2]*a_m[1]+a[i3]*a_m[2];
				if (ada > 1.0) ada = 1.0;
                                if (ada < -1.0) ada = -1.0;
				beta[i] = acos(ada);
                                sb = sin(beta[i]);
                                if (beta[i] > 1e-10){
                                        f1 = (a_m[0]*b[i1]+a_m[1]*b[i2]+a_m[2]*b[i3])/sb;
                                }else{
                                        f1 = (b_m[0]*b[i1]+b_m[1]*b[i2]+b_m[2]*b[i3]);
                                }
				if (f1 > 1.0) f1 = 1.0;
                                if (f1 < -1.0) f1 = -1.0;
				Ac = acos(f1);
                                f2 = a_m[0]*c[i1] + a_m[1]*c[i2] + a_m[2]*c[i3];
                                if (f2 >= 0){
                                        alpha[i] = Ac;
                                }else{
                                        alpha[i] = -Ac;
                                }

                                f1 = (b[i1]*b_m[0]+b[i2]*b_m[1]+b[i3]*b_m[2]+c[i1]*c[if1]+c[i2]*c[if2]+c[i3]*c[if3])/(1.0 + ada);
                                if (f1 > 1.0) f1 = 1.0;
                                if (f1 < -1.0) f1 = -1.0;
                                apg = acos(f1);
                                f2 = (c[i1]*b_m[0]+c[i2]*b_m[1]+c[i3]*b_m[2]-(b[i1]*c[if1]+b[i2]*c[if2]+b[i3]*c[if3]))/(1.0 + ada);
                                if (f2 >= 0.0){
                                        gamma[i] = apg - alpha[i];
                                }else{
                                        gamma[i] = -apg - alpha[i];
                                }
			}
		}else{
			ic = 0;
			for (count = 0; count <= i; count++){
				if(type[count]==1) ic=ic+1;
			}
			ic = ic-1;
			ic1 = 3*ic;
			ic2 = ic1+1;
			ic3 = ic2+1;

			a_dna[ic1] = (r[if1] - (r[i1]-ro*b[i1]-d1*c[i1]))/length[i];
			a_dna[ic2] = (r[if2] - (r[i2]-ro*b[i2]-d1*c[i2]))/length[i];
			a_dna[ic3] = (r[if3] - (r[i3]-ro*b[i3]-d1*c[i3]))/length[i];

			cb = a[i1]*a_dna[ic1] + a[i2]*a_dna[ic2] + a[i3]*a_dna[ic3];
			if (cb > 1.0) cb = 1.0;
			if (cb < -1.0) cb = -1.0;
			beta_p[ic] = acos(cb);
			sb = sin(beta_p[ic]);
			if (beta_p[ic] >= 1e-10){
				b_m[0] = (a_dna[ic1]-cb*a[i1])/sb;
				b_m[1] = (a_dna[ic2]-cb*a[i2])/sb;
				b_m[2] = (a_dna[ic3]-cb*a[i3])/sb;
				ca = b_m[0]*b[i1]+b_m[1]*b[i2]+b_m[2]*b[i3];
				if (ca > 1.0) ca = 1.0;
                       		if (ca < -1.0) ca = -1.0;
				Ac = acos(ca);
				f1 = a_dna[ic1]*c[i1]+a_dna[ic2]*c[i2]+a_dna[ic3]*c[i3];
				if (f1 >= 0){
					alpha_p[ic] = Ac;
				}else{
					alpha_p[ic] = -Ac;
				}
				gamma_p[ic] = -alpha_p[ic];
				sa = sin(alpha_p[ic]);
				sg = sin(gamma_p[ic]);
				cg = cos(gamma_p[ic]);
				R21 = -cg*sb;
				R22 = cg*cb*ca-sg*sa;
				R23 = cg*cb*sa+sg*ca;

				b_dna[ic1] = R21*a[i1] + R22*b[i1] + R23*c[i1];
				b_dna[ic2] = R21*a[i2] + R22*b[i2] + R23*c[i2];
				b_dna[ic3] = R21*a[i3] + R22*b[i3] + R23*c[i3];

				R31 = sg*sb;
				R32 = -sg*cb*ca-cg*sa;
				R33 = -sg*cb*sa+cg*ca;

				c_dna[ic1] = R31*a[i1]+R32*b[i1]+R33*c[i1];
				c_dna[ic2] = R31*a[i2]+R32*b[i2]+R33*c[i2];
				c_dna[ic3] = R31*a[i3]+R32*b[i3]+R33*c[i3];
	

			}else{
				b_dna[ic1] = b[i1];
				b_dna[ic2] = b[i2];
				b_dna[ic3] = b[i3];
				c_dna[ic1] = c[i1];
				c_dna[ic2] = c[i2];
				c_dna[ic3] = c[i3];
			}
			ada = a_dna[ic1]*a[if1]+a_dna[ic2]*a[if2]+a_dna[ic3]*a[if3];
			if (ada > 1.0) ada = 1.0;
			if (ada < -1.0) ada = -1.0;
			beta[i] = acos(ada);
			sb  = sin(beta[i]);
			if (beta[i] >= 1e-10){
				f1 = (a[if1]*b_dna[ic1]+a[if2]*b_dna[ic2]+a[if3]*b_dna[ic3])/sb;
			}else{
				f1 = (b[if1]*b_dna[ic1]+b[if2]*b_dna[ic2]+b[if3]*b_dna[ic3]);
			}
			if (f1 > 1.0) f1 =1.0;
			if (f1 < -1.0) f1 = -1.0;
			Ac = acos(f1);
			f2 = a[if1]*c_dna[ic1]+a[if2]*c_dna[ic2]+a[if3]*c_dna[ic3];
			if (f2 >= 0){
				alpha[i] = Ac;
			}else{
				alpha[i] = -Ac;
			}

			f1 = (b_dna[ic1]*b[if1]+b_dna[ic2]*b[if2]+b_dna[ic3]*b[if3]+c_dna[ic1]*c[if1]+c_dna[ic2]*c[if2]+c_dna[ic3]*c[if3])/(1.0+ada);
			if (f1 > 1.0) f1 = 1.0;
			if (f1 < -1.0) f1 = -1.0;
			apg = acos(f1);
			f2 = (c_dna[ic1]*b[if1]+c_dna[ic2]*b[if2]+c_dna[ic3]*b[if3]-(b_dna[ic1]*c[if1]+b_dna[ic2]*c[if2]+b_dna[ic3]*c[if3]))/(1.0+ada);
			if (f2 >= 0.0){
				gamma[i] = apg - alpha[i];
			}else{
				gamma[i] = -apg - alpha[i];
			}

		}

	}



}

void first_coord(int t, double* r, double* a, double* b, double* c, double* r_f){

	double b_m[3];
	double si, co;

	si = sin(theta);
	co = cos(theta);

	if (t==0){
		r_f[0] = r[0];
		r_f[1] = r[1];
		r_f[2] = r[2];
	}else{
		b_m[0] = -si*a[0]+co*b[0];
		b_m[1] = -si*a[1]+co*b[1];
		b_m[2] = -si*a[2]+co*b[2];

		r_f[0] = r[0] - ro*b_m[0]+d1*c[0];
		r_f[1] = r[1] - ro*b_m[1]+d1*c[1];
		r_f[2] = r[2] - ro*b_m[2]+d1*c[2];
	}

}

void second_coord(int t, double* r, double* a, double* b, double* c, double* r_s){

	if (t==0){
		r_s[0] = r[0];
                r_s[1] = r[1];
                r_s[2] = r[2];
	}else{
		r_s[0] = r[0] - (ro*b[0]+d1*c[0]);
		r_s[1] = r[1] - (ro*b[1]+d1*c[1]);
		r_s[2] = r[2] - (ro*b[2]+d1*c[2]);
	}

}




void Diffusion_Tensor(int n, int n3, double* r, double a1, double a2, double* rad, std::vector<std::vector<double>> &D, std::vector<std::vector<double>> &Chol){

        double sij[3];
        double ssq, s, f, f1, f2;
        int i, j, k, l;
        int index1, index2, index3, index4;

	//Diffusion Tensor
	
	for (i=0;i<n;i++){
                for (j=i;j<n;j++){
                        if(i == j){
                                for (k=0;k<3;k++){
                                        for (l=0;l<3;l++){
                                                index1 = 3*i+k;
                                                index2 = 3*j+l;
                                                if(k==l){
                                                        D[index1][index2] = a1/rad[i];
                                                }else{
                                                        D[index1][index2] = 0.0;
                                                        D[index2][index1] = 0.0;
                                                }
                                        }
                                }
                        }else{
                                sij[0] = r[3*i] - r[3*j];
                                sij[1] = r[3*i+1] - r[3*j+1];
                                sij[2] = r[3*i+2] - r[3*j+2];
                                ssq = sij[0]*sij[0] + sij[1]*sij[1] + sij[2]*sij[2];
                                s = sqrt(ssq);
                               
				if (s >= rad[i]+rad[j]){
					for (k=0;k<3;k++){
						for(l=0;l<3;l++){
							index1 = 3*i+k;
							index2 = 3*j+l;
							index3 = 3*i+l;
							index4 = 3*j+k;
							f = (rad[i]*rad[i]+rad[j]*rad[j])/ssq;
							f1 = 1.0 + f/3;
							f2 = 1.0 -f;
							if(k==l){
								D[index1][index2] = (a2/s)*(f1+f2*sij[k]*sij[k]/ssq);
								D[index2][index1] = D[index1][index2];
							}else{
								D[index1][index2] = (a2/s)*(f2*sij[k]*sij[l]/ssq);
								D[index2][index1] = D[index1][index2];
								D[index3][index4] = D[index1][index2];
								D[index4][index3] = D[index1][index2];
							}
						}
					}
				}else{
					s = pow((rad[i]*rad[i]*rad[i]+rad[j]*rad[j]*rad[j])/2.0, 1.0/3);
					for (k=0;k<3;k++){
                                                for(l=0;l<3;l++){
                                                        index1 = 3*i+k;
                                                        index2 = 3*j+l;
                                                        index3 = 3*i+l;
                                                        index4 = 3*j+k;
                                                        if(k==l){
                                                                D[index1][index2] = (a1/s)*(1.0-9.0*sqrt(ssq)/(32*s) + 3.0*sij[k]*sij[k]/(32.0*s*sqrt(ssq)));
                                                                D[index2][index1] = D[index1][index2];
                                                        }else{
                                                                D[index1][index2] = (a1/s)*(3.0*sij[k]*sij[l]/(32.0*s*sqrt(ssq)));
                                                                D[index2][index1] = D[index1][index2];
                                                                D[index3][index4] = D[index1][index2];
                                                                D[index4][index3] = D[index1][index2];
                                                        }
                                                }
                                        }
				}
				
                                
                        }
                }
        }

	//Cholesky decomposition

	Chol[0][0] = sqrt(D[0][0]);

        for (i=1;i<n*3;i++){
                Chol[i][0] = D[i][0]/Chol[0][0];

                for (j=1;j<i;j++){
                        ssq = 0.0;
                        for (k=0;k<j;k++){
                                ssq = ssq+Chol[i][k]*Chol[j][k];
                        }
                        Chol[i][j] = (D[i][j]-ssq)/Chol[j][j];
                }
                ssq = 0.0;
                for (k=0;k<i;k++){
                        ssq = ssq + Chol[i][k]*Chol[i][k];
                }

                f = D[i][i] - ssq;

		if (f < 0){
			ofstream wrongDfile;
			wrongDfile.open("wrongD.txt", std::ios_base::app);
			wrongDfile << "wrong D here." << endl;
			wrongDfile.close();
			f=-f;
		}

                Chol[i][i] = sqrt(f);
        }



}

int Diffusion_Tensor_check(int n, int n3, double* r, double a1, double a2, double* rad, std::vector<std::vector<double>> &D, std::vector<std::vector<double>> &Chol){

        double sij[3];
        double ssq, s, f, f1, f2;
        int i, j, k, l;
	int D_recal = 1;
        int index1, index2, index3, index4;

        //Diffusion Tensor

        for (i=0;i<n;i++){
                for (j=i;j<n;j++){
                        if(i == j){
                                for (k=0;k<3;k++){
                                        for (l=0;l<3;l++){
                                                index1 = 3*i+k;
                                                index2 = 3*j+l;
                                                if(k==l){
                                                        D[index1][index2] = a1/rad[i];
                                                }else{
                                                        D[index1][index2] = 0.0;
                                                        D[index2][index1] = 0.0;
                                                }
                                        }
                                }
                        }else{
                                sij[0] = r[3*i] - r[3*j];
                                sij[1] = r[3*i+1] - r[3*j+1];
                                sij[2] = r[3*i+2] - r[3*j+2];
                                ssq = sij[0]*sij[0] + sij[1]*sij[1] + sij[2]*sij[2];
                                s = sqrt(ssq);

                                if (s >= rad[i]+rad[j]){
                                        for (k=0;k<3;k++){
                                                for(l=0;l<3;l++){
                                                        index1 = 3*i+k;
                                                        index2 = 3*j+l;
                                                        index3 = 3*i+l;
                                                        index4 = 3*j+k;
                                                        f = (rad[i]*rad[i]+rad[j]*rad[j])/ssq;
                                                        f1 = 1.0 + f/3;
                                                        f2 = 1.0 -f;
                                                        if(k==l){
                                                                D[index1][index2] = (a2/s)*(f1+f2*sij[k]*sij[k]/ssq);
                                                                D[index2][index1] = D[index1][index2];
                                                        }else{
                                                                D[index1][index2] = (a2/s)*(f2*sij[k]*sij[l]/ssq);
                                                                D[index2][index1] = D[index1][index2];
                                                                D[index3][index4] = D[index1][index2];
                                                                D[index4][index3] = D[index1][index2];
                                                        }
                                                }
                                        }
                                }else{
					s = pow((rad[i]*rad[i]*rad[i]+rad[j]*rad[j]*rad[j])/2.0, 1.0/3);
                                        for (k=0;k<3;k++){
                                                for(l=0;l<3;l++){
                                                        index1 = 3*i+k;
                                                        index2 = 3*j+l;
                                                        index3 = 3*i+l;
                                                        index4 = 3*j+k;
                                                        if(k==l){
                                                                D[index1][index2] = (a1/s)*(1.0-9.0*sqrt(ssq)/(32*s) + 3.0*sij[k]*sij[k]/(32.0*s*sqrt(ssq)));
                                                                D[index2][index1] = D[index1][index2];
                                                        }else{
                                                                D[index1][index2] = (a1/s)*(3.0*sij[k]*sij[l]/(32.0*s*sqrt(ssq)));
                                                                D[index2][index1] = D[index1][index2];
                                                                D[index3][index4] = D[index1][index2];
                                                                D[index4][index3] = D[index1][index2];
                                                        }
                                                }
                                        }
                                }


                        }
                }
        }

        //Cholesky decomposition

        Chol[0][0] = sqrt(D[0][0]);

        for (i=1;i<n*3;i++){
                Chol[i][0] = D[i][0]/Chol[0][0];

                for (j=1;j<i;j++){
                        ssq = 0.0;
                        for (k=0;k<j;k++){
                                ssq = ssq+Chol[i][k]*Chol[j][k];
                        }
                        Chol[i][j] = (D[i][j]-ssq)/Chol[j][j];
                }
                ssq = 0.0;
                for (k=0;k<i;k++){
                        ssq = ssq + Chol[i][k]*Chol[i][k];
                }

                f = D[i][i] - ssq;

                if (f<0) D_recal=0;

                Chol[i][i] = sqrt(f);
        }

	return D_recal;

}

void mechanical_force_and_torque(int n_c, int nc3, int n, int n3, int* type, double* r, double* a, double* b, double* c, double* alpha, double* beta, double* gamma, double* length, double* a_dna, double* b_dna, double* c_dna, double* alpha_p, double* beta_p, double* gamma_p, double h, double g, double s, double* phi_o, double* force, double* torque, double& Energy, int rank, int comm_size){

	int i, j, k, l, ic, count;
	int i1, i2, i3, j1, j2, j3, k1, k2, k3, l1, l2, l3, ib1, ib2, ib3, if1, if2, if3, im1, im2, im3, ic1, ic2, ic3;
	double c1, c2, s1, s2, g1, g2, dist, si, co;
	double a_m[3], b_m[3];
	double alpha_b, beta_b, gamma_b, phi_o_b, gamma_n;
	double a_f[3], a_b[3], a_o[3], b_o[3], c_o[3];
	double df[3];
	double Stri[3], Strim1[3], Ai[3], Aim1[3], Bi[3], Bim1[3], Chi[3], Chim1[3], Zhi[3], Zhim1[3];
	double mag, fa, fb, fc;
	double ada, adb, adc, cda, bda, cdb, bdb, cdc, bdc;

	double force_Stretching, force_Bending, force_Twisting;

	double torque1[3], torque2[3], torque3[3];
	double force_projection1[3], force_projection2[3], force_projection3[3];
	double r_f[3], r_s[3], r_f2[3];

	double r_tmp1[3], r_tmp2[3], a_tmp1[3], a_tmp2[3], b_tmp1[3], b_tmp2[3], c_tmp1[3], c_tmp2[3];


	si = sin(theta);
	co = cos(theta);

	//Stretching

        for (i=0;i<n-1;i++){

		if (i%comm_size != rank) continue;

		for (int xi = 0; xi <3; xi++){
			r_tmp1[xi] = r[i*3+3+xi];
			r_tmp2[xi] = r[i*3+xi];
			a_tmp1[xi] = a[i*3+3+xi];
			a_tmp2[xi] = a[i*3+xi];
			b_tmp1[xi] = b[i*3+3+xi];
                        b_tmp2[xi] = b[i*3+xi];
			c_tmp1[xi] = c[i*3+3+xi];
                        c_tmp2[xi] = c[i*3+xi];
		}

		first_coord(type[i+1], r_tmp1, a_tmp1, b_tmp1, c_tmp1, r_f);
		second_coord(type[i], r_tmp2, a_tmp2, b_tmp2, c_tmp2, r_s);

		Force_Stretching(h,r_s,r_f,lo,force_projection1,force_projection2, Energy);

		force[i*3] = force[i*3] + force_projection1[0];
                force[i*3+1] = force[i*3+1] + force_projection1[1];
                force[i*3+2] = force[i*3+2] + force_projection1[2];
                force[i*3+3] = force[i*3+3] + force_projection2[0];
                force[i*3+4] = force[i*3+4] + force_projection2[1];
                force[i*3+5] = force[i*3+5] + force_projection2[2];

//		printf("host: %d %f %f %f %f %f %f \n", i, force_projection1[0], force_projection1[1], force_projection1[2], force_projection2[0], force_projection2[1], force_projection2[2]);

		torque_due_to_force(force_projection1, r_s, r_tmp2, a_tmp2, b_tmp2, c_tmp2, torque1);
		torque_due_to_force(force_projection2, r_f, r_tmp1, a_tmp1, b_tmp1, c_tmp1, torque2);

		torque[i*3] = torque[i*3] + torque1[0];
                torque[i*3+1] = torque[i*3+1] + torque1[1];
                torque[i*3+2] = torque[i*3+2] + torque1[2];

		torque[i*3+3] = torque[i*3+3] + torque2[0];
                torque[i*3+4] = torque[i*3+4] + torque2[1];
                torque[i*3+5] = torque[i*3+5] + torque2[2];

        }


	//Bending
	for (i = 0; i < n-1; i++){

		if (i%comm_size != rank) continue;

		ic = 0;
                for (count = 0; count <= i; count++){
                        if(type[count]==1) ic=ic+1;
                }
                ic = ic-1;
                ic1 = 3*ic;
                ic2 = ic1+1;
                ic3 = ic2+1;

		for (int xi = 0; xi <3; xi++){
                        r_tmp1[xi] = r[i*3+3+xi];
                        r_tmp2[xi] = r[i*3+xi];
                        a_tmp1[xi] = a[i*3+3+xi];
                        a_tmp2[xi] = a[i*3+xi];
                        b_tmp1[xi] = b[i*3+3+xi];
                        b_tmp2[xi] = b[i*3+xi];
                        c_tmp1[xi] = c[i*3+3+xi];
                        c_tmp2[xi] = c[i*3+xi];
                }

		first_coord(type[i+1], r_tmp1, a_tmp1, b_tmp1, c_tmp1, r_f);
                second_coord(type[i], r_tmp2, a_tmp2, b_tmp2, c_tmp2, r_s);
		
		if (type[i]==0){
                        if (type[i+1]==0){
				for (int xi = 0; xi<3; xi++){
					a_f[xi] = a[i*3+3+xi];
				}
                        }else{
                                a_m[0] = co*a[i*3+3] + si*b[i*3+3];
                                a_m[1] = co*a[i*3+4] + si*b[i*3+4];
                                a_m[2] = co*a[i*3+5] + si*b[i*3+5];
				for (int xi = 0; xi<3; xi++){
                                        a_f[xi] = a_m[xi];
                                }
                        }
                        if (type[i-1]==0){
				for (int xi = 0; xi<3; xi++){
                                        a_b[xi] = a[i*3-3+xi];
                                }
                        }else{
				for (int xi = 0; xi<3; xi++){
                                        a_b[xi] = a_dna[ic1+xi];
                                }
                        }
			for (int xi = 0; xi<3; xi++){
                                a_o[xi] = a[i*3+xi];
                        }
			beta_b = beta[i-1];
                                       
                }else{
			for (int xi = 0; xi<3; xi++){
				a_f[xi] = a[i*3+3+xi];
				a_b[xi] = a[i*3+xi];
                                a_o[xi] = a_dna[ic1+xi];
                        }
			beta_b = beta_p[ic];
                }
		
		Bending_force_projection(g, beta[i], beta_b, length[i], a_f, a_b, a_o, force_projection1, force_projection2, Energy);
		
		force[i*3] = force[i*3] + force_projection1[0];
                force[i*3+1] = force[i*3+1] + force_projection1[1];
                force[i*3+2] = force[i*3+2] + force_projection1[2];
                force[i*3+3] = force[i*3+3] + force_projection2[0];
                force[i*3+4] = force[i*3+4] + force_projection2[1];
                force[i*3+5] = force[i*3+5] + force_projection2[2];

		torque_due_to_force(force_projection1, r_s, r_tmp2, a_tmp2, b_tmp2, c_tmp2, torque1);
                torque_due_to_force(force_projection2, r_f, r_tmp1, a_tmp1, b_tmp1, c_tmp1, torque2);

                torque[i*3] = torque[i*3] + torque1[0];
                torque[i*3+1] = torque[i*3+1] + torque1[1];
                torque[i*3+2] = torque[i*3+2] + torque1[2];

                torque[i*3+3] = torque[i*3+3] + torque2[0];
                torque[i*3+4] = torque[i*3+4] + torque2[1];
                torque[i*3+5] = torque[i*3+5] + torque2[2];

	}


	//Twisting
	
	for (i = 0; i < n-1; i++){

		if (i%comm_size != rank) continue;

                ic = 0;
                for (count = 0; count <= i; count++){
                        if(type[count]==1) ic=ic+1;
                }
                ic = ic-1;
                ic1 = 3*ic;
                ic2 = ic1+1;
                ic3 = ic2+1;

		for (int xi = 0; xi <3; xi++){
                        r_tmp1[xi] = r[i*3+3+xi];
                        r_tmp2[xi] = r[i*3+xi];
                        a_tmp1[xi] = a[i*3+3+xi];
                        a_tmp2[xi] = a[i*3+xi];
                        b_tmp1[xi] = b[i*3+3+xi];
                        b_tmp2[xi] = b[i*3+xi];
                        c_tmp1[xi] = c[i*3+3+xi];
                        c_tmp2[xi] = c[i*3+xi];
                }

                first_coord(type[i+1], r_tmp1, a_tmp1, b_tmp1, c_tmp1, r_f);
                second_coord(type[i], r_tmp2, a_tmp2, b_tmp2, c_tmp2, r_s);

                if (type[i]==0){
			alpha_b = alpha[i-1];
			beta_b = beta[i-1];
			gamma_b = gamma[i-1];
			phi_o_b = phi_o[i-1];
			gamma_n = gamma[i-1];
			for (int xi = 0; xi <3; xi++){
				b_o[xi] = b[i*3+xi];
				c_o[xi] = c[i*3+xi];
			}
                }else{
			alpha_b = alpha[i];
			beta_b = beta_p[ic];
                        gamma_b = gamma[i];
                        phi_o_b = phi_o[i];
                        gamma_n = gamma_p[ic];
			for (int xi = 0; xi <3; xi++){
                                b_o[xi] = b_dna[ic1+xi];
                                c_o[xi] = c_dna[ic1+xi];
                        }
                }

                Twisting_force_projection(s, alpha[i], beta[i], gamma[i], phi_o[i],  length[i], alpha_b, beta_b, gamma_b, phi_o_b, gamma_n, b_o, c_o, force_projection1, force_projection2, Energy);

                force[i*3] = force[i*3] + force_projection1[0];
                force[i*3+1] = force[i*3+1] + force_projection1[1];
                force[i*3+2] = force[i*3+2] + force_projection1[2];
                force[i*3+3] = force[i*3+3] + force_projection2[0];
                force[i*3+4] = force[i*3+4] + force_projection2[1];
                force[i*3+5] = force[i*3+5] + force_projection2[2];

		torque_due_to_force(force_projection1, r_s, r_tmp2, a_tmp2, b_tmp2, c_tmp2, torque1);
                torque_due_to_force(force_projection2, r_f, r_tmp1, a_tmp1, b_tmp1, c_tmp1, torque2);

                torque[i*3] = torque[i*3] + torque1[0];
                torque[i*3+1] = torque[i*3+1] + torque1[1];
                torque[i*3+2] = torque[i*3+2] + torque1[2];

                torque[i*3+3] = torque[i*3+3] + torque2[0];
                torque[i*3+4] = torque[i*3+4] + torque2[1];
                torque[i*3+5] = torque[i*3+5] + torque2[2];

        }



	//Mechanical Torques

	for (i = 0; i < n-1; i++){

		if (i%comm_size != rank) continue;

		im1 = i-1;
		i1 = i*3;
		i2 = i1+1;
		i3 = i2+1;
		ib1 = i1-3;
		ib2 = ib1+1;
		ib3 = ib2+1;

		ic = 0;
                for (count = 0; count <= i; count++){
                        if(type[count]==1) ic=ic+1;
                }
                ic = ic-1;
                ic1 = 3*ic;
                ic2 = ic1+1;
                ic3 = ic2+1;

		if (type[i]==0){
			torque[i1] = s*(alpha[i]+gamma[i]+phi_o[i]-alpha[im1]-gamma[im1]-phi_o[im1]);
			torque[i2] = 0.0;
			torque[i3] = 0.0;
		}else{
			ada = a_dna[ic1]*a[i1] + a_dna[ic2]*a[i2] + a_dna[ic3]*a[i3];
			adb = a_dna[ic1]*b[i1] + a_dna[ic2]*b[i2] + a_dna[ic3]*b[i3];
			adc = a_dna[ic1]*c[i1] + a_dna[ic2]*c[i2] + a_dna[ic3]*c[i3];
			
			mag = s*(alpha[i]+gamma[i]-phi_o[i]);
			torque[i1] = torque[i1] + mag*ada;
			torque[i2] = torque[i2] + mag*adb;
			torque[i3] = torque[i3] + mag*adc;
			if (i > 0){
				mag = -s*(alpha[im1]+gamma[im1]-phi_o[im1]);
				torque[i1] = torque[i1] + mag*co;
				torque[i2] = torque[i2] + mag*si;
				torque[i3] = torque[i3] + 0.0;
			}

			//Extra Bending torque

			torque[i2] = torque[i2] - g*beta_p[ic]*adc/sin(beta_p[ic]);
			torque[i3] = torque[i3] + g*beta_p[ic]*adb/sin(beta_p[ic]);
		
			if (i > 0){
				ada = a[ib1]*a[i1] + a[ib2]*a[i2] + a[ib3]*a[i3];
				adb = a[ib1]*b[i1] + a[ib2]*b[i2] + a[ib3]*b[i3];
				adc = a[ib1]*c[i1] + a[ib2]*c[i2] + a[ib3]*c[i3];
				
				torque[i1] = torque[i1] + g*beta[im1]*(si*adc)/sin(beta[im1]);
				torque[i2] = torque[i2] - g*beta[im1]*(co*adc)/sin(beta[im1]);
				torque[i3] = torque[i3] + g*beta[im1]*(co*adb-si*ada)/sin(beta[im1]);
			}

			//Extra Twisting torque

			s1 = sin(alpha_p[ic]);
			c1 = cos(alpha_p[ic]);
			mag = s*(alpha[i]+gamma[i]-phi_o[i])*tan(0.5*beta_p[ic]);
			cda = c_dna[ic1]*a[i1] + c_dna[ic2]*a[i2] + c_dna[ic3]*a[i3];
			bda = b_dna[ic1]*a[i1] + b_dna[ic2]*a[i2] + b_dna[ic3]*a[i3];
			cdb = c_dna[ic1]*b[i1] + c_dna[ic2]*b[i2] + c_dna[ic3]*b[i3];
			bdb = b_dna[ic1]*b[i1] + b_dna[ic2]*b[i2] + b_dna[ic3]*b[i3];
			cdc = c_dna[ic1]*c[i1] + c_dna[ic2]*c[i2] + c_dna[ic3]*c[i3];
			bdc = b_dna[ic1]*c[i1] + b_dna[ic2]*c[i2] + b_dna[ic3]*c[i3];

			torque[i1] = torque[i1] - mag*(s1*cda + c1*bda);
			torque[i2] = torque[i2] - mag*(s1*cdb + c1*bdb);
			torque[i3] = torque[i3] - mag*(s1*cdc + c1*bdc);

			if (i > 0){
				s1 = sin(gamma[im1]);
				c1 = cos(gamma[im1]);
				mag = s*(alpha[im1]+gamma[im1]-phi_o[im1])*tan(0.5*beta[im1]);
				cda = 0.0;
				bda = -si;
				cdb = 0.0;
				bdb = co;
				cdc = 1.0;
				bdc = 0.0;
				
				torque[i1] = torque[i1] - mag*(s1*cda - c1*bda);
                        	torque[i2] = torque[i2] - mag*(s1*cdb - c1*bdb);
                        	torque[i3] = torque[i3] - mag*(s1*cdc - c1*bdc);
			}

		}


	
	}

	//Additional torque for last bead

	if (n3%comm_size == rank){

		torque[n3-3] = -s*(alpha[n-2]+gamma[n-2]-phi_o[n-2]);
		torque[n3-2] = 0.0;
		torque[n3-1] = 0.0;

	}

}

void Electrostatic_and_Excluded_volume_force(int n, int n3, int n_c, int nc3, int* type, double* r, double* a, double* b, double* c, double debyell, double debye, double q_l, double k_e, double k_ex, double k_h1, double sigma_DNA_DNA, double sigma_DNA_Core, double sigma_Core_Core, int Nq, int Nq3, std::vector<double> core_pos, std::vector<double> core_q, double* force, double* torque, double& Energy, int rank, int comm_size){

	double ql_ql, dist;
	int i, j, k, l, ch;
	int i1, i2, i3, j1, j2, j3, k1, k2, k3, l1, l2, l3;
	double mi, Rcut, l_h1;
	double z[3];
	double fa, fb, fc;
	double g1, g2, s1, s2, mag;

	Rcut = 25.0;
	ql_ql = q_l*q_l;
	
	for (i = 0; i < n-1; i++){
		i1 = 3*i;
		i2 = i1+1;
		i3 = i2+1;
		for (j = i+1; j < n; j++){
			if ((i+j)%comm_size != rank) continue;
			j1 = j*3;
			j2 = j1+1;
			j3 = j2+1;

			dist = sqrt((r[j1]-r[i1])*(r[j1]-r[i1])+(r[j2]-r[i2])*(r[j2]-r[i2])+(r[j3]-r[i3])*(r[j3]-r[i3]));

			ch = 1;
			if (dist > Rcut) ch=0;
			if (type[i+1] != 0) ch=0;
			if (i > 0 and type[i-1] != 0) ch=0;
			if (j < n-1){
				if (type[j+1] != 0) ch=0;
			}
			if (type[j-1] != 0) ch=0;

			if (ch == 1){
				if (type[i] == 0){
					if (type[j] == 0){
						//DNA-DNA interaction
						if (abs(i-j) > 1){
							mi = 1.0/dist;
							z[0] = mi*(r[i1]-r[j1]);
							z[1] = mi*(r[i2]-r[j2]);
							z[2] = mi*(r[i3]-r[j3]);
			
							g1 = k_e*ql_ql*exp(-debyell*dist)*(debyell*dist+1.0)/(dist*dist);
							Energy = Energy + (k_e*ql_ql/dist)*exp(-debyell*dist);


							if (dist <= 8){ 
								s1 = sigma_DNA_DNA;
								s2 = sigma_DNA_DNA;
								g1 = g1 + k_ex*((12.0/s1)*pow((s1/dist),13)-(6.0/s2)*pow((s2/dist),7));
								Energy = Energy + k_ex*(pow(s1/dist,12) - pow(s2/dist,6));
							}

							force[i1] = force[i1] + g1*z[0];
							force[i2] = force[i2] + g1*z[1];
							force[i3] = force[i3] + g1*z[2];
							force[j1] = force[j1] - g1*z[0];
							force[j2] = force[j2] - g1*z[1];
							force[j3] = force[j3] - g1*z[2];

//							printf("host %d %d %f %f %f \n", i, j, g1*z[0], g1*z[1], g1*z[2]);
						}
					}else{
						//DNA-Core interaction

//						printf("host %d %d \n", i, j);

						for (k=0;k<Nq;k++){
							k1 = 3*k;
							k2 = k1+1;
							k3 = k2+1;
							z[0] = (r[i1]-(r[j1]+a[j1]*core_pos[k1]+b[j1]*core_pos[k2]+c[j1]*core_pos[k3]));
                                                        z[1] = (r[i2]-(r[j2]+a[j2]*core_pos[k1]+b[j2]*core_pos[k2]+c[j2]*core_pos[k3]));
                                                        z[2] = (r[i3]-(r[j3]+a[j3]*core_pos[k1]+b[j3]*core_pos[k2]+c[j3]*core_pos[k3]));
                                                        dist = sqrt(z[0]*z[0]+z[1]*z[1]+z[2]*z[2]);
							mi = 1.0/dist;
							z[0] = mi*z[0];
							z[1] = mi*z[1];
							z[2] = mi*z[2];

//							printf("host %d %d %d %f %f %f \n", i, j, k, z[0], z[1], z[2]);

							if (abs(i-j) > 1){
								g1 = k_e*q_l*core_q[k]*exp(-debye*dist)*(debye*dist+1.0)/(dist*dist);
								Energy = Energy + (k_e*q_l*core_q[k]/dist)*exp(-debye*dist);
								
								force[i1] = force[i1] + g1*z[0];
                                                        	force[i2] = force[i2] + g1*z[1];
                                                        	force[i3] = force[i3] + g1*z[2];
                                                        	force[j1] = force[j1] - g1*z[0];
                                                        	force[j2] = force[j2] - g1*z[1];
                                                        	force[j3] = force[j3] - g1*z[2];

//								printf("host %d %d %d %f %f %f %f \n", i, j, k, core_q[k], dist, exp(-debye*dist), g1);

								//torque due to force
								fa = -g1*(a[j1]*z[0]+a[j2]*z[1]+a[j3]*z[2]);
								fb = -g1*(b[j1]*z[0]+b[j2]*z[1]+b[j3]*z[2]);
								fc = -g1*(c[j1]*z[0]+c[j2]*z[1]+c[j3]*z[2]);
								torque[j1] = torque[j1] + fc*core_pos[k2] - fb*core_pos[k3];
								torque[j2] = torque[j2] + fa*core_pos[k3] - fc*core_pos[k1];
								torque[j3] = torque[j3] + fb*core_pos[k1] - fa*core_pos[k2];

//								printf("host %d %d %d %f %f %f \n", i, j, k, fc*core_pos[k2] - fb*core_pos[k3], fa*core_pos[k3] - fc*core_pos[k1], fb*core_pos[k1] - fa*core_pos[k2]);
							}

							// Excluded Volume force

							if (dist <= 8.0 and core_q[k]>0){
								s1 = sigma_DNA_Core;
								s2 = sigma_DNA_Core;
								g1 = k_ex*((12.0/s1)*pow((s1/dist),13)-(6.0/s2)*pow((s2/dist),7));
								Energy = Energy + k_ex*(pow(s1/dist,12) - pow(s2/dist,6));

								force[i1] = force[i1] + g1*z[0];
                                                                force[i2] = force[i2] + g1*z[1];
                                                                force[i3] = force[i3] + g1*z[2];
                                                                force[j1] = force[j1] - g1*z[0];
                                                                force[j2] = force[j2] - g1*z[1];
                                                                force[j3] = force[j3] - g1*z[2];

								//torque due to force
								fa = -g1*(a[j1]*z[0]+a[j2]*z[1]+a[j3]*z[2]);
                                                                fb = -g1*(b[j1]*z[0]+b[j2]*z[1]+b[j3]*z[2]);
                                                                fc = -g1*(c[j1]*z[0]+c[j2]*z[1]+c[j3]*z[2]);
                                                                torque[j1] = torque[j1] + fc*core_pos[k2] - fb*core_pos[k3];
                                                                torque[j2] = torque[j2] + fa*core_pos[k3] - fc*core_pos[k1];
                                                                torque[j3] = torque[j3] + fb*core_pos[k1] - fa*core_pos[k2];

							}


						}

					}
				}else{

					//Core-DNA interaction
					if (type[j] == 0){
						for (k=0;k<Nq;k++){
							k1 = 3*k;
							k2 = k1+1;
							k3 = k2+1;
							z[0] = (-r[j1]+(r[i1]+a[i1]*core_pos[k1]+b[i1]*core_pos[k2]+c[i1]*core_pos[k3]));
                                                        z[1] = (-r[j2]+(r[i2]+a[i2]*core_pos[k1]+b[i2]*core_pos[k2]+c[i2]*core_pos[k3]));
                                                        z[2] = (-r[j3]+(r[i3]+a[i3]*core_pos[k1]+b[i3]*core_pos[k2]+c[i3]*core_pos[k3]));
                                                        dist = sqrt(z[0]*z[0]+z[1]*z[1]+z[2]*z[2]);

							mi = 1.0/dist;
							z[0] = mi*z[0];
                                                        z[1] = mi*z[1];
                                                        z[2] = mi*z[2];
							if (abs(i-j) > 1){
								g1 = k_e*q_l*core_q[k]*exp(-debye*dist)*(debye*dist+1.0)/(dist*dist);
								Energy = Energy + (k_e*q_l*core_q[k]/dist)*exp(-debye*dist);
								
								force[i1] = force[i1] + g1*z[0];
                                                                force[i2] = force[i2] + g1*z[1];
                                                                force[i3] = force[i3] + g1*z[2];
                                                                force[j1] = force[j1] - g1*z[0];
                                                                force[j2] = force[j2] - g1*z[1];
                                                                force[j3] = force[j3] - g1*z[2];

								//torque due to force
								fa = g1*(a[i1]*z[0]+a[i2]*z[1]+a[i3]*z[2]);
                                                                fb = g1*(b[i1]*z[0]+b[i2]*z[1]+b[i3]*z[2]);
                                                                fc = g1*(c[i1]*z[0]+c[i2]*z[1]+c[i3]*z[2]);
                                                                torque[i1] = torque[i1] + fc*core_pos[k2] - fb*core_pos[k3];
                                                                torque[i2] = torque[i2] + fa*core_pos[k3] - fc*core_pos[k1];
                                                                torque[i3] = torque[i3] + fb*core_pos[k1] - fa*core_pos[k2];
							}

							//Excluded Volume force
							if (dist <= 8.0 and core_q[k]>0){
								s1 = sigma_DNA_Core;
								s2 = sigma_DNA_Core;
								Energy = Energy + k_ex*(pow(s1/dist,12) - pow(s2/dist,6));
								g1 = k_ex*((12.0/s1)*pow((s1/dist),13)-(6.0/s2)*pow((s2/dist),7));

                                                                force[i1] = force[i1] + g1*z[0];
                                                                force[i2] = force[i2] + g1*z[1];
                                                                force[i3] = force[i3] + g1*z[2];
                                                                force[j1] = force[j1] - g1*z[0];
                                                                force[j2] = force[j2] - g1*z[1];
                                                                force[j3] = force[j3] - g1*z[2];

								//torque due to force
								fa = g1*(a[i1]*z[0]+a[i2]*z[1]+a[i3]*z[2]);
                                                                fb = g1*(b[i1]*z[0]+b[i2]*z[1]+b[i3]*z[2]);
                                                                fc = g1*(c[i1]*z[0]+c[i2]*z[1]+c[i3]*z[2]);
                                                                torque[i1] = torque[i1] + fc*core_pos[k2] - fb*core_pos[k3];
                                                                torque[i2] = torque[i2] + fa*core_pos[k3] - fc*core_pos[k1];
                                                                torque[i3] = torque[i3] + fb*core_pos[k1] - fa*core_pos[k2];
							}
						}

					}else{//Core-Core interaction

						for (k = 0; k < Nq; k++){
							k1 = 3*k;
							k2 = k1+1;
							k3 = k2+1;
							for (l = 0; l < Nq; l++){
								l1 = 3*l;
								l2 = l1+1;
								l3 = l2+1;
								z[0] = (r[i1]+a[i1]*core_pos[k1]+b[i1]*core_pos[k2]+c[i1]*core_pos[k3] - (r[j1]+a[j1]*core_pos[l1]+b[j1]*core_pos[l2]+c[j1]*core_pos[l3]));
                                                                z[1] = (r[i2]+a[i2]*core_pos[k1]+b[i2]*core_pos[k2]+c[i2]*core_pos[k3] - (r[j2]+a[j2]*core_pos[l1]+b[j2]*core_pos[l2]+c[j2]*core_pos[l3]));
                                                                z[2] = (r[i3]+a[i3]*core_pos[k1]+b[i3]*core_pos[k2]+c[i3]*core_pos[k3] - (r[j3]+a[j3]*core_pos[l1]+b[j3]*core_pos[l2]+c[j3]*core_pos[l3]));
                                                                dist = sqrt(z[0]*z[0]+z[1]*z[1]+z[2]*z[2]);

								mi = 1.0/dist;
								z[0] = mi*z[0];
                                                        	z[1] = mi*z[1];
                                                        	z[2] = mi*z[2];

								g1 = k_e*core_q[k]*core_q[l]*exp(-debye*dist)*(debye*dist+1.0)/(dist*dist);
								Energy = Energy + (k_e*core_q[l]*core_q[k]/dist)*exp(-debye*dist);

								force[i1] = force[i1] + g1*z[0];
                                                                force[i2] = force[i2] + g1*z[1];
                                                                force[i3] = force[i3] + g1*z[2];
                                                                force[j1] = force[j1] - g1*z[0];
                                                                force[j2] = force[j2] - g1*z[1];
                                                                force[j3] = force[j3] - g1*z[2];

								//torque due to force
								fa = g1*(a[i1]*z[0]+a[i2]*z[1]+a[i3]*z[2]);
                                                                fb = g1*(b[i1]*z[0]+b[i2]*z[1]+b[i3]*z[2]);
                                                                fc = g1*(c[i1]*z[0]+c[i2]*z[1]+c[i3]*z[2]);
                                                                torque[i1] = torque[i1] + fc*core_pos[k2] - fb*core_pos[k3];
                                                                torque[i2] = torque[i2] + fa*core_pos[k3] - fc*core_pos[k1];
                                                                torque[i3] = torque[i3] + fb*core_pos[k1] - fa*core_pos[k2];
								fa = -g1*(a[j1]*z[0]+a[j2]*z[1]+a[j3]*z[2]);
                                                                fb = -g1*(b[j1]*z[0]+b[j2]*z[1]+b[j3]*z[2]);
                                                                fc = -g1*(c[j1]*z[0]+c[j2]*z[1]+c[j3]*z[2]);
                                                                torque[j1] = torque[j1] + fc*core_pos[l2] - fb*core_pos[l3];
                                                                torque[j2] = torque[j2] + fa*core_pos[l3] - fc*core_pos[l1];
                                                                torque[j3] = torque[j3] + fb*core_pos[l1] - fa*core_pos[l2];

								
								//Excluded Volume force
								if (dist <= 8.0){
									s1 = sigma_Core_Core;
									s2 = sigma_Core_Core;
									g1 = k_ex*((12.0/s1)*pow((s1/dist),13)-(6.0/s2)*pow((s2/dist),7));
									Energy = Energy + k_ex*(pow(s1/dist,12) - pow(s2/dist,6));

									force[i1] = force[i1] + g1*z[0];
                                                                	force[i2] = force[i2] + g1*z[1];
                                                                	force[i3] = force[i3] + g1*z[2];
                                                                	force[j1] = force[j1] - g1*z[0];
                                                                	force[j2] = force[j2] - g1*z[1];
                                                                	force[j3] = force[j3] - g1*z[2];
				
									//torque due to force
									fa = g1*(a[i1]*z[0]+a[i2]*z[1]+a[i3]*z[2]);
                                                                	fb = g1*(b[i1]*z[0]+b[i2]*z[1]+b[i3]*z[2]);
                                                                	fc = g1*(c[i1]*z[0]+c[i2]*z[1]+c[i3]*z[2]);
                                                                	torque[i1] = torque[i1] + fc*core_pos[k2] - fb*core_pos[k3];
                                                                	torque[i2] = torque[i2] + fa*core_pos[k3] - fc*core_pos[k1];
                                                                	torque[i3] = torque[i3] + fb*core_pos[k1] - fa*core_pos[k2];
                                                                	fa = -g1*(a[j1]*z[0]+a[j2]*z[1]+a[j3]*z[2]);
                                                                	fb = -g1*(b[j1]*z[0]+b[j2]*z[1]+b[j3]*z[2]);
                                                                	fc = -g1*(c[j1]*z[0]+c[j2]*z[1]+c[j3]*z[2]);
                                                                	torque[j1] = torque[j1] + fc*core_pos[l2] - fb*core_pos[l3];
                                                                	torque[j2] = torque[j2] + fa*core_pos[l3] - fc*core_pos[l1];
 	                                                              	torque[j3] = torque[j3] + fb*core_pos[l1] - fa*core_pos[l2];
								}

							}
						}

					}


				}


			}



		}

	
	}
/*
	l_h1 = 3.0;

	for (i = 1; i < n; i++){
		if (i%comm_size != rank) continue;

                if (type[i]==1){
			z[0] = r[(i-2)*3]-r[(i+2)*3];
			z[1] = r[(i-2)*3+1]-r[(i+2)*3+1];
			z[2] = r[(i-2)*3+2]-r[(i+2)*3+2];
                        dist = sqrt((r[(i-2)*3]-r[(i+2)*3])*(r[(i-2)*3]-r[(i+2)*3]) + (r[(i-2)*3+1]-r[(i+2)*3+1])*(r[(i-2)*3+1]-r[(i+2)*3+1]) + (r[(i-2)*3+2]-r[(i+2)*3+2])*(r[(i-2)*3+2]-r[(i+2)*3+2]));
                        mi = 1.0/dist;
			z[0] = mi*z[0];
			z[1] = mi*z[1];
			z[2] = mi*z[2];
			mag = 2.0*k_h1*(dist - l_h1);
			force[(i-2)*3]=  force[(i-2)*3] - mag*z[0];
			force[(i-2)*3+1]=  force[(i-2)*3+1] - mag*z[1];
			force[(i-2)*3+2]=  force[(i-2)*3+2] - mag*z[2];
			force[(i+2)*3]=  force[(i+2)*3] + mag*z[0];
			force[(i+2)*3+1]=  force[(i+2)*3+1] + mag*z[1];
			force[(i+2)*3+2]=  force[(i+2)*3+2] + mag*z[2];
                }
        }
*/
}

void rotate(int n, int n3, double* a, double* b, double* c, double* a_n, double* b_n, double* c_n, double* d_theta, double dt){

	double theta_o, si, co;
	double temp[3];
	double wa, wb, wc, g1, g2, g3;
	double z, z2, wa2, wb2, wc2, czt, Omczt, szt;
	int i, i1, i2, i3;

	for (i = 0; i < n; i++){
		i1 = 3*i;
		i2 = i1+1;
		i3 = i2+1;
		wa = d_theta[i1];
		wb = d_theta[i2];
		wc = d_theta[i3];
		wa2 = wa*wa;
		wb2 = wb*wb;
		wc2 = wc*wc;
		z2 = wa2 + wb2 + wc2;
		z = sqrt(z2);
		czt = cos(z*dt);
		szt = sin(z*dt);
		Omczt = 1.0 - czt;

		if (z2 > 0.0){
			//rotation of a
			g1 = ((wb2+wc2)*czt+wa2)/z2;
			g2 = wa*wb*Omczt/z2 + wc*szt/z;
			g3 = wa*wc*Omczt/z2 - wb*szt/z;
			a_n[i1] = g1*a[i1] + g2*b[i1] + g3*c[i1];
			a_n[i2] = g1*a[i2] + g2*b[i2] + g3*c[i2];
			a_n[i3] = g1*a[i3] + g2*b[i3] + g3*c[i3];

			//rotation of b
			g1 = wa*wb*Omczt/z2 - wc*szt/z;
                        g2 = ((wa2+wc2)*czt+wb2)/z2;
                        g3 = wb*wc*Omczt/z2 + wa*szt/z;
                        b_n[i1] = g1*a[i1] + g2*b[i1] + g3*c[i1];
                        b_n[i2] = g1*a[i2] + g2*b[i2] + g3*c[i2];
                        b_n[i3] = g1*a[i3] + g2*b[i3] + g3*c[i3];

			//rotation of c
			g1 = wa*wc*Omczt/z2 + wb*szt/z;
			g2 = wb*wc*Omczt/z2 - wa*szt/z;
                        g3 = ((wa2+wb2)*czt+wc2)/z2;
                        c_n[i1] = g1*a[i1] + g2*b[i1] + g3*c[i1];
                        c_n[i2] = g1*a[i2] + g2*b[i2] + g3*c[i2];
                        c_n[i3] = g1*a[i3] + g2*b[i3] + g3*c[i3];

		}else{
			a_n[i1] = a[i1];
			a_n[i2] = a[i2];
			a_n[i3] = a[i3];
			b_n[i1] = b[i1];
			b_n[i2] = b[i2];
			b_n[i3] = b[i3];
			c_n[i1] = c[i1];
			c_n[i2] = c[i2];
			c_n[i3] = c[i3];
		}

	}

}

void build_tail(int n, int* type, int n_t, std::vector<double> tail_pos, std::vector<double> tail_q, std::vector<int> tail_grp, std::vector<int> tail_fix, std::vector<double> tail_rad, std::vector<double> tail_bond_c, std::vector<double> tail_bond_v, std::vector<double> tail_angle_c, std::vector<double> tail_angle_v, double* r, double* a, double* b, double* c, int n_tail, int* nc_t_flag, double* r_t, double* h_t, double* g_t, double* lo_t, double* beta_o_t, double* t_q, double* t_rad, int* t_grp, int* t_fix){

	int i,j,k,k1;

	int n_grp;

	n_grp=1;
	
	for (i=1;i<n_t;i++){
		if(tail_grp[i] != tail_grp[i-1]){
			n_grp = n_grp+1;
		}
	}

	k1 = 0;
	k = 0;

	for (i=0;i<n;i++){
		if (type[i]==1){

			if (nc_t_flag[k1]==1){
				for (j=0;j<n_t;j++){
					r_t[k*n_t*3+j*3] = r[i*3] + a[i*3]*tail_pos[j*3] + b[i*3]*tail_pos[j*3+1] + c[i*3]*tail_pos[j*3+2];
					r_t[k*n_t*3+j*3+1] = r[i*3+1] + a[i*3+1]*tail_pos[j*3] + b[i*3+1]*tail_pos[j*3+1] + c[i*3+1]*tail_pos[j*3+2];
					r_t[k*n_t*3+j*3+2] = r[i*3+2] + a[i*3+2]*tail_pos[j*3] + b[i*3+2]*tail_pos[j*3+1] + c[i*3+2]*tail_pos[j*3+2];

					t_q[k*n_t+j] = tail_q[j];
					t_rad[k*n_t+j] = tail_rad[j];
					t_fix[k*n_t+j] = tail_fix[j];
					t_grp[k*n_t+j] = tail_grp[j] + k*n_grp;
					h_t[k*n_t+j] = tail_bond_c[j];
					lo_t[k*n_t+j] = tail_bond_v[j];
					g_t[k*n_t+j] = tail_angle_c[j];
					beta_o_t[k*n_t+j] = tail_angle_v[j];
				}
				k = k+1;
			}
			k1 = k1+1;
		}
	}

}

void tail_local(int n, int* type, int n_t, double* r, double* a, double* b, double* c, int* nc_t_flag, int n_tail, double* r_t, double* r_t_local){

	int i,j,k,k1;

	k=0;
	k1=0;

	for (i=0;i<n;i++){
                if (type[i]==1){
                        if (nc_t_flag[k1]==1){
                                for (j=0;j<n_t;j++){
					r_t_local[k*n_t*3+j*3] = (r_t[k*n_t*3+j*3] - r[i*3])*a[i*3] + (r_t[k*n_t*3+j*3+1] - r[i*3+1])*a[i*3+1] + (r_t[k*n_t*3+j*3+2] - r[i*3+2])*a[i*3+2];
					r_t_local[k*n_t*3+j*3+1] = (r_t[k*n_t*3+j*3] - r[i*3])*b[i*3] + (r_t[k*n_t*3+j*3+1] - r[i*3+1])*b[i*3+1] + (r_t[k*n_t*3+j*3+2] - r[i*3+2])*b[i*3+2];
					r_t_local[k*n_t*3+j*3+2] = (r_t[k*n_t*3+j*3] - r[i*3])*c[i*3] + (r_t[k*n_t*3+j*3+1] - r[i*3+1])*c[i*3+1] + (r_t[k*n_t*3+j*3+2] - r[i*3+2])*c[i*3+2];
                                }
                                k = k+1;
                        }
                        k1 = k1+1;
                }
        }

}

void tail_grow(int n, int* type, int n_t, double* r, double* a, double* b, double* c, int* nc_t_flag, int n_tail, double* r_t, double* r_t_local){

	int i,j,k,k1;

        k=0;
        k1=0;

        for (i=0;i<n;i++){
                if (type[i]==1){
                        if (nc_t_flag[k1]==1){
                                for (j=0;j<n_t;j++){
					r_t[k*n_t*3+j*3] = r[i*3] + a[i*3]*r_t_local[k*n_t*3+j*3] + b[i*3]*r_t_local[k*n_t*3+j*3+1] + c[i*3]*r_t_local[k*n_t*3+j*3+2];
                                        r_t[k*n_t*3+j*3+1] = r[i*3+1] + a[i*3+1]*r_t_local[k*n_t*3+j*3] + b[i*3+1]*r_t_local[k*n_t*3+j*3+1] + c[i*3+1]*r_t_local[k*n_t*3+j*3+2];
                                        r_t[k*n_t*3+j*3+2] = r[i*3+2] + a[i*3+2]*r_t_local[k*n_t*3+j*3] + b[i*3+2]*r_t_local[k*n_t*3+j*3+1] + c[i*3+2]*r_t_local[k*n_t*3+j*3+2];
                                }
                                k = k+1;
                        }
                        k1 = k1+1;
                }
        }	

}



void update_tail_beta(int n_tail, double* r_t, double* beta_t, int* t_grp){

	double x12, y12, z12, r12, x32, y32, z32, r32, p123;

	for (int i = 0; i < n_tail; i++){
		if((t_grp[i] == t_grp[i+1]) and (t_grp[i] == t_grp[i+2]) and (i < n_tail-2)){
			x12 = r_t[i*3] - r_t[i*3+3];
			y12 = r_t[i*3+1] - r_t[i*3+4];
			z12 = r_t[i*3+2] - r_t[i*3+5];
	
			r12 = sqrt(x12*x12 + y12*y12 + z12*z12);

			x32 = r_t[i*3+6] - r_t[i*3+3];
			y32 = r_t[i*3+7] - r_t[i*3+4];
			z32 = r_t[i*3+8] - r_t[i*3+5];

			r32 = sqrt(x32*x32 + y32*y32 + z32*z32);

			p123 = (x12*x32 + y12*y32 + z12*z32)/(r12*r32);

			beta_t[i] = acos(p123);
		}else{
			beta_t[i] = 0.0;
		}

	}

}


void tail_force(int n, int n_t, int* type, int n_tail, int n_tail3, std::vector<double> tail_pos, std::vector<int> tail_fix, int* nc_t_flag, double* r_t, double* beta_t, double* h_t, double* g_t, double* lo_t, double* beta_o_t, double h, double k_e, double debye, double k_ex, double* t_q, double* t_rad, int* t_grp, int* t_fix, double* r, double* a, double* b, double* c, double q_l, int Nq, int Nq3, std::vector<double> core_pos, std::vector<double> core_q, double* t_force, double* force, double* torque, double sigma_Tail_Tail, double sigma_Tail_Linker, double sigma_Tail_Core, double& Energy, int rank, int comm_size){

	double Stri[3], Strim1[3], Ai[3], Bi[3], Aim1[3], Bim1[3];
	double distance, fa, fb, fc;
	double a_f[3], a_b[3], a_o[3];
	double beta_b, dis;
	double distancef, distanceb;
	int i, j, k, count, n_ct_cnt, k1, index;
	double r_to[3], z[3];
	double g1, g2, c1, c2;
	double torque1[3], torque2[3], torque3[3];
	double force_projection1[3], force_projection2[3],force_projection3[3];
	double torque_tc[3];

	double r_t_tmp1[3], r_t_tmp2[3], r_t_tmp3[3], r_tmp[3], core_pos_tmp[3], tail_pos_tmp[3], a_tmp[3], b_tmp[3], c_tmp[3], h_stri_tmp[3];

	//Stretching between tail beads

	for (i=0; i<n_tail-1; i++){
		if (i%comm_size != rank) continue;
		if (t_grp[i]==t_grp[i+1]){
			for (int ix = 0; ix<3; ix++){
				r_t_tmp1[ix] = r_t[i*3+ix];
				r_t_tmp2[ix] = r_t[i*3+3+ix];
			}
			Force_Stretching(h_t[i], r_t_tmp1, r_t_tmp2, lo_t[i], force_projection1, force_projection2, Energy);

			t_force[i*3] = t_force[i*3] + force_projection1[0];
                	t_force[i*3+1] = t_force[i*3+1] + force_projection1[1];
                	t_force[i*3+2] = t_force[i*3+2] + force_projection1[2];
                	t_force[i*3+3] = t_force[i*3+3] + force_projection2[0];
                	t_force[i*3+4] = t_force[i*3+4] + force_projection2[1];
                	t_force[i*3+5] = t_force[i*3+5] + force_projection2[2];
		}
	}


	//Bending between tail beads

	for (i=0; i<n_tail-2;i++){
		if (i%comm_size != rank) continue;
		if ((t_grp[i]==t_grp[i+1]) and (t_grp[i]==t_grp[i+2])){
			for (int ix = 0; ix<3; ix++){
                                r_t_tmp1[ix] = r_t[i*3+ix];
                                r_t_tmp2[ix] = r_t[i*3+3+ix];
				r_t_tmp3[ix] = r_t[i*3+6+ix];
                        }
			Force_Bending(g_t[i],beta_t[i],beta_o_t[i], r_t_tmp1, r_t_tmp2, r_t_tmp3, force_projection1, force_projection2, force_projection3, Energy);

	//		printf("host %d %f %f %f \n", i, g_t[i], beta_t[i], beta_o_t[i]);
//			printf("host %d %f %f %f \n", i, force_projection1[0], force_projection1[1], force_projection1[2]);

			t_force[i*3] = t_force[i*3] - force_projection1[0];
			t_force[i*3+1] = t_force[i*3+1] - force_projection1[1];
			t_force[i*3+2] = t_force[i*3+2] - force_projection1[2];

			t_force[i*3+3] = t_force[i*3+3] - force_projection2[0];
			t_force[i*3+4] = t_force[i*3+4] - force_projection2[1];
			t_force[i*3+5] = t_force[i*3+5] - force_projection2[2];

			t_force[i*3+6] = t_force[i*3+6] - force_projection3[0];
			t_force[i*3+7] = t_force[i*3+7] - force_projection3[1];
			t_force[i*3+8] = t_force[i*3+8] - force_projection3[2];

		}
	}


	//Electrostatics and Excluded volume between tail beads

	for (i=0;i<n_tail-1;i++){
		for (j=i+1;j<n_tail;j++){
//			printf("host %d %d \n", i, j);
			if ((i+j)%comm_size != rank) continue;
			if ((t_grp[i]!=t_grp[j]) or (j>i+2)){
				for (int ix = 0; ix<3; ix++){
                                	r_t_tmp1[ix] = r_t[i*3+ix];
                                	r_t_tmp2[ix] = r_t[j*3+ix];
                        	}
//				Force_Ele_Vdw(t_q[i], t_q[j], 1/(4*PI*k_e), debye, k_ex, sigma_Tail_Tail, r_t_tmp1, r_t_tmp2, force_projection1, force_projection2, Energy);

				Force_Electrostatics(t_q[i], t_q[j], 1/(4*PI*k_e), debye, r_t_tmp1, r_t_tmp2, force_projection1, force_projection2, Energy);
/*
				if(i==132){
				printf("host %d %f %f %f %f %f %f \n", j, force_projection1[0], force_projection1[1], force_projection1[2], force_projection2[0], force_projection2[1], force_projection2[2]);
				}
*/
				
				t_force[i*3] = t_force[i*3] + force_projection1[0];
	                        t_force[i*3+1] = t_force[i*3+1] + force_projection1[1];
                        	t_force[i*3+2] = t_force[i*3+2] + force_projection1[2];

				t_force[j*3] = t_force[j*3] + force_projection2[0];
                                t_force[j*3+1] = t_force[j*3+1] + force_projection2[1];
                                t_force[j*3+2] = t_force[j*3+2] + force_projection2[2];

			}
		}
	}




	//Electrostatics and Excluded volume between tail and DNA beads

	for (i=0;i<n_tail;i++){
		for (j=0;j<n;j++){
			if ((i+j)%comm_size != rank) continue;
			if(type[j]==0){
			        for (int ix = 0; ix<3; ix++){
                                        r_t_tmp1[ix] = r_t[i*3+ix];
                                        r_tmp[ix] = r[j*3+ix];
                                }
				Force_Ele_Vdw(t_q[i], q_l, 1/(4*PI*k_e), debye,k_ex, sigma_Tail_Linker, r_t_tmp1, r_tmp, force_projection1, force_projection2, Energy);

				t_force[i*3] = t_force[i*3] + force_projection1[0];
                                t_force[i*3+1] = t_force[i*3+1] + force_projection1[1];
                                t_force[i*3+2] = t_force[i*3+2] + force_projection1[2];

				force[j*3] = force[j*3] + force_projection2[0];
                                force[j*3+1] = force[j*3+1] + force_projection2[1];
                                force[j*3+2] = force[j*3+2] + force_projection2[2];
			}
		}
	}

	//Electrostatics and Excluded volume between tail and core
/*
	count=0;

	n_ct_cnt=0;

	for (i=0;i<n;i++){
		if (type[i]==1){
			if (nc_t_flag[count]==1){
				n_ct_cnt=n_ct_cnt+1;
			}
			count = count+1;
			for (j=0;j<Nq;j++){
				for (k=0;k<n_tail;k++){
					if ((i+j+k)%comm_size != rank) continue;
					z[0] = r[i*3]+a[i*3]*core_pos[j*3]+b[i*3]*core_pos[j*3+1]+c[i*3]*core_pos[j*3+2];
                                        z[1] = r[i*3+1]+a[i*3+1]*core_pos[j*3]+b[i*3+1]*core_pos[j*3+1]+c[i*3+1]*core_pos[j*3+2];
                                        z[2] = r[i*3+2]+a[i*3+2]*core_pos[j*3]+b[i*3+2]*core_pos[j*3+1]+c[i*3+2]*core_pos[j*3+2];

					if (k/n_t != n_ct_cnt-1){
						for (int ix = 0; ix<3; ix++){
                                        		r_t_tmp1[ix] = r_t[k*3+ix];
							core_pos_tmp[ix] =  core_pos[j*3+ix];
							a_tmp[ix] = a[i*3+ix];
							b_tmp[ix] = b[i*3+ix];
							c_tmp[ix] = c[i*3+ix];
                                		}	
						Force_Ele_Vdw(core_q[j], t_q[k], 1/(4*PI*k_e), debye,k_ex, sigma_Tail_Core, z, r_t_tmp1,force_projection1, force_projection2, Energy);

//						Force_Electrostatics(core_q[j], t_q[k], 1/(4*PI*k_e), debye, z, r_t_tmp1, force_projection1, force_projection2, Energy);

						torque_due_to_force_relative(force_projection1, core_pos_tmp, a_tmp, b_tmp, c_tmp, torque_tc);
						force[i*3] = force[i*3] + force_projection1[0];
                                		force[i*3+1] = force[i*3+1] + force_projection1[1];
	                                	force[i*3+2] = force[i*3+2] + force_projection1[2];

						t_force[k*3] = t_force[k*3] + force_projection2[0];
                	                	t_force[k*3+1] = t_force[k*3+1] + force_projection2[1];
                        	        	t_force[k*3+2] = t_force[k*3+2] + force_projection2[2];
						torque[i*3] = torque[i*3] + torque_tc[0];
						torque[i*3+1] = torque[i*3+1] + torque_tc[1];
						torque[i*3+2] = torque[i*3+2] + torque_tc[2]; 
					}
				}
			}
		}
	}

*/

	//Stretching between tail and core

	k = 0;

	for (i = 0; i < n; i++){
		if (type[i]==1){
			if(nc_t_flag[k]==1){
				for (j = 0; j < n_t; j++){
					if ((i+j)%comm_size != rank) continue;
					if (tail_fix[j]==1){
						
						r_to[0] = r[i*3] + a[i*3]*tail_pos[j*3] + b[i*3]*tail_pos[j*3+1] + c[i*3]*tail_pos[j*3+2];
						r_to[1] = r[i*3+1] + a[i*3+1]*tail_pos[j*3] + b[i*3+1]*tail_pos[j*3+1] + c[i*3+1]*tail_pos[j*3+2];
						r_to[2] = r[i*3+2] + a[i*3+2]*tail_pos[j*3] + b[i*3+2]*tail_pos[j*3+1] + c[i*3+2]*tail_pos[j*3+2];

						Stri[0] = r_t[k*n_t*3+j*3] - r_to[0];
						Stri[1] = r_t[k*n_t*3+j*3+1] - r_to[1];
						Stri[2] = r_t[k*n_t*3+j*3+2] - r_to[2];
			
						Energy = Energy + h*(Stri[0]*Stri[0]+Stri[1]*Stri[1]+Stri[2]*Stri[2])/2;

//						printf("host %d %d %d %f %f %f %f \n", i, j, k, h, r_to[0], r_to[1], r_to[2]);

//						printf("host %d %d %d %f %f %f %f %f %f %f %f %f %f \n", i, j, k, h, r_to[0], r_to[1], r_to[2], r_t[k*n_t*3+j*3], r_t[k*n_t*3+j*3+1], r_t[k*n_t*3+j*3+2], Stri[0], Stri[1], Stri[2]);

						t_force[k*n_t*3+j*3] = t_force[k*n_t*3+j*3] - h*Stri[0];
						t_force[k*n_t*3+j*3+1] = t_force[k*n_t*3+j*3+1] - h*Stri[1];
						t_force[k*n_t*3+j*3+2] = t_force[k*n_t*3+j*3+2] - h*Stri[2];

						force[i*3] = force[i*3] + h*Stri[0];
						force[i*3+1] = force[i*3+1] + h*Stri[1];
						force[i*3+2] = force[i*3+2] + h*Stri[2];

						// torque due to force

						for (int ix = 0; ix<3; ix++){
							tail_pos_tmp[ix] =  tail_pos[j*3+ix];
							a_tmp[ix] = a[i*3+ix];
							b_tmp[ix] = b[i*3+ix];
							c_tmp[ix] = c[i*3+ix];
							h_stri_tmp[ix] = h*Stri[ix];
						}

						torque_due_to_force_relative(h_stri_tmp, tail_pos_tmp, a_tmp, b_tmp, c_tmp, torque_tc);
						torque[i*3] = torque[i*3] + torque_tc[0];
						torque[i*3+1] = torque[i*3+1] + torque_tc[1];
						torque[i*3+2] = torque[i*3+2] + torque_tc[2];

					}
				}
			}
			k = k+1;
		}
	}	



}


void build_LH(int n, int* type, int n_LH, std::vector<double> LH_n_pos, std::vector<double> LH_c_pos, std::vector<double> LH_radius, double* r, double* a, double* b, double* c, int n_lh_n, int n_lh_c, int* nc_lh_flag, double* r_lh, double* lh_rad){

        int i,j,k,k1;

	int nlh = n_lh_n + n_lh_c;

	k1 = 0;
        k = 0;	

        for (i=0;i<n;i++){
                if (type[i]==1){
			if (nc_lh_flag[k1]==1){
        	                for (j=0;j<nlh;j++){
					if (j < n_lh_n){
						r_lh[k*nlh*3+j*3] = r[i*3] + a[i*3]*LH_n_pos[j*3] + b[i*3]*LH_n_pos[j*3+1] + c[i*3]*LH_n_pos[j*3+2];
                                	        r_lh[k*nlh*3+j*3+1] = r[i*3+1] + a[i*3+1]*LH_n_pos[j*3] + b[i*3+1]*LH_n_pos[j*3+1] + c[i*3+1]*LH_n_pos[j*3+2];
                        	                r_lh[k*nlh*3+j*3+2] = r[i*3+2] + a[i*3+2]*LH_n_pos[j*3] + b[i*3+2]*LH_n_pos[j*3+1] + c[i*3+2]*LH_n_pos[j*3+2];
					}else{
        	                        	r_lh[k*nlh*3+j*3] = r[i*3] + a[i*3]*LH_c_pos[j*3] + b[i*3]*LH_c_pos[j*3+1] + c[i*3]*LH_c_pos[j*3+2];
	                                	r_lh[k*nlh*3+j*3+1] = r[i*3+1] + a[i*3+1]*LH_c_pos[j*3] + b[i*3+1]*LH_c_pos[j*3+1] + c[i*3+1]*LH_c_pos[j*3+2];
                                		r_lh[k*nlh*3+j*3+2] = r[i*3+2] + a[i*3+2]*LH_c_pos[j*3] + b[i*3+2]*LH_c_pos[j*3+1] + c[i*3+2]*LH_c_pos[j*3+2];
					}
                                	lh_rad[k*nlh+j] = LH_radius[j];
                        	}
				k = k+1;
			}
                        k1 = k1+1;
                }
        }

}


void update_LH_beta(int n_LH, int n_lh_c, double* r_lh, double* beta_lh){

        double x12, y12, z12, r12, x32, y32, z32, r32, p123;

        for (int i = 0; i < n_LH; i++){
                if((i%n_lh_c!=0) and (i%n_lh_c != n_lh_c-1)){
                        x12 = r_lh[i*3-3] - r_lh[i*3];
                        y12 = r_lh[i*3-2] - r_lh[i*3+1];
                        z12 = r_lh[i*3-1] - r_lh[i*3+2];

                        r12 = sqrt(x12*x12 + y12*y12 + z12*z12);

                        x32 = r_lh[i*3+3] - r_lh[i*3];
                        y32 = r_lh[i*3+4] - r_lh[i*3+1];
                        z32 = r_lh[i*3+5] - r_lh[i*3+2];

                        r32 = sqrt(x32*x32 + y32*y32 + z32*z32);

                        p123 = (x12*x32 + y12*y32 + z12*z32)/(r12*r32);

                        beta_lh[i] = acos(p123);
                }else{
                        beta_lh[i] = 0.0;
                }

        }

}



void Linker_Histone_force(int n, int n_tail, int n_lh_n, int n_lh_g, int n_lh_c, int* t_grp, int* t_fix, int* type, int n_LH, int n_LH3, std::vector<double> LH_g_pos, std::vector<int> LH_conn, int* nc_lh_flag, double* beta_lh, double* r_lh, std::vector<double> LH_q, double k_e, double debye, double k_ex, std::vector<double> LH_vdw_hh, std::vector<double> LH_vdw_hc, std::vector<double> LH_vdw_hl, std::vector<double> LH_vdw_ht, std::vector<double> LH_kstr, std::vector<double> LH_kben, std::vector<double> LH_streq, std::vector<double> LH_betaeq, double* r, double* a, double* b, double* c, double q_l, int Nq, int Nq3, std::vector<double> core_pos, std::vector<double> core_q, double* r_t, double* t_q, double* t_force, double* LH_force, double* force, double* torque, double& Energy, int rank, int comm_size){

	double Stri[3], Strim1[3];
        double distance, fa, fb, fc;
        int i, j, k, k2, m, m1, xi, count, n_clh_cnt;
        double r_lho[3], r_lho2[3], z[3];
	int index, index1, index2, index_c;
	int nlh = n_lh_n+n_lh_c;
	double torque_hgl[3], torque_lhc[3], f_p[3], torque1[3], torque2[3];
	double force_projection1[3], force_projection2[3],force_projection3[3];

	double r_lh_tmp1[3], r_lh_tmp2[3], r_lh_tmp3[3], r_t_tmp[3], r_tmp[3], a_tmp[3], b_tmp[3], c_tmp[3], core_pos_tmp[3], LH_g_pos_tmp[3];

	ofstream wf;
	wf.open("wf.txt", std::ios_base::app);

	//Stretching force between LH beads

	for (i=0; i<n_LH; i++){
		if (i%comm_size != rank) continue;

		index = i%(n_lh_n+n_lh_c);
                if (index < n_lh_n){
                        index = index;
                }else{
                        index = index + n_lh_g;
                }

		if (LH_conn[index]==1){
			for (xi = 0; xi<3; xi++){
				r_lh_tmp1[xi] = r_lh[i*3+xi];
				r_lh_tmp2[xi] = r_lh[i*3+3+xi];
			}

			Force_Stretching(LH_kstr[index], r_lh_tmp1, r_lh_tmp2, LH_streq[index],force_projection1, force_projection2, Energy);
//			printf("host %d %d %f %f %f %f \n", i, index, LH_kstr[index], r_lh_tmp1[0], r_lh_tmp1[1], r_lh_tmp1[2]);

			LH_force[i*3] = LH_force[i*3] + force_projection1[0];
                        LH_force[i*3+1] = LH_force[i*3+1] + force_projection1[1];
                        LH_force[i*3+2] = LH_force[i*3+2] + force_projection1[2];
 
 			LH_force[i*3+3] = LH_force[i*3+3] + force_projection2[0];
                        LH_force[i*3+4] = LH_force[i*3+4] + force_projection2[1];
                        LH_force[i*3+5] = LH_force[i*3+5] + force_projection2[2];

			if (abs(force_projection1[0])>1 || abs(force_projection1[1])>1 || abs(force_projection1[2])>1 || abs(force_projection2[0])>1 || abs(force_projection2[1])>1 || abs(force_projection2[2])>1){
				wf << "stretching " << endl;
			}

		}
	}



	//Bending between LH beads

	for (i=0; i<n_LH;i++){
		if (i%comm_size != rank) continue;

		index = i%(n_lh_n+n_lh_c);
                if (index < n_lh_n){
                        index = index;
                }else{
                        index = index + n_lh_g;
                }
                if ((i%n_lh_c!=0) and (i%n_lh_c!=n_lh_c-1)){
			for (xi = 0; xi<3; xi++){
                                r_lh_tmp1[xi] = r_lh[i*3-3+xi];
                                r_lh_tmp2[xi] = r_lh[i*3+xi];
				r_lh_tmp3[xi] = r_lh[i*3+3+xi];
                        }
			Force_Bending(LH_kben[index],beta_lh[i],LH_betaeq[index],r_lh_tmp1, r_lh_tmp2, r_lh_tmp3,force_projection1,force_projection2,force_projection3, Energy);
                        LH_force[i*3-3] = LH_force[i*3-3] - force_projection1[0];
                        LH_force[i*3-2] = LH_force[i*3-2] - force_projection1[1];
                        LH_force[i*3-1] = LH_force[i*3-1] - force_projection1[2];
                        LH_force[i*3] = LH_force[i*3] - force_projection2[0];
                        LH_force[i*3+1] = LH_force[i*3+1] - force_projection2[1];
                        LH_force[i*3+2] = LH_force[i*3+2] - force_projection2[2];
                        LH_force[i*3+3] = LH_force[i*3+3] - force_projection3[0];
                        LH_force[i*3+4] = LH_force[i*3+4] - force_projection3[1];
                        LH_force[i*3+5] = LH_force[i*3+5] - force_projection3[2];

			if (abs(force_projection1[0])>1 || abs(force_projection1[1])>1 || abs(force_projection1[2])>1 || abs(force_projection2[0])>1 || abs(force_projection2[1])>1 || abs(force_projection2[2])>1 || abs(force_projection3[0])>1 || abs(force_projection3[1])>1 || abs(force_projection3[2])>1){
                                wf << "bending " << endl;
                        }


                }
        }






	//Electrostatic and Excluded Volume between LH	
	
	for (i=0;i<n_LH;i++){
		index1 = i%(n_lh_n+n_lh_c);
                if (index1 < n_lh_n){
                        index1 = index1;
                }else{
                        index1 = index1 + n_lh_g;
                }
		for (j=i;j<n_LH;j++){
//			printf("host %d %d \n", i, j);
			if ((i+j)%comm_size != rank) continue;
			index2 = j%(n_lh_n+n_lh_c);
                	if (index2 < n_lh_n){
                        	index2 = index2;
                	}else{
				index2 = index2 + n_lh_g;
                        }
			
			if ((i/(n_lh_n+n_lh_c) == j/(n_lh_n+n_lh_c)) and (j > i+1)){

				for (xi = 0; xi<3; xi++){
                                        r_lh_tmp1[xi] = r_lh[i*3+xi];
                                        r_lh_tmp2[xi] = r_lh[j*3+xi];
                                }

				Force_Exclude_Volume(k_ex, (LH_vdw_hh[index1]+LH_vdw_hh[index2])/2, r_lh_tmp1, r_lh_tmp2,force_projection1, force_projection2, Energy);
/*				if(i==0){
                                printf("host %d %f %f %f %f %f %f \n", j, force_projection1[0], force_projection1[1], force_projection1[2], force_projection2[0], force_projection2[1], force_projection2[2]);
                                }
*/
				LH_force[i*3] = LH_force[i*3] + force_projection1[0];
	                        LH_force[i*3+1] = LH_force[i*3+1] + force_projection1[1];
        	                LH_force[i*3+2] = LH_force[i*3+2] + force_projection1[2];
				LH_force[j*3] = LH_force[j*3] + force_projection2[0];
                	        LH_force[j*3+1] = LH_force[j*3+1] + force_projection2[1];
                        	LH_force[j*3+2] = LH_force[j*3+2] + force_projection2[2];

			}else if ((i/(n_lh_n+n_lh_c) != j/(n_lh_n+n_lh_c))){
				for (xi = 0; xi<3; xi++){
                                	r_lh_tmp1[xi] = r_lh[i*3+xi];
                                	r_lh_tmp2[xi] = r_lh[j*3+xi];
                        	}


				Force_Ele_Vdw(LH_q[index1], LH_q[index2], 1/(4*PI*k_e), debye,k_ex, (LH_vdw_hh[index1]+LH_vdw_hh[index2])/2, r_lh_tmp1, r_lh_tmp2,force_projection1, force_projection2, Energy);
				LH_force[i*3] = LH_force[i*3] + force_projection1[0]; 
                                LH_force[i*3+1] = LH_force[i*3+1] + force_projection1[1];
                                LH_force[i*3+2] = LH_force[i*3+2] + force_projection1[2];
                                LH_force[j*3] = LH_force[j*3] + force_projection2[0];
                                LH_force[j*3+1] = LH_force[j*3+1] + force_projection2[1];
                                LH_force[j*3+2] = LH_force[j*3+2] + force_projection2[2];

			}

		}
	}


	//Electrostatic and Excluded Volume between LH and DNA

	for (i=0;i<n_LH;i++){
		index = i%(n_lh_n+n_lh_c);
                if (index < n_lh_n){
                        index = index;
                }else{
                        index = index + n_lh_g;
                }
		for (j=0;j<n;j++){
			if ((i+j)%comm_size != rank) continue;
			if (type[j]==0){
				for (xi = 0; xi<3; xi++){
                                	r_lh_tmp1[xi] = r_lh[i*3+xi];
                                	r_tmp[xi] = r[j*3+xi];
                        	}
				Force_Ele_Vdw(LH_q[index], q_l, 1/(4*PI*k_e), debye,k_ex, LH_vdw_hl[index], r_lh_tmp1, r_tmp,force_projection1, force_projection2, Energy);
				LH_force[i*3] = LH_force[i*3] + force_projection1[0];
                                LH_force[i*3+1] = LH_force[i*3+1] + force_projection1[1];
                                LH_force[i*3+2] = LH_force[i*3+2] + force_projection1[2];

				force[j*3] = force[j*3] + force_projection2[0];
                                force[j*3+1] = force[j*3+1] + force_projection2[1];
                                force[j*3+2] = force[j*3+2] + force_projection2[2];

				if (abs(force_projection1[0])>1 || abs(force_projection1[1])>1 || abs(force_projection1[2])>1 || abs(force_projection2[0])>1 || abs(force_projection2[1])>1 || abs(force_projection2[2])>1){
                                	wf << "LH-DNA " << endl;
                        	}


			}
		}
	}


	//Electrostatic and Excluded Volume between LH and core
/*

	for (i=0;i<n_LH;i++){
		index = i%(n_lh_n+n_lh_c);
                if (index < n_lh_n){
                        index = index;
                }else{
                        index = index + n_lh_g;
                }
		count = 0;
		k2 = 0;
		for (j=0; j<n; j++){
			if (type[j]==1){
				if (i/(n_lh_n+n_lh_c)!=k2){
					for (k=0;k<Nq;k++){
						z[0] = r[j*3]+a[j*3]*core_pos[k*3]+b[j*3]*core_pos[k*3+1]+c[j*3]*core_pos[k*3+2];
						z[1] = r[j*3+1]+a[j*3+1]*core_pos[k*3]+b[j*3+1]*core_pos[k*3+1]+c[j*3+1]*core_pos[k*3+2];
						z[2] = r[j*3+2]+a[j*3+2]*core_pos[k*3]+b[j*3+2]*core_pos[k*3+1]+c[j*3+2]*core_pos[k*3+2];
						for (xi = 0; xi<3; xi++){
							r_lh_tmp1[xi] = r_lh[i*3+xi];
							core_pos_tmp[xi] = core_pos[k*3+xi];
							a_tmp[xi] = a[j*3+xi];
							b_tmp[xi] = b[j*3+xi];
							c_tmp[xi] = c[j*3+xi];
						}

//						Force_Electrostatics(core_q[k], LH_q[index], 1/(4*PI*k_e), debye, z, r_lh_tmp1, force_projection1, force_projection2, Energy);
						Force_Ele_Vdw(core_q[k], LH_q[index], 1/(4*PI*k_e), debye,k_ex, LH_vdw_hc[index], z, r_lh_tmp1, force_projection1, force_projection2, Energy);
						torque_due_to_force_relative(force_projection1, core_pos_tmp, a_tmp, b_tmp, c_tmp,torque_lhc);

						force[j*3] = force[j*3] + force_projection1[0];
						force[j*3+1] = force[j*3+1] + force_projection1[1];
						force[j*3+2] = force[j*3+2] + force_projection1[2];
						LH_force[i*3] = LH_force[i*3] + force_projection2[0];
						LH_force[i*3+1] = LH_force[i*3+1] + force_projection2[1];
						LH_force[i*3+2] = LH_force[i*3+2] + force_projection2[2];
						torque[j*3] = torque[j*3] + torque_lhc[0];
						torque[j*3+1] = torque[j*3+1] + torque_lhc[1];
						torque[j*3+2] = torque[j*3+2] + torque_lhc[2];

						if (abs(force_projection1[0])>1 || abs(force_projection1[1])>1 || abs(force_projection1[2])>1 || abs(force_projection2[0])>1 || abs(force_projection2[1])>1 || abs(force_projection2[2])>1){
                                			wf << "LH-core " << endl;
                        			}

					}
				}
				if (nc_lh_flag[count]==1){
					k2=k2+1;
				}
	
				count=count+1;
			}
		}

	}
*/

	
	//Electrostatic and Excluded Volume between LH and Tails

	for (i=0;i<n_LH;i++){
		index = i%(n_lh_n+n_lh_c);
                if (index < n_lh_n){
                        index = index;
                }else{
                        index = index + n_lh_g;
                }
		for (j=0;j<n_tail;j++){
			if ((i+j)%comm_size != rank) continue;
			for (xi = 0; xi<3; xi++){
                                r_lh_tmp1[xi] = r_lh[i*3+xi];
                                r_t_tmp[xi] = r_t[j*3+xi];
                        }
			Force_Ele_Vdw(LH_q[index], t_q[j], 1/(4*PI*k_e), debye,k_ex, LH_vdw_ht[index], r_lh_tmp1, r_t_tmp,force_projection1, force_projection2, Energy);
			LH_force[i*3] = LH_force[i*3] + force_projection1[0];
                        LH_force[i*3+1] = LH_force[i*3+1] + force_projection1[1];
                        LH_force[i*3+2] = LH_force[i*3+2] + force_projection1[2];

			if (t_fix[j]==0){
				t_force[j*3] = t_force[j*3] + force_projection2[0];
        	                t_force[j*3+1] = t_force[j*3+1] + force_projection2[1];
                	        t_force[j*3+2] = t_force[j*3+2] + force_projection2[2];
			}else{
                                k=0;
                                for (int ccnt=0; ccnt<n; ccnt++){
                                        if (k==(t_grp[j]-1)/10 and type[ccnt]==1){
                                                index_c = ccnt;
                                                break;
                                        }
                                        if(type[ccnt]==1){
                                                k=k+1;
                                        }
                                }

                                for (int ix = 0; ix<3; ix++){
                                        r_tmp[ix] = r[index_c*3+ix];
                                        a_tmp[ix] = a[index_c*3+ix];
                                        b_tmp[ix] = b[index_c*3+ix];
                                        c_tmp[ix] = c[index_c*3+ix];
                                }
                                torque_due_to_force(force_projection2, r_t_tmp, r_tmp, a_tmp, b_tmp, c_tmp, torque2);
//				printf("host %d %d %d %f %f %f %f %f %f \n", i, j, index_c, force_projection2[0], force_projection2[1], force_projection2[2], torque2[0], torque2[1], torque2[2]);
				force[index_c*3] = force[index_c*3] + force_projection2[0];
                                force[index_c*3+1] = force[index_c*3+1] + force_projection2[1];
                                force[index_c*3+2] = force[index_c*3+2] + force_projection2[2];
                                torque[index_c*3] = torque[index_c*3] + torque2[0];
                                torque[index_c*3+1] = torque[index_c*3+1] + torque2[1];
                                torque[index_c*3+2] = torque[index_c*3+2] + torque2[2];
                        }

			if (abs(force_projection1[0])>1 || abs(force_projection1[1])>1 || abs(force_projection1[2])>1 || abs(force_projection2[0])>1 || abs(force_projection2[1])>1 || abs(force_projection2[2])>1){
                                wf << "LH-tail " << endl;
                        }

		}
	}


	//Electrostatic and Excluded Volume between LH global head and Tails

	k2=0;

	for (i=0;i<n;i++){
		if (type[i]==1){ 
			if (nc_lh_flag[k2]==1){
				for (m=0;m<n_lh_g;m++){
					r_lho[0] = r[i*3] + a[i*3]*LH_g_pos[m*3] + b[i*3]*LH_g_pos[m*3+1] + c[i*3]*LH_g_pos[m*3+2];
					r_lho[1] = r[i*3+1] + a[i*3+1]*LH_g_pos[m*3] + b[i*3+1]*LH_g_pos[m*3+1] + c[i*3+1]*LH_g_pos[m*3+2];
					r_lho[2] = r[i*3+2] + a[i*3+2]*LH_g_pos[m*3] + b[i*3+2]*LH_g_pos[m*3+1] + c[i*3+2]*LH_g_pos[m*3+2];

					for (j=0;j<n_tail;j++){
						if ((i+j)%comm_size != rank) continue;
						for (xi = 0; xi<3; xi++){
							r_t_tmp[xi] = r_t[j*3+xi];
							LH_g_pos_tmp[xi] = LH_g_pos[m*3+xi];
							a_tmp[xi] = a[i*3+xi];
							b_tmp[xi] = b[i*3+xi];
							c_tmp[xi] = c[i*3+xi];
						}
						Force_Ele_Vdw(LH_q[m], t_q[j], 1/(4*PI*k_e), debye,k_ex, LH_vdw_ht[m], r_lho, r_t_tmp,force_projection1, force_projection2, Energy);

//						printf("host %d %d %d %f %f %f %f %f %f \n", i, m, j, force_projection1[0], force_projection1[1], force_projection1[2], LH_g_pos[m*3], LH_g_pos[m*3+1], LH_g_pos[m*3+2]);

						force[i*3] = force[i*3] + force_projection1[0];
						force[i*3+1] = force[i*3+1] + force_projection1[1];
						force[i*3+2] = force[i*3+2] + force_projection1[2];

						torque_due_to_force_relative(force_projection1, LH_g_pos_tmp, a_tmp, b_tmp, c_tmp, torque_hgl);
						torque[i*3] = torque[i*3] + torque_hgl[0];
						torque[i*3+1] = torque[i*3+1] + torque_hgl[1];
						torque[i*3+2] = torque[i*3+2] + torque_hgl[2];


						if (t_fix[j]==0){
							t_force[j*3] = t_force[j*3] + force_projection2[0];
							t_force[j*3+1] = t_force[j*3+1] + force_projection2[1];
							t_force[j*3+2] = t_force[j*3+2] + force_projection2[2];
						}else{
							k=0;
							for (int ccnt=0; ccnt<n; ccnt++){
								if (k==(t_grp[j]-1)/10 and type[ccnt]==1){
									index_c = ccnt;
									break;
								}
								if(type[ccnt]==1){
									k=k+1;
								}
							}

							for (int ix = 0; ix<3; ix++){
								r_tmp[ix] = r[index_c*3+ix];
								a_tmp[ix] = a[index_c*3+ix];
								b_tmp[ix] = b[index_c*3+ix];
								c_tmp[ix] = c[index_c*3+ix];
							}
							torque_due_to_force(force_projection2, r_t_tmp, r_tmp, a_tmp, b_tmp, c_tmp, torque2);
							force[index_c*3] = force[index_c*3] + force_projection2[0];
							force[index_c*3+1] = force[index_c*3+1] + force_projection2[1];
							force[index_c*3+2] = force[index_c*3+2] + force_projection2[2];
							torque[index_c*3] = torque[index_c*3] + torque2[0];
							torque[index_c*3+1] = torque[index_c*3+1] + torque2[1];
							torque[index_c*3+2] = torque[index_c*3+2] + torque2[2];
						}

						if (abs(force_projection1[0])>1 || abs(force_projection1[1])>1 || abs(force_projection1[2])>1 || abs(force_projection2[0])>1 || abs(force_projection2[1])>1 || abs(force_projection2[2])>1){
                                			wf << "LHG-tail " << endl;
                        			}

					}
				}
			}
			k2=k2+1;
        	}		
	}

	
	//Electrostatic and Excluded Volum between LH global head and Linker DNA

	k = 0;

	for (i=0;i<n;i++){
		if (type[i]==1){ 
			if (nc_lh_flag[k]==1){
				for (j=0;j<n_lh_g;j++){
					if ((i+j)%comm_size != rank) continue;
					r_lho[0] = r[i*3] + a[i*3]*LH_g_pos[j*3] + b[i*3]*LH_g_pos[j*3+1] + c[i*3]*LH_g_pos[j*3+2];
					r_lho[1] = r[i*3+1] + a[i*3+1]*LH_g_pos[j*3] + b[i*3+1]*LH_g_pos[j*3+1] + c[i*3+1]*LH_g_pos[j*3+2];
					r_lho[2] = r[i*3+2] + a[i*3+2]*LH_g_pos[j*3] + b[i*3+2]*LH_g_pos[j*3+1] + c[i*3+2]*LH_g_pos[j*3+2];
					for (int l=0;l<n;l++){
						if (type[l]==0){
							for (xi = 0; xi<3; xi++){
								r_tmp[xi] = r[l*3+xi];
								LH_g_pos_tmp[xi] = LH_g_pos[j*3+xi];
								a_tmp[xi] = a[i*3+xi];
								b_tmp[xi] = b[i*3+xi];
								c_tmp[xi] = c[i*3+xi];
							}
							Force_Ele_Vdw(LH_q[j], q_l, 1/(4*PI*k_e), debye, k_ex, LH_vdw_hl[j], r_lho, r_tmp,force_projection1, force_projection2, Energy);
							force[i*3] = force[i*3] + force_projection1[0];
							force[i*3+1] = force[i*3+1] + force_projection1[1];
							force[i*3+2] = force[i*3+2] + force_projection1[2];
							force[l*3] = force[l*3] + force_projection2[0];
							force[l*3+1] = force[l*3+1] + force_projection2[1];
							force[l*3+2] = force[l*3+2] + force_projection2[2];

							torque_due_to_force_relative(force_projection1, LH_g_pos_tmp, a_tmp, b_tmp, c_tmp, torque_hgl);
							torque[i*3] = torque[i*3] + torque_hgl[0];
							torque[i*3+1] = torque[i*3+1] + torque_hgl[1];
							torque[i*3+2] = torque[i*3+2] + torque_hgl[2];
						
							if (abs(force_projection1[0])>1 || abs(force_projection1[1])>1 || abs(force_projection1[2])>1 || abs(force_projection2[0])>1 || abs(force_projection2[1])>1 || abs(force_projection2[2])>1){
                                				wf << "LHG-DNA " << endl;
                        				}
						}
					}
				}
			}
			k = k+1;
		}
	}	

	//Electrostatic and Excluded Volum between LH global head and core

	k = 0;

        for (i=0;i<n;i++){
                if (type[i]==1){ 
			if (nc_lh_flag[k]==1){
				for (j=0;j<n_lh_g;j++){
					if ((i+j)%comm_size != rank) continue;
					r_lho[0] = r[i*3] + a[i*3]*LH_g_pos[j*3] + b[i*3]*LH_g_pos[j*3+1] + c[i*3]*LH_g_pos[j*3+2];
					r_lho[1] = r[i*3+1] + a[i*3+1]*LH_g_pos[j*3] + b[i*3+1]*LH_g_pos[j*3+1] + c[i*3+1]*LH_g_pos[j*3+2];
					r_lho[2] = r[i*3+2] + a[i*3+2]*LH_g_pos[j*3] + b[i*3+2]*LH_g_pos[j*3+1] + c[i*3+2]*LH_g_pos[j*3+2];
					for (int l=0;l<n;l++){
						if (i!=l and type[l]==1){
							for (m=0; m<Nq; m++){
								z[0] = r[l*3]+a[l*3]*core_pos[m*3]+b[l*3]*core_pos[m*3+1]+c[l*3]*core_pos[m*3+2];
								z[1] = r[l*3+1]+a[l*3+1]*core_pos[m*3]+b[l*3+1]*core_pos[m*3+1]+c[l*3+1]*core_pos[m*3+2];
								z[2] = r[l*3+2]+a[l*3+2]*core_pos[m*3]+b[l*3+2]*core_pos[m*3+1]+c[l*3+2]*core_pos[m*3+2];
						
								Force_Ele_Vdw(LH_q[j], core_q[m], 1/(4*PI*k_e), debye, k_ex, LH_vdw_hc[j], r_lho, z, force_projection1, force_projection2, Energy);
			
								for (xi = 0; xi<3; xi++){
									LH_g_pos_tmp[xi] = LH_g_pos[j*3+xi];
									a_tmp[xi] = a[i*3+xi];
									b_tmp[xi] = b[i*3+xi];
									c_tmp[xi] = c[i*3+xi];
								}

//								printf("host %d %d %d %d %f %f %f %f %f %f \n", i, j, l, m, force_projection1[0], force_projection1[1], force_projection1[2], r_lho[0], r_lho[1], r_lho[2]);
								force[i*3] = force[i*3] + force_projection1[0];
								force[i*3+1] = force[i*3+1] + force_projection1[1];
								force[i*3+2] = force[i*3+2] + force_projection1[2];
								force[l*3] = force[l*3] + force_projection2[0];
								force[l*3+1] = force[l*3+1] + force_projection2[1];
								force[l*3+2] = force[l*3+2] + force_projection2[2];

								if (abs(force_projection1[0])>1 || abs(force_projection1[1])>1 || abs(force_projection1[2])>1 || abs(force_projection2[0])>1 || abs(force_projection2[1])>1 || abs(force_projection2[2])>1){
                                					wf << "LHG-core " << endl;
                        					}

								torque_due_to_force_relative(force_projection1, LH_g_pos_tmp, a_tmp, b_tmp, c_tmp, torque1);
								torque[i*3] = torque[i*3] + torque1[0];
								torque[i*3+1] = torque[i*3+1] + torque1[1];
								torque[i*3+2] = torque[i*3+2] + torque1[2];
								for (xi = 0; xi<3; xi++){
									core_pos_tmp[xi] = core_pos[m*3+xi];
									a_tmp[xi] = a[l*3+xi];
									b_tmp[xi] = b[l*3+xi];
									c_tmp[xi] = c[l*3+xi];
								}

								torque_due_to_force_relative(force_projection2, core_pos_tmp, a_tmp, b_tmp, c_tmp, torque2);
								torque[l*3] = torque[l*3] + torque2[0];
								torque[l*3+1] = torque[l*3+1] + torque2[1];
								torque[l*3+2] = torque[l*3+2] + torque2[2];
							}
						}
					}
				}
			}
                        k = k+1;
                }
        }

	//Electrostatic and Excluded Volume between LH global head and LH

	k2=0;
	count=0;

        for (i=0;i<n;i++){
                if (type[i]==1){
			 if (nc_lh_flag[k2]==1){
				for (m=0;m<n_lh_g;m++){
					r_lho[0] = r[i*3] + a[i*3]*LH_g_pos[m*3] + b[i*3]*LH_g_pos[m*3+1] + c[i*3]*LH_g_pos[m*3+2];
					r_lho[1] = r[i*3+1] + a[i*3+1]*LH_g_pos[m*3] + b[i*3+1]*LH_g_pos[m*3+1] + c[i*3+1]*LH_g_pos[m*3+2];
					r_lho[2] = r[i*3+2] + a[i*3+2]*LH_g_pos[m*3] + b[i*3+2]*LH_g_pos[m*3+1] + c[i*3+2]*LH_g_pos[m*3+2];

					for (j=0;j<n_LH;j++){
						if ((i+j)%comm_size != rank) continue;
						if (j/(n_lh_n +n_lh_c) != count){
							index2 = j%(n_lh_n+n_lh_c);
							if (index2 < n_lh_n){
								index2 = index2;
							}else{
								index2 = index2 + n_lh_g;
							}

							for (xi = 0; xi<3; xi++){
								r_lh_tmp2[xi] = r_lh[j*3+xi];
								LH_g_pos_tmp[xi] = LH_g_pos[m*3+xi];
								a_tmp[xi] = a[i*3+xi];
								b_tmp[xi] = b[i*3+xi];
								c_tmp[xi] = c[i*3+xi];
							}

							Force_Ele_Vdw(LH_q[m], LH_q[index2], 1/(4*PI*k_e), debye,k_ex, (LH_vdw_hh[m]+LH_vdw_hh[index2])/2, r_lho, r_lh_tmp2,force_projection1, force_projection2, Energy);

							force[i*3] = force[i*3] + force_projection1[0];
							force[i*3+1] = force[i*3+1] + force_projection1[1];
							force[i*3+2] = force[i*3+2] + force_projection1[2];

							LH_force[j*3] = LH_force[j*3] + force_projection2[0];
                                			LH_force[j*3+1] = LH_force[j*3+1] + force_projection2[1];
                                			LH_force[j*3+2] = LH_force[j*3+2] + force_projection2[2];

							torque_due_to_force_relative(force_projection1, LH_g_pos_tmp, a_tmp, b_tmp, c_tmp, torque_hgl);
							torque[i*3] = torque[i*3] + torque_hgl[0];
							torque[i*3+1] = torque[i*3+1] + torque_hgl[1];
							torque[i*3+2] = torque[i*3+2] + torque_hgl[2];

							if (abs(force_projection1[0])>1 || abs(force_projection1[1])>1 || abs(force_projection1[2])>1 || abs(force_projection2[0])>1 || abs(force_projection2[1])>1 || abs(force_projection2[2])>1){
                                				wf << "LHG-LH " << endl;
                        				}
						}
					}
				}
				count = count+1;
			}
                	k2=k2+1;
                }
        }


	//Electrostatic and Excluded Volume between LH global head and LH global head

	k2=0;

        for (i=0;i<n;i++){
                if (type[i]==1){
                         if (nc_lh_flag[k2]==1){
                                for (m=0;m<n_lh_g;m++){
                                        r_lho[0] = r[i*3] + a[i*3]*LH_g_pos[m*3] + b[i*3]*LH_g_pos[m*3+1] + c[i*3]*LH_g_pos[m*3+2];
                                        r_lho[1] = r[i*3+1] + a[i*3+1]*LH_g_pos[m*3] + b[i*3+1]*LH_g_pos[m*3+1] + c[i*3+1]*LH_g_pos[m*3+2];
                                        r_lho[2] = r[i*3+2] + a[i*3+2]*LH_g_pos[m*3] + b[i*3+2]*LH_g_pos[m*3+1] + c[i*3+2]*LH_g_pos[m*3+2];

					k = k2+1;
                                        for (j=i+1;j<n;j++){
						if (type[j]==1){
							if ((i+j)%comm_size != rank) continue;
							if (nc_lh_flag[k]==1){
								for (m1=0;m1<n_lh_g;m1++){
									r_lho2[0] = r[j*3] + a[j*3]*LH_g_pos[m1*3] + b[j*3]*LH_g_pos[m1*3+1] + c[j*3]*LH_g_pos[m1*3+2];
									r_lho2[1] = r[j*3+1] + a[j*3+1]*LH_g_pos[m1*3] + b[j*3+1]*LH_g_pos[m1*3+1] + c[j*3+1]*LH_g_pos[m1*3+2];
									r_lho2[2] = r[j*3+2] + a[j*3+2]*LH_g_pos[m1*3] + b[j*3+2]*LH_g_pos[m1*3+1] + c[j*3+2]*LH_g_pos[m1*3+2];

									Force_Ele_Vdw(LH_q[m], LH_q[m1], 1/(4*PI*k_e), debye,k_ex, (LH_vdw_hh[m]+LH_vdw_hh[m1])/2, r_lho, r_lho2,force_projection1, force_projection2, Energy);

									force[i*3] = force[i*3] + force_projection1[0];
									force[i*3+1] = force[i*3+1] + force_projection1[1];
									force[i*3+2] = force[i*3+2] + force_projection1[2];

									force[j*3] = force[j*3] + force_projection2[0];
									force[j*3+1] = force[j*3+1] + force_projection2[1];
									force[j*3+2] = force[j*3+2] + force_projection2[2];



									if (abs(force_projection1[0])>1 || abs(force_projection1[1])>1 || abs(force_projection1[2])>1 || abs(force_projection2[0])>1 || abs(force_projection2[1])>1 || abs(force_projection2[2])>1){
                                						wf << "LHG-LHG " << endl;
                        						}

									for (xi = 0; xi<3; xi++){
										LH_g_pos_tmp[xi] = LH_g_pos[m*3+xi];
										a_tmp[xi] = a[i*3+xi];
										b_tmp[xi] = b[i*3+xi];
										c_tmp[xi] = c[i*3+xi];
									}

									torque_due_to_force_relative(force_projection1, LH_g_pos_tmp, a_tmp, b_tmp, c_tmp, torque_hgl);
									torque[i*3] = torque[i*3] + torque_hgl[0];
									torque[i*3+1] = torque[i*3+1] + torque_hgl[1];
									torque[i*3+2] = torque[i*3+2] + torque_hgl[2];

									for (xi = 0; xi<3; xi++){
										LH_g_pos_tmp[xi] = LH_g_pos[m1*3+xi];
										a_tmp[xi] = a[j*3+xi];
										b_tmp[xi] = b[j*3+xi];
										c_tmp[xi] = c[j*3+xi];
									}

									torque_due_to_force_relative(force_projection2, LH_g_pos_tmp, a_tmp, b_tmp, c_tmp, torque_hgl);
									torque[j*3] = torque[j*3] + torque_hgl[0];
									torque[j*3+1] = torque[j*3+1] + torque_hgl[1];
									torque[j*3+2] = torque[j*3+2] + torque_hgl[2];

								}
							}
							k=k+1;
						}
                                        }
                                }
                       
                        }
                        k2=k2+1;
                }
        }





	//Stretching between LH global head with C-term

        k = 0;
	k2 = 0;

        for (i = 0; i < n; i++){
                if (type[i]==1){ 
			if (nc_lh_flag[k2]==1){
				for (j = 0; j < n_lh_g; j++){
					if ((i+j)%comm_size != rank) continue;
					if (LH_conn[j+n_lh_n]==1){
						r_lho[0] = r[i*3] + a[i*3]*LH_g_pos[j*3] + b[i*3]*LH_g_pos[j*3+1] + c[i*3]*LH_g_pos[j*3+2];
						r_lho[1] = r[i*3+1] + a[i*3+1]*LH_g_pos[j*3] + b[i*3+1]*LH_g_pos[j*3+1] + c[i*3+1]*LH_g_pos[j*3+2];
						r_lho[2] = r[i*3+2] + a[i*3+2]*LH_g_pos[j*3] + b[i*3+2]*LH_g_pos[j*3+1] + c[i*3+2]*LH_g_pos[j*3+2];

						distance = sqrt((r_lh[k*nlh*3]-r_lho[0])*(r_lh[k*nlh*3]-r_lho[0])+(r_lh[k*nlh*3+1]-r_lho[1])*(r_lh[k*nlh*3+1]-r_lho[1])+(r_lh[k*nlh*3+2]-r_lho[2])*(r_lh[k*nlh*3+2]-r_lho[2]));
						
						Energy = Energy + LH_kstr[j+n_lh_n]*(distance-LH_streq[j+n_lh_n])*(distance-LH_streq[j+n_lh_n])/2;
			

						Stri[0] = (distance - LH_streq[j+n_lh_n])*(r_lh[k*nlh*3]-r_lho[0])/distance;
						Stri[1] = (distance - LH_streq[j+n_lh_n])*(r_lh[k*nlh*3+1]-r_lho[1])/distance;
						Stri[2] = (distance - LH_streq[j+n_lh_n])*(r_lh[k*nlh*3+2]-r_lho[2])/distance;


						LH_force[k*nlh*3] = LH_force[k*nlh*3] - LH_kstr[j+n_lh_n]*Stri[0]*10;
						LH_force[k*nlh*3+1] = LH_force[k*nlh*3+1] - LH_kstr[j+n_lh_n]*Stri[1]*10;
						LH_force[k*nlh*3+2] = LH_force[k*nlh*3+2] - LH_kstr[j+n_lh_n]*Stri[2]*10;

						force[i*3] = force[i*3] + LH_kstr[j+n_lh_n]*Stri[0]*10;
						force[i*3+1] = force[i*3+1] + LH_kstr[j+n_lh_n]*Stri[1]*10;
						force[i*3+2] = force[i*3+2] + LH_kstr[j+n_lh_n]*Stri[2]*10;


						if (abs(LH_kstr[j+n_lh_n]*Stri[0]*10)>1 || abs(LH_kstr[j+n_lh_n]*Stri[1]*10)>1 || abs(LH_kstr[j+n_lh_n]*Stri[2]*10)>1){
                                			wf << "LH-LHG " << endl;
                        			}


						fa = LH_kstr[j+n_lh_n]*(a[i*3]*Stri[0]+a[i*3+1]*Stri[1]+a[i*3+2]*Stri[2]);
						fb = LH_kstr[j+n_lh_n]*(b[i*3]*Stri[0]+b[i*3+1]*Stri[1]+b[i*3+2]*Stri[2]);
						fc = LH_kstr[j+n_lh_n]*(c[i*3]*Stri[0]+c[i*3+1]*Stri[1]+c[i*3+2]*Stri[2]);
						torque[i*3] = torque[i*3] + fc*LH_g_pos[j*3+1] - fb*LH_g_pos[j*3+2];
						torque[i*3+1] = torque[i*3+1] + fa*LH_g_pos[j*3+2] - fc*LH_g_pos[j*3];
						torque[i*3+2] = torque[i*3+2] + fb*LH_g_pos[j*3] - fa*LH_g_pos[j*3+1];

					}
				}
				k=k+1;
			}
                        k2 = k2+1;
                }
        }

	wf.close();

}
