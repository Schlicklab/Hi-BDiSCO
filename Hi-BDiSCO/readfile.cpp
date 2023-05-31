// **************************************************************************************//
//											 //
//		 		Brownian Dynamics Simulation Algorithm			 //
// 		     Copyright Zilong Li, Tamar Schlick and New York University 	 //
//					    April 2020					 //
//                                                                                       //
// **************************************************************************************//


#include "readfile.h"

using namespace std;

std::vector<double> Read_const(string filename){

	std::vector<double> input_const;


	ifstream ReadFile(filename);

	int i=0;
	string line;
	while (getline (ReadFile, line)){
		if (i!=0){
			std::istringstream iss(line);
                        std::vector<std::string> result((std::istream_iterator<std::string>(iss)),std::istream_iterator<std::string>());

			stringstream  consts(result[1]);
			double tmp_const;
			consts >> tmp_const;
			input_const.push_back(tmp_const);
		}
		i=i+1;
	}

	return input_const;

}

void Read_initial_struc(int &n, int &n_c, int &n3, int &nc3, std::vector<double> &tmp_r, std::vector<double> &tmp_a, std::vector<double> &tmp_b, std::vector<double> &tmp_c, std::vector<double> &tmp_rad, std::vector<double> &tmp_Er, std::vector<int> &tmp_type, string filename){

	ifstream ReadFile(filename);


	int i = 0;
	n = 0;
	n_c = 0;

	string line;
	while (getline (ReadFile, line)){
		std::istringstream iss(line);
                std::vector<std::string> result((std::istream_iterator<std::string>(iss)),std::istream_iterator<std::string>());

		if (i%4 == 0){
			stringstream xs(result[0]);
			stringstream ys(result[1]);
			stringstream zs(result[2]);
			stringstream ts(result[3]);

			double x,y,z;
			int t;

			xs >> x;
			ys >> y;
			zs >> z;
			ts >> t;
			
			tmp_r.push_back(x);
			tmp_r.push_back(y);
			tmp_r.push_back(z);
			tmp_type.push_back(t);
			if (t==1){
				n_c = n_c+1;
				tmp_rad.push_back(5.0);
				tmp_Er.push_back(8*PI*eta*125.0);
				tmp_Er.push_back(8*PI*eta*125.0);
				tmp_Er.push_back(8*PI*eta*125.0);
			}else{
				tmp_rad.push_back(1.5);
                                tmp_Er.push_back(4*PI*eta*r_h*r_h*lo);
                                tmp_Er.push_back(0.0);
                                tmp_Er.push_back(0.0);
			}
			n = n+1;		
		}else if (i%4==1){
                        stringstream xs(result[0]);
                        stringstream ys(result[1]);
                        stringstream zs(result[2]);
                
                        double x,y,z;

                        xs >> x;
                        ys >> y;
                        zs >> z;

			tmp_a.push_back(x);
			tmp_a.push_back(y);
			tmp_a.push_back(z);
		}else if (i%4==2){
                        stringstream xs(result[0]);
                        stringstream ys(result[1]);
                        stringstream zs(result[2]);

                        double x,y,z;

                        xs >> x;
                        ys >> y;
                        zs >> z;

                        tmp_b.push_back(x);
                        tmp_b.push_back(y);
                        tmp_b.push_back(z);
                }else if (i%4==3){
                        stringstream xs(result[0]);
                        stringstream ys(result[1]);
                        stringstream zs(result[2]);

                        double x,y,z;

                        xs >> x;
                        ys >> y;
                        zs >> z;

                        tmp_c.push_back(x);
                        tmp_c.push_back(y);
                        tmp_c.push_back(z);
                }
		i = i+1;
	}

	n3 = n*3;
	nc3 = n_c*3;
	
}


void Read_core(int &Nq, int &Nq3, std::vector<double> &core_pos, std::vector<double> &core_q, string filename){

        ifstream ReadFile(filename);

        string line;

        Nq=0;

        while (getline (ReadFile, line)) {
                std::istringstream iss(line);
                std::vector<std::string> result((std::istream_iterator<std::string>(iss)),std::istream_iterator<std::string>());
                
                if (result.size() == 3){
                        stringstream xs(result[0]);
                        stringstream ys(result[1]);
                        stringstream zs(result[2]);

                        double x, y, z;

                        xs >> x;
                        ys >> y;
                        zs >> z;
                       
			core_pos.push_back(x/10);
			core_pos.push_back(y/10);
			core_pos.push_back(z/10);
                        Nq = Nq+1;
                }else{
                        stringstream charges(result[0]);
                        double charge;
                        charges >> charge;
                        core_q.push_back(charge);
                }
        }
	Nq3 = Nq*3;


}

void Read_extra(int &ex_n, std::vector<int> &ex_m, string filename){

        ifstream ReadFile(filename);

        string line;

        ex_n=0;

        while (getline (ReadFile, line)) {
                std::istringstream iss(line);
                std::vector<std::string> result((std::istream_iterator<std::string>(iss)),std::istream_iterator<std::string>());

		stringstream xs(result[0]);
		stringstream ys(result[1]);

		int x,y;

		xs >> x;
		ys >> y;
	
		ex_m.push_back(x);
		ex_m.push_back(y);

		ex_n+=1;

        }


}


void Read_tail(int &n_t, std::vector<double> &tail_pos, std::vector<double> &tail_q, std::vector<int> &tail_grp, std::vector<int> &tail_fix, std::vector<double> &tail_rad, std::vector<double> &tail_bond_c, std::vector<double> &tail_bond_v, std::vector<double> &tail_angle_c, std::vector<double> &tail_angle_v, string filename){

        ifstream ReadFile(filename);

        string line;

        n_t=0;
        int t_grp, t_fix;
        double t_chg, t_rad, t_mass, x, y, z, t_bond_v, t_bond_c, t_angle_v, t_angle_c; 


        while (getline (ReadFile, line)) {
                std::istringstream iss(line);
                std::vector<std::string> result((std::istream_iterator<std::string>(iss)),std::istream_iterator<std::string>());

                stringstream t_grps(result[1]);
                stringstream t_fixs(result[2]);
                stringstream t_chgs(result[3]);
                stringstream t_rads(result[4]);
                stringstream t_masss(result[5]);
                stringstream xs(result[6]);
                stringstream ys(result[7]);
                stringstream zs(result[8]);
                stringstream t_bond_vs(result[9]);
                stringstream t_bond_cs(result[10]);
                stringstream t_angle_vs(result[11]);
                stringstream t_angle_cs(result[12]);

                t_grps >> t_grp;
                t_fixs >> t_fix;
                t_chgs >> t_chg;
                t_rads >> t_rad;
                t_masss >> t_mass;
                xs >> x;
                ys >> y;
                zs >> z;
                t_bond_vs >> t_bond_v;
                t_bond_cs >> t_bond_c;
                t_angle_vs >> t_angle_v;
                t_angle_cs >> t_angle_c;

		tail_pos.push_back(x/10);
		tail_pos.push_back(y/10);
		tail_pos.push_back(z/10);

		tail_q.push_back(t_chg*1.12);
		tail_grp.push_back(t_grp);
		tail_fix.push_back(t_fix);
		tail_rad.push_back(t_rad/10);
		tail_bond_c.push_back(t_bond_c*100*4.142e-3/0.5962);
		tail_bond_v.push_back(t_bond_v/10);
		tail_angle_c.push_back(t_angle_c*4.142e-3/0.5962);
		tail_angle_v.push_back(t_angle_v*PI/180);

//		cout << n_t << endl;
                n_t = n_t+1;
        }
 

}


void Read_LH(int &n_lh_g, int &n_lh_n, int &n_lh_c, std::vector<double> &LH_g_pos, std::vector<double> &LH_n_pos, std::vector<double> &LH_c_pos, std::vector<int> &LH_conn, std::vector<double> &LH_q, std::vector<double> &LH_vdw_hh, std::vector<double> &LH_vdw_hc, std::vector<double> &LH_vdw_hl, std::vector<double> &LH_vdw_ht, std::vector<double> &LH_kstr, std::vector<double> &LH_kben, std::vector<double> &LH_streq, std::vector<double> &LH_betaeq, std::vector<double> &LH_radius, string filename){

        ifstream ReadFile(filename);

        string line;
        string LH_type;

	n_lh_g = 0;
	n_lh_n = 0;
	n_lh_c = 0;

        while (getline (ReadFile, line)) {
        
                std::istringstream iss(line);
                std::vector<std::string> result((std::istream_iterator<std::string>(iss)),std::istream_iterator<std::string>());
                
                double x, y, z, charge, Vdw_hh, Vdw_hc, Vdw_hl, Vdw_ht, kstr, kben, ktor, rstreq, betaeq;
                int seq, conn;


                if (result[1] == "N-term"){
                        LH_type = "LH_N";
                }else if (result[1] == "G-head"){
                        LH_type = "LH_G";
                }else if (result[1] == "C-term"){
                        LH_type = "LH_C";
                }else if (result[0] != "#"){
                        stringstream xs(result[0]);
                        stringstream ys(result[1]);
                        stringstream zs(result[2]);
        
                        stringstream charges(result[3]);

                        stringstream Vdw_hhs(result[4]);
                        stringstream Vdw_hcs(result[5]);
                        stringstream Vdw_hls(result[6]);
                        stringstream Vdw_hts(result[7]);

                        stringstream kstrs(result[8]);
                        stringstream kbens(result[9]);
                        stringstream ktors(result[10]);
        
                        stringstream seqs(result[11]);
                        stringstream conns(result[12]);

                        stringstream rstreqs(result[13]);
                        stringstream betaeqs(result[14]);

                        xs >> x;
                        ys >> y;
                        zs >> z;
                        charges >> charge;
                        Vdw_hhs >> Vdw_hh;
                        Vdw_hcs >> Vdw_hc;
                        Vdw_hls >> Vdw_hl;
                        Vdw_hts >> Vdw_ht;
                        kstrs >> kstr;
                        kbens >> kben;
                        ktors >> ktor;
                        seqs >> seq;
                        conns >> conn;
                        rstreqs >> rstreq;
                        betaeqs >> betaeq;

			if (LH_type == "LH_N"){
				n_lh_n = n_lh_n+1;
				LH_n_pos.push_back(x);
				LH_n_pos.push_back(y);
				LH_n_pos.push_back(z);
			}else if (LH_type == "LH_G"){
				n_lh_g = n_lh_g+1;
				LH_g_pos.push_back(x);
                                LH_g_pos.push_back(y);
                                LH_g_pos.push_back(z);
			}else if (LH_type == "LH_C"){
				n_lh_c = n_lh_c+1;
                                LH_c_pos.push_back(x);
                                LH_c_pos.push_back(y);
                                LH_c_pos.push_back(z);
			}
	

			LH_q.push_back(charge*0.66);
			LH_vdw_hh.push_back(Vdw_hh);
			LH_vdw_hc.push_back(Vdw_hc);
			LH_vdw_hl.push_back(Vdw_hl);
			LH_vdw_ht.push_back(Vdw_ht);
			LH_kstr.push_back(kstr*4.142e-3/0.5962);
			LH_kben.push_back(kben*4.142e-3/0.5962);
			LH_streq.push_back(rstreq);
			LH_betaeq.push_back(betaeq*PI/180);
			LH_conn.push_back(conn);
			LH_radius.push_back(0.5);


                }
                
        }

}








void write_xyz_append(int n, int n3, int n_c, int Nq, int Nq3, int* type, double* r, double* a, double* b, double* c, std::vector<double> core_pos, string filename){

        std::ofstream outfile;
        outfile.open(filename, std::ios_base::app);

	double z[3];
	int k = 0;

	outfile << n*4+n_c*Nq << endl;
        outfile << "comments zl3765" << endl;

	for (int i = 0; i < n; i++){
		if (type[i]==0){
			outfile << "CA" << " " << r[i*3] << " " << r[i*3+1] << " " << r[i*3+2] << endl;
			outfile << "H1" << " " << a[i*3] << " " << a[i*3+1] << " " << a[i*3+2] << endl;
			outfile << "H2" << " " << b[i*3] << " " << b[i*3+1] << " " << b[i*3+2] << endl;
			outfile << "H3" << " " << c[i*3] << " " << c[i*3+1] << " " << c[i*3+2] << endl;
		}else{
			outfile << "OC" << " " << r[i*3] << " " << r[i*3+1] << " " << r[i*3+2] << endl;
			outfile << "H1" << " " << a[i*3] << " " << a[i*3+1] << " " << a[i*3+2] << endl;
                        outfile << "H2" << " " << b[i*3] << " " << b[i*3+1] << " " << b[i*3+2] << endl;
                        outfile << "H3" << " " << c[i*3] << " " << c[i*3+1] << " " << c[i*3+2] << endl;
			for (int j = 0; j < Nq; j++){
				z[0] = r[i*3] + a[i*3]*core_pos[j*3] + b[i*3]*core_pos[j*3+1] + c[i*3]*core_pos[j*3+2];
				z[1] = r[i*3+1] + a[i*3+1]*core_pos[j*3] + b[i*3+1]*core_pos[j*3+1] + c[i*3+1]*core_pos[j*3+2];
				z[2] = r[i*3+2] + a[i*3+2]*core_pos[j*3] + b[i*3+2]*core_pos[j*3+1] + c[i*3+2]*core_pos[j*3+2];
				outfile << "NC" << " " << z[0] << " " << z[1] << " " << z[2] << endl;
			}
		}

	}

	outfile.close();


}

void write_restart(int n3, double* r, double* a, double* b, double* c){

	ofstream restart_r ("restart_r.txt", std::ofstream::trunc);
	ofstream restart_a ("restart_a.txt", std::ofstream::trunc);
	ofstream restart_b ("restart_b.txt", std::ofstream::trunc);
	ofstream restart_c ("restart_c.txt", std::ofstream::trunc);

	for (int i=0; i<n3; i++){
		restart_r << std::fixed << std::setprecision(15) << r[i] << endl;
		restart_a << std::fixed << std::setprecision(15) << a[i] << endl;
		restart_b << std::fixed << std::setprecision(15) << b[i] << endl;
		restart_c << std::fixed << std::setprecision(15) << c[i] << endl;
	}



	restart_r.close();
	restart_a.close();
	restart_b.close();
	restart_c.close();

}

void Read_restart_ini(int n3, double* r, double* a, double* b, double* c){

        ifstream ReadFile("restart_r.txt");

        int i=0;
        string line;
        while (getline (ReadFile, line)){
		std::istringstream iss(line);
                std::vector<std::string> result((std::istream_iterator<std::string>(iss)),std::istream_iterator<std::string>());
		if (i<n3){
			stringstream  consts(result[0]);
                        double tmp_const;
                        consts >> tmp_const;
			r[i] = tmp_const;
		}
                i=i+1;
        }

	ifstream ReadFile_a("restart_a.txt");

        i=0;
        while (getline (ReadFile_a, line)){
                std::istringstream iss(line);
                std::vector<std::string> result((std::istream_iterator<std::string>(iss)),std::istream_iterator<std::string>());
                if (i<n3){
                        stringstream  consts(result[0]);
                        double tmp_const;
                        consts >> tmp_const;
                        a[i] = tmp_const;
                }
                i=i+1;
        }

	ifstream ReadFile_b("restart_b.txt");

        i=0;
        while (getline (ReadFile_b, line)){
                std::istringstream iss(line);
                std::vector<std::string> result((std::istream_iterator<std::string>(iss)),std::istream_iterator<std::string>());
                if (i<n3){
                        stringstream  consts(result[0]);
                        double tmp_const;
                        consts >> tmp_const;
                        b[i] = tmp_const;
                }
                i=i+1;
        }

	ifstream ReadFile_c("restart_c.txt");

        i=0;
        while (getline (ReadFile_c, line)){
                std::istringstream iss(line);
                std::vector<std::string> result((std::istream_iterator<std::string>(iss)),std::istream_iterator<std::string>());
                if (i<n3){
                        stringstream  consts(result[0]);
                        double tmp_const;
                        consts >> tmp_const;
                        c[i] = tmp_const;
                }
                i=i+1;
        }

}

void Read_restart_tail(int n_tail3, double* r_t){

	ifstream ReadFile("restart_tail.txt");

        int i=0;
        string line;
        while (getline (ReadFile, line)){
                std::istringstream iss(line);
                std::vector<std::string> result((std::istream_iterator<std::string>(iss)),std::istream_iterator<std::string>());
                if (i<n_tail3){
                        stringstream  consts(result[0]);
                        double tmp_const;
                        consts >> tmp_const;
                        r_t[i] = tmp_const;
                }
                i=i+1;
        }

}

void Read_restart_LH(int n_LH3, double* r_lh){

        ifstream ReadFile("restart_lh.txt");

        int i=0;
        string line;
        while (getline (ReadFile, line)){
                std::istringstream iss(line);
                std::vector<std::string> result((std::istream_iterator<std::string>(iss)),std::istream_iterator<std::string>());
                if (i<n_LH3){
                        stringstream  consts(result[0]);
                        double tmp_const;
                        consts >> tmp_const;
                        r_lh[i] = tmp_const;
                }
                i=i+1;
        }

}

void Read_LH_core(int n_c, int &nc_lh, int* nc_lh_flag){

	ifstream ReadFile("LH.in");

	int i=0;
	string line;
	while (getline (ReadFile, line)){
                std::istringstream iss(line);
                std::vector<std::string> result((std::istream_iterator<std::string>(iss)),std::istream_iterator<std::string>());
                if (i<n_c){
                        stringstream  consts(result[0]);
                        int tmp_const;
                        consts >> tmp_const;
                        nc_lh_flag[i] = tmp_const;
			nc_lh = nc_lh + nc_lh_flag[i];
                }
                i=i+1;
        }



}






