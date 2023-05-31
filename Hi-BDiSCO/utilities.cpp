
#include "utilities.h"

using namespace std;


void Message(){
    
    cout << endl << endl;
    cout << "***********************************************************************************" << endl;
    cout << "*                                                                                 *" << endl;
    cout << "*                                 FOLDING ENGINE                                  *" << endl;
    cout << "*                                                                                 *" << endl;
    cout << "*                      (a chromatin folding simulation program)                   *" << endl;
    cout << "*                                                                                 *" << endl;
    cout << "*                                  Version 2.05b                                  *" << endl;
    cout << "*                                                                                 *" << endl;
    cout << "*     (c) Ognjen Perisic, Tamar Schlick and New York University, 2017 - 2019.     *" << endl;
    cout << "*                                                                                 *" << endl;
    cout << "*               email: ognjen.perisic@gmail.com, schlick@nyu.edu                  *" << endl;
    cout << "*                                                                                 *" << endl;
    cout << "*                                 January 29 2019.                               *" << endl;
    cout << "*                                                                                 *" << endl;
    cout << "***********************************************************************************" << endl;
    cout << endl << endl << endl;

    cout << "    References:" << endl << endl;

    cout << "   [1] Schlick, T.; Li, B.; Olson, W.K. \"The Influence of Salt on the" << endl;
    cout << "       Structure and Energetics of Supercoiled DNA\". " << endl;
    cout << "       Biophys. J. (1994), 67, 2146 - 2166." << endl << endl;

    cout << "   [2] Beard, D.; Schlick, T. \"Modeling Salt - Mediated Electrostatics of" << endl;
    cout << "       Macromolecules : The Discrete Surface Charge Optimization Algorithm" << endl;
    cout << "       and Its Application to the Nucleosome\"." << endl;
    cout << "       Biopolymers(2001), 58, 106 - 115." << endl << endl;

    cout << "   [3] Zhang, Q.; Beard, D.; Schlick, T. \"Constructing Irregular Surfaces" << endl;
    cout << "       to Enclose Macromolecular Complexes for Mesoscale Modeling Using" << endl;
    cout << "       the Discrete Surface Charge Optimization(DISCO) Algorithm\"." << endl;
    cout << "       J. Comput. Chem. (2003), 24, 2063 - 2074." << endl << endl;

    cout << "   [4] Arya, G.; Schlick, T. \"Role of Histone Tails in Chromatin Folding" << endl;
    cout << "       Revealed by a New Mesoscopic Oligonucleosome Model\"." << endl;
    cout << "       Proc. Natl. Acad.Sci. USA (2006), 103, 16236 - 16241." << endl << endl;

    cout << "   [5] Arya, G.; Schlick, T. \"A Tale of Tails : How Histone Tails" << endl;
    cout << "       Mediate Chromatin Compaction in Different Salt and Linker Histone" << endl;
    cout << "       Environments\". J. Phys. Chem. A (2009), 113, 4045 - 4059." << endl << endl;

    cout << "   [6] Schlick, T.; Perisic, O. \"Mesoscale Simulations of Two" << endl;
    cout << "       Nucleosome - Repeat Length Oligonucleosomes\"." << endl;
    cout << "       Phys. Chem. Chem. Phys. (2009), 11, 10729 - 10737." << endl << endl;

    cout << "   [7] Perisic, O.; Collepardo - Guevara, R.; Schlick, T. \"Modeling Studies" << endl;
    cout << "       of Chromatin Fiber Structure as a Function of DNA Linker Length\"." << endl;
    cout << "       J. Mol. Biol. (2010), 403, 777 - 802." << endl << endl;

    cout << "   [8] Collepardo - Guevara, R.; Schlick, T. \"Crucial Role of Dynamic" << endl;
    cout << "       Linker Histone Binding and Divalent Ions for DNA Accessibility and" << endl;
    cout << "       Gene Regulation Revealed by Mesoscale Modeling of Oligonucleosomes\"." << endl;
    cout << "       Nucleic Acids Res. (2012), 40, 8803 - 8817." << endl << endl;

    cout << "   [9] Collepardo - Guevara, R.; Schlick, T. \"Chromatin fiber polymorphism" << endl;
    cout << "       triggered by variations of DNA linker lengths\"." << endl;
    cout << "       Proc. Natl. Acad. Sci. USA (2014), 111, 8061 - 8066." << endl << endl;

    cout << "  [10] Luque, A.; Collepardo - Guevara, R.; Grigoryev, S.; Schlick, T." << endl;
    cout << "       \"Dynamic Condensation of Linker Histone C - terminal Domain" << endl;
    cout << "       Regulates Chromatin Structure\". Nucleic Acids Res." << endl;
    cout << "       (2014), 42, 7553 - 7560." << endl << endl;

    cout << "  [11] Luque, A.; Ozer, G.; Schlick, T. \"Correlation among DNA Linker" << endl;
    cout << "       Length, Linker Histone Concentration, and Histone Tails in" << endl;
    cout << "       Chromatin\". Biophys.J. (2016), 110, 2309 - 2319." << endl << endl;

    cout << "  [12] Grigoryev, S.A.; Bascom, G.; Buckwalter, J.M.; Schubert, M.B.;" << endl;
    cout << "       Woodcock, C.L.; Schlick, T. \"Hierarchical Looping of Zigzag" << endl;
    cout << "       Nucleosome Chains in Metaphase Chromosomes\"." << endl;
    cout << "       Proc. Natl. Acad. Sci. USA (2016), 113, 1238 - 1243." << endl << endl;

    cout << "  [13] Bascom, G.D.; Sanbonmatsu, K.Y.; Schlick, T. \"Mesoscale Modeling" << endl;
    cout << "       Reveals Hierarchical Looping of Chromatin Fibers near Gene" << endl;
    cout << "       Regulatory Elements\". J. Phys. Chem. B (2016), 120, 8642 - 8653." << endl << endl;

    cout << "  [14] Bascom, G.; Schlick, T. \"Linking Chromatin Fibers to Gene" << endl;
    cout << "       Folding by Hierarchical Looping. Biophys. J. (2017), 112, 434 - 445." << endl << endl;

    cout << "  [15] Perisic O.; Schlick T. Dependence of the Linker Histone" << endl;
    cout << "       and Chromatin Condensation on the Nucleosome Environment\"." << endl;
    cout << "       J. Phys. Chem. B (2017), 121:7823 - 7832" << endl << endl;

    cout << "  [16] Bascom G.; Schlick T. \"Chromatin Fiber Folding Directed by" << endl;
    cout << "       Cooperative Histone Tail Acetylation and Linker Histone Binding\"." << endl;
    cout << "       Biophys J. (2018), 114(10): 2376-2385" << endl << endl;


    cout << endl;
    cout << endl;
}

double norm(vector<double> &r1, vector<double> &r2){
    int i;
    double norma = 0;
    for (i = 0; i<3; i++){
        norma = norma + (r1[i] - r2[i])*(r1[i] - r2[i]);
    }
    return sqrt(norma);
}



double norm(double r1[3], vector<double> &r2){
    int i;
    double norma = 0;
    for (i = 0; i<3; i++){
        norma = norma + (r1[i] - r2[i])*(r1[i] - r2[i]);
    }
    return sqrt(norma);
}

double norm(double r1[3], double r2[3]){
    int i;
    double norma = 0;
    for (i = 0; i<3; i++){
        norma = norma + (r1[i] - r2[i])*(r1[i] - r2[i]);
    }
    return sqrt(norma);
}

double norm(vector<double> r){
    double norma = 0;
    norma = norma + r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
    
    return sqrt(norma);
}


vector<double> cross_product(vector<double> &r1, vector<double> &r2){
    vector<double> product;
    product.resize(3, 0.0);

    product[0] = r1[1] * r2[2] - r1[2] * r2[1];
    product[1] = -1 * (r1[0] * r2[2] - r1[2] * r2[0]);
    product[2] = r1[0] * r2[1] - r1[1] * r2[0];

    return product;
}


vector <double> cross_product(double r1[3], vector < double > &r2){
    vector<double> product;
    product.resize(3, 0.0);

    product[0] = r1[1] * r2[2] - r1[2] * r2[1];
    product[1] = -1 * (r1[0] * r2[2] - r1[2] * r2[0]);
    product[2] = r1[0] * r2[1] - r1[1] * r2[0];

    return product;
}

double dot_product(vector<double> &r1, vector<double> &r2){
    return r1[0] * r2[0] + r1[1] * r2[1] + r1[2] * r2[2];
}

double dot_product(double r1[3], double r2[3]){
    return r1[0] * r2[0] + r1[1] * r2[1] + r1[2] * r2[2];
}

double dot_product(vector<double> &r1, double r2[3]){
    return r1[0] * r2[0] + r1[1] * r2[1] + r1[2] * r2[2];
}


double dot_product(double r1[3], vector < double > &r2){
    return r1[0] * r2[0] + r1[1] * r2[1] + r1[2] * r2[2];
}

vector <double> multiply_vector_by_matrix(vector<double> v, vector< vector <double > > m){
    vector <double> product;
    int i, j;
    product.resize(3, 0.0);

    for (i = 0; i < 3; i++){
        for (j = 0; j < 3; j++){
            product[i] = product[i] + m[i][j] * v[j];
        }
    }
    return product;
}




void multiply_array_by_matrix(double v[3], vector< vector <double > > m){
    double product[3];
    int i, j;
    product[0] = 0.0; product[1] = 0.0; product[2] = 0.0;

    for (i = 0; i < 3; i++){
        for (j = 0; j < 3; j++){
            product[i] = product[i] + m[i][j] * v[j];
        }
    }
    for (i = 0; i < 3; i++)
        v[i] = product[i];    
}



int scopyn(char *t, char *s, int n){
    int i;
    if (n <= strlen(s)){
        for (i = 0; i<n; i++){
            *(t + i) = *(s + i);
        }
        t[i++] = '\0';
        return 0;
    }
    return -1;
}


int scopyn2n(char *t, char *s, int n1, int n2){
    int i, j;
    if ((n2 <= strlen(s))&(n1 >= 0)&(n1 <= n2)) {
        j = 0;
        for (i = n1 - 1; i<n2; i++){
            *(t + j++) = *(s + i);
        }
        t[++j] = '\0';
        return 0;
    }
    return -1;
}


double normalize_angle(double x){
    x = x - 2 * PI*floor((x + PI) / (2 * PI));
    
    return x;
}



void multiply_array_by_inverse_matrix(double v[3], double a[3], double b[3], double c[3]){
    double n_v[3];
    double inv_mat[3][3];
    double det, invdet;

    det = a[0] * (b[1] * c[2] - c[1] * b[2]) -
        a[1] * (b[0] * c[2] - b[2] * c[0]) +
        a[2] * (b[0] * c[1] - b[1] * c[0]);

    invdet = 1 / det;

    inv_mat[0][0] = (b[1]*c[2]-c[1]*b[2])*invdet;
    inv_mat[0][1] = (a[2]*c[1]-a[1]*c[2])*invdet;
    inv_mat[0][2] = (a[1]*b[2]-a[2]*b[1])*invdet;

    inv_mat[1][0] = (b[2]*c[0]-b[0]*c[2])*invdet;
    inv_mat[1][1] = (a[0]*c[2]-a[2]*c[0])*invdet;
    inv_mat[1][2] = (b[0]*a[2]-a[0]*b[2])*invdet;

    inv_mat[2][0] = (b[0]*c[1]-c[0]*b[1])*invdet;
    inv_mat[2][1] = (c[0]*a[1]-a[0]*c[1])*invdet;
    inv_mat[2][2] = (a[0]*b[1]-b[0]*a[1])*invdet;

    n_v[0] = v[0] * inv_mat[0][0] + v[1] * inv_mat[1][0] + v[2] * inv_mat[2][0];
    n_v[1] = v[0] * inv_mat[0][1] + v[1] * inv_mat[1][1] + v[2] * inv_mat[2][1];
    n_v[2] = v[0] * inv_mat[0][2] + v[1] * inv_mat[1][2] + v[2] * inv_mat[2][2];

    v[0] = n_v[0]; v[1] = n_v[1]; v[2] = n_v[2];


}


void calculate_inverse_matrix(double inv_mat[][3], double a[3], double b[3], double c[3]){
    for (int i = 0; i < 3; ++i){        
        inv_mat[i][0] = a[i];
        inv_mat[i][1] = b[i];
        inv_mat[i][2] = c[i];
    }       
}

void multiply_array_by_matrix(double v[3], double m[3][3]){
    double product[3];
    int i, j;
    product[0] = 0.0; product[1] = 0.0; product[2] = 0.0;

    for (i = 0; i < 3; i++){
        for (j = 0; j < 3; j++){
            //product[i] = product[i] + m[i][j] * v[j];
            product[i] = product[i] + m[j][i] * v[j];
        }
        
    }
    for (i = 0; i < 3; i++)
        v[i] = product[i];
}


double rand_normal(double mean, double stddev){
    
    //Box muller method

    static double n2 = 0.0;
    static int n2_cached = 0;


    if (!n2_cached)
    {
        double x, y, r;
        do
        {
            x = 2.0*genrand_real1() - 1;
            y = 2.0*genrand_real1() - 1;

            r = x*x + y*y;
        } while (r == 0.0 || r > 1.0);
        {
            double d = sqrt(-2.0*log(r) / r);
            double n1 = x*d;
            n2 = y*d;
            double result = n1*stddev + mean;
            n2_cached = 1;
            return result;
        }
    }
    else
    {
        n2_cached = 0;
        return n2*stddev + mean;
    }
}



string remove_ws(const std::string& str){
    std::string str_no_ws;
    for (char c : str) if (!std::isspace(c)) str_no_ws += c;
    return str_no_ws;
}
