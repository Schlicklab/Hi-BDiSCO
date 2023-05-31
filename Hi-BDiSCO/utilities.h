#ifndef UTILITIES
#define UTILITIES


#include <vector>
#include <iostream>
#include <algorithm> 
#include <functional>
#include <list>
#include <iterator>
#include <cmath>
#include <string.h>
#include <cctype>
#include "constants.h"
#include "mt.h"

using namespace std;


double norm(vector < double > &r1, vector < double > &r2);
double norm(vector < double > &r1, vector < double > &r2);
double norm(double r1[3], vector < double > &r2);
double norm(double r1[3], double r2[3]);
double norm(vector < double > r);
void transpose(vector < vector < double > > &mat);
vector <double> cross_product(vector < double > &r1, vector < double > &r2);
vector <double> cross_product(double r1[3], vector < double > &r2);
double dot_product(vector < double > &r1, vector < double > &r2);
double dot_product(vector < double > &r1, double r2[3]);
double dot_product(double r1[3], double r2[3]);
double dot_product(double r1[3], vector < double > &r2);
vector <double> multiply_vector_by_matrix(vector<double> v, vector< vector <double > > m);

void multiply_array_by_matrix(double v[3], vector< vector <double > > m);
void multiply_array_by_matrix(double v[3],  double m[3][3]);

void multiply_array_by_inverse_matrix(double v[3], double a[3], double b[3], double c[3]);
void calculate_inverse_matrix(double inv_mat[][3], double a[3], double b[3], double c[3]);

double normalize_angle(double angle);

double rand_normal(double mean, double stddev);

int scopyn(char *t, char *s, int n);
int scopyn2n(char *t, char *s, int n1, int n2);

string remove_ws(const std::string& str);

void Message();

#endif 
