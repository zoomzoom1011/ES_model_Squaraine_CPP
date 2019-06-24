#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <cstring>
#include <time.h>
#include <cstdlib>
#include <math.h>
#include <vector>
#include <iomanip>
#include <unistd.h>
#include <algorithm>    // std::min
#include "mkl_lapacke.h"
#include <sys/stat.h>

using namespace std;

// defining the global variable
string task_title;
int vibmax, sys_vibmax;
double hw;
bool periodic;
string xyz_file[2];
double lambda_n, lambda_z1, lambda_z2, lambda_a, lambda_c;
double abs_lw;
double nz, tz, nct, tct;
int spec_step;
double spec_start_ab, spec_end_ab;
double spec_start_pl, spec_end_pl;
bool lorentzian;
double agg_angle;
double dielectric;
bool no_frenkel;
bool calc_pl, nearest_neighbor;

// constant
double pi = atan(1.00) * 4;

// system parameters
const int es_zmax = 3;
const int da = 3;
int nmax = 2;
const int es_n = 1, es_z1 = 2, es_z2 = 3;
int mon_kount, sys_kount;
double *sys_h = NULL;

// geometry
double es_lambda_z1[3] = { 0.0,1.0,0.0 }; //N, Z1, Z2 for arm z1
double es_lambda_z2[3] = { 0.0,0.0,1.0 }; //N, Z1, Z2 for arm z2, shift if in the specified electronic state

// z_state parts
const int leftdonor = 1;
const int acceptor = 2;
const int rightdonor = 3;

// class method
class basis
{
public:
    int es_state;
    int vib1;
    int vib2;
};

basis *mon_state = NULL;
basis **sys_state = NULL;
// tdm
double* ux = NULL;
double* uy = NULL;
double* uz = NULL;

// absorption spectra
double* ab_sys_eval; 
bool abs_freq_dep; 
double* ab_osc_x = NULL;
double* ab_osc_y = NULL;
double* ab_osc_z = NULL;      
double* ab_x = NULL;   
double* ab_y = NULL;   
double* ab_z = NULL;    

// pl spectra
double* pl_sys_eval; 
double* pl_osc_x = NULL;
double* pl_osc_y = NULL;
double* pl_osc_z = NULL;
double* pl_x = NULL;  
double* pl_y = NULL;  
double* pl_z = NULL;   
int pl_start; 

// geometry
double* coulomb_coupling = NULL;
double* mol1pos = NULL;
int anum = 3;   // SQ dye