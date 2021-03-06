#ifndef __UTIL_H__
#define __UTIL_H__

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include "global.h"
using namespace std;

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)
#define AT __FILE__ ":" TOSTRING(__LINE__)
#define report_exit(a) _report_exit(AT,a)
#define fzero(a) (fabs(a)<FLOAT_ZERO)
#define feqn(a,b) (fzero((a)-(b)))

void _report_exit(const char *location, const char *msg);
double ** new_square_matrix(int n);
void delete_matrix(double **,int);
void output_matrix(double **,int);
string get_basename(char * filename);

void open_logfile(const char * logname);
void close_logfile();

DIRECTION get_opposite_dir(DIRECTION dir);

// given a vector, copy its element to a basic array
template<class T>
void vector_to_array(vector<T> v, T * arr){
	copy(v.begin(), v.end(), arr);
}

#endif
