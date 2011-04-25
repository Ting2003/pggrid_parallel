// ----------------------------------------------------------------//
// Filename : util.cpp
// Author : Xiao Zigang <zxiao2@illinois.edu>
//
// Source code for some utility functions
// ----------------------------------------------------------------//

#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>
#include "global.h"
#include "util.h"
using namespace std;

// report an error string and exit the program
// Note that it is wrapped by a macro `report_exit'
// Do NOT use this directly.
void _report_exit(const char *location, const char *msg){
#ifdef DEBUG
	  fprintf(stderr,"Error at %s: %s", location, msg);
#else
	  fprintf(stderr,"%s", msg);
#endif
	  exit(1);
}

// dynamically allocate a 2D square array of size n
double ** new_square_matrix(int n){
	double **m = new double *[n];
	int i,j;
	for(i=0;i<n;i++){
		m[i] = new double [n];
		for(j=0;j<n;j++) // initialize as 0
			m[i][j]=0.0;
	}
	return m;
}

// delete the allocated space from new_square_matrix
void delete_matrix(double ** m, int n){
	for(int i=0;i<n;i++) delete m[i];
	delete [] m;
}

// give a 2D array, output it
void output_matrix(double ** m,int n){
	int precision = 4;
	int width = precision+8;
	cout.precision(precision);
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			cout<<scientific<<setw(width)<<m[i][j];
		}
		cout<<endl;
	}
	cout<<endl;
}

string get_basename(char * filename){
	string fn(filename);
	size_t pos = fn.find(".");
	return fn.substr(0,pos);
}

streambuf * clog_save=NULL;
ofstream logstream;
void open_logfile(const char * logname){
	if( logname == NULL ) return;
	clog_save = clog.rdbuf();
	logstream.open(logname);
	clog.rdbuf(logstream.rdbuf());
}

void close_logfile(){
	if(clog_save == NULL) return;
	clog.rdbuf(clog_save);
	logstream.close();
}

DIRECTION get_opposite_dir(DIRECTION dir){
	DIRECTION ret;
	switch(dir){
	case NORTH: ret = SOUTH; break;
	case SOUTH: ret = NORTH; break;
	case EAST: ret = WEST; break;
	case WEST: ret = EAST; break;
	case TOP: ret = BOTTOM; break;
	case BOTTOM: ret = TOP; break;
	case UNDEFINED: report_exit("Invalid usage of get_opposite_dir"); break;
	}
	return ret;
}
