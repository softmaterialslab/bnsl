#ifndef _FUNCTIONS_H
#define _FUNCTIONS_H

#include<iostream>
#include<iomanip>
#include<fstream>
#include<vector>
#include<cmath>
#include "vector3d.h"
#include "particle.h"

using namespace std;

//  Overloaded << to print 3d vectors
ostream& operator<<(ostream&, VECTOR3D);

//  Fuction used for sorting purposes and select output as a function of particle types in 'generate_inst_dump_files()'
bool return_Lesser_Type(string, string);

//  Generate from a time-series movie file the instantaneous files for g(r) calculation
void generate_inst_dump_files(int, int, int, int);

#endif
