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

// general functions
// -----------------

//  Overloaded << to print 3d vectors:
ostream& operator<<(ostream&, VECTOR3D);

//  Make movie:
void make_movie(int, vector<PARTICLE>&, double, double, double, ofstream&);

//  Function telling if the second column (1st number after 1st space) in the less than for input 'line1' than 'line2':
    //  Used for sorting purposes and selective output as a function of particle types in 'generate_inst_dump_files()'.
bool return_Lesser_Type(string, string);

//  Generate from a time-series movie file the instantaneous files for g(r) calculation.
void generate_inst_dump_files(int, int, int, int);

// initialize particle velocities
//void initialize_particle_velocities(vector<PARTICLE>&, double temperature);

#endif
