// This is the particle class
// It provides features to the particle and also lists how their positions and velocities are updated (in case this is needed)

#ifndef _PARTICLE_H
#define _PARTICLE_H

#include "vector3d.h"

class PARTICLE 
{
  public:

  // members
  int id;			// id of the particle
  double diameter;		// diameter of the particle
  double m; 			// mass of the particle
  double charge;
  VECTOR3D posvec;		// position vector of the particle
  VECTOR3D velvec;		// velocity vector of the particle
  VECTOR3D forvec;		// force vector on the particle
  VECTOR3D initPosVec; 		// the initial position of the particle (for MSD code only)
  double pe;			// potential energy
  long double ke;		// kinetic energy
  double energy;		// energy
  long double lx, ly, lz;	// box lengths-- useful information for particles to know, hence here

  // member functions
  
  // make a particle		// this defines the constructor
  PARTICLE(int id = 0, double diameter = 0, double mass = 0, double q = 1, VECTOR3D position = VECTOR3D(0,0,0), VECTOR3D initPosition = VECTOR3D(0,0,0), long double lx = 0, long double ly = 0, long double lz = 0): id(id), diameter(diameter), m(mass), charge(q), posvec(position), initPosVec(initPosition), lx(lx), ly(ly), lz(lz)
  {
  }
  
  // the next two functions are central to the velocity-Verlet algorithm
  // update position of the particle
  void update_position(double dt)		// dt is taken from outside; it is the time-step that the user supplies (based on steps)
  {
    posvec = ( posvec + (velvec ^ dt) );	// position updated to a full time-step
    // periodic boundary
    if (posvec.x > lx/2.0)
      posvec.x = posvec.x - lx;
    if (posvec.x < -lx/2.0)
      posvec.x = posvec.x + lx;
    if (posvec.y > ly/2.0)
      posvec.y = posvec.y - ly;
    if (posvec.y < -ly/2.0)
      posvec.y = posvec.y + ly;
    //if (posvec.z > lz/2.0)
      //posvec.z = posvec.z - lz;
    //if (posvec.z < -lz/2.0)
      //posvec.z = posvec.z + lz;
    return;
  }
  
  // update velocity of the particle
  void update_velocity(double dt)	
  {
    velvec = ( velvec + ( (forvec) ^ ( 0.5 * dt / m ) ) );	// notice the half time-step
    return;
  }
  
  // calculate kinetic energy of a particle
  void kinetic_energy()				
  {
    ke = 0.5 * m * velvec.GetMagnitude() * velvec.GetMagnitude();
    return;
  }
};

#endif
