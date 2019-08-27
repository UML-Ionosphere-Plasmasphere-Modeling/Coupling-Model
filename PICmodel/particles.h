#ifndef _PARTICLES_H_
#define _PARTICLES_H_
/// Particle information
#include <iostream>
#include "mathutil.h"
#include "parameters.h"
#include "vector3.h"
#include "structdef.h"
#include "fieldsgrids.h"

class Particles
{
public:
    friend class GridsPoints;
    friend class Vector3;


//************************************************************************
//************************************************************************
// Transfore uint_64 posuInt to F I J K of the grids, from which we can latey determin the total
// 8 grids-points of the cell in which the particles is.

//************************************************************************
//************************************************************************
inline structg InttoStrp1()
{
    struct structg strg = {0,0,0,0,0,0,0, 0.0, 0.0, 0.0, 0.0};
    strg.face = posUint >> 61;
//std::cout << posUint << " "<< strg.face <<" "<< strg.ig<<" "<< strg.jg<<" "<< strg.kg << std::endl;
    for( int i = 0; i < fieldsGridsLevel; i++) 
    {
        strg.ig = (strg.ig << 1) + ((posUint >> 60   - i*3) & 1);
        strg.jg = (strg.jg << 1) + ((posUint >> 60-1 - i*3) & 1);
        strg.kg = (strg.kg << 1) + ((posUint >> 60-2 - i*3) & 1);
    }
//std::cout << posUint << " "<< strg.face <<" "<< strg.ig<<" "<< strg.jg<<" "<< strg.kg << std::endl;

    for( int i = fieldsGridsLevel; i < particlesGridsLevel; i++)
    {
        strg.iw = (strg.iw << 1) + ((posUint >> 60   - i*3) & 1);
        strg.jw = (strg.jw << 1) + ((posUint >> 60-1 - i*3) & 1);
        strg.kw = (strg.kw << 1) + ((posUint >> 60-2 - i*3) & 1);
    }
    strg.vx = vp.x(); strg.vy = vp.y(); strg.vz = vp.z();
    strg.mass = mass1;
    return strg;
}
//************************************************************************
//************************************************************************
// Transfor uint_64 posuInt to i j k of the particle in the cell, which can help to determin the 
// weighting of particle density and velocity on the grids-points of the cell in which the particles
// is.
//************************************************************************
//************************************************************************
inline structg InttoStrp2()
{
    struct structg strg;
    strg.face = posUint >> 61;
    for( int i = 0; i < fieldsGridsSize; i++) 
    {
        strg.ig = (strg.ig << 1) + ((posUint >> 60   - i*3) & 1);
        strg.jg = (strg.jg << 1) + ((posUint >> 60-1 - i*3) & 1);
        strg.kg = (strg.kg << 1) + ((posUint >> 60-2 - i*3) & 1);
    }
    for( int i = fieldsGridsSize+1; i < particlesGridsSize; i++)
    {
        strg.iw = (strg.iw << 1) + ((posUint >> 60   - i*3) & 1);
        strg.jw = (strg.jw << 1) + ((posUint >> 60-1 - i*3) & 1);
        strg.kw = (strg.kw << 1) + ((posUint >> 60-2 - i*3) & 1);
    }
    strg.vx = vp.x(); strg.vy = vp.y(); strg.vz = vp.z();
    strg.mass = mass2;
    return strg;
}


//************************************************************************
//************************************************************************
// Calculate the local E and B of the particles, return vector3 E or B
//
//************************************************************************
//************************************************************************
Vector3 LocalE( struct structg *strg_in, GridsPoints***** ptrArray_in);
Vector3 LocalB( struct structg *strg_in, GridsPoints***** ptrArray_in);

//************************************************************************
//************************************************************************
// Calculate the new vp by Boris' Method, update vector3 vp
// And return a int "0" means in the main domain
// "1" means out of the main domain
//
//************************************************************************
//************************************************************************
int BorisMethod( struct structg *strg_in, GridsPoints***** ptrArray_in);

//************************************************************************
//************************************************************************
// Calculate the new position of the particles in uint_64
// And return a int "0" means in the main domain
// "1" means out of the main domain
//
//************************************************************************
//************************************************************************
int UpdateUint_64();

inline Vector3 VelParticles()
{
    return vp;
}

inline double WeightMi()
{
    return mi;
}

//////////////////////////////////Constructor//////////////////////////////

    Particles( uint_64 posInt_in, Vector3 vx_in, double mi_in);
    Particles();

private:
    uint_64 posUint; // single unsigned int for position
    Vector3 vp; // velocity of particles
    double mi; // weight
    static constexpr double mass1 = 16;
    static constexpr double mass2 = 1;
//    double px; double py; double pz;
//    double vx; double vy; double vz;
//    int face; int ri; int rj; int rk; //face, (i,j) in fieldsgrids, radial
};
#endif