#include<iostream>
#include "vector3.h"
#include "fieldsgrids.h"
#include "parameters.h"

GridsPoints::GridsPoints( double px_in, double py_in, double pz_in,
                          double ex_in, double ey_in, double ez_in,
                          double bx_in, double by_in, double bz_in,
                          double vx_in, double vy_in, double vz_in,
                          double density_in)
{
//    ex = ex_in; ey = ey_in; ez = ez_in;
//    density = density_in;
//    FIJKtoXYZ( face_in, gi_in, gj_in, gk_in); // calculate px py pz
//    XYZtoB( pos3.x(), pos3.y(), pos3.z()); // calculate bx, by, bz

    pos3 = Vector3( px_in, py_in, pz_in);
    e3 =   Vector3( ex_in, ey_in, ez_in);
    b3 =   Vector3( bx_in, by_in, bz_in);
    v3 =   Vector3( vx_in, vy_in, vz_in);
    density = density_in;
//    face = 0; gi = 0; gj = 0; gk =0;
}

GridsPoints::GridsPoints(const GridsPoints& other)
{
    
    density = other.density;
    pos3.SetVector3( other.pos3);
    pos3 = Vector3( other.pos3);
    e3 = Vector3( other.e3);
    b3 = Vector3( other.b3);
    v3 = Vector3( other.v3);
    density = other.density;
}

GridsPoints::GridsPoints()
{
    pos3 = Vector3( 0.0, 0.0, 0.0);
    e3 =   Vector3( 0.0, 0.0, 0.0);
    b3 =   Vector3( 0.0, 0.0, 0.0);
    v3 =   Vector3( 0.0, 0.0, 0.0);
    density = 0.0;
//    face = 0; gi = 0; gj = 0; gk =0;
}