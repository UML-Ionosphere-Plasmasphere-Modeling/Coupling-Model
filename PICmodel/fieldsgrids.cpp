#include<iostream>
#include "vector3.h"
#include "fieldsgrids.h"
#include "parameters.h"

GridsPoints::GridsPoints( double px_in, double py_in, double pz_in,
                          double ex_in, double ey_in, double ez_in,
                          double bx_in, double by_in, double bz_in,
<<<<<<< HEAD:PICmodel/fieldsgrids.cpp
                          double vx_in, double vy_in, double vz_in,
                          double density_in)
=======
                          double dBx_in,double dBy_in,double dBz_in,
                          double vx_in, double vy_in, double vz_in,
                          double vex_in, double vey_in, double vez_in,
                          double gradBx_in, double gradBy_in, double gradBz_in,
                          double density_H_in, double density_He_in, double density_O_in,
                          double temperature_in,
                          int stopSign_in)
>>>>>>> 865083259b6dad8f4208c771b239bab419e7fcb4:PICmodel/fieldsgrids.cpp
{
//    ex = ex_in; ey = ey_in; ez = ez_in;
//    density = density_in;
//    FIJKtoXYZ( face_in, gi_in, gj_in, gk_in); // calculate px py pz
//    XYZtoB( pos3.x(), pos3.y(), pos3.z()); // calculate bx, by, bz

    pos3 = Vector3( px_in, py_in, pz_in);
    e3 =   Vector3( ex_in, ey_in, ez_in);
    b3 =   Vector3( bx_in, by_in, bz_in);
<<<<<<< HEAD:PICmodel/fieldsgrids.cpp
    v3 =   Vector3( vx_in, vy_in, vz_in);
    density = density_in;
=======
    dB3=   Vector3( dBx_in, dBy_in, dBz_in);
    v3 =   Vector3( vx_in, vy_in, vz_in);
    ve3=   Vector3( vex_in, vey_in, vez_in);
    gradB3 = Vector3( gradBx_in, gradBy_in, gradBz_in);

    density_H = density_H_in;
    density_He= density_He_in;
    density_O = density_O_in;
    
    temperature = temperature_in;
    stopSign = stopSign_in;
>>>>>>> 865083259b6dad8f4208c771b239bab419e7fcb4:PICmodel/fieldsgrids.cpp
//    face = 0; gi = 0; gj = 0; gk =0;
}

GridsPoints::GridsPoints(const GridsPoints& other)
{
<<<<<<< HEAD:PICmodel/fieldsgrids.cpp
    
    density = other.density;
=======
>>>>>>> 865083259b6dad8f4208c771b239bab419e7fcb4:PICmodel/fieldsgrids.cpp
    pos3.SetVector3( other.pos3);
    pos3 = Vector3( other.pos3);
    e3 = Vector3( other.e3);
    b3 = Vector3( other.b3);
<<<<<<< HEAD:PICmodel/fieldsgrids.cpp
    v3 = Vector3( other.v3);
    density = other.density;
=======
    dB3= Vector3( other.dB3);
    v3 = Vector3( other.v3);
    ve3= Vector3( other.ve3);
    gradB3 = Vector3( other.gradB3);

    density_H = other.density_H;
    density_He = other.density_He;
    density_O = other.density_O;

    temperature = other.temperature;
    stopSign = other.stopSign;
>>>>>>> 865083259b6dad8f4208c771b239bab419e7fcb4:PICmodel/fieldsgrids.cpp
}

GridsPoints::GridsPoints()
{
    pos3 = Vector3( 0.0, 0.0, 0.0);
    e3 =   Vector3( 0.0, 0.0, 0.0);
    b3 =   Vector3( 0.0, 0.0, 0.0);
<<<<<<< HEAD:PICmodel/fieldsgrids.cpp
    v3 =   Vector3( 0.0, 0.0, 0.0);
    density = 0.0;
=======
    dB3=   Vector3( 0.0, 0.0, 0.0);
    v3 =   Vector3( 0.0, 0.0, 0.0);
    ve3=   Vector3( 0.0, 0.0, 0.0);
    gradB3= Vector3( 0.0, 0.0, 0.0);

    density_H = 0.0;
    density_He= 0.0;
    density_O = 0.0;

    temperature = 0.0;
    stopSign= 0;
>>>>>>> 865083259b6dad8f4208c771b239bab419e7fcb4:PICmodel/fieldsgrids.cpp
//    face = 0; gi = 0; gj = 0; gk =0;
}