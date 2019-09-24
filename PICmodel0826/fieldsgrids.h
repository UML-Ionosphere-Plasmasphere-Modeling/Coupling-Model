#ifndef _FIELDSGRIDS_H_
#define _FIELDSGRIDS_H_
/// Field grids information
#include <iostream>
#include "parameters.h"
#include "vector3.h"

// set up grids container to contain class GridsPoints, the position 
// can be calculated but not stored
class GridsPoints
{
    public:
//************************************************************************
//************************************************************************
// initialize bx, by, bz using dipole field
//
//************************************************************************
//************************************************************************
inline void XYZtoB( Vector3 const& v)
{
    double r = sqrt(pow(v.x(),2.0) + pow(v.y(),2.0) + pow(v.z(),2.0));
    b3.Setx( 3 * dMoment * v.x() * v.z() / pow(r,5.0));
    b3.Sety( 3 * dMoment * v.y() * v.z() / pow(r,5.0));
    b3.Setz( dMoment * (3 * pow(v.z(),2.0) - pow(r,2.0)) / pow(r,5.0));    
}
//************************************************************************
//************************************************************************
// initialize vx, vy, vz using corotation assumption
// with pos3
//************************************************************************
//************************************************************************
inline void XYZtoVel( )
{
    Vector3 tempPos = pos3;
    Vector3 tempOmega = Vector3( 0.0, 0.0, omega_earth);
    tempPos.Setz( 0.0);
    v3 = tempOmega.CrossProduct(tempPos);
}
//************************************************************************
//************************************************************************
// initialize Ex, Ey, Ez using electron momentum equation
// with gradient Pe ( in r direction) and vel X B
// - grad(Pe) / Ni = mi *g , quasi-neutrality He = Hi 
//************************************************************************
//************************************************************************
inline void XYZtoE()
{
    Vector3 temp; // for calculating grad(Pe) term
    if( pos3.norm() > 0)
    {
    temp = pos3.NormalizedVector().ScaleProduct( mi0 * gravity);
    e3 = temp.PlusProduct( b3.CrossProduct(v3));
    }
}
//************************************************************************
//************************************************************************
// initialize density
//************************************************************************
//************************************************************************
inline void XYZtoDensity( )
{
    double scaleHeight = ikT / mi0 / gravity;
    if( pos3.norm() > 0)
    density = N0_i * mi0 * exp(-1 * (pos3.norm() - radius) / scaleHeight);              
}

//************************************************************************
//************************************************************************
// Initialization the pos3 for gridspoints
// giving face, i, j, k, all are int not int64
//************************************************************************
//************************************************************************

inline void InttoPos3( int face, int i, int j, int k)
{
    double px, py, pz;
    double temp[2];
    // 2.1 radial
    double L = LMin * pow(10, logRatio *  k );

    // 2.2 IgJg to ST note 0<ST<1
    temp[0] = (1.0 / fieldsGridsSize) * (i - 1);
    temp[1] = (1.0 / fieldsGridsSize) * (j - 1);
    // 2.3 ST to UV note -1<UV<1
    for ( int i=0; i<=1; i++)
    {
        if (temp[i] >= 0.5) 
        {   
            temp[i]= (1/3.) * (4*temp[i]*temp[i] - 1);
        }
        else
        {   
            temp[i]= (1/3.) * (1 - 4*(1-temp[i])*(1-temp[i]));
        }
    }
    // 2.4 UV to xyz 
    double kk = L * radius / sqrt(pow(1.0,2.0) + pow(temp[0],2.0) + pow(temp[1],2.0));
    switch (face)
    {
        case 0: px=1.0;           py=temp[0];     pz=temp[1]; break;
        case 1: px=-1.0*temp[0];  py=1.0;         pz=temp[1]; break;
        case 2: px=-1.0*temp[1];  py=temp[0];     pz=1.0;     break;
        case 3: px=-1.0;          py=-1.0*temp[0];pz=temp[1]; break;
        case 4: px=temp[0];       py=-1.0;        pz=temp[1]; break;
        default:px=temp[1];       py=temp[0];     pz=-1.0;    break;
    }
    px *= kk; py *= kk; pz *= kk;
    pos3 = Vector3( px, py, pz);

//    std::cout << L << " " << px << " " << py << " " << pz << std::endl;
}

//************************************************************************
//************************************************************************
// Reset parameters
//
//************************************************************************
//************************************************************************
inline void ResetParameters()
{
    density = 0.0;
    v3 = Vector3( 0.0, 0.0, 0.0);
}

//************************************************************************
//************************************************************************
// Calculate weighting of density on grids as well as velocity
// weight: iw, jw, kw
// double mass_in: weight of each simualtion particle
// Vector3 vp_in: velocity of each simulation particle
//************************************************************************
//************************************************************************
inline void UpdateDueToWgt( int iw, int jw, int kw, double mass_in, Vector3 vp_in)
{
    density += mass_in * iw * jw * kw / cellSize3; // acutally is mass not density
    v3 = v3.PlusProduct( Vector3( mass_in * vp_in.x() * iw * jw * kw / cellSize3,
                                  mass_in * vp_in.y() * iw * jw * kw / cellSize3,
                                  mass_in * vp_in.z() * iw * jw * kw / cellSize3
                                ));
}
// After all simulation are calculated once, the density means the total 
// mass at each grid points, and the v3 means the total momentum. 
inline void UpdateDueToWgt( GridsPoints***** ptrArray_in, double volume_in)
{
    if ( density != 0.0)
    {
    v3 = v3.ScaleProduct(1/density);
    density = density / volume_in;
    }
}
//************************************************************************
//************************************************************************
// Calculate curl B and curl E
//************************************************************************
//************************************************************************
inline Vector3 CurlB()
{
    return( Vector3(1.0, 1.0, 1.0));
}
inline Vector3 CurlE()
{
    return( Vector3(1.0, 1.0, 1.0));
}

// Calculate the velocity of electrons from Ampere's Law
inline Vector3 Velocity3e()
{
    return v3.MinusProduct( CurlB().ScaleProduct(1 / ( alpha * Ni)));
}

// Calculate div of Pe
inline Vector3 divPe()
{
    return( Vector3(1.0,1.0,1.0));
}

// update E from electron's momentum equation
// E = - Ve X B - grad Pe / N
// input (Ve) and (grad Pe), both are Vector3
inline void updateE( Vector3 Ve_in, Vector3 GradPe_in)
{
    e3 = b3.CrossProduct(Ve_in).MinusProduct(GradPe_in.ScaleProduct(1/density));
}

// update B from Faraday's Law
// B = B - tstep * (curl E)
// input (curl E)
void updateB( Vector3 E_in)
{
    b3 = b3.PlusProduct( E_in.ScaleProduct(-1 * tstep));
}

// return Vector3 e3
inline Vector3 E3()
{
    return e3;
}

// return Vector3 B3
inline Vector3 B3()
{ 
    return b3;
}

// return Vector3 pos3
inline Vector3 Pos3()
{
    return pos3;
}

// return density
inline double Density()
{
    return density;
}

// return velocity
inline Vector3 Vel3()
{
    return v3;
}

// return stopSign
inline int StopSign()
{
    return stopSign;
}

// set temperature
inline void SetTemperature( double temperature_in)
{
    temperature = temperature_in;
}

// set stopSign
inline void SetStopSign( int stopSign_in)
{
    stopSign = stopSign_in;
}

//////////////////////////    
    // Constructors
    GridsPoints( double px_in, double py_in, double pz_in,
                 double ex_in, double ey_in, double ez_in,
                 double bx_in, double by_in, double bz_in,
                 double vx_in, double vy_in, double vz_in,
                 double density_in, double temperature_in,
                 int stopSign_in);
                 
    GridsPoints( const GridsPoints& other);

    GridsPoints();

    
    private:
    Vector3 pos3;            // pos3: position for vector 3
    Vector3 e3;             // e3
    Vector3 b3;             // b3
    Vector3 v3;             // velocity3

    double density;
    double temperature;
    int stopSign;
//    int face; int gi; int gj; int gk; //face, i, j, in fieldsgrids, radial
};
#endif