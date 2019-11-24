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
// initialize vx, vy, vz of top boundary
//************************************************************************
//************************************************************************
inline void SetVel_topBoundary( const Vector3& vel_in)
{
    v3 = vel_in;
    ve3= vel_in;
}

//************************************************************************
//************************************************************************
// initialize Ex, Ey, Ez using electron momentum equation 
// - vel X B
//************************************************************************
//************************************************************************
inline void XYZtoE()
{
    Vector3 temp; // for calculating grad(Pe) term
    if( pos3.norm() > 0)
    {
//    temp = pos3.NormalizedVector().ScaleProduct( mi0_O * gravity);
//    e3 = temp.PlusProduct( b3.CrossProduct(v3));
    e3 = b3.CrossProduct(v3);
    }
}

//************************************************************************
//************************************************************************
// initialize density
//************************************************************************
//************************************************************************
inline void VeltoE_topBoundary()
{
    e3 = b3.CrossProduct(v3);
}
//************************************************************************
//************************************************************************
// initialize density
//************************************************************************
//************************************************************************
inline void XYZtoDensity( )
{
    double scaleHeight = ikT / mi0_H / gravity;
    if( pos3.norm() > 0)
    density_H = N0_H * exp(-1 * (pos3.norm() - radius) / scaleHeight);

    scaleHeight = ikT / mi0_He / gravity;
    if( pos3.norm() > 0)
    density_He= N0_He * exp(-1 * (pos3.norm() - radius) / scaleHeight);
    
    scaleHeight = ikT / mi0_O / gravity;
    if( pos3.norm() > 0)
    density_O= N0_O * exp(-1 * (pos3.norm() - radius) / scaleHeight);  
}
//************************************************************************
//************************************************************************
// initialize GradBNorm
//************************************************************************
//************************************************************************
inline void XYZtoGradBNorm()
{
    double r_xyz, r_xy;
    double sintheta, costheta, sinphi, cosphi;
    double fact_1, fact_2;

    r_xyz = pos3.norm();
    r_xy = pow( pow( pos3.x(), 2.0) + pow( pos3.y(), 2.0), 0.5);
    sintheta = r_xy / r_xyz;
    costheta = pos3.z() / r_xyz;
    sinphi = pos3.y() / r_xy;
    cosphi = pos3.x() / r_xy;

    fact_1 = - 3.0 * pow( 1+ 3* pow(costheta,2.0), 0.5) / pow( r_xyz, 4.0);
    fact_2 = - 3.0 * sintheta * costheta / pow( r_xyz, 4.0) / pow( 1+ 3* pow(costheta,2.0), 0.5);
        
    gradB3.Setx( fact_1 * sintheta * cosphi + fact_2 * costheta * cosphi);
    gradB3.Sety( fact_1 * sintheta * sinphi + fact_2 * costheta * sinphi);
    gradB3.Setz( fact_1 * costheta - fact_2 * sintheta);

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
    density_H = 0.0;
    density_He = 0.0;
    density_O = 0.0;
    v3 = Vector3( 0.0, 0.0, 0.0);
    vH3 = Vector3( 0.0, 0.0, 0.0);
    vHe3 = Vector3( 0.0, 0.0, 0.0);
    vO3 = Vector3( 0.0, 0.0, 0.0);
}

//************************************************************************
//************************************************************************
// Calculate weighting of density on grids as well as velocity
// weight: iw, jw, kw
// double number_in: number weight of each simualtion particle
// Vector3 vp_in: velocity of each simulation particle
//************************************************************************
//************************************************************************
inline void UpdateDueToWgt( int iw, int jw, int kw, double number_in, Vector3 vp_in, int ionType_in)
{
    switch (ionType_in)
    {
    case 1: {
    density_H  += number_in * iw * jw * kw / cellSize3; // acutally is number not number density
    vH3 = vH3.PlusProduct( Vector3( number_in * vp_in.x() * iw * jw * kw / cellSize3,
                                     number_in * vp_in.y() * iw * jw * kw / cellSize3,
                                    number_in * vp_in.z() * iw * jw * kw / cellSize3 ));
    break;}
    case 4:{
    density_He += number_in * iw * jw * kw / cellSize3; // acutally is mass not density
    vHe3 = vHe3.PlusProduct( Vector3( number_in * vp_in.x() * iw * jw * kw / cellSize3,
                                     number_in * vp_in.y() * iw * jw * kw / cellSize3,
                                    number_in * vp_in.z() * iw * jw * kw / cellSize3 ));
    break;}
    default:{
    density_O  += number_in * iw * jw * kw / cellSize3; // acutally is mass not density
    vO3 = vO3.PlusProduct( Vector3( number_in * vp_in.x() * iw * jw * kw / cellSize3,
                                     number_in * vp_in.y() * iw * jw * kw / cellSize3,
                                    number_in * vp_in.z() * iw * jw * kw / cellSize3 ));
    break;}
    }
    
    v3 = v3.PlusProduct( Vector3( number_in * vp_in.x() * iw * jw * kw / cellSize3,
                                  number_in * vp_in.y() * iw * jw * kw / cellSize3,
                                  number_in * vp_in.z() * iw * jw * kw / cellSize3
                                ));
}
// After all simulation are calculated once, the density means the total 
// mass at each grid points, and the v3 means the total momentum. 
inline void UpdateDueToWgt( GridsPoints***** ptrArray_in, double volume_in)
{
    if ( density_H != 0 || density_He != 0 || density_O != 0)
    {
        v3 = v3.ScaleProduct(1.0 / (density_H + density_He + density_O) / updateInfoPeriod);
       
        if( density_H != 0)
        {
            vH3= vH3.ScaleProduct(1.0 / density_H / updateInfoPeriod);
        } else
        {
            vH3 = Vector3( 0.0, 0.0, 0.0);
        }

        if( density_He !=0)
        {
        vHe3= vHe3.ScaleProduct(1.0 / density_He / updateInfoPeriod);
        } else
        {
            vHe3 = Vector3( 0.0, 0.0, 0.0);
        }

        if( density_O !=0)
        {
            vO3= vO3.ScaleProduct(1.0 / density_O / updateInfoPeriod);
        } else
        {
            vO3 = Vector3( 0.0, 0.0, 0.0);
        }
        
        density_H = density_H / volume_in / updateInfoPeriod;
        density_He = density_He / volume_in / updateInfoPeriod;
        density_O = density_O / volume_in / updateInfoPeriod;
    } else
    {
        v3 = Vector3( 0.0, 0.0, 0.0);
    }
}

// update E from electron's momentum equation
// E = - Ve X B - grad Pe / N / qi0
// input grad Pe at grids

inline void updateE( Vector3 GradPe_in)
{
    if( update_type == 0)
    {
    e3 = b3.CrossProduct(v3).MinusProduct(GradPe_in.ScaleProduct(1.0 / (density_H + density_He + density_O ) / qi0));
    } else
    {
    e3 = b3.CrossProduct(ve3).MinusProduct(GradPe_in.ScaleProduct(1.0 / (density_H + density_He + density_O ) / qi0));
    }
    
}

// Update ve from ampere's law
// ve = vi - curl B / mu0 / qi0 / total number density
// input curl B
inline void updateve3( Vector3 curlB_in)
{
    ve3 = v3.MinusProduct( curlB_in.ScaleProduct( 1/ (density_H + density_He + density_O ) / qi0 / mu0));
}

// update B from Faraday's Law
// B = B - tstep * (curl E)
// input (curl E)
void updatedB( Vector3 E_in)
{

    dB3 = dB3.PlusProduct( E_in.ScaleProduct(-1 * tstep * updateInfoPeriod));
}





// return Vector3 gradB3
inline Vector3 GradB3()
{
    return gradB3;
}

// return Vector3 e3
inline Vector3 E3()
{
    return e3;
}
// return Vector3 b3 only
inline Vector3 B3_base()
{
    return b3;
}
// return Vector3 B3
inline Vector3 B3()
{ 
    return b3.PlusProduct(dB3);
}

// return Vector3 dB3
inline Vector3 DB3()
{
    return dB3;
}

// return Vector3 pos3
inline Vector3 Pos3()
{
    return pos3;
}

// return density
inline double Density()
{
    return density_H + density_He + density_O;
}

// return density of H
inline double Density_H()
{
    return density_H;
}

// return density of H
inline double Density_He()
{
    return density_He;
}

// return density of H
inline double Density_O()
{
    return density_O;
}
// set density 
inline double Density_H( double den_in)
{
    density_H = den_in;
}
inline double Density_He( double den_in)
{
    density_He = den_in;
}
inline double Density_O( double den_in)
{
    density_O = den_in;
}
// return velocity
inline Vector3 Vel3()
{
    return v3;
}

inline Vector3 VelH3()
{
    return vH3;
}

inline Vector3 VelHe3()
{
    return vHe3;
}

inline Vector3 VelO3()
{
    return vO3;
}

// return velocity e3
inline Vector3 Vel_e3()
{
    return ve3;
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
inline double Temperature( )
{
    return temperature;
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
                 double dBx_in,double dBy_in,double dBz_in,
                 double vx_in, double vy_in, double vz_in,
                 double vHx_in, double vHy_in, double vHz_in,
                 double vHex_in, double vHey_in, double vHez_in,
                 double vOx_in, double vOy_in, double vOz_in,
                 double vex_in, double vey_in, double vez_in,
                 double gradBx_in, double gradBy_in, double gradBz_in,
                 double density_H_in, double density_He_in, double density_O_in,
                 double temperature_in,
                 int stopSign_in);
                 
    GridsPoints( const GridsPoints& other);

    GridsPoints();

    
    private:
    Vector3 pos3;            // pos3: position for vector 3
    Vector3 e3;             // e3

    Vector3 b3;             // b3
    Vector3 dB3;            // perturbation of b3
    
    Vector3 v3;             // average ions velocity3
    Vector3 vH3;
    Vector3 vHe3;
    Vector3 vO3;

    Vector3 ve3;            // e velocity

    Vector3 gradB3;          // gradB in which B is norm of vector3 B

    double density_H;    // mass density
    double density_He;   // mass density
    double density_O;    // mass density

    double temperature;     // electron temperature: Te
    int stopSign;
//    int face; int gi; int gj; int gk; //face, i, j, in fieldsgrids, radial
};
#endif