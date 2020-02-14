#ifndef _VECTOR3_H_
#define _VECTOR3_H_

#include<iostream>
#include "mathutil.h"
#include "parameters.h"
#include "structdef.h"

class Vector3
{
    public:
    inline void printv3()
    {
        std::cout << v_x << " " << v_y << " " << v_z << std:: endl;
    }
    
    inline void Setx( const double &v) {v_x = v;}
    
    inline void Sety( const double &v) {v_y = v;}
    
    inline void Setz( const double &v) {v_z = v;}

    inline double x() const { return v_x;}

    inline double y() const { return v_y;}

    inline double z() const { return v_z;}

// Calculate V1 X V2 , return a Vector3
    inline Vector3 CrossProduct(const Vector3& v3b)
    {
        return( Vector3(v_y * v3b.v_z - v_z * v3b.v_y,
                        v_z * v3b.v_x - v_x * v3b.v_z,
                        v_x * v3b.v_y - v_y * v3b.v_x));
    }

// Calculate V1 dot V2, return a double
    inline double DotProduct( const Vector3& v3b)
    {
        return( v_x * v3b.v_x + v_y * v3b.v_y + v_z * v3b.v_z);
    }

// Calculate a times V, return a Vector3
    inline Vector3 ScaleProduct( const double& scale)
    {
        return( Vector3( v_x* scale, v_y* scale, v_z* scale));
    }
// Calculate mix product a dot b cross c
    inline double MixProduct( const Vector3& v3b, const Vector3& v3c)
    {
        return  (v_y * v3b.v_z - v_z * v3b.v_y) * v3c.v_x +
                (v_z * v3b.v_x - v_x * v3b.v_z) * v3c.v_y +
                (v_x * v3b.v_y - v_y * v3b.v_x) * v3c.v_z;
    }
// Calculate V1 + V2, return a Vector3
    inline Vector3 PlusProduct( const Vector3& v3b)
    {
        double a = v_x+v3b.x();
        double b = v_y+v3b.y();
        double c = v_z+v3b.z();
        return( Vector3( a, b, c));
    }

// Calculate V1 - V2, return a Vector3
    inline Vector3 MinusProduct( const Vector3& v3b)
    {
        double a = v_x-v3b.x();
        double b = v_y-v3b.y();
        double c = v_z-v3b.z();
        return( Vector3( a, b, c));
    }

// Calculate the |v|
inline double norm()
{
    return sqrt(v_x * v_x + v_y * v_y + v_z * v_z);   
}

// Calculate the v^2
inline double norm2()
{
    return v_x * v_x + v_y * v_y + v_z * v_z;   
}

// Calculate unit vector of self
inline Vector3 NormalizedVector()
{
    double r = norm();
    return ScaleProduct(1/r);
}

// Calculate unit Vector orthogonal to a plane constructed by two vectors
// Take care about the direction
inline Vector3 OrthoUnitVector( const Vector3& v3b)
{
    return CrossProduct( v3b).NormalizedVector();
}

// Set value

inline void SetVector3( const Vector3& v3b)
{
    v_x = v3b.v_x;
    v_y = v3b.v_y;
    v_z = v3b.v_z;
}

// Show the Vector3
inline Vector3 V3()
{
    return Vector3(v_x, v_y, v_z);
}

// return the Uint for vector3 assume it is a position vector
inline uint_64 Uint_64_Trans()
{   
    // fixed from particles.cpp Uintupdate_64()
    double px = v_x;
    double py = v_y;
    double pz = v_z;
    double temp[2];
    // 4.1 radial kp
    double L = sqrt(pow(px,2.0)+pow(py,2.0)+pow(pz,2.0))/radius;
    uint_64 kp = static_cast<uint_64>(floor(log10(L / LMin) / logRatio *cellSize1) );
    // 4.2 XYZtoUV, note that -1<UV<1 
    uint_64 face = Getface(px, py, pz);
    switch (face)
    {
        case 0: temp[0] = py/px; temp[1] = pz/px; break;
        case 1: temp[0] =-px/py; temp[1] = pz/py; break;
        case 2: temp[0] = py/pz; temp[1] =-px/pz; break;
        case 3: temp[0] = py/px; temp[1] =-pz/px; break;
        case 4: temp[0] =-px/py; temp[1] =-pz/py; break;
        default:temp[0] =-py/pz; temp[1] =-px/pz; break;
    }
    // 4.3 UVtoST, note that 0<ST<1
    for (int i=0; i<=1; i++)
    {
        if (temp[i] >= 0) temp[i] = 0.5 * std::sqrt(1 + 3*temp[i]);
        else            temp[i] = 1 - 0.5 * std::sqrt(1 - 3*temp[i]);
    }
    // 4.4 STtoIpJp // Notice the structure of grids, main domain not from zero
    uint_64 ip= static_cast<unsigned int>(floor(temp[0] * particlesGridsSize) + 1);
    uint_64 jp= static_cast<unsigned int>(floor(temp[1] * particlesGridsSize) + 1);
    
    // 5. F ip jp kp to Uint_64
    uint_64 posUint = face << 61 ;
    for( int i = 0; i < particlesGridsLevel; i++)
    {
        posUint += (((ip >> particlesGridsLevel-1-i) & 1 )<< 60 - i *3) 
                    + (((jp >> particlesGridsLevel-1-i) & 1 )<< 60-1 - i *3)
                    + (((kp >> particlesGridsLevel-1-i) & 1 )<< 60-2 - i *3) ;
    }
    return posUint;
}

// Calculate the face B vector from a linear equation
// The equation is like this
// (vector v1)    (x)      (c1)
// (vector v2) *  (y)  =   (c2)
// (vector v3)    (z)      (c3)
// Solve for ( x, y ,z) and the v is the face area vector, and c1 is the integration
// of E along the completed loop for outside direction
inline Vector3 FaceBSolver( const Vector3& v1, const Vector3& v2, const Vector3& v3,
                            const double& c1, const double& c2, const double& c3)
                            {
                              
                                // computes the inverse of a matrix m
                                double det = v1.x() * (v2.y() * v3.z() - v3.y() * v2.z()) -
                                             v1.y() * (v2.x() * v3.z() - v2.z() * v3.x()) +
                                             v1.z() * (v2.x() * v3.y() - v2.y() * v3.x());
                                if(det == 0)
                                std::cout << " Face B solver DET can't be zero" << std::endl;
                                double invdet = 1.0 / det;

                                double minv[3][3]; // inverse of matrix m
                                minv[0][0] = (v2.y() * v3.z() - v3.y() * v2.z()) * invdet;
                                minv[0][1] = (v1.z() * v3.y() - v1.y() * v3.z()) * invdet;
                                minv[0][2] = (v1.y() * v2.z() - v1.z() * v2.y()) * invdet;
                                minv[1][0] = (v2.z() * v3.x() - v2.x() * v3.z()) * invdet;
                                minv[1][1] = (v1.x() * v3.z() - v1.z() * v3.x()) * invdet;
                                minv[1][2] = (v2.x() * v1.z() - v1.x() * v2.z()) * invdet;
                                minv[2][0] = (v2.x() * v3.y() - v3.x() * v2.y()) * invdet;
                                minv[2][1] = (v3.x() * v1.y() - v1.x() * v3.y()) * invdet;
                                minv[2][2] = (v1.x() * v2.y() - v2.x() * v1.y()) * invdet;

                                v_x = minv[0][0]*c1 + minv[0][1]*c2 + minv[0][2]*c3;
                                v_y = minv[1][0]*c1 + minv[1][1]*c2 + minv[1][2]*c3;
                                v_z = minv[2][0]*c1 + minv[2][1]*c2 + minv[2][2]*c3;
                            }
// constructor
    Vector3()
    {
        v_x=0.0;
        v_y=0.0;
        v_z=0.0;
    }
    Vector3(double x, double y, double z) 
    {
        v_x = x;
        v_y = y;
        v_z = z;
    };
    Vector3( const Vector3& v3b)
    {
        v_x = v3b.v_x;
        v_y = v3b.v_y;
        v_z = v3b.v_z;
    }

    private:
    double v_x;
    double v_y;
    double v_z;
};

#endif


//***************************************************************************
