#ifndef _MODULE_H_1_
#define _MODULE_H_1_
#include<iostream>
#include "parameters.h"
#include "module_1.h"
#include <memory>
#include "particles.h"
#include <vector>
#include "fieldsgrids.h"
#include <cmath>
#include <limits>
#include <bitset>

using std::vector;

//************************************************************************
//************************************************************************
// FUNCTION
// Initialization the particles for velocity and position

vector<Particles>* ParticlesLists( GridsPoints***** ptrArray_in, double*** ptrVolumeCellArray_in, double mi0, double N0);

//************************************************************************
//************************************************************************
// FUNCTION
// Initialization the particles for velocity and position 
// Generate lists of particles for bot and top region temp

vector<Particles>* ParticlesListsTemp( GridsPoints***** ptrArray_in, double*** ptrVolumeCellArray_in, double mi0, int ionType_in);



//************************************************************************
//************************************************************************
// FUNCTION
// Transform the ubit64 to vector of position

Vector3 Uint64ToVector3( uint_64 intPos_in);


//************************************************************************
//************************************************************************
// FUNCTION rand (0 - 1)
inline double dRand()
{
    return (double) rand() / (RAND_MAX);
}

//************************************************************************
//************************************************************************
// FUNCTION return uniform distrubtion random Vector3
// in cell (face, i, j, k).
// Remind that there are limited position for particles to locate, we just
// need to random the position in smaller cells
inline uint_64 UniDisInCell( GridsPoints***** ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    // 1. Randomly get three ints in 3 direction
    uint_64 ipar = static_cast<uint_64>( floor(cellSize1 * dRand())) ;
    uint_64 jpar = static_cast<uint_64>( floor(cellSize1 * dRand()));   
    uint_64 kpar = static_cast<uint_64>( floor(cellSize1 * dRand()));
    // 2. Transfer the base point from grids location for posUint calculation
    uint_64 face = face_in ;
    uint_64 ig = i_in -1 ;
    uint_64 jg = j_in -1 ;
    uint_64 kg = k_in;
   
//    std::cout << ipar << " " << jpar << " " << kpar << " " << std::endl;

    // 3. Transfer the location in int to the location in Uint_64 
    uint_64 posUint = face << 61;
    for( int i = 0; i < fieldsGridsLevel; i++)
    {
        posUint += (((ig >> fieldsGridsLevel-1-i) & 1 )<< 60 - i *3) 
                    + (((jg >> fieldsGridsLevel-1-i) & 1 )<< 60-1 - i *3)
                    + (((kg >> fieldsGridsLevel-1-i) & 1 )<< 60-2 - i *3) ;
    }


//    std::cout << "1    " << std::bitset<64>(posUint) << std::endl; 
    int cellLevel = particlesGridsLevel - fieldsGridsLevel;
    for( int i = fieldsGridsLevel; i < particlesGridsLevel; i++)
    {  
        int j = i - fieldsGridsLevel;
        posUint += (((ipar >> cellLevel-1- j) & 1 )<< 60 - i *3) 
                    + (((jpar >> cellLevel-1- j) & 1 )<< 60-1 - i *3)
                    + (((kpar >> cellLevel-1- j) & 1 )<< 60-2 - i *3) ;
    }
/*
    std::cout << "2    " << std::bitset<64>(posUint) << std::endl; 

    struct structg strg = {0,0,0,0,0,0,0, 0.0, 0.0, 0.0, 0.0};
    strg.face = posUint >> 61;
    for( int i = 0; i < fieldsGridsLevel; i++) 
    {
        strg.ig = (strg.ig << 1) + ((posUint >> 60   - i*3) & 1);
        strg.jg = (strg.jg << 1) + ((posUint >> 60-1 - i*3) & 1);
        strg.kg = (strg.kg << 1) + ((posUint >> 60-2 - i*3) & 1);
    }

    std::cout << "3    " << std::bitset<64>(posUint) << std::endl; 

    for( int i = fieldsGridsLevel; i < particlesGridsLevel; i++)
    {
        strg.iw = (strg.iw << 1) + ((posUint >> 60   - i*3) & 1);
        strg.jw = (strg.jw << 1) + ((posUint >> 60-1 - i*3) & 1);
        strg.kw = (strg.kw << 1) + ((posUint >> 60-2 - i*3) & 1);
    }

    
   std::cout <<"face " << std::bitset<64>(strg.face) << std::endl;
    std::cout << "ig   " <<std::bitset<64>(strg.ig) << std::endl;
    std::cout <<"jg   " << std::bitset<64>(strg.jg) << std::endl;
    std::cout << "kg   " <<std::bitset<64>(strg.kg) << std::endl << std::endl;
    std::cout << "iw   " <<std::bitset<64>(strg.iw) << std::endl;
    std::cout <<"jw   " << std::bitset<64>(strg.jw) << std::endl;
    std::cout << "kw   " <<std::bitset<64>(strg.kw) << std::endl;

int pause ;
std::cin >> pause;
*/

    return posUint;
}

//************************************************************************
//************************************************************************

// This function computes the inverse of the error function `erf` in the C math
// library. The implementation is based on the rational approximation of Normal
// quantile function available from http://www.jstor.org/stable/2347330
//
// MIT License
//
// Copyright (c) 2017 Lakshay Garg <lakshayg@outlook.in>
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

template <typename T>
T erfinv(T x) {

  if (x < -1 || x > 1) {
    return std::numeric_limits<T>::quiet_NaN();
  } else if (x == 1.0) {
    return std::numeric_limits<T>::infinity();
  } else if (x == -1.0) {
    return -std::numeric_limits<T>::infinity();
  }

  const T LN2 = 6.931471805599453094172321214581e-1;

  const T A0 = 1.1975323115670912564578e0;
  const T A1 = 4.7072688112383978012285e1;
  const T A2 = 6.9706266534389598238465e2;
  const T A3 = 4.8548868893843886794648e3;
  const T A4 = 1.6235862515167575384252e4;
  const T A5 = 2.3782041382114385731252e4;
  const T A6 = 1.1819493347062294404278e4;
  const T A7 = 8.8709406962545514830200e2;

  const T B0 = 1.0000000000000000000e0;
  const T B1 = 4.2313330701600911252e1;
  const T B2 = 6.8718700749205790830e2;
  const T B3 = 5.3941960214247511077e3;
  const T B4 = 2.1213794301586595867e4;
  const T B5 = 3.9307895800092710610e4;
  const T B6 = 2.8729085735721942674e4;
  const T B7 = 5.2264952788528545610e3;

  const T C0 = 1.42343711074968357734e0;
  const T C1 = 4.63033784615654529590e0;
  const T C2 = 5.76949722146069140550e0;
  const T C3 = 3.64784832476320460504e0;
  const T C4 = 1.27045825245236838258e0;
  const T C5 = 2.41780725177450611770e-1;
  const T C6 = 2.27238449892691845833e-2;
  const T C7 = 7.74545014278341407640e-4;

  const T D0 = 1.4142135623730950488016887e0;
  const T D1 = 2.9036514445419946173133295e0;
  const T D2 = 2.3707661626024532365971225e0;
  const T D3 = 9.7547832001787427186894837e-1;
  const T D4 = 2.0945065210512749128288442e-1;
  const T D5 = 2.1494160384252876777097297e-2;
  const T D6 = 7.7441459065157709165577218e-4;
  const T D7 = 1.4859850019840355905497876e-9;

  const T E0 = 6.65790464350110377720e0;
  const T E1 = 5.46378491116411436990e0;
  const T E2 = 1.78482653991729133580e0;
  const T E3 = 2.96560571828504891230e-1;
  const T E4 = 2.65321895265761230930e-2;
  const T E5 = 1.24266094738807843860e-3;
  const T E6 = 2.71155556874348757815e-5;
  const T E7 = 2.01033439929228813265e-7;

  const T F0 = 1.414213562373095048801689e0;
  const T F1 = 8.482908416595164588112026e-1;
  const T F2 = 1.936480946950659106176712e-1;
  const T F3 = 2.103693768272068968719679e-2;
  const T F4 = 1.112800997078859844711555e-3;
  const T F5 = 2.611088405080593625138020e-5;
  const T F6 = 2.010321207683943062279931e-7;
  const T F7 = 2.891024605872965461538222e-15;

  T abs_x = abs(x);

  if (abs_x <= 0.85) {
    T r =  0.180625 - 0.25 * x * x;
    T num = (((((((A7 * r + A6) * r + A5) * r + A4) * r + A3) * r + A2) * r + A1) * r + A0);
    T den = (((((((B7 * r + B6) * r + B5) * r + B4) * r + B3) * r + B2) * r + B1) * r + B0);
    return x * num / den; 
  }

  T r = sqrt(LN2 - log(1.0 - abs_x));

  T num, den;
  if (r <= 5.0) {
    r = r - 1.6;
    num = (((((((C7 * r + C6) * r + C5) * r + C4) * r + C3) * r + C2) * r + C1) * r + C0);
    den = (((((((D7 * r + D6) * r + D5) * r + D4) * r + D3) * r + D2) * r + D1) * r + D0);
  } else {
    r = r - 5.0;
    num = (((((((E7 * r + E6) * r + E5) * r + E4) * r + E3) * r + E2) * r + E1) * r + E0);
    den = (((((((F7 * r + F6) * r + F5) * r + F4) * r + F3) * r + F2) * r + F1) * r + F0);
  }

  if (x < 0) {
    return -num / den;
  } else {
    return num / den;
  }
}

//************************************************************************
//************************************************************************
// FUNCTION return maxwell distrubtion random double number 
// Only select 0.25-0.75
inline double MaxwellDisV( double iKT, double bulkV, double mi0_in)
{
    double sigma = sqrt(iKT / mi0_in);
    double temp = 0.0;
    while ( temp < 0.25 || temp > 0.75)
    {
        temp = dRand();
    }
    return erfinv( 2.0* temp - 1.0) * sqrt(2.0) * sigma + bulkV;
}


//************************************************************************
//************************************************************************
// FUNCTION return maxwell distrubtion random double number of particle 
// energy for magnetic moment/ adiabatic invarient 0.1~10 eV
// the range for random number is -1 ~ 1 and then negetive it to make the 
// major particles is between 0.1~1 eV. should make sigma_in < 0.2
// Only select 0.25-0.75
inline double MaxwellDisEnergy( GridsPoints***** ptrArray_in, uint_64 posUint)
{
    double temp = 0.0;
    double sigma_in = sigma_MaxwellDis;
    double mu_in = mu_MaxwellDis;
    while( temp < 0.15 || temp > 0.60)
    {
        temp = dRand();
    }
    double x_temp = pow( 10.0, erfinv( 2.0* temp -1.0 ) * sqrt(2.0) * sigma_in + mu_in);


//************************************************************************
    struct structg strg_in = {0,0,0,0,0,0,0, 0.0, 0.0, 0.0, 0.0};
    strg_in.face = posUint >> 61;
    for( int i = 0; i < fieldsGridsLevel; i++) 
    {
        strg_in.ig = (strg_in.ig << 1) + ((posUint >> 60   - i*3) & 1);
        strg_in.jg = (strg_in.jg << 1) + ((posUint >> 60-1 - i*3) & 1);
        strg_in.kg = (strg_in.kg << 1) + ((posUint >> 60-2 - i*3) & 1);
    }

    for( int i = fieldsGridsLevel; i < particlesGridsLevel; i++)
    {
        strg_in.iw = (strg_in.iw << 1) + ((posUint >> 60   - i*3) & 1);
        strg_in.jw = (strg_in.jw << 1) + ((posUint >> 60-1 - i*3) & 1);
        strg_in.kw = (strg_in.kw << 1) + ((posUint >> 60-2 - i*3) & 1);
    }
    strg_in.vx = 0.0; strg_in.vy = 0.0; strg_in.vz = 0.0; // not need
    strg_in.mass = 0.0; // not need
//************************************************************************

    Vector3 tempb1=ptrArray_in[strg_in.face][strg_in.ig+1][strg_in.jg+1][strg_in.kg]->B3();
    Vector3 tempb2=ptrArray_in[strg_in.face][strg_in.ig+2][strg_in.jg+1][strg_in.kg]->B3();
    Vector3 tempb3=ptrArray_in[strg_in.face][strg_in.ig+2][strg_in.jg+2][strg_in.kg]->B3();
    Vector3 tempb4=ptrArray_in[strg_in.face][strg_in.ig+1][strg_in.jg+2][strg_in.kg]->B3();
    Vector3 tempb5=ptrArray_in[strg_in.face][strg_in.ig+1][strg_in.jg+1][strg_in.kg+1]->B3();
    Vector3 tempb6=ptrArray_in[strg_in.face][strg_in.ig+2][strg_in.jg+1][strg_in.kg+1]->B3();
    Vector3 tempb7=ptrArray_in[strg_in.face][strg_in.ig+2][strg_in.jg+2][strg_in.kg+1]->B3();
    Vector3 tempb8=ptrArray_in[strg_in.face][strg_in.ig+1][strg_in.jg+2][strg_in.kg+1]->B3();

    double w1 = 1- (strg_in.iw +1) * (strg_in.jw +1) * (strg_in.kw +1) / cellSize3;
    double w2 = 1- (cellSize1- strg_in.iw)* (strg_in.jw +1) * (strg_in.kw +1) / cellSize3;
    double w3 = 1- (cellSize1- strg_in.iw) * (cellSize1- strg_in.jw) * (strg_in.kw +1) / cellSize3;
    double w4 = 1- (strg_in.iw +1) * (cellSize1- strg_in.jw )* (strg_in.kw +1)/ cellSize3;
    double w5 = 1- (strg_in.iw +1) * (strg_in.jw +1) * (cellSize1- strg_in.kw) / cellSize3;
    double w6 = 1- (cellSize1- strg_in.iw)* (strg_in.jw +1) * (cellSize1- strg_in.kw) / cellSize3;
    double w7 = 1- (cellSize1- strg_in.iw) * (cellSize1- strg_in.jw) * (cellSize1- strg_in.kw) / cellSize3;
    double w8 = 1- (strg_in.iw +1) * (cellSize1- strg_in.jw )* (cellSize1- strg_in.kw) / cellSize3;

    Vector3 tempb;
    tempb.Setx(tempb1.x()*w1 + tempb2.x()*w2 + tempb3.x()*w3 + tempb4.x()*w4 
                + tempb5.x()*w5 + tempb6.x()*w6 + tempb7.x()*w7 + tempb8.x()*w8);
                
    tempb.Sety(tempb1.y()*w1 + tempb2.y()*w2 + tempb3.y()*w3 + tempb4.y()*w4 
                + tempb5.y()*w5 + tempb6.y()*w6 + tempb7.y()*w7 + tempb8.y()*w8);
                
    tempb.Setz(tempb1.z()*w1 + tempb2.z()*w2 + tempb3.z()*w3 + tempb4.z()*w4 
                + tempb5.z()*w5 + tempb6.z()*w6 + tempb7.z()*w7 + tempb8.z()*w8);

    return x_temp / tempb.norm();
}


//************************************************************************
//************************************************************************
// FUNCTION return uniform distrubtion random double number
inline double UniformDisX( double x1, double x2)
{
    if ( x1 < x2)
    return dRand() * ( x2 - x1) + x1;
    else 
    return dRand() * ( x1 - x2) + x2; 
}

//************************************************************************
//************************************************************************
// FUNCTION return uniform distrubtion random Vector3
// return a Vector3 starting at the end of v1 and pointing to v2
inline Vector3 UniformDisVector3( Vector3 v1, Vector3 v2)
{   
    Vector3 temp = v2.MinusProduct( v1).ScaleProduct( dRand());
    return temp;
}

//************************************************************************
//************************************************************************
// FUNCTION
// Update info in the grids due to the info of particles and related 
// weighting
void UpdateInfoGrids( GridsPoints***** ptrArray_in, 
                      vector<Particles>* ptrParticlesList_H_in, 
                      vector<Particles>* ptrParticlesList_He_in,
                      vector<Particles>* ptrParticlesList_O_in,
                      vector<Particles>* ptrParticlesListTemp_H_in,
                      vector<Particles>* ptrParticlesListTemp_He_in,
                      vector<Particles>* ptrParticlesListTemp_O_in,
                      double*** ptrVolumeGridArray_in,
                      int timeline_in, int updateInfoPeriod_in);

#endif