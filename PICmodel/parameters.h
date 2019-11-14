#ifndef _PARAMETERS_H_
#define _PARAMETERS_H_
#include<iostream>
#include<cmath>
    using uint_64 = unsigned long long;
    // Parameters constants of the earth
    // Length of face in 6
    const double length = 9220000.0;
    // radius of earth
    const double radius = 6371000.0;
    // diople moment for earth, T*m^3
    const double dMoment = -8e15;

    // Simulation parameters
    const double tstep = 0.1 ;
    // alpha
    const double alpha = 1.0;
    // Ni
    const double Ni = 1.0;
    // charge (C)
    const double qi0 = 1.602e-19;
    // mass (kg)
    const double mi0 = 1.66e-27;
    // electron kT
    const double ekT = 1.0;
    // ion kT k =  1.38 x 10−23 J·K−1, T = 1000 K
    const double ikT = 1.38e-20;
    // initial particle numbers per cell ( count)
    const int iniParticleNumberPerCell = 1;
    // g (m / s2)
    const double gravity = 9.8;
    // number density at base level ( / m^3)
    const double N0_i = 100000000000.0;
    // angular velocity of Earth ( rad/s)
    const double omega_earth = 7.292e-5;

    // These two levels are between 1 and 20, and particlesgridslevel is greater than the fieldgridslevel.
    // The radialgridslevel is calculated from the fieldgridslevel to make the fieldgrids-cell similar to 
    // a cubic. Why do we want a cubic similar cell? The reason is that, the first, we want to use Z-curve 
    // to represent the location, therefore we need the same number grids on each three coordinates; the second,
    // it is easy to calculate the weighting of each particle on the gridspoints, otherwise, with a cuboid cell, 
    // it may be hard to judge which particles in the cell will have weighting on the grids points.
    const int fieldsGridsLevel = 2;
    const int particlesGridsLevel = 3;
    const int cellSize1 = 1 << (particlesGridsLevel - fieldsGridsLevel);
    const int cellSize3 = 1 << (particlesGridsLevel - fieldsGridsLevel + 3);
    const int fieldsGridsSize = 1 << fieldsGridsLevel;
    const int particlesGridsSize = 1 << particlesGridsLevel;
  //  extern uint_64 fieldsGridsSize;
   // extern uint_64 particlesGridsSize;
    //extern uint_64 radialGridsSize;

    // Determin the bottom L: LMin; Calculate the top L: LMax
    const int totalFace = 6;
    const double AltitudeMin = 600000.0; // unit: m 
    const double LMin = 1 + AltitudeMin / radius;

    // Calculate the array in the radial direction
    const double ratio = 1.0+(length/fieldsGridsSize)/radius;
    const double logRatio = log10(1.0+(length/fieldsGridsSize)/radius);
    // L(m) = radio ** m, L(0) = 1
    // Construct an array from L(m) to L(n), the number is
 //   const int numMin = round(log10(LMin)/logRatio);
//    const int numMax = numMin + radialGridsSize;
    const double LMax = LMin * pow(ratio, fieldsGridsSize);

    
    static int h5FileCheck = 0; // 0: create a new h5, 1: open exist h5


#endif