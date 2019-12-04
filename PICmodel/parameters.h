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

    // charge (C)
    const double qi0 = 1.602e-19;
    // mass of ion in unit of (kg) for Hydrogen
    const double mi0_H = 1.66e-27; 
    // mass of ion in unit of kg for Helium
    const double mi0_He = 3.32e-27;
    // mass of ion in unit of kg for Oxygen
    const double mi0_O = 2.656e-26;
    // Vacuum permeability mu0
    const double mu0 = 1.256637e-6;

    // ion kT k =  1.38 x 10−23 J·K−1, T = 1000 K
    const double ikT = 1.38e-20;
    // Boltzmann_k constant in unit J·K−1
    const double boltzmann_k = 1.38e-23;
    // initial particle numbers per cell ( count)
    const int iniParticleNumberPerCell = 100;
    // g (m / s2)
    const double gravity = 9.8;
    // number density at base level ( / m^3) for H
    const double N0_H =  100000000000.0;
    // number density at base level for He
    const double N0_He = 100000000000.0;
    // number density at base level for O
    const double N0_O =  100000000000.0;
    // angular velocity of Earth ( rad/s)
    const double omega_earth = 7.292e-5;

    // These two levels are between 1 and 20, and particlesgridslevel is greater than the fieldgridslevel.
    // The radialgridslevel is calculated from the fieldgridslevel to make the fieldgrids-cell similar to 
    // a cubic. Why do we want a cubic similar cell? The reason is that, the first, we want to use Z-curve 
    // to represent the location, therefore we need the same number grids on each three coordinates; the second,
    // it is easy to calculate the weighting of each particle on the gridspoints, otherwise, with a cuboid cell, 
    // it may be hard to judge which particles in the cell will have weighting on the grids points.
    // fieldGridsLevel max shouble better be 9 
    // particlesGridsLevel max should be 10 greater than fieldGridsLevel
    
    const int fieldsGridsLevel = 5;
    const int particlesGridsLevel = fieldsGridsLevel + 10;
    const int cellSize1 = 1 << (particlesGridsLevel - fieldsGridsLevel);
    const int cellSize3 = 1 << (particlesGridsLevel - fieldsGridsLevel) <<  (particlesGridsLevel - fieldsGridsLevel) <<  (particlesGridsLevel - fieldsGridsLevel);
    const int fieldsGridsSize = 1 << fieldsGridsLevel;
    const int particlesGridsSize = 1 << particlesGridsLevel;
    const int tempGridsCellLevel = 1;
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

    
    const double LMin_maindomain = LMin * pow(ratio, 1.0);
    const double LMax_maindomain = LMin * pow(ratio, fieldsGridsSize-1);
//************************************************************************
//************************************************************************
// run time control
//************************************************************************
//************************************************************************

    // Simulation parameters (unit) s
    const double tstep = 0.01 ;
    static int timeLineLimit = 720000;
    static int printTimePeriod = 60000;
    static int updateInfoPeriod = 10;

//************************************************************************
//************************************************************************
// For printout
//************************************************************************
//************************************************************************
    
    static int h5FileCheck = 0; // 0: create a new h5, 1: open exist h5


//************************************************************************
//************************************************************************
// For some parameters of some functions
//************************************************************************
//************************************************************************           
    const double mu_MaxwellDis = 0.0;
    const double sigma_MaxwellDis = 0.15;

    
//************************************************************************
//************************************************************************
// For some control
//************************************************************************
//************************************************************************ 
    const int update_type = 1; // 0- no current 1- with current
    
//************************************************************************
//************************************************************************
// For top boundary initialization, in degree
// c0_latitude > r0_latitude
//************************************************************************
//************************************************************************ 
    const double r0_latitude = 63.0;
    const double c0_latitude = 75.0;
    const double t0_convection = 1800.0;
//************************************************************************
//************************************************************************
// For bpt boundary initialization, in degree
//************************************************************************
//************************************************************************ 
    const double rho_max = 7.0e-18;
    const double rho_min = 1.0e-18;
    const double ratioH = 1.0 / 21.0;
    const double ratioHe = 4.0 /21.0;
    const double ratioO = 16.0 / 21.0;

#endif