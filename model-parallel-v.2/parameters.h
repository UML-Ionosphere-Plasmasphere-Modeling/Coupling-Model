#ifndef _PARAMETERS_H_
#define _PARAMETERS_H_
#include<iostream>
#include<cmath>
    using uint_64 = unsigned long long;
//************************************************************************
//************************************************************************
// earth physics parameters control
//************************************************************************
//************************************************************************

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
    // g (m / s2)
    const double gravity = 9.8;
    // angular velocity of Earth ( rad/s)
    const double omega_earth = 7.292e-5;

//************************************************************************
//************************************************************************
// grids size control
//************************************************************************
//************************************************************************

    // These two levels are between 1 and 20, and particlesgridslevel is greater than the fieldgridslevel.
    // The radialgridslevel is calculated from the fieldgridslevel to make the fieldgrids-cell similar to 
    // a cubic. Why do we want a cubic similar cell? The reason is that, the first, we want to use Z-curve 
    // to represent the location, therefore we need the same number grids on each three coordinates; the second,
    // it is easy to calculate the weighting of each particle on the gridspoints, otherwise, with a cuboid cell, 
    // it may be hard to judge which particles in the cell will have weighting on the grids points.
    // fieldGridsLevel max shouble better be 9 
    // particlesGridsLevel max should be 10 greater than fieldGridsLevel
    
    const int fieldsGridsLevel = 4;
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
    const double tstep = 3 ;
    static int timeLineLimit = 240000;
    static int printTimePeriod = 10; //60000;
    static int updateInfoPeriod = 5; //10;

//************************************************************************
//************************************************************************
// For printout
//************************************************************************
//************************************************************************
    
    static int h5FileCheck = 0; // 0: create a new h5, 1: open exist h5


//************************************************************************
//************************************************************************
// Initial particles control
//************************************************************************
//************************************************************************           
    const double mu_MaxwellDis = 0.0;
    const double sigma_MaxwellDis = 0.15;
    // number density at base level ( / m^3) for H
    const double N0_H =  100000000000.0;
    // number density at base level for He
    const double N0_He = 100000000000.0;
    // number density at base level for O
    const double N0_O =  100000000000.0;
    // initial particle numbers per cell ( count)
    const int iniParticleNumberPerCell = 10;
    // temp cell particle number per cell ( count)
    const int tempParticleNumberPerCell = 100;

    
//************************************************************************
//************************************************************************
// For some control
//************************************************************************
//************************************************************************ 
    // 0- no current 1- with current
    const int update_type = 1;  

    // work only with update_type = 1
    // 0- zero velocity 
    //    increased inner boundary velocity to normal
    //    normal density for all
    // 1- normal velocity for all ( not applied)
    //    increased outer boundary velocity
    //    normal density for all
    const int initial_bot_type = 0;
    // 0- no convectional top boundary
    // 1- with convectional top boundary
    const int initial_top_type = 1;

    const double botBoundaryInitialTimeStart = 0.0;
    const double topBoundaryInitialTimeStart = 180.0;
    const double botBoundaryInitialTime = 60.0;
    const double topBoundaryInitialTime = 60.0;
    
//************************************************************************
//************************************************************************
// For top boundary initialization, in degree
// c0_latitude > r0_latitude
//************************************************************************
//************************************************************************ 
    const double r0_latitude = 70.0;
    const double c0_latitude = 82.0;
    const double t0_convection = 3600.0;
//************************************************************************
//************************************************************************
// For bot boundary initialization, in degree
//************************************************************************
//************************************************************************ 
    const double rho_max = 7.0e-18;
    const double rho_min = 1.0e-18;
    const double ratioH = 1.0 / 3.0;
    const double ratioHe = 1.0 /3.0;
    const double ratioO = 1.0 / 3.0;
    const double ratioH_bot = 0.05;
    const double ratioHe_bot = 0.1;
    const double ratioO_bot = 0.85;
    const double ratioH_top = 0.85;
    const double ratioHe_top = 0.1;
    const double ratioO_top = 0.05;

#endif