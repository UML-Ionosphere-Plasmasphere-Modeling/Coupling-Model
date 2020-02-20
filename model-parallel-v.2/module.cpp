#include<iostream>
#include <list>
#include <vector>
#include <memory>
#include <string>
#include "parameters.h"
#include "fieldsgrids.h"
#include "particles.h"
#include "vector3.h"
#include "module.h"
#include "module_0.h"
#include "module_1.h"
#include <cmath>
#include "H5Cpp.h"
#include <bitset>

using std::cout;
using std::endl;
using std::vector;



//************************************************************************
//************************************************************************
// Test FUNCTION // Put at end of module.cpp
// This test fuction only apply for fixed global E and B, particle moving would not 
// affect the E and B at gridspoints.
// This function only apply for testing the moving particles.
// Assume we have a list of Particles
//************************************************************************
//************************************************************************
void ProcessFunc()
{
    // Prerun 1.0 // Create Grids, including B and Pos. And then Velocity (corotation) and N (exponential )
    // Prerun 1.1 // And then E (electron momentum equation).
    cout << " Create Grids" << endl;
    GridsPoints***** ptrArray =  GridsCreation();

    Titheridge_Te( ptrArray); // initial Temprature of electron
    
//    SetTopBoundary( ptrArray);
//    SetBotBoundary( ptrArray);
    
    cout << LMin << " " << LMax << endl;
    // Prerun 1.2 // Create Cell centered field array for nesseary calculation for one face of six
    // The size is [fsize+2][fsize+2][fsize+2]
    Vector3*** ptrVectorCellArray = VectorCellField();  
    Vector3*** ptrVelVectorCellArray = VectorCellField();
    Vector3*** ptrGradVectorCellArray= VectorCellField();

    // Prerun 1.3 // Create grids field array of volum for one face of six
    // The size is [fsize+2][fsize+2][fsize+2]
    cout << " Create array of Volume at cells and at grids" << endl;
    double*** ptrVolumeCellArray = VolumeCellsField( ptrArray);
    
    cout << " Create array of Volume at cells and at grids" << endl;
    // The size if [fsize+1][fsize+1][fsize+1]
    double*** ptrVolumeGridArray = VolumeGridsField( ptrVolumeCellArray);

    // Presun 1.4 // Create Cell centered field array for E
    // [totalface * fsize+2 * fsize+2 * fsize+2] 
    Vector3***** ptrEVectorCellArray = EVectorCellArray( ptrArray);
    // Presun 1.5 // Create Face centered field array for B
    // [direction * face * (fsize+1) * (fsize+1) * (fsize+1)]
    Vector3***** ptrBVectorFaceArray = BVectorFaceArray( ptrArray);
    
    // Initialize condition
    if( update_type == 0)
    {
        SetInitialCondition( ptrArray, ptrVectorCellArray, ptrVolumeCellArray);
    }
    // Prerun 1.4 // Create particles list, initialize the velocity and position of each particles
    cout << " Create particles list of main domain" << endl;
    
    vector<Particles> ptrParticlesList_H;
    vector<Particles> ptrParticlesList_He;
    vector<Particles> ptrParticlesList_O;
    vector<int> ptrParticlesList_out_H;
    vector<int> ptrParticlesList_out_He; 
    vector<int> ptrParticlesList_out_O;
    ptrParticlesList_H.reserve(500000000);
    ptrParticlesList_He.reserve(500000000);   
    ptrParticlesList_O.reserve(500000000);
    ptrParticlesList_out_H.reserve(50000000);
    ptrParticlesList_out_He.reserve(50000000);
    ptrParticlesList_out_O.reserve(50000000);

    vector<Particles> ptrParticlesListTemp_H; 
    vector<Particles> ptrParticlesListTemp_He;
    vector<Particles> ptrParticlesListTemp_O; 
    ptrParticlesListTemp_H.reserve(5000000);
    ptrParticlesListTemp_He.reserve(5000000);
    ptrParticlesListTemp_O.reserve(5000000);

    #pragma omp parallel
    {
        #pragma omp sections
        {
            #pragma omp section
            {
                ParticlesLists( ptrParticlesList_H, 
                                ptrArray, 
                                ptrVolumeCellArray, 
                                mi0_H, 
                                N0_H);
            }
            #pragma omp section
            {
                ParticlesLists( ptrParticlesList_He, 
                                ptrArray, 
                                ptrVolumeCellArray, 
                                mi0_He, 
                                N0_He);
            }
            #pragma omp section
            {
                ParticlesLists( ptrParticlesList_O, 
                                ptrArray, 
                                ptrVolumeCellArray, 
                                mi0_O, 
                                N0_O);
            }
        }
        #pragma omp barrier
    }
   
    
    // Run 2.0
    
    cout << " Start" << endl;


    for( int timeline = 1; timeline <= timeLineLimit; timeline++)   // timeline start with 1
    {
        // set boundary condition: 60s initial time interval
        // can be developed for info exchanging boundary
        if( update_type == 1)
        { 
            if( timeline >= botBoundaryInitialTimeStart/tstep && 
                timeline <= (botBoundaryInitialTimeStart+botBoundaryInitialTime)/tstep &&
                initial_bot_type ==0)
            {
                SetRotationalVelBotBoundary( ptrArray, timeline);
            }
            if( timeline >= topBoundaryInitialTimeStart/tstep && 
                timeline <= (topBoundaryInitialTimeStart+topBoundaryInitialTime)/tstep && 
                initial_top_type ==1)
            {
                SetConvectionVelTopBoundary( ptrArray, timeline);
            }
        }
        
        // const info print out
        if( timeline == 1)
        {
            std::cout << " PrintOut Const" << std::endl;
            PrintOutHdf5( ptrArray, timeline, h5FileCheck);
        }
        // average pho, v, update grids info B, E & reset pho, v
        if( timeline % updateInfoPeriod ==0)
        {
            cout << timeline << " H " << ptrParticlesList_H.size() << " " << ptrParticlesList_out_H.size();
            cout << " He " << ptrParticlesList_He.size() << " " << ptrParticlesList_out_He.size();
            cout << " O " << ptrParticlesList_O.size() << " " << ptrParticlesList_out_O.size() << endl;
        
            // average pho and v
            CalculatingAveragedPhoVatGrids( ptrArray, 
                                            ptrVolumeGridArray,
                                            updateInfoPeriod);
            // Run 2.5.2 
            for( int face = 0; face < 6; face++)
            {
                if( update_type == 0)
                {
                // Update E without current
                // Calculate curl dB update ve3, ve3 = v3
                ptrVectorCellArray = ValueCurlField(ptrVectorCellArray, 
                                                    ptrVolumeCellArray, 
                                                    ptrArray, 
                                                    face, 
                                                    'D');
                UpdateVe3(  ptrVectorCellArray, 
                            ptrArray,   
                            face);
                // Calculate the gradient of Pe
                ptrVectorCellArray = ValueGradient( ptrVectorCellArray, 
                                                    ptrVolumeCellArray, 
                                                    ptrArray, 
                                                    face, 
                                                    'P');
                UpdateE3(   ptrVectorCellArray, 
                            ptrArray, 
                            face); // update E
                } else
                {
                // 1. Calculate the curl B, need ptrBFaceArray 
                // With the B on the faces and area vectors
                ptrVectorCellArray = CurlBCellArray(ptrArray, 
                                                    ptrVectorCellArray,
                                                    ptrBVectorFaceArray,
                                                    ptrVolumeCellArray,
                                                    face);
                // 2. Calculate the gradient of Pe
                // ( fsize+2 * fsize+2 * fsize)
                ptrGradVectorCellArray = ValueGradient( ptrVectorCellArray, 
                                                    ptrVolumeCellArray, 
                                                    ptrArray, 
                                                    face, 
                                                    'P');
                // 3. Calculate the B at the center of cells
                // 4. Update E at the center of cells and at the grids
                UpdateECellArray(  ptrArray, 
                                   ptrEVectorCellArray,
                                   ptrVectorCellArray,
                                   ptrGradVectorCellArray,
                                   face);
                // 5. Update B at center of cells and at the grids
                BVectorFaceArrayUpdate( ptrArray, ptrBVectorFaceArray);
                BVectorGridsArrayUpdate( ptrArray, ptrBVectorFaceArray);
                

                // Update gradient norm B
                ptrVectorCellArray = ValueGradient( ptrVectorCellArray, 
                                                    ptrVolumeCellArray, 
                                                    ptrArray, 
                                                    face, 
                                                    'B');
                UpdateGradBNorm(ptrVectorCellArray, 
                                ptrArray, 
                                face);
                }
            }
            // printout 
            if( timeline % printTimePeriod == 0)
            {
                std::cout << " PrintOut  " << timeline << std::endl;
                PrintOutHdf5( ptrArray, timeline, h5FileCheck);
            }
            // reset pho and v
            ResetPhoVatGrids( ptrArray);
        }

        // iterate particles
        IterateParticlesMain(   ptrArray, 
                                ptrParticlesList_H, 
                                ptrParticlesList_out_H,
                                mi0_H);
        IterateParticlesMain(   ptrArray, 
                                ptrParticlesList_He, 
                                ptrParticlesList_out_He,
                                mi0_He);
        IterateParticlesMain(   ptrArray, 
                                ptrParticlesList_O, 
                                ptrParticlesList_out_O,
                                mi0_O);
                                
        #pragma omp parallel
        {
            #pragma omp sections
            {
                #pragma omp section
                {
                    ParticlesListsTemp( ptrParticlesListTemp_H, 
                                        ptrArray, 
                                        ptrVolumeCellArray, 
                                        mi0_H, 
                                        1);
                }
                #pragma omp section
                {
                    ParticlesListsTemp( ptrParticlesListTemp_He, 
                                        ptrArray, 
                                        ptrVolumeCellArray, 
                                        mi0_He, 
                                        4);
                }
                #pragma omp section
                {
                    ParticlesListsTemp( ptrParticlesListTemp_O, 
                                        ptrArray, 
                                        ptrVolumeCellArray, 
                                        mi0_O, 
                                        16);
                }
            }
            #pragma omp barrier
        }
            
        // temp particles
        IterateParticlesTemp(   ptrArray, 
                                ptrParticlesList_H, 
                                ptrParticlesListTemp_H,
                                ptrParticlesList_out_H,
                                mi0_H);
        IterateParticlesTemp(   ptrArray, 
                                ptrParticlesList_He, 
                                ptrParticlesListTemp_He,
                                ptrParticlesList_out_He,
                                mi0_He);
        IterateParticlesTemp(   ptrArray, 
                                ptrParticlesList_O, 
                                ptrParticlesListTemp_O,
                                ptrParticlesList_out_O,
                                mi0_O);         

        ptrParticlesListTemp_H.clear();
        ptrParticlesListTemp_He.clear();
        ptrParticlesListTemp_O.clear();

        }

    delete ptrVectorCellArray;
    delete ptrArray;
    ptrParticlesList_H.clear();
    ptrParticlesList_He.clear();
    ptrParticlesList_O.clear();
    ptrParticlesList_out_H.clear();
    ptrParticlesList_out_He.clear();
    ptrParticlesList_out_O.clear();
    ptrParticlesList_H.shrink_to_fit();
    ptrParticlesList_He.shrink_to_fit();
    ptrParticlesList_O.shrink_to_fit();
    ptrParticlesList_out_H.shrink_to_fit();
    ptrParticlesList_out_He.shrink_to_fit();
    ptrParticlesList_out_O.shrink_to_fit();

}







//************************************************************************
//************************************************************************
// FUNCTION 
// Process control function, basicly it has two parts and repeat between them.
// 1 the particles parts that is mainly for moving the particles, and 2 the
// grid parts that is mainly for updating E B or other general variables on 
// grid nodes. 
//
// Procedures:
// Assume E and B on gridspoints are know before the first loop
// Particles Part I: Only need to go through the particles list once, for 
// each particle class:
// 1. Get local E and B, details in particles.cpp, update locations 
// and velocity of each particles, and return a structp of 
// (ig, jg, kg, iw, jw, kw, vp) for each particles.
// 2. For each strucp of each particles, calculate related information of 
// density and vi on each gridspoints. Old density and vi can be covered
// directly.
// Gridspoints Part II:
// Follwed Part I, assume density and vi are updated at gridspoints.
// Matrix of curl E and B need to be applied in this part to store curl
// E and B.
// The information update face by face. 
// The calculation of curl E and B would be demonstrated elsewhere.
// 1. Ampere's Law: Calculate curl B at each gridspoints, and then calculate
// ve using Ampere's Law.
// 2. Electron's momentum equation: Calculate gradient of Pe, and then 
// calculate new E at each gridspoints.
// 3. Faraday's Law: Calculate curl E, and then update B.
//************************************************************************
//************************************************************************

