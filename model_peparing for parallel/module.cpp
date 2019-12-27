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
    Vector3*** ptrVectorCellArray = VectorCellField();  

    // Prerun 1.3 // Create grids field array of volum for one face of six
    cout << " Create array of Volume at cells and at grids" << endl;
    double*** ptrVolumeCellArray = VolumeCellsField( ptrArray);
    
    cout << " Create array of Volume at cells and at grids" << endl;
    double*** ptrVolumeGridArray = VolumeGridsField( ptrVolumeCellArray);
    
    // Initialize condition
    SetInitialCondition( ptrArray, ptrVectorCellArray, ptrVolumeCellArray);

    // Prerun 1.4 // Create particles list, initialize the velocity and position of each particles
    cout << " Create particles list of main domain" << endl;
    
    vector<Particles>* ptrParticlesList_H;
    vector<Particles>* ptrParticlesList_He;
    vector<Particles>* ptrParticlesList_O;
    vector<int>* ptrParticlesList_He_out = new vector<int>;
    vector<int>* ptrParticlesList_H_out = new vector<int>; 
    vector<int>* ptrParticlesList_O_out = new vector<int>;
#pragma omp parallel
{
    #pragma omp sections
    {
        #pragma omp section
        {
        ptrParticlesList_H = ParticlesLists( ptrArray, ptrVolumeCellArray, mi0_H, N0_H);
        }
        #pragma omp section
        {
        ptrParticlesList_He = ParticlesLists( ptrArray, ptrVolumeCellArray, mi0_He, N0_He);    
        }
        #pragma omp section
        {
        ptrParticlesList_O = ParticlesLists( ptrArray, ptrVolumeCellArray, mi0_O, N0_O);    
        }
    }
    #pragma omp barrier
}
   
    
    // Run 2.0
    
    cout << " Start" << endl;
    for( int timeline = 1; timeline <= timeLineLimit; timeline++)   // timeline start with 1
    {
    if( timeline==1){
    // Printout the initial condition 
    PrintOutHdf5( ptrArray, timeline, h5FileCheck);
    }

    std::cout << "timeline" << timeline << std::endl;

#pragma omp parallel
{
    #pragma omp sections
    {
        #pragma omp section
        {

        // Run 2.1 // Particles in main domain
        for( auto iteratorM = ptrParticlesList_H->begin(); iteratorM != ptrParticlesList_H->end(); ++iteratorM)
        {
            Particles temp = *iteratorM;
            struct structg tempStr = temp.InttoStrp1();
            
            int check; // check whether in the main domain or not, "0" means in "1" means out
            // update velocity // update position
            check = temp.BorisMethod( &tempStr, ptrArray, mi0_H);
            // check if still in the main domain
            if( check == 1) // out of the domain
            {
                iteratorM = ptrParticlesList_H->erase( iteratorM);
                int tempint = iteratorM - ptrParticlesList_H->begin();
                ptrParticlesList_H_out->push_back( tempint);
                
            }
        }
        
        cout << "Particles H " << ptrParticlesList_H->size() << endl;
        }

        #pragma omp section
        {
        for( auto iteratorM = ptrParticlesList_He->begin(); iteratorM != ptrParticlesList_He->end(); ++iteratorM)
        {
            Particles temp = *iteratorM;
            struct structg tempStr = temp.InttoStrp1();
            
            int check; // check whether in the main domain or not, "0" means in "1" means out
            //test function//  double xxx= temp.VelParticles().x();    //cout << xxx << " ";
            // update velocity // update position
            check = temp.BorisMethod( &tempStr, ptrArray, mi0_He);
            // check if still in the main domain
            if( check == 1) // out of the domain
            {
                iteratorM = ptrParticlesList_He->erase( iteratorM);
                int tempint = iteratorM - ptrParticlesList_He->begin();
                ptrParticlesList_He_out->push_back( tempint);
            }
        }
        cout << "Particles He " << ptrParticlesList_He->size() << endl;
        }
        #pragma omp section
        {
        for( auto iteratorM = ptrParticlesList_O->begin(); iteratorM != ptrParticlesList_O->end(); ++iteratorM)
        {
            Particles temp = *iteratorM;
            struct structg tempStr = temp.InttoStrp1();
            
            int check; // check whether in the main domain or not, "0" means in "1" means out
            // update velocity // update position
            check = temp.BorisMethod( &tempStr, ptrArray, mi0_O);
            // check if still in the main domain
            if( check == 1) // out of the domain
            {
                iteratorM = ptrParticlesList_O->erase( iteratorM);
                int tempint = iteratorM - ptrParticlesList_O->begin();
                ptrParticlesList_O_out->push_back( tempint);
            }
        }
        cout << "Particles O " << ptrParticlesList_O->size() << endl;
        }     
    }   
    #pragma omp barrier
}

        // Run 2.2 // Create temp particle lists    
        vector<Particles>* ptrParticlesListTemp_H; // = ParticlesListsTemp( ptrArray, ptrVolumeCellArray, mi0_H , 1);
        vector<Particles>* ptrParticlesListTemp_He; // = ParticlesListsTemp( ptrArray, ptrVolumeCellArray, mi0_He, 4);
        vector<Particles>* ptrParticlesListTemp_O; // = ParticlesListsTemp( ptrArray, ptrVolumeCellArray, mi0_O, 16);

    #pragma omp parallel
    {
    #pragma omp sections
    {
    #pragma omp section
    {   

        ptrParticlesListTemp_H = ParticlesListsTemp( ptrArray, ptrVolumeCellArray, mi0_H , 1);     
        // Run 2.3 // Particles in temp domain
        for( auto iterator = ptrParticlesListTemp_H->begin(); iterator != ptrParticlesListTemp_H->end(); ++iterator)
        {
            Particles temp = *iterator;
            struct structg tempStr = temp.InttoStrp1();

       

            int check; // check whether in the main domain or not, "0" means in "1" means out
            // update velocity  // update position
            check = temp.BorisMethod( &tempStr, ptrArray, mi0_H);
   

            // check if still in the main domain
            if( check == 0) // in the domain
            {    
        
/*        std::cout << std::endl << std::bitset<64>(temp.PosUint()) << " info " <<
        tempStr.face << " " << tempStr.ig << " " << tempStr.jg << " " << tempStr.kg <<" vel " << tempStr.vx << " " << tempStr.vy << " " << tempStr.vz << std::endl;

        tempStr = temp.InttoStrp1();
        
        std::cout << std::bitset<64>(temp.PosUint()) << " info " << 
        tempStr.face << " " << tempStr.ig << " " << tempStr.jg << " " << tempStr.kg << " vel " << tempStr.vx << " " << tempStr.vy << " " << tempStr.vz << std::endl;
*/
                if( ptrParticlesList_H_out->size() > 0)
                {
                    auto temp_pos_out = ptrParticlesList_H_out->end() - 1;
                    (*ptrParticlesList_H)[ *temp_pos_out] = temp;
                    temp_pos_out = ptrParticlesList_H_out->erase(temp_pos_out);
                }
                else
                {
                    ptrParticlesList_H->push_back( temp);                   
                }
                       //        iterator = ptrParticlesListTemp_H->erase( iterator);
            }
        }
    }
    #pragma omp section
    {  
        ptrParticlesListTemp_He = ParticlesListsTemp( ptrArray, ptrVolumeCellArray, mi0_He, 4);
        for( auto iterator = ptrParticlesListTemp_He->begin(); iterator != ptrParticlesListTemp_He->end(); ++iterator)
        {
            Particles temp = *iterator;
            struct structg tempStr = temp.InttoStrp1();
            int check; // check whether in the main domain or not, "0" means in "1" means out
            // update velocity  // update position
            check = temp.BorisMethod( &tempStr, ptrArray, mi0_He);
            // check if still in the main domain
            if( check == 0) // in the domain
            {   
                if( ptrParticlesList_He_out->size() > 0)
                {
                    auto temp_pos_out = ptrParticlesList_He_out->end() - 1;
                    (*ptrParticlesList_He)[ *temp_pos_out] = temp;
                    temp_pos_out = ptrParticlesList_He_out->erase(temp_pos_out);
                }
                else
                {
                    ptrParticlesList_He->push_back( temp);
                }
                 //    iterator = ptrParticlesListTemp_He->erase( iterator);
            }
        }
    }
    #pragma omp section
    {    
        ptrParticlesListTemp_O = ParticlesListsTemp( ptrArray, ptrVolumeCellArray, mi0_O, 16);
        for( auto iterator = ptrParticlesListTemp_O->begin(); iterator != ptrParticlesListTemp_O->end(); ++iterator)
        {
            Particles temp = *iterator;
            struct structg tempStr = temp.InttoStrp1();
            int check; // check whether in the main domain or not, "0" means in "1" means out
            // update velocity  // update position
            check = temp.BorisMethod( &tempStr, ptrArray, mi0_O);
            // check if still in the main domain
            if( check == 0) // in the domain
            {   
                if( ptrParticlesList_O_out->size() > 0)
                {
                    auto temp_pos_out = ptrParticlesList_O_out->end() - 1;
                    (*ptrParticlesList_O)[ *temp_pos_out] = temp;
                    temp_pos_out = ptrParticlesList_O_out->erase(temp_pos_out);
                }
                else
                {
                    ptrParticlesList_O->push_back( temp);
                }
            }
        }
        
    }
    }
    #pragma omp barrier

    }

    // Run 2.5 // Update info in grids 
    // Run 2.5.1 // Accumulate density and velocity per timestep and average them to get the
    // value of density and velocity per updatetimeperiod
   UpdateInfoGrids( ptrArray, 
                     ptrParticlesList_H,
                     ptrParticlesList_He,
                     ptrParticlesList_O, 
                     ptrParticlesListTemp_H, 
                     ptrParticlesListTemp_He,
                     ptrParticlesListTemp_O,
                     ptrVolumeGridArray, 
                     timeline, 
                     updateInfoPeriod);
    

    // Run 2.4 // delete temp particlesLists
    delete ptrParticlesListTemp_H;
    delete ptrParticlesListTemp_He;
    delete ptrParticlesListTemp_O;



    if( timeline % updateInfoPeriod ==0)
    {
    std::cout << " Update gridspoints info " << timeline << std::endl;
        // Run 2.5.2 
        for( int face = 0; face < 6; face++)
        {
            if( update_type == 0)
            {
            // Update E without current
            // Calculate curl dB update ve3, ve3 = v3
            ptrVectorCellArray = ValueCurlField( ptrVectorCellArray, ptrVolumeCellArray, ptrArray, face, 'D');
            UpdateVe3( ptrVectorCellArray, ptrArray, face);
            // Calculate the gradient of Pe
            ptrVectorCellArray = ValueGradient( ptrVectorCellArray, ptrVolumeCellArray, ptrArray, face, 'P');
            UpdateE3( ptrVectorCellArray, ptrArray, face); // update E
            } else
            {
            // Update grids info
            // Calculate curl dB update ve3
            ptrVectorCellArray = ValueCurlField( ptrVectorCellArray, ptrVolumeCellArray, ptrArray, face, 'D');
            UpdateVe3( ptrVectorCellArray, ptrArray, face);
            // Calculate the gradient of Pe update E
            ptrVectorCellArray = ValueGradient( ptrVectorCellArray, ptrVolumeCellArray, ptrArray, face, 'P');
            UpdateE3( ptrVectorCellArray, ptrArray, face);
            // Calculate the curl E update B
            ptrVectorCellArray = ValueCurlField( ptrVectorCellArray, ptrVolumeCellArray, ptrArray, face, 'E');
            UpdateB3( ptrVectorCellArray, ptrArray, face);
            // Update gradient norm B
            ptrVectorCellArray = ValueGradient( ptrVectorCellArray, ptrVolumeCellArray, ptrArray, face, 'B');
            UpdateGradBNorm( ptrVectorCellArray, ptrArray, face);
            }
        }
    }

    // Postrun 3.0 // Printout info in grids points
    if( timeline % printTimePeriod ==0)
    {
        std::cout << " PrintOut " << std::endl;
        PrintOutHdf5( ptrArray, timeline, h5FileCheck);
    }
    }
    delete ptrVectorCellArray;
    delete ptrArray;
    delete ptrParticlesList_H;
    delete ptrParticlesList_He;
    delete ptrParticlesList_O;

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

