#include<iostream>
#include <list>
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
using std::list;
using std::shared_ptr;
using std::make_shared;

//************************************************************************
//************************************************************************
// FUNCTION
// Initialization the particles for velocity and position 
// Generate lists of particles
//************************************************************************
//************************************************************************
list<Particles>* ParticlesLists( GridsPoints***** ptrArray_in, double*** ptrVolumeCellArray_in, double mi0, double N0)
{
   //list<Particles> listP;
//    auto listsPtr = make_shared<list<Particles>>();
    list<Particles>* listsPtr = new list<Particles>;
//   shared_ptr<Particles> p1 = make_shared<Particles> (); // test
    double scaleHeight = ikT / mi0 / gravity; // assume the T is const for initiallization

    for( int i = 0; i < totalFace; i++)
    {
        for( int j = 1; j <= fieldsGridsSize; j++)
        {
            for( int k = 1; k <= fieldsGridsSize; k++)
            {
                for( int s = 1; s <= fieldsGridsSize-2; s++)
                {
                    // number of real particles ( notice the unit is number density)
                    double N = N0 * exp(-1.0 * (ptrArray_in[i][j][k][s]->Pos3().norm() - radius) / scaleHeight);
                    // mass of each simulation particle
                    double mi_simu = N / iniParticleNumberPerCell *mi0 * ptrVolumeCellArray_in[j][k][s];

                    for ( int t = 1; t <= iniParticleNumberPerCell; t++)
                    {
                    // calculate random position
                    //test 
                    /*
                    Vector3 temp1 = ptrArray_in[i][j][k][s]->Pos3();
                    Vector3 temp2 = ptrArray_in[i][j][k][s+1]->Pos3();
                    Vector3 temp = UniformDisVector3( ptrArray_in[i][j][k][s]->Pos3(),ptrArray_in[i][j][k][s+1]->Pos3());

                    Vector3 vPos = ptrArray_in[i][j][k][s]->Pos3().PlusProduct(UniformDisVector3( ptrArray_in[i][j][k][s]->Pos3(),ptrArray_in[i][j][k][s+1]->Pos3()));
                    vPos = vPos.PlusProduct( UniformDisVector3( ptrArray_in[i][j][k][s]->Pos3(),ptrArray_in[i][j][k+1][s]->Pos3()));
                    vPos = vPos.PlusProduct( UniformDisVector3( ptrArray_in[i][j][k][s]->Pos3(),ptrArray_in[i][j+1][k][s]->Pos3()));
                    uint_64 intPos = vPos.Uint_64_Trans();
                    */
                    uint_64 intPos = UniDisInCell( ptrArray_in, i, j, k, s);
                    // calculate random velocity
         //           std::cout<< MaxwellDisV( ikT, ptrArray_in[i][j][k][s]->Vel3().x()) << std::endl;
         //           std:: cout << ptrArray_in[i][j][k][s]->B3().norm()<< " B "<< std::endl;                        
                    

                    Vector3 vVel = Vector3( MaxwellDisV( ikT, ptrArray_in[i][j][k][s]->Vel3().x(), mi0),
                                            MaxwellDisV( ikT, ptrArray_in[i][j][k][s]->Vel3().y(), mi0),
                                            MaxwellDisV( ikT, ptrArray_in[i][j][k][s]->Vel3().z(), mi0));
                    double mu_simu = MaxwellDisEnergy( ptrArray_in, intPos);
                    // put the particles at the end of list
                    Particles tempP= Particles(intPos, vVel, mi_simu, mu_simu);
                    listsPtr->push_back( tempP);       
                    }                
                }
            }
        }
    }    
    cout << "Particles initial size" << listsPtr->size() << endl;
//    listsP.emplace_back(radius*3, radius*5, radius*1, 0.0, 0.0, 0.0);

//    for( auto ptrL1 = listsPtr->begin(); ptrL1 != listsPtr->end(); ptrL1++)
//   {  ptrL1->PrintP();  }
   return listsPtr; 
}



//************************************************************************
//************************************************************************
// FUNCTION
// Initialization the particles for velocity and position 
// Generate lists of particles for bot and top region temp
//************************************************************************
//************************************************************************
list<Particles>* ParticlesListsTemp( GridsPoints***** ptrArray_in, double*** ptrVolumeCellArray_in, double mi0)
{
    list<Particles>* listsPtrTemp = new list<Particles>;
    double N, mi_simu;
   
    
    // for bottom temp domain
    for( int i = 0; i < totalFace; i++)
    {
        for( int j = 1; j <= fieldsGridsSize; j++)
        {
            for( int k =1; k<= fieldsGridsSize; k++)
            {
                int s = 0;
                // number density
                N = ptrArray_in[i][j][k][s]->Density() / mi0;
                // mass of each simulation particle 
                mi_simu = N / iniParticleNumberPerCell *mi0 * ptrVolumeCellArray_in[j][k][s];
                for ( int t = 1; t <= iniParticleNumberPerCell; t++)
                {
                // calculate random position
                Vector3 temp1 = ptrArray_in[i][j][k][s]->Pos3();
                Vector3 temp2 = ptrArray_in[i][j][k][s+1]->Pos3();
                Vector3 temp = UniformDisVector3( ptrArray_in[i][j][k][s]->Pos3(),ptrArray_in[i][j][k][s+1]->Pos3());

                Vector3 vPos = ptrArray_in[i][j][k][s]->Pos3().PlusProduct(UniformDisVector3( ptrArray_in[i][j][k][s]->Pos3(),ptrArray_in[i][j][k][s+1]->Pos3()));
                vPos = vPos.PlusProduct( UniformDisVector3( ptrArray_in[i][j][k][s]->Pos3(),ptrArray_in[i][j][k+1][s]->Pos3()));
                vPos = vPos.PlusProduct( UniformDisVector3( ptrArray_in[i][j][k][s]->Pos3(),ptrArray_in[i][j+1][k][s]->Pos3()));

                uint_64 intPos = vPos.Uint_64_Trans();
                // calculate random velocity
                Vector3 vVel = Vector3( MaxwellDisV( ikT, ptrArray_in[i][j][k][s]->Vel3().x(), mi0),
                                            MaxwellDisV( ikT, ptrArray_in[i][j][k][s]->Vel3().y(), mi0),
                                            MaxwellDisV( ikT, ptrArray_in[i][j][k][s]->Vel3().z(), mi0));
                double mu_simu = MaxwellDisEnergy( ptrArray_in, intPos);
                // put the particles at the end of list
                Particles tempP= Particles(intPos, vVel, mi_simu, mu_simu);
                listsPtrTemp->push_back( tempP);           
                }
            }
        }
    }


    // for top temp domain
    for( int i = 0; i < totalFace; i++)
    {
        for( int j = 1; j <= fieldsGridsSize; j++)
        {
            for( int k =1; k<= fieldsGridsSize; k++)
            {
                int s = fieldsGridsSize - 1;
                // number density
                N = ptrArray_in[i][j][k][s]->Density() / mi0;
                // mass of each simulation particle 
                mi_simu = N / iniParticleNumberPerCell *mi0 * ptrVolumeCellArray_in[j][k][s];
                for ( int t = 1; t <= iniParticleNumberPerCell; t++)
                {
                // calculate random position
                Vector3 temp1 = ptrArray_in[i][j][k][s]->Pos3();
                Vector3 temp2 = ptrArray_in[i][j][k][s+1]->Pos3();
                Vector3 temp = UniformDisVector3( ptrArray_in[i][j][k][s]->Pos3(),ptrArray_in[i][j][k][s+1]->Pos3());

                Vector3 vPos = ptrArray_in[i][j][k][s]->Pos3().PlusProduct(UniformDisVector3( ptrArray_in[i][j][k][s]->Pos3(),ptrArray_in[i][j][k][s+1]->Pos3()));
                vPos = vPos.PlusProduct( UniformDisVector3( ptrArray_in[i][j][k][s]->Pos3(),ptrArray_in[i][j][k+1][s]->Pos3()));
                vPos = vPos.PlusProduct( UniformDisVector3( ptrArray_in[i][j][k][s]->Pos3(),ptrArray_in[i][j+1][k][s]->Pos3()));

                uint_64 intPos = vPos.Uint_64_Trans();
                // calculate random velocity
                Vector3 vVel = Vector3( MaxwellDisV( ikT, ptrArray_in[i][j][k][s]->Vel3().x(), mi0),
                                            MaxwellDisV( ikT, ptrArray_in[i][j][k][s]->Vel3().y(), mi0),
                                            MaxwellDisV( ikT, ptrArray_in[i][j][k][s]->Vel3().z(), mi0));
                double mu_simu = MaxwellDisEnergy( ptrArray_in, intPos);
                // put the particles at the end of list
                Particles tempP= Particles(intPos, vVel, mi_simu, mu_simu);
                listsPtrTemp->push_back( tempP);           
                }
            }
        }
    }


    return listsPtrTemp;
}

//************************************************************************
//************************************************************************
// Generate a grids of pointers which point to the 
// objects of GridsPoints in heap.
//
//************************************************************************
//************************************************************************

GridsPoints***** GridsCreation()
{ 
//    GridsPoints **a = new GridsPoints* [totalFace * (fieldsGridsSize+1) * (fieldsGridsSize+1) * (fieldsGridsSize+1)];
//    GridsPoints *ptrArray[totalFace][fieldsGridsSize+1][fieldsGridsSize+1][radialGridsSize+1];
    GridsPoints *****ptrArray;
    ptrArray = new GridsPoints****[totalFace];
    for( int face = 0; face < totalFace; face++)
    {
        ptrArray[face] = new GridsPoints***[fieldsGridsSize+3];
        for( int i = 0; i <= fieldsGridsSize+2; i++)
        {
            ptrArray[face][i] = new GridsPoints**[fieldsGridsSize+3];
            for( int j = 0; j <= fieldsGridsSize+2; j++)
            {
                ptrArray[face][i][j] = new GridsPoints*[fieldsGridsSize+1];
            }
        }
    }
    // face 0 (to us)
    for (int i = 1; i <= fieldsGridsSize+1; i++)
    {
        for(int j = 1; j <= fieldsGridsSize+1; j++)
        {
            for(int k = 0; k <= fieldsGridsSize; k++)
            {
                ptrArray[0][i][j][k] = new GridsPoints( 0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0);
                ptrArray[0][i][j][k]->InttoPos3( 0, i, j, k);                                        
                ptrArray[0][i][j][k]->XYZtoB(ptrArray[0][i][j][k]->Pos3());
                ptrArray[0][i][j][k]->XYZtoVel();
                ptrArray[0][i][j][k]->XYZtoE();
                ptrArray[0][i][j][k]->XYZtoDensity();
                ptrArray[0][i][j][k]->SetStopSign(0);
                ptrArray[0][i][j][k]->SetTemperature(0.0);
            }
        }
    }
    
    // face 1 (on the right)
    // share len with face 0
    for( int j = 1; j <= fieldsGridsSize+1; j++)
    {
        for( int k =0; k <= fieldsGridsSize; k++)
        {
            ptrArray[1][1][j][k] = ptrArray[0][fieldsGridsSize+1][j][k];
        }
    }

    for( int i = 2; i <= fieldsGridsSize+1; i++)
    {
        for( int j = 1; j <= fieldsGridsSize+1; j++)
        {
            for( int k = 0; k <= fieldsGridsSize; k++)
            {
                ptrArray[1][i][j][k] = new GridsPoints( 0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0);
                ptrArray[1][i][j][k]->InttoPos3( 1, i, j, k);                                        
                ptrArray[1][i][j][k]->XYZtoB(ptrArray[1][i][j][k]->Pos3());
                ptrArray[1][i][j][k]->XYZtoVel();
                ptrArray[1][i][j][k]->XYZtoE();     
                ptrArray[1][i][j][k]->XYZtoDensity();
                ptrArray[1][i][j][k]->SetStopSign(0);
                ptrArray[1][i][j][k]->SetTemperature(0.0);                                
            }
        }
    }
    // face 2 (on the top)
    // share len with face 0
    for( int i = 1; i <= fieldsGridsSize+1; i++)
    {
        for( int k = 0; k <= fieldsGridsSize; k++)
        {
            ptrArray[2][i][1][k] = ptrArray[0][i][fieldsGridsSize+1][k];
        }
    }
    // share len with face 1
    for( int j = 1; j <= fieldsGridsSize+1; j++)
    {
        for( int k = 0; k <= fieldsGridsSize; k++)
        {
            ptrArray[2][fieldsGridsSize+1][j][k] = ptrArray[1][j][fieldsGridsSize+1][k];
        }
    }
    for( int i = 1; i <= fieldsGridsSize; i++)
    {   
        for( int j = 2; j <= fieldsGridsSize+1; j++)
        {
            for( int k = 0; k <= fieldsGridsSize; k++)
            {
                ptrArray[2][i][j][k] = new GridsPoints( 0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0);
                ptrArray[2][i][j][k]->InttoPos3( 2, i, j, k);                                        
                ptrArray[2][i][j][k]->XYZtoB(ptrArray[2][i][j][k]->Pos3());
                ptrArray[2][i][j][k]->XYZtoVel();
                ptrArray[2][i][j][k]->XYZtoE();
                ptrArray[2][i][j][k]->XYZtoDensity();
                ptrArray[2][i][j][k]->SetStopSign(0);
                ptrArray[2][i][j][k]->SetTemperature(0.0);
            }
        }
    }

    // face 4 (on the left)
    // share len with face 0
    for( int j = 1; j <= fieldsGridsSize+1; j++)
    {
        for( int k = 0; k <= fieldsGridsSize; k++)
        {
            ptrArray[4][fieldsGridsSize+1][j][k] = ptrArray[0][1][j][k];
        }
    }
    // share len with face 2
    for( int i = 1; i <= fieldsGridsSize+1; i++)
    {
        for( int k = 0; k <= fieldsGridsSize; k++)
        {
            ptrArray[4][i][fieldsGridsSize+1][k] = ptrArray[2][1][fieldsGridsSize+2-i][k];
        }
    }

    for( int i = 1; i <= fieldsGridsSize; i++)
    {
        for( int j = 1; j <= fieldsGridsSize; j++)
        {
            for( int k = 0; k <= fieldsGridsSize; k++)
            {
                ptrArray[4][i][j][k] = new GridsPoints( 0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0);
                ptrArray[4][i][j][k]->InttoPos3( 4, i, j, k);                                        
                ptrArray[4][i][j][k]->XYZtoB(ptrArray[4][i][j][k]->Pos3());
                ptrArray[4][i][j][k]->XYZtoVel();
                ptrArray[4][i][j][k]->XYZtoE();
                ptrArray[4][i][j][k]->XYZtoDensity();
                ptrArray[4][i][j][k]->SetStopSign(0);
                ptrArray[4][i][j][k]->SetTemperature(0.0);
            }
        }
    }
    // face 5 (on the bottom)
    // share len with face 0
    for( int i = 1; i <= fieldsGridsSize+1; i++)
    {
        for( int k = 0; k <= fieldsGridsSize; k++)
        {
            ptrArray[5][i][fieldsGridsSize+1][k] = ptrArray[0][i][1][k];
        }
    }
    // share len with face 1
    for( int j = 1; j <= fieldsGridsSize+1; j++)
    {
        for( int k = 0; k <= fieldsGridsSize; k++)
        {
            ptrArray[5][fieldsGridsSize+1][j][k] = ptrArray[1][fieldsGridsSize+2-j][1][k];
        }
    }
    // share len with face 4
    for( int j = 1; j <= fieldsGridsSize+1; j++)
    {
        for( int k = 0; k <= fieldsGridsSize; k++)
        {
            ptrArray[5][1][j][k] = ptrArray[4][j][1][k];
        }
    }
    for( int i = 2; i <= fieldsGridsSize; i++)
    {
        for( int j = 1; j <= fieldsGridsSize; j++)
        {
            for( int k = 0; k <= fieldsGridsSize; k++)
            {
                ptrArray[5][i][j][k] = new GridsPoints( 0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0);
                ptrArray[5][i][j][k]->InttoPos3( 5, i, j, k);                                        
                ptrArray[5][i][j][k]->XYZtoB(ptrArray[5][i][j][k]->Pos3());
                ptrArray[5][i][j][k]->XYZtoVel();
                ptrArray[5][i][j][k]->XYZtoE();
                ptrArray[5][i][j][k]->XYZtoDensity();
                ptrArray[5][i][j][k]->SetStopSign(0);
                ptrArray[5][i][j][k]->SetTemperature(0.0);
            }
        }
    }
    // face 3
    // share len with face 1
    for( int j = 1; j <= fieldsGridsSize+1; j++)
    {
        for( int k = 0; k <= fieldsGridsSize; k++)
        {
            ptrArray[3][1][j][k] = ptrArray[1][fieldsGridsSize+1][j][k];
        }
    }
    // share len with face 2
    for( int i = 1; i <= fieldsGridsSize+1; i++)
    {
        for( int k = 0; k <= fieldsGridsSize; k++)
        {
            ptrArray[3][i][fieldsGridsSize+1][k] = ptrArray[2][fieldsGridsSize+2-i][fieldsGridsSize+1][k];
        }
    }
    // share len with face 4
    for( int j = 1; j <= fieldsGridsSize+1; j++)
    {
        for( int k = 0; k <= fieldsGridsSize; k++)
        {
            ptrArray[3][fieldsGridsSize+1][j][k] = ptrArray[4][1][j][k];
        }
    }

    // share len with face 5
    for( int i = 1; i <= fieldsGridsSize+1; i++)
    {
        for( int k = 0; k <= fieldsGridsSize; k++)
        {
            ptrArray[3][i][1][k] = ptrArray[5][fieldsGridsSize+2-i][1][k];
        }
    }
    for( int i = 2; i <= fieldsGridsSize; i++)
    {
        for( int j = 2; j <= fieldsGridsSize; j++)
        {
            for( int k = 0; k <= fieldsGridsSize; k++)
            {
                ptrArray[3][i][j][k] = new GridsPoints( 0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0);
                ptrArray[3][i][j][k]->InttoPos3( 3, i, j, k);                                        
                ptrArray[3][i][j][k]->XYZtoB(ptrArray[3][i][j][k]->Pos3());
                ptrArray[3][i][j][k]->XYZtoVel();
                ptrArray[3][i][j][k]->XYZtoE();
                ptrArray[3][i][j][k]->XYZtoDensity();
                ptrArray[3][i][j][k]->SetStopSign(0);
                ptrArray[3][i][j][k]->SetTemperature(0.0);
            }
        }
    }

    // info lens for adjant face
    //face 0 
    for( int i = 1; i <= fieldsGridsSize+1; i++)
    {
        for( int k = 0; k <= fieldsGridsSize; k++)
        {
        ptrArray[0][i][0][k] = ptrArray[5][i][fieldsGridsSize][k]; // bot
        ptrArray[0][fieldsGridsSize+2][i][k] = ptrArray[1][2][i][k]; // right
        ptrArray[0][i][fieldsGridsSize+2][k] = ptrArray[2][i][2][k]; // top
        ptrArray[0][0][i][k] = ptrArray[4][fieldsGridsSize][i][k];  // left
        }
    }
    //face 1
    for( int i = 1; i <= fieldsGridsSize+1; i++)
    {
        for( int k = 0; k <= fieldsGridsSize; k++)
        {
        ptrArray[1][i][0][k] = ptrArray[5][fieldsGridsSize][fieldsGridsSize+2-i][k]; // bot
        ptrArray[1][fieldsGridsSize+2][i][k] = ptrArray[3][2][i][k]; // right
        ptrArray[1][i][fieldsGridsSize+2][k] = ptrArray[2][fieldsGridsSize][i][k]; // top
        ptrArray[1][0][i][k] = ptrArray[0][fieldsGridsSize][i][k];  // left
        }
    }
    //face 2
    for( int i = 1; i <= fieldsGridsSize+1; i++)
    {
        for( int k = 0; k <= fieldsGridsSize; k++)
        {
        ptrArray[2][i][0][k] = ptrArray[0][i][fieldsGridsSize][k]; // bot
        ptrArray[2][fieldsGridsSize+2][i][k] = ptrArray[1][i][fieldsGridsSize][k]; // right
        ptrArray[2][i][fieldsGridsSize+2][k] = ptrArray[3][fieldsGridsSize+2-i][fieldsGridsSize][k]; // top
        ptrArray[2][0][i][k] = ptrArray[4][fieldsGridsSize+2-i][fieldsGridsSize][k];  // left
        }
    }
    //face 3
    for( int i = 1; i <= fieldsGridsSize+1; i++)
    {
        for( int k = 0; k <= fieldsGridsSize; k++)
        {
        ptrArray[3][i][0][k] = ptrArray[5][fieldsGridsSize+2-i][2][k]; // bot
        ptrArray[3][fieldsGridsSize+2][i][k] = ptrArray[4][2][i][k]; // right
        ptrArray[3][i][fieldsGridsSize+2][k] = ptrArray[2][fieldsGridsSize+2-i][fieldsGridsSize][k]; // top
        ptrArray[3][0][i][k] = ptrArray[1][fieldsGridsSize][i][k];  // left
        }
    }
    //face 4
    for( int i = 1; i <= fieldsGridsSize+1; i++)
    {
        for( int k = 0; k <= fieldsGridsSize; k++)
        {
        ptrArray[4][i][0][k] = ptrArray[5][2][i][k]; // bot
        ptrArray[4][fieldsGridsSize+2][i][k] = ptrArray[0][2][i][k]; // right
        ptrArray[4][i][fieldsGridsSize+2][k] = ptrArray[2][2][fieldsGridsSize+2-i][k]; // top
        ptrArray[4][0][i][k] = ptrArray[3][fieldsGridsSize][i][k];  // left
        }
    }   
    //face 5
    for( int i = 1; i <= fieldsGridsSize+1; i++)
    {
        for( int k = 0; k <= fieldsGridsSize; k++)
        {
        ptrArray[5][i][0][k] = ptrArray[3][fieldsGridsSize+2-i][2][k]; // bot
        ptrArray[5][fieldsGridsSize+2][i][k] = ptrArray[1][fieldsGridsSize+2-i][2][k]; // right
        ptrArray[5][i][fieldsGridsSize+2][k] = ptrArray[0][i][2][k]; // top
        ptrArray[5][0][i][k] = ptrArray[4][i][2][k];  // left
        }
    }    
    
    // Info for not used pointers at four corners of each face
    for( int face = 0; face < totalFace; face++)
    {
        for( int k = 0; k <= fieldsGridsSize; k++)
        {
            ptrArray[face][0][0][k] = new GridsPoints( 0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0);
            ptrArray[face][0][fieldsGridsSize+2][k] = new GridsPoints( 0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0);
            ptrArray[face][fieldsGridsSize+2][0][k] = new GridsPoints( 0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0);
            ptrArray[face][fieldsGridsSize+2][fieldsGridsSize+2][k] = new GridsPoints( 0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0);
        }
    }

    //////////// test of location for peak points 

    cout << "fieldsGridsSize " << fieldsGridsSize  << endl;
//    cout << "sizeof " << sizeof(ptrArray[0])/sizeof(ptrArray[0][0][0][0]) << endl;
    
/*    // face 0, bottom left
    cout <<"0 "<< ptrArray[0][1][1][0]<<" 5 "<< ptrArray[5][1][fieldsGridsSize+1][0] <<" 4 "<< ptrArray[4][fieldsGridsSize+1][1][0] << endl;
    // face 0, bottom right
    cout <<"0 "<< ptrArray[0][fieldsGridsSize+1][1][0]<<" 5 "<< ptrArray[5][fieldsGridsSize+1][fieldsGridsSize+1][0] <<" 1 "<< ptrArray[1][1][1][0] << endl;
    // face 0, top left
    cout <<"0 "<< ptrArray[0][1][fieldsGridsSize+1][0]<<" 2 "<< ptrArray[2][1][1][0] <<" 4 "<< ptrArray[4][fieldsGridsSize+1][fieldsGridsSize+1][0] << endl;
    // face 0, top right
    cout <<"0 "<< ptrArray[0][fieldsGridsSize+1][fieldsGridsSize+1][0]<<" 1 "<< ptrArray[1][1][fieldsGridsSize+1][0] <<" 2 "<< ptrArray[2][fieldsGridsSize+1][1][0] << endl;
    // face 3, bottom left
    cout <<"3 "<< ptrArray[3][1][1][0]<<" 1 "<< ptrArray[1][fieldsGridsSize+1][1][0] <<" 5 "<< ptrArray[5][fieldsGridsSize+1][1][0] << endl;
    // face 3, bottom right
    cout <<"3 "<< ptrArray[3][fieldsGridsSize+1][1][0]<<" 4 "<< ptrArray[4][1][1][0] <<" 5 "<< ptrArray[5][1][1][0] << endl;
    // face 3, top left
    cout <<"3 "<< ptrArray[3][1][fieldsGridsSize+1][0]<<" 1 "<< ptrArray[1][fieldsGridsSize+1][fieldsGridsSize+1][0] <<" 2 "<< ptrArray[2][fieldsGridsSize+1][fieldsGridsSize+1][0] << endl;
    // face 3, top right
    cout <<"3 "<< ptrArray[3][fieldsGridsSize+1][fieldsGridsSize+1][0]<<" 4 "<< ptrArray[4][1][fieldsGridsSize+1][0] <<" 2 "<< ptrArray[2][1][fieldsGridsSize+1][0] << endl;
    //face 0 
    for( int i =0; i < 6; i++)
    {   
        cout << "face " << i << endl;
        cout << "botleft " << ptrArray[i][0][1][0] << " " << ptrArray[i][1][0][0] << endl;
        cout << "botright" << ptrArray[i][fieldsGridsSize+2][1][0] << " " << ptrArray[i][fieldsGridsSize+1][0][0] << endl;
        cout << "topleft " << ptrArray[i][0][fieldsGridsSize+1][0] << " " << ptrArray[i][1][fieldsGridsSize+2][0] << endl;
        cout << "topright" << ptrArray[i][fieldsGridsSize+2][fieldsGridsSize+1][0] << " " << ptrArray[i][fieldsGridsSize+1][fieldsGridsSize+2][0] << endl << endl;        
    }
*/   
    return ptrArray;
}


//************************************************************************
//************************************************************************
// FUNCTION // Set up a matrix to store the curl E or B for Faraday's Law
// and for Ampere's Law, or the gradient of Pe. 
// The size of the matrix should be 1 smaller than 
// the size of gridspoints in main doman which is a cubic, which is [fsize+2].
// Therefore, it is [fsize+2][fsize+2][fsize]
// For each face, 8 corner cell should be excluded. ?
// Notice that the curl E or B is at the center of each cell.
// The data structure is array of Vector3, which is created in heap. Return
// a pointer(may not need to be a smart pointer), and would not need to 
// delete, or would be deleted as a smart pointer.
//************************************************************************
//************************************************************************
Vector3*** VectorCellField()
{
    static Vector3* mem_VectorCellField = new Vector3[ (fieldsGridsSize+2)*(fieldsGridsSize+2)*fieldsGridsSize];
    Vector3*** cellArray = new Vector3**[fieldsGridsSize+2];
    for( int i = 0; i < fieldsGridsSize+2; i++)
    {
        cellArray[i] = new Vector3*[fieldsGridsSize+2];
        for( int j = 0; j < fieldsGridsSize+2; j++)
        {
            cellArray[i][j] = new Vector3[fieldsGridsSize];
            cellArray[i][j] = mem_VectorCellField + i*(fieldsGridsSize+2)*(fieldsGridsSize)
                            + j* fieldsGridsSize;
            for( int k = 0; k< fieldsGridsSize; k++)
            {
                cellArray[i][j][k] = Vector3(0.0, 0.0, 0.0);
            }
        }
    }
    return cellArray;
/*
    Vector3**** cellArray;
    cellArray = new Vector3***[fieldsGridsSize+2];
    for( int i = 0; i < fieldsGridsSize+2; i++)
    {
        cellArray[i] = new Vector3**[fieldsGridsSize+2];
        for( int j=0; j < fieldsGridsSize+2; j++)
        {
            cellArray[i][j] = new Vector3*[fieldsGridsSize];
            for( int k = 0; k < fieldsGridsSize; k++)
            {
                cellArray[i][j][k] = new Vector3(0.0, 0.0, 0.0);
            }
        }
    }
    return cellArray;
    */
}

//************************************************************************
//************************************************************************
// FUNCTION 
// Value the matrix field using finite volume method, put in the pointer 
// of the MatrixField, value it, and return the pointer.
// Notice that the cell at corners should be absent in calculation.
//************************************************************************
//************************************************************************
Vector3*** ValueCurlField( Vector3*** curlArray_in, double*** ptrVolumeCellArray_in, GridsPoints***** ptrArray_in, int face_in, char field_in)
{
    for( int i = 0; i < fieldsGridsSize+2; i++)
    {
        for( int j = 0; j < fieldsGridsSize+2; j++)
        {
            for( int k = 0; k < fieldsGridsSize; k++)
            {
                if((i==0&&j==0)||
                            (i==0&&j==fieldsGridsSize+1)||
                            (i==fieldsGridsSize+1&&j==fieldsGridsSize+1)||
                            (i==fieldsGridsSize+1&&j==0))  
                            continue;


        /*        Vector3 test = AreaVectorL( ptrArray_in, face_in, i, j, k);
                std::cout << test.x() << " " << test.y() << " " << test.z() << " " << test.norm();
                int pause;
                std::cin >> pause;
        */
                // for each cell, calculate sum(n X B( on face)) and devided by
                // Volume to get the curl B at the center of cell
                Vector3 temp = AreaVectorL( ptrArray_in, face_in, i, j, k).CrossProduct(
                               FaceFieldVectorL(ptrArray_in, face_in, i, j, k, field_in));
                        temp = temp.PlusProduct(
                               AreaVectorR( ptrArray_in, face_in, i, j, k).CrossProduct(
                               FaceFieldVectorR(ptrArray_in, face_in, i, j, k, field_in)));
                        temp = temp.PlusProduct(
                               AreaVectorT( ptrArray_in, face_in, i, j, k).CrossProduct(
                               FaceFieldVectorT(ptrArray_in, face_in, i, j, k, field_in)));
                        temp = temp.PlusProduct(
                               AreaVectorBot( ptrArray_in, face_in, i, j, k).CrossProduct(
                               FaceFieldVectorBot(ptrArray_in, face_in, i, j, k, field_in)));
                        temp = temp.PlusProduct(
                               AreaVectorF( ptrArray_in, face_in, i, j, k).CrossProduct(
                               FaceFieldVectorF(ptrArray_in, face_in, i, j, k, field_in)));
                        temp = temp.PlusProduct(
                               AreaVectorBack( ptrArray_in, face_in, i, j, k).CrossProduct(
                               FaceFieldVectorBack(ptrArray_in, face_in, i, j, k, field_in)));
                double volumetemp = ptrVolumeCellArray_in[i][j][k];

                temp = temp.ScaleProduct( 1.0 / volumetemp);
                curlArray_in[i][j][k].SetVector3( temp); 
/*
                std::cout << face_in << field_in << i << j << k << " " << temp.norm() << std::endl;;
                int pause;
                if( temp.norm() !=0 ){
     //           std::cin >> pause;
                    }
 */           }
        }
    }
    return curlArray_in;
}

//************************************************************************
//************************************************************************
// FUNCTION 
// Value gradient field of Pe.
// gradientArray_in is in size of ( fsize+2 * fsize+2 * fsize) with vector3
// ptrVolumeCellArray is in size of ( fsize+2 * fsize+2 * fsize) with double
// Pe = n k T, in which n is the number density, k is the boltzmann constant, and T is the Te
//************************************************************************
//************************************************************************
Vector3*** ValueGradient(Vector3*** gradientArray_in, double*** ptrVolumeCellArray_in, GridsPoints***** ptrArray_in, int face_in, char char_in)
{
    for( int i = 0; i < fieldsGridsSize+2; i++)
    {
        for( int j = 0; j < fieldsGridsSize+2; j++)
        {
            for( int k = 0; k < fieldsGridsSize; k++)
            {
                if((i==0&&j==0)||
                            (i==0&&j==fieldsGridsSize+1)||
                            (i==fieldsGridsSize+1&&j==fieldsGridsSize+1)||
                            (i==fieldsGridsSize+1&&j==0))
                continue;

                Vector3 temp;

                if( char_in == 'P')
                {
                // for each cell, calculate sum of n(face vector) * densities and devided by
                // Volume to get the gradient at the center of cell  
                temp = AreaVectorL( ptrArray_in, face_in, i, j, k).ScaleProduct(
                           FaceNumberDensityL(ptrArray_in, face_in, i, j, k) * FaceTemperatureL(ptrArray_in, face_in, i, j, k) * boltzmann_k);
                temp = temp.PlusProduct(
                               AreaVectorR( ptrArray_in, face_in, i, j, k).ScaleProduct(
                               FaceNumberDensityR(ptrArray_in, face_in, i, j, k) * FaceTemperatureR(ptrArray_in, face_in, i, j, k) * boltzmann_k));
                temp = temp.PlusProduct(
                               AreaVectorT( ptrArray_in, face_in, i, j, k).ScaleProduct(
                               FaceNumberDensityT(ptrArray_in, face_in, i, j, k) * FaceTemperatureT(ptrArray_in, face_in, i, j, k) * boltzmann_k));
                temp = temp.PlusProduct(
                               AreaVectorBot( ptrArray_in, face_in, i, j, k).ScaleProduct(
                               FaceNumberDensityBot(ptrArray_in, face_in, i, j, k) * FaceTemperatureBot(ptrArray_in, face_in, i, j, k) * boltzmann_k));
                temp = temp.PlusProduct(
                               AreaVectorF( ptrArray_in, face_in, i, j, k).ScaleProduct(
                               FaceNumberDensityF(ptrArray_in, face_in, i, j, k) * FaceTemperatureF(ptrArray_in, face_in, i, j, k) * boltzmann_k));
                temp = temp.PlusProduct(
                               AreaVectorBack( ptrArray_in, face_in, i, j, k).ScaleProduct(
                               FaceNumberDensityBack(ptrArray_in, face_in, i, j, k) * FaceTemperatureBack(ptrArray_in, face_in, i, j, k) * boltzmann_k));
                double volumetemp = ptrVolumeCellArray_in[i][j][k];
                temp = temp.ScaleProduct(1/volumetemp);
   //             std::cout << volumetemp << " " ;
     //           std::cout << " --> " << temp.x() << " " << temp.y() << " " << temp.z() << std::endl;
                gradientArray_in[i][j][k].SetVector3(temp); 
                } 
                else
                {
                temp = AreaVectorL( ptrArray_in, face_in, i, j, k).ScaleProduct( FaceNormBL(ptrArray_in, face_in, i, j, k) );
                temp = temp.PlusProduct(
                               AreaVectorR( ptrArray_in, face_in, i, j, k).ScaleProduct( FaceNormBR(ptrArray_in, face_in, i, j, k) ));
                temp = temp.PlusProduct(
                               AreaVectorT( ptrArray_in, face_in, i, j, k).ScaleProduct( FaceNormBT(ptrArray_in, face_in, i, j, k) ));
                temp = temp.PlusProduct(
                               AreaVectorBot( ptrArray_in, face_in, i, j, k).ScaleProduct( FaceNormBBot(ptrArray_in, face_in, i, j, k) ));
                temp = temp.PlusProduct(
                               AreaVectorF( ptrArray_in, face_in, i, j, k).ScaleProduct( FaceNormBF(ptrArray_in, face_in, i, j, k) ));
                temp = temp.PlusProduct(
                               AreaVectorBack( ptrArray_in, face_in, i, j, k).ScaleProduct( FaceNormBBack(ptrArray_in, face_in, i, j, k) ));
                double volumetemp = ptrVolumeCellArray_in[i][j][k];
                temp = temp.ScaleProduct(1/volumetemp);
                gradientArray_in[i][j][k].SetVector3(temp); 
                }
                
            }
        }
    }
    return gradientArray_in;
}

//************************************************************************
//************************************************************************
// FUNCTION 
// UpdateVe3
//************************************************************************
//************************************************************************
void UpdateVe3( Vector3*** curlField_in, GridsPoints***** ptrArray_in, int face_in)
{
    for( int i = 1; i < fieldsGridsSize+2; i++)
    {
        for( int j = 1; j < fieldsGridsSize+2; j++)
        {
            for( int k = 1; k < fieldsGridsSize; k++)
            {
                Vector3 curl_field = Vector3(0.0, 0.0, 0.0);
                    
                if(i==1&&j==1) 
                {
                // gradPe at gridspoints
                curl_field = curlField_in[1][0][k-1].PlusProduct(
                                    curlField_in[1][1][k-1]);
                curl_field = curl_field.PlusProduct(           
                                    curlField_in[0][1][k-1]);
                curl_field = curl_field.PlusProduct(
                                    curlField_in[1][0][k]);
                curl_field = curl_field.PlusProduct(
                                    curlField_in[1][1][k]);
                curl_field = curl_field.PlusProduct(
                                    curlField_in[0][1][k]).ScaleProduct(1.0/6.0);            
                }
                else if(i==1&&j==fieldsGridsSize+1) 
                {
                
                curl_field = curlField_in[1][fieldsGridsSize+1][k-1].PlusProduct(
                                    curlField_in[1][fieldsGridsSize][k-1]);
                curl_field = curl_field.PlusProduct(           
                                    curlField_in[0][fieldsGridsSize][k-1]);
                curl_field = curl_field.PlusProduct(
                                    curlField_in[1][fieldsGridsSize+1][k]);
                curl_field = curl_field.PlusProduct(
                                    curlField_in[1][fieldsGridsSize][k]);
                curl_field = curl_field.PlusProduct(
                                    curlField_in[0][fieldsGridsSize][k]).ScaleProduct(1.0/6.0);
                }
                else if(i==fieldsGridsSize+1&&j==fieldsGridsSize+1)
                {
                curl_field = curlField_in[fieldsGridsSize][fieldsGridsSize+1][k-1].PlusProduct(
                                    curlField_in[fieldsGridsSize][fieldsGridsSize][k-1]);
                curl_field = curl_field.PlusProduct(           
                                    curlField_in[fieldsGridsSize+1][fieldsGridsSize][k-1]);
                curl_field = curl_field.PlusProduct(
                                    curlField_in[fieldsGridsSize][fieldsGridsSize+1][k]);
                curl_field = curl_field.PlusProduct(
                                    curlField_in[fieldsGridsSize][fieldsGridsSize][k]);
                curl_field = curl_field.PlusProduct(
                                    curlField_in[fieldsGridsSize+1][fieldsGridsSize][k]).ScaleProduct(1.0/6.0);
                }
                else if(i==fieldsGridsSize+1&&j==1) 
                {
                curl_field = curlField_in[fieldsGridsSize][0][k-1].PlusProduct(
                                    curlField_in[fieldsGridsSize][1][k-1]);
                curl_field = curl_field.PlusProduct(           
                                    curlField_in[fieldsGridsSize+1][1][k-1]);
                curl_field = curl_field.PlusProduct(
                                    curlField_in[fieldsGridsSize][0][k]);
                curl_field = curl_field.PlusProduct(
                                    curlField_in[fieldsGridsSize][1][k]);
                curl_field = curl_field.PlusProduct(
                                    curlField_in[fieldsGridsSize+1][1][k]).ScaleProduct(1.0/6.0);                
                } else
                {
                curl_field = curlField_in[i-1][j-1][k-1].PlusProduct(
                                    curlField_in[i][j-1][k-1]);
                curl_field = curl_field.PlusProduct(           
                                    curlField_in[i-1][j][k-1]);
                curl_field = curl_field.PlusProduct(
                                    curlField_in[i][j][k-1]);
                curl_field = curl_field.PlusProduct(
                                    curlField_in[i-1][j-1][k]);
                curl_field = curl_field.PlusProduct(
                                    curlField_in[i][j-1][k]);
                curl_field = curl_field.PlusProduct(
                                    curlField_in[i-1][j][k]);    
                curl_field = curl_field.PlusProduct(
                                    curlField_in[i][j][k]).ScaleProduct(1.0/8.0);
                }
                ptrArray_in[face_in][i][j][k]->updateve3(curl_field);
            }
        }
    }
}


//************************************************************************
//************************************************************************
// FUNCTION
// As in the updating curlField and gradientPe array, some variables are
// repeating calculating, it is suitable to put them in one function.
// Therefore, we need three matrix of curlB, curlE, and gradientPe. 
// Assume they are curlB, curlE and gradPe, respectively.
//************************************************************************
//************************************************************************
void updateCellMatrix(Vector3**** curlB_in, Vector3**** curlE_in,
                      Vector3**** gradPe_in, GridsPoints***** ptrArray_in, int face_in)
{
     for( int i = 0; i < fieldsGridsSize+2; i++)
    {
        for( int j = 0; j < fieldsGridsSize+2; j++)
        {
            for( int k = 0; k < fieldsGridsSize+2; k++)
            {
                if((i==0&&j==0)||
                            (i==0&&j==fieldsGridsSize+1)||
                            (i==fieldsGridsSize+1&&j==fieldsGridsSize+1)||
                            (i==fieldsGridsSize+1&&j==0) )
                continue;
                // Volume
                double volumetemp = CellVolume(ptrArray_in, face_in, i, j, k);
                // AreaVectors
                Vector3 nL = AreaVectorL( ptrArray_in, face_in, i, j, k);
                Vector3 nR = AreaVectorR( ptrArray_in, face_in, i, j, k);
                Vector3 nT = AreaVectorT( ptrArray_in, face_in, i, j, k);
                Vector3 nBot= AreaVectorBot( ptrArray_in, face_in, i, j, k);
                Vector3 nF = AreaVectorF( ptrArray_in, face_in, i, j, k);
                Vector3 nBack = AreaVectorBack( ptrArray_in, face_in, i, j, k);

                // for each cell, calculate sum of n(vector) * densities and devided by
                // Volume to get the gradient at the center of cell
                
                Vector3 tempGrad = nL.ScaleProduct(
                               FaceDensityL(ptrArray_in, face_in, i, j, k));
                        tempGrad = tempGrad.PlusProduct(
                               nR.ScaleProduct(
                               FaceDensityR(ptrArray_in, face_in, i, j, k)));
                        tempGrad = tempGrad.PlusProduct(
                               nT.ScaleProduct(
                               FaceDensityT(ptrArray_in, face_in, i, j, k)));
                        tempGrad = tempGrad.PlusProduct(
                               nBot.ScaleProduct(
                               FaceDensityBot(ptrArray_in, face_in, i, j, k)));
                        tempGrad = tempGrad.PlusProduct(
                               nF.ScaleProduct(
                               FaceDensityF(ptrArray_in, face_in, i, j, k)));
                        tempGrad = tempGrad.PlusProduct(
                               nBack.ScaleProduct(
                               FaceDensityBack(ptrArray_in, face_in, i, j, k)));

                
                tempGrad = tempGrad.ScaleProduct(1/volumetemp);
                gradPe_in[i][j][k]->SetVector3(tempGrad); 

                Vector3 tempCurlB = nL.CrossProduct(
                               FaceFieldVectorL(ptrArray_in, face_in, i, j, k, 'B'));
                        tempCurlB = tempCurlB.PlusProduct(
                               nR.CrossProduct(
                               FaceFieldVectorR(ptrArray_in, face_in, i, j, k, 'B')));
                        tempCurlB = tempCurlB.PlusProduct(
                               nT.CrossProduct(
                               FaceFieldVectorT(ptrArray_in, face_in, i, j, k, 'B')));
                        tempCurlB = tempCurlB.PlusProduct(
                               nBot.CrossProduct(
                               FaceFieldVectorBot(ptrArray_in, face_in, i, j, k, 'B')));
                        tempCurlB = tempCurlB.PlusProduct(
                               nF.CrossProduct(
                               FaceFieldVectorF(ptrArray_in, face_in, i, j, k, 'B')));
                        tempCurlB = tempCurlB.PlusProduct(
                               nBack.CrossProduct(
                               FaceFieldVectorBack(ptrArray_in, face_in, i, j, k, 'B')));

                tempCurlB = tempCurlB.ScaleProduct(1/volumetemp);
                curlB_in[i][j][k]->SetVector3(tempCurlB); 

                Vector3 tempCurlE = nL.CrossProduct(
                               FaceFieldVectorL(ptrArray_in, face_in, i, j, k, 'E'));
                        tempCurlE = tempCurlE.PlusProduct(
                               nR.CrossProduct(
                               FaceFieldVectorR(ptrArray_in, face_in, i, j, k, 'E')));
                        tempCurlE = tempCurlE.PlusProduct(
                               nT.CrossProduct(
                               FaceFieldVectorT(ptrArray_in, face_in, i, j, k, 'E')));
                        tempCurlE = tempCurlE.PlusProduct(
                               nBot.CrossProduct(
                               FaceFieldVectorBot(ptrArray_in, face_in, i, j, k, 'E')));
                        tempCurlE = tempCurlE.PlusProduct(
                               nF.CrossProduct(
                               FaceFieldVectorF(ptrArray_in, face_in, i, j, k, 'E')));
                        tempCurlE = tempCurlE.PlusProduct(
                               nBack.CrossProduct(
                               FaceFieldVectorBack(ptrArray_in, face_in, i, j, k, 'E')));

                tempCurlE = tempCurlE.ScaleProduct(1/volumetemp);
                curlE_in[i][j][k]->SetVector3(tempCurlE); 
            }
        }
    }
}


//************************************************************************
//************************************************************************
// FUNCTION
// Each calculation are on the girds.
// Update E at grids for ve = vi ( no current)
// Update E at grids for ve ( with current)
//************************************************************************
//************************************************************************
void UpdateE3( Vector3*** gradPe_in, GridsPoints***** ptrArray_in, int face_in)
{


    for( int i = 1; i < fieldsGridsSize+2; i++)
    {
        for( int j = 1; j < fieldsGridsSize+2; j++)
        {
            for( int k = 1; k < fieldsGridsSize; k++)
            {
                Vector3 tempGradPe = Vector3(0.0, 0.0, 0.0);
                    
                if(i==1&&j==1) 
                {
                // gradPe at gridspoints
                tempGradPe = gradPe_in[1][0][k-1].PlusProduct(
                                    gradPe_in[1][1][k-1]);
                tempGradPe = tempGradPe.PlusProduct(           
                                    gradPe_in[0][1][k-1]);
                tempGradPe = tempGradPe.PlusProduct(
                                    gradPe_in[1][0][k]);
                tempGradPe = tempGradPe.PlusProduct(
                                    gradPe_in[1][1][k]);
                tempGradPe = tempGradPe.PlusProduct(
                                    gradPe_in[0][1][k]).ScaleProduct(1.0/6.0);            
                }
                else if(i==1&&j==fieldsGridsSize+1) 
                {
                
                tempGradPe = gradPe_in[1][fieldsGridsSize+1][k-1].PlusProduct(
                                    gradPe_in[1][fieldsGridsSize][k-1]);
                tempGradPe = tempGradPe.PlusProduct(           
                                    gradPe_in[0][fieldsGridsSize][k-1]);
                tempGradPe = tempGradPe.PlusProduct(
                                    gradPe_in[1][fieldsGridsSize+1][k]);
                tempGradPe = tempGradPe.PlusProduct(
                                    gradPe_in[1][fieldsGridsSize][k]);
                tempGradPe = tempGradPe.PlusProduct(
                                    gradPe_in[0][fieldsGridsSize][k]).ScaleProduct(1.0/6.0);
                }
                else if(i==fieldsGridsSize+1&&j==fieldsGridsSize+1)
                {
                tempGradPe = gradPe_in[fieldsGridsSize][fieldsGridsSize+1][k-1].PlusProduct(
                                    gradPe_in[fieldsGridsSize][fieldsGridsSize][k-1]);
                tempGradPe = tempGradPe.PlusProduct(           
                                    gradPe_in[fieldsGridsSize+1][fieldsGridsSize][k-1]);
                tempGradPe = tempGradPe.PlusProduct(
                                    gradPe_in[fieldsGridsSize][fieldsGridsSize+1][k]);
                tempGradPe = tempGradPe.PlusProduct(
                                    gradPe_in[fieldsGridsSize][fieldsGridsSize][k]);
                tempGradPe = tempGradPe.PlusProduct(
                                    gradPe_in[fieldsGridsSize+1][fieldsGridsSize][k]).ScaleProduct(1.0/6.0);
                }
                else if(i==fieldsGridsSize+1&&j==1) 
                {
                tempGradPe = gradPe_in[fieldsGridsSize][0][k-1].PlusProduct(
                                    gradPe_in[fieldsGridsSize][1][k-1]);
                tempGradPe = tempGradPe.PlusProduct(           
                                    gradPe_in[fieldsGridsSize+1][1][k-1]);
                tempGradPe = tempGradPe.PlusProduct(
                                    gradPe_in[fieldsGridsSize][0][k]);
                tempGradPe = tempGradPe.PlusProduct(
                                    gradPe_in[fieldsGridsSize][1][k]);
                tempGradPe = tempGradPe.PlusProduct(
                                    gradPe_in[fieldsGridsSize+1][1][k]).ScaleProduct(1.0/6.0);                
                } else
                {
                tempGradPe = gradPe_in[i-1][j-1][k-1].PlusProduct(
                                    gradPe_in[i][j-1][k-1]);
                tempGradPe = tempGradPe.PlusProduct(           
                                    gradPe_in[i-1][j][k-1]);
                tempGradPe = tempGradPe.PlusProduct(
                                    gradPe_in[i][j][k-1]);
                tempGradPe = tempGradPe.PlusProduct(
                                    gradPe_in[i-1][j-1][k]);
                tempGradPe = tempGradPe.PlusProduct(
                                    gradPe_in[i][j-1][k]);
                tempGradPe = tempGradPe.PlusProduct(
                                    gradPe_in[i-1][j][k]);    
                tempGradPe = tempGradPe.PlusProduct(
                                    gradPe_in[i][j][k]).ScaleProduct(1.0/8.0);
                }
                // update E
                ptrArray_in[face_in][i][j][k]->updateE(tempGradPe);
            }
        }
    }        
}


//************************************************************************
//************************************************************************
// FUNCTION
// Each calculation are on the girds.
// Step 1: Setup Ve for each gridpoint
// Step 2: Update E at each gridpoints
//************************************************************************
//************************************************************************
void update__( Vector3**** curlB_in, Vector3**** gradPe_in, GridsPoints***** ptrArray_in, int face_in)
{
    Vector3 tempCurlB;
    Vector3 tempGradPe;

    for( int i = 1; i < fieldsGridsSize+2; i++)
    {
        for( int j = 1; j < fieldsGridsSize+2; j++)
        {
            for( int k = 1; k < fieldsGridsSize; k++)
            {
                // curl B at gridspoint
                if(i==1&&j==1) 
                {
                tempCurlB = curlB_in[1][0][k-1]->PlusProduct(
                                    *curlB_in[1][1][k-1]);
                tempCurlB = tempCurlB.PlusProduct(           
                                    *curlB_in[0][1][k-1]);
                tempCurlB = tempCurlB.PlusProduct(
                                    *curlB_in[1][0][k]);
                tempCurlB = tempCurlB.PlusProduct(
                                    *curlB_in[1][1][k]);
                tempCurlB = tempCurlB.PlusProduct(
                                    *curlB_in[0][1][k]).ScaleProduct(1/6);
                // gradPe at gridspoints
                tempGradPe = gradPe_in[1][0][k-1]->PlusProduct(
                                    *gradPe_in[1][1][k-1]);
                tempGradPe = tempGradPe.PlusProduct(           
                                    *gradPe_in[0][1][k-1]);
                tempGradPe = tempGradPe.PlusProduct(
                                    *gradPe_in[1][0][k]);
                tempGradPe = tempGradPe.PlusProduct(
                                    *gradPe_in[1][1][k]);
                tempGradPe = tempGradPe.PlusProduct(
                                    *gradPe_in[0][1][k]).ScaleProduct(1/6);            
                }
                else if(i==1&&j==fieldsGridsSize+1) 
                {
                tempCurlB = curlB_in[1][fieldsGridsSize+1][k-1]->PlusProduct(
                                    *curlB_in[1][fieldsGridsSize][k-1]);
                tempCurlB = tempCurlB.PlusProduct(           
                                    *curlB_in[0][fieldsGridsSize][k-1]);
                tempCurlB = tempCurlB.PlusProduct(
                                    *curlB_in[1][fieldsGridsSize+1][k]);
                tempCurlB = tempCurlB.PlusProduct(
                                    *curlB_in[1][fieldsGridsSize][k]);
                tempCurlB = tempCurlB.PlusProduct(
                                    *curlB_in[0][fieldsGridsSize][k]).ScaleProduct(1/6);
                
                tempGradPe = gradPe_in[1][fieldsGridsSize+1][k-1]->PlusProduct(
                                    *gradPe_in[1][fieldsGridsSize][k-1]);
                tempGradPe = tempGradPe.PlusProduct(           
                                    *gradPe_in[0][fieldsGridsSize][k-1]);
                tempGradPe = tempGradPe.PlusProduct(
                                    *gradPe_in[1][fieldsGridsSize+1][k]);
                tempGradPe = tempGradPe.PlusProduct(
                                    *gradPe_in[1][fieldsGridsSize][k]);
                tempGradPe = tempGradPe.PlusProduct(
                                    *gradPe_in[0][fieldsGridsSize][k]).ScaleProduct(1/6);
                }
                else if(i==fieldsGridsSize+1&&j==fieldsGridsSize+1)
                {
                tempCurlB = curlB_in[fieldsGridsSize][fieldsGridsSize+1][k-1]->PlusProduct(
                                    *curlB_in[fieldsGridsSize][fieldsGridsSize][k-1]);
                tempCurlB = tempCurlB.PlusProduct(           
                                    *curlB_in[fieldsGridsSize+1][fieldsGridsSize][k-1]);
                tempCurlB = tempCurlB.PlusProduct(
                                    *curlB_in[fieldsGridsSize][fieldsGridsSize+1][k]);
                tempCurlB = tempCurlB.PlusProduct(
                                    *curlB_in[fieldsGridsSize][fieldsGridsSize][k]);
                tempCurlB = tempCurlB.PlusProduct(
                                    *curlB_in[fieldsGridsSize+1][fieldsGridsSize][k]).ScaleProduct(1/6);
                
                tempGradPe = gradPe_in[fieldsGridsSize][fieldsGridsSize+1][k-1]->PlusProduct(
                                    *gradPe_in[fieldsGridsSize][fieldsGridsSize][k-1]);
                tempGradPe = tempGradPe.PlusProduct(           
                                    *gradPe_in[fieldsGridsSize+1][fieldsGridsSize][k-1]);
                tempGradPe = tempGradPe.PlusProduct(
                                    *gradPe_in[fieldsGridsSize][fieldsGridsSize+1][k]);
                tempGradPe = tempGradPe.PlusProduct(
                                    *gradPe_in[fieldsGridsSize][fieldsGridsSize][k]);
                tempGradPe = tempGradPe.PlusProduct(
                                    *gradPe_in[fieldsGridsSize+1][fieldsGridsSize][k]).ScaleProduct(1/6);
                }
                else if(i==fieldsGridsSize+1&&j==1) 
                {
                tempCurlB = curlB_in[fieldsGridsSize][0][k-1]->PlusProduct(
                                    *curlB_in[fieldsGridsSize][1][k-1]);
                tempCurlB = tempCurlB.PlusProduct(           
                                    *curlB_in[fieldsGridsSize+1][1][k-1]);
                tempCurlB = tempCurlB.PlusProduct(
                                    *curlB_in[fieldsGridsSize][0][k]);
                tempCurlB = tempCurlB.PlusProduct(
                                    *curlB_in[fieldsGridsSize][1][k]);
                tempCurlB = tempCurlB.PlusProduct(
                                    *curlB_in[fieldsGridsSize+1][1][k]).ScaleProduct(1/6);
                
                tempGradPe = gradPe_in[fieldsGridsSize][0][k-1]->PlusProduct(
                                    *gradPe_in[fieldsGridsSize][1][k-1]);
                tempGradPe = tempGradPe.PlusProduct(           
                                    *gradPe_in[fieldsGridsSize+1][1][k-1]);
                tempGradPe = tempGradPe.PlusProduct(
                                    *gradPe_in[fieldsGridsSize][0][k]);
                tempGradPe = tempGradPe.PlusProduct(
                                    *gradPe_in[fieldsGridsSize][1][k]);
                tempGradPe = tempGradPe.PlusProduct(
                                    *gradPe_in[fieldsGridsSize+1][1][k]).ScaleProduct(1/6);                
                } else
                {
                tempCurlB = curlB_in[i-1][j-1][k-1]->PlusProduct(
                                    *curlB_in[i][j-1][k-1]);
                tempCurlB = tempCurlB.PlusProduct(           
                                    *curlB_in[i-1][j][k-1]);
                tempCurlB = tempCurlB.PlusProduct(
                                    *curlB_in[i][j][k-1]);
                tempCurlB = tempCurlB.PlusProduct(
                                    *curlB_in[i-1][j-1][k]);
                tempCurlB = tempCurlB.PlusProduct(
                                    *curlB_in[i][j-1][k]);
                tempCurlB = tempCurlB.PlusProduct(
                                    *curlB_in[i-1][j][k]);    
                tempCurlB = tempCurlB.PlusProduct(
                                    *curlB_in[i][j][k]).ScaleProduct(1/8);

                tempGradPe = gradPe_in[i-1][j-1][k-1]->PlusProduct(
                                    *gradPe_in[i][j-1][k-1]);
                tempGradPe = tempGradPe.PlusProduct(           
                                    *gradPe_in[i-1][j][k-1]);
                tempGradPe = tempGradPe.PlusProduct(
                                    *gradPe_in[i][j][k-1]);
                tempGradPe = tempGradPe.PlusProduct(
                                    *gradPe_in[i-1][j-1][k]);
                tempGradPe = tempGradPe.PlusProduct(
                                    *gradPe_in[i][j-1][k]);
                tempGradPe = tempGradPe.PlusProduct(
                                    *gradPe_in[i-1][j][k]);    
                tempGradPe = tempGradPe.PlusProduct(
                                    *gradPe_in[i][j][k]).ScaleProduct(1/8);
                }
                Vector3 Ve = ptrArray_in[face_in][i][j][k]->Vel3().MinusProduct(
                                    tempCurlB.ScaleProduct(1/alpha/
                                    ptrArray_in[face_in][i][j][k]->Density()));
                // update E
                ptrArray_in[face_in][i][j][k]->updateE(tempGradPe);
            }
        }
    }        
}


//************************************************************************
//************************************************************************
// FUNCTION
// UpdateB3 vased on faraday's law
//************************************************************************
//************************************************************************
void UpdateB3( Vector3*** curlField_in, GridsPoints***** ptrArray_in, int face_in)
{


    for( int i = 1; i < fieldsGridsSize+2; i++)
    {
        for( int j = 1; j < fieldsGridsSize+2; j++)
        {
            for( int k = 1; k < fieldsGridsSize; k++)
            {
                Vector3 curl_field = Vector3(0.0, 0.0, 0.0);
                    
                if(i==1&&j==1) 
                {
                // gradPe at gridspoints
                curl_field = curlField_in[1][0][k-1].PlusProduct(
                                    curlField_in[1][1][k-1]);
                curl_field = curl_field.PlusProduct(           
                                    curlField_in[0][1][k-1]);
                curl_field = curl_field.PlusProduct(
                                    curlField_in[1][0][k]);
                curl_field = curl_field.PlusProduct(
                                    curlField_in[1][1][k]);
                curl_field = curl_field.PlusProduct(
                                    curlField_in[0][1][k]).ScaleProduct(1.0/6.0);            
                }
                else if(i==1&&j==fieldsGridsSize+1) 
                {
                
                curl_field = curlField_in[1][fieldsGridsSize+1][k-1].PlusProduct(
                                    curlField_in[1][fieldsGridsSize][k-1]);
                curl_field = curl_field.PlusProduct(           
                                    curlField_in[0][fieldsGridsSize][k-1]);
                curl_field = curl_field.PlusProduct(
                                    curlField_in[1][fieldsGridsSize+1][k]);
                curl_field = curl_field.PlusProduct(
                                    curlField_in[1][fieldsGridsSize][k]);
                curl_field = curl_field.PlusProduct(
                                    curlField_in[0][fieldsGridsSize][k]).ScaleProduct(1.0/6.0);
                }
                else if(i==fieldsGridsSize+1&&j==fieldsGridsSize+1)
                {
                curl_field = curlField_in[fieldsGridsSize][fieldsGridsSize+1][k-1].PlusProduct(
                                    curlField_in[fieldsGridsSize][fieldsGridsSize][k-1]);
                curl_field = curl_field.PlusProduct(           
                                    curlField_in[fieldsGridsSize+1][fieldsGridsSize][k-1]);
                curl_field = curl_field.PlusProduct(
                                    curlField_in[fieldsGridsSize][fieldsGridsSize+1][k]);
                curl_field = curl_field.PlusProduct(
                                    curlField_in[fieldsGridsSize][fieldsGridsSize][k]);
                curl_field = curl_field.PlusProduct(
                                    curlField_in[fieldsGridsSize+1][fieldsGridsSize][k]).ScaleProduct(1.0/6.0);
                }
                else if(i==fieldsGridsSize+1&&j==1) 
                {
                curl_field = curlField_in[fieldsGridsSize][0][k-1].PlusProduct(
                                    curlField_in[fieldsGridsSize][1][k-1]);
                curl_field = curl_field.PlusProduct(           
                                    curlField_in[fieldsGridsSize+1][1][k-1]);
                curl_field = curl_field.PlusProduct(
                                    curlField_in[fieldsGridsSize][0][k]);
                curl_field = curl_field.PlusProduct(
                                    curlField_in[fieldsGridsSize][1][k]);
                curl_field = curl_field.PlusProduct(
                                    curlField_in[fieldsGridsSize+1][1][k]).ScaleProduct(1.0/6.0);                
                } else
                {
                curl_field = curlField_in[i-1][j-1][k-1].PlusProduct(
                                    curlField_in[i][j-1][k-1]);
                curl_field = curl_field.PlusProduct(           
                                    curlField_in[i-1][j][k-1]);
                curl_field = curl_field.PlusProduct(
                                    curlField_in[i][j][k-1]);
                curl_field = curl_field.PlusProduct(
                                    curlField_in[i-1][j-1][k]);
                curl_field = curl_field.PlusProduct(
                                    curlField_in[i][j-1][k]);
                curl_field = curl_field.PlusProduct(
                                    curlField_in[i-1][j][k]);    
                curl_field = curl_field.PlusProduct(
                                    curlField_in[i][j][k]).ScaleProduct(1.0/8.0);
                }
                // update B
                ptrArray_in[face_in][i][j][k]->updatedB(curl_field);
            }
        }
    }        
}


//************************************************************************
//************************************************************************
// FUNCTION
// Update info in the grids due to the info of particles and related 
// weighting
// input: ptrArray_in for the grids address
//        ptrParticlesLists for the particles list address
//************************************************************************
//************************************************************************
void UpdateInfoGrids( GridsPoints***** ptrArray_in, 
                      list<Particles>* ptrParticlesList_H_in, 
                      list<Particles>* ptrParticlesList_He_in,
                      list<Particles>* ptrParticlesList_O_in,
                      list<Particles>* ptrParticlesListTemp_H_in,
                      list<Particles>* ptrParticlesListTemp_He_in,
                      list<Particles>* ptrParticlesListTemp_O_in,
                      double*** ptrVolumeGridArray_in,
                      int timeline_in, int updateInfoPeriod_in)
{
    if( timeline_in == 0 || (timeline_in - 1) % updateInfoPeriod_in == 0)   // timeline_in should start with 1
    {
        // reset, clear the previous value: density and velocity
        for( int face = 0; face < totalFace; face++)
        {
            for( int i = 1; i < fieldsGridsSize+2; i++)
            {
                for( int j = 1; j < fieldsGridsSize+2; j++)
                {
                    for( int k = 1; k < fieldsGridsSize; k++)
                    {  
                        ptrArray_in[face][i][j][k]->ResetParameters();
                    }
                }
            }
        }
    }

// For H particles in main domain    
    for( list<Particles>::iterator iter= ptrParticlesList_H_in->begin(); iter!=ptrParticlesList_H_in->end(); ++iter)
    {
        // locate the particle
        Particles temp = *iter;
        // get the weighting info of this particle
        struct structg tempStr = temp.InttoStrp1();
        // get the info of mass(weight of each simulation particle)
        double tempMass = temp.WeightMi();

/*      // std::cout << temp.VelParticles().x() << " " << temp.VelParticles().y() << " " << temp.VelParticles().z() << std::endl;
        // std::cout << temp.WeightMi() << " <== " << std::endl;
        // int pause;
        // std::cin >>pause;
        // std::cout << "checkUp "<< tempStr.face << " " << tempStr.ig << " " << tempStr.jg << " " << tempStr.kg << " " << tempStr.iw << " " << tempStr.jw << " " << tempStr.kw << " " ;
        // std::cout << "str.v "<< tempStr.vx << " " << tempStr.vy << " " << tempStr.vz << " " ;
        // std::cout << "tempVel "<< tempVel.x() << " " << tempVel.y() << " " << tempVel.z() << " ";
        // std::cout << "tempMass "<< tempMass << std::endl;
        // int pause ;
        // std::cin >> pause;
        // std::cout << tempStr.face <<  tempStr.ig+1 << tempStr.jg+1 << tempStr.kg << " " << tempStr.iw << tempStr.jw << tempStr.kw << " mass " << tempMass << " xxxx "; ;
        // get the info of velocity
*/        
        Vector3 tempVel = temp.VelParticles();
        // update density and velocity in the grids
        ptrArray_in[tempStr.face][tempStr.ig+1][tempStr.jg+1][tempStr.kg]->UpdateDueToWgt( cellSize1 - tempStr.iw, cellSize1 - tempStr.jw, cellSize1 - tempStr.kw, tempMass, tempVel, 1);
        ptrArray_in[tempStr.face][tempStr.ig+2][tempStr.jg+1][tempStr.kg]->UpdateDueToWgt( tempStr.iw + 1,         cellSize1 - tempStr.jw, cellSize1 - tempStr.kw, tempMass, tempVel, 1);
        ptrArray_in[tempStr.face][tempStr.ig+1][tempStr.jg+2][tempStr.kg]->UpdateDueToWgt( cellSize1 - tempStr.iw, tempStr.jw + 1,         cellSize1 - tempStr.kw, tempMass, tempVel, 1);
        ptrArray_in[tempStr.face][tempStr.ig+2][tempStr.jg+2][tempStr.kg]->UpdateDueToWgt( tempStr.iw + 1,         tempStr.jw + 1,         cellSize1 - tempStr.kw, tempMass, tempVel, 1);
        
        ptrArray_in[tempStr.face][tempStr.ig+1][tempStr.jg+1][tempStr.kg+1]->UpdateDueToWgt( cellSize1 - tempStr.iw, cellSize1 - tempStr.jw, tempStr.kw + 1,         tempMass, tempVel, 1);
        ptrArray_in[tempStr.face][tempStr.ig+2][tempStr.jg+1][tempStr.kg+1]->UpdateDueToWgt( tempStr.iw + 1,         cellSize1 - tempStr.jw, tempStr.kw + 1,         tempMass, tempVel, 1);
        ptrArray_in[tempStr.face][tempStr.ig+1][tempStr.jg+2][tempStr.kg+1]->UpdateDueToWgt( cellSize1 - tempStr.iw, tempStr.jw + 1,         tempStr.kw + 1,         tempMass, tempVel, 1);
        ptrArray_in[tempStr.face][tempStr.ig+2][tempStr.jg+2][tempStr.kg+1]->UpdateDueToWgt( tempStr.iw + 1,         tempStr.jw + 1,         tempStr.kw + 1,         tempMass, tempVel, 1);
/*     
        std::cout << ptrArray_in[tempStr.face][tempStr.ig+1][tempStr.jg+1][tempStr.kg]->Density() << " "
        << ptrArray_in[tempStr.face][tempStr.ig+2][tempStr.jg+1][tempStr.kg]->Density() << " "
        << ptrArray_in[tempStr.face][tempStr.ig+1][tempStr.jg+2][tempStr.kg]->Density() << " "
        << ptrArray_in[tempStr.face][tempStr.ig+2][tempStr.jg+2][tempStr.kg]->Density() << " nnnn "
        << ptrArray_in[tempStr.face][tempStr.ig+1][tempStr.jg+1][tempStr.kg+1]->Density() << "  "
        << ptrArray_in[tempStr.face][tempStr.ig+2][tempStr.jg+1][tempStr.kg+1]->Density() << " "
        << ptrArray_in[tempStr.face][tempStr.ig+1][tempStr.jg+2][tempStr.kg+1]->Density() << " "
        << ptrArray_in[tempStr.face][tempStr.ig+2][tempStr.jg+2][tempStr.kg+1]->Density() << std::endl;
*/
        // int pause ;
        // std::cin >> pause;
    }

// For He particles in main domain    
    for( list<Particles>::iterator iter= ptrParticlesList_He_in->begin(); iter!=ptrParticlesList_He_in->end(); ++iter)
    {
        // locate the particle
        Particles temp = *iter;
        // get the weighting info of this particle
        struct structg tempStr = temp.InttoStrp1();
        // get the info of mass(weight of each simulation particle)
        double tempMass = temp.WeightMi();
        // get the info of velocity
        Vector3 tempVel = temp.VelParticles();
        // update density and velocity in the grids
        ptrArray_in[tempStr.face][tempStr.ig+1][tempStr.jg+1][tempStr.kg]->UpdateDueToWgt( cellSize1 - tempStr.iw, cellSize1 - tempStr.jw, cellSize1 - tempStr.kw, tempMass, tempVel, 4);
        ptrArray_in[tempStr.face][tempStr.ig+2][tempStr.jg+1][tempStr.kg]->UpdateDueToWgt( tempStr.iw + 1,         cellSize1 - tempStr.jw, cellSize1 - tempStr.kw, tempMass, tempVel, 4);
        ptrArray_in[tempStr.face][tempStr.ig+1][tempStr.jg+2][tempStr.kg]->UpdateDueToWgt( cellSize1 - tempStr.iw, tempStr.jw + 1,         cellSize1 - tempStr.kw, tempMass, tempVel, 4);
        ptrArray_in[tempStr.face][tempStr.ig+2][tempStr.jg+2][tempStr.kg]->UpdateDueToWgt( tempStr.iw + 1,         tempStr.jw + 1,         cellSize1 - tempStr.kw, tempMass, tempVel, 4);
        
        ptrArray_in[tempStr.face][tempStr.ig+1][tempStr.jg+1][tempStr.kg+1]->UpdateDueToWgt( cellSize1 - tempStr.iw, cellSize1 - tempStr.jw, tempStr.kw + 1,         tempMass, tempVel, 4);
        ptrArray_in[tempStr.face][tempStr.ig+2][tempStr.jg+1][tempStr.kg+1]->UpdateDueToWgt( tempStr.iw + 1,         cellSize1 - tempStr.jw, tempStr.kw + 1,         tempMass, tempVel, 4);
        ptrArray_in[tempStr.face][tempStr.ig+1][tempStr.jg+2][tempStr.kg+1]->UpdateDueToWgt( cellSize1 - tempStr.iw, tempStr.jw + 1,         tempStr.kw + 1,         tempMass, tempVel, 4);
        ptrArray_in[tempStr.face][tempStr.ig+2][tempStr.jg+2][tempStr.kg+1]->UpdateDueToWgt( tempStr.iw + 1,         tempStr.jw + 1,         tempStr.kw + 1,         tempMass, tempVel, 4);
    }


// For O particles in main domain    
    for( list<Particles>::iterator iter= ptrParticlesList_O_in->begin(); iter!=ptrParticlesList_O_in->end(); ++iter)
    {
        // locate the particle
        Particles temp = *iter;
        // get the weighting info of this particle
        struct structg tempStr = temp.InttoStrp1();
        // get the info of mass(weight of each simulation particle)
        double tempMass = temp.WeightMi();
        // get the info of velocity
        Vector3 tempVel = temp.VelParticles();
        // update density and velocity in the grids
        ptrArray_in[tempStr.face][tempStr.ig+1][tempStr.jg+1][tempStr.kg]->UpdateDueToWgt( cellSize1 - tempStr.iw, cellSize1 - tempStr.jw, cellSize1 - tempStr.kw, tempMass, tempVel, 16);
        ptrArray_in[tempStr.face][tempStr.ig+2][tempStr.jg+1][tempStr.kg]->UpdateDueToWgt( tempStr.iw + 1,         cellSize1 - tempStr.jw, cellSize1 - tempStr.kw, tempMass, tempVel, 16);
        ptrArray_in[tempStr.face][tempStr.ig+1][tempStr.jg+2][tempStr.kg]->UpdateDueToWgt( cellSize1 - tempStr.iw, tempStr.jw + 1,         cellSize1 - tempStr.kw, tempMass, tempVel, 16);
        ptrArray_in[tempStr.face][tempStr.ig+2][tempStr.jg+2][tempStr.kg]->UpdateDueToWgt( tempStr.iw + 1,         tempStr.jw + 1,         cellSize1 - tempStr.kw, tempMass, tempVel, 16);
        
        ptrArray_in[tempStr.face][tempStr.ig+1][tempStr.jg+1][tempStr.kg+1]->UpdateDueToWgt( cellSize1 - tempStr.iw, cellSize1 - tempStr.jw, tempStr.kw + 1,         tempMass, tempVel, 16);
        ptrArray_in[tempStr.face][tempStr.ig+2][tempStr.jg+1][tempStr.kg+1]->UpdateDueToWgt( tempStr.iw + 1,         cellSize1 - tempStr.jw, tempStr.kw + 1,         tempMass, tempVel, 16);
        ptrArray_in[tempStr.face][tempStr.ig+1][tempStr.jg+2][tempStr.kg+1]->UpdateDueToWgt( cellSize1 - tempStr.iw, tempStr.jw + 1,         tempStr.kw + 1,         tempMass, tempVel, 16);
        ptrArray_in[tempStr.face][tempStr.ig+2][tempStr.jg+2][tempStr.kg+1]->UpdateDueToWgt( tempStr.iw + 1,         tempStr.jw + 1,         tempStr.kw + 1,         tempMass, tempVel, 16);
    }


// For H particles in temp domain    
    for( list<Particles>::iterator iter= ptrParticlesListTemp_H_in->begin(); iter!=ptrParticlesListTemp_H_in->end(); ++iter)
    {
        // locate the particle
        Particles temp = *iter;
        // get the weighting info of this particle
        struct structg tempStr = temp.InttoStrp1();
        // get the info of mass(weight of each simulation particle)
        double tempMass = temp.WeightMi();

        Vector3 tempVel = temp.VelParticles();
   
        if( tempStr.kg == 0)
        {    
        ptrArray_in[tempStr.face][tempStr.ig+1][tempStr.jg+1][tempStr.kg+1]->UpdateDueToWgt( cellSize1 - tempStr.iw, cellSize1 - tempStr.jw, tempStr.kw + 1,         tempMass, tempVel, 1);
        ptrArray_in[tempStr.face][tempStr.ig+2][tempStr.jg+1][tempStr.kg+1]->UpdateDueToWgt( tempStr.iw + 1,         cellSize1 - tempStr.jw, tempStr.kw + 1,         tempMass, tempVel, 1);
        ptrArray_in[tempStr.face][tempStr.ig+1][tempStr.jg+2][tempStr.kg+1]->UpdateDueToWgt( cellSize1 - tempStr.iw, tempStr.jw + 1,         tempStr.kw + 1,         tempMass, tempVel, 1);
        ptrArray_in[tempStr.face][tempStr.ig+2][tempStr.jg+2][tempStr.kg+1]->UpdateDueToWgt( tempStr.iw + 1,         tempStr.jw + 1,         tempStr.kw + 1,         tempMass, tempVel, 1);
        }
        else if( tempStr.kg == fieldsGridsSize)
        {
        ptrArray_in[tempStr.face][tempStr.ig+1][tempStr.jg+1][tempStr.kg]->UpdateDueToWgt( cellSize1 - tempStr.iw, cellSize1 - tempStr.jw, cellSize1 - tempStr.kw, tempMass, tempVel, 1);
        ptrArray_in[tempStr.face][tempStr.ig+2][tempStr.jg+1][tempStr.kg]->UpdateDueToWgt( tempStr.iw + 1,         cellSize1 - tempStr.jw, cellSize1 - tempStr.kw, tempMass, tempVel, 1);
        ptrArray_in[tempStr.face][tempStr.ig+1][tempStr.jg+2][tempStr.kg]->UpdateDueToWgt( cellSize1 - tempStr.iw, tempStr.jw + 1,         cellSize1 - tempStr.kw, tempMass, tempVel, 1);
        ptrArray_in[tempStr.face][tempStr.ig+2][tempStr.jg+2][tempStr.kg]->UpdateDueToWgt( tempStr.iw + 1,         tempStr.jw + 1,         cellSize1 - tempStr.kw, tempMass, tempVel, 1);
        }

    }

// For He particles in temp domain    
    for( list<Particles>::iterator iter= ptrParticlesListTemp_He_in->begin(); iter!=ptrParticlesListTemp_He_in->end(); ++iter)
    {
        // locate the particle
        Particles temp = *iter;
        // get the weighting info of this particle
        struct structg tempStr = temp.InttoStrp1();
        // get the info of mass(weight of each simulation particle)
        double tempMass = temp.WeightMi();

        Vector3 tempVel = temp.VelParticles();
   
        if( tempStr.kg == 0)
        {    
        ptrArray_in[tempStr.face][tempStr.ig+1][tempStr.jg+1][tempStr.kg+1]->UpdateDueToWgt( cellSize1 - tempStr.iw, cellSize1 - tempStr.jw, tempStr.kw + 1,         tempMass, tempVel, 4);
        ptrArray_in[tempStr.face][tempStr.ig+2][tempStr.jg+1][tempStr.kg+1]->UpdateDueToWgt( tempStr.iw + 1,         cellSize1 - tempStr.jw, tempStr.kw + 1,         tempMass, tempVel, 4);
        ptrArray_in[tempStr.face][tempStr.ig+1][tempStr.jg+2][tempStr.kg+1]->UpdateDueToWgt( cellSize1 - tempStr.iw, tempStr.jw + 1,         tempStr.kw + 1,         tempMass, tempVel, 4);
        ptrArray_in[tempStr.face][tempStr.ig+2][tempStr.jg+2][tempStr.kg+1]->UpdateDueToWgt( tempStr.iw + 1,         tempStr.jw + 1,         tempStr.kw + 1,         tempMass, tempVel, 4);
        }
        else if( tempStr.kg == fieldsGridsSize)
        {
        ptrArray_in[tempStr.face][tempStr.ig+1][tempStr.jg+1][tempStr.kg]->UpdateDueToWgt( cellSize1 - tempStr.iw, cellSize1 - tempStr.jw, cellSize1 - tempStr.kw, tempMass, tempVel, 4);
        ptrArray_in[tempStr.face][tempStr.ig+2][tempStr.jg+1][tempStr.kg]->UpdateDueToWgt( tempStr.iw + 1,         cellSize1 - tempStr.jw, cellSize1 - tempStr.kw, tempMass, tempVel, 4);
        ptrArray_in[tempStr.face][tempStr.ig+1][tempStr.jg+2][tempStr.kg]->UpdateDueToWgt( cellSize1 - tempStr.iw, tempStr.jw + 1,         cellSize1 - tempStr.kw, tempMass, tempVel, 4);
        ptrArray_in[tempStr.face][tempStr.ig+2][tempStr.jg+2][tempStr.kg]->UpdateDueToWgt( tempStr.iw + 1,         tempStr.jw + 1,         cellSize1 - tempStr.kw, tempMass, tempVel, 4);
        }

    }

// For O particles in temp domain    
    for( list<Particles>::iterator iter= ptrParticlesListTemp_O_in->begin(); iter!=ptrParticlesListTemp_O_in->end(); ++iter)
    {
        // locate the particle
        Particles temp = *iter;
        // get the weighting info of this particle
        struct structg tempStr = temp.InttoStrp1();
        // get the info of mass(weight of each simulation particle)
        double tempMass = temp.WeightMi();

        Vector3 tempVel = temp.VelParticles();

        if( tempStr.kg == 0)
        {    
        ptrArray_in[tempStr.face][tempStr.ig+1][tempStr.jg+1][tempStr.kg+1]->UpdateDueToWgt( cellSize1 - tempStr.iw, cellSize1 - tempStr.jw, tempStr.kw + 1,         tempMass, tempVel, 16);
        ptrArray_in[tempStr.face][tempStr.ig+2][tempStr.jg+1][tempStr.kg+1]->UpdateDueToWgt( tempStr.iw + 1,         cellSize1 - tempStr.jw, tempStr.kw + 1,         tempMass, tempVel, 16);
        ptrArray_in[tempStr.face][tempStr.ig+1][tempStr.jg+2][tempStr.kg+1]->UpdateDueToWgt( cellSize1 - tempStr.iw, tempStr.jw + 1,         tempStr.kw + 1,         tempMass, tempVel, 16);
        ptrArray_in[tempStr.face][tempStr.ig+2][tempStr.jg+2][tempStr.kg+1]->UpdateDueToWgt( tempStr.iw + 1,         tempStr.jw + 1,         tempStr.kw + 1,         tempMass, tempVel, 16);
        }
        else if( tempStr.kg == fieldsGridsSize)
        {
        ptrArray_in[tempStr.face][tempStr.ig+1][tempStr.jg+1][tempStr.kg]->UpdateDueToWgt( cellSize1 - tempStr.iw, cellSize1 - tempStr.jw, cellSize1 - tempStr.kw, tempMass, tempVel, 16);
        ptrArray_in[tempStr.face][tempStr.ig+2][tempStr.jg+1][tempStr.kg]->UpdateDueToWgt( tempStr.iw + 1,         cellSize1 - tempStr.jw, cellSize1 - tempStr.kw, tempMass, tempVel, 16);
        ptrArray_in[tempStr.face][tempStr.ig+1][tempStr.jg+2][tempStr.kg]->UpdateDueToWgt( cellSize1 - tempStr.iw, tempStr.jw + 1,         cellSize1 - tempStr.kw, tempMass, tempVel, 16);
        ptrArray_in[tempStr.face][tempStr.ig+2][tempStr.jg+2][tempStr.kg]->UpdateDueToWgt( tempStr.iw + 1,         tempStr.jw + 1,         cellSize1 - tempStr.kw, tempMass, tempVel, 16);
        }
    }

// finish culmulating and average the density and velocity
    if( timeline_in % updateInfoPeriod_in == 0)
    {
        for( int face = 0; face < totalFace; face++)
        {
            for( int i = 1; i < fieldsGridsSize+2; i++)
            {
                for( int j = 1; j < fieldsGridsSize+2; j++)
                {
                    for( int k = 1; k < fieldsGridsSize; k++)
                    {           
                    //check stopsign
                    if( ptrArray_in[face][i][j][k]->StopSign() == 1) continue;
                    // set volume 
                    double volume = ptrVolumeGridArray_in[i-1][j-1][k]; // face of ptrArray is greater than that of ptrVolumeGridArray

//                    std::cout << face << i << j << k << " " ;
//                    std::cout << " volume " << volume << " density " << ptrArray_in[face][i][j][k]->Density() << " ==> " ;
                              
                    ptrArray_in[face][i][j][k]->UpdateDueToWgt(ptrArray_in, volume);
                    // set stopSign
                    ptrArray_in[face][i][j][k]->SetStopSign(1);

//                    std::cout << ptrArray_in[face][i][j][k]->Density() << std::endl;
                    }
                }
            }
        }   

        // reset stopsign
        for( int face = 0; face < totalFace; face++)
        {
            for( int i = 1; i < fieldsGridsSize+2; i++)
            {
                for( int j = 1; j < fieldsGridsSize+2; j++)
                {
                    for( int k = 1; k < fieldsGridsSize; k++)
                    {           
                    // set stopSign
                    ptrArray_in[face][i][j][k]->SetStopSign(0);
                    }
                }
            }
        }  
    }
}




//************************************************************************
//************************************************************************
// FUNCTION
// Printout the gridpoints on the girds as hdf5 format
// Step 1: Generate a new matrix fulling with gridspoints class
// Step 2: Print out it as .h5
//************************************************************************
//************************************************************************
void PrintOutHdf5( GridsPoints***** ptrArray_in, int i_in, int h5FileCheck_in)
{
    using namespace H5;
    char filename[80];
    sprintf( filename, "ArrayOfGrids_%d", i_in);

    H5std_string FILE_NAME( "GridsData.h5");
    H5std_string DATASET_NAME( filename);
    H5std_string DATASET_CONST_NAME( "ArrayOfGrids_const");

    H5std_string MEMBERx( "x");
    H5std_string MEMBERy( "y");
    H5std_string MEMBERz( "z");

    H5std_string MEMBER_pos3( "pos3");
    H5std_string MEMBER_e3( "e3");

    H5std_string MEMBER_b3( "b3");
    H5std_string MEMBER_dB3( "dB3");

    H5std_string MEMBER_ve3( "ve3");
    H5std_string MEMBER_v3( "v3");
    H5std_string MEMBER_vH3( "vH3");
    H5std_string MEMBER_vHe3( "vHe3");
    H5std_string MEMBER_vO3( "vO3");

    H5std_string MEMBER_density_H( "densityH");
    H5std_string MEMBER_density_He( "densityHe");
    H5std_string MEMBER_density_O( "densityO");
    H5std_string MEMBER_density( "density");
    H5std_string MEMBER_te( "temperature");

    const int RANK = 4;

    typedef struct Vector3_h5{
        double v_x;
        double v_y;
        double v_z;
    } Vector3_h5;

    typedef struct GridsPoints_h5{
        Vector3_h5 e3;
        Vector3_h5 dB3;

        Vector3_h5 ve3;
        Vector3_h5 v3;
        Vector3_h5 vH3;
        Vector3_h5 vHe3;
        Vector3_h5 vO3;

        double density;
        double density_H;
        double density_He;
        double density_O;
    } GridsPoints_h5;   // varies with timeline

    typedef struct GridsPoints_const_h5{
        Vector3_h5 pos3;
        Vector3_h5 b3;

        double temperature;
    } GridsPoints_const_h5;     // dont varies with timeline


    // Apply for continus memory 
    GridsPoints_h5* data_mem = new GridsPoints_h5[totalFace* (1+fieldsGridsSize)* (1+fieldsGridsSize)* (1+fieldsGridsSize)];
    GridsPoints_h5**** array_data = new GridsPoints_h5***[totalFace];
    for( int face=0; face<totalFace; face++)
    {
        array_data[face] = new GridsPoints_h5**[1+fieldsGridsSize];
        for( int i = 0; i< 1+fieldsGridsSize; i++)
        {
            array_data[face][i] = new GridsPoints_h5*[1+fieldsGridsSize];
            for( int j=0; j< 1+fieldsGridsSize; j++)
            {
                array_data[face][i][j] = new GridsPoints_h5[1+fieldsGridsSize];
                array_data[face][i][j] = data_mem + face*(1+fieldsGridsSize)*(1+fieldsGridsSize)*(1+fieldsGridsSize)+
                                         i*(1+fieldsGridsSize)*(1+fieldsGridsSize)+
                                         j*(1+fieldsGridsSize);

                for( int k = 0; k < 1 + fieldsGridsSize; k++)
                {
                    array_data[face][i][j][k].e3 
                        = {ptrArray_in[face][i+1][j+1][k]->E3().x(), ptrArray_in[face][i+1][j+1][k]->E3().y(), ptrArray_in[face][i+1][j+1][k]->E3().z()} ;

                    array_data[face][i][j][k].dB3 
                        = {ptrArray_in[face][i+1][j+1][k]->DB3().x(), ptrArray_in[face][i+1][j+1][k]->DB3().y(), ptrArray_in[face][i+1][j+1][k]->DB3().z()} ;
                    
                    array_data[face][i][j][k].ve3 
                        = {ptrArray_in[face][i+1][j+1][k]->Vel_e3().x(), ptrArray_in[face][i+1][j+1][k]->Vel_e3().y(), ptrArray_in[face][i+1][j+1][k]->Vel_e3().z()} ;
                    
                    array_data[face][i][j][k].v3 
                        = {ptrArray_in[face][i+1][j+1][k]->Vel3().x(), ptrArray_in[face][i+1][j+1][k]->Vel3().y(), ptrArray_in[face][i+1][j+1][k]->Vel3().z()} ;
                    array_data[face][i][j][k].vH3 
                        = {ptrArray_in[face][i+1][j+1][k]->VelH3().x(), ptrArray_in[face][i+1][j+1][k]->VelH3().y(), ptrArray_in[face][i+1][j+1][k]->VelH3().z()} ;
                    array_data[face][i][j][k].vHe3 
                        = {ptrArray_in[face][i+1][j+1][k]->VelHe3().x(), ptrArray_in[face][i+1][j+1][k]->VelHe3().y(), ptrArray_in[face][i+1][j+1][k]->VelHe3().z()} ;
                    array_data[face][i][j][k].vO3 
                        = {ptrArray_in[face][i+1][j+1][k]->VelO3().x(), ptrArray_in[face][i+1][j+1][k]->VelO3().y(), ptrArray_in[face][i+1][j+1][k]->VelO3().z()} ;

                    array_data[face][i][j][k].density
                        = ptrArray_in[face][i+1][j+1][k]->Density();
                    array_data[face][i][j][k].density
                        = ptrArray_in[face][i+1][j+1][k]->Density_H();
                    array_data[face][i][j][k].density
                        = ptrArray_in[face][i+1][j+1][k]->Density_He();
                    array_data[face][i][j][k].density
                        = ptrArray_in[face][i+1][j+1][k]->Density_O();
                    //    std::cout << array_data[face][i][j][k].b3.v_x << std::endl;
                }
            }
        }
    }
    
    GridsPoints_const_h5* data_mem_const = new GridsPoints_const_h5[totalFace* (1+fieldsGridsSize)* (1+fieldsGridsSize)* (1+fieldsGridsSize)];
    GridsPoints_const_h5**** array_data_const = new GridsPoints_const_h5***[totalFace];
    for( int face=0; face<totalFace; face++)
    {
        array_data_const[face] = new GridsPoints_const_h5**[1+fieldsGridsSize];
        for( int i = 0; i< 1+fieldsGridsSize; i++)
        {
            array_data_const[face][i] = new GridsPoints_const_h5*[1+fieldsGridsSize];
            for( int j=0; j< 1+fieldsGridsSize; j++)
            {
                array_data_const[face][i][j] = new GridsPoints_const_h5[1+fieldsGridsSize];
                array_data_const[face][i][j] = data_mem_const + face*(1+fieldsGridsSize)*(1+fieldsGridsSize)*(1+fieldsGridsSize)+
                                         i*(1+fieldsGridsSize)*(1+fieldsGridsSize)+
                                         j*(1+fieldsGridsSize);

                for( int k = 0; k < 1 + fieldsGridsSize; k++)
                {       
                    array_data_const[face][i][j][k].pos3 
                        = {ptrArray_in[face][i+1][j+1][k]->Pos3().x(), ptrArray_in[face][i+1][j+1][k]->Pos3().y(), ptrArray_in[face][i+1][j+1][k]->Pos3().z()} ;
                    array_data_const[face][i][j][k].b3 
                        = {ptrArray_in[face][i+1][j+1][k]->B3_base().x(), ptrArray_in[face][i+1][j+1][k]->B3_base().y(), ptrArray_in[face][i+1][j+1][k]->B3_base().z()} ;
                    array_data_const[face][i][j][k].temperature = ptrArray_in[face][i+1][j+1][k]->Temperature();
                }
            }
        }
    }            



    Exception::dontPrint();



    hsize_t dim[] = {totalFace,fieldsGridsSize+1,fieldsGridsSize+1,fieldsGridsSize+1};
    DataSpace space( RANK, dim);


        // non-const group variables
        CompType mtype_vector3( sizeof( Vector3_h5));
        mtype_vector3.insertMember( MEMBERx, HOFFSET(Vector3_h5,v_x), PredType::NATIVE_DOUBLE);
        mtype_vector3.insertMember( MEMBERy, HOFFSET(Vector3_h5,v_y), PredType::NATIVE_DOUBLE);
        mtype_vector3.insertMember( MEMBERz, HOFFSET(Vector3_h5,v_z), PredType::NATIVE_DOUBLE);

        CompType mtype_grids( sizeof( GridsPoints_h5));
        mtype_grids.insertMember( MEMBER_e3, HOFFSET(GridsPoints_h5,e3), mtype_vector3);

        mtype_grids.insertMember( MEMBER_dB3, HOFFSET(GridsPoints_h5,dB3), mtype_vector3);
        
        mtype_grids.insertMember( MEMBER_ve3, HOFFSET(GridsPoints_h5,ve3), mtype_vector3);
        mtype_grids.insertMember( MEMBER_v3, HOFFSET(GridsPoints_h5,v3), mtype_vector3);
        mtype_grids.insertMember( MEMBER_vH3, HOFFSET(GridsPoints_h5,vH3), mtype_vector3);
        mtype_grids.insertMember( MEMBER_vHe3, HOFFSET(GridsPoints_h5,vHe3), mtype_vector3);
        mtype_grids.insertMember( MEMBER_vO3, HOFFSET(GridsPoints_h5,vO3), mtype_vector3);
        
        mtype_grids.insertMember( MEMBER_density, HOFFSET(GridsPoints_h5,density), PredType::NATIVE_DOUBLE);

        mtype_grids.insertMember( MEMBER_density_H, HOFFSET(GridsPoints_h5,density_H), PredType::NATIVE_DOUBLE);
        mtype_grids.insertMember( MEMBER_density_He, HOFFSET(GridsPoints_h5,density_He), PredType::NATIVE_DOUBLE);
        mtype_grids.insertMember( MEMBER_density_O, HOFFSET(GridsPoints_h5,density_O), PredType::NATIVE_DOUBLE);

        // const group variables
        CompType mtype_grids_const( sizeof( GridsPoints_const_h5));
        mtype_grids_const.insertMember( MEMBER_pos3, HOFFSET(GridsPoints_const_h5,pos3), mtype_vector3);
        mtype_grids_const.insertMember( MEMBER_b3, HOFFSET(GridsPoints_const_h5,b3), mtype_vector3);

        mtype_grids_const.insertMember( MEMBER_te, HOFFSET(GridsPoints_const_h5,temperature), PredType::NATIVE_DOUBLE);



    if(h5FileCheck_in == 0)
    {
        H5File* file = new H5File( FILE_NAME, H5F_ACC_TRUNC);
        h5FileCheck = 1;
        DataSet* dataset;
        dataset = new DataSet(file->createDataSet(DATASET_NAME, mtype_grids, space));
        dataset->write( array_data[0][0][0], mtype_grids);

        DataSet* dataset_const;
        dataset_const = new DataSet(file->createDataSet( DATASET_CONST_NAME, mtype_grids_const, space));
        dataset_const->write( array_data_const[0][0][0], mtype_grids_const);

        delete dataset;
        delete dataset_const;
        delete data_mem;
        delete data_mem_const;

        delete file;

    }
    else
    {
        H5File* file = new H5File( FILE_NAME, H5F_ACC_RDWR);
  
        DataSet* dataset;
        dataset = new DataSet(file->createDataSet(DATASET_NAME, mtype_grids, space));
        dataset->write( array_data[0][0][0], mtype_grids);

        delete dataset;
        delete file;
        delete data_mem;
    }
    
}

//************************************************************************
//************************************************************************
// FUNCTION
// Create a array to store volume of cells needed to calculate density 
// at each grids points 
// The volume means the sum volume of adjacent 6 or 8 cells
// Only consider the main domain so that the size is (fieldsGridsSize+1)*(fieldsGridsSize+1)*(fieldsGridsSize+1)
//************************************************************************
//************************************************************************
double*** VolumeGridsField( double*** ptrVolumeCellArray_in)
{
    static double* mem_VolumeGridsArray = new double[(fieldsGridsSize+1)*(fieldsGridsSize+1)*(fieldsGridsSize+1)]; 
    double*** VolumeGridsArray = new double **[fieldsGridsSize+1];

    for( int i =0; i < fieldsGridsSize+1; i++)
    {
        VolumeGridsArray[i] = new double*[fieldsGridsSize+1];
        for( int j =0; j < fieldsGridsSize+1; j++)
        {
            VolumeGridsArray[i][j]= new double [fieldsGridsSize+1];
            VolumeGridsArray[i][j]= mem_VolumeGridsArray + i* (fieldsGridsSize+1)*(fieldsGridsSize+1) + j*(fieldsGridsSize+1);
            for( int k = 0; k < fieldsGridsSize+1; k++)
            {
                if( k == 0 || k == fieldsGridsSize)
                {
                    VolumeGridsArray[i][j][k] = 999.99; // should not be used
                }
                else
                {
                    if( i ==0 && j == 0) { VolumeGridsArray[i][j][k] = 3* ptrVolumeCellArray_in[i+1][j+1][k] + 3* ptrVolumeCellArray_in[i+1][j+1][k+1]; }

                    else if( i ==0 && j ==fieldsGridsSize) { VolumeGridsArray[i][j][k] = 3* ptrVolumeCellArray_in[i+1][j][k] + 3* ptrVolumeCellArray_in[i+1][j][k+1];}
  
                    else if( i ==fieldsGridsSize && j ==0) { VolumeGridsArray[i][j][k] = 3* ptrVolumeCellArray_in[i][j+1][k] + 3* ptrVolumeCellArray_in[i][j+1][k+1];}
     
                    else if( i ==fieldsGridsSize && j ==fieldsGridsSize) { VolumeGridsArray[i][j][k] = 3* ptrVolumeCellArray_in[i][j][k] + 3* ptrVolumeCellArray_in[i][j][k+1];}
              
                    else if( i ==0 && j !=0 && j != fieldsGridsSize) { VolumeGridsArray[i][j][k] = 2* ptrVolumeCellArray_in[i+1][j+1][k] + 2* ptrVolumeCellArray_in[i+1][j][k]
                                                                        + 2* ptrVolumeCellArray_in[i+1][j+1][k+1] + 2* ptrVolumeCellArray_in[i+1][j][k+1];}

                    else if( i ==fieldsGridsSize && j !=0 && j!= fieldsGridsSize) { VolumeGridsArray[i][j][k] = 2* ptrVolumeCellArray_in[i][j+1][k] + 2* ptrVolumeCellArray_in[i][j][k]
                                                                                     + 2* ptrVolumeCellArray_in[i][j+1][k+1] + 2* ptrVolumeCellArray_in[i][j][k+1];}

                    else if( j ==0 && i != 0 && i != fieldsGridsSize) { VolumeGridsArray[i][j][k] = 2* ptrVolumeCellArray_in[i+1][j+1][k] + 2* ptrVolumeCellArray_in[i][j+1][k]
                                                                                                    + 2* ptrVolumeCellArray_in[i+1][j+1][k+1] + 2* ptrVolumeCellArray_in[i][j+1][k];}

                    else if( j ==fieldsGridsSize && i != 0 && i != fieldsGridsSize) { VolumeGridsArray[i][j][k] = 2* ptrVolumeCellArray_in[i+1][j][k] + 2* ptrVolumeCellArray_in[i][j][k]
                                                                                                    + 2* ptrVolumeCellArray_in[i+1][j][k+1] + 2* ptrVolumeCellArray_in[i][j][k+1];}

                    else 
                    {
                        VolumeGridsArray[i][j][k] = ptrVolumeCellArray_in[i+1][j+1][k] + ptrVolumeCellArray_in[i][j+1][k] + ptrVolumeCellArray_in[i+1][j][k] + ptrVolumeCellArray_in[i][j][k] 
                                                  + ptrVolumeCellArray_in[i+1][j+1][k+1] + ptrVolumeCellArray_in[i][j+1][k+1] + ptrVolumeCellArray_in[i+1][j][k+1] + ptrVolumeCellArray_in[i][j][k+1];                        
                    }
                }
            }
        }
    }
    return VolumeGridsArray;
}


//************************************************************************
//************************************************************************
// FUNCTION
// Create a array to store volume of cells needed to calculate density 
// at each cell [fieldsize+2][fieldsize+2][fieldsize]
//************************************************************************
//************************************************************************
double*** VolumeCellsField( GridsPoints***** ptrArray_in)
{
    static double* mem_VolumeCellsArray = new double[(fieldsGridsSize+2) * (fieldsGridsSize+2) * (fieldsGridsSize)];
    double*** VolumeCellsArray = new double**[fieldsGridsSize+2];
    for( int i = 0; i < fieldsGridsSize+2; i++)
    {
        VolumeCellsArray[i] = new double*[fieldsGridsSize+2];
        for( int j = 0; j < fieldsGridsSize+2; j++)
        {         
            VolumeCellsArray[i][j] = new double[fieldsGridsSize];
            VolumeCellsArray[i][j] = mem_VolumeCellsArray + i* (fieldsGridsSize)* (fieldsGridsSize+2) 
                                     + j* (fieldsGridsSize);
            for( int k = 0; k< fieldsGridsSize; k++)
            {
                VolumeCellsArray[i][j][k] = CellVolume( ptrArray_in, 0, i, j, k);
            /*  if( VolumeCellsArray[i][j][k] == 0){
                std::cout << i << " " << j << " " << k << std::endl;
                std::cout << VolumeCellsArray[i][j][k] << std::endl;
                }
            */
            }
        }
    }
    return VolumeCellsArray;
}

//************************************************************************
//************************************************************************
// FUNCTION
// Calculation the temprature
//************************************************************************
//************************************************************************
void Titheridge_Te( GridsPoints***** ptrArray_in) 
//  by Yuzhang Ma AUG 16,2019
{
//inputs
double  SEC = 18*3600; // second of day
double  DOY = 31;      // day of year

double  Kp   = 10;
double  F107 = 100;

//  consts
const double rad2deg = 57.2957763671875;
const double pi = 3.1415927410125732;
const double PLAT=1.37846; //geographic latitude (rad) of Earth's north magnetic pole
const double PLON=5.04557; //geographic longitude (rad) of Earth's north magnetic pole
const double h0 = 400.0e0,     R0 = 1.0627825e0,R02 = 1.1295067e0;
const double por = 0.2857142e0,Re = 6371.2e0,   x1 = 1.0470869e0;
const double a0 = 1.23e3,a1 = 2.2e3,  a2 = -3.29e3,a3 = -2.6e-1,a4 = -6.8e-1;
const double b0 = 4.65e0,b1 = -8.55e0,b2 = 4.14e0, b3 = -2.16e0,b4 = 1.45e0;
const double aa0 = 9.85e2, aa1 = -9.63e2,aa2 = 1.125e3,aa3 = -0.6e0,aa4 = 0.10e0; 
const double bb0 = 7.56e-1,bb1 = -8.8e-1,bb2 = 2.9e-1,bb3 = -2.63e0,bb4 = 1.84e0;

double  GLONR=15.0*(12.0-SEC/3600.0)/rad2deg;//GLONR in rad

        if(GLONR < 0.0) 
		{  
	    GLONR += 2.0*pi; // GLONR = GLONR + 2*pi
		}
		if(GLONR > 2.0*pi) 
		{  
	    GLONR -= 2.0*pi; // GLONR = GLONR - 2*pi
		}

double  MLON12=atan(sin(PLAT)*tan(GLONR-PLON));

        if(MLON12 < 0.0)
		{
	    MLON12+=2.0*pi; // MLON12 = MLON12 + 2*pi
	    }		  

//		  
//  caution: Px_in py_in pz_in should be transformed to SM Coordinate before using. 
//

   
    for( int face = 0; face < totalFace; face++)
    {
        for( int i = 1; i < fieldsGridsSize+2; i++)
        {
            for( int j = 1; j < fieldsGridsSize+2; j++)
            {
                for( int k = 0; k < fieldsGridsSize+1; k++)
                { 
                    
                    //check stopsign
                    if( ptrArray_in[face][i][j][k]->StopSign() == 1) continue;

                    double px_in = ptrArray_in[face][i][j][k]->Pos3().x();
                    double py_in = ptrArray_in[face][i][j][k]->Pos3().y();
                    double pz_in = ptrArray_in[face][i][j][k]->Pos3().z();

    //    std::cout << face << i<<j<< k<< " --> ";
    //    std::cout << px_in << " " << py_in << " " << pz_in << " A--> ";
        
                    double  rr=ptrArray_in[face][i][j][k]->Pos3().norm(); // in unit of m
                    double  Lat= asin( pz_in / rr);	
                    double  phi=atan(py_in/px_in);    
                    double  L= rr/ pow(cos( Lat), 2.0);  
                    double  heq=(L-1.0e0)*Re;	
                    double  SL=1.0e0-x1/L;
                    double  BLATD=asin(sqrt(SL))*rad2deg; // BLATD in degree geomagnetic latitude at 300km       
                    double  BLOND=(phi+MLON12)*rad2deg;// BLOND in degree // NaN for polar points

        // std::cout << phi << " " << BLATD << " " << BLOND << " B==> ";

                    //GLAT calculate
                    double  BLONR=BLOND/rad2deg;
                    double  BLATR=BLATD/rad2deg;
                    double  XM=cos(BLATR)*cos(BLONR);
                    double  ZM=sin(BLATR);
                    double  XG=XM*sin(PLAT)+ZM*cos(PLAT);
                    double  ZG=-XM*cos(PLAT)+ZM*sin(PLAT);
                    double	GLATD=asin(ZG)*rad2deg;// GLATD in degree
                    double  GLATR=GLATD/rad2deg;
                    //  day time T0 and G0  Eq.(19) in Titheridge, JGR1998
                    double  T0d=(a0+a1*SL+a2*SL*SL)/(1.0e0+a3*SL+a4*SL*SL);
                    double  G0d=(b0+b1*SL+b2*SL*SL)/(1.0e0+b3*SL+b4*SL*SL);
    
                    //  night time T0 and G0
                    double  T0n=(aa0+aa1*SL+aa2*SL*SL)/(1.0e0+aa3*SL+aa4*SL*SL);
                    double  G0n=(bb0+bb1*SL+bb2*SL*SL)/(1.0e0+bb3*SL+bb4*SL*SL);

                    //  change in solar activity
                    double  delta_T0=3.4*(F107-120.0);
                    G0d=G0d*T0d/(T0d+delta_T0);
                    T0d=T0d+delta_T0;
                    G0n=G0n*T0n/(T0n+delta_T0);
                    T0n=T0n+delta_T0;
                    //  change in magnetic activity
                    if(BLATD > 46.0e0)		  
                    {
                        double  z=0.135*(BLATD-46.0)/rad2deg;
                        double  Dk=37.0+1.33*(BLATD-46.0)-37.0*cos(z);
                        delta_T0=Dk*(Kp-1.5);
                        G0d=G0d*pow(T0d/(T0d+delta_T0),2.5);
                        T0d=T0d+delta_T0;
                        G0n=G0n*pow(T0n/(T0n+delta_T0),2.5);
                        T0n=T0n+delta_T0;
		            }	  

                    double  G0T0d=G0d/T0d;
                    double  G0T0n=G0n/T0n;

                    double alt=(rr / 1000 / Re -1)*Re;
                    double alg=log(alt/h0);
                    double R2=pow(1.0e0+alt/Re,2.0);
                    double Bh=0.05e0/(2.0e0*L-R0)*(88.0e0+alg*(10.5e0-alg));

       // std::cout << T0d << " " << Bh << " " << G0T0d << " " << G0T0n << " heq " << heq << " alt " << alt << " C==> ";
                    //  mean day values from Eq.(13) in Titheridge, 1998 JGR 
                    double Tday0=T0d*pow(1.0e0+Bh*G0T0d*((heq-h0)/R02-(heq-alt)/R2),por);
                    //  mean night values from Eq.(13) in Titheridge, 1998 JGR 
                    double Tnig0=T0n*pow(1.0e0+Bh*G0T0n*((heq-h0)/R02-(heq-alt)/R2),por);


      //  std::cout << Tnig0 << " " << Tday0 << " D==> ";
                    //  local solar time (SAT) inputs: day of year,GLONR,local time
                    //  https://pvcdrom.pveducation.org/SUNLIGHT/SOLART.HTM
                    double B=300/365*(DOY-81)/rad2deg; //B is degree; trans to rad
                    double EOT=9.87*sin(2*B)-7.53*cos(B)-1.5*sin(B);
                    // get the time difference form GLONR 
                    // https://blog.csdn.net/fct2001140269/article/details/86513925
                    int	currentLon=GLONR*rad2deg;
                    if (currentLon>180)
                    {
	                    currentLon=-1*(360-currentLon);
                    }
                    int timeZone;
                    int shangValue = (int) (currentLon / 15);
                    double yushuValue = abs(currentLon % 15);
                    if (yushuValue <= 7.5) 
		            {
                        timeZone = shangValue;
                    } else 
		            {
                        timeZone = shangValue + (currentLon > 0 ? 1 : -1);
                    }
                    double LSTM=15*timeZone;
                    double TC=4*(LSTM-GLONR*rad2deg)+EOT;
                    double SAT=SEC/3600+TC/60;
                    if (SAT < 0) {SAT+=24.0;}
	                if (SAT > 24.0) {SAT-=24.0;}
                    //DEC -- solar declination angle in radians
                    //solar declination angle inputs: day of year, output is in dgree, trans to rad)
                    //https://www.sciencedirect.com/topics/engineering/solar-declination
                    double DEC=23.45*sin(360/365*(284+DOY))/rad2deg;
                    double S_LONG=-1.0e0*tan(GLATR)*tan(DEC);
                    double s_r;
                    if(S_LONG >= -1.0e0 && S_LONG <= 1.0e0)
	                {
		                s_r=12.0e0-3.82*acos(S_LONG);
		            }
	                else
		            { 
                        if(S_LONG<=-1.0e0){s_r=0.0e0;}
                        else{s_r=12.0e0;}
	                }


                    //duration of the day-night transition 
                    double Delta_tr=1.2+0.5*SL;
                    double Delta_ts=Delta_tr+0.9;
                    //time of the center of the transition in hours after ground sunrise or sunset
                    double  t_s=0.15*s_r;
                    if(t_s >= 1.0e0){t_s=1.0;}
                    double  t_r=-0.5*t_s;
                    double  D;      
                    if(SAT<=12.0e0)
		            {D=(s_r-0.5*t_s-SAT)/Delta_tr;}
                    else
		            {D=(SAT-t_s+s_r-24.0e0)/Delta_ts;}
                    phi=1.0/(1.0+exp(3.2*D));
                    double  y=s_r-12.5;
                    if(y<=-4.0){y=-4.0;}
                    double  TDvar=0.97+0.22*pow((SAT-12.5)/y,2);
                    double  Tday;        
                    if(TDvar<=1.2)
		            {Tday=TDvar*Tday0;}
                    else
		            {Tday=1.2*Tday0;}
        
                    if(SAT<=12.0e0)
		            {y=SAT;}
                    else
                    {y=SAT-24.0;}
                    double  TNvar=0.83+0.15*exp(-0.2*y*(1.4-0.8*SL));
                    double  Tnig;       
                    if(TNvar<=1.2)
		            {Tnig=TNvar*Tnig0;}
                    else
                    {Tnig=1.2*Tnig0;}
         
                    double  Te=Tnig+(Tday-Tnig)*phi;// should be the avrage value of 2 Te points,it could be done out of the function

       // std::cout << Te << std::endl;

                    ptrArray_in[face][i][j][k]->SetTemperature(Te);
                    
                    // set stopSign
                    ptrArray_in[face][i][j][k]->SetStopSign(1);
                }
            }
        }
    }

    for( int face = 0; face < totalFace; face++)
    {
        for( int i = 1; i < fieldsGridsSize+2; i++)
        {
            for( int j = 1; j < fieldsGridsSize+2; j++)
            {
                for( int k = 0; k < fieldsGridsSize+1; k++)
                {           
                    // set stopSign  
                    ptrArray_in[face][i][j][k]->SetStopSign(0);
                }
            }
        }
    }    

}

//************************************************************************
//************************************************************************
// Function
// initial the bot boundary for the  velocity of magnetic field line
//************************************************************************
//************************************************************************
void SetBotBoundary( GridsPoints***** ptrArray_in)
{
    double PI = 3.1415926535897;
/*  double rho_max; 
    double rho_min;
    double ratioH = 0.05;
    double ratioHe = 0.05;
    double ratioO = 0.9;
*/
    // Set the dawn and equatorial point as zero point
    
    std::cout << " Set Bot boundary " << std::endl;
    for( int face = 0; face < totalFace; face++)
    {
        for( int i = 1; i < fieldsGridsSize+2; i++)
        {
            for( int j = 1; j < fieldsGridsSize+2; j++)
            {
                int k = 0;
    if( ptrArray_in[face][i][j][k]->StopSign() == 1) continue;

    double x = ptrArray_in[face][i][j][k]->Pos3().x();
    double y = ptrArray_in[face][i][j][k]->Pos3().y();
    double z = ptrArray_in[face][i][j][k]->Pos3().z();


    double longtitude;
    double latitude;

    double A = 0.5 * ( rho_max - rho_min);
    double A_average = 0.5 * ( rho_max + rho_min);
    double rho;

    if( x == 0 && y == 0)
    { rho = A_average;}
    else if( x == 0 && y > 0)
    { longtitude = PI / 2.0;}
    else if( x == 0 && y < 0)
    { longtitude = PI / 2.0 * 3.0;}
    else if( x!= 0)
    {

    longtitude = atan( y / x);
    
    if( x<0) { longtitude = longtitude + PI;}
    
    }
    
    latitude = PI / 2.0 - acos( z / sqrt( x*x + y*y + z*z));   

    rho = ( A - 2.0 * A / PI * abs( latitude)) * sin( longtitude) + A_average;

    ptrArray_in[face][i][j][k]->Density_H( rho * ratioH);
    ptrArray_in[face][i][j][k]->Density_He( rho * ratioHe);
    ptrArray_in[face][i][j][k]->Density_O( rho * ratioO);

    ptrArray_in[face][i][j][k]->SetStopSign(1);
            }
        }
    }

    for( int face = 0; face < totalFace; face++)
    {
        for( int i = 1; i < fieldsGridsSize+2; i++)
        {
            for( int j = 1; j < fieldsGridsSize+2; j++)
            {
                int k = 0;

                // set stopSign  
                ptrArray_in[face][i][j][k]->SetStopSign(0);
            }
        }
    } 

}


//************************************************************************
//************************************************************************
// Function
// initial the top boundary for the  velocity of magnetic field line
//************************************************************************
//************************************************************************
void SetTopBoundary( GridsPoints***** ptrArray_in)
{
    double PI = 3.1415926535897;
    // input two const
    double r0 = radius * cos( r0_latitude);
    double c0 = radius * cos( c0_latitude);
    double t0 = t0_convection;

    std::cout << " Set Top boundary " << std::endl;
    //non-circle
    for( int face = 0; face < totalFace; face++)
    {
        for( int i = 1; i < fieldsGridsSize+2; i++)
        {
            for( int j = 1; j < fieldsGridsSize+2; j++)
            {
                int k = fieldsGridsSize;
                
    if( ptrArray_in[face][i][j][k]->StopSign() == 1) continue;

    // point at the top shell
    double x = ptrArray_in[face][i][j][k]->Pos3().x();
    double y = ptrArray_in[face][i][j][k]->Pos3().y();
    double z = ptrArray_in[face][i][j][k]->Pos3().z();
    
    double r_top = sqrt( x * x + y * y + z * z);
    double theta_top = acos( z / r_top); 
    double phi_top;
    if( x != 0.0)
    {
        phi_top = atan( y/x );
        if( x < 0) 
        {
            phi_top = phi_top + PI;
        }
    }
    else if ( x == 0.0 && y > 0.0)
    {
        phi_top = PI / 2.0;
    }
    else if ( x == 0.0 && y < 0.0)
    {
        phi_top = PI / 2.0 * ( - 1.0);
    }
    else if ( x ==0.0 && y == 0.0)
    {
        phi_top = 0.0; // special points
    }

   
    // Step 1
    // find related point on the earth shell in the north
    double r_earth;
    double theta_earth = asin( sin( theta_top) * sqrt( r_earth / r_top));
    double phi_earth = phi_top;
    double x_earth = r_earth * sin( theta_earth) * cos(phi_earth);
    double y_earth = r_earth * sin( theta_earth) * sin(phi_earth);
    double z_earth = r_earth * cos( theta_earth);

    if( x==0.0 && y == 0.0)
    {
        x_earth = 0.0;
        y_earth = 0.0;
    }

    // Step 2
    // find the velocity of the point ( |x_earth|, y_earth) on the x-y plane
    double xx = x_earth;
    double yy = y_earth;
    double vx_earth, vy_earth;
    double L;

    double x_prime, y_prime;    // used for region 2
    if( x_earth < 0.0) 
    {
        xx = -1.0 * x_earth;
    }

    if( yy <= r0 - r0 *xx / c0 && yy >= -1.0 * r0 + xx * r0 / c0)       // region 1
    {
        L = 2.0 * r0 * ( 1 - xx / c0);
        vy_earth = -1.0 * PI * L / t0 * sqrt( 0.25 - yy*yy / L / L);
        vx_earth = 0.0;
    }
    else if( yy > r0 - r0 *xx / c0 && yy < -1.0 * r0 + xx * r0 / c0 && xx *xx + yy* yy <= r0) // region 2
    {
        x_prime = -1.0 * ( xx - r0*r0/c0 + sqrt( (xx-r0*r0/c0)*(xx-r0*r0/c0) - (r0*r0/c0/c0 -1.0)*(r0*r0-xx*xx-yy*yy))) / (r0*r0/c0/c0 - 1.0);
        y_prime = r0 * ( 1.0- x_prime / c0);
        
        L = 2.0 * y_prime;  // y direction
        vy_earth = PI * L / t0 * sqrt( 0.25 - yy * yy / L / L);
        L = y_prime; // x direction
        vx_earth = PI * L / t0 * sqrt( (xx - x_prime) / L * ( 1.0 - ( xx - x_prime) / L));
        if( yy < 0.0)
        {
            vx_earth = -1.0 * vx_earth;
        }
    }
    else // other places
    {
        vx_earth = 0.0;
        vy_earth = 0.0;
    }

    // Step 3
    // find the realted velocity on the earth ( x_earth, y_earth, z_earth) or ( r_earth, theta_earth, phi_earth)
    // the velocity on the x and y direction is known as ( vx_earth, vy_earth)
    double vtheta_earth = (vx_earth * cos( phi_earth) + vy_earth * sin( phi_earth)) / cos( theta_earth);
    double vphi_earth = vy_earth * cos( phi_earth) - vx_earth * sin( phi_earth);

    // Step 4
    // find the related velocity on the arbitrary shell as we want using the equation A21 and A22 of 
    // Rasmussen et.al 1992
    double sinchi_top = 2.0 * cos( theta_top) / sqrt( 1.0+ 3.0* cos( theta_top)* cos( theta_top));
    double coschi_top = sin( theta_top) / sqrt( 1.0 + 3.0* cos( theta_top)* cos( theta_top));

    double vr_top = r_top * coschi_top * coschi_top * vtheta_earth / r_earth * 2.0 * cos(theta_earth) / sin( theta_earth);
    double vtheta_top = r_top * sinchi_top * coschi_top * vtheta_earth / r_earth * 2.0 * cos(theta_earth) / sin( theta_earth);
    double vphi_top = r_top * sin( theta_top) * vphi_earth;

    double vx_top = vr_top * sin( theta_top) * cos( phi_top) + 
                    vtheta_top * cos(theta_top) * cos( phi_top) -
                    vphi_top * sin(phi_top);
    double vy_top = vr_top * sin( theta_top) * sin(phi_top) +
                    vtheta_top * cos(theta_top) * sin(phi_top) +
                    vphi_top * cos(phi_top);
    double vz_top = vr_top * cos( theta_top) -
                    vtheta_top* sin(theta_top);

    Vector3 temp = Vector3( vx_top, vy_top, vz_top);
    Vector3 original_vel = ptrArray_in[face][i][j][k]->Vel3();
    ptrArray_in[face][i][j][k]->SetVel_topBoundary( original_vel.PlusProduct( temp));

    Vector3 temp_gradPe = Vector3( 0.0, 0.0, 0.0);
    ptrArray_in[face][i][j][k]->updateE( temp_gradPe);

    double N_H = N0_H * exp(-1.0 * (ptrArray_in[face][i][j][k]->Pos3().norm() - radius) / (ikT / mi0_H / gravity));
    ptrArray_in[face][i][j][k]->Density_H( N_H * mi0_H);
         
    double N_He = N0_He * exp(-1.0 * (ptrArray_in[face][i][j][k]->Pos3().norm() - radius) / (ikT / mi0_He / gravity));
    ptrArray_in[face][i][j][k]->Density_H( N_He * mi0_He);

    double N_O = N0_O * exp(-1.0 * (ptrArray_in[face][i][j][k]->Pos3().norm() - radius) / (ikT / mi0_O / gravity));
    ptrArray_in[face][i][j][k]->Density_H( N_O * mi0_O);   

    // set stopSign
    ptrArray_in[face][i][j][k]->SetStopSign(1);
        
            }
        }
    }

    for( int face = 0; face < totalFace; face++)
    {
        for( int i = 1; i < fieldsGridsSize+2; i++)
        {
            for( int j = 1; j < fieldsGridsSize+2; j++)
            {
                int k = fieldsGridsSize;

                // set stopSign  
                ptrArray_in[face][i][j][k]->SetStopSign(0);
            }
        }
    } 

}

//************************************************************************
//************************************************************************
// FUNCTION
// Calculate the gradient of normal of B
//************************************************************************
//************************************************************************
void GradBNorm( GridsPoints***** ptrArray_in)
{

    for( int face = 0; face < totalFace; face++)
    {
        for( int i = 1; i < fieldsGridsSize+2; i++)
        {
            for( int j = 1; j < fieldsGridsSize+2; j++)
            {
                for( int k = 0; k < fieldsGridsSize+1; k++)
                { 

                    //check stopsign
                    if( ptrArray_in[face][i][j][k]->StopSign() == 1) continue;
                    
                    ptrArray_in[face][i][j][k]->XYZtoGradBNorm();

                    // set stopSign
                    ptrArray_in[face][i][j][k]->SetStopSign(1);
                }
            }
        }
    }

    for( int face = 0; face < totalFace; face++)
    {
        for( int i = 1; i < fieldsGridsSize+2; i++)
        {
            for( int j = 1; j < fieldsGridsSize+2; j++)
            {
                for( int k = 1; k < fieldsGridsSize; k++)
                {           
                    // set stopSign
                    ptrArray_in[face][i][j][k]->SetStopSign(0);
                }
            }
        }
    }    
}



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
    SetTopBoundary( ptrArray);
    SetBotBoundary( ptrArray);

    // Prerun 1.2 // Create Cell centered field array for nesseary calculation for one face of six
    Vector3*** ptrVectorCellArray = VectorCellField();  

    // Prerun 1.3 // Create grids field array of volum for one face of six
    cout << " Create array of Volume at cells and at grids" << endl;
    double*** ptrVolumeCellArray = VolumeCellsField( ptrArray);
    
    cout << " Create array of Volume at cells and at grids" << endl;
    double*** ptrVolumeGridArray = VolumeGridsField( ptrVolumeCellArray);
    
    // Prerun 1.4 // Create particles list, initialize the velocity and position of each particles
    cout << " Create particles list of main domain" << endl;
    list<Particles>* ptrParticlesList_H = ParticlesLists( ptrArray, ptrVolumeCellArray, mi0_H, N0_H);

    list<Particles>* ptrParticlesList_He = ParticlesLists( ptrArray, ptrVolumeCellArray, mi0_He, N0_He);
    
    list<Particles>* ptrParticlesList_O = ParticlesLists( ptrArray, ptrVolumeCellArray, mi0_O, N0_O);
    // Run 2.0
    
    cout << " Start" << endl;
    for( int timeline = 1; timeline <= timeLineLimit; timeline++)   // timeline start with 1
    {
    std::cout << "timeline" << timeline << std::endl;
#pragma omp parallel
{
    #pragma omp sections
    {
        #pragma omp section
        {
    // Run 2.1 // Particles in main domain
        for( list<Particles>::iterator iteratorM = ptrParticlesList_H->begin(); iteratorM != ptrParticlesList_H->end(); ++iteratorM)
        {
//std::cout << " check 0" << std::endl;
            Particles temp = *iteratorM;
            struct structg tempStruct = temp.InttoStrp1();
            
            int check; // check whether in the main domain or not, "0" means in "1" means out
            //test function//  double xxx= temp.VelParticles().x();    //cout << xxx << " ";
            // update velocity // update position
            check = temp.BorisMethod( &tempStruct, ptrArray, mi0_H);
            // check if still in the main domain
            if( check == 1) // out of the domain
            {
                iteratorM = ptrParticlesList_H->erase( iteratorM);
                
            }
        }
        cout << "Particles H " << ptrParticlesList_H->size() << endl;
        }
        #pragma omp section
        {
        for( list<Particles>::iterator iteratorM = ptrParticlesList_He->begin(); iteratorM != ptrParticlesList_He->end(); ++iteratorM)
        {
//std::cout << " check 0" << std::endl;
            Particles temp = *iteratorM;
            struct structg tempStruct = temp.InttoStrp1();
            
            int check; // check whether in the main domain or not, "0" means in "1" means out
            //test function//  double xxx= temp.VelParticles().x();    //cout << xxx << " ";
            // update velocity // update position
            check = temp.BorisMethod( &tempStruct, ptrArray, mi0_He);
            // check if still in the main domain
            if( check == 1) // out of the domain
            {
                iteratorM = ptrParticlesList_He->erase( iteratorM);
            }
        }
        cout << "Particles He " << ptrParticlesList_He->size() << endl;
        }
        #pragma omp section
        {
        for( list<Particles>::iterator iteratorM = ptrParticlesList_O->begin(); iteratorM != ptrParticlesList_O->end(); ++iteratorM)
        {
//std::cout << " check 0" << std::endl;
            Particles temp = *iteratorM;
            struct structg tempStruct = temp.InttoStrp1();
            
            int check; // check whether in the main domain or not, "0" means in "1" means out
            //test function//  double xxx= temp.VelParticles().x();    //cout << xxx << " ";
            // update velocity // update position
            check = temp.BorisMethod( &tempStruct, ptrArray, mi0_O);
            // check if still in the main domain
            if( check == 1) // out of the domain
            {
                iteratorM = ptrParticlesList_O->erase( iteratorM);
            }
        }
        cout << "Particles O " << ptrParticlesList_O->size() << endl;
        }     
    }   
    #pragma omp barrier
}
        // Run 2.2 // Create temp particle lists    
        list<Particles>* ptrParticlesListTemp_H = ParticlesListsTemp( ptrArray, ptrVolumeCellArray, mi0_H);
        list<Particles>* ptrParticlesListTemp_He = ParticlesListsTemp( ptrArray, ptrVolumeCellArray, mi0_He);
        list<Particles>* ptrParticlesListTemp_O = ParticlesListsTemp( ptrArray, ptrVolumeCellArray, mi0_O);

    #pragma omp parallel
    {
    #pragma omp sections
    {
    #pragma omp section
    {
        // Run 2.3 // Particles in temp domain
        for( list<Particles>::iterator iterator = ptrParticlesListTemp_H->begin(); iterator != ptrParticlesListTemp_H->end(); ++iterator)
        {
            Particles temp = *iterator;
            struct structg tempStruct = temp.InttoStrp1();
            int check; // check whether in the main domain or not, "0" means in "1" means out
            // update velocity  // update position
            check = temp.BorisMethod( &tempStruct, ptrArray, mi0_H);
            // check if still in the main domain
            if( check == 0) // in the domain
            {   
                ptrParticlesList_H->push_back( temp);
                iterator = ptrParticlesListTemp_H->erase( iterator);
            }
        }
    }
    #pragma omp section
    {  
   
        for( list<Particles>::iterator iterator = ptrParticlesListTemp_He->begin(); iterator != ptrParticlesListTemp_He->end(); ++iterator)
        {
            Particles temp = *iterator;
            struct structg tempStruct = temp.InttoStrp1();
            int check; // check whether in the main domain or not, "0" means in "1" means out
            // update velocity  // update position
            check = temp.BorisMethod( &tempStruct, ptrArray, mi0_He);
            // check if still in the main domain
            if( check == 0) // in the domain
            {   
                ptrParticlesList_He->push_back( temp);
                iterator = ptrParticlesListTemp_He->erase( iterator);
            }
        }
    }
    #pragma omp section
    {    
    
        for( list<Particles>::iterator iterator = ptrParticlesListTemp_O->begin(); iterator != ptrParticlesListTemp_O->end(); ++iterator)
        {
            Particles temp = *iterator;
            struct structg tempStruct = temp.InttoStrp1();
            int check; // check whether in the main domain or not, "0" means in "1" means out
            // update velocity  // update position
            check = temp.BorisMethod( &tempStruct, ptrArray, mi0_O);
            // check if still in the main domain
            if( check == 0) // in the domain
            {   
                ptrParticlesList_O->push_back( temp);
                iterator = ptrParticlesListTemp_O->erase( iterator);
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
    std::cout << " Update gridspoints info" << std::endl;
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

            }
        }
    }
    std::cout << " PrintOut " << std::endl;
    // Postrun 3.0 // Printout info in grids points
    if( timeline % printTimePeriod ==0)
    {
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

