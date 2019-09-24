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
list<Particles>* ParticlesLists( GridsPoints***** ptrArray_in, double*** ptrVolumeCellArray_in)
{
   //list<Particles> listP;
//    auto listsPtr = make_shared<list<Particles>>();
    list<Particles>* listsPtr = new list<Particles>;
//   shared_ptr<Particles> p1 = make_shared<Particles> (); // test
    double scaleHeight = ikT / mi0 / gravity;

    for( int i = 0; i < totalFace; i++)
    {
        for( int j = 1; j <= fieldsGridsSize; j++)
        {
            for( int k = 1; k <= fieldsGridsSize; k++)
            {
                for( int s = 1; s <= fieldsGridsSize-2; s++)
                {
                    // number of real particles ( notice the unit is number density)
                    double N = N0_i * exp(-1 * (ptrArray_in[i][j][k][s]->Pos3().norm() - radius) / scaleHeight);
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
                    

                    Vector3 vVel = Vector3( MaxwellDisV( ikT, ptrArray_in[i][j][k][s]->Vel3().x()),
                                            MaxwellDisV( ikT, ptrArray_in[i][j][k][s]->Vel3().y()),
                                            MaxwellDisV( ikT, ptrArray_in[i][j][k][s]->Vel3().z()));
                    // put the particles at the end of list
                    Particles tempP= Particles(intPos, vVel, mi_simu);
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
// x = "B" for bot temp
// x = "T" for top temp
//************************************************************************
//************************************************************************
list<Particles>* ParticlesListsTemp( GridsPoints***** ptrArray_in, const char* x)
{
    list<Particles>* listsPtrTemp = new list<Particles>;
    double scaleHeight = ikT / mi0 / gravity;
    double N, mi_simu;
    int s;
    if( x = "B")
    {
        for( int i = 0; i < totalFace; i++)
        {
            for( int j = 1; j <= fieldsGridsSize; j++)
            {
                for( int k =1; k<= fieldsGridsSize; k++)
                {
                    // number of real particles
                    if( x = "B")
                    { s = 0;}
                    else
                    { s = 3+fieldsGridsSize;}
                    
                    N = N0_i * exp(-1 * (ptrArray_in[i][j][k][s]->Pos3().norm() - radius) / scaleHeight);
                    
                    // mass of each simulation particle
                    mi_simu = N / iniParticleNumberPerCell *mi0;
                    for ( int t = 1; t <= iniParticleNumberPerCell; t++)
                    {
                    // calculate random position
                    //test 
                    Vector3 temp1 = ptrArray_in[i][j][k][s]->Pos3();
                    Vector3 temp2 = ptrArray_in[i][j][k][s+1]->Pos3();
                    Vector3 temp = UniformDisVector3( ptrArray_in[i][j][k][s]->Pos3(),ptrArray_in[i][j][k][s+1]->Pos3());

                    Vector3 vPos = ptrArray_in[i][j][k][s]->Pos3().PlusProduct(UniformDisVector3( ptrArray_in[i][j][k][s]->Pos3(),ptrArray_in[i][j][k][s+1]->Pos3()));
                    vPos = vPos.PlusProduct( UniformDisVector3( ptrArray_in[i][j][k][s]->Pos3(),ptrArray_in[i][j][k+1][s]->Pos3()));
                    vPos = vPos.PlusProduct( UniformDisVector3( ptrArray_in[i][j][k][s]->Pos3(),ptrArray_in[i][j+1][k][s]->Pos3()));

                    uint_64 intPos = vPos.Uint_64_Trans();
                    // calculate random velocity
                    Vector3 vVel = Vector3( MaxwellDisV( ikT, ptrArray_in[i][j][k][s]->Vel3().x()),
                                            MaxwellDisV( ikT, ptrArray_in[i][j][k][s]->Vel3().y()),
                                            MaxwellDisV( ikT, ptrArray_in[i][j][k][s]->Vel3().z()));
                    // put the particles at the end of list
                    Particles tempP= Particles(intPos, vVel, mi_simu);
                    listsPtrTemp->push_back( tempP);           
                    }
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
                                                        0.0, 0.0, 0);
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
                                                        0.0, 0.0, 0);
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
                                                        0.0, 0.0, 0);
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
                                                        0.0, 0.0, 0);
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
                                                        0.0, 0.0, 0);
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
                                                        0.0, 0.0, 0);
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
                                                        0.0, 0.0, 0);
            ptrArray[face][0][fieldsGridsSize+2][k] = new GridsPoints( 0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0);
            ptrArray[face][fieldsGridsSize+2][0][k] = new GridsPoints( 0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0);
            ptrArray[face][fieldsGridsSize+2][fieldsGridsSize+2][k] = new GridsPoints( 0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0);
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
Vector3**** ValueCurlField( Vector3**** curlArray_in, GridsPoints***** ptrArray_in, int face_in, char field_in)
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
                            (i==fieldsGridsSize+1&&j==0))  
                            continue;

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
                double volumetemp = CellVolume(ptrArray_in, face_in, i, j, k);
                temp = temp.ScaleProduct(1/volumetemp);
                curlArray_in[i][j][k]->SetVector3(temp); 
            }
        }
    }
    return curlArray_in;
}

//************************************************************************
//************************************************************************
// FUNCTION 
// Value gradient field of Pe.
//************************************************************************
//************************************************************************
Vector3*** ValueGradientPe(Vector3*** gradientArray_in, double*** ptrVolumeCellArray_in, GridsPoints***** ptrArray_in, int face_in)
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

                // for each cell, calculate sum of n(face vector) * densities and devided by
                // Volume to get the gradient at the center of cell  
                Vector3 temp = AreaVectorL( ptrArray_in, face_in, i, j, k).ScaleProduct(
                               FaceDensityL(ptrArray_in, face_in, i, j, k));
                        temp = temp.PlusProduct(
                               AreaVectorR( ptrArray_in, face_in, i, j, k).ScaleProduct(
                               FaceDensityR(ptrArray_in, face_in, i, j, k)));
                        temp = temp.PlusProduct(
                               AreaVectorT( ptrArray_in, face_in, i, j, k).ScaleProduct(
                               FaceDensityT(ptrArray_in, face_in, i, j, k)));
                        temp = temp.PlusProduct(
                               AreaVectorBot( ptrArray_in, face_in, i, j, k).ScaleProduct(
                               FaceDensityBot(ptrArray_in, face_in, i, j, k)));
                        temp = temp.PlusProduct(
                               AreaVectorF( ptrArray_in, face_in, i, j, k).ScaleProduct(
                               FaceDensityF(ptrArray_in, face_in, i, j, k)));
                        temp = temp.PlusProduct(
                               AreaVectorBack( ptrArray_in, face_in, i, j, k).ScaleProduct(
                               FaceDensityBack(ptrArray_in, face_in, i, j, k)));
                double volumetemp = ptrVolumeCellArray_in[i][j][k];
                temp = temp.ScaleProduct(1/volumetemp);
   //             std::cout << volumetemp << " " ;
     //           std::cout << " --> " << temp.x() << " " << temp.y() << " " << temp.z() << std::endl;
                gradientArray_in[i][j][k].SetVector3(temp); 
            }
        }
    }
    return gradientArray_in;
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
//************************************************************************
//************************************************************************
void updateE_nocurrent( Vector3*** gradPe_in, GridsPoints***** ptrArray_in, int face_in)
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
//                std::cout << tempGradPe.x() << " " << tempGradPe.y() << " " << tempGradPe.z() << std::endl;
                // update E
                ptrArray_in[face_in][i][j][k]->updateE(ptrArray_in[face_in][i][j][k]->Vel3(), tempGradPe);
            
//                std::cout << ptrArray_in[face_in][i][j][k]->E3().x() << " " << ptrArray_in[face_in][i][j][k]->E3().y() 
//                          << " " << ptrArray_in[face_in][i][j][k]->E3().z() << std::endl;
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
void updateE( Vector3**** curlB_in, Vector3**** gradPe_in, GridsPoints***** ptrArray_in, int face_in)
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
                ptrArray_in[face_in][i][j][k]->updateE(Ve, tempGradPe);
            }
        }
    }        
}


//************************************************************************
//************************************************************************
// FUNCTION
// Each calculation are on the girds.
// Step 1: Setup curl E for each gridpoint
// Step 2: Update B at each gridpoints
//************************************************************************
//************************************************************************
void updateB( Vector3**** curlE_in, GridsPoints***** ptrArray_in, int face_in)
{
    Vector3 tempCurlE;
    for( int i = 1; i < fieldsGridsSize+2; i++)
    {
        for( int j = 1; j < fieldsGridsSize+2; j++)
        {
            for( int k = 1; k < fieldsGridsSize; k++)
            {
                // curl E at gridspoint
                if(i==1&&j==1) 
                {
                tempCurlE = curlE_in[1][0][k-1]->PlusProduct(
                                    *curlE_in[1][1][k-1]);
                tempCurlE = tempCurlE.PlusProduct(           
                                    *curlE_in[0][1][k-1]);
                tempCurlE = tempCurlE.PlusProduct(
                                    *curlE_in[1][0][k]);
                tempCurlE = tempCurlE.PlusProduct(
                                    *curlE_in[1][1][k]);
                tempCurlE = tempCurlE.PlusProduct(
                                    *curlE_in[0][1][k]).ScaleProduct(1/6);           
                }
                else if(i==1&&j==fieldsGridsSize+1) 
                {
                tempCurlE = curlE_in[1][fieldsGridsSize+1][k-1]->PlusProduct(
                                    *curlE_in[1][fieldsGridsSize][k-1]);
                tempCurlE = tempCurlE.PlusProduct(           
                                    *curlE_in[0][fieldsGridsSize][k-1]);
                tempCurlE = tempCurlE.PlusProduct(
                                    *curlE_in[1][fieldsGridsSize+1][k]);
                tempCurlE = tempCurlE.PlusProduct(
                                    *curlE_in[1][fieldsGridsSize][k]);
                tempCurlE = tempCurlE.PlusProduct(
                                    *curlE_in[0][fieldsGridsSize][k]).ScaleProduct(1/6);
                }
                else if(i==fieldsGridsSize+1&&j==fieldsGridsSize+1)
                {
                tempCurlE = curlE_in[fieldsGridsSize][fieldsGridsSize+1][k-1]->PlusProduct(
                                    *curlE_in[fieldsGridsSize][fieldsGridsSize][k-1]);
                tempCurlE = tempCurlE.PlusProduct(           
                                    *curlE_in[fieldsGridsSize+1][fieldsGridsSize][k-1]);
                tempCurlE = tempCurlE.PlusProduct(
                                    *curlE_in[fieldsGridsSize][fieldsGridsSize+1][k]);
                tempCurlE = tempCurlE.PlusProduct(
                                    *curlE_in[fieldsGridsSize][fieldsGridsSize][k]);
                tempCurlE = tempCurlE.PlusProduct(
                                    *curlE_in[fieldsGridsSize+1][fieldsGridsSize][k]).ScaleProduct(1/6);
                }
                else if(i==fieldsGridsSize+1&&j==1) 
                {
                tempCurlE = curlE_in[fieldsGridsSize][0][k-1]->PlusProduct(
                                    *curlE_in[fieldsGridsSize][1][k-1]);
                tempCurlE = tempCurlE.PlusProduct(           
                                    *curlE_in[fieldsGridsSize+1][1][k-1]);
                tempCurlE = tempCurlE.PlusProduct(
                                    *curlE_in[fieldsGridsSize][0][k]);
                tempCurlE = tempCurlE.PlusProduct(
                                    *curlE_in[fieldsGridsSize][1][k]);
                tempCurlE = tempCurlE.PlusProduct(
                                    *curlE_in[fieldsGridsSize+1][1][k]).ScaleProduct(1/6);          
                } else
                {
                tempCurlE = curlE_in[i-1][j-1][k-1]->PlusProduct(
                                    *curlE_in[i][j-1][k-1]);
                tempCurlE = tempCurlE.PlusProduct(           
                                    *curlE_in[i-1][j][k-1]);
                tempCurlE = tempCurlE.PlusProduct(
                                    *curlE_in[i][j][k-1]);
                tempCurlE = tempCurlE.PlusProduct(
                                    *curlE_in[i-1][j-1][k]);
                tempCurlE = tempCurlE.PlusProduct(
                                    *curlE_in[i][j-1][k]);
                tempCurlE = tempCurlE.PlusProduct(
                                    *curlE_in[i-1][j][k]);    
                tempCurlE = tempCurlE.PlusProduct(
                                    *curlE_in[i][j][k]).ScaleProduct(1/8);
                }
                // update B
                ptrArray_in[face_in][i][j][k]->updateB(tempCurlE);
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
                      list<Particles>* ptrParticlesList_in, 
                      list<Particles>* ptrParticlesListBottom_in,
                      list<Particles>* ptrParticlesListTop_in,
                      double*** ptrVolumeGridArray_in,
                      int timeline_in, int updateInfoPeriod_in)
{
    if( timeline_in == 0)
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

    for( list<Particles>::iterator iter= ptrParticlesList_in->begin(); iter!=ptrParticlesList_in->end(); ++iter)
    {
        // locate the particle
        Particles temp = *iter;
        // get the weighting info of this particle
        struct structg tempStr = temp.InttoStrp1();
        // get the info of mass(weight of each simulation particle)
        double tempMass = temp.WeightMi();

 //       std::cout << temp.VelParticles().x() << " " << temp.VelParticles().y() << " " << temp.VelParticles().z() << std::endl;
 //       std::cout << temp.WeightMi() << " <== " << std::endl;
 //       int pause;
  //      std::cin >>pause;
//        std::cout << "checkUp "<< tempStr.face << " " << tempStr.ig << " " << tempStr.jg << " " << tempStr.kg << " " << tempStr.iw << " " << tempStr.jw << " " << tempStr.kw << " " ;
//        std::cout << "str.v "<< tempStr.vx << " " << tempStr.vy << " " << tempStr.vz << " " ;
//        std::cout << "tempVel "<< tempVel.x() << " " << tempVel.y() << " " << tempVel.z() << " ";
//        std::cout << "tempMass "<< tempMass << std::endl;
//        int pause ;
  //      std::cin >> pause;
//    std::cout << tempStr.face <<  tempStr.ig+1 << tempStr.jg+1 << tempStr.kg << " " << tempStr.iw << tempStr.jw << tempStr.kw << " mass " << tempMass << " xxxx "; ;
        // get the info of velocity
        Vector3 tempVel = temp.VelParticles();
        // update density and velocity in the grids
        ptrArray_in[tempStr.face][tempStr.ig+1][tempStr.jg+1][tempStr.kg]->UpdateDueToWgt( cellSize1 - tempStr.iw, cellSize1 - tempStr.jw, cellSize1 - tempStr.kw, tempMass, tempVel);
        ptrArray_in[tempStr.face][tempStr.ig+2][tempStr.jg+1][tempStr.kg]->UpdateDueToWgt( tempStr.iw + 1,         cellSize1 - tempStr.jw, cellSize1 - tempStr.kw, tempMass, tempVel);
        ptrArray_in[tempStr.face][tempStr.ig+1][tempStr.jg+2][tempStr.kg]->UpdateDueToWgt( cellSize1 - tempStr.iw, tempStr.jw + 1,         cellSize1 - tempStr.kw, tempMass, tempVel);
        ptrArray_in[tempStr.face][tempStr.ig+2][tempStr.jg+2][tempStr.kg]->UpdateDueToWgt( tempStr.iw + 1,         tempStr.jw + 1,         cellSize1 - tempStr.kw, tempMass, tempVel);
        
        ptrArray_in[tempStr.face][tempStr.ig+1][tempStr.jg+1][tempStr.kg+1]->UpdateDueToWgt( cellSize1 - tempStr.iw, cellSize1 - tempStr.jw, tempStr.kw + 1,         tempMass, tempVel);
        ptrArray_in[tempStr.face][tempStr.ig+2][tempStr.jg+1][tempStr.kg+1]->UpdateDueToWgt( tempStr.iw + 1,         cellSize1 - tempStr.jw, tempStr.kw + 1,         tempMass, tempVel);
        ptrArray_in[tempStr.face][tempStr.ig+1][tempStr.jg+2][tempStr.kg+1]->UpdateDueToWgt( cellSize1 - tempStr.iw, tempStr.jw + 1,         tempStr.kw + 1,         tempMass, tempVel);
        ptrArray_in[tempStr.face][tempStr.ig+2][tempStr.jg+2][tempStr.kg+1]->UpdateDueToWgt( tempStr.iw + 1,         tempStr.jw + 1,         tempStr.kw + 1,         tempMass, tempVel);
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
//    int pause ;
  //  std::cin >> pause;
    }

    // go through temp particles list Bottom
    for( list<Particles>::iterator iter= ptrParticlesListBottom_in->begin(); iter!=ptrParticlesListBottom_in->end(); ++iter)
    {
        // locate the particle
        Particles temp = *iter;
        // get the weighting info of this particle
        struct structg tempStr = temp.InttoStrp1();
        // get the info of mass(weight of each simulation particle)
        double tempMass = temp.WeightMi();

        Vector3 tempVel = temp.VelParticles();
   
        ptrArray_in[tempStr.face][tempStr.ig+1][tempStr.jg+1][tempStr.kg+1]->UpdateDueToWgt( cellSize1 - tempStr.iw, cellSize1 - tempStr.jw, tempStr.kw + 1,         tempMass, tempVel);
        ptrArray_in[tempStr.face][tempStr.ig+2][tempStr.jg+1][tempStr.kg+1]->UpdateDueToWgt( tempStr.iw + 1,         cellSize1 - tempStr.jw, tempStr.kw + 1,         tempMass, tempVel);
        ptrArray_in[tempStr.face][tempStr.ig+1][tempStr.jg+2][tempStr.kg+1]->UpdateDueToWgt( cellSize1 - tempStr.iw, tempStr.jw + 1,         tempStr.kw + 1,         tempMass, tempVel);
        ptrArray_in[tempStr.face][tempStr.ig+2][tempStr.jg+2][tempStr.kg+1]->UpdateDueToWgt( tempStr.iw + 1,         tempStr.jw + 1,         tempStr.kw + 1,         tempMass, tempVel);

    }

    // go through temp particles list top
    for( list<Particles>::iterator iter= ptrParticlesListTop_in->begin(); iter!=ptrParticlesListTop_in->end(); ++iter)
    {
        // locate the particle
        Particles temp = *iter;
        // get the weighting info of this particle
        struct structg tempStr = temp.InttoStrp1();
        // get the info of mass(weight of each simulation particle)
        double tempMass = temp.WeightMi();

        Vector3 tempVel = temp.VelParticles();

        ptrArray_in[tempStr.face][tempStr.ig+1][tempStr.jg+1][tempStr.kg]->UpdateDueToWgt( cellSize1 - tempStr.iw, cellSize1 - tempStr.jw, cellSize1 - tempStr.kw, tempMass, tempVel);
        ptrArray_in[tempStr.face][tempStr.ig+2][tempStr.jg+1][tempStr.kg]->UpdateDueToWgt( tempStr.iw + 1,         cellSize1 - tempStr.jw, cellSize1 - tempStr.kw, tempMass, tempVel);
        ptrArray_in[tempStr.face][tempStr.ig+1][tempStr.jg+2][tempStr.kg]->UpdateDueToWgt( cellSize1 - tempStr.iw, tempStr.jw + 1,         cellSize1 - tempStr.kw, tempMass, tempVel);
        ptrArray_in[tempStr.face][tempStr.ig+2][tempStr.jg+2][tempStr.kg]->UpdateDueToWgt( tempStr.iw + 1,         tempStr.jw + 1,         cellSize1 - tempStr.kw, tempMass, tempVel);
    
    }



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
    H5std_string MEMBERx( "x_name");
    H5std_string MEMBERy( "y_name");
    H5std_string MEMBERz( "z_name");

    H5std_string MEMBER1( "pos3_name");
    H5std_string MEMBER2( "e3_name");
    H5std_string MEMBER3( "b3_name");
    H5std_string MEMBER4( "v3_name");
    H5std_string MEMBER5( "density_name");

    const int RANK = 4;

    typedef struct Vector3_h5{
        double v_x;
        double v_y;
        double v_z;
    } Vector3_h5;

    typedef struct GridsPoints_h5{
        Vector3_h5 pos3;
        Vector3_h5 e3;
        Vector3_h5 b3;
        Vector3_h5 v3;
        double density;
    } GridsPoints_h5;


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
                    array_data[face][i][j][k].pos3 
                        = {ptrArray_in[face][i+1][j+1][k]->Pos3().x(), ptrArray_in[face][i+1][j+1][k]->Pos3().y(), ptrArray_in[face][i+1][j+1][k]->Pos3().z()} ;

                    array_data[face][i][j][k].e3 
                        = {ptrArray_in[face][i+1][j+1][k]->E3().x(), ptrArray_in[face][i+1][j+1][k]->E3().y(), ptrArray_in[face][i+1][j+1][k]->E3().z()} ;

                    array_data[face][i][j][k].b3 
                        = {ptrArray_in[face][i+1][j+1][k]->B3().x(), ptrArray_in[face][i+1][j+1][k]->B3().y(), ptrArray_in[face][i+1][j+1][k]->B3().z()} ;

                    array_data[face][i][j][k].v3 
                        = {ptrArray_in[face][i+1][j+1][k]->Vel3().x(), ptrArray_in[face][i+1][j+1][k]->Vel3().y(), ptrArray_in[face][i+1][j+1][k]->Vel3().z()} ;
                    array_data[face][i][j][k].density
                        = ptrArray_in[face][i+1][j+1][k]->Density();

                 //   std::cout << array_data[face][i][j][k].b3.v_x << std::endl;
                }
            }
        }
    }
    

    Exception::dontPrint();


    hsize_t dim[] = {totalFace,fieldsGridsSize+1,fieldsGridsSize+1,fieldsGridsSize+1};
    DataSpace space( RANK, dim);

    if(h5FileCheck_in == 0)
    {
        H5File* file = new H5File( FILE_NAME, H5F_ACC_TRUNC);
        h5FileCheck = 1;

        CompType mtype_vector3( sizeof( Vector3_h5));
        mtype_vector3.insertMember( MEMBERx, HOFFSET(Vector3_h5,v_x), PredType::NATIVE_DOUBLE);
        mtype_vector3.insertMember( MEMBERy, HOFFSET(Vector3_h5,v_y), PredType::NATIVE_DOUBLE);
        mtype_vector3.insertMember( MEMBERz, HOFFSET(Vector3_h5,v_z), PredType::NATIVE_DOUBLE);

        CompType mtype_grids( sizeof( GridsPoints_h5));
        mtype_grids.insertMember( MEMBER1, HOFFSET(GridsPoints_h5,pos3), mtype_vector3);
        mtype_grids.insertMember( MEMBER2, HOFFSET(GridsPoints_h5,e3), mtype_vector3);
        mtype_grids.insertMember( MEMBER3, HOFFSET(GridsPoints_h5,b3), mtype_vector3);
        mtype_grids.insertMember( MEMBER4, HOFFSET(GridsPoints_h5,v3), mtype_vector3);
        mtype_grids.insertMember( MEMBER5, HOFFSET(GridsPoints_h5,density), PredType::NATIVE_DOUBLE);

        DataSet* dataset;
        dataset = new DataSet(file->createDataSet(DATASET_NAME, mtype_grids, space));
        dataset->write( array_data[0][0][0], mtype_grids);

        delete dataset;
        delete file;
        delete data_mem;
    }
    else
    {
        H5File* file = new H5File( FILE_NAME, H5F_ACC_RDWR);

        CompType mtype_vector3( sizeof( Vector3_h5));
        mtype_vector3.insertMember( MEMBERx, HOFFSET(Vector3_h5,v_x), PredType::NATIVE_DOUBLE);
        mtype_vector3.insertMember( MEMBERy, HOFFSET(Vector3_h5,v_y), PredType::NATIVE_DOUBLE);
        mtype_vector3.insertMember( MEMBERz, HOFFSET(Vector3_h5,v_z), PredType::NATIVE_DOUBLE);

        CompType mtype_grids( sizeof( GridsPoints_h5));
        mtype_grids.insertMember( MEMBER1, HOFFSET(GridsPoints_h5,pos3), mtype_vector3);
        mtype_grids.insertMember( MEMBER2, HOFFSET(GridsPoints_h5,e3), mtype_vector3);
        mtype_grids.insertMember( MEMBER3, HOFFSET(GridsPoints_h5,b3), mtype_vector3);
        mtype_grids.insertMember( MEMBER4, HOFFSET(GridsPoints_h5,v3), mtype_vector3);
        mtype_grids.insertMember( MEMBER5, HOFFSET(GridsPoints_h5,density), PredType::NATIVE_DOUBLE);

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
// at each cell
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
    //            std::cout << s<< " " << i << j << k << std::endl;
            }
        }
    }
    return VolumeCellsArray;
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
    static int timeLineLimit = 1;
    static int printTimePeriod = 1;
    static int updateInfoPeriod = 1;
    // Prerun 1.0 // Create Grids, including B and Pos. And then Velocity (corotation) and N (exponential )
    // Prerun 1.1 // And then E (electron momentum equation).
    cout << " Create Grids" << endl;
    GridsPoints***** ptrArray =  GridsCreation();

    // Prerun 1.3 // Create Cell centered field array for nesseary calculation
    Vector3*** ptrVectorCellArray = VectorCellField();  

    // Prerun 1.4 // Create grids field array of volum for one face
    cout << " Create array of Volume at cells and at grids" << endl;
    double*** ptrVolumeCellArray = VolumeCellsField( ptrArray);
    
    cout << " Create array of Volume at cells and at grids" << endl;
    double*** ptrVolumeGridArray = VolumeGridsField( ptrVolumeCellArray);
    
    // Prerun 1.2 // Create particles list, initialize the velocity and position of each particles
    cout << " Create particles list of main domain" << endl;
    list<Particles>* ptrParticlesList = ParticlesLists( ptrArray, ptrVolumeCellArray);
    // Run 2.0
    
    cout << " Start" << endl;
    for( int timeline = 0; timeline <= timeLineLimit; timeline++)
    {
    std::cout << "timeline" << timeline << std::endl;

    // Run 2.1 // Particles in main domain
        for( list<Particles>::iterator iteratorM = ptrParticlesList->begin(); iteratorM != ptrParticlesList->end(); ++iteratorM)
        {
//std::cout << " check 0" << std::endl;
            Particles temp = *iteratorM;
            struct structg tempStruct = temp.InttoStrp1();
            
            int check; // check whether in the main domain or not, "0" means in "1" means out
            //test function//  double xxx= temp.VelParticles().x();    //cout << xxx << " ";
            // update velocity // update position
            check = temp.BorisMethod( &tempStruct, ptrArray);
            // check if still in the main domain
            if( check == 1) // out of the domain
            {
                ptrParticlesList->erase( iteratorM);
                iteratorM--;
            }
        }

    // Run 2.2 // Create temp particle lists    
    list<Particles>* ptrParticlesListTop = ParticlesListsTemp( ptrArray, "T");
    list<Particles>* ptrParticlesListBottom = ParticlesListsTemp( ptrArray, "B");

    // Run 2.3 // Particles in temp domain
        for( list<Particles>::iterator iteratorT = ptrParticlesListTop->begin(); iteratorT != ptrParticlesListTop->end(); ++iteratorT)
        {
            Particles temp = *iteratorT;
            struct structg tempStruct = temp.InttoStrp1();
            int check; // check whether in the main domain or not, "0" means in "1" means out
            // update velocity  // update position
            check = temp.BorisMethod( &tempStruct, ptrArray);
            // check if still in the main domain
            if( check == 0) // in the domain
            {   
                ptrParticlesList->push_back( temp);
                ptrParticlesListTop->erase( iteratorT);
                iteratorT--;
            }
        }

        for( list<Particles>::iterator iteratorB = ptrParticlesListBottom->begin(); iteratorB != ptrParticlesListBottom->end(); ++iteratorB)
        {
            Particles temp = *iteratorB;
            struct structg tempStruct = temp.InttoStrp1();
            int check; // check whether in the main domain or not, "0" means in "1" means out
            // update velocity  // update position
            check = temp.BorisMethod( &tempStruct, ptrArray);
            // check if still in the main domain
            if( check == 0) // in the domain
            {   
                ptrParticlesList->push_back( temp);
                ptrParticlesListBottom->erase( iteratorB);
                iteratorB--;
            }
        }       

                



    std::cout << " Update gridspoints info" << std::endl;
    // Run 2.5 // Update info in grids 
    // Run 2.5.1 // Accumulate density and velocity per timestep and average them to get the
    // value of density and velocity per updatetimeperiod
    UpdateInfoGrids( ptrArray, ptrParticlesList, ptrParticlesListBottom, ptrParticlesListTop,
                    ptrVolumeGridArray, timeline, updateInfoPeriod);
    
    // Run 2.4 // delete temp particlesLists
    delete ptrParticlesListTop;
    delete ptrParticlesListBottom;





    if( timeline % updateInfoPeriod ==0)
    {
        // Run 2.5.2 // Update E
        for( int face = 0; face < 6; face++)
        {
            ptrVectorCellArray = ValueGradientPe( ptrVectorCellArray, ptrVolumeCellArray, ptrArray, face);
            updateE_nocurrent( ptrVectorCellArray, ptrArray, face);
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
    delete ptrParticlesList;

}


//************************************************************************
//************************************************************************
// FUNCTION
// Calculation the temprature
//************************************************************************
//************************************************************************
void Titheridge_Te(GridsPoints***** ptrArray_in)
{
    
    const double h0 = 400.0, R0 = 1.0627825, R02 = 1.1295067, por = 0.2857142, Re = 6371.2,
                 x1 = 1.0470869;
    const double a0 = 1.23e3, a1 = 2.2e3, a2 = -3.29e3, a3 = -2.6e-1, a4 = -6.8e-1,
                 b0 = 4.65e0, b1 = -8.55e0, b2 = 4.14e0, b3 = -2.16e0, b4 = 1.45e0;
    const double aa0 = 9.85e2, aa1 = -9.63e2, aa2 = 1.125e3, aa3 = -0.6e0, aa4 = 0.10e0,
                 bb0 = 7.56e-1, bb1 = -8.8e-1, bb2 = 2.9e-1, bb3 = -2.63e0, bb4 = 1.84e0;


    //assume DT = 0, DOY = 0;
    int DOY= 0;
    double DT = 0.0;
    double F107 = 0.0;
    double Kp = 0.0;

    
    for( int face = 0; face < totalFace; face++)
    {
        for( int i = 1; i < fieldsGridsSize+2; i++)
        {
            for( int j = 1; j < fieldsGridsSize+2; j++)
            {
                for( int k = 1; k < fieldsGridsSize; k++)
                { 
                    double alt = ptrArray_in[face][i][j][k]->Pos3().norm() / 1000.0;
                    double L = alt / Re;
                    double GLOND = atan( ptrArray_in[face][i][j][k]->Pos3().y() / ptrArray_in[face][i][j][k]->Pos3().x()) * 57.2957795;
                    double GLATD = asin( ptrArray_in[face][i][j][k]->Pos3().z() / ptrArray_in[face][i][j][k]->Pos3().norm()) * 57.2957795;
                    double LT = GLOND /15;
    // part 1
                    double heq = (L-1.0e0)*Re;
                    double SL = 1.0e0 - x1 / L;
                    //geomagnetic latitude an an altitude of 300km
                    double Lat = asin( sqrt( SL)) * 57.29578; // in degree
                    // day time T0 and G0 Eq.(19) in Titheridge, JGR1998
                    double T0d=(a0+a1*SL+a2* pow (SL,2))/(1.0e0+a3*SL+a4* pow(SL,2));
                    double G0d=(b0+b1*SL+b2*pow(SL,2))/(1.0e0+b3*SL+b4*pow(SL,2));
                    // night time T0 and G0
                    double T0n=(aa0+aa1*SL+aa2 * pow(SL,2))/(1.0e0+aa3*SL+aa4* pow(SL,2));
                    double G0n=(bb0+bb1*SL+bb2 * pow(SL,2))/(1.0e0+bb3*SL+bb4* pow(SL,2));

    // part 2
                    // change in solar activity
                    double delta_T0 = 3.4 * (F107 - 120.0);
                    G0d=G0d*T0d/(T0d+delta_T0);
                    T0d=T0d+delta_T0;
                    G0n=G0n*T0n/(T0n+delta_T0);
                    T0n=T0n+delta_T0;
                    // change in magnetic activity
                    if(Lat < 46.0e0) 
                    {
                        double z=0.135*(Lat-46.0)/57.29578;
                        double Dk=37.0+1.33*(Lat-46.0)-37.0*cos(z);
                        // change in magnetic activity
                        delta_T0=Dk*(Kp-1.5);
                        G0d=G0d * pow(T0d/(T0d+delta_T0),2.5);
                        T0d=T0d+delta_T0;
                        G0n=G0n * pow(T0n/(T0n+delta_T0),2.5);
                        T0n=T0n+delta_T0;        
                    }

    // PART3 
                    double G0T0d = G0d / T0d;
                    double G0T0n = G0n / T0n;


                    double alg=log(alt/h0);
                    double R2=pow(1.0e0+alt/Re,2.0);
                    double Bh=0.05e0/(2.0e0*L-R0)*(88.0e0+alg*(10.5e0-alg));
                    // mean day values from Eq.(13) in Titheridge, 1998 JGR 
                    double Tday0=T0d * pow(1.0e0+Bh*G0T0d*((heq-h0)/R02-(heq-alt)/R2),por);
                    // mean night values from Eq.(13) in Titheridge, 1998 JGR 
                    double Tnig0=T0n * pow(1.0e0+Bh*G0T0n*((heq-h0)/R02-(heq-alt)/R2),por);

    // PART 4: 
    //local solar time (SAT) SAT need: day of year，local time，GLOND, GLATD and local time(LT), and different time with Greenwich Time (DT)，
	//see https://pvcdrom.pveducation.org/SUNLIGHT/SOLART.HTM
	                double B=300/365*(DOY-81);
	                double EOT=9.87*sin(2*B)-7.53*cos(B)-1.5*sin(B);
	                double LSTM=15*DT;   // 15 degree per hour
	                double TC=4*(LSTM-GLOND)+EOT;
                    double SAT=LT+TC/60;
	                if(SAT < 0.0e0) SAT=SAT+24.0e0;
                    if(SAT > 24.0e0) SAT=SAT-24.0e0;
                    // degrees to radians
                    double GLATR=GLATD/57.2957795;
    // DEC -- solar declination angle in radians
	// solar declination angle need: day of year in unit of dgree, notice: dgree to radians((180/2pi)~57.2957795)
	//see https://www.sciencedirect.com/topics/engineering/solar-declination
	                double DEC=23.45*sin(360/365*(284+DOY))/57.2957795;
                    double S_Long=-1.0e0*tan(GLATR)*tan(DEC);
                    double s_r;
                    if(S_Long > -1.0e0 && S_Long < 1.0e0)
                    {
                    // Local sunrise
                        s_r=12.0e0-3.82*acos(S_Long);
                    } else 
                    {
                        if(S_Long < -1.0e0)
                    {
                        s_r = 0.0e0;
                    } else
                    {   
                        s_r=12.0e0;
                    }
                    }
                    //duration of the day-night transition 
                    double Delta_tr=1.2+0.5*SL;
                    double Delta_ts=Delta_tr+0.9;
                    // time of the center of the transition in hours after ground sunrise or sunset
                    double t_s=0.15*s_r;
                    if(t_s > 1.0e0) t_s=1.0;
                    double t_r=-0.5*t_s;      
                    double D;  
                    if(SAT < 12.0e0){
                        D=(s_r-0.5*t_s-SAT)/Delta_tr;
                    } else {
                        D=(SAT-t_s+s_r-24.0e0)/Delta_ts;
                    }
                    double phi=1.0/(1.0+exp(3.2*D));
                    double y=s_r-12.5;
                    if(y < -4.0) y=-4.0;
                    double TDvar=0.97+0.22* pow((SAT-12.5)/y, 2.0);

                    double Tday;    
                    if(TDvar < 1.2) {
                        Tday=TDvar*Tday0;
                    } else {
                        Tday=1.2*Tday0;
                    }
       
        
                    if(SAT < 12.0e0) {
                        y=SAT;
                    } else {
                        y=SAT-24.0;
                    }
                    double TNvar=0.83+0.15*exp(-0.2*y*(1.4-0.8*SL));
       
                    double Tnig;   
                    if(TNvar < 1.2){
                        Tnig=TNvar*Tnig0;
                    } else {
                        Tnig=1.2*Tnig0;
                    }
                    
                    double Tes=Tnig+(Tday-Tnig)*phi;
                    ptrArray_in[face][i][j][k]->SetTemperature(Tes);

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

