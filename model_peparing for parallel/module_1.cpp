#include<iostream>
#include <list>
#include <memory>
#include <string>
#include "parameters.h"
#include "fieldsgrids.h"
#include "particles.h"
#include "vector3.h"
#include "module_1.h"
#include <cmath>
#include "H5Cpp.h"
#include <bitset>

using std::cout;
using std::endl;
using std::shared_ptr;
using std::make_shared;


//************************************************************************
//************************************************************************
// FUNCTION
// Initialization the particles for velocity and position 
// Generate lists of particles
//************************************************************************
//************************************************************************
vector<Particles>* ParticlesLists( GridsPoints***** ptrArray_in, double*** ptrVolumeCellArray_in, double mi0, double N0)
{
   //list<Particles> listP;
//    auto listsPtr = make_shared<list<Particles>>();
    vector<Particles>* listsPtr = new vector<Particles>;
//   shared_ptr<Particles> p1 = make_shared<Particles> (); // test
//    double scaleHeight = ikT / mi0 / gravity; // assume the T is const for initiallization

    for( int i = 0; i < totalFace; i++)
    {
        for( int j = 1; j <= fieldsGridsSize; j++)
        {
            for( int k = 1; k <= fieldsGridsSize; k++)
            {
                for( int s = 1; s <= fieldsGridsSize-2; s++)
                {
                    // number of real particles ( notice the unit is number density)
                    //double N = N0 * exp(-1.0 * (ptrArray_in[i][j][k][s]->Pos3().norm() - radius) / scaleHeight);
                    double r = ptrArray_in[i][j][k][s]->Pos3().norm() / radius;
                    double N = N0 / r * ( 1.0 - tanh( r - 6.5));
                    // weightNi of each simulation particle
                    double Ni_simu = N / iniParticleNumberPerCell * ptrVolumeCellArray_in[j][k][s];

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
                    Vector3 tempVector3 = Uint64ToVector3 ( intPos);        

                    // calculate random velocity
         //           std::cout<< MaxwellDisV( ikT, ptrArray_in[i][j][k][s]->Vel3().x()) << std::endl;
         //           std:: cout << ptrArray_in[i][j][k][s]->B3().norm()<< " B "<< std::endl;                        
                    

                    Vector3 vVel = Vector3( MaxwellDisV( ikT, ptrArray_in[i][j][k][s]->Vel3().x(), mi0),
                                            MaxwellDisV( ikT, ptrArray_in[i][j][k][s]->Vel3().y(), mi0),
                                            MaxwellDisV( ikT, ptrArray_in[i][j][k][s]->Vel3().z(), mi0));
                    double mu_simu = MaxwellDisEnergy( ptrArray_in, intPos);
                    // put the particles at the end of list
                    Particles tempP= Particles(intPos, tempVector3, vVel, Ni_simu, mu_simu);
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
// Transform the ubit64 to vector of position
//************************************************************************
//************************************************************************
Vector3 Uint64ToVector3( uint_64 intPos_in)
{
    double px, py, pz;
    double temp[2];
    uint_64 posUint = intPos_in;
    uint_64 face = 0, ip = 0, jp = 0, kp = 0;

    face = posUint >> 61;
    for( int i = 0; i < particlesGridsLevel; i++) 
    {
        ip = (ip << 1) + ((posUint >> 60   - i*3) & 1);
        jp = (jp << 1) + ((posUint >> 60-1 - i*3) & 1);
        kp = (kp << 1) + ((posUint >> 60-2 - i*3) & 1);
    }
    // 2. transfor to double x y z
    // 2.1 radial

    // cellSize1 = /particlesGridsSize * fieldsGridsSize
    // particles located at center of cell
    double L = LMin * pow(10, logRatio *  ( (kp +0.5)/ cellSize1 )); 
    
    // 2.2 IgJg to ST note 0<ST<1
    temp[0] = (1.0 / particlesGridsSize) * ip;
    temp[1] = (1.0 / particlesGridsSize) * jp;
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
    double k = L * radius / sqrt(pow(1.0,2.0) + pow(temp[0],2.0) + pow(temp[1],2.0));
    switch (face)
    {
        case 0: px=1.0;           py=temp[0];     pz=temp[1]; break;
        case 1: px=-1.0*temp[0];  py=1.0;         pz=temp[1]; break;
        case 2: px=-1.0*temp[1];  py=temp[0];     pz=1.0;     break;
        case 3: px=-1.0;          py=-1.0*temp[0];pz=temp[1]; break;
        case 4: px=temp[0];       py=-1.0;        pz=temp[1]; break;
        default:px=temp[1];       py=temp[0];     pz=-1.0;    break;
    }
    
    px *= k; py *= k; pz *= k;
    Vector3 tempV = Vector3( px, py, pz);
    return tempV;
}


//************************************************************************
//************************************************************************
// FUNCTION
// Initialization the particles for velocity and position 
// Generate lists of particles for bot and top region temp
//************************************************************************
//************************************************************************
vector<Particles>* ParticlesListsTemp( GridsPoints***** ptrArray_in, double*** ptrVolumeCellArray_in, double mi0, int ionType_in)
{
    vector<Particles>* listsPtrTemp = new vector<Particles>;
    double N, Ni_simu;
   
    
    // for bottom temp domain
    for( int i = 0; i < totalFace; i++)
    {
        for( int j = 1; j <= fieldsGridsSize; j++)
        {
            for( int k =1; k<= fieldsGridsSize; k++)
            {
                int s = 0;
                // number density
                switch (ionType_in)
                {
                case 1:
                    N = ptrArray_in[i][j][k][s]->Density_H() * 3.0 * ratioH_bot;
                    break;
                case 4:
                    N = ptrArray_in[i][j][k][s]->Density_He()* 3.0 * ratioHe_bot;
                    break;                
                default:
                    N = ptrArray_in[i][j][k][s]->Density_O()* 3.0 * ratioO_bot;
                    break;
                }
                // mass of each simulation particle 
                Ni_simu = N / tempParticleNumberPerCell * ptrVolumeCellArray_in[j][k][s];
                for ( int t = 1; t <= tempParticleNumberPerCell; t++)
                {
                // calculate random position
                uint_64 intPos = UniDisInCell( ptrArray_in, i, j, k, s);
                Vector3 tempVector3 = Uint64ToVector3 ( intPos);   

                // calculate random velocity
                Vector3 vVel = Vector3( MaxwellDisV( ikT, ptrArray_in[i][j][k][s]->Vel3().x(), mi0),
                                            MaxwellDisV( ikT, ptrArray_in[i][j][k][s]->Vel3().y(), mi0),
                                            MaxwellDisV( ikT, ptrArray_in[i][j][k][s]->Vel3().z(), mi0));
                double mu_simu = MaxwellDisEnergy( ptrArray_in, intPos);
                // put the particles at the end of list
                Particles tempP= Particles(intPos, tempVector3, vVel, Ni_simu, mu_simu);
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
                switch (ionType_in)
                {
                case 1:
                    N = ptrArray_in[i][j][k][s]->Density_H()* 3.0 * ratioH_top;
                    break;
                case 4:
                    N = ptrArray_in[i][j][k][s]->Density_He()* 3.0 * ratioHe_top;
                    break;                
                default:
                    N = ptrArray_in[i][j][k][s]->Density_O()* 3.0 * ratioO_top;
                    break;
                }
                

                // mass of each simulation particle 
                Ni_simu = N / tempParticleNumberPerCell *mi0 * ptrVolumeCellArray_in[j][k][s];
                for ( int t = 1; t <= tempParticleNumberPerCell; t++)
                {
                // calculate random position
                uint_64 intPos = UniDisInCell( ptrArray_in, i, j, k, s);
                Vector3 tempVector3 = Uint64ToVector3 ( intPos);   

                // calculate random velocity
                Vector3 vVel = Vector3( MaxwellDisV( ikT, ptrArray_in[i][j][k][s]->Vel3().x(), mi0),
                                            MaxwellDisV( ikT, ptrArray_in[i][j][k][s]->Vel3().y(), mi0),
                                            MaxwellDisV( ikT, ptrArray_in[i][j][k][s]->Vel3().z(), mi0));
                double mu_simu = MaxwellDisEnergy( ptrArray_in, intPos);
                // put the particles at the end of list
                Particles tempP= Particles(intPos, tempVector3, vVel, Ni_simu, mu_simu);
                listsPtrTemp->push_back( tempP);           
                }
            }
        }
    }


    return listsPtrTemp;
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
                      vector<Particles>* ptrParticlesList_H_in, 
                      vector<Particles>* ptrParticlesList_He_in,
                      vector<Particles>* ptrParticlesList_O_in,
                      vector<Particles>* ptrParticlesListTemp_H_in,
                      vector<Particles>* ptrParticlesListTemp_He_in,
                      vector<Particles>* ptrParticlesListTemp_O_in,
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
                    for( int k = 1 + tempGridsCellLevel; k < fieldsGridsSize - tempGridsCellLevel; k++)
                    {  
                        ptrArray_in[face][i][j][k]->ResetParameters();
                    }
                }
            }
        }
    }

// For H particles in main domain    
    for( auto iter= ptrParticlesList_H_in->begin(); iter!=ptrParticlesList_H_in->end(); ++iter)
    {
        // locate the particle
        Particles temp = *iter;
        // get the weighting info of this particle
        struct structg tempStr = temp.InttoStrp1();
        // get the info of mass(weight of each simulation particle)
        double tempNumber = temp.WeightNi();

/*      // std::cout << temp.VelParticles().x() << " " << temp.VelParticles().y() << " " << temp.VelParticles().z() << std::endl;
        // std::cout << temp.WeightMi() << " <== " << std::endl;
        // int pause;
        // std::cin >>pause;
        // std::cout << "checkUp "<< tempStr.face << " " << tempStr.ig << " " << tempStr.jg << " " << tempStr.kg << " " << tempStr.iw << " " << tempStr.jw << " " << tempStr.kw << " " ;
        // std::cout << "str.v "<< tempStr.vx << " " << tempStr.vy << " " << tempStr.vz << " " ;
        // std::cout << "tempVel "<< tempVel.x() << " " << tempVel.y() << " " << tempVel.z() << " ";
        // std::cout << "tempNumber "<< tempNumber << std::endl;
        // int pause ;
        // std::cin >> pause;
        // std::cout << tempStr.face <<  tempStr.ig+1 << tempStr.jg+1 << tempStr.kg << " " << tempStr.iw << tempStr.jw << tempStr.kw << " mass " << tempNumber << " xxxx "; ;
        // get the info of velocity
*/        
        Vector3 tempVel = temp.VelParticles();
        // update density and velocity in the grids

        if( tempStr.kg > tempGridsCellLevel ){
        ptrArray_in[tempStr.face][tempStr.ig+1][tempStr.jg+1][tempStr.kg]->UpdateDueToWgt( cellSize1 - tempStr.iw, cellSize1 - tempStr.jw, cellSize1 - tempStr.kw, tempNumber, tempVel, 1);
        ptrArray_in[tempStr.face][tempStr.ig+2][tempStr.jg+1][tempStr.kg]->UpdateDueToWgt( tempStr.iw + 1,         cellSize1 - tempStr.jw, cellSize1 - tempStr.kw, tempNumber, tempVel, 1);
        ptrArray_in[tempStr.face][tempStr.ig+1][tempStr.jg+2][tempStr.kg]->UpdateDueToWgt( cellSize1 - tempStr.iw, tempStr.jw + 1,         cellSize1 - tempStr.kw, tempNumber, tempVel, 1);
        ptrArray_in[tempStr.face][tempStr.ig+2][tempStr.jg+2][tempStr.kg]->UpdateDueToWgt( tempStr.iw + 1,         tempStr.jw + 1,         cellSize1 - tempStr.kw, tempNumber, tempVel, 1);
        }
        if( tempStr.kg < fieldsGridsSize - 1 - tempGridsCellLevel ){
        ptrArray_in[tempStr.face][tempStr.ig+1][tempStr.jg+1][tempStr.kg+1]->UpdateDueToWgt( cellSize1 - tempStr.iw, cellSize1 - tempStr.jw, tempStr.kw + 1,         tempNumber, tempVel, 1);
        ptrArray_in[tempStr.face][tempStr.ig+2][tempStr.jg+1][tempStr.kg+1]->UpdateDueToWgt( tempStr.iw + 1,         cellSize1 - tempStr.jw, tempStr.kw + 1,         tempNumber, tempVel, 1);
        ptrArray_in[tempStr.face][tempStr.ig+1][tempStr.jg+2][tempStr.kg+1]->UpdateDueToWgt( cellSize1 - tempStr.iw, tempStr.jw + 1,         tempStr.kw + 1,         tempNumber, tempVel, 1);
        ptrArray_in[tempStr.face][tempStr.ig+2][tempStr.jg+2][tempStr.kg+1]->UpdateDueToWgt( tempStr.iw + 1,         tempStr.jw + 1,         tempStr.kw + 1,         tempNumber, tempVel, 1);
        }
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
    for( auto iter= ptrParticlesList_He_in->begin(); iter!=ptrParticlesList_He_in->end(); ++iter)
    {
        // locate the particle
        Particles temp = *iter;
        // get the weighting info of this particle
        struct structg tempStr = temp.InttoStrp1();
        // get the info of mass(weight of each simulation particle)
        double tempNumber = temp.WeightNi();
        // get the info of velocity
        Vector3 tempVel = temp.VelParticles();
        // update density and velocity in the grids
        if( tempStr.kg > tempGridsCellLevel ){
        ptrArray_in[tempStr.face][tempStr.ig+1][tempStr.jg+1][tempStr.kg]->UpdateDueToWgt( cellSize1 - tempStr.iw, cellSize1 - tempStr.jw, cellSize1 - tempStr.kw, tempNumber, tempVel, 4);
        ptrArray_in[tempStr.face][tempStr.ig+2][tempStr.jg+1][tempStr.kg]->UpdateDueToWgt( tempStr.iw + 1,         cellSize1 - tempStr.jw, cellSize1 - tempStr.kw, tempNumber, tempVel, 4);
        ptrArray_in[tempStr.face][tempStr.ig+1][tempStr.jg+2][tempStr.kg]->UpdateDueToWgt( cellSize1 - tempStr.iw, tempStr.jw + 1,         cellSize1 - tempStr.kw, tempNumber, tempVel, 4);
        ptrArray_in[tempStr.face][tempStr.ig+2][tempStr.jg+2][tempStr.kg]->UpdateDueToWgt( tempStr.iw + 1,         tempStr.jw + 1,         cellSize1 - tempStr.kw, tempNumber, tempVel, 4);
        }
        if( tempStr.kg < fieldsGridsSize - 1 - tempGridsCellLevel ){
        ptrArray_in[tempStr.face][tempStr.ig+1][tempStr.jg+1][tempStr.kg+1]->UpdateDueToWgt( cellSize1 - tempStr.iw, cellSize1 - tempStr.jw, tempStr.kw + 1,         tempNumber, tempVel, 4);
        ptrArray_in[tempStr.face][tempStr.ig+2][tempStr.jg+1][tempStr.kg+1]->UpdateDueToWgt( tempStr.iw + 1,         cellSize1 - tempStr.jw, tempStr.kw + 1,         tempNumber, tempVel, 4);
        ptrArray_in[tempStr.face][tempStr.ig+1][tempStr.jg+2][tempStr.kg+1]->UpdateDueToWgt( cellSize1 - tempStr.iw, tempStr.jw + 1,         tempStr.kw + 1,         tempNumber, tempVel, 4);
        ptrArray_in[tempStr.face][tempStr.ig+2][tempStr.jg+2][tempStr.kg+1]->UpdateDueToWgt( tempStr.iw + 1,         tempStr.jw + 1,         tempStr.kw + 1,         tempNumber, tempVel, 4);
        }
    }


// For O particles in main domain    
    for( auto iter= ptrParticlesList_O_in->begin(); iter!=ptrParticlesList_O_in->end(); ++iter)
    {
        // locate the particle
        Particles temp = *iter;
        // get the weighting info of this particle
        struct structg tempStr = temp.InttoStrp1();
        // get the info of mass(weight of each simulation particle)
        double tempNumber = temp.WeightNi();
        // get the info of velocity
        Vector3 tempVel = temp.VelParticles();
        // update density and velocity in the grids
        if( tempStr.kg > tempGridsCellLevel ){
        ptrArray_in[tempStr.face][tempStr.ig+1][tempStr.jg+1][tempStr.kg]->UpdateDueToWgt( cellSize1 - tempStr.iw, cellSize1 - tempStr.jw, cellSize1 - tempStr.kw, tempNumber, tempVel, 16);
        ptrArray_in[tempStr.face][tempStr.ig+2][tempStr.jg+1][tempStr.kg]->UpdateDueToWgt( tempStr.iw + 1,         cellSize1 - tempStr.jw, cellSize1 - tempStr.kw, tempNumber, tempVel, 16);
        ptrArray_in[tempStr.face][tempStr.ig+1][tempStr.jg+2][tempStr.kg]->UpdateDueToWgt( cellSize1 - tempStr.iw, tempStr.jw + 1,         cellSize1 - tempStr.kw, tempNumber, tempVel, 16);
        ptrArray_in[tempStr.face][tempStr.ig+2][tempStr.jg+2][tempStr.kg]->UpdateDueToWgt( tempStr.iw + 1,         tempStr.jw + 1,         cellSize1 - tempStr.kw, tempNumber, tempVel, 16);
        }
        if( tempStr.kg < fieldsGridsSize - 1 - tempGridsCellLevel ){
        ptrArray_in[tempStr.face][tempStr.ig+1][tempStr.jg+1][tempStr.kg+1]->UpdateDueToWgt( cellSize1 - tempStr.iw, cellSize1 - tempStr.jw, tempStr.kw + 1,         tempNumber, tempVel, 16);
        ptrArray_in[tempStr.face][tempStr.ig+2][tempStr.jg+1][tempStr.kg+1]->UpdateDueToWgt( tempStr.iw + 1,         cellSize1 - tempStr.jw, tempStr.kw + 1,         tempNumber, tempVel, 16);
        ptrArray_in[tempStr.face][tempStr.ig+1][tempStr.jg+2][tempStr.kg+1]->UpdateDueToWgt( cellSize1 - tempStr.iw, tempStr.jw + 1,         tempStr.kw + 1,         tempNumber, tempVel, 16);
        ptrArray_in[tempStr.face][tempStr.ig+2][tempStr.jg+2][tempStr.kg+1]->UpdateDueToWgt( tempStr.iw + 1,         tempStr.jw + 1,         tempStr.kw + 1,         tempNumber, tempVel, 16);
        }
    }

/*
// For H particles in temp domain    
    for( list<Particles>::iterator iter= ptrParticlesListTemp_H_in->begin(); iter!=ptrParticlesListTemp_H_in->end(); ++iter)
    {
        // locate the particle
        Particles temp = *iter;
        // get the weighting info of this particle
        struct structg tempStr = temp.InttoStrp1();
        // get the info of mass(weight of each simulation particle)
        double tempNumber = temp.WeightNi();

        cout << tempNumber << endl;

        Vector3 tempVel = temp.VelParticles();
   
        if( tempStr.kg == 0)
        {    
        ptrArray_in[tempStr.face][tempStr.ig+1][tempStr.jg+1][tempStr.kg+1]->UpdateDueToWgt( cellSize1 - tempStr.iw, cellSize1 - tempStr.jw, tempStr.kw + 1,         tempNumber, tempVel, 1);
        ptrArray_in[tempStr.face][tempStr.ig+2][tempStr.jg+1][tempStr.kg+1]->UpdateDueToWgt( tempStr.iw + 1,         cellSize1 - tempStr.jw, tempStr.kw + 1,         tempNumber, tempVel, 1);
        ptrArray_in[tempStr.face][tempStr.ig+1][tempStr.jg+2][tempStr.kg+1]->UpdateDueToWgt( cellSize1 - tempStr.iw, tempStr.jw + 1,         tempStr.kw + 1,         tempNumber, tempVel, 1);
        ptrArray_in[tempStr.face][tempStr.ig+2][tempStr.jg+2][tempStr.kg+1]->UpdateDueToWgt( tempStr.iw + 1,         tempStr.jw + 1,         tempStr.kw + 1,         tempNumber, tempVel, 1);
        }
        else if( tempStr.kg == fieldsGridsSize - 1 )
        {
        ptrArray_in[tempStr.face][tempStr.ig+1][tempStr.jg+1][tempStr.kg]->UpdateDueToWgt( cellSize1 - tempStr.iw, cellSize1 - tempStr.jw, cellSize1 - tempStr.kw, tempNumber, tempVel, 1);
        ptrArray_in[tempStr.face][tempStr.ig+2][tempStr.jg+1][tempStr.kg]->UpdateDueToWgt( tempStr.iw + 1,         cellSize1 - tempStr.jw, cellSize1 - tempStr.kw, tempNumber, tempVel, 1);
        ptrArray_in[tempStr.face][tempStr.ig+1][tempStr.jg+2][tempStr.kg]->UpdateDueToWgt( cellSize1 - tempStr.iw, tempStr.jw + 1,         cellSize1 - tempStr.kw, tempNumber, tempVel, 1);
        ptrArray_in[tempStr.face][tempStr.ig+2][tempStr.jg+2][tempStr.kg]->UpdateDueToWgt( tempStr.iw + 1,         tempStr.jw + 1,         cellSize1 - tempStr.kw, tempNumber, tempVel, 1);
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
        double tempNumber = temp.WeightNi();

        Vector3 tempVel = temp.VelParticles();
   
        if( tempStr.kg == 0)
        {    
        ptrArray_in[tempStr.face][tempStr.ig+1][tempStr.jg+1][tempStr.kg+1]->UpdateDueToWgt( cellSize1 - tempStr.iw, cellSize1 - tempStr.jw, tempStr.kw + 1,         tempNumber, tempVel, 4);
        ptrArray_in[tempStr.face][tempStr.ig+2][tempStr.jg+1][tempStr.kg+1]->UpdateDueToWgt( tempStr.iw + 1,         cellSize1 - tempStr.jw, tempStr.kw + 1,         tempNumber, tempVel, 4);
        ptrArray_in[tempStr.face][tempStr.ig+1][tempStr.jg+2][tempStr.kg+1]->UpdateDueToWgt( cellSize1 - tempStr.iw, tempStr.jw + 1,         tempStr.kw + 1,         tempNumber, tempVel, 4);
        ptrArray_in[tempStr.face][tempStr.ig+2][tempStr.jg+2][tempStr.kg+1]->UpdateDueToWgt( tempStr.iw + 1,         tempStr.jw + 1,         tempStr.kw + 1,         tempNumber, tempVel, 4);
        }
        else if( tempStr.kg == fieldsGridsSize - 1)
        {
        ptrArray_in[tempStr.face][tempStr.ig+1][tempStr.jg+1][tempStr.kg]->UpdateDueToWgt( cellSize1 - tempStr.iw, cellSize1 - tempStr.jw, cellSize1 - tempStr.kw, tempNumber, tempVel, 4);
        ptrArray_in[tempStr.face][tempStr.ig+2][tempStr.jg+1][tempStr.kg]->UpdateDueToWgt( tempStr.iw + 1,         cellSize1 - tempStr.jw, cellSize1 - tempStr.kw, tempNumber, tempVel, 4);
        ptrArray_in[tempStr.face][tempStr.ig+1][tempStr.jg+2][tempStr.kg]->UpdateDueToWgt( cellSize1 - tempStr.iw, tempStr.jw + 1,         cellSize1 - tempStr.kw, tempNumber, tempVel, 4);
        ptrArray_in[tempStr.face][tempStr.ig+2][tempStr.jg+2][tempStr.kg]->UpdateDueToWgt( tempStr.iw + 1,         tempStr.jw + 1,         cellSize1 - tempStr.kw, tempNumber, tempVel, 4);
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
        double tempNumber = temp.WeightNi();

        Vector3 tempVel = temp.VelParticles();

        if( tempStr.kg == 0)
        {    
        ptrArray_in[tempStr.face][tempStr.ig+1][tempStr.jg+1][tempStr.kg+1]->UpdateDueToWgt( cellSize1 - tempStr.iw, cellSize1 - tempStr.jw, tempStr.kw + 1,         tempNumber, tempVel, 16);
        ptrArray_in[tempStr.face][tempStr.ig+2][tempStr.jg+1][tempStr.kg+1]->UpdateDueToWgt( tempStr.iw + 1,         cellSize1 - tempStr.jw, tempStr.kw + 1,         tempNumber, tempVel, 16);
        ptrArray_in[tempStr.face][tempStr.ig+1][tempStr.jg+2][tempStr.kg+1]->UpdateDueToWgt( cellSize1 - tempStr.iw, tempStr.jw + 1,         tempStr.kw + 1,         tempNumber, tempVel, 16);
        ptrArray_in[tempStr.face][tempStr.ig+2][tempStr.jg+2][tempStr.kg+1]->UpdateDueToWgt( tempStr.iw + 1,         tempStr.jw + 1,         tempStr.kw + 1,         tempNumber, tempVel, 16);
        }
        else if( tempStr.kg == fieldsGridsSize - 1)
        {
        ptrArray_in[tempStr.face][tempStr.ig+1][tempStr.jg+1][tempStr.kg]->UpdateDueToWgt( cellSize1 - tempStr.iw, cellSize1 - tempStr.jw, cellSize1 - tempStr.kw, tempNumber, tempVel, 16);
        ptrArray_in[tempStr.face][tempStr.ig+2][tempStr.jg+1][tempStr.kg]->UpdateDueToWgt( tempStr.iw + 1,         cellSize1 - tempStr.jw, cellSize1 - tempStr.kw, tempNumber, tempVel, 16);
        ptrArray_in[tempStr.face][tempStr.ig+1][tempStr.jg+2][tempStr.kg]->UpdateDueToWgt( cellSize1 - tempStr.iw, tempStr.jw + 1,         cellSize1 - tempStr.kw, tempNumber, tempVel, 16);
        ptrArray_in[tempStr.face][tempStr.ig+2][tempStr.jg+2][tempStr.kg]->UpdateDueToWgt( tempStr.iw + 1,         tempStr.jw + 1,         cellSize1 - tempStr.kw, tempNumber, tempVel, 16);
        }
    }

*/
// finish culmulating and average the density and velocity
    if( timeline_in % updateInfoPeriod_in == 0)
    {
        for( int face = 0; face < totalFace; face++)
        {
            for( int i = 1; i < fieldsGridsSize+2; i++)
            {
                for( int j = 1; j < fieldsGridsSize+2; j++)
                {
                    for( int k = 1 + tempGridsCellLevel; k < fieldsGridsSize - tempGridsCellLevel; k++)
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

#pragma omp parallel for collapse(4)
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
