#ifndef _MODULE_H_0_
#define _MODULE_H_0_
#include<iostream>
#include "parameters.h"
#include <memory>
#include "particles.h"
#include <vector>
#include "fieldsgrids.h"
#include "module_0.h"
#include <cmath>
#include <limits>
#include <bitset>


//************************************************************************
//************************************************************************
// FUNCTION // Return a authogonal vector to a face when provied a location
// (face, i, j, k). The norm is 1.
//  We want to point out all the six face orthogonal unit vector.

inline Vector3 UnitVectorL(GridsPoints***** ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    return ptrArray_in[face_in][i_in][j_in][k_in]->Pos3().OrthoUnitVector( ptrArray_in[face_in][i_in][j_in+1][k_in]->Pos3());
}

inline Vector3 UnitVectorT(GridsPoints***** ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    return ptrArray_in[face_in][i_in][j_in+1][k_in]->Pos3().OrthoUnitVector( ptrArray_in[face_in][i_in+1][j_in+1][k_in]->Pos3());
}

inline Vector3 UnitVectorR(GridsPoints***** ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    return ptrArray_in[face_in][i_in+1][j_in+1][k_in]->Pos3().OrthoUnitVector( ptrArray_in[face_in][i_in+1][j_in][k_in]->Pos3());
}

inline Vector3 UnitVectorBot(GridsPoints***** ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    return ptrArray_in[face_in][i_in+1][j_in][k_in]->Pos3().OrthoUnitVector( ptrArray_in[face_in][i_in][j_in][k_in]->Pos3());
}

inline Vector3 UnitVectorF(GridsPoints***** ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    Vector3 temp1 = ptrArray_in[face_in][i_in][j_in][k_in]->Pos3().PlusProduct( ptrArray_in[face_in][i_in][j_in+1][k_in]->Pos3());
    Vector3 temp2 = ptrArray_in[face_in][i_in+1][j_in][k_in]->Pos3().PlusProduct( ptrArray_in[face_in][i_in+1][j_in+1][k_in]->Pos3());
    return temp1.PlusProduct(temp2).NormalizedVector();
}

inline Vector3 UnitVectorBack(GridsPoints***** ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    return UnitVectorF(ptrArray_in, face_in, i_in, j_in, k_in).ScaleProduct(-1.0);
}


//************************************************************************
//************************************************************************
// FUNCTION // Return a average field vector on a face when provied a 
// location (face, i, j, k). The norm is of the average of four vector at the 
// nearest gridspoints. 
// We want to point out all the six face field vector for E and B and P

inline Vector3 FaceFieldVectorL(GridsPoints***** ptrArray_in, int face_in, int i_in, int j_in, int k_in, char field_in)
{
    Vector3 temp1, temp2;
    switch(field_in)
    {
    case 'E': {
    temp1 = ptrArray_in[face_in][i_in][j_in][k_in]->E3().PlusProduct( ptrArray_in[face_in][i_in][j_in+1][k_in]->E3());
    temp2 = ptrArray_in[face_in][i_in][j_in][k_in+1]->E3().PlusProduct( ptrArray_in[face_in][i_in][j_in+1][k_in+1]->E3());
    break;}
    case 'B': {
    temp1 = ptrArray_in[face_in][i_in][j_in][k_in]->B3().PlusProduct( ptrArray_in[face_in][i_in][j_in+1][k_in]->B3());
    temp2 = ptrArray_in[face_in][i_in][j_in][k_in+1]->B3().PlusProduct( ptrArray_in[face_in][i_in][j_in+1][k_in+1]->B3());
    break;}
    case 'D': {
    temp1 = ptrArray_in[face_in][i_in][j_in][k_in]->DB3().PlusProduct( ptrArray_in[face_in][i_in][j_in+1][k_in]->DB3());
    temp2 = ptrArray_in[face_in][i_in][j_in][k_in+1]->DB3().PlusProduct( ptrArray_in[face_in][i_in][j_in+1][k_in+1]->DB3());
    break;}
    case 'P': {
    temp1 = ptrArray_in[face_in][i_in][j_in][k_in]->Pos3().PlusProduct( ptrArray_in[face_in][i_in][j_in+1][k_in]->Pos3());
    temp2 = ptrArray_in[face_in][i_in][j_in][k_in+1]->Pos3().PlusProduct( ptrArray_in[face_in][i_in][j_in+1][k_in+1]->Pos3());
    break;}
    }
    return temp1.PlusProduct(temp2).ScaleProduct(0.25);
};

inline Vector3 FaceFieldVectorT(GridsPoints***** ptrArray_in, int face_in, int i_in, int j_in, int k_in, char field_in)
{
    Vector3 temp1, temp2;
    switch(field_in)
    {
    case 'E': {
    temp1 = ptrArray_in[face_in][i_in][j_in+1][k_in]->E3().PlusProduct( ptrArray_in[face_in][i_in+1][j_in+1][k_in]->E3());
    temp2 = ptrArray_in[face_in][i_in][j_in+1][k_in+1]->E3().PlusProduct( ptrArray_in[face_in][i_in+1][j_in+1][k_in+1]->E3());
    break;}
    case 'B': {
    temp1 = ptrArray_in[face_in][i_in][j_in+1][k_in]->B3().PlusProduct( ptrArray_in[face_in][i_in+1][j_in+1][k_in]->B3());
    temp2 = ptrArray_in[face_in][i_in][j_in+1][k_in+1]->B3().PlusProduct( ptrArray_in[face_in][i_in+1][j_in+1][k_in+1]->B3());
    break;}
    case 'P': {
    temp1 = ptrArray_in[face_in][i_in][j_in+1][k_in]->Pos3().PlusProduct( ptrArray_in[face_in][i_in+1][j_in+1][k_in]->Pos3());
    temp2 = ptrArray_in[face_in][i_in][j_in+1][k_in+1]->Pos3().PlusProduct( ptrArray_in[face_in][i_in+1][j_in+1][k_in+1]->Pos3());
    break;}
    }
    return temp1.PlusProduct(temp2).ScaleProduct(0.25);    
};

inline Vector3 FaceFieldVectorR(GridsPoints***** ptrArray_in, int face_in, int i_in, int j_in, int k_in, char field_in)
{
    Vector3 temp1, temp2;
    switch(field_in)
    {
    case 'E': { 
    temp1 = ptrArray_in[face_in][i_in+1][j_in][k_in]->E3().PlusProduct( ptrArray_in[face_in][i_in+1][j_in+1][k_in]->E3());
    temp2 = ptrArray_in[face_in][i_in+1][j_in][k_in+1]->E3().PlusProduct( ptrArray_in[face_in][i_in+1][j_in+1][k_in+1]->E3());
    break;}
    case 'B': {
    temp1 = ptrArray_in[face_in][i_in+1][j_in][k_in]->B3().PlusProduct( ptrArray_in[face_in][i_in+1][j_in+1][k_in]->B3());
    temp2 = ptrArray_in[face_in][i_in+1][j_in][k_in+1]->B3().PlusProduct( ptrArray_in[face_in][i_in+1][j_in+1][k_in+1]->B3());
    break;}
    case 'P': {
    temp1 = ptrArray_in[face_in][i_in+1][j_in][k_in]->Pos3().PlusProduct( ptrArray_in[face_in][i_in+1][j_in+1][k_in]->Pos3());
    temp2 = ptrArray_in[face_in][i_in+1][j_in][k_in+1]->Pos3().PlusProduct( ptrArray_in[face_in][i_in+1][j_in+1][k_in+1]->Pos3());
    break;}
    }
    return temp1.PlusProduct(temp2).ScaleProduct(0.25);
};

inline Vector3 FaceFieldVectorBot(GridsPoints***** ptrArray_in, int face_in, int i_in, int j_in, int k_in, char field_in)
{
    Vector3 temp1, temp2;
    switch(field_in)
    {
    case 'E': { 
    temp1 = ptrArray_in[face_in][i_in][j_in][k_in]->E3().PlusProduct( ptrArray_in[face_in][i_in+1][j_in][k_in]->E3());
    temp2 = ptrArray_in[face_in][i_in][j_in][k_in+1]->E3().PlusProduct( ptrArray_in[face_in][i_in+1][j_in][k_in+1]->E3());
    break;}
    case 'B': {
    temp1 = ptrArray_in[face_in][i_in][j_in][k_in]->B3().PlusProduct( ptrArray_in[face_in][i_in+1][j_in+1][k_in]->B3());
    temp2 = ptrArray_in[face_in][i_in][j_in][k_in+1]->B3().PlusProduct( ptrArray_in[face_in][i_in+1][j_in][k_in+1]->B3());
    break;}
    case 'P': {
    temp1 = ptrArray_in[face_in][i_in][j_in][k_in]->Pos3().PlusProduct( ptrArray_in[face_in][i_in+1][j_in+1][k_in]->Pos3());
    temp2 = ptrArray_in[face_in][i_in][j_in][k_in+1]->Pos3().PlusProduct( ptrArray_in[face_in][i_in+1][j_in][k_in+1]->Pos3());
    break;}
    }
    return temp1.PlusProduct(temp2).ScaleProduct(0.25); // 1/4 is not acceptable !
};

inline Vector3 FaceFieldVectorF(GridsPoints***** ptrArray_in, int face_in, int i_in, int j_in, int k_in, char field_in)
{
    Vector3 temp1, temp2;
    switch(field_in)
    {
    case 'E': { 
    temp1 = ptrArray_in[face_in][i_in][j_in][k_in+1]->E3().PlusProduct( ptrArray_in[face_in][i_in][j_in+1][k_in+1]->E3());
    temp2 = ptrArray_in[face_in][i_in+1][j_in][k_in+1]->E3().PlusProduct( ptrArray_in[face_in][i_in+1][j_in+1][k_in+1]->E3());
    break;}
    case 'B': {
    temp1 = ptrArray_in[face_in][i_in][j_in][k_in+1]->B3().PlusProduct( ptrArray_in[face_in][i_in][j_in+1][k_in+1]->B3());
    temp2 = ptrArray_in[face_in][i_in+1][j_in][k_in+1]->B3().PlusProduct( ptrArray_in[face_in][i_in+1][j_in+1][k_in+1]->B3());
    break;}
    case 'P': {
    temp1 = ptrArray_in[face_in][i_in][j_in][k_in+1]->Pos3().PlusProduct( ptrArray_in[face_in][i_in][j_in+1][k_in+1]->Pos3());
    temp2 = ptrArray_in[face_in][i_in+1][j_in][k_in+1]->Pos3().PlusProduct( ptrArray_in[face_in][i_in+1][j_in+1][k_in+1]->Pos3());
    break;}
    }
//    std::cout << " PPP" << temp1.PlusProduct(temp2).ScaleProduct(0.25).x() << " " << temp1.PlusProduct(temp2).y() << " " << temp1.PlusProduct(temp2).z() << " " ;
  
    return temp1.PlusProduct(temp2).ScaleProduct(0.25);
};

inline Vector3 FaceFieldVectorBack(GridsPoints***** ptrArray_in, int face_in, int i_in, int j_in, int k_in, char field_in)
{
    Vector3 temp1, temp2;
    switch(field_in)
    {
    case 'E': { 
    temp1 = ptrArray_in[face_in][i_in][j_in][k_in]->E3().PlusProduct( ptrArray_in[face_in][i_in][j_in+1][k_in]->E3());
    temp2 = ptrArray_in[face_in][i_in+1][j_in][k_in]->E3().PlusProduct( ptrArray_in[face_in][i_in+1][j_in+1][k_in]->E3());
    break;}
    case 'B': {
    temp1 = ptrArray_in[face_in][i_in][j_in][k_in]->B3().PlusProduct( ptrArray_in[face_in][i_in][j_in+1][k_in]->B3());
    temp2 = ptrArray_in[face_in][i_in+1][j_in][k_in]->B3().PlusProduct( ptrArray_in[face_in][i_in+1][j_in+1][k_in]->B3());
    break;}
    case 'P': {
    temp1 = ptrArray_in[face_in][i_in][j_in][k_in]->Pos3().PlusProduct( ptrArray_in[face_in][i_in][j_in+1][k_in]->Pos3());
    temp2 = ptrArray_in[face_in][i_in+1][j_in][k_in]->Pos3().PlusProduct( ptrArray_in[face_in][i_in+1][j_in+1][k_in]->Pos3());
    break;}
    }
    return temp1.PlusProduct(temp2).ScaleProduct(0.25);
};


//************************************************************************
//************************************************************************
// FUNCTION 
// Return a average density of ions( electrons) on a face when provied a 
// location (face, i, j, k). The norm is of the average of four densities
// at the nearest gridspoints. 
// We want to point out all the six face field vector for E and B and P

inline double FaceDensityL(GridsPoints***** ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    double temp1, temp2;
    temp1 = ptrArray_in[face_in][i_in][j_in][k_in]->Density() + ptrArray_in[face_in][i_in][j_in+1][k_in]->Density();
    temp2 = ptrArray_in[face_in][i_in][j_in][k_in+1]->Density() + ptrArray_in[face_in][i_in][j_in+1][k_in+1]->Density();

    return (temp1 + temp2) * (1.0/4.0);
};

inline double FaceDensityT(GridsPoints***** ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    double temp1, temp2;
    temp1 = ptrArray_in[face_in][i_in][j_in+1][k_in]->Density() + ptrArray_in[face_in][i_in+1][j_in+1][k_in]->Density();
    temp2 = ptrArray_in[face_in][i_in][j_in+1][k_in+1]->Density() + ptrArray_in[face_in][i_in+1][j_in+1][k_in+1]->Density();
    return (temp1 + temp2) * (1.0/4.0);
};

inline double FaceDensityR(GridsPoints***** ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    double temp1, temp2;
    temp1 = ptrArray_in[face_in][i_in+1][j_in][k_in]->Density() +  ptrArray_in[face_in][i_in+1][j_in+1][k_in]->Density();
    temp2 = ptrArray_in[face_in][i_in+1][j_in][k_in+1]->Density() + ptrArray_in[face_in][i_in+1][j_in+1][k_in+1]->Density();
    return (temp1 + temp2) * (1.0/4.0);
};

inline double FaceDensityBot(GridsPoints***** ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    double temp1, temp2;
    temp1 = ptrArray_in[face_in][i_in][j_in][k_in]->Density() + ptrArray_in[face_in][i_in+1][j_in][k_in]->Density();
    temp2 = ptrArray_in[face_in][i_in][j_in][k_in+1]->Density() + ptrArray_in[face_in][i_in+1][j_in][k_in+1]->Density();
    return (temp1 + temp2) * (1.0/4.0);
};

inline double FaceDensityF(GridsPoints***** ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    double temp1, temp2;
    temp1 = ptrArray_in[face_in][i_in][j_in][k_in+1]->Density() + ptrArray_in[face_in][i_in][j_in+1][k_in+1]->Density();
    temp2 = ptrArray_in[face_in][i_in+1][j_in][k_in+1]->Density() + ptrArray_in[face_in][i_in+1][j_in+1][k_in+1]->Density();
    return (temp1 + temp2) * (1.0/4.0);
};

inline double FaceDensityBack(GridsPoints***** ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    double temp1, temp2;
    temp1 = ptrArray_in[face_in][i_in][j_in][k_in]->Density() + ptrArray_in[face_in][i_in][j_in+1][k_in]->Density();
    temp2 = ptrArray_in[face_in][i_in+1][j_in][k_in]->Density() + ptrArray_in[face_in][i_in+1][j_in+1][k_in]->Density();
    return (temp1 + temp2) * (1.0/4.0);
};


//************************************************************************
//************************************************************************
// FUNCTION 
// Return a average number density of ions( electrons) on a face when provied a 
// location (face, i, j, k). The norm is of the average of four densities
// at the nearest gridspoints. 

inline double FaceNumberDensityL(GridsPoints***** ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    double temp1, temp2;
    temp1 = ptrArray_in[face_in][i_in][j_in][k_in]->Density_H() + ptrArray_in[face_in][i_in][j_in+1][k_in]->Density_H() +
            ptrArray_in[face_in][i_in][j_in][k_in]->Density_He() + ptrArray_in[face_in][i_in][j_in+1][k_in]->Density_He() +
            ptrArray_in[face_in][i_in][j_in][k_in]->Density_O() + ptrArray_in[face_in][i_in][j_in+1][k_in]->Density_O();
    temp2 = ptrArray_in[face_in][i_in][j_in][k_in+1]->Density_H() + ptrArray_in[face_in][i_in][j_in+1][k_in+1]->Density_H() +
            ptrArray_in[face_in][i_in][j_in][k_in+1]->Density_He() + ptrArray_in[face_in][i_in][j_in+1][k_in+1]->Density_He() +
            ptrArray_in[face_in][i_in][j_in][k_in+1]->Density_O() + ptrArray_in[face_in][i_in][j_in+1][k_in+1]->Density_O();
    return (temp1 + temp2) * (1.0/4.0);
};

inline double FaceNumberDensityT(GridsPoints***** ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    double temp1, temp2;
    temp1 = ptrArray_in[face_in][i_in][j_in+1][k_in]->Density_H() + ptrArray_in[face_in][i_in+1][j_in+1][k_in]->Density_H() +
            ptrArray_in[face_in][i_in][j_in+1][k_in]->Density_He() + ptrArray_in[face_in][i_in+1][j_in+1][k_in]->Density_He() +
            ptrArray_in[face_in][i_in][j_in+1][k_in]->Density_O() + ptrArray_in[face_in][i_in+1][j_in+1][k_in]->Density_O();
    temp2 = ptrArray_in[face_in][i_in][j_in+1][k_in+1]->Density_H() + ptrArray_in[face_in][i_in+1][j_in+1][k_in+1]->Density_H() +
            ptrArray_in[face_in][i_in][j_in+1][k_in+1]->Density_He() + ptrArray_in[face_in][i_in+1][j_in+1][k_in+1]->Density_He() +
            ptrArray_in[face_in][i_in][j_in+1][k_in+1]->Density_O() + ptrArray_in[face_in][i_in+1][j_in+1][k_in+1]->Density_O();
    return (temp1 + temp2) * (1.0/4.0);
};

inline double FaceNumberDensityR(GridsPoints***** ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    double temp1, temp2;
    temp1 = ptrArray_in[face_in][i_in+1][j_in][k_in]->Density_H() +  ptrArray_in[face_in][i_in+1][j_in+1][k_in]->Density_H() + 
            ptrArray_in[face_in][i_in+1][j_in][k_in]->Density_He() +  ptrArray_in[face_in][i_in+1][j_in+1][k_in]->Density_He() +
            ptrArray_in[face_in][i_in+1][j_in][k_in]->Density_O() +  ptrArray_in[face_in][i_in+1][j_in+1][k_in]->Density_O();
    
    temp2 = ptrArray_in[face_in][i_in+1][j_in][k_in+1]->Density_H() + ptrArray_in[face_in][i_in+1][j_in+1][k_in+1]->Density_H() +
            ptrArray_in[face_in][i_in+1][j_in][k_in+1]->Density_He() + ptrArray_in[face_in][i_in+1][j_in+1][k_in+1]->Density_He() +
            ptrArray_in[face_in][i_in+1][j_in][k_in+1]->Density_O() + ptrArray_in[face_in][i_in+1][j_in+1][k_in+1]->Density_O();
    return (temp1 + temp2) * (1.0/4.0);
};

inline double FaceNumberDensityBot(GridsPoints***** ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    double temp1, temp2;
    temp1 = ptrArray_in[face_in][i_in][j_in][k_in]->Density_H() + ptrArray_in[face_in][i_in+1][j_in][k_in]->Density_H() +
            ptrArray_in[face_in][i_in][j_in][k_in]->Density_He() + ptrArray_in[face_in][i_in+1][j_in][k_in]->Density_He() +
            ptrArray_in[face_in][i_in][j_in][k_in]->Density_O() + ptrArray_in[face_in][i_in+1][j_in][k_in]->Density_O();

    temp2 = ptrArray_in[face_in][i_in][j_in][k_in+1]->Density_H() + ptrArray_in[face_in][i_in+1][j_in][k_in+1]->Density_H() + 
            ptrArray_in[face_in][i_in][j_in][k_in+1]->Density_He() + ptrArray_in[face_in][i_in+1][j_in][k_in+1]->Density_He() +
            ptrArray_in[face_in][i_in][j_in][k_in+1]->Density_O() + ptrArray_in[face_in][i_in+1][j_in][k_in+1]->Density_O();
    return (temp1 + temp2) * (1.0/4.0);
};

inline double FaceNumberDensityF(GridsPoints***** ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    double temp1, temp2;
    temp1 = ptrArray_in[face_in][i_in][j_in][k_in+1]->Density_H() + ptrArray_in[face_in][i_in][j_in+1][k_in+1]->Density_H() + 
            ptrArray_in[face_in][i_in][j_in][k_in+1]->Density_He() + ptrArray_in[face_in][i_in][j_in+1][k_in+1]->Density_He() + 
            ptrArray_in[face_in][i_in][j_in][k_in+1]->Density_O() + ptrArray_in[face_in][i_in][j_in+1][k_in+1]->Density_O();
    temp2 = ptrArray_in[face_in][i_in+1][j_in][k_in+1]->Density_H() + ptrArray_in[face_in][i_in+1][j_in+1][k_in+1]->Density_H() +
            ptrArray_in[face_in][i_in+1][j_in][k_in+1]->Density_He() + ptrArray_in[face_in][i_in+1][j_in+1][k_in+1]->Density_He() +
            ptrArray_in[face_in][i_in+1][j_in][k_in+1]->Density_O() + ptrArray_in[face_in][i_in+1][j_in+1][k_in+1]->Density_O();
    return (temp1 + temp2) * (1.0/4.0);
};

inline double FaceNumberDensityBack(GridsPoints***** ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    double temp1, temp2;
    temp1 = ptrArray_in[face_in][i_in][j_in][k_in]->Density_H() + ptrArray_in[face_in][i_in][j_in+1][k_in]->Density_H() +
            ptrArray_in[face_in][i_in][j_in][k_in]->Density_He() + ptrArray_in[face_in][i_in][j_in+1][k_in]->Density_He() +
            ptrArray_in[face_in][i_in][j_in][k_in]->Density_O() + ptrArray_in[face_in][i_in][j_in+1][k_in]->Density_O();
    temp2 = ptrArray_in[face_in][i_in+1][j_in][k_in]->Density_H() + ptrArray_in[face_in][i_in+1][j_in+1][k_in]->Density_H() +
            ptrArray_in[face_in][i_in+1][j_in][k_in]->Density_He() + ptrArray_in[face_in][i_in+1][j_in+1][k_in]->Density_He() +
            ptrArray_in[face_in][i_in+1][j_in][k_in]->Density_O() + ptrArray_in[face_in][i_in+1][j_in+1][k_in]->Density_O() ;
    return (temp1 + temp2) * (1.0/4.0);
};


//************************************************************************
//************************************************************************
// FUNCTION 
// Return a average Temperature of electrons on a face when provied a 
// location (face, i, j, k). The norm is of the average of four temperature
// at the nearest gridspoints. 
// We want to point out all the six face field vector for E and B and P

inline double FaceTemperatureL(GridsPoints***** ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    double temp1, temp2;
    temp1 = ptrArray_in[face_in][i_in][j_in][k_in]->Temperature() + ptrArray_in[face_in][i_in][j_in+1][k_in]->Temperature();
    temp2 = ptrArray_in[face_in][i_in][j_in][k_in+1]->Temperature() + ptrArray_in[face_in][i_in][j_in+1][k_in+1]->Temperature();

    return (temp1 + temp2) * (1.0/4.0);
};

inline double FaceTemperatureT(GridsPoints***** ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    double temp1, temp2;
    temp1 = ptrArray_in[face_in][i_in][j_in+1][k_in]->Temperature() + ptrArray_in[face_in][i_in+1][j_in+1][k_in]->Temperature();
    temp2 = ptrArray_in[face_in][i_in][j_in+1][k_in+1]->Temperature() + ptrArray_in[face_in][i_in+1][j_in+1][k_in+1]->Temperature();
    return (temp1 + temp2) * (1.0/4.0);
};

inline double FaceTemperatureR(GridsPoints***** ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    double temp1, temp2;
    temp1 = ptrArray_in[face_in][i_in+1][j_in][k_in]->Temperature() +  ptrArray_in[face_in][i_in+1][j_in+1][k_in]->Temperature();
    temp2 = ptrArray_in[face_in][i_in+1][j_in][k_in+1]->Temperature() + ptrArray_in[face_in][i_in+1][j_in+1][k_in+1]->Temperature();
    return (temp1 + temp2) * (1.0/4.0);
};

inline double FaceTemperatureBot(GridsPoints***** ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    double temp1, temp2;
    temp1 = ptrArray_in[face_in][i_in][j_in][k_in]->Temperature() + ptrArray_in[face_in][i_in+1][j_in][k_in]->Temperature();
    temp2 = ptrArray_in[face_in][i_in][j_in][k_in+1]->Temperature() + ptrArray_in[face_in][i_in+1][j_in][k_in+1]->Temperature();
    return (temp1 + temp2) * (1.0/4.0);
};

inline double FaceTemperatureF(GridsPoints***** ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    double temp1, temp2;
    temp1 = ptrArray_in[face_in][i_in][j_in][k_in+1]->Temperature() + ptrArray_in[face_in][i_in][j_in+1][k_in+1]->Temperature();
    temp2 = ptrArray_in[face_in][i_in+1][j_in][k_in+1]->Temperature() + ptrArray_in[face_in][i_in+1][j_in+1][k_in+1]->Temperature();
    return (temp1 + temp2) * (1.0/4.0);
};

inline double FaceTemperatureBack(GridsPoints***** ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    double temp1, temp2;
    temp1 = ptrArray_in[face_in][i_in][j_in][k_in]->Temperature() + ptrArray_in[face_in][i_in][j_in+1][k_in]->Temperature();
    temp2 = ptrArray_in[face_in][i_in+1][j_in][k_in]->Temperature() + ptrArray_in[face_in][i_in+1][j_in+1][k_in]->Temperature();
    return (temp1 + temp2) * (1.0/4.0);
};


//************************************************************************
//************************************************************************
// FUNCTION 
// Return a average>B3().norm(of electrons on a face when provied a 
// location (face, i, j, k). The norm is of the average of four temperature
// at the nearest gridspoints. 
// We want to point out all the six face field vector for E and B and P

inline double FaceNormBL(GridsPoints***** ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    double temp1, temp2;
    temp1 = ptrArray_in[face_in][i_in][j_in][k_in]->B3().norm() + ptrArray_in[face_in][i_in][j_in+1][k_in]->B3().norm();
    temp2 = ptrArray_in[face_in][i_in][j_in][k_in+1]->B3().norm() + ptrArray_in[face_in][i_in][j_in+1][k_in+1]->B3().norm();

    return (temp1 + temp2) * (1.0/4.0);
};

inline double FaceNormBT(GridsPoints***** ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    double temp1, temp2;
    temp1 = ptrArray_in[face_in][i_in][j_in+1][k_in]->B3().norm() + ptrArray_in[face_in][i_in+1][j_in+1][k_in]->B3().norm();
    temp2 = ptrArray_in[face_in][i_in][j_in+1][k_in+1]->B3().norm() + ptrArray_in[face_in][i_in+1][j_in+1][k_in+1]->B3().norm();
    return (temp1 + temp2) * (1.0/4.0);
};

inline double FaceNormBR(GridsPoints***** ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    double temp1, temp2;
    temp1 = ptrArray_in[face_in][i_in+1][j_in][k_in]->B3().norm() +  ptrArray_in[face_in][i_in+1][j_in+1][k_in]->B3().norm();
    temp2 = ptrArray_in[face_in][i_in+1][j_in][k_in+1]->B3().norm() + ptrArray_in[face_in][i_in+1][j_in+1][k_in+1]->B3().norm();
    return (temp1 + temp2) * (1.0/4.0);
};

inline double FaceNormBBot(GridsPoints***** ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    double temp1, temp2;
    temp1 = ptrArray_in[face_in][i_in][j_in][k_in]->B3().norm() + ptrArray_in[face_in][i_in+1][j_in][k_in]->B3().norm();
    temp2 = ptrArray_in[face_in][i_in][j_in][k_in+1]->B3().norm() + ptrArray_in[face_in][i_in+1][j_in][k_in+1]->B3().norm();
    return (temp1 + temp2) * (1.0/4.0);
};

inline double FaceNormBF(GridsPoints***** ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    double temp1, temp2;
    temp1 = ptrArray_in[face_in][i_in][j_in][k_in+1]->B3().norm() + ptrArray_in[face_in][i_in][j_in+1][k_in+1]->B3().norm();
    temp2 = ptrArray_in[face_in][i_in+1][j_in][k_in+1]->B3().norm() + ptrArray_in[face_in][i_in+1][j_in+1][k_in+1]->B3().norm();
    return (temp1 + temp2) * (1.0/4.0);
};

inline double FaceNormBBack(GridsPoints***** ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    double temp1, temp2;
    temp1 = ptrArray_in[face_in][i_in][j_in][k_in]->B3().norm() + ptrArray_in[face_in][i_in][j_in+1][k_in]->B3().norm();
    temp2 = ptrArray_in[face_in][i_in+1][j_in][k_in]->B3().norm() + ptrArray_in[face_in][i_in+1][j_in+1][k_in]->B3().norm();
    return (temp1 + temp2) * (1.0/4.0);
};


//************************************************************************
//************************************************************************
// FUNCTION // Return a vector athogonal to a face when provied a 
// location (face, i, j, k). 
// We want to point out all the six vectors whose norm is the face area.
// For the side area, use bigger tri-angle minus smaller tri-angle to represent 
// the norm.
// For the front and back area, use two small tri-angle sum to represent the 
// norm, and use average vector, as UnitVectorF and UnitVectorBack do, to 
// represent the direction.

inline Vector3 AreaVectorL(GridsPoints***** ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    Vector3 temp1 = ptrArray_in[face_in][i_in][j_in][k_in]->Pos3().CrossProduct( ptrArray_in[face_in][i_in][j_in+1][k_in]->Pos3());
    Vector3 temp2 = ptrArray_in[face_in][i_in][j_in][k_in+1]->Pos3().CrossProduct( ptrArray_in[face_in][i_in][j_in+1][k_in+1]->Pos3());
    return temp2.MinusProduct(temp1).ScaleProduct(0.5);
}

inline Vector3 AreaVectorT(GridsPoints***** ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    Vector3 temp1 = ptrArray_in[face_in][i_in][j_in+1][k_in]->Pos3().CrossProduct( ptrArray_in[face_in][i_in+1][j_in+1][k_in]->Pos3());
    Vector3 temp2 = ptrArray_in[face_in][i_in][j_in+1][k_in+1]->Pos3().CrossProduct( ptrArray_in[face_in][i_in+1][j_in+1][k_in+1]->Pos3());
    return temp2.MinusProduct(temp1).ScaleProduct(0.5);
};

inline Vector3 AreaVectorR(GridsPoints***** ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    Vector3 temp1 = ptrArray_in[face_in][i_in+1][j_in+1][k_in]->Pos3().CrossProduct( ptrArray_in[face_in][i_in+1][j_in][k_in]->Pos3());
    Vector3 temp2 = ptrArray_in[face_in][i_in+1][j_in+1][k_in+1]->Pos3().CrossProduct( ptrArray_in[face_in][i_in+1][j_in][k_in+1]->Pos3());
    return temp2.MinusProduct(temp1).ScaleProduct(0.5);
};

inline Vector3 AreaVectorBot(GridsPoints***** ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    Vector3 temp1 = ptrArray_in[face_in][i_in+1][j_in][k_in]->Pos3().CrossProduct( ptrArray_in[face_in][i_in][j_in][k_in]->Pos3());
    Vector3 temp2 = ptrArray_in[face_in][i_in+1][j_in][k_in+1]->Pos3().CrossProduct( ptrArray_in[face_in][i_in][j_in][k_in+1]->Pos3());
    return temp2.MinusProduct(temp1).ScaleProduct(0.5);
};

inline Vector3 AreaVectorF(GridsPoints***** ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    Vector3 temp1 = ptrArray_in[face_in][i_in][j_in][k_in+1]->Pos3().MinusProduct( ptrArray_in[face_in][i_in+1][j_in][k_in+1]->Pos3());
    Vector3 temp2 = ptrArray_in[face_in][i_in+1][j_in+1][k_in+1]->Pos3().MinusProduct( ptrArray_in[face_in][i_in+1][j_in][k_in+1]->Pos3());
    
    Vector3 temp3 = ptrArray_in[face_in][i_in][j_in][k_in+1]->Pos3().MinusProduct( ptrArray_in[face_in][i_in][j_in+1][k_in+1]->Pos3());
    Vector3 temp4 = ptrArray_in[face_in][i_in+1][j_in+1][k_in+1]->Pos3().MinusProduct( ptrArray_in[face_in][i_in][j_in+1][k_in+1]->Pos3());
/*    std::cout << " check " << ptrArray_in[face_in][i_in][j_in][k_in+1]->Pos3().x() << 
                      " - " << ptrArray_in[face_in][i_in+1][j_in+1][k_in+1]->Pos3().x() << 
                      " - " << ptrArray_in[face_in][i_in][j_in][k_in+1]->Pos3().x() << 
                      " - " << ptrArray_in[face_in][i_in+1][j_in+1][k_in+1]->Pos3().x() << " - ";
*/    return temp2.CrossProduct(temp1).PlusProduct(temp3.CrossProduct(temp4));
};

inline Vector3 AreaVectorBack(GridsPoints***** ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    Vector3 temp1 = ptrArray_in[face_in][i_in][j_in][k_in]->Pos3().MinusProduct( ptrArray_in[face_in][i_in+1][j_in][k_in]->Pos3());
    Vector3 temp2 = ptrArray_in[face_in][i_in+1][j_in+1][k_in]->Pos3().MinusProduct( ptrArray_in[face_in][i_in+1][j_in][k_in]->Pos3());
    
    Vector3 temp3 = ptrArray_in[face_in][i_in][j_in][k_in]->Pos3().MinusProduct( ptrArray_in[face_in][i_in][j_in+1][k_in]->Pos3());
    Vector3 temp4 = ptrArray_in[face_in][i_in+1][j_in+1][k_in]->Pos3().MinusProduct( ptrArray_in[face_in][i_in][j_in+1][k_in]->Pos3());
    
    return temp2.CrossProduct(temp1).PlusProduct(temp3.CrossProduct(temp4)).ScaleProduct(-1.0);
};


//************************************************************************
//************************************************************************
// FUNCTION // Return a volume for a cell. Input (face, i, j, k) and return 
// a double. The value is (average of back and front area) * (norm of vector
// between back and front face vector)

inline double CellVolume(GridsPoints***** ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    double normV = FaceFieldVectorF(ptrArray_in, face_in, i_in, j_in, k_in, 'P').PlusProduct(
                     FaceFieldVectorBack(ptrArray_in, face_in, i_in, j_in, k_in, 'P')).norm()/2.0;
   
    double avgArea = AreaVectorBack(ptrArray_in, face_in, i_in, j_in, k_in).PlusProduct(
                     AreaVectorF(ptrArray_in, face_in, i_in, j_in, k_in).ScaleProduct(-1.0)).norm()/2.0;
/*    std::cout << " test" << AreaVectorF(ptrArray_in, face_in, i_in, j_in, k_in).x() << " " 
                << AreaVectorF(ptrArray_in, face_in, i_in, j_in, k_in).y() << " "
                << AreaVectorF(ptrArray_in, face_in, i_in, j_in, k_in).z() << " "
                << normV * avgArea;
    int pause;
    std::cin >> pause;
*/    return normV * avgArea;
}


//************************************************************************
//************************************************************************
// FUNCTION
// As in the updating curlField and gradientPe array, some variables are
// repeating calculating, it is suitable to put them in one function.
// Therefore, we need three matrix of curlB, curlE, and gradientPe. 
// Assume they are curlB, curlE and gradPe, respectively.

void updateCellMatrix(Vector3**** curlB_in, Vector3**** curlE_in,
                      Vector3**** gradPe_in, GridsPoints***** ptrArray_in, int face_in);

                      
//************************************************************************
//************************************************************************
// FUNCTION 
// Value gradient field of Pe.

Vector3*** ValueGradientPe(Vector3*** gradientArray_in, double*** ptrVolumeCellArray_in, GridsPoints***** ptrArray_in, int face_in);

//************************************************************************
//************************************************************************
// FUNCTION 
// Value gradient field of Pe.
// gradientArray_in is in size of ( fsize+2 * fsize+2 * fsize) with vector3
// ptrVolumeCellArray is in size of ( fsize+2 * fsize+2 * fsize) with double
// Pe = n k T, in which n is the number density, k is the boltzmann constant, and T is the Te
Vector3*** ValueGradient(Vector3*** gradientArray_in, double*** ptrVolumeCellArray_in, GridsPoints***** ptrArray_in, int face_in, char char_in);


//************************************************************************
//************************************************************************
// FUNCTION 
// Value the matrix field using finite volume method, put in the pointer 
// of the MatrixField, value it, and return the pointer.
// Notice that the cell at corners should be absent in calculation.
Vector3*** ValueCurlField( Vector3*** curlArray_in, double*** ptrVolumeCellArray_in, GridsPoints***** ptrArray_in, int face_in, char field_in);


//************************************************************************
//************************************************************************
// FUNCTION
// Create a array to store volume of cells needed to calculate density 
// at each grids points
// (fieldGridsSize+1)*(fieldGridsSize+1)*(fieldGridsSize+1) for "double" 

double*** VolumeGridsField( double*** ptrVolumeCellArray_in);

//************************************************************************
//************************************************************************
// FUNCTION
// Create a array to store volume of cells needed to calculate density 
// at each cell
// (fieldGridsSize * fieldGridsSize * fieldGridsSize) for "double" type 

double*** VolumeCellsField( GridsPoints***** ptrArray_in);

//************************************************************************
//************************************************************************
// FUNCTION
// Each calculation are on the girds.
// Update E at grids for ve = vi ( no current)

void updateGrids_nocurrent( Vector3*** gradPe_in, GridsPoints***** ptrArray_in, int face_in);

//************************************************************************
//************************************************************************
// FUNCTION 
// UpdateVe3
void UpdateVe3( Vector3*** curlField_in, GridsPoints***** ptrArray_in, int face_in);

//************************************************************************
//************************************************************************
// FUNCTION
// Each calculation are on the girds.
// Update E at grids for ve = vi ( no current)
// Update E at grids for ve ( with current)
void UpdateE3( Vector3*** gradPe_in, GridsPoints***** ptrArray_in, int face_in);

//************************************************************************
//************************************************************************
// FUNCTION
// UpdateB3 vased on faraday's law

void UpdateB3( Vector3*** curlField_in, GridsPoints***** ptrArray_in, int face_in);

//************************************************************************
//************************************************************************
// FUNCTION
// Each calculation are on the girds.

void updateGrids_withcurrent( Vector3*** ptrVectorCellArray_in, GridsPoints***** ptrArray_in, int face_in);

//************************************************************************
//************************************************************************
// FUNCTION
// Calculation the temprature

void Titheridge_Te(GridsPoints***** ptrArray_in);


//************************************************************************
//************************************************************************
// FUNCTION
// Calculate the gradient of normal of B

void GradBNorm( GridsPoints***** ptrArray_in);


//************************************************************************
//************************************************************************
// FUNCTION
// Sec convection velocity due to a ideal convectional cell and a dipole
// magnetic field

void SetConvectionVel( GridsPoints***** ptrArray_in, int face_in, int i_in, int k_in, int j_in);

//************************************************************************
//************************************************************************
// Function
// Set initial condition

void SetInitialCondition( GridsPoints***** ptrArray_in, Vector3*** ptrVectorCellArray_in, double*** ptrVolumeCellArray_in);

//************************************************************************
//************************************************************************
// FUNCTION
// Each calculation are on the girds.
// Calculate the grad|B| on the gridspoints

void UpdateGradBNorm( Vector3*** gradBNorm_in, GridsPoints***** ptrArray_in, int face_in);


//************************************************************************
//************************************************************************
// FUNCTION
// Each calculation are on the girds.
// Step 1: Generate a new matrix fulling with gridspoints class
// Step 2: Print out it as .h5
void PrintOutHdf5( GridsPoints***** ptrArray_in, int i_in, int h5FileCheck_in);

//************************************************************************
//************************************************************************
// FUNCTION
// Create grids

GridsPoints***** GridsCreation();

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
Vector3*** VectorCellField();

//************************************************************************
//************************************************************************
// FUNCTION
// finish culmulating and average the density and velocity
void CalculatingAveragedPhoVatGrids(GridsPoints***** ptrArray_in, 
                                    double*** ptrVolumeGridArray_in,
                                    int updateInfoPeriod_in);
                                    
//************************************************************************
//************************************************************************
// FUNCTION
// Set zero for pho and v at each points
void  ResetPhoVatGrids( GridsPoints***** ptrArray_in);

#endif