#ifndef _MODULE_H_
#define _MODULE_H_
#include<iostream>
#include "parameters.h"
#include <memory>
#include "particles.h"
#include <vector>
#include "fieldsgrids.h"
#include <cmath>
#include <limits>
#include <bitset>
using std::shared_ptr;
using std::vector;

//************************************************************************
//************************************************************************
// FUNCTION // Return a authogonal vector to a face when provied a location
// (face, i, j, k). The norm is 1.
//  We want to point out all the six face orthogonal unit vector.
//************************************************************************
//************************************************************************
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
//************************************************************************
//************************************************************************
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
//************************************************************************
//************************************************************************
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
//************************************************************************
//************************************************************************
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
//************************************************************************
//************************************************************************
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
//************************************************************************
//************************************************************************
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
//************************************************************************
//************************************************************************
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
//************************************************************************
//************************************************************************
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
//************************************************************************
//************************************************************************
void updateCellMatrix(Vector3**** curlB_in, Vector3**** curlE_in,
                      Vector3**** gradPe_in, GridsPoints***** ptrArray_in, int face_in);

                      
//************************************************************************
//************************************************************************
// FUNCTION 
// Value gradient field of Pe.
//************************************************************************
//************************************************************************
Vector3*** ValueGradientPe(Vector3*** gradientArray_in, double*** ptrVolumeCellArray_in, GridsPoints***** ptrArray_in, int face_in);


//************************************************************************
//************************************************************************
// FUNCTION 
// Value the matrix field using finite volume method, put in the pointer 
// of the MatrixField, value it, and return the pointer.
// Notice that the cell at corners should be absent in calculation.
//************************************************************************
//************************************************************************
Vector3**** ValueCurlField( Vector3**** curlArray_in, GridsPoints***** ptrArray_in, int face_in, char field_in);

//************************************************************************
//************************************************************************
// FUNCTION rand (0 - 1)
//************************************************************************
//************************************************************************
inline double dRand()
{
    return (double) rand() / (RAND_MAX);
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
//************************************************************************
//************************************************************************
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
//************************************************************************
//************************************************************************
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
//************************************************************************
//************************************************************************
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
//************************************************************************
//************************************************************************
inline Vector3 UniformDisVector3( Vector3 v1, Vector3 v2)
{   
    Vector3 temp = v2.MinusProduct( v1).ScaleProduct( dRand());
    return temp;
}

//************************************************************************
//************************************************************************
// FUNCTION return uniform distrubtion random Vector3
// in cell (face, i, j, k).
// Remind that there are limited position for particles to locate, we just
// need to random the position in smaller cells
//************************************************************************
//************************************************************************
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
// FUNCTION
// Initialization the particles for velocity and position
//************************************************************************
//************************************************************************
vector<Particles>* ParticlesLists( GridsPoints***** ptrArray_in, double*** ptrVolumeCellArray_in);

//************************************************************************
//************************************************************************
// FUNCTION
// Each calculation are on the girds.
// Step 1: Generate a new matrix fulling with gridspoints class
// Step 2: Print out it as .h5
//************************************************************************
//************************************************************************
void PrintOutHdf5( GridsPoints***** ptrArray_in, int i_in, int h5FileCheck_in);

//************************************************************************
//************************************************************************
// FUNCTION
// Create grids
//************************************************************************
//************************************************************************
GridsPoints***** GridsCreation();

//************************************************************************
//************************************************************************
// FUNCTION
// Main general control funtion
//************************************************************************
//************************************************************************
void ProcessFunc();


//************************************************************************
//************************************************************************
// FUNCTION
// Update info in the grids due to the info of particles and related 
// weighting
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
                      int timeline_in, int updateInfoPeriod_in);

//************************************************************************
//************************************************************************
// FUNCTION
// Create a array to store volume of cells needed to calculate density 
// at each grids points
// (fieldGridsSize+1)*(fieldGridsSize+1)*(fieldGridsSize+1) for "double" 
//************************************************************************
//************************************************************************
double*** VolumeGridsField( double*** ptrVolumeCellArray_in);

//************************************************************************
//************************************************************************
// FUNCTION
// Create a array to store volume of cells needed to calculate density 
// at each cell
// (fieldGridsSize * fieldGridsSize * fieldGridsSize) for "double" type 
//************************************************************************
//************************************************************************
double*** VolumeCellsField( GridsPoints***** ptrArray_in);

//************************************************************************
//************************************************************************
// FUNCTION
// Each calculation are on the girds.
// Update E at grids for ve = vi ( no current)
//************************************************************************
//************************************************************************
void updateGrids_nocurrent( Vector3*** gradPe_in, GridsPoints***** ptrArray_in, int face_in);

//************************************************************************
//************************************************************************
// FUNCTION
// Each calculation are on the girds.
//************************************************************************
//************************************************************************
void updateGrids_withcurrent( Vector3*** ptrVectorCellArray_in, GridsPoints***** ptrArray_in, int face_in);

//************************************************************************
//************************************************************************
// FUNCTION
// Calculation the temprature
//************************************************************************
//************************************************************************
void Titheridge_Te(GridsPoints***** ptrArray_in);


//************************************************************************
//************************************************************************
// FUNCTION
// Calculate the gradient of normal of B
//************************************************************************
//************************************************************************
void GradBNorm( GridsPoints***** ptrArray_in);

//************************************************************************
//************************************************************************
// FUNCTION
// Initialization the particles for velocity and position 
// Generate lists of particles for bot and top region temp
//************************************************************************
//************************************************************************
vector<Particles>* ParticlesListsTemp( GridsPoints***** ptrArray_in, double*** ptrVolumeCellArray_in, double mi0, int ionType_in);


//************************************************************************
//************************************************************************
// FUNCTION
// Sec convection velocity due to a ideal convectional cell and a dipole
// magnetic field
//************************************************************************
//************************************************************************
void SetConvectionVel( GridsPoints***** ptrArray_in, int face_in, int i_in, int k_in, int j_in);

//************************************************************************
//************************************************************************
// Function
// Set initial condition
//************************************************************************
//************************************************************************
void SetInitialCondition( GridsPoints***** ptrArray_in, Vector3*** ptrVectorCellArray_in, double*** ptrVolumeCellArray_in);

//************************************************************************
//************************************************************************
// FUNCTION
// Each calculation are on the girds.
// Calculate the grad|B| on the gridspoints
//************************************************************************
//************************************************************************
void UpdateGradBNorm( Vector3*** gradBNorm_in, GridsPoints***** ptrArray_in, int face_in);



//************************************************************************
//************************************************************************
// FUNCTION
// Transform the ubit64 to vector of position
//************************************************************************
//************************************************************************
Vector3 Uint64ToVector3( uint_64 intPos_in);
//************************************************************************
//************************************************************************
// FUNCTION
// Initialization the particles for the magnetic moment / adiabatic invarient
// 
//************************************************************************
//************************************************************************



#endif