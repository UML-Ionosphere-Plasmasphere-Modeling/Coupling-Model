//**********************************
// IMPLEMENTATION FILE (spheregrids.cpp)
// This file implements the spheregrids member functions
//**********************************

#include "spheregrids.h"
#include <iostream>
#include <cmath>
using namespace std;

// Private members of class:
// double cx,cy,cz;
// SphereGrid* neighbor[6];
// int k;

//******************************************************
//constructor
SphereGrids::SphereGrids(/* in */ double inicor[3])
{
    for (int i = 0; i < 3; i++)
        cor[i] = inicor[i];
}

//******************************************************
//constructor
SphereGrids::SphereGrids()
{
    cor[0] = 0; cor[1] = 0; cor[2] = 0;
}
//*****************************************************
void SphereGrids::SetCorrdinate( /* in */ double inicor[3])
{
    for (int i = 0; i < 3; i++)
        cor[i] = inicor[i];
}
//******************************************************
void SphereGrids::Write() const
{
    cout << cor[0] << " " << cor[1] << " " << cor[2] << endl;
}
//******************************************************
double SphereGrids::cx() const
{
    double temp=cor[0];
    return temp;
}
double SphereGrids::cy() const
{
    double temp=cor[1];
    return temp;
}
double SphereGrids::cz() const
{
    double temp=cor[2];
    return temp;
}
//******************************************************
void SphereGrids::CenterProjectionOnSphere()
{
    double ratio = sqrt( 1.0 / (pow(cor[0],2) +
                                pow(cor[1],2) + 
                                pow(cor[2],2)));
    cor[0] *= ratio;
    cor[1] *= ratio;
    cor[2] *= ratio;
}
//******************************************************
double SphereGrids::DistanceBetweenPoints( /* in */ SphereGrids otherPoint) const
{
    double distance;
    distance =  sqrt(pow(cor[0]-otherPoint.cor[0],2) + 
                      pow(cor[1]-otherPoint.cor[1],2) +
                      pow(cor[2]-otherPoint.cor[2],2));
    return distance;
}
//******************************************************
