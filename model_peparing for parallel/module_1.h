#ifndef _MODULE_H_1_
#define _MODULE_H_1_
#include<iostream>
#include "parameters.h"
#include <memory>
#include "particles.h"
#include <vector>
#include "fieldsgrids.h"
#include <cmath>
#include <limits>
#include <bitset>

//************************************************************************
//************************************************************************
// FUNCTION
// Initialization the particles for velocity and position

vector<Particles>* ParticlesLists( GridsPoints***** ptrArray_in, double*** ptrVolumeCellArray_in);

//************************************************************************
//************************************************************************
// FUNCTION
// Initialization the particles for velocity and position 
// Generate lists of particles for bot and top region temp

vector<Particles>* ParticlesListsTemp( GridsPoints***** ptrArray_in, double*** ptrVolumeCellArray_in, double mi0, int ionType_in);



//************************************************************************
//************************************************************************
// FUNCTION
// Transform the ubit64 to vector of position

Vector3 Uint64ToVector3( uint_64 intPos_in);