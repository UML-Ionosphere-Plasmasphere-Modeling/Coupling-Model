#include <iostream>
#include <string>
#include "H5Cpp.h"
#include "vector3.h"
#include "parameters.h"
#include 
using namespace H5;

const H5std_string FLIE_NAME( "GridsData.h5");
const H5std_string DATASET_NAME( "ArrayOfStructures");
const H5std_string MEMBER1( "pos3_name");
const H5std_string MEMBER2( "e3_name");
const H5std_string MEMBER3( "b3_name");
const H5std_string MEMBER4( "v3_name");
const H5std_string MEMBER5( "density_name");

const int RANK = 4;

