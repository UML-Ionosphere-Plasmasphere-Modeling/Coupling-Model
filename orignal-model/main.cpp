#include<iostream>
#include<cmath>
#include "parameters.h"
#include "fieldsgrids.h"
#include "particles.h"
#include "vector3.h"
#include "module.h"
 using std::cout;
 using std::endl;
#include "vector"
#include <omp.h>
    // levels of grids for fields and particles


int main()
{
    ProcessFunc();
    
    #pragma omp parallel num_threads(4)
    	{
    		cout << "Test" << endl;
    	}
    
    std::cout << "OKOKOK" << std::endl<< std::endl;
    // t= 0. put in initial condition and create particles initial condition

    // t = n. 
    // Go through the particles list, 
    // 1. calculate the local E and B
    // 2. calculate the new velocity using Boris' Method
    // 3. update the location of the particles
    // 4. update the weighting of N and V on grids points
//    ProcessFunc();
    // t = n+1 , repeat t =n

    return 0;

}