#ifndef _STRUCTDEF_H_
#define _STRUCTDEF_H_
#include <iostream>
#include "parameters.h"


struct structg 
    {
        int face;
        int ig;  int jg;  int kg;
        int iw;  int jw;  int kw;
        double vx; double vy; double vz;
//        Vector3 vp;
//        double mass;
    };



// Getface     
inline uint_64 Getface( double px_in, double py_in, double pz_in)
{
    return abs(px_in) > abs(py_in) ?
                abs(px_in) > abs(pz_in) ? 
                    px_in > 0 ? 0 : 3 :
                    pz_in > 0 ? 2 : 5 :
                abs(py_in) > abs(pz_in) ?
                    py_in > 0 ? 1 : 4 :
                    pz_in > 0 ? 2 : 5 ;
}

#endif