#include<iostream>
#include "mathutil.h"
#include "parameters.h"
#include "particles.h"
#include "vector3.h"
#include "structdef.h"
#include <bitset>

// private:
// uint_64 posUint;
// Vector3 vp;

// FUNCTION //Constructor
Particles::Particles(   uint_64 posUint_in,
                        Vector3 posP_in, 
                        Vector3 vp_in,
                        double weightNi_in,
                        double mu_in)
{
    posUint = posUint_in;
    posP = posP_in;
    vp = vp_in; 
    weightNi = weightNi_in;
    mu = mu_in;
}

// FUNCTION //Default Constructor
Particles::Particles()
{
    posUint = 0;
    posP = Vector3( 0.0, 0.0, 0.0);
    vp= Vector3(0.0, 0.0, 0.0);
    weightNi = 0.0;
    mu = 0.0;
}

//************************************************************************
//************************************************************************
// FUNCTION // Update uint_64 IN 
// And return a int "0" means in the main domain
// "1" means out of the main domain
//************************************************************************
//************************************************************************
int Particles::UpdateUint_64()
{
    uint_64 face = 0, ip = 0, jp = 0, kp = 0;
    double px, py, pz;
    double temp[2];
    int check=0;

    px = posP.x(); py = posP.y(); pz = posP.z();
    // 3. update double x y z
    px += vp.x() * tstep;
    py += vp.y() * tstep;
    pz += vp.z() * tstep;

    // 4. transfor to face ip kp jp
    // 4.1 radial kp
    double L = sqrt( px*px + py*py + pz*pz )/radius;

    // check if in the main domain
    if( L > LMax_maindomain || L < LMin_maindomain) 
    {
        check = 1;
        return check;
    }    
    else
    {
        kp = static_cast<uint_64>( floor( log10(L / LMin)/logRatio *cellSize1)); 
        // 4.2 XYZtoUV, note that -1<UV<1 
        face = Getface(px, py, pz);
        switch (face)
        {
        case 0: temp[0] = py/px; temp[1] = pz/px; break;
        case 1: temp[0] =-px/py; temp[1] = pz/py; break;
        case 2: temp[0] = py/pz; temp[1] =-px/pz; break;
        case 3: temp[0] = py/px; temp[1] =-pz/px; break;
        case 4: temp[0] =-px/py; temp[1] =-pz/py; break;
        default:temp[0] =-py/pz; temp[1] =-px/pz; break;
        }
        // 4.3 UVtoST, note that 0<ST<1
        for (int i=0; i<=1; i++)
        {
        if (temp[i] >= 0) temp[i] = 0.5 * std::sqrt(1 + 3*temp[i]);
        else            temp[i] = 1 - 0.5 * std::sqrt(1 - 3*temp[i]);
        }
        // 4.4 STtoIpJp
        ip= static_cast<unsigned int>(floor(temp[0] * particlesGridsSize ));
        jp= static_cast<unsigned int>(floor(temp[1] * particlesGridsSize ));

        // 5. F ip jp kp to Uint_64
        posUint = face << 61 ;
        for( int i = 0; i < particlesGridsLevel; i++)
        {
        posUint += (((ip >> particlesGridsLevel-1-i) & 1 )<< 60 - i *3) 
                    + (((jp >> particlesGridsLevel-1-i) & 1 )<< 60-1 - i *3)
                    + (((kp >> particlesGridsLevel-1-i) & 1 )<< 60-2 - i *3) ;
        }
/*        std::cout << " next " << std::endl;
    std::cout << std::bitset<64>(face) << " " << face << std::endl;
    std::cout << std::bitset<64>(ip) << " " << ip << std::endl;
    std::cout << std::bitset<64>(jp) << " " << jp << std::endl;
    std::cout << std::bitset<64>(kp) << " " << kp << std::endl << std::endl;
    int pause;
    std::cin >>pause;
*/
    return check;        
    }
}

//************************************************************************
//************************************************************************
// FUNCTION // BorisMethod update velocity
// And return a int "0" means in the main domain
// "1" means out of the main domain
//
//************************************************************************
//************************************************************************
int Particles::BorisMethod( struct structg *strg_in, GridsPoints***** ptrArray_in, double mi0_in)
{
    // 1. pre-set some variables
    // 1.1 qtm = q * dt / m /2 
    double qtm = qi0 * tstep / mi0_in / 2;
    // 1.2 local B and E
//    std::cout << std::endl;
//    std::cout << strg_in->face << strg_in->ig << strg_in->jg << strg_in->kg << std::endl;

    Vector3 tempb1=ptrArray_in[strg_in->face][strg_in->ig+1][strg_in->jg+1][strg_in->kg]->B3();
    Vector3 tempb2=ptrArray_in[strg_in->face][strg_in->ig+2][strg_in->jg+1][strg_in->kg]->B3();
    Vector3 tempb3=ptrArray_in[strg_in->face][strg_in->ig+2][strg_in->jg+2][strg_in->kg]->B3();
    Vector3 tempb4=ptrArray_in[strg_in->face][strg_in->ig+1][strg_in->jg+2][strg_in->kg]->B3();
    Vector3 tempb5=ptrArray_in[strg_in->face][strg_in->ig+1][strg_in->jg+1][strg_in->kg+1]->B3();
    Vector3 tempb6=ptrArray_in[strg_in->face][strg_in->ig+2][strg_in->jg+1][strg_in->kg+1]->B3();
    Vector3 tempb7=ptrArray_in[strg_in->face][strg_in->ig+2][strg_in->jg+2][strg_in->kg+1]->B3();
    Vector3 tempb8=ptrArray_in[strg_in->face][strg_in->ig+1][strg_in->jg+2][strg_in->kg+1]->B3();
    Vector3 tempe1=ptrArray_in[strg_in->face][strg_in->ig+1][strg_in->jg+1][strg_in->kg]->E3();
    Vector3 tempe2=ptrArray_in[strg_in->face][strg_in->ig+2][strg_in->jg+1][strg_in->kg]->E3();
    Vector3 tempe3=ptrArray_in[strg_in->face][strg_in->ig+2][strg_in->jg+2][strg_in->kg]->E3();
    Vector3 tempe4=ptrArray_in[strg_in->face][strg_in->ig+1][strg_in->jg+2][strg_in->kg]->E3();
    Vector3 tempe5=ptrArray_in[strg_in->face][strg_in->ig+1][strg_in->jg+1][strg_in->kg+1]->E3();
    Vector3 tempe6=ptrArray_in[strg_in->face][strg_in->ig+2][strg_in->jg+1][strg_in->kg+1]->E3();
    Vector3 tempe7=ptrArray_in[strg_in->face][strg_in->ig+2][strg_in->jg+2][strg_in->kg+1]->E3();
    Vector3 tempe8=ptrArray_in[strg_in->face][strg_in->ig+1][strg_in->jg+2][strg_in->kg+1]->E3();
    Vector3 gradB1=ptrArray_in[strg_in->face][strg_in->ig+1][strg_in->jg+1][strg_in->kg]->GradB3();
    Vector3 gradB2=ptrArray_in[strg_in->face][strg_in->ig+2][strg_in->jg+1][strg_in->kg]->GradB3();
    Vector3 gradB3=ptrArray_in[strg_in->face][strg_in->ig+2][strg_in->jg+2][strg_in->kg]->GradB3();
    Vector3 gradB4=ptrArray_in[strg_in->face][strg_in->ig+1][strg_in->jg+2][strg_in->kg]->GradB3();
    Vector3 gradB5=ptrArray_in[strg_in->face][strg_in->ig+1][strg_in->jg+1][strg_in->kg+1]->GradB3();
    Vector3 gradB6=ptrArray_in[strg_in->face][strg_in->ig+2][strg_in->jg+1][strg_in->kg+1]->GradB3();
    Vector3 gradB7=ptrArray_in[strg_in->face][strg_in->ig+2][strg_in->jg+2][strg_in->kg+1]->GradB3();
    Vector3 gradB8=ptrArray_in[strg_in->face][strg_in->ig+1][strg_in->jg+2][strg_in->kg+1]->GradB3();
    Vector3 tempPos1=ptrArray_in[strg_in->face][strg_in->ig+1][strg_in->jg+1][strg_in->kg]->Pos3();
    Vector3 tempPos2=ptrArray_in[strg_in->face][strg_in->ig+2][strg_in->jg+1][strg_in->kg]->Pos3();
    Vector3 tempPos3=ptrArray_in[strg_in->face][strg_in->ig+2][strg_in->jg+2][strg_in->kg]->Pos3();
    Vector3 tempPos4=ptrArray_in[strg_in->face][strg_in->ig+1][strg_in->jg+2][strg_in->kg]->Pos3();
    Vector3 tempPos5=ptrArray_in[strg_in->face][strg_in->ig+1][strg_in->jg+1][strg_in->kg+1]->Pos3();
    Vector3 tempPos6=ptrArray_in[strg_in->face][strg_in->ig+2][strg_in->jg+1][strg_in->kg+1]->Pos3();
    Vector3 tempPos7=ptrArray_in[strg_in->face][strg_in->ig+2][strg_in->jg+2][strg_in->kg+1]->Pos3();
    Vector3 tempPos8=ptrArray_in[strg_in->face][strg_in->ig+1][strg_in->jg+2][strg_in->kg+1]->Pos3();
    // iw jw kw are 0, 1, 2, ... , cellSize-1
    double w1 = 1- (strg_in->iw +1) * (strg_in->jw +1) * (strg_in->kw +1) / cellSize3;
    double w2 = 1- (cellSize1- strg_in->iw)* (strg_in->jw +1) * (strg_in->kw +1) / cellSize3;
    double w3 = 1- (cellSize1- strg_in->iw) * (cellSize1- strg_in->jw) * (strg_in->kw +1) / cellSize3;
    double w4 = 1- (strg_in->iw +1) * (cellSize1- strg_in->jw )* (strg_in->kw +1)/ cellSize3;
    double w5 = 1- (strg_in->iw +1) * (strg_in->jw +1) * (cellSize1- strg_in->kw) / cellSize3;
    double w6 = 1- (cellSize1- strg_in->iw)* (strg_in->jw +1) * (cellSize1- strg_in->kw) / cellSize3;
    double w7 = 1- (cellSize1- strg_in->iw) * (cellSize1- strg_in->jw) * (cellSize1- strg_in->kw) / cellSize3;
    double w8 = 1- (strg_in->iw +1) * (cellSize1- strg_in->jw )* (cellSize1- strg_in->kw) / cellSize3;
    Vector3 tempb;
    tempb.Setx(tempb1.x()*w1 + tempb2.x()*w2 + tempb3.x()*w3 + tempb4.x()*w4 
                + tempb5.x()*w5 + tempb6.x()*w6 + tempb7.x()*w7 + tempb8.x()*w8);
                
    tempb.Sety(tempb1.y()*w1 + tempb2.y()*w2 + tempb3.y()*w3 + tempb4.y()*w4 
                + tempb5.y()*w5 + tempb6.y()*w6 + tempb7.y()*w7 + tempb8.y()*w8);
                
    tempb.Setz(tempb1.z()*w1 + tempb2.z()*w2 + tempb3.z()*w3 + tempb4.z()*w4 
                + tempb5.z()*w5 + tempb6.z()*w6 + tempb7.z()*w7 + tempb8.z()*w8);
//    std::cout << " tempb " << tempb.x() << " " << tempb.y() << " " << tempb.z() << " " << std::endl;
    
    Vector3 tempe;
    tempe.Setx(tempe1.x()*w1 + tempe2.x()*w2 + tempe3.x()*w3 + tempe4.x()*w4 
                + tempe5.x()*w5 + tempe6.x()*w6 + tempe7.x()*w7 + tempe8.x()*w8);
                
    tempe.Sety(tempe1.y()*w1 + tempe2.y()*w2 + tempe3.y()*w3 + tempe4.y()*w4 
                + tempe5.y()*w5 + tempe6.y()*w6 + tempe7.y()*w7 + tempe8.y()*w8);
                
    tempe.Setz(tempe1.z()*w1 + tempe2.z()*w2 + tempe3.z()*w3 + tempe4.z()*w4 
                + tempe5.z()*w5 + tempe6.z()*w6 + tempe7.z()*w7 + tempe8.z()*w8);
    
//    std::cout << " tempe1 " << tempe.norm() <<" xyz "<< tempe.x() << " " << tempe.y() << " " << tempe.z() << " " << std::endl;
    
    // 1 Revised E ( including gravity which is not accurate to the loation of unique particles)
    Vector3 tempPos;

    tempPos.Setx(tempPos1.x()*w1 + tempPos2.x()*w2 + tempPos3.x()*w3 + tempPos4.x()*w4 
                + tempPos5.x()*w5 + tempPos6.x()*w6 + tempPos7.x()*w7 + tempPos8.x()*w8);
                
    tempPos.Sety(tempPos1.y()*w1 + tempPos2.y()*w2 + tempPos3.y()*w3 + tempPos4.y()*w4 
                + tempPos5.y()*w5 + tempPos6.y()*w6 + tempPos7.y()*w7 + tempPos8.y()*w8);
                
    tempPos.Setz(tempPos1.z()*w1 + tempPos2.z()*w2 + tempPos3.z()*w3 + tempPos4.z()*w4 
                + tempPos5.z()*w5 + tempPos6.z()*w6 + tempPos7.z()*w7 + tempPos8.z()*w8);

//    std::cout << " tempPos1 " << tempPos.norm() << " " << tempPos.norm2() << " xyz " << tempPos.x() << " " << tempPos.y() << " " << tempPos.z() << " " << std::endl;
//    std::cout << std::endl << " L^2 " << radius*radius/tempPos.norm2() << " qtm " << qtm << std::endl;
    
    tempPos = tempPos.NormalizedVector().ScaleProduct(-1*gravity*radius*radius/tempPos.norm2() *mi0_in/qi0);
    
    // 2 Revised E ( including gradB at the location of unique particles)
    Vector3 tempGradB;
    tempGradB.Setx(gradB1.x()*w1 + gradB2.x()*w2 + gradB3.x()*w3 + gradB4.x()*w4 
                + gradB5.x()*w5 + gradB6.x()*w6 + gradB7.x()*w7 + gradB8.x()*w8);
    tempGradB.Sety(gradB1.x()*w1 + gradB2.x()*w2 + gradB3.x()*w3 + gradB4.x()*w4 
                + gradB5.x()*w5 + gradB6.x()*w6 + gradB7.x()*w7 + gradB8.x()*w8);
    tempGradB.Setz(gradB1.x()*w1 + gradB2.x()*w2 + gradB3.x()*w3 + gradB4.x()*w4 
                + gradB5.x()*w5 + gradB6.x()*w6 + gradB7.x()*w7 + gradB8.x()*w8);

//    std::cout << std::endl << "e " << tempe.x() << " b " << tempb.x() << " pos " << tempPos.x() << " gradB " << tempGradB.x() << std::endl; 

    // mu / e , mu = m * v^2 / e / B
    tempGradB = tempGradB.ScaleProduct( mu / qi0 );

    // revised E including gravity and gradB
    tempe = tempe.PlusProduct(tempPos).PlusProduct(tempGradB.ScaleProduct( -1.0));
    
//    std::cout << " tempPos2 " << tempPos.norm() << " xyz " << tempPos.x() << " " << tempPos.y() << " " << tempPos.z() << " " << std::endl;
//    std::cout << " scale " << -1*gravity*radius*radius/tempPos.norm2() *mi0/qi0 << std::endl;
//    std::cout << " tempe2 " << tempe.norm() <<" xyz " << tempe.x() << " " << tempe.y() << " " << tempe.z() << " " << std::endl;
    // 1.3 Vector3 t ; Vector3 s
    Vector3 t = tempb.ScaleProduct(qtm);
    Vector3 s = t.ScaleProduct(2/(1+t.norm2()));
    // 2. Boris equations
    
    // 2.1 equation I: v1 = v + E * qtm
    Vector3 v1 = vp.PlusProduct(tempe.ScaleProduct(qtm));
    // 2.2 equation II: v2 = v1 + v1 X t
    Vector3 v2 = v1.PlusProduct(v1.CrossProduct(t));
    // 2.3 equation III: v3 = v1 + v2 X s
    Vector3 v3 = v1.PlusProduct(v2.CrossProduct(s));
    // 2.4 equation IV: v = v3 + E * qtm
    vp = v3.PlusProduct(tempe.ScaleProduct(qtm));
    // 3. update the postion 
    return UpdateUint_64();

}


// FUNCTION // output structg
//void Particles::OutputStrg()
//{}