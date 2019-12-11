#include <iostream>
#include <vector>
#include "spheregrids.h"

using namespace std;

void gridcreate(int );
#define MINS 0.0000001

int main ()
{
    int initial=0;

    cout << "Enter Int for grid resolution" << endl;
 //   cint >> initial >> endl; 
    
    gridcreate(3);
    
    return 0;
}

void gridcreate(int resolution)
{
    int n=1;
    double dist=1.1, tempDist=1.0;
    double a[12][3]=
    {
    {  0.00000000000000000000000000000000,  0.00000000000000000000000000000000,  1.00000000000000000000000000000000},
    {  0.00000000000000000000000000000000,  0.00000000000000000000000000000000, -1.00000000000000000000000000000000},
    {  0.72360679774997893609622678923188, -0.52573111211913359230862852200516,  0.44721359549995792770360480972158},
    {  0.89442719099991585540720961944317,  0.00000000000000000000000000000001, -0.44721359549995792770360480972158},
    {  0.72360679774997893609622678923188,  0.52573111211913359230862852200516,  0.44721359549995792770360480972158},
    {  0.27639320225002100839262197951030,  0.85065080835203987774661982257385, -0.44721359549995792770360480972158},
    { -0.27639320225002100839262197951030,  0.85065080835203987774661982257385,  0.44721359549995792770360480972158},
    { -0.72360679774997893609622678923188,  0.52573111211913359230862852200516, -0.44721359549995792770360480972158},
    { -0.89442719099991585540720961944317, -0.00000000000000000000000000000001,  0.44721359549995792770360480972158},
    { -0.72360679774997893609622678923188, -0.52573111211913359230862852200516, -0.44721359549995792770360480972158},
    { -0.27639320225002100839262197951030, -0.85065080835203987774661982257385,  0.44721359549995792770360480972158},   
    {  0.27639320225002100839262197951030, -0.85065080835203987774661982257385, -0.44721359549995792770360480972158}
    };
    vector<SphereGrids> gridList(12);
    gridList.reserve(50000);
    double (*p)[3]=&a[0];
    int j=1;

    std::vector<SphereGrids>::iterator iter;
    // Initial first 12 points
    for( iter = gridList.begin(); iter != gridList.end(); iter++)
    {
        iter->SetCorrdinate(*p);
        p++;
        cout << j++ << " " << iter->cx() << " " << iter->cy() << " " << iter->cz() <<endl;
    }
    iter=iter-1;

    while (n<=resolution)
    {
    // Create next level resolution points
        for(auto iter1 = gridList.begin(); iter1 <= iter; iter1++)
        {
            for(auto iter2 = iter1+1; iter2 <= iter; iter2++)
            {
      //      cout << n << endl;
                double distBtwPoints = iter1->DistanceBetweenPoints(*iter2);
                if(distBtwPoints <= dist + MINS)
                {
                    double centerPoint[3]={(iter1->cx()+iter2->cx())*0.50,
                                           (iter1->cy()+iter2->cy())*0.50,
                                           (iter1->cz()+iter2->cz())*0.50};
                    gridList.emplace_back(centerPoint);
                    (gridList.end()-1)->CenterProjectionOnSphere();

                    if(tempDist > iter1->DistanceBetweenPoints(*(gridList.end()-1)) - MINS)
                    {
                        tempDist = iter1->DistanceBetweenPoints(*(gridList.end()-1));
                    }

                    cout << j++ << " "<< (gridList.end()-1)->cx()  << " " 
                    << (gridList.end()-1)->cy() << " " << (gridList.end()-1)->cz() << " ";
                    cout << iter1-gridList.begin() << "  ";
                    cout << dist << " " << tempDist << endl;
                    
                }   
            }
        }
    dist=tempDist*1.2;
    iter=gridList.end()-1;
    n++;
    }

}

