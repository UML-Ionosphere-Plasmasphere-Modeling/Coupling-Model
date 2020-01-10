#include <iostream>
#include <omp.h>
#include <vector>
using std::vector;
// go through ptrVector and change ptrCheck
void subdomain( vector<int>& ptrVector, vector<int>& ptrCheck, int istart, int ipoints)
{
    int i;
    vector<int>::iterator iter = ptrVector.begin();
    for (i = 0; i < ipoints; i++)
    {
        if( *(iter+istart+i)%3==0 )
        {
            int item = iter+istart+i - ptrVector.begin();
            #pragma omp critical
            ptrCheck.push_back(item);
        }
    }
}
// go through ptrCheck and change ptrVector
void subdomain2( vector<int>& ptrVector, vector<int>& ptrVector2, vector<int>& ptrCheck, int istart, int ipoints)
{
    int i;
    vector<int>::iterator iter = ptrVector2.begin();
    vector<int>::iterator iterCheck;
    for (i = 0; i < ipoints; i++)
    {
        int temp = *(iter+istart+i);
        
        #pragma omp critical
        if( ptrCheck.size()!=0)
        {
        ptrVector[*(ptrCheck.end()-1)]=temp;
        iterCheck = ptrCheck.erase(ptrCheck.end()-1);
        iterCheck--;
        } else
        {ptrVector.push_back(temp);
        }
        
    }
}
void sub( vector<int>& ptrVector, vector<int>& ptrVector2, vector<int>& ptrCheck)
{
    int npoints = ptrVector.size();
    int cpoints = ptrVector2.size();
    int iam, nt, ipoints, istart;
    #pragma omp parallel default(shared) private(iam,nt,ipoints,istart)
    {
        iam = omp_get_thread_num();
        nt =  omp_get_num_threads();
        ipoints = npoints / nt;    /*size of partition*/
        istart = iam*ipoints;  /*starting array index*/
        if (iam == nt-1)     /*last thread may do more*/
        ipoints = npoints - istart;
        subdomain(ptrVector, ptrCheck, istart, ipoints);
        #pragma omp critical
        std::cout << iam << " " << istart << " " << ipoints << std::endl;

        #pragma omp barrier
        std::cout << " step 2 " << std::endl;

        ipoints = cpoints /nt;
        istart = iam * ipoints;
        if(iam == nt-1)
        ipoints = cpoints - istart;
        subdomain2(ptrVector, ptrVector2, ptrCheck, istart, ipoints);
        #pragma omp critical
        std::cout << iam << " " << istart << " " << ipoints << std::endl;


    }
}

int main()
{
    vector<int> ptrVector={1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20};
    vector<int> ptrVector2(20,0);
    vector<int> ptrCheck;

    sub(ptrVector, ptrVector2, ptrCheck);   
    
    for( auto iterator = ptrVector.begin(); iterator != ptrVector.end(); iterator++)
    std::cout << *iterator << " ";    
    std::cout << std::endl ;
    for( auto iterator = ptrCheck.begin(); iterator != ptrCheck.end(); iterator++)
    std::cout << *iterator << " ";  
    std::cout << std::endl ;
    return 0;
}