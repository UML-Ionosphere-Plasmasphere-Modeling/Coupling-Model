#include <iostream>
#include <vector>
using std::vector;

int main()
{
    vector<int>* ptrVector;
    ptrVector = new vector<int>{ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13,14,15,16,17,18,19,20};
    vector<int>* ptrCheck;
    ptrCheck = new vector<int>;

    for( auto iterator = ptrVector->begin(); iterator != ptrVector->end(); iterator++)
    {
        std::cout << *iterator << std::endl; 
        if( *iterator % 3 == 0)
        {
            int temp = iterator - ptrVector->begin();
            ptrCheck->push_back( temp);  
        }
    } 

    for( auto iterator = ptrVector->begin(); iterator != ptrVector->end(); iterator++)
    std::cout << *iterator << " ";

    std::cout << std::endl ;
    for( auto iterator = ptrCheck->begin(); iterator != ptrCheck->end(); iterator++)
    std::cout << *iterator << " ";

    std::cout << std::endl ;
    

    for( auto iterator = ptrCheck->rbegin(); iterator != ptrCheck->rend(); iterator++)
    {
        (*ptrVector)[ *iterator] = 0;
        ptrCheck->erase((++iterator).base());
    }

    
    for( auto iterator = ptrVector->begin(); iterator != ptrVector->end(); iterator++)
    std::cout << *iterator << " ";
    
    std::cout << std::endl ;
    for( auto iterator = ptrCheck->begin(); iterator != ptrCheck->end(); iterator++)
    std::cout << *iterator << " ";

    std::cout << std::endl ;
    return 0;
}