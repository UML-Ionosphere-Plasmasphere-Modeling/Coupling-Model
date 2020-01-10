#include <iostream>
#include <vector>
using std::vector;

int main()
{
    vector<int> ptrVector = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13,14,15,16,17,18,19,20};
    vector<int> ptrCheck;

    for( auto iterator = ptrVector.begin(); iterator != ptrVector.end(); iterator++)
    {
        if( *iterator % 4 == 0)
        {
            int temp = iterator - ptrVector.begin();
            ptrCheck.push_back( temp);  
        }
    } 

    for( auto iterator = ptrVector.begin(); iterator != ptrVector.end(); iterator++)
    std::cout << *iterator << " ";

    std::cout << std::endl ;
    for( auto iterator = ptrCheck.begin(); iterator != ptrCheck.end(); iterator++)
    std::cout << *iterator << " ";

    std::cout << std::endl ;
    
 
    for( auto iterator = ptrCheck.end()-1; iterator >= ptrCheck.begin(); iterator--)
    {
        (ptrVector)[ *iterator] = 0;
        std::cout << " test " << *iterator ; 

        iterator = ptrCheck.erase(iterator); // if erase the last element, I think the iterator won't go to NULL as end().
        std::cout << " test " << *iterator << std::endl; 
    }

    std::cout << " end test " << std::endl;

    for( auto iterator = ptrVector.begin(); iterator != ptrVector.end(); iterator++)
    std::cout << *iterator << " ";
    
    std::cout << std::endl ;
    for( auto iterator = ptrCheck.begin(); iterator != ptrCheck.end(); iterator++)
    std::cout << *iterator << " ";

    std::cout << std::endl ;
    return 0;
}