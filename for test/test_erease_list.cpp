#include <iostream>
#include <list>

int main()
{
    std::list<int> llist;
    for (size_t ix = 0; ix != 20; ++ix)
    llist.push_back(ix);

    for(auto x = llist.begin(); x!= llist.end(); x++)
    {
        std::cout << *x << " " ;
    }
    std::cout << std::endl;

    int a[4] = {4,5,6,7};
    for(auto x = llist.begin(); x!= llist.end(); x++)
    {
        
    
    llist.erase(x);
    std::cout << x << " " ;

    }
    for(auto x = llist.begin(); x!= llist.end(); x++)
    {
        std::cout << *x << " " ;
    }
    std::cout << std::endl;

    return 0;
}