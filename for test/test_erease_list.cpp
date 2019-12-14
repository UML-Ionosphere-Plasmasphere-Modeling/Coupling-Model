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

    // int a[4] = {4,5,6,7};
    for(auto x = llist.begin(); x!= llist.end(); x++)
    {
    std::cout << " test " << std::endl;
    if( *x == 4 || *x ==5 || *x ==6)    
    {
    x= llist.erase(x);
    
    }

    }
    for(auto x = llist.begin(); x!= llist.end(); x++)
    {
        std::cout << *x << " " ;
    }
    std::cout << std::endl;

    return 0;
}