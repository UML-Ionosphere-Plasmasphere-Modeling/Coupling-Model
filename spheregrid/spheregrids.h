//*****************************************
//SPECIFICATION FILE (spheregrids.h)
//This file gives the specification of grid-points
//
//*****************************************
class SphereGrids
{
    public:

    SphereGrids( double cor[3]);    //Constructor
    SphereGrids();                  //Constructor

    void SetCorrdinate( /* in */ double cor[3]);
        // Precondition:
        // Corrdinate 
        // Postcondition:
        // Set/Reset corrdinate
    void Write() const;
        // Write function
    double cx() const;
        // return cx
    double cy() const;
        // return cy
    double cz() const;
        // return cz        
    double DistanceBetweenPoints(/* in */ SphereGrids otherPoint ) const;
        // Precondition;
        // One another grid point on sphere which should be checked 
        // Postcondition:
        // Return the distance square between two points

    void CenterProjectionOnSphere() ;
        // Precondition:
        // One another grid point on sphere which should be checked
        // Postcondition:
        // Calculate the centerpoint and return the on-sphere cor of it

//    void SetNeighbor();
        // Precondition:
        // 
        // Postcondition:
        // Set/Reset neighbor[6], wrap k++

    private:
    double cor[3];
//    SphereGrids* neighbor[6];
//    int k;  // for the assignment of neighbor
};