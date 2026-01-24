#ifndef csvArray_H
#define csvArray_H

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class HLLCS Declaration
\*---------------------------------------------------------------------------*/

class csvArray
{
    // Private data
    const word& csvName_;
    const bool hasHeader_;
    const scalar coef_;
    const List<scalarList> array_;
        
            
public:

    // Constructors
    csvArray
    (
        const word& csvName,
        const bool hasHeader = false,
        const scalar coef = 1
    );


    //- Destructor
    ~csvArray()
    {}


    // Member Functions          
    List<scalarList> readCSV();
        
    // Access
    const List<scalarList>& array() const
    {
        return array_;
    }
    
    const scalar& element
    (
        const label x,
        const label y
    ) const
    {
        return array_[y][x];
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

