#ifndef arrayToField_H
#define arrayToField_H

#include "fvCFD.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class HLLCS Declaration
\*---------------------------------------------------------------------------*/

class arrayToField
{
    // Private data
    volScalarField& field_;
    const List<scalarList>& array_;
    const scalar& Xmin_;
    const scalar& Xmax_;
    const scalar& Ymin_;
    const scalar& Ymax_;
    const scalar& delta_;  
            
public:

    // Constructors
    arrayToField
    (
        volScalarField& field,
        const List<scalarList>& array,
        const scalar& Xmin,
        const scalar& Xmax,
        const scalar& Ymin,
        const scalar& Ymax,
        const scalar& delta
    );


    //- Destructor
    ~arrayToField()
    {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

