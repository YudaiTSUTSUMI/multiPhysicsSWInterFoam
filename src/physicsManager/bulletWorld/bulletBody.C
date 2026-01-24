#include "bulletBody.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::bulletBody::bulletBody
(
    const dictionary& bodyDict
)
:
    bodyDict_(bodyDict),
    dynamic_(false),
    P_(Zero),
    P0_(Zero),
    Q_(I),
    Q0_(I),
    initialP_(Zero),
    localCoM_(Zero),
    localCoM0_(Zero),
    v_(Zero),
    v0_(Zero),
    ev_(Zero),
    ev0_(Zero),
    a_(Zero),
    pi_(Zero),
    pi0_(Zero),
    tau_(Zero),
    rhoSolid_(0.0),
    massSolid_(0.0),
    inertiaSolid_(Zero),
    CoMSolid_(Zero),
    massFluid_(0.0),
    inertiaFluid_(Zero),
    inertiaTotal_(Zero),
    inertiaTotal0_(Zero),
    CoMFluid_(Zero)
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //




// ************************************************************************* //
