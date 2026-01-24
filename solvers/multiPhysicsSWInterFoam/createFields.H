#include "readGravitationalAcceleration.H"
const vector gValue = g.value();
const dimensionedScalar magg = mag(g);

#include "createOIDFFields.H"
#include "createSWFields.H"
#include "createOlaFields.H"
#include "createLaplacianFields.H"

PtrList<fvMesh> allVOFRegions(OIDFNames.size() + olaNames.size());
wordList allVOFNames(OIDFNames.size() + olaNames.size());

forAll(OIDFRegions, i)
{
    allVOFRegions.set(i, &OIDFRegions[i]);  
    allVOFNames[i] = OIDFNames[i];
}

forAll(olaRegions, i)
{
    allVOFRegions.set(OIDFRegions.size() + i, &olaRegions[i]);  
    allVOFNames[OIDFRegions.size() + i] = olaNames[i];
}

regionConnectors regionConnector 
(
    runTime,
    g.value(),
    SWRegions,
    allVOFRegions,
    laplacianRegions,
    SWNames,
    allVOFNames,
    laplacianNames
);
//FatalErrorInFunction <<  abort(FatalError);
