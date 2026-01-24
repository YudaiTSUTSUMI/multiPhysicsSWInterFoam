scalar maxAlphaCo
(
    runTime.controlDict().get<scalar>("maxAlphaCo")
);

scalar alphaCoNum = 0.0;
scalar meanAlphaCoNum = 0.0;

scalarList alphaCoNumOIDF(OIDFRegions.size(), 0.0);
scalarList meanAlphaCoNumOIDF(OIDFRegions.size(), 0.0);

forAll(OIDFNames, i)
{
    dynamicFvMesh& mesh = OIDFRegions[i];
    surfaceScalarField& phi = phiOIDF[i];
    volScalarField& cellMask = cellMaskOIDF[i];
    immiscibleIncompressibleTwoPhaseMixture& mixture = mixtureOIDF[i];
    
    
    if (mesh.nInternalFaces())
    {
        surfaceScalarField phiMask(localMin<scalar>(mesh).interpolate(cellMask));

        scalarField sumPhi
        (
            mixture.nearInterface()().internalField()
           *fvc::surfaceSum(mag(phiMask*phi))().internalField()
        );

        alphaCoNum = 0.5*gMax(sumPhi/mesh.V().field())*runTime.deltaTValue();

        meanAlphaCoNum =
            0.5*(gSum(sumPhi)/gSum(mesh.V().field()))*runTime.deltaTValue();
    }
    
}

alphaCoNum = max(alphaCoNumOIDF);
meanAlphaCoNum = sum(meanAlphaCoNumOIDF)/OIDFRegions.size();

Info<< "Interface Courant Number mean: " << meanAlphaCoNum
    << " max: " << alphaCoNum << endl;

// ************************************************************************* //
