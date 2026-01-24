scalar maxAlphaCo
(
    runTime.controlDict().get<scalar>("maxAlphaCo")
);
scalar alphaCoNum = 0.0;
scalar meanAlphaCoNum = 0.0;

scalarList alphaCoNumAllVOF(allVOFRegions.size(), 0.0);
scalarList meanAlphaCoNumAllVOF(allVOFRegions.size(), 0.0);

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

        alphaCoNumAllVOF[i] = 0.5*gMax(sumPhi/mesh.V().field())*runTime.deltaTValue();

        meanAlphaCoNumAllVOF[i] =
            0.5*(gSum(sumPhi)/gSum(mesh.V().field()))*runTime.deltaTValue();
    }
    
}

forAll(olaNames, i)
{
    dynamicFvMesh& mesh = olaRegions[i];
    surfaceScalarField& phi = phiOla[i];
    immiscibleIncompressibleTwoPhaseMixture& mixture = mixtureOla[i];
    
    
    if (mesh.nInternalFaces())
    {
        scalarField sumPhi
        (
            mixture.nearInterface()().internalField()
           *fvc::surfaceSum(mag(phi))().internalField()
        );

        alphaCoNumAllVOF[OIDFNames.size() + i] = 0.5*gMax(sumPhi/mesh.V().field())*runTime.deltaTValue();

        meanAlphaCoNumAllVOF[OIDFNames.size() + i]  =
            0.5*(gSum(sumPhi)/gSum(mesh.V().field()))*runTime.deltaTValue();
    }
    
}


if(allVOFRegions.size())
{
    alphaCoNum = max(alphaCoNumAllVOF);
    meanAlphaCoNum = sum(meanAlphaCoNumAllVOF)/allVOFRegions.size();
}


Info<< "Interface Courant Number mean: " << meanAlphaCoNum
    << " max: " << alphaCoNum << endl;

// ************************************************************************* //
