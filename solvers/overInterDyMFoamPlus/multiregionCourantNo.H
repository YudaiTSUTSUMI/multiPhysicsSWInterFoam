
scalar CoNum = 0.0;
scalar meanCoNum = 0.0;

scalarList CoNumOIDF(OIDFRegions.size(), 0.0);
scalarList meanCoNumOIDF(OIDFRegions.size(), 0.0);

forAll(OIDFNames, i)
{
    dynamicFvMesh& mesh = OIDFRegions[i];
    surfaceScalarField& phi = phiOIDF[i];
    volScalarField& cellMask = cellMaskOIDF[i];
    volScalarField& interpolatedCells = interpolatedCellsOIDF[i];
    
    
    if (mesh.nInternalFaces())
    {
        surfaceScalarField phiMask
        (
            localMin<scalar>(mesh).interpolate(cellMask + interpolatedCells)
        );
        scalarField sumPhi
        (
            fvc::surfaceSum(mag(phiMask*phi))().internalField()
        );
        CoNumOIDF[i] = 0.5*gMax(sumPhi/mesh.V().field())*runTime.deltaTValue();
        meanCoNumOIDF[i] =
            0.5*(gSum(sumPhi)/gSum(mesh.V().field()))*runTime.deltaTValue();
    }
    
}

CoNum = max(CoNumOIDF);
meanCoNum = sum(meanCoNumOIDF)/OIDFRegions.size();

Info<< "Courant Number mean: " << meanCoNum
    << " max: " << CoNum << endl;

