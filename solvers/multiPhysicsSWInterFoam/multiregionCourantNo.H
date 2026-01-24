scalar CoNum = 0.0;
scalar meanCoNum = 0.0;

scalarList CoNumAllVOF(allVOFRegions.size(), 0.0);
scalarList meanCoNumAllVOF(allVOFRegions.size(), 0.0);

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
        CoNumAllVOF[i] = 0.5*gMax(sumPhi/mesh.V().field())*runTime.deltaTValue();
        meanCoNumAllVOF[i] =
            0.5*(gSum(sumPhi)/gSum(mesh.V().field()))*runTime.deltaTValue();
    }
}

forAll(olaNames, i)
{
    dynamicFvMesh& mesh = olaRegions[i];
    surfaceScalarField& phi = phiOla[i];
    
    if (mesh.nInternalFaces())
    {
        scalarField sumPhi
        (
            fvc::surfaceSum(mag(phi))().internalField()
        );
        CoNumAllVOF[OIDFNames.size() + i] = 0.5*gMax(sumPhi/mesh.V().field())*runTime.deltaTValue();
        meanCoNumAllVOF[OIDFNames.size() + i] =
            0.5*(gSum(sumPhi)/gSum(mesh.V().field()))*runTime.deltaTValue();
    }
}

if(allVOFRegions.size())
{
    CoNum = max(CoNumAllVOF);
    meanCoNum = sum(meanCoNumAllVOF)/allVOFRegions.size();
}

Info << "Courant Number mean: " << meanCoNum
    << " max: " << CoNum << endl;

