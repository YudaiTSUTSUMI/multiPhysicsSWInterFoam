//#include "postProcess.H"

#include "setRootCaseLists.H"
#include "createTime.H"

regionProperties rp(runTime);
const wordList OIDFNames(rp["OIDF"]);

Info << "Registered field names: " << OIDFNames << endl;

#include "createOIDFMeshes.H"

PtrList<pimpleControl> pimpleOIDF(OIDFRegions.size()); 
PtrList<scalar> cumulativeContErrOIDF(OIDFRegions.size());
forAll(OIDFNames, i)
{
    dynamicFvMesh& mesh = OIDFRegions[i];
    
    cumulativeContErrOIDF.set
    (
        i, 
        new scalar
        (
            uniformDimensionedScalarField
            (
                IOobject
                (
                    "cumulativeContErr",
                    runTime.timeName(),
                    "uniform",
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                dimensionedScalar(dimless, Zero)
            ).value()
        )
    );
    pimpleOIDF.set(i, new pimpleControl(mesh));
}

#include "createTimeControls.H"
#include "createOIDFControls.H"

#include "createOIDFFields.H"
#include "createOIDFAlphaFluxes.H"

PtrList<fv::options> fvOptionsOIDF(OIDFRegions.size());
forAll(OIDFNames, i)
{
    dynamicFvMesh& mesh = OIDFRegions[i];
    fvOptionsOIDF.set(i, new fv::options(mesh));
}

//#include "createUf.H"
PtrList<volScalarField> rAUOIDF(OIDFRegions.size());
PtrList<surfaceVectorField> UfOIDF(OIDFRegions.size());
forAll(OIDFNames, i)
{
    dynamicFvMesh& mesh = OIDFRegions[i];
    volScalarField& rho = rhoOIDF[i];
    volVectorField& U = UOIDF[i];
    
    rAUOIDF.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "rAU",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("rAUf", dimTime/rho.dimensions(), 1.0)
        )
    );
    
    UfOIDF.set
    (
        i,
        new surfaceVectorField
        (
            IOobject
            (
                "Uf",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            fvc::interpolate(U)
        )
    );
}


forAll(OIDFNames, i)
{
    dynamicFvMesh& mesh = OIDFRegions[i];
    volScalarField& cellMask = cellMaskOIDF[i];
    surfaceScalarField& faceMask = faceMaskOIDF[i];
    volScalarField& interpolatedCells = interpolatedCellsOIDF[i];
    
    #include "setCellMask.H"
    #include "setInterpolatedCells.H"
}

#include "multiregionCourantNo.H"
#include "setInitialDeltaT.H"
