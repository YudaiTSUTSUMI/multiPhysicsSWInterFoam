forAll(OIDFNames, i)
{
    Info<< "Create OIDF mesh for region " << OIDFNames[i]
        << " for time = " << runTime.timeName() << nl << endl;
        
    OIDFRegions.set
    (
        i,
        new bulletOversetFvMesh
        (
            IOobject
            (
                OIDFNames[i],
                runTime.timeName(),
                runTime,
                IOobject::MUST_READ
            )
        )
    );
    
    Info << endl;
}

PtrList<pimpleControl> pimpleOIDF(OIDFRegions.size()); 
PtrList<scalar> cumulativeContErrOIDF(OIDFRegions.size());
forAll(OIDFNames, i)
{
    dynamicFvMesh& mesh = OIDFRegions[i];
    
    //#include "initContinuityErrs.H"
    //#include "createControl.H"
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