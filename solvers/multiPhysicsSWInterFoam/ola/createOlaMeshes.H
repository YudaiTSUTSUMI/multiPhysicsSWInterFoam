forAll(olaNames, i)
{
    Info<< "Create ola mesh for region " << olaNames[i]
        << " for time = " << runTime.timeName() << nl << endl;

    olaRegions.set
    (
        i,
        dynamicFvMesh::New
        (
            IOobject
            (
                olaNames[i],
                runTime.timeName(),
                runTime,
                IOobject::MUST_READ
            )
        ).ptr()
    );
    
    Info << endl;
}

PtrList<pimpleControl> pimpleOla(olaRegions.size()); 
PtrList<scalar> cumulativeContErrOla(olaRegions.size());
forAll(olaNames, i)
{
    dynamicFvMesh& mesh = olaRegions[i];
    
    //#include "initContinuityErrs.H"
    //#include "createControl.H"
    cumulativeContErrOla.set
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
    pimpleOla.set(i, new pimpleControl(mesh));
}
