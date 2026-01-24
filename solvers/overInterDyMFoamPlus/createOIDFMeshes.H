PtrList<bulletOversetFvMesh> OIDFRegions(OIDFNames.size());

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
