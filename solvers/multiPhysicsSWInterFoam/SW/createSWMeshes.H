forAll(SWNames, i)
{
    Info<< "Create SW mesh for region " << SWNames[i]
        << " for time = " << runTime.timeName() << nl << endl;

    SWRegions.set
    (
        i,
        new fvMesh
        (
            IOobject
            (
                SWNames[i],
                runTime.timeName(),
                runTime,
                IOobject::MUST_READ
            )
        )
    );
}
