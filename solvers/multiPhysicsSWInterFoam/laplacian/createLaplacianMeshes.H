forAll(laplacianNames, i)
{
    Info<< "Create laplacian mesh for region " << laplacianNames[i]
        << " for time = " << runTime.timeName() << nl << endl;

    laplacianRegions.set
    (
        i,
        new fvMesh
        (
            IOobject
            (
                laplacianNames[i],
                runTime.timeName(),
                runTime,
                IOobject::MUST_READ
            )
        )
    );
}

PtrList<simpleControl> simpleLaplacian(laplacianRegions.size());

forAll(laplacianNames, i)
{
    fvMesh& mesh = laplacianRegions[i];
    
    simpleLaplacian.set(i, new simpleControl(mesh));
}
