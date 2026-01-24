PtrList<volScalarField> pLaplacian(laplacianRegions.size());

forAll(laplacianRegions, i)
{
    fvMesh& mesh = laplacianRegions[i];
    
    pLaplacian.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "p",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )
    );
}

