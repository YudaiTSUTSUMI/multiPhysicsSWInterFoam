forAll(laplacianRegions, i)
{
	volScalarField& p = pLaplacian[i];
	
	simpleControl& simple = simpleLaplacian[i];
	
	while (simple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(p)
        );

        pEqn.solve();
    }
}

