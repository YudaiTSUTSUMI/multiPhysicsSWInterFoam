PtrList<IOobject> alphaPhi10HeaderOla(olaRegions.size());
PtrList<surfaceScalarField> alphaPhi10Ola(olaRegions.size());
PtrList<surfaceScalarField> talphaPhi1Corr0Ola(olaRegions.size());
PtrList<surfaceScalarField> alphaPhiUnOla(olaRegions.size());

PtrList<bool> alphaRestartOla(olaRegions.size());

forAll(olaRegions, i)
{   
    dynamicFvMesh& mesh = olaRegions[i];
    immiscibleIncompressibleTwoPhaseMixture& mixture = mixtureOla[i];
    volScalarField& alpha1(mixture.alpha1());
    surfaceScalarField& phi = phiOla[i];
    
    alphaPhi10HeaderOla.set
    (
        i,
        new IOobject
        (
            IOobject::groupName("alphaPhi0", alpha1.group()),
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    );
    IOobject& alphaPhi10Header = alphaPhi10HeaderOla[i];
    
    alphaRestartOla.set
    (
        i,
        new bool
        (
            alphaPhi10Header.typeHeaderOk<surfaceScalarField>(true)
        )
    );
    const bool& alphaRestart = alphaRestartOla[i];
     
    if (alphaRestart)
    {
        Info << "Restarting alpha" << endl;
    }
     
    // MULES flux from previous time-step
    alphaPhi10Ola.set
    (
        i,
        new surfaceScalarField
        (
            alphaPhi10Header,
            phi*fvc::interpolate(alpha1)
        )
    );
     
    // MULES Correction
    //tmp<surfaceScalarField> talphaPhi1Corr0; //how to deal??
    talphaPhi1Corr0Ola.set(i, new surfaceScalarField(0*phi*fvc::interpolate(alpha1)));
     
    // MULES compressed flux is registered in case scalarTransport FO needs it.
    alphaPhiUnOla.set
    (
        i,
        new surfaceScalarField
        (
            IOobject
            (
                "alphaPhiUn",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar(phi.dimensions(), Zero)
        )
    );
}
