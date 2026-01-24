PtrList<IOobject> alphaPhi10HeaderOIDF(OIDFRegions.size());
PtrList<surfaceScalarField> alphaPhi10OIDF(OIDFRegions.size());
PtrList<surfaceScalarField> talphaPhi1Corr0OIDF(OIDFRegions.size());
PtrList<surfaceScalarField> alphaPhiUnOIDF(OIDFRegions.size());

PtrList<bool> alphaRestartOIDF(OIDFRegions.size());

forAll(OIDFRegions, i)
{   
    dynamicFvMesh& mesh = OIDFRegions[i];
    immiscibleIncompressibleTwoPhaseMixture& mixture = mixtureOIDF[i];
    volScalarField& alpha1(mixture.alpha1());
    surfaceScalarField& phi = phiOIDF[i];
    
    alphaPhi10HeaderOIDF.set
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
    IOobject& alphaPhi10Header = alphaPhi10HeaderOIDF[i];
    
    alphaRestartOIDF.set
    (
        i,
        new bool
        (
            alphaPhi10Header.typeHeaderOk<surfaceScalarField>(true)
        )
    );
    const bool& alphaRestart = alphaRestartOIDF[i];
     
    if (alphaRestart)
    {
        Info << "Restarting alpha" << endl;
    }
     
    // MULES flux from previous time-step
    alphaPhi10OIDF.set
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
    talphaPhi1Corr0OIDF.set(i, new surfaceScalarField(0*phi*fvc::interpolate(alpha1)));
     
    // MULES compressed flux is registered in case scalarTransport FO needs it.
    alphaPhiUnOIDF.set
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
