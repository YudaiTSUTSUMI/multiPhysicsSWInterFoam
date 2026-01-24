#include "createOlaPorosities.H"


PtrList<volScalarField> p_rghOla(olaRegions.size());
PtrList<volVectorField> UOla(olaRegions.size());
PtrList<surfaceScalarField> phiOla(olaRegions.size());
PtrList<volScalarField> rhoOla(olaRegions.size());
PtrList<surfaceScalarField> rhoPhiOla(olaRegions.size());

PtrList<volScalarField> ghOla(olaRegions.size());
PtrList<surfaceScalarField> ghfOla(olaRegions.size());

PtrList<volScalarField> pOla(olaRegions.size());

PtrList<immiscibleIncompressibleTwoPhaseMixture> mixtureOla(olaRegions.size());

PtrList<dimensionedScalar> ghRefOla(olaRegions.size());

PtrList<label> pRefCellOla(olaRegions.size());
PtrList<scalar> pRefValueOla(olaRegions.size());

PtrList<incompressible::turbulenceModel> turbulenceOla(olaRegions.size());

PtrList<IOMRFZoneList> MRFOla(olaRegions.size());


forAll(olaNames, i)
{
    Info<< "*** Region " << olaRegions[i].name() << nl << endl;
    dynamicFvMesh& mesh = olaRegions[i];
    pimpleControl& pimple = pimpleOla[i];
    
    #include "createRDeltaT.H"

    Info<< "Reading field p_rgh\n" << endl;
    p_rghOla.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "p_rgh",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )
    );
    volScalarField& p_rgh = p_rghOla[i];

    Info<< "Reading field U\n" << endl;
    UOla.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "U",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )
    );
    volVectorField& U = UOla[i];

    //#include "createPhi.H"
    Info<< "Reading/calculating face flux field phi\n" << endl;
    phiOla.set
    (
        i,
        new surfaceScalarField
        (
            IOobject
            (
                "phi",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fvc::flux(U)
        )
    );
    surfaceScalarField& phi = phiOla[i];


    Info<< "Reading transportProperties\n" << endl;
    mixtureOla.set
    (
        i,
        new immiscibleIncompressibleTwoPhaseMixture
        (
            U,
            phi
        )
    );
    //immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);
    immiscibleIncompressibleTwoPhaseMixture& mixture = mixtureOla[i];

    volScalarField& alpha1(mixture.alpha1());
    volScalarField& alpha2(mixture.alpha2());

    const dimensionedScalar& rho1 = mixture.rho1();
    const dimensionedScalar& rho2 = mixture.rho2();


    // Need to store rho for ddt(rho, U)
    rhoOla.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "rho",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT
            ),
            alpha1*rho1 + alpha2*rho2
        )
    );
    volScalarField& rho = rhoOla[i];
    rho.oldTime();
   

    // Mass flux
    rhoPhiOla.set
    (
        i,
        new surfaceScalarField
        (
            IOobject
            (
                "rhoPhi",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            fvc::interpolate(rho)*phi
        )
    );


    // Construct incompressible turbulence model
    turbulenceOla.set
    (
        i,
        incompressible::turbulenceModel::New
        (
            U,
            phi,
            mixture
        )
    );
    
    turbulenceOla[i].validate();

    
    //#include "readhRef.H"
    Info<< "\nReading hRef" << endl;
    uniformDimensionedScalarField hRef
    (
        IOobject
        (
            "hRef",
            runTime.constant(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        dimensionedScalar(dimLength, Zero)
    );
    //#include "gh.H"
    Info<< "Calculating field g.h\n" << endl;
    ghRefOla.set
    (
        i,
        new dimensionedScalar
        (
            mag(g.value()) > SMALL
          ? g & (cmptMag(g.value())/mag(g.value()))*hRef
          : dimensionedScalar("ghRef", g.dimensions()*dimLength, 0)
        )
    );
    dimensionedScalar& ghRef = ghRefOla[i];
    
    //volScalarField gh("gh", (g & mesh.C()) - ghRef);
    ghOla.set
    (
        i,
        new volScalarField((g & mesh.C()) - ghRef)
    );
    //surfaceScalarField ghf("ghf", (g & mesh.Cf()) - ghRef);
    ghfOla.set
    (
        i,
        new surfaceScalarField((g & mesh.Cf()) - ghRef)
    );
    
    volScalarField& gh = ghOla[i];
    
    pOla.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "p",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            p_rgh + rho*gh
        )
    );
    volScalarField& p = pOla[i];

    pRefCellOla.set(i, new label(0));
    pRefValueOla.set(i, new scalar(0.0));
    
    label& pRefCell = pRefCellOla[i];
    scalar& pRefValue = pRefValueOla[i];
    
    setRefCell
    (
        p,
        p_rgh,
        pimple.dict(),
        pRefCell,
        pRefValue
    );


    if (p_rgh.needReference())
    {
        Info << "p_rgh need Reference" << endl;
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - rho*gh;
    }


    mesh.setFluxRequired(p_rgh.name());
    mesh.setFluxRequired(alpha1.name());

    //#include "createMRF.H"
    MRFOla.set(i, new IOMRFZoneList(mesh));
    
}


#include "createOlaAlphaFluxes.H"


PtrList<fv::options> fvOptionsOla(olaRegions.size());
forAll(olaNames, i)
{
    dynamicFvMesh& mesh = olaRegions[i];
    fvOptionsOla.set(i, new fv::options(mesh));
}

PtrList<volScalarField> rAUOla(olaRegions.size());
List<autoPtr<surfaceVectorField>> UfOla(olaRegions.size());

forAll(olaNames, i)
{
    dynamicFvMesh& mesh = olaRegions[i];
    volScalarField& rho = rhoOla[i];
    volVectorField& U = UOla[i];
    
    rAUOla.set
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

    if (mesh.dynamic())
    {
        Info<< "Constructing face velocity Uf\n" << endl;

        UfOla[i].reset
        (
            new surfaceVectorField
            (
                IOobject
                (
                    "Uf",
                    runTime.timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                fvc::interpolate(U)
            )
        );
    }
}









