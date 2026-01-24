PtrList<volScalarField> p_rghOIDF(OIDFRegions.size());
PtrList<volVectorField> UOIDF(OIDFRegions.size());
PtrList<surfaceScalarField> phiOIDF(OIDFRegions.size());

PtrList<volScalarField> cellMaskOIDF(OIDFRegions.size());
PtrList<surfaceScalarField> faceMaskOIDF(OIDFRegions.size());
PtrList<volScalarField> interpolatedCellsOIDF(OIDFRegions.size());

PtrList<volScalarField> rhoOIDF(OIDFRegions.size());
PtrList<surfaceScalarField> rhoPhiOIDF(OIDFRegions.size());

PtrList<scalar> massLiquidOIDF(OIDFRegions.size());

PtrList<volScalarField> ghOIDF(OIDFRegions.size());
PtrList<surfaceScalarField> ghfOIDF(OIDFRegions.size());

PtrList<volScalarField> pOIDF(OIDFRegions.size());

PtrList<immiscibleIncompressibleTwoPhaseMixture> mixtureOIDF(OIDFRegions.size());

PtrList<dimensionedScalar> ghRefOIDF(OIDFRegions.size());

PtrList<label> pRefCellOIDF(OIDFRegions.size());
PtrList<scalar> pRefValueOIDF(OIDFRegions.size());

PtrList<incompressible::turbulenceModel> turbulenceOIDF(OIDFRegions.size());

PtrList<IOMRFZoneList> MRFOIDF(OIDFRegions.size());

PtrList<volScalarField> boundaryCellsOIDF(OIDFRegions.size());
PtrList<labelList> boundaryCellsListOIDF(OIDFRegions.size());

forAll(OIDFNames, i)
{
    Info<< "*** Region " << OIDFRegions[i].name() << nl << endl;
    dynamicFvMesh& mesh = OIDFRegions[i];
    pimpleControl& pimple = pimpleOIDF[i];
    
    #include "createRDeltaT.H"

    Info<< "Reading field p_rgh\n" << endl;
    p_rghOIDF.set
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
    volScalarField& p_rgh = p_rghOIDF[i];

    Info<< "Reading field U\n" << endl;
    UOIDF.set
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
    volVectorField& U = UOIDF[i];

    //#include "createPhi.H"
    Info<< "Reading/calculating face flux field phi\n" << endl;
    phiOIDF.set
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
    surfaceScalarField& phi = phiOIDF[i];
    

    //- Overset specific

    // Add solver-specific interpolations
    {
        wordHashSet& nonInt =
            const_cast<wordHashSet&>(Stencil::New(mesh).nonInterpolatedFields());

        nonInt.insert("HbyA");
        nonInt.insert("grad(p_rgh)");
        nonInt.insert("nHat");
        nonInt.insert("surfaceIntegrate(nHatf)");
        nonInt.insert("surfaceIntegrate(phi+meshPhi)");
        nonInt.insert("surfaceIntegrate(phi)");
        nonInt.insert("surfaceIntegrate(phiHbyA)");
        nonInt.insert("surfaceSum(((S|magSf)*S)");
        nonInt.insert("surfaceIntegrate(((rAUf*magSf)*snGradCorr(p_rgh)))");
        nonInt.insert("cellMask");
        nonInt.insert("cellDisplacement");
        nonInt.insert("interpolatedCells");
        nonInt.insert("cellInterpolationWeight");
        nonInt.insert("pcorr");
        
        nonInt.insert("boundaryCells");
    }


    // Mask field for zeroing out contributions on hole cells
    //#include "createCellMask.H"
    cellMaskOIDF.set
    (
        i, 
        new volScalarField
        (
            IOobject
            (
                "cellMask",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("one", dimless, 1.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    volScalarField& cellMask = cellMaskOIDF[i];
    
    faceMaskOIDF.set
    (
        i, 
        new surfaceScalarField
        (
            localMin<scalar>(mesh).interpolate(cellMask)
        )
    );
    surfaceScalarField& faceMask = faceMaskOIDF[i];
    
    #include "setCellMask.H"
    
    // Create bool field with interpolated cells
    //#include "createInterpolatedCells.H"
    interpolatedCellsOIDF.set
    (
        i, 
        new volScalarField
        (
            IOobject
            (
                "interpolatedCells",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("one", dimless, 1.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    volScalarField& interpolatedCells = interpolatedCellsOIDF[i];
    
    #include "setInterpolatedCells.H"
    

    Info<< "Reading transportProperties\n" << endl;
    mixtureOIDF.set
    (
        i,
        new immiscibleIncompressibleTwoPhaseMixture
        (
            U,
            phi
        )
    );
    immiscibleIncompressibleTwoPhaseMixture& mixture = mixtureOIDF[i];

    volScalarField& alpha1(mixture.alpha1());
    volScalarField& alpha2(mixture.alpha2());

    const dimensionedScalar& rho1 = mixture.rho1();
    const dimensionedScalar& rho2 = mixture.rho2();


    // Need to store rho for ddt(rho, U)
    rhoOIDF.set
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
    volScalarField& rho = rhoOIDF[i];
    rho.oldTime();
    
    massLiquidOIDF.set(i, new scalar(0.0));
    scalar& massLiquid = massLiquidOIDF[i];
    
    forAll(mesh.cells(), celli)
    {
        scalar massElement = rho[celli] * mesh.V()[celli];
        
        massLiquid += massElement;
    }
    reduce(massLiquid, sumOp<scalar>());

    // Mass flux
    rhoPhiOIDF.set
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
    turbulenceOIDF.set
    (
        i,
        incompressible::turbulenceModel::New
        (
            U,
            phi,
            mixture
        )
    );
    
    turbulenceOIDF[i].validate();

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
    ghRefOIDF.set
    (
        i,
        new dimensionedScalar
        (
            mag(g.value()) > SMALL
          ? g & (cmptMag(g.value())/mag(g.value()))*hRef
          : dimensionedScalar("ghRef", g.dimensions()*dimLength, 0)
        )
    );
    dimensionedScalar& ghRef = ghRefOIDF[i];
    
    ghOIDF.set
    (
        i,
        new volScalarField((g & mesh.C()) - ghRef)
    );
    
    ghfOIDF.set
    (
        i,
        new surfaceScalarField((g & mesh.Cf()) - ghRef)
    );
    
    volScalarField& gh = ghOIDF[i];
    
    pOIDF.set
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
    volScalarField& p = pOIDF[i];

    pRefCellOIDF.set(i, new label(0));
    pRefValueOIDF.set(i, new scalar(0.0));
    
    label& pRefCell = pRefCellOIDF[i];
    scalar& pRefValue = pRefValueOIDF[i];
    
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
    MRFOIDF.set(i, new IOMRFZoneList(mesh));
    
}


#include "createOIDFAlphaFluxes.H"


PtrList<fv::options> fvOptionsOIDF(OIDFRegions.size());
forAll(OIDFNames, i)
{
    dynamicFvMesh& mesh = OIDFRegions[i];
    fvOptionsOIDF.set(i, new fv::options(mesh));
}

PtrList<volScalarField> rAUOIDF(OIDFRegions.size());
PtrList<surfaceVectorField> UfOIDF(OIDFRegions.size());
forAll(OIDFNames, i)
{
    dynamicFvMesh& mesh = OIDFRegions[i];
    volScalarField& rho = rhoOIDF[i];
    volVectorField& U = UOIDF[i];
    
    rAUOIDF.set
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
    
    UfOIDF.set
    (
        i,
        new surfaceVectorField
        (
            IOobject
            (
                "Uf",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            fvc::interpolate(U)
        )
    );
}

forAll(OIDFNames, i)
{
    dynamicFvMesh& mesh = OIDFRegions[i];
    volScalarField& cellMask = cellMaskOIDF[i];
    surfaceScalarField& faceMask = faceMaskOIDF[i];
    volScalarField& interpolatedCells = interpolatedCellsOIDF[i];
    
    #include "setCellMask.H"
    #include "setInterpolatedCells.H"
}
