PtrList<volScalarField> porosityIndexOla(olaRegions.size());
PtrList<bool> activePorosityOla(olaRegions.size());
PtrList<volScalarField> aPorFieldOla(olaRegions.size());
PtrList<volScalarField> bPorFieldOla(olaRegions.size());
PtrList<volScalarField> cPorFieldOla(olaRegions.size());
PtrList<volScalarField> useTransMaskOla(olaRegions.size());
PtrList<volScalarField> KCPorFieldOla(olaRegions.size());
PtrList<volScalarField> D50FieldOla(olaRegions.size());
PtrList<volScalarField> porosityOla(olaRegions.size());
PtrList<surfaceScalarField> porosityFOla(olaRegions.size());

forAll(olaNames, i)
{
    Info<< "*** Region " << olaRegions[i].name()
        << " for time = " << runTime.timeName() << nl << endl;

    dynamicFvMesh& mesh = olaRegions[i];

    // -------------------------------------------------------------------------
    // Read porosityIndex
    // -------------------------------------------------------------------------
    Info<< "Reading field porosityIndex\n" << endl;
    
    porosityIndexOla.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "porosityIndex",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("porosityIndex", dimless, 0.0)
        )
    );

    volScalarField& porosityIndex = porosityIndexOla[i];

    scalar nPor = gMax(porosityIndex);  // Number of porous materials
    
    activePorosityOla.set(i, new bool(nPor > 0.0));

    // -------------------------------------------------------------------------
    // Create porosity field
    // -------------------------------------------------------------------------
    porosityOla.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "porosity",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("porosity", dimless, 1.0)
        )
    );

    volScalarField& porosity = porosityOla[i];
    
    porosityFOla.set
    (
        i,
        new surfaceScalarField
        (
            IOobject
            (
                "porosityF",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            fvc::interpolate(porosity)
        )
    );

    surfaceScalarField& porosityF = porosityFOla[i];

    // -------------------------------------------------------------------------
    // Allocate local working fields
    // -------------------------------------------------------------------------
    aPorFieldOla.set(i, new volScalarField("aPorField", 0.0*porosityIndex));
    bPorFieldOla.set(i, new volScalarField("bPorField", 0.0*porosityIndex));
    cPorFieldOla.set(i, new volScalarField("cPorField", 0.0*porosityIndex));
    useTransMaskOla.set(i, new volScalarField("useTransMask", 0.0*porosityIndex));
    KCPorFieldOla.set(i, new volScalarField("KCPorField", 1.0 + 0.0*porosityIndex));

    aPorFieldOla[i].dimensions().reset(dimless);
    bPorFieldOla[i].dimensions().reset(dimless);
    cPorFieldOla[i].dimensions().reset(dimless);

    D50FieldOla.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "D50Field",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("D50", dimLength, 1.0)
        )
    );

    volScalarField& aPorField = aPorFieldOla[i];
    volScalarField& bPorField = bPorFieldOla[i];
    volScalarField& cPorField = cPorFieldOla[i];
    volScalarField& useTransMask = useTransMaskOla[i];
    volScalarField& KCPorField = KCPorFieldOla[i];
    volScalarField& D50Field   = D50FieldOla[i];

    // -------------------------------------------------------------------------
    // If porosity is active, read porosityDict and fill fields
    // -------------------------------------------------------------------------
    if (activePorosityOla[i])
    {
        Info<< "Porosity activated\nReading porosity variables\n" << endl;

        IOdictionary porosityDict
        (
            IOobject
            (
                "porosityDict",
                runTime.constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        // Read parameters
        List<scalar> aPor  = porosityDict.lookupOrDefault("a",     List<scalar>(1,0.0));
        List<scalar> bPor  = porosityDict.lookupOrDefault("b",     List<scalar>(1,0.0));
        List<scalar> cPor  = porosityDict.lookupOrDefault("c",     List<scalar>(1,0.0));
        List<scalar> D50Por= porosityDict.lookupOrDefault("D50",   List<scalar>(1,1.0));
        List<scalar> phiPor= porosityDict.lookupOrDefault("porosity",List<scalar>(1,1.0));

        bool debugPor      = porosityDict.lookupOrDefault<bool>("debugPor", false);
        bool useTransient  = porosityDict.lookupOrDefault<bool>("useTransient", false);
        scalar KC          = porosityDict.lookupOrDefault<scalar>("KC", 1.0);

        Info<< "Number of materials: " << nPor+1 << endl;
        Info<< "a        = " << aPor << endl;
        Info<< "b        = " << bPor << endl;
        Info<< "c        = " << cPor << endl;
        Info<< "D50      = " << D50Por << endl;
        Info<< "porosity = " << phiPor << endl;
        Info<< "Transient formulation: " << useTransient << endl;
        if (useTransient) Info<< "KC = " << KC << endl;

        // Consistency check
        if
        (
            aPor.size() != bPor.size() ||
            aPor.size() != cPor.size() ||
            aPor.size() != phiPor.size() ||
            aPor.size() != D50Por.size() ||
            aPor.size() != (nPor+1)
        )
        {
            FatalError
                << "Mismatch in size of aPor/bPor/cPor/phiPor/D50Por\n"
                << "Check porosityDict and porosityIndex"
                << exit(FatalError);
        }

        // Fill fields based on porosityIndex
        forAll(porosityIndex, cellI)
        {
            label idx = label(porosityIndex[cellI]);

            if (idx > 0)
            {
                aPorField[cellI] = aPor[idx];
                bPorField[cellI] = bPor[idx];
                cPorField[cellI] = cPor[idx];

                D50Field[cellI]  = D50Por[idx];
                porosity[cellI]  = phiPor[idx];

                KCPorField[cellI] = KC;

                if (useTransient) useTransMask[cellI] = 1.0;
            }
        }

        porosityF = fvc::interpolate(porosity);

        porosity.write();

        if (debugPor)
        {
            aPorField.write();
            bPorField.write();
            cPorField.write();
            D50Field.write();
        }
    }
    else
    {
        Info << "Porosity NOT activated\n" << endl;
    }
}
