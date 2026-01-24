forAll(olaNames, i)
{   
    Info<< "Calculate region: " << olaNames[i] << nl;
    #include "setOlaFields.H"         

    while (pimple.loop())
    {
        if (pimple.firstIter())
        {
            mesh.update();

            if (mesh.changing())
            {
                // Do not apply previous time-step mesh compression flux
                // if the mesh topology changed
                if (mesh.topoChanging())
                {
                    talphaPhi1Corr0.clear();
                }

                gh = (g & mesh.C()) - ghRef;
                ghf = (g & mesh.Cf()) - ghRef;

                MRF.update();

                if (correctPhi)
                {
                    // Calculate absolute flux
                    // from the mapped surface velocity
                    phi = mesh.Sf() & UfOla[i]();

                    #include "correctPhiOla.H"

                    // Make the flux relative to the mesh motion
                    fvc::makeRelative(phi, U);

                    mixture.correct();
                }

                if (checkMeshCourantNo)
                {
                    #include "meshCourantNo.H"
                }
            }
        }

        #include "alphaControls.H"
        #include "alphaEqnSubCycleOla.H"

        mixture.correct();

        if (pimple.frozenFlow())
        {
            continue;
        }

        #include "UEqnOla.H"

        // --- Pressure corrector loop
        while (pimple.correct())
        {
            #include "pEqnOla.H"
        }

        if (pimple.turbCorr())
        {
            turbulence.correct();
        }
    }

    runTime.write();
    // Write Porous Variables
    if( activePorosity && runTime.outputTime() ) 
    {
        porosity.write();
        porosityIndex.write();
    }
        
}
