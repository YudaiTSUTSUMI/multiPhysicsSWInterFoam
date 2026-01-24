label iCorrector = 0;
        
while(iCorrector < PhysicsManager.nOuterCorrectors())
{
    if(iCorrector > 0)
    {
        PhysicsManager.restoreStates();
    }
    
    PhysicsManager.update();
    
    forAll(OIDFNames, i)
    {   
        Info<< "Calculate region: " << OIDFNames[i] << nl;
        #include "setOIDFFields.H"            
    
        while (pimple.loop())
        {
                
            if (pimple.firstIter())
            {
                scalar timeBeforeMeshUpdate = runTime.elapsedCpuTime();

                mesh.update();
                

                if (mesh.changing())
                {
                    Info<< "Execution time for mesh.update() = "
                        << runTime.elapsedCpuTime() - timeBeforeMeshUpdate
                        << " s" << endl;

                    // Do not apply previous time-step mesh compression flux
                    // if the mesh topology changed
                    if (mesh.topoChanging())
                    {
                        //talphaPhi1Corr0.clear();
                        talphaPhi1Corr0 = 0.0*talphaPhi1Corr0;
                    }
                    
                    // Update cellMask field for blocking out hole cells
                    #include "setCellMask.H"
                    #include "setInterpolatedCells.H"
                    #include "correctPhiFaceMask.H"
                    
                    gh = (g & mesh.C()) - ghRef;
                    ghf = (g & mesh.Cf()) - ghRef;

                    // Make the flux relative to the mesh motion
                    fvc::makeRelative(phi, U);
                }

                if (mesh.changing() && checkMeshCourantNo)
                {
                    #include "meshCourantNo.H"
                }
            }
            
            
            if (adjustFringe)
            {
                oversetAdjustPhi(phi, U, zoneIdMass);
            }
            
            if(updateVOFOIDF[i])
            {
                #include "alphaControls.H"
                #include "alphaEqnSubCycleOIDF.H"
                            
                rhoPhi *= faceMask;

                mixture.correct();
                
                #include "UEqnOIDF.H"

                // --- Pressure corrector loop
                while (pimple.correct())
                {
                    #include "pEqnOIDF.H"
                }
                
            
                if (pimple.turbCorr())
                {
                    turbulence.correct();
                }
            }
        }
    }
    
    iCorrector++;
}
