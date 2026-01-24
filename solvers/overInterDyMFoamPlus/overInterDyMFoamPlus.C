/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2014 OpenFOAM Foundation
    Copyright (C) 2016-2021 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    overInterDyMFoamPlus

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "fvcSmooth.H"
#include "cellCellStencilObject.H"
#include "localMin.H"
#include "localMax.H"
#include "oversetAdjustPhi.H"
#include "oversetPatchPhiErr.H"
#include "regionProperties.H"
#include "bulletOversetFvMesh.H"
#include "bulletMotionSolver.H"

#include "physicsManager.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Solver for two incompressible, isothermal immiscible fluids using"
        " VOF phase-fraction based interface capturing\n"
        "With optional mesh motion and mesh topology changes including"
        " adaptive re-meshing."
    );
    
    #include "initOIDF.H"
    #include "initPhysics.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;
        
    while (runTime.run())
    {                
        #include "multiregionCourantNo.H"
        #include "multiregionAlphaCourantNo.H"
        #include "setDeltaT.H"
        
        
        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;
        
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
                                talphaPhi1Corr0.clear();
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
                        #include "alphaEqnSubCycle.H"
                                    
                        rhoPhi *= faceMask;

                        mixture.correct();
                        
                        #include "UEqn.H"

                        // --- Pressure corrector loop
                        while (pimple.correct())
                        {
                            #include "pEqn.H"
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
        
        runTime.write();
        
        PhysicsManager.writeAndStoreStates();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
