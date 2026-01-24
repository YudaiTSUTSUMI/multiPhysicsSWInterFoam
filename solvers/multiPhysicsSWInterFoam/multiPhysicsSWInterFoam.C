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
    multiPhysicsSWInterFoam

Description
    Multi-physics solver extended from interFoam, incorporating multiple models
    to simulate free-surface flows, wave-structure interactions, and porous
    media dynamics.

    The solver integrates the following components:

    1) 2D Shallow Water Equation (SW):
       - Calculates wave propagation in shallow regions

    2) 2D Laplace Equation (laplacian):
       - Solves pressure distribution in narrow micro gaps

    3) overInterDyMFoam+ (OIDF):
       - Performs FSI simulations using overset mesh method
       - Models jointed and moored structures with internal fluid

    4) olaFlow (based on Higuera et al., 2014) (ola):
       - Supports advanced wave generation
       - Simulates flow through porous structures such as breakwaters

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
#include "simpleControl.H"
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

#include "SWRiemannSolver.H"
#include "regionConnector.H"
#include "SWWaveSource.H"

#include "physicsManager.H"

#include "CorrectPhi.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    //#include "postProcess.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMeshes.H"
    #include "createTimeControls.H"
    #include "createControls.H"
    #include "createFields.H"
    
    #include "multiregionCourantNo.H"
    #include "setInitialDeltaT.H"
        
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
        
        regionConnector.updateSWVOF();
        
        #include "solveSW.H"
        
        regionConnector.updateSWVOF();
                
        #include "solveOIDF.H"
        
        #include "solveOla.H"        
        
        regionConnector.updateVOFLaplacian();
        
        #include "solveLaplacian.H"  
        
        runTime.write();
        
        PhysicsManager.writeAndStoreStates();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
