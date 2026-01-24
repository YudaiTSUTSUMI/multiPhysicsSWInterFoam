/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2017 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "bulletMotionSolver.H"
#include "addToRunTimeSelectionTable.H"
#include "polyMesh.H"
#include "pointPatchDist.H"
#include "pointConstraints.H"
#include "uniformDimensionedFields.H"
#include "forces.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(bulletMotionSolver, 0);

    addToRunTimeSelectionTable
    (
        motionSolver,
        bulletMotionSolver,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::bulletMotionSolver::bulletMotionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict
)
:
    displacementMotionSolver(mesh, dict, typeName),
    motion_
    (
        coeffDict(),
        IOobject
        (
            "bulletMotionState",
            mesh.time().timeName(),
            "uniform",
            mesh
        ).typeHeaderOk<IOdictionary>(true)
      ? IOdictionary
        (
            IOobject
            (
                "bulletMotionState",
                mesh.time().timeName(),
                "uniform",
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            )
        )
      : coeffDict(),
        mesh.time()
    ),
    patches_(coeffDict().get<wordRes>("patches")),
    patchSet_(mesh.boundaryMesh().patchSet(patches_)),
    di_(coeffDict().get<scalar>("innerDistance")),
    do_(coeffDict().get<scalar>("outerDistance")),
    test_(coeffDict().getOrDefault("test", false)),
    rhoInf_(1.0),
    rhoName_(coeffDict().getOrDefault<word>("rho", "rho")),
    scale_
    (
        IOobject
        (
            "motionScale",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        pointMesh::New(mesh),
        dimensionedScalar(dimless, Zero)
    ),
    multiOversetMesh_(false),
    oversetScale_
    (
        IOobject
        (
            "oversetScale",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        pointMesh::New(mesh),
        dimensionedScalar(dimless, Zero)
    ),
    curTimeIndex_(-1),
    cOfGdisplacement_(coeffDict().getOrDefault<word>("cOfGdisplacement", "none")),
    solverDict_(dict),
    force_(Zero),
    moment_(Zero)
{
    if (rhoName_ == "rhoInf")
    {
        coeffDict().readEntry("rhoInf", rhoInf_);
    }

    // Calculate scaling factor everywhere
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    {
        const pointMesh& pMesh = pointMesh::New(mesh);

        pointPatchDist pDist(pMesh, patchSet_, points0());

        // Scaling: 1 up to di then linear down to 0 at do away from patches
        scale_.primitiveFieldRef() =
            min
            (
                max
                (
                    (do_ - pDist.primitiveField())/(do_ - di_),
                    scalar(0)
                ),
                scalar(1)
            );

        // Convert the scale function to a cosine
        scale_.primitiveFieldRef() =
            min
            (
                max
                (
                    0.5
                  - 0.5
                   *cos(scale_.primitiveField()
                   *Foam::constant::mathematical::pi),
                    scalar(0)
                ),
                scalar(1)
            );
        
        // when multi overset, need multi scale
        //scale_.write();
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


const Foam::bulletMotion&
Foam::bulletMotionSolver::motion() const
{
    return motion_;
}


Foam::tmp<Foam::pointField>
Foam::bulletMotionSolver::curPoints() const
{
    tmp<pointField> newPoints
    (
        points0() + pointDisplacement_.primitiveField()
    );
    
    if (!moveAllCells())
    {
        tmp<pointField> ttransformedPts(new pointField(mesh().points()));
        pointField& transformedPts = ttransformedPts.ref();

        UIndirectList<point>(transformedPts, pointIDs()) =
            pointField(newPoints.ref(), pointIDs());

        return ttransformedPts;
    }

    return newPoints;
}


Foam::tmp<Foam::pointField> Foam::bulletMotionSolver::newPoints()
{    
    pointDisplacement_.primitiveFieldRef() =
    motion_.transform(points0(), scale_) - points0();
    
    // Displacement has changed. Update boundary conditions
    pointConstraints::New
    (
        pointDisplacement_.mesh()
    ).constrainDisplacement(pointDisplacement_);
    
    return curPoints();
}



// This is for eccentric object 
void Foam::bulletMotionSolver::calcForceMoment
(
    const vector& worldCoM
)
{
    if (mesh().nPoints() != points0().size())
    {
        FatalErrorInFunction
            << "The number of points in the mesh seems to have changed." << endl
            << "In constant/polyMesh there are " << points0().size()
            << " points; in the current mesh there are " << mesh().nPoints()
            << " points." << exit(FatalError);
    }
    
    
    // Store the motion state at the beginning of the time-step
    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        curTimeIndex_ = this->db().time().timeIndex();
    }

    dimensionedVector g("g", dimAcceleration, Zero);

    if (db().time().foundObject<uniformDimensionedVectorField>("g"))
    {
        g = db().time().lookupObject<uniformDimensionedVectorField>("g");
    }
    else
    {
        coeffDict().readIfPresent("g", g);
    }

    dictionary forcesDict;
    
    point worldCoMPoint(worldCoM);

    forcesDict.add("type", functionObjects::forces::typeName);
    forcesDict.add("patches", patches_);
    forcesDict.add("rhoInf", rhoInf_);
    forcesDict.add("rho", rhoName_);
    forcesDict.add("CofR", worldCoMPoint);//change!!motion_.centreOfRotation()

    functionObjects::forces f("forces", db(), forcesDict);

    f.calcForcesMoments();//
    
    force_ = f.forceEff();
    moment_ = f.momentEff();
}


void Foam::bulletMotionSolver::setOversetScale
(
    const pointScalarField& oversetScale
)
{
    Info  << "setting oveset scale" << nl;
    multiOversetMesh_ = true;
    oversetScale_ = oversetScale;
    const pointMesh& pMesh = pointMesh::New(mesh());
    pointConstraints::New(pMesh).constrain(oversetScale_);
}


void Foam::bulletMotionSolver::initializeMotionState
(
    const vector initialCentreOfRotation,
    const tensor initialQ
)
{
    motion_.initialCentreOfRotation() = initialCentreOfRotation;
    motion_.initialQ() = initialQ;
}


void Foam::bulletMotionSolver::updateMotionState
(
    const vector P,
    const tensor Q
)
{
    motion_.centreOfRotation() = P;
    motion_.Q() = Q;
}


const Foam::vector& Foam::bulletMotionSolver::force() const
{
    return force_;
}


const Foam::vector& Foam::bulletMotionSolver::moment() const
{
    return moment_;
}



// ************************************************************************* //
