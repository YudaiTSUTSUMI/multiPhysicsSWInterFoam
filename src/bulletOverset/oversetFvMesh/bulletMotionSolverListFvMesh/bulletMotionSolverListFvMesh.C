/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
    Copyright (C) 2019-2022 OpenCFD Ltd.
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

#include "bulletMotionSolverListFvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "bulletMotionSolver.H"
#include "pointMesh.H"
#include "pointConstraints.H"
#include "volFields.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(bulletMotionSolverListFvMesh, 0);
    addToRunTimeSelectionTable
    (
        dynamicFvMesh,
        bulletMotionSolverListFvMesh,
        IOobject
    );
    addToRunTimeSelectionTable
    (
        dynamicFvMesh,
        bulletMotionSolverListFvMesh,
        doInit
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::bulletMotionSolverListFvMesh::bulletMotionSolverListFvMesh
(
    const IOobject& io,
    const bool doInit
)
:
    dynamicFvMesh(io, doInit),
    bulletMotionSolvers_()
{
    if (doInit)
    {
        init(false);    // do not initialise lower levels
    }
}


bool Foam::bulletMotionSolverListFvMesh::init
(
    const bool doInit,
    const bool mandatory
)
{
    if (doInit)
    {
        dynamicFvMesh::init(doInit);
    }

    IOobject ioDict
    (
        "dynamicMeshDict",
        time().constant(),
        *this,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        IOobject::NO_REGISTER
    );

    IOdictionary dict(ioDict);

    label i = 0;

    const auto* dictptr = dict.findDict("solvers");
    if (dictptr)
    {
        const dictionary& solverDict = *dictptr;

        bulletMotionSolvers_.setSize(solverDict.size());
        oversetScaleList_.setSize(solverDict.size());
        
        bool multiOversetMesh = false;
        if(solverDict.size() >= 2)
        {
            multiOversetMesh = true;
            Info << "multiOversetMesh" << nl;
        }

        for (const entry& dEntry : solverDict)
        {
            if (dEntry.isDict())
            {
                IOobject io(ioDict);
                io.readOpt(IOobject::NO_READ);
                io.writeOpt(IOobject::AUTO_WRITE);
                io.rename(dEntry.dict().dictName());

                IOdictionary IOsolverDict
                (
                    io,
                    dEntry.dict()
                );
                
                
                bulletMotionSolver* solver = new bulletMotionSolver(*this, IOsolverDict); // create pointer
                
                if(multiOversetMesh)
                {
                    pointScalarField oversetScale = setOversetScale(IOsolverDict, i);
                }
                
                bulletMotionSolvers_.set
                (
                    i,
                    solver //function for bulletMotionSolver
                );
                
                solverDictList_.append(dEntry.dict());
                
                i++;
            }
        }
        bulletMotionSolvers_.setSize(i);
    }
    else if (mandatory)
    {
        bulletMotionSolvers_.setSize(1);
        bulletMotionSolver* solver = new bulletMotionSolver(*this, dict);
        bulletMotionSolvers_.set(i++, solver);
    }
    
    // Assume something changed
    return true;
}



bool Foam::bulletMotionSolverListFvMesh::init(const bool doInit)
{
    // Fall-back to always constructing bulletMotionSolver
    return init(doInit, true);
}


Foam::pointScalarField Foam::bulletMotionSolverListFvMesh::setOversetScale
(
    const IOdictionary solverDict,
    const int i
)
{
    fvMesh& mesh = *this;
    
    volScalarField zoneIDList
    (
        IOobject
        (
            "zoneID",
            time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );
    
    const label zoneID = solverDict.get<label>("zoneID");
    
    oversetScaleList_.set
    (
        i,
        new pointScalarField
        (
            IOobject
            (
                "oversetScaleForZone" + Foam::name(zoneID),
                time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            pointMesh::New(mesh),
            dimensionedScalar(dimless, Zero)
        )
    );
    pointScalarField& oversetScale = oversetScaleList_[i];
    
    forAll(mesh.cells(), celli)
    {
        if(zoneIDList[celli] == zoneID)
        {
            forAll(mesh.cellPoints()[celli], pointi)
            {
                label cellpID = mesh.cellPoints()[celli][pointi];
                oversetScale[cellpID] = 1;
            }
        }
    }
    
    oversetScale.write();
    
    return oversetScaleList_[i];
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::bulletMotionSolverListFvMesh::~bulletMotionSolverListFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::bulletMotionSolverListFvMesh::mapFields
(
    const mapPolyMesh& mpm
)
{
    dynamicFvMesh::mapFields(mpm);

    // Update the bulletMotionSolvers for any topo change ...
    for (auto& ms : bulletMotionSolvers_)
    {
        ms.updateMesh(mpm);
    }
}


bool Foam::bulletMotionSolverListFvMesh::update()
{
    Info << nl << "bulletMotionSolverListFvMesh::update()" << nl;
    if (bulletMotionSolvers_.size())
    {
        // Accumulated displacement
        pointField disp(0*fvMesh::points()); 
        
        bool multiOversetMesh = false;
        if(bulletMotionSolvers_.size() >= 2)
        {
            multiOversetMesh = true;
        }
        
        for (label i = 0; i < bulletMotionSolvers_.size(); i++)
        {   
            if(multiOversetMesh)
            {
                Info << "multiOversetMesh, No." << i << nl;
                tmp<pointScalarField> scale
                (
                    oversetScaleList_[i]
                );
                
                scalar movingPoints = gSum(oversetScaleList_[i]);
                Info << "moving " << movingPoints << " points" << nl;
            
                disp += scale*(bulletMotionSolvers_[i].newPoints() - fvMesh::points());
            }
            else
            {
                disp += bulletMotionSolvers_[i].newPoints() - fvMesh::points();
            }
        }
        
        fvMesh::movePoints(points() + disp);
    }
    Info << "finish move points" << endl;
    return true;
}


// ************************************************************************* //
