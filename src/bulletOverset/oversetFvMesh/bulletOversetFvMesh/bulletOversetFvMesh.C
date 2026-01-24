/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify i
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

#include "bulletOversetFvMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(bulletOversetFvMesh, 0);
    addToRunTimeSelectionTable(dynamicFvMesh, bulletOversetFvMesh, IOobject);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::bulletOversetFvMesh::bulletOversetFvMesh
(
    const IOobject& io,
    const bool doInit
)
:
    bulletMotionSolverListFvMesh(io, doInit),
    oversetFvMeshBase(static_cast<const fvMesh&>(*this))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::bulletOversetFvMesh::~bulletOversetFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::bulletOversetFvMesh::update()
{
    if (!bulletMotionSolverListFvMesh::update())
    {
        return false;
    }

    oversetFvMeshBase::update();

    return true;
}


bool Foam::bulletOversetFvMesh::writeObject
(
    IOstreamOption streamOpt,
    const bool writeOnProc
) const
{
    bool ok =
        bulletMotionSolverListFvMesh::writeObject(streamOpt, writeOnProc);
    ok = oversetFvMeshBase::writeObject(streamOpt, writeOnProc) && ok;
    return ok;
}


// ************************************************************************* //
