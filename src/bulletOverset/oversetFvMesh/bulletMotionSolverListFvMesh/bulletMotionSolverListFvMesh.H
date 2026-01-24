/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
    Copyright (C) 2020,2022 OpenCFD Ltd.
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

Class
    Foam::bulletMotionSolverListFvMesh

SourceFiles
    bulletMotionSolverListFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef bulletMotionSolverListFvMesh_H
#define bulletMotionSolverListFvMesh_H

#include "dynamicFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class bulletMotionSolver;

/*---------------------------------------------------------------------------*\
                Class bulletMotionSolverListFvMesh Declaration
\*---------------------------------------------------------------------------*/

class bulletMotionSolverListFvMesh
:
    public dynamicFvMesh
{
    // Private data

        PtrList<bulletMotionSolver> bulletMotionSolvers_;
        
        //- points to move per zoneID for oversetMesh
        PtrList<pointScalarField> oversetScaleList_;
        
        List<dictionary> solverDictList_;

    // Private Member Functions

        //- No copy construct
        bulletMotionSolverListFvMesh
        (
            const bulletMotionSolverListFvMesh&
        ) = delete;

        //- No copy assignment
        void operator=(const bulletMotionSolverListFvMesh&) = delete;


public:

    //- Runtime type information
    TypeName("bulletMotionSolverListFvMesh");


    // Constructors

        //- Construct from IOobject
        bulletMotionSolverListFvMesh
        (
            const IOobject& io,
            const bool doInit=true
        );


    //- Destructor
    ~bulletMotionSolverListFvMesh();


    // Member Functions

        //- Initialise all non-demand-driven data
        virtual bool init(const bool doInit);

        //- Initialise all non-demand-driven data. Used for initialising
        //  without default bulletMotionSolver (used by dynamicRefineFvMesh)
        virtual bool init(const bool doInit, const bool mandatory);
        
        pointScalarField setOversetScale(const IOdictionary solverDict, const int i);

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();

        //- Map all fields in time using given map. Triggered by topo change
        virtual void mapFields(const mapPolyMesh& mpm);
        
        void updateState();
        
        dictionary getSolverDict(int i)
        {
            return solverDictList_[i];
        }
        
    // Edit
        PtrList<bulletMotionSolver>& bulletMotionSolvers()
        {
            return bulletMotionSolvers_;
        }
        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
