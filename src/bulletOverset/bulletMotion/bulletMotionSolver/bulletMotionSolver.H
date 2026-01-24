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

Class
    Foam::bulletMotionSolver
\*---------------------------------------------------------------------------*/

#ifndef bulletMotionSolver_H
#define bulletMotionSolver_H

#include "displacementMotionSolver.H"
#include "bulletMotion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class bulletMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class bulletMotionSolver
:
    public displacementMotionSolver
{
    // Private data

        //- Six DoF motion object
        bulletMotion motion_;

        wordRes patches_;

        //- Patches to integrate forces
        const labelHashSet patchSet_;

        //- Inner morphing distance (limit of solid-body region)
        const scalar di_;

        //- Outer morphing distance (limit of linear interpolation region)
        const scalar do_;

        //- Test-mode in which only the gravitational body-force is applied
        bool test_;

        //- Reference density required by the forces object for
        //  incompressible calculations, required if rho == rhoInf
        scalar rhoInf_;

        //- Name of density field, optional unless used for an
        //  incompressible simulation, when this needs to be specified
        //  as rhoInf
        word rhoName_;

        //- Current interpolation scale (1 at patches, 0 at distance_)
        pointScalarField scale_;
        
        bool multiOversetMesh_;
        pointScalarField oversetScale_;

        //- Current time index (used for updating)
        label curTimeIndex_;

        //- Name of the uniformVectorField for CofG displacement
        word cOfGdisplacement_;
        
        dictionary solverDict_;
        
        vector force_;
        
        vector moment_;


    // Private Member Functions

        //- No copy construct
        bulletMotionSolver
        (
            const bulletMotionSolver&
        ) = delete;

        //- No copy assignment
        void operator=(const bulletMotionSolver&) = delete;


public:

    //- Runtime type information
    TypeName("bulletMotion");


    // Constructors

        //- Construct from polyMesh and IOdictionary
        bulletMotionSolver
        (
            const polyMesh&,
            const IOdictionary& dict
        );


    //- Destructor
    ~bulletMotionSolver() = default;


    // Member Functions

        //- Return the six DoF motion object
        const bulletMotion& motion() const;

        //- Return point location obtained from the current motion field
        virtual tmp<pointField> curPoints() const;
        
        virtual tmp<pointField> newPoints();

        //- Solve for motion
        virtual void solve(){}
       
        // new function
        void calcForceMoment  
        (
            const vector& worldCoM
        );  
        
        void setOversetScale
        (
            const pointScalarField& oversetScale
        );
        
        //- Initialize motion state due to compund body shape
        void initializeMotionState
        (
            const vector initialCentreOfRotation,
            const tensor initialQ
        );
        
        //- Update motion state from Bullet
        void updateMotionState
        (
            const vector P,
            const tensor Q
        );
        

        //- Return access
            //- Return force
            const vector& force() const;

            //- Return moment
            const vector& moment() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
