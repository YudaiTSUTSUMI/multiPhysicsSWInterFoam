/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2016-2023 OpenCFD Ltd.
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
    Foam::bulletMotion
    
\*---------------------------------------------------------------------------*/

#ifndef bulletMotion_H
#define bulletMotion_H

#include "pointField.H"
#include "Tuple2.H"

#include "Time.H"
#include "dictionary.H"
#include "autoPtr.H"
#include "vector.H"
#include "runTimeSelectionTables.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class bulletSolver;

/*---------------------------------------------------------------------------*\
                      Class bulletMotion Declaration
\*---------------------------------------------------------------------------*/

class bulletMotion
{
    friend class bulletSolver;

    // Private data

        //- Reference to time database
        const Time& time_;

        //- Centre of rotation of initial state
        point initialCentreOfRotation_;
        
        point centreOfRotation_;

        //- Orientation of initial state
        tensor initialQ_;
        
        tensor Q_;

        //- Mass of the body#include "Time.H"
#include "dictionary.H"
#include "autoPtr.H"
#include "vector.H"
#include "runTimeSelectionTables.H"

    // Private Member Functions

        
        // Access functions retained as private because of the risk of
        // confusion over what is a body local frame vector and what is global

        // Access

            //- Return the initial centre of rotation
            inline const point& initialCentreOfRotation() const;
            
            inline const point& centreOfRotation() const;

            //- Return the initial orientation
            inline const tensor& initialQ() const;

            //- Return the orientation
            inline const tensor& Q() const;
            
public:

    // Constructors

        //- Construct null
        bulletMotion(const Time&);

        //- Construct from constant and state dictionaries
        bulletMotion
        (
            const dictionary& dict,
            const dictionary& stateDict,
            const Time& time
        );

        //- Construct as copy
        bulletMotion(const bulletMotion&);


    //- Destructor
    ~bulletMotion();


    // Member Functions

        // Access

            //- Return time
            inline const Time& time() const;

        // Transformations

            //- Transform the given initial state point by the current motion
            //  state
            inline point transform(const point& initialPoints) const;

            //- Transform the given initial state pointField by the current
            //  motion state
            tmp<pointField> transform(const pointField& initialPoints) const;

            //- Transform the given initial state pointField by the current
            //  motion state scaled by the given scale
            tmp<pointField> transform
            (
                const pointField& initialPoints,
                const scalarField& scale
            ) const;

        
        // Edit

            //- Return the centre of rotation
            inline point& initialCentreOfRotation();
            
            //- Return non-const access to the centre of rotation
            inline point& centreOfRotation();

            //- Return initial orientation
            inline tensor& initialQ();

            //- Return non-const access to the orientation
            inline tensor& Q();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "bulletMotionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
