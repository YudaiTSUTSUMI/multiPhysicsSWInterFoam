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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline const Foam::point&
Foam::bulletMotion::initialCentreOfRotation() const
{
    return initialCentreOfRotation_;
}

inline const Foam::point&
Foam::bulletMotion::centreOfRotation() const
{
    return centreOfRotation_;
}


inline const Foam::tensor&
Foam::bulletMotion::initialQ() const
{
    return initialQ_;
}

inline const Foam::tensor&
Foam::bulletMotion::Q() const
{
    return Q_;
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

inline Foam::point& Foam::bulletMotion::initialCentreOfRotation()
{
    return initialCentreOfRotation_;
}


inline Foam::point& Foam::bulletMotion::centreOfRotation()
{
    return centreOfRotation_;
}


inline Foam::tensor& Foam::bulletMotion::initialQ()
{
    return initialQ_;
}


inline Foam::tensor& Foam::bulletMotion::Q()
{
    return Q_;
}

inline const Foam::Time& Foam::bulletMotion::time() const
{
    return time_;
}

inline Foam::point Foam::bulletMotion::transform
(
    const point& initialPoint
) const
{
    return
    (
        centreOfRotation()
      + (Q() & initialQ_.T() & (initialPoint - initialCentreOfRotation_))
    );
}

// ************************************************************************* //
