#include "readGravitationalAcceleration.H"
const vector gValue = g.value();
const dimensionedScalar magg = mag(g);

IOdictionary SWProperties
(
    IOobject
    (
        "SWProperties",
        runTime.constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    )
);

bool adjustTimeStep =
     runTime.controlDict().lookupOrDefault("adjustTimeStep", false);
  
 scalar maxCo =
     runTime.controlDict().lookupOrDefault<scalar>("maxCo", 1.0);
  
 scalar maxDeltaT =
     runTime.controlDict().lookupOrDefault<scalar>("maxDeltaT", GREAT);

bool limitFrontVelocity = SWProperties.lookupOrDefault("limitFrontVelocity", false);

dimensionedScalar hDry
(
    "hDry",
    dimLength,
    SWProperties.lookupOrDefault("hDry", SMALL)
);

volScalarField h
(
    IOobject
    (
        "h",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField wetdry
(
     IOobject
     (
         "wetdry",
         runTime.timeName(),
         mesh,
         IOobject::READ_IF_PRESENT,
         IOobject::AUTO_WRITE
     ),
     pos(h-hDry)
);

volScalarField h0
(
    IOobject
    (
        "h0",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volVectorField hU
(
    IOobject
    (
        "hU",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("U", dimensionSet(0,1,-1,0,0,0,0), vector::zero)
);

volScalarField hTotal
(
    IOobject
    (
        "hTotal",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    h+h0
);

volScalarField N
(
    IOobject
    (
        "N",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("N", dimensionSet(0,0,0,0,0,0,0),0.0)
);
    
SWRiemannSolver riemannSolver
(
    h,
    h0,
    hU,
    N,
    wetdry,
    gValue
);

SWWaveSource waveSource
(
    runTime,
    h,
    hU		    
);

