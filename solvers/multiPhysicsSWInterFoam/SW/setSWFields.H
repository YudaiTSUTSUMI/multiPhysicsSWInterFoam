fvMesh& mesh = SWRegions[i];

bool limitFrontVelocity = limitFrontVelocitySW[i];
dimensionedScalar& hDry = hDrySW[i];

volScalarField& h = hSW[i];
volScalarField& h0 = h0SW[i];
volScalarField& hTotal = hTotalSW[i];

volVectorField& hU = hUSW[i];
volVectorField& U = USW[i];

volScalarField& wetdry = wetdrySW[i];

SWRiemannSolver& riemannSolver = riemannSolverSW[i];

SWWaveSource& waveSource = waveSourceSW[i];

