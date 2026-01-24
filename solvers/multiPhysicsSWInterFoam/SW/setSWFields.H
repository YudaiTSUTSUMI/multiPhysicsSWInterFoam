fvMesh& mesh = SWRegions[i];

volScalarField& h          = hSW[i];
volScalarField& h0         = h0SW[i];
volScalarField& hTotal     = hTotalSW[i];

volVectorField& hU         = hUSW[i];
volVectorField& U          = USW[i];

volScalarField& wetdry     = wetdrySW[i];
surfaceScalarField& wetdryf= wetdryfSW[i];

surfaceScalarField& lsHeight = lsHeightSW[i];
surfaceScalarField& lsLabel  = lsLabelSW[i];

SWRiemannSolver& SWRiemannSolver          = SWRiemannSolverSW[i];

SWWaveSource& waveSource   = waveSourceSW[i];

