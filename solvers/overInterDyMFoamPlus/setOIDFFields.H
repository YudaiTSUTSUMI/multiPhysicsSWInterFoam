dynamicFvMesh& mesh = OIDFRegions[i];
pimpleControl& pimple = pimpleOIDF[i];

scalar& cumulativeContErr = cumulativeContErrOIDF[i];

volScalarField& p_rgh = p_rghOIDF[i];
volVectorField& U = UOIDF[i];
surfaceScalarField& phi = phiOIDF[i];
volScalarField& rho = rhoOIDF[i];
surfaceScalarField& rhoPhi = rhoPhiOIDF[i];

incompressible::turbulenceModel& turbulence = turbulenceOIDF[i];

dimensionedScalar& ghRef = ghRefOIDF[i];
volScalarField& gh = ghOIDF[i];
surfaceScalarField& ghf = ghfOIDF[i];

label& pRefCell = pRefCellOIDF[i];
scalar& pRefValue = pRefValueOIDF[i];

volScalarField& p = pOIDF[i];
volScalarField& rAU = rAUOIDF[i];
surfaceVectorField& Uf = UfOIDF[i];

immiscibleIncompressibleTwoPhaseMixture& mixture = mixtureOIDF[i];
const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();
volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

surfaceScalarField& alphaPhi10 = alphaPhi10OIDF[i];
surfaceScalarField& talphaPhi1Corr0 = talphaPhi1Corr0OIDF[i];
surfaceScalarField& alphaPhiUn = alphaPhiUnOIDF[i];
const bool& alphaRestart = alphaRestartOIDF[i];

volScalarField& cellMask = cellMaskOIDF[i];
surfaceScalarField& faceMask = faceMaskOIDF[i];
volScalarField& interpolatedCells = interpolatedCellsOIDF[i];

bool& checkMeshCourantNo = checkMeshCourantNoOIDF[i];

bool& adjustFringe = adjustFringeOIDF[i];
label& zoneIdMass = zoneIdMassOIDF[i];
bool& oversetPatchErrOutput = oversetPatchErrOutputOIDF[i];

IOMRFZoneList& MRF = MRFOIDF[i];

fv::options& fvOptions = fvOptionsOIDF[i];
