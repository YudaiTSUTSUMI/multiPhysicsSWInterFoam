dynamicFvMesh& mesh = olaRegions[i];
pimpleControl& pimple = pimpleOla[i];

scalar& cumulativeContErr = cumulativeContErrOla[i];

volScalarField& p_rgh = p_rghOla[i];

volVectorField& U = UOla[i];
surfaceScalarField& phi = phiOla[i];
volScalarField& rho = rhoOla[i];
surfaceScalarField& rhoPhi = rhoPhiOla[i];

incompressible::turbulenceModel& turbulence = turbulenceOla[i];

dimensionedScalar& ghRef = ghRefOla[i];
volScalarField& gh = ghOla[i];
surfaceScalarField& ghf = ghfOla[i];

label& pRefCell = pRefCellOla[i];
scalar& pRefValue = pRefValueOla[i];

volScalarField& p = pOla[i];
volScalarField& rAU = rAUOla[i];

immiscibleIncompressibleTwoPhaseMixture& mixture = mixtureOla[i];
const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();
volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

//IOobject& alphaPhi10Header = alphaPhi10HeaderOla[i];
surfaceScalarField& alphaPhi10 = alphaPhi10Ola[i];
surfaceScalarField& talphaPhi1Corr0 = talphaPhi1Corr0Ola[i];
surfaceScalarField& alphaPhiUn = alphaPhiUnOla[i];
const bool& alphaRestart = alphaRestartOla[i];

bool& correctPhi = correctPhiOla[i];
bool& checkMeshCourantNo = checkMeshCourantNoOla[i];

IOMRFZoneList& MRF = MRFOla[i];

fv::options& fvOptions = fvOptionsOla[i];


volScalarField& porosityIndex = porosityIndexOla[i];
bool& activePorosity = activePorosityOla[i];
volScalarField& aPorField = aPorFieldOla[i];
volScalarField& bPorField = bPorFieldOla[i];
volScalarField& cPorField = cPorFieldOla[i];
volScalarField& useTransMask = useTransMaskOla[i];
volScalarField& KCPorField = KCPorFieldOla[i];
volScalarField& D50Field = D50FieldOla[i];
volScalarField& porosity = porosityOla[i];
surfaceScalarField& porosityF = porosityFOla[i];


