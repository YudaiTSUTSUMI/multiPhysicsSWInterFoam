regionProperties rp(runTime);

const wordList OIDFNames(rp["OIDF"]);
PtrList<bulletOversetFvMesh> OIDFRegions(OIDFNames.size());

Info << "Registered OIDFRegions names: " << OIDFNames << endl;

const wordList SWNames(rp["SW"]);
PtrList<fvMesh> SWRegions(SWNames.size());

Info << "Registered SWRegion names: " << SWNames << endl;

const wordList olaNames(rp["ola"]);
PtrList<dynamicFvMesh> olaRegions(olaNames.size());

Info << "Registered olaRegion names: " << olaNames << endl;

const wordList laplacianNames(rp["laplacian"]);
PtrList<fvMesh> laplacianRegions(laplacianNames.size());

Info << "Registered laplacianRegion names: " << laplacianNames << endl;

#include "createOIDFMeshes.H"
#include "createSWMeshes.H"
#include "createOlaMeshes.H"
#include "createLaplacianMeshes.H"


