PtrList<bool> correctPhiOla(olaRegions.size());
PtrList<bool> checkMeshCourantNoOla(olaRegions.size());

forAll(olaNames, i)
{
    dynamicFvMesh& mesh = olaRegions[i];
    pimpleControl& pimple = pimpleOla[i];
    
    correctPhiOla.set
    (
        i,
        new bool
        (
            pimple.dict().getOrDefault("correctPhi", mesh.dynamic())
        )
    );
    
    checkMeshCourantNoOla.set
    (
        i,
        new bool
        (
            pimple.dict().getOrDefault("checkMeshCourantNo", false)
        )
    );
}
