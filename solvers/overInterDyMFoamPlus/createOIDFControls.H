PtrList<bool> correctPhiOIDF(OIDFRegions.size());
PtrList<bool> checkMeshCourantNoOIDF(OIDFRegions.size());

forAll(OIDFNames, i)
{
    dynamicFvMesh& mesh = OIDFRegions[i];
    pimpleControl& pimple = pimpleOIDF[i];
    
    correctPhiOIDF.set
    (
        i,
        new bool
        (
            pimple.dict().getOrDefault("correctPhi", mesh.dynamic())
        )
    );
    
    checkMeshCourantNoOIDF.set
    (
        i,
        new bool
        (
            pimple.dict().getOrDefault("checkMeshCourantNo", false)
        )
    );
}


PtrList<bool> adjustFringeOIDF(OIDFRegions.size());
PtrList<label> zoneIdMassOIDF(OIDFRegions.size());
PtrList<bool> oversetPatchErrOutputOIDF(OIDFRegions.size());

forAll(OIDFNames, i)
{
    pimpleControl& pimple = pimpleOIDF[i];
    
    adjustFringeOIDF.set
    (
        i,
        new bool
        (
            pimple.dict().getOrDefault("oversetAdjustPhi", false)
        )
    );
    
    zoneIdMassOIDF.set
    (
        i,
        new label
        (
            pimple.dict().getOrDefault("adjustZone", -1)
        )
    );
    
    oversetPatchErrOutputOIDF.set
    (
        i,
        new bool
        (
            pimple.dict().getOrDefault("oversetPatchErrOutput", false)
        )
    );
}

