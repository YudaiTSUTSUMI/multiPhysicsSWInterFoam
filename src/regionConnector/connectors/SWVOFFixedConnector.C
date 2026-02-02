#include "SWVOFFixedConnector.H"
#include "IOdictionary.H"
#include "Time.H"
#include "processorFvPatch.H"
#include "syncTools.H"
#include "triSurface.H"
#include "triSurfaceSearch.H"
#include "cellClassification.H"
#include "meshSearch.H"

// * * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //
template<class templateList>
templateList Foam::SWVOFFixedConnector::generateGlobalList
(
    const templateList& list
)
{
    // Number of processors
    label nProcs = Pstream::nProcs();

    // My Number of processor
    label proci = Pstream::myProcNo(); 
        
    // Creating a placeholder list
    templateList globalList;
    List<templateList> shareList(nProcs);
    
    forAll(shareList, i)
    {
        shareList[i] = templateList(0);
    }
    shareList[proci] = list;
         
    Pstream::gatherList(shareList);
    Pstream::scatterList(shareList);
    
    forAll(shareList, listI)
    {
        globalList.append(shareList[listI]);
    }
    
    return globalList;
}


Foam::labelList Foam::SWVOFFixedConnector::generateGlobalIDs
(
    label localSize
)
{
    label nProcs = Pstream::nProcs();
    label myRank = Pstream::myProcNo();

    List<label> allLocalSizes(nProcs, 0);
    allLocalSizes[myRank] = localSize;
    Pstream::gatherList(allLocalSizes);
    Pstream::scatterList(allLocalSizes);

    List<label> offset(nProcs, 0);
    for (label i = 1; i < nProcs; ++i)
    {
        offset[i] = offset[i-1] + allLocalSizes[i-1];
    }

    label start = offset[myRank];

    labelList globalIDs(localSize);
    for (label i = 0; i < localSize; ++i)
    {
        globalIDs[i] = start + i;
    }

    return globalIDs;
}


template<class FieldType>
List<FieldType> Foam::SWVOFFixedConnector::extractFieldSubset
(
    const labelList& donorIDs,
    const UList<FieldType>& globalField
)
{
    List<FieldType> result(donorIDs.size());

    forAll(donorIDs, i)
    {
        result[i] = globalField[donorIDs[i]];
    }

    return result;
}


void Foam::SWVOFFixedConnector::updateConnection
(
    scalar& h,
    scalar& hFraction,
   
    vector& hU,
    scalar& hUFraction,
    
    const scalar& h0,
    
    scalarField& alpha,
    scalarField& alphaFraction,
    
    vectorField& U,
    scalarField& UFraction,
    
    const vectorField& SfVOF,
    const vectorField& CVOF,
    const scalarField& deltaVOF,
    
    const scalar hVOF,
    const vector hUVOF,
    
    const scalar Fr
)
{
    //valueFraction = 0.0 : zeroGradient
    scalar hTotal = h + h0;
    vector USW = hU/stabilise(h, SMALL);
    
    if (Fr >= 1)
    {
        hFraction = 0;
        hUFraction = 0;
        
        forAll(alpha, VOFI)
        {
            alphaFraction[VOFI] = 1;  
            UFraction[VOFI] = 1;
        }            
        
        forAll(alpha, VOFI)
        {            
            scalar z = CVOF[VOFI].z();
            scalar delta = deltaVOF[VOFI];
            
            if(hTotal < z - 0.5*delta)
            {
                alpha[VOFI] = 0;
            }
            else if(hTotal >= z + 0.5*delta)
            {
                alpha[VOFI] = 1;
            }
            else
            {
                scalar alpha2 = (hTotal - z)/delta + 0.5;
                alpha[VOFI] = max(min(alpha2, 1), 0);
            }
            
            if(alpha[VOFI] < 0.01)
            {
                UFraction[VOFI] = 0;
                U[VOFI] *= 0;
            }
            else
            {
                U[VOFI] = alpha[VOFI]*USW;                    
            }
            
        }   
    }
    else if((Fr > -1 && Fr < 1))
    {
        hFraction = 1;
        hUFraction = 0;
        
        forAll(alpha, VOFI)
        {
            alphaFraction[VOFI] = 0;
            UFraction[VOFI] = 1;
        }            
        
        h = hVOF;
        
        forAll(alpha, VOFI)
        {
            if(alpha[VOFI] < 0.01)
            {
                UFraction[VOFI] = 0;
                U[VOFI] *= 0;
            }
            else
            {
               U[VOFI] = alpha[VOFI]*USW;
            }
        }   
    }
    else
    {
        hFraction = 1;
        hUFraction = 1;
        
        forAll(alpha, VOFI)
        {
            alphaFraction[VOFI] = 0;
            UFraction[VOFI] = 0;
        }
        
        h = hVOF;
        hU = hUVOF;
    }
}


void Foam::SWVOFFixedConnector::storeBasicData
(
    const dictionary& regionConnectorEntry,
    const fvMesh& SWRegion,
    const fvMesh& VOFRegion
)
{
    SWRegionData_.SWMesh_ = &SWRegion;
    VOFRegionData_.VOFMesh_ = &VOFRegion;
    
    const fvMesh* SWMesh = SWRegionData_.SWMesh_;
    const fvMesh* VOFMesh = VOFRegionData_.VOFMesh_;
    
    SWRegionData& SWData = SWRegionData_;
    VOFRegionData& VOFData = VOFRegionData_;    
    
    //SWData
    {
        SWData.h_ = const_cast<volScalarField*>(&SWMesh->lookupObject<volScalarField>("h"));
        SWData.hU_ = const_cast<volVectorField*>(&SWMesh->lookupObject<volVectorField>("hU"));
        SWData.h0_ = const_cast<volScalarField*>(&SWMesh->lookupObject<volScalarField>("h0"));
        
        SWData.cellType_ = labelList(SWMesh->nCells(), -1);
        SWData.faceType_ = labelList(SWMesh->nFaces(), 0);    
        
        SWData.hBound_ = const_cast<surfaceScalarField*>(&SWMesh->lookupObject<surfaceScalarField>("hBound"));
        SWData.hUBound_ = const_cast<surfaceVectorField*>(&SWMesh->lookupObject<surfaceVectorField>("hUBound"));        
    }
    
    
    //VOFData
    {
        VOFData.patchName_ = regionConnectorEntry.get<word>("VOFPatch");
        label patchID = VOFMesh->boundaryMesh().findPatchID(VOFData.patchName_);
        if (patchID < 0)
        {
            FatalErrorInFunction
                << "Patch " << VOFData.patchName_ << " not found in SWRegion " << exit(FatalError);
        }
        else
        {
            VOFData.patchID_ = patchID;
        }
        
        volVectorField* Uptr = VOFMesh->getObjectPtr<volVectorField>("U");
        VOFData.U_ = Uptr;
        VOFData.UBound_ = &refCast<mixedFvPatchField<vector>>(Uptr->boundaryFieldRef()[patchID]);
        
        surfaceScalarField* phiptr = VOFMesh->getObjectPtr<surfaceScalarField>("phi");
        VOFData.phiBound_ = &phiptr->boundaryFieldRef()[patchID];
        
        const IOdictionary transDict
        (
            IOobject
            (
                "transportProperties",
                Uptr->time().constant(),
                Uptr->db(), 
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
            )
        );
        const word alphaName = transDict.get<wordList>("phases")[0];
        
        const dictionary alphaDict = transDict.subDict(alphaName);
        
        volScalarField* alphaptr = VOFMesh->getObjectPtr<volScalarField>("alpha." + alphaName);
        VOFData.alpha_ = alphaptr;
        VOFData.alphaBound_ = &refCast<mixedFvPatchField<scalar>>(alphaptr->boundaryFieldRef()[patchID]);
        
        volScalarField* rhoptr = VOFMesh->getObjectPtr<volScalarField>("rho");
        VOFData.rho_ = rhoptr;
        VOFData.rhoBound_ = &rhoptr->boundaryFieldRef()[patchID];
        
        VOFData.globalIDs_ = generateGlobalIDs(VOFRegionData_.UBound_->size());
    }
}


// * * * * * * * * * * * * * * * * Update Functions  * * * * * * * * * * * * * * //


void Foam::SWVOFFixedConnector::createConnection()
{
    const fvMesh& SWMesh = *SWRegionData_.SWMesh_;
    const fvMesh& VOFMesh = *VOFRegionData_.VOFMesh_;
    
    SWRegionData& SWData = SWRegionData_;
    VOFRegionData& VOFData = VOFRegionData_; 
    
    labelList& SWCellType = SWData.cellType_;        
    labelList& SWFaceType = SWData.faceType_;
    
    vectorField VOFCBoundLocal = VOFMesh.boundaryMesh()[VOFData.patchID_].faceCentres();
    vectorField VOFCBoundGlobal = generateGlobalList(VOFCBoundLocal);
    
    //generate donor
    {        
        boundBox SWBB(SWMesh.points(), true);
        
        scalar SWmidZ = (SWBB.min().z() + SWBB.max().z())/2;
        
        vectorField VOFCExBoundLocal = VOFCBoundLocal;
        forAll(VOFCExBoundLocal, faceI)
        {
            vector Sf = VOFMesh.Sf().boundaryField()[VOFData.patchID_][faceI];
            vector normal = Sf/mag(Sf);
            scalar delta = 0.01*sqrt(mag(Sf));
            
            VOFCExBoundLocal[faceI] += delta*normal;
            
            VOFCExBoundLocal[faceI].z() = SWmidZ;
        }
        vectorField VOFCExBoundGlobal = generateGlobalList(VOFCExBoundLocal);
        
        labelList VOFProcIDLocal = labelList(VOFCBoundLocal.size(), Pstream::myProcNo());
        labelList VOFProcIDGlobal = generateGlobalList(VOFProcIDLocal);
        
        labelList VOFFaceIDLocal = labelList(VOFCBoundLocal.size());
        forAll(VOFFaceIDLocal, faceI)
        {
            VOFFaceIDLocal[faceI] = faceI;
        }
        labelList VOFFaceIDGlobal = generateGlobalList(VOFFaceIDLocal);
        
        List<labelList> VOFGlobalDonorLocal;        
        
        meshSearch meshSearch(SWMesh);
        
        forAll(VOFCBoundGlobal, fbi)
        {
            vector VOFC = VOFCBoundGlobal[fbi];
            vector VOFCEx = VOFCExBoundGlobal[fbi];
                        
            label SWFaceI = meshSearch.findNearestFace(VOFCEx);
            
            //label SWCellI = SWMesh.findCell(VOFCEx);   
            label SWCellI = SWMesh.findNearestCell(VOFCEx);    
            {
                scalar dist = mag(VOFCEx - SWMesh.C()[SWCellI]);
                scalar minDist = dist;
                reduce(minDist, minOp<scalar>());

                if(dist != minDist)
                {
                    SWCellI = -1;
                }
            }  
            
            if(SWCellI >= 0)
            {
                labelList VOFDPC(2);
                VOFDPC[0] = VOFProcIDGlobal[fbi];
                VOFDPC[1] = VOFFaceIDGlobal[fbi];
                
                label localFaceIndex = -1;
                forAll(SWData.localFaceIDs_, i)
                {
                    if(SWData.localFaceIDs_[i] == SWFaceI)
                    {
                        localFaceIndex = i;
                        break;
                    }
                }
                
                if(localFaceIndex  == -1)
                {
                    SWData.localFaceIDs_.append(SWFaceI);
                    
                    SWFaceType[SWFaceI] = 1;
                    SWCellType[SWCellI] = 1; //donor for VOF
                    
                    VOFGlobalDonorLocal.append(labelList(1,fbi));        
                }
                else
                {
                    VOFGlobalDonorLocal[localFaceIndex].append(fbi); 
                }
            }
        }
        
        SWData.VOFGlobalDonor_ = generateGlobalList(VOFGlobalDonorLocal);
        
        if(SWData.VOFGlobalDonor_.size() == 0)
        {
            FatalErrorInFunction << "failed to create connection" << exit(FatalError);
        }
        
        SWData.globalIDs_ = generateGlobalIDs(SWData.localFaceIDs_.size());
    }
        
    // SW cell label
    {
        //boundBox faceBB(VOFCBoundGlobal, true);
        boundBox faceBB(VOFMesh.points(), true);
        scalar xMin = faceBB.min().x();
        scalar xMax = faceBB.max().x();
        scalar yMin = faceBB.min().y();
        scalar yMax = faceBB.max().y();
        
        reduce(xMin, minOp<scalar>());
        reduce(xMax, maxOp<scalar>());
        reduce(yMin, minOp<scalar>());
        reduce(yMax, maxOp<scalar>());
        
        faceBB.min().x() = xMin;
        faceBB.max().x() = xMax;
        faceBB.min().y() = yMin;
        faceBB.max().y() = yMax;
        faceBB.min().z() = -GREAT;
        faceBB.max().z() =  GREAT;
        
        DynamicList<label> stack;
         
        forAll(SWMesh.C(), celli) 
        { 
            if(!faceBB.contains(SWMesh.C()[celli]) && SWCellType[celli] == -1)
            { 
                 SWCellType[celli] = 0; 
                 stack.append(celli);
            } 
        }

        while (!stack.empty())
        {
            label cell = stack.last();
            stack.setSize(stack.size() - 1);

            const labelList& faces = SWMesh.cells()[cell];

            forAll(faces, fi)
            {
                label faceI = faces[fi];

                if (faceI < SWMesh.nInternalFaces())
                {
                    label own = SWMesh.owner()[faceI];
                    label nei = SWMesh.neighbour()[faceI];

                    label otherCell = (cell == own ? nei : own);

                    if (SWCellType[otherCell] == -1) 
                    {
                        SWCellType[otherCell] = 0;
                        stack.append(otherCell); 
                    }
                }
            }
        }
            
        forAll(SWFaceType, SWfi)
        {
            if(SWfi < SWMesh.nInternalFaces())
            {
                const label own = SWMesh.owner()[SWfi];
                const label nei = SWMesh.neighbour()[SWfi];
                 
                if(SWCellType[own] == 1 && SWCellType[nei] == -1)
                {
                    SWCellType[nei] = 2;
                }
                else if(SWCellType[own] == -1 && SWCellType[nei] == 1)
                {
                    SWCellType[own] = 2;
                }     
            }                                   
        }
        
        labelList SWnbrType;
        syncTools::swapBoundaryCellList
        (
            SWMesh,
            SWCellType,
            SWnbrType
        );
        
        for (label faceI = SWMesh.nInternalFaces(); faceI < SWMesh.nFaces(); faceI++)
        {
            const label own = SWMesh.faceOwner()[faceI];
            
            label ownType = SWCellType[own];
            label neiType = SWnbrType[faceI-SWMesh.nInternalFaces()];
        
            label patchID = SWMesh.boundaryMesh().whichPatch(faceI);
            const polyPatch& patch = SWMesh.boundaryMesh()[patchID];
            
            if (patch.coupled())
            {
                if(ownType == -1 && neiType == 1)
                {
                    SWCellType[own] = 2;
                }
            }
        }
        
        
        forAll(SWCellType, SWci)
        {
            label cellType = SWCellType[SWci];
            
            if(cellType == -1)
            {
                SWCellType[SWci] = 3;
            }
        }
    }
}


void Foam::SWVOFFixedConnector::updateConnectBoundary()
{
    SWRegionData& SWData = SWRegionData_;
    VOFRegionData& VOFData = VOFRegionData_;
    
    const fvMesh& SWMesh = *SWRegionData_.SWMesh_;
    const fvMesh& VOFMesh = *VOFRegionData_.VOFMesh_;
    
    scalarField alphaLocal = VOFData.alphaBound_->patchInternalField();
    scalarField alphaGlobal = generateGlobalList(alphaLocal);
    scalarField alphaFractionGlobal = scalarField(alphaGlobal.size(), -1);
   
    vectorField ULocal = VOFData.UBound_->patchInternalField();
    vectorField UGlobal = generateGlobalList(ULocal);
    scalarField UFractionGlobal = scalarField(UGlobal.size(), -1);
    
    
    vectorField SfVOFLocal = VOFMesh.boundaryMesh()[VOFData.patchID_].faceAreas();
    vectorField SfVOFGlobal = generateGlobalList(SfVOFLocal);
    
    vectorField CLocal = VOFMesh.boundaryMesh()[VOFData.patchID_].faceCentres();
    vectorField CGlobal = generateGlobalList(CLocal);
    
    scalarField FrGlobal(CGlobal.size());
    
    scalarField minZLocal(CLocal.size());
    scalarField maxZLocal(CLocal.size());
    {
        const polyPatch& patch = VOFMesh.boundaryMesh()[VOFData.patchID_];
        const pointField& points = VOFMesh.points();
        
        forAll(minZLocal, faceI)
        {
            scalar minZ = GREAT;
            scalar maxZ = -GREAT;
            
            const face& f = patch[faceI];                            

            forAll(f, pi)
            {
                scalar z = points[f[pi]].z();

                minZ = min(minZ, z);
                maxZ = max(maxZ, z);
            }
            
            minZLocal[faceI] = minZ;
            maxZLocal[faceI] = maxZ;
        }
    }
    
    scalarField minZGlobal = generateGlobalList(minZLocal);
    scalarField maxZGlobal = generateGlobalList(maxZLocal);
    
    scalarField hLocal = scalarField(SWData.localFaceIDs_.size(), 0.0);
    scalarField h0Local = scalarField(SWData.localFaceIDs_.size(), 0.0);
    vectorField hULocal = vectorField(SWData.localFaceIDs_.size(), vector::zero);
    vectorField SfSWLocal = vectorField(SWData.localFaceIDs_.size(), vector::zero);
    
    labelList& SWCellType = SWData.cellType_;        
    
    forAll(SWData.localFaceIDs_, SWi)//local loop
    {            
        label facei = SWData.localFaceIDs_[SWi];
        
        if(facei < SWMesh.nInternalFaces())
        {
            const label own = SWMesh.owner()[facei];
            const label nei = SWMesh.neighbour()[facei];
            
            label boundCelli = -1;
             
            if(SWCellType[own] == 1)
            {
                boundCelli = own;
                SfSWLocal[SWi] = SWMesh.Sf()[facei];
            }
            else if(SWCellType[nei] == 1)
            {
                boundCelli = nei;
                SfSWLocal[SWi] = -SWMesh.Sf()[facei];
            }
            
            hLocal[SWi] = (*SWData.h_)[boundCelli];
            h0Local[SWi] = (*SWData.h0_)[boundCelli];
            hULocal[SWi] = (*SWData.hU_)[boundCelli];
        }
        else
        {
            label patchID = SWMesh.boundaryMesh().whichPatch(facei);
            const polyPatch& patch = SWMesh.boundaryMesh()[patchID];
            
            label localFaceID = facei - SWMesh.boundaryMesh()[patchID].start();
            
            const labelUList& faceCells = patch.faceCells();
            const label fCelli = faceCells[localFaceID]; 
            
            hLocal[SWi] = (*SWData.h_)[fCelli];
            h0Local[SWi] = (*SWData.h0_)[fCelli];
            hULocal[SWi] = (*SWData.hU_)[fCelli];
            SfSWLocal[SWi] = SWMesh.Sf().boundaryField()[patchID][localFaceID];
        }
    }
    
    scalarField hGlobal = generateGlobalList(hLocal);
    scalarField hFractionGlobal = scalarField(hGlobal.size(), -1);
    
    vectorField hUGlobal = generateGlobalList(hULocal);
    scalarField hUFractionGlobal = scalarField(hUGlobal.size(), -1);
    
    scalarField h0Global = generateGlobalList(h0Local);
    
    vectorField SfSWGlobal = generateGlobalList(SfSWLocal);
    
    
    forAll(SWData.VOFGlobalDonor_, SWI)
    {
        labelList VOFDonor = SWData.VOFGlobalDonor_[SWI];
        
        scalarField alphaList = extractFieldSubset(VOFDonor, alphaGlobal);
        scalarField alphaFractionList = extractFieldSubset(VOFDonor, alphaFractionGlobal);
        vectorField UList = extractFieldSubset(VOFDonor, UGlobal);
        scalarField UFractionList = extractFieldSubset(VOFDonor, UFractionGlobal);
        vectorField CList = extractFieldSubset(VOFDonor, CGlobal);
        vectorField SfVOFList = extractFieldSubset(VOFDonor, SfVOFGlobal);
        scalarField minZList = extractFieldSubset(VOFDonor, minZGlobal);
        scalarField maxZList = extractFieldSubset(VOFDonor, maxZGlobal);
        
        scalarField deltaList = maxZList - minZList;
        scalar donorHeight = max(maxZList) - min(minZList);
        
        scalar hVOF = 0.0;
        vector hUVOF = vector::zero;
        scalar totalArea = 0.0;
        
        forAll(alphaList, i)
        {
            scalar alpha = alphaList[i];
            vector U = UList[i];
            scalar Area = mag(SfVOFList[i]);

            hVOF += alpha * Area;
            hUVOF += alpha * U * Area;
            totalArea += Area;
        }
        
        scalar bottomLength = totalArea / donorHeight;
        
        hVOF /= bottomLength;
        hUVOF /= bottomLength;
        
        vector normal = SfSWGlobal[SWI]/mag(SfSWGlobal[SWI]);
        
        scalar VSW = (hUGlobal[SWI] & normal)/stabilise(hGlobal[SWI], SMALL);
        scalar Fr = VSW/stabilise(Foam::sqrt(mag(g_)*hGlobal[SWI]), SMALL);
        
        updateConnection
        (
            hGlobal[SWI],
            hFractionGlobal[SWI],
           
            hUGlobal[SWI],
            hUFractionGlobal[SWI],
            
            h0Global[SWI],
            
            alphaList,
            alphaFractionList,
            
            UList,
            UFractionList,
            
            SfVOFList,
            CList,
            deltaList,
            
            hVOF,
            hUVOF,
            
            Fr
        );
        
        forAll(VOFDonor, fi)
        {
            label globalID = VOFDonor[fi];
            
            alphaGlobal[globalID] = alphaList[fi];
            alphaFractionGlobal[globalID] = alphaFractionList[fi];
            
            UGlobal[globalID] = UList[fi];
            UFractionGlobal[globalID] = UFractionList[fi];
            
            FrGlobal[globalID] = Fr;
        }        
    }   
        
    
    // update boundary
    {        
        forAll(VOFRegionData_.globalIDs_, VOFI) //local loop
        {
            label globalID =  VOFData.globalIDs_[VOFI];
            
            VOFData.alphaBound_->refValue()[VOFI] = alphaGlobal[globalID];
            VOFData.alphaBound_->valueFraction()[VOFI] = alphaFractionGlobal[globalID];     
            VOFData.alphaBound_->updateCoeffs();   
                        
            VOFData.UBound_->refValue()[VOFI] = UGlobal[globalID];
            VOFData.UBound_->valueFraction()[VOFI] = UFractionGlobal[globalID];      
            VOFData.UBound_->updateCoeffs();
            
            if(FrGlobal[globalID] >= 1 && VOFData.UBound_->valueFraction()[VOFI] == 1)
            {
                (*VOFData.phiBound_)[VOFI] = UGlobal[globalID] & SfVOFGlobal[globalID];
            }
        }
        
        surfaceScalarField& hBound = *(SWData.hBound_);
        surfaceVectorField& hUBound = *(SWData.hUBound_);
        
        scalarField hTrans = (*SWData.h_).internalField();
        vectorField hUTrans = (*SWData.hU_).internalField();
                
        forAll(SWData.localFaceIDs_, SWi) //local loop
        {            
            label globalID = SWData.globalIDs_[SWi];
            
            label facei = SWData.localFaceIDs_[SWi];
            
            if(facei < SWMesh.nInternalFaces())
            {
                const label own = SWMesh.owner()[facei];
                const label nei = SWMesh.neighbour()[facei];
                
                label ghostCelli = -1;
                 
                if(SWCellType[own] == 1 && SWCellType[nei] == 2)
                {
                    ghostCelli = nei;
                }
                else if(SWCellType[own] == 2 && SWCellType[nei] == 1)
                {
                    ghostCelli = own;
                }
                
                if(ghostCelli >= 0)
                {
                    SWData.h_->internalFieldRef()[ghostCelli] = hGlobal[globalID];
                    SWData.hU_->internalFieldRef()[ghostCelli] = hUGlobal[globalID];   
                }
                
                hBound.internalFieldRef()[facei] = hGlobal[globalID];
                hUBound.internalFieldRef()[facei] = hUGlobal[globalID];
            }
            else
            {
                label patchID = SWMesh.boundaryMesh().whichPatch(facei);
                
                label localFaceID = facei - SWMesh.boundaryMesh()[patchID].start();
                
                const label fCelli = SWMesh.faceOwner()[facei];
                
                hBound.boundaryFieldRef()[patchID][localFaceID] = hGlobal[globalID];
                hUBound.boundaryFieldRef()[patchID][localFaceID] = hUGlobal[globalID];
                                
                hTrans[fCelli] = hGlobal[globalID];
                hUTrans[fCelli] = hUGlobal[globalID];
            }
        }

        labelList nbrType;
        scalarField nbrH;
        vectorField nbrHU;

        syncTools::swapBoundaryCellList
        (
            SWMesh,
            SWCellType,
            nbrType
        );
        
        syncTools::swapBoundaryCellList
        (
            SWMesh,
            hTrans,
            nbrH
        );
        
        syncTools::swapBoundaryCellList
        (
            SWMesh,
            hUTrans,
            nbrHU
        );
        
        for (label faceI = SWMesh.nInternalFaces(); faceI < SWMesh.nFaces(); faceI++)
        {
            const label fCelli = SWMesh.faceOwner()[faceI];
        
            label patchID = SWMesh.boundaryMesh().whichPatch(faceI);
            const polyPatch& patch = SWMesh.boundaryMesh()[patchID];
            
            if (patch.coupled())
            {
                if(SWData.cellType_[fCelli] == 2 && nbrType[faceI-SWMesh.nInternalFaces()] == 1)
                {
                    SWData.h_->internalFieldRef()[fCelli] = nbrH[faceI-SWMesh.nInternalFaces()];
                    SWData.hU_->internalFieldRef()[fCelli] = nbrHU[faceI-SWMesh.nInternalFaces()];
                }
            }
        }
    }
    
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::SWVOFFixedConnector::SWVOFFixedConnector
(
    const vector g,
    const dictionary& regionConnectorEntry,
    const fvMesh& SWRegion,
    const fvMesh& VOFRegion
)
:
    g_(g)
{
    storeBasicData
    (
        regionConnectorEntry,
        SWRegion,
        VOFRegion
    );
    
    createConnection();
    
    updateConnectBoundary();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::SWVOFFixedConnector::update()
{    
    updateConnectBoundary();
}


// ************************************************************************* //
