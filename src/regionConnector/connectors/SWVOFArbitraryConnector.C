#include "SWVOFArbitraryConnector.H"
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
templateList Foam::SWVOFArbitraryConnector::generateGlobalList
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


Foam::labelList Foam::SWVOFArbitraryConnector::generateGlobalIDs
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
List<FieldType> Foam::SWVOFArbitraryConnector::extractFieldSubset
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

void Foam::SWVOFArbitraryConnector::storeBasicData  // store data
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
        
        VOFData.donorCell_ = labelList(VOFMesh->nCells(), 0);
            
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
        
        VOFData.SfBound_ = VOFMesh->boundaryMesh()[patchID].faceAreas();
        VOFData.CBound_ = VOFMesh->boundaryMesh()[patchID].faceCentres();
        
        VOFData.globalIDs_ = generateGlobalIDs(VOFRegionData_.UBound_->size());
        
        VOFData.SWDonorProcCell_ = List<labelList>(generateGlobalList(VOFData.SfBound_).size(), labelList{-1, -1});
    }
}


// * * * * * * * * * * * * * * * * Update Functions  * * * * * * * * * * * * * * //


void Foam::SWVOFArbitraryConnector::createConnection()
{
    const fvMesh& SWMesh = *SWRegionData_.SWMesh_;
    const fvMesh& VOFMesh = *VOFRegionData_.VOFMesh_;
    
    SWRegionData& SWData = SWRegionData_;
    VOFRegionData& VOFData = VOFRegionData_; 
        
    labelList& SWCellType = SWData.cellType_;        
    labelList& SWFaceType = SWData.faceType_;
    SWCellType = labelList(SWMesh.nCells(), -1);
    SWFaceType = labelList(SWMesh.nFaces(), 0);    
    
    SWData.localFaceIDs_.clear();
    SWData.globalIDs_.clear();
    
    labelList& VOFDonorCell = VOFData.donorCell_;
    VOFDonorCell = labelList(VOFMesh.nCells(), -1);
    
    vectorField VOFCBoundGlobal = generateGlobalList(VOFData.CBound_);
    vectorField VOFSfBoundGlobal = generateGlobalList(VOFData.SfBound_);
    
    VOFData.SWDonorProcCell_ = List<labelList>(VOFCBoundGlobal.size());    
    
    {
        boundBox SWBB(SWMesh.points(), true);
        scalar SWmidZ = (SWBB.min().z() + SWBB.max().z())/2;
        
        forAll(VOFCBoundGlobal, fbi)
        {
            vector VOFC = VOFCBoundGlobal[fbi];
            VOFC.z() = SWmidZ;
            
            //label SWI = SWMesh.findCell(VOFC);
            label SWI = SWMesh.findNearestCell(VOFC);    
            {
                scalar dist = mag(VOFC - SWMesh.C()[SWI]);
                scalar minDist = dist;
                reduce(minDist, minOp<scalar>());

                if(dist != minDist)
                {
                    SWI = -1;
                }
            }  

            label procI = -1;
            label cellI = -1;
            
            if(SWI >= 0)
            {
                SWCellType[SWI] = 1; //donor for VOF
                
                //store SW donor for VOF
                procI = Pstream::myProcNo();
                cellI = SWI;
            }
            
            reduce(procI, maxOp<label>());
            reduce(cellI, maxOp<label>());
            
            if(procI >= 0 && cellI >= 0)
            {
                labelList SWDPC(2);
                SWDPC[0] = procI;
                SWDPC[1] = cellI;

                VOFData.SWDonorProcCell_[fbi] = SWDPC;
            }
            else
            {
                FatalErrorInFunction << "Donor for VOF was not created, No. " << fbi << exit(FatalError);
            }
        }
        
        // label for outside cell
        {
            boundBox faceBB(VOFCBoundGlobal, true);
            
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
        }        
        
        forAll(SWFaceType, SWfi)
        {
            if(SWfi < SWMesh.nInternalFaces())
            {
                const label own = SWMesh.owner()[SWfi];
                const label nei = SWMesh.neighbour()[SWfi];
                 
                if(SWCellType[own] == 1 && SWCellType[nei] == -1)
                {
                    SWFaceType[SWfi] = 1;
                    SWCellType[nei] = 2;
                    
                    SWData.localFaceIDs_.append(SWfi);
                }
                else if(SWCellType[own] == -1 && SWCellType[nei] == 1)
                {
                    SWFaceType[SWfi] = 1;
                    SWCellType[own] = 2;
                    
                    SWData.localFaceIDs_.append(SWfi);
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
                if(ownType == 1 && neiType == -1)
                {
                    SWFaceType[faceI] = 1;
                    
                    SWData.localFaceIDs_.append(faceI);
                }
                else if(ownType == -1 && neiType == 1)
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
        
        SWData.globalIDs_ = generateGlobalIDs(SWData.localFaceIDs_.size());
        labelList allIDs = generateGlobalList(SWData.globalIDs_);
        SWData.VOFDonorProcCells_ = List<List<labelList>>(allIDs.size());
    }
    
    
    {
        const pointField& points = SWMesh.points();
        List<pointList> localSWFacePoints(SWData.localFaceIDs_.size());
        
        forAll(localSWFacePoints, i)
        {
            const label SWfi = SWData.localFaceIDs_[i];
            const face& f = SWMesh.faces()[SWfi]; 
            
            forAll(f, pi)
            {
                localSWFacePoints[i].append(points[f[pi]]);
            }
        }
        
        List<pointList> globalSWFacePoints(generateGlobalList(localSWFacePoints));
                
        forAll(SWData.VOFDonorProcCells_, SWI)
        {
            const List<point>& facePoints = globalSWFacePoints[SWI];

            boundBox faceBB(facePoints, true);
            
            vector lengths = faceBB.span();  
            
            scalar mergin = sqrt(sqr(lengths.x())+sqr(lengths.y()));
            
            vector A = vector(faceBB.min().x(), faceBB.min().y(), 0);
            vector B = vector(faceBB.max().x(), faceBB.max().y(), 0);

            point newMin
            (
                faceBB.min().x() - mergin,
                faceBB.min().y() - mergin,
                faceBB.min().z() - mergin
            );

            point newMax
            (
                faceBB.max().x() + mergin,
                faceBB.max().y() + mergin,
                faceBB.max().z() + mergin
            );
            
            boundBox expandedXYBox(newMin, newMax);
            
            List<label> candidateCells;
            
            forAll(VOFMesh.C(), celli) 
            { 
                if(expandedXYBox.contains(VOFMesh.C()[celli]))
                { 
                    candidateCells.append(celli); 
                } 
            }
            
            List<labelList> SWDonorProcCells(0);
            
            forAll(candidateCells, i)
            {
                label fci = candidateCells[i];
                
                bool intersect = false;
                
                {
                    const labelList& cFaces = VOFMesh.cells()[fci];
                    labelHashSet edgeSet;

                    for (const label faceI : cFaces)
                    {
                        const labelList& fEdges = VOFMesh.faceEdges()[faceI];
                        for (const label edgeI : fEdges)
                        {
                            edgeSet.insert(edgeI);
                        }
                    }

                    for (const label edgeI : edgeSet)
                    {
                        const edge& e = VOFMesh.edges()[edgeI];
                        
                        point p1 = VOFMesh.points()[e.start()];
                        point p2 = VOFMesh.points()[e.end()];
                        
                        vector C = vector(p1.x(), p1.y(), 0);
                        vector D = vector(p2.x(), p2.y(), 0);
                        
                        if 
                        (
                            max(A.x(), B.x()) < min(C.x(), D.x()) - SMALL ||
                            min(A.x(), B.x()) > max(C.x(), D.x()) + SMALL ||
                            max(A.y(), B.y()) < min(C.y(), D.y()) - SMALL ||
                            min(A.y(), B.y()) > max(C.y(), D.y()) + SMALL
                        )
                        {
                            continue;
                        }
                        
                        vector AB = B - A;
                        vector AC = C - A;
                        vector AD = D - A;
                        
                        vector CD = D - C;
                        vector CA = A - C;
                        vector CB = B - C;
                        
                        scalar cp1 = AB.x()*AC.y() - AB.y()*AC.x();
                        scalar cp2 = AB.x()*AD.y() - AB.y()*AD.x();
                        scalar cp3 = CD.x()*CA.y() - CD.y()*CA.x();
                        scalar cp4 = CD.x()*CB.y() - CD.y()*CB.x();
                        
                        scalar tol = 1e-6;

                        if (mag(cp1) < tol && mag(cp2) < tol && mag(cp3) < tol && mag(cp4) < tol)
                        {
                            intersect = true; 
                            break;
                        }
                        else if ((cp1 * cp2 <= 0) && (cp3 * cp4 <= 0))
                        {
                            intersect = true;
                            break;
                        }
                    }
                }
                
                if(intersect)
                {
                    const label procI = Pstream::myProcNo();
                    SWDonorProcCells.append(labelList{procI, fci});
                    VOFDonorCell[fci] = 1;
                }
            }
            
            SWData.VOFDonorProcCells_[SWI] = SWDonorProcCells;
        }
        
        forAll(globalSWFacePoints, SWI)
        {
            List<labelList> SWDonorProcCells = SWData.VOFDonorProcCells_[SWI];
            SWData.VOFDonorProcCells_[SWI] = generateGlobalList(SWDonorProcCells);
            
            if(SWData.VOFDonorProcCells_[SWI].size() == 0)
            {                
                Info << "SW points: " << globalSWFacePoints[SWI] << endl;
                
                FatalErrorInFunction << "Donor for SW was not created, No. " << SWI << exit(FatalError);
            }
        }
    }
    
}


void Foam::SWVOFArbitraryConnector::updateConnectBoundary()
{
    SWRegionData& SWData = SWRegionData_;
    VOFRegionData& VOFData = VOFRegionData_;
    
    const fvMesh& SWMesh = *SWRegionData_.SWMesh_;
    const fvMesh& VOFMesh = *VOFRegionData_.VOFMesh_;
    
    // update flude boundary
    {
        scalarField alphaLocal = VOFData.alphaBound_->patchInternalField();
        scalarField alphaGlobal = generateGlobalList(alphaLocal);
        scalarField alphaFractionGlobal = scalarField(alphaGlobal.size(), -1);
       
        vectorField ULocal = VOFData.UBound_->patchInternalField();
        vectorField UGlobal = generateGlobalList(ULocal);
        scalarField UFractionGlobal = scalarField(UGlobal.size(), -1);
        
        //generate donor
        scalarField hDonor(VOFData.SWDonorProcCell_.size(), 0);
        vectorField hUDonor(VOFData.SWDonorProcCell_.size(), vector::zero);
        
        forAll(VOFData.SWDonorProcCell_, fi)
        {
            label proci = VOFData.SWDonorProcCell_[fi][0];
            label celli = VOFData.SWDonorProcCell_[fi][1];
            
            if(proci == Pstream::myProcNo())
            {
                hDonor[fi] = (*SWData.h_)[celli];
                hUDonor[fi] = (*SWData.hU_)[celli];
            }
            reduce(hDonor[fi], sumOp<scalar>());
            reduce(hUDonor[fi], sumOp<vector>());
        }
        
        forAll(UGlobal, fi)
        {
            alphaFractionGlobal[fi] = 0;
            UFractionGlobal[fi] = 1;
            
            if(alphaGlobal[fi] < 0.01)
            {
                UFractionGlobal[fi] = 0;
                UGlobal[fi] *= 0;
                
                alphaFractionGlobal[fi] = 1;
                alphaGlobal[fi] = 0;
            }
            else
            {
               vector UDonor = hUDonor[fi]/stabilise(hDonor[fi], SMALL);            
               UGlobal[fi] = alphaGlobal[fi]*UDonor;
            }
        }   
        
        forAll(VOFRegionData_.globalIDs_, VOFI)
        {
            label globalID =  VOFData.globalIDs_[VOFI];
            
            VOFData.alphaBound_->refValue()[VOFI] = alphaGlobal[globalID];
            VOFData.alphaBound_->valueFraction()[VOFI] = alphaFractionGlobal[globalID];     
            VOFData.alphaBound_->updateCoeffs();   
                        
            VOFData.UBound_->refValue()[VOFI] = UGlobal[globalID];
            VOFData.UBound_->valueFraction()[VOFI] = UFractionGlobal[globalID];      
            VOFData.UBound_->updateCoeffs();
        }
    }
    // update SW boundary
    {
        surfaceScalarField& hBound = *(SWData.hBound_);
        surfaceVectorField& hUBound = *(SWData.hUBound_);

        const labelList& cellType = SWData.cellType_;
        const labelList& faceType = SWData.faceType_;
        
        scalarField hLocal = scalarField(SWData.localFaceIDs_.size(), 0.0);
        vectorField hULocal = vectorField(SWData.localFaceIDs_.size(), vector::zero);
        // set initial as zeroGradient
        forAll(SWData.localFaceIDs_, SWi)//local loop
        {            
            label facei = SWData.localFaceIDs_[SWi];
            
            if(facei < SWMesh.nInternalFaces())
            {
                const label own = SWMesh.owner()[facei];
                const label nei = SWMesh.neighbour()[facei];
                
                label boundCelli = -1;
                 
                if(cellType[own] == 1)
                {
                    boundCelli = own;
                }
                else if(cellType[nei] == 1)
                {
                    boundCelli = nei;
                }
                
                hLocal[SWi] = (*SWData.h_)[boundCelli];
                hULocal[SWi] = (*SWData.hU_)[boundCelli];
            }
            else
            {
                label patchID = SWMesh.boundaryMesh().whichPatch(facei);
                const polyPatch& patch = SWMesh.boundaryMesh()[patchID];
                
                label localFaceID = facei - SWMesh.boundaryMesh()[patchID].start();

                if (patch.coupled())
                {
                    const labelUList& faceCells = patch.faceCells();
                    const label fCelli = faceCells[localFaceID]; 
                    
                    if(cellType[fCelli] == 1)
                    {
                        hLocal[SWi] = (*SWData.h_)[fCelli];
                        hULocal[SWi] = (*SWData.hU_)[fCelli];
                    }
                    else if(cellType[fCelli] == 2)
                    {
                        hLocal[SWi] = (*SWData.h_).boundaryField()[patchID][localFaceID];
                        hULocal[SWi] = (*SWData.hU_).boundaryField()[patchID][localFaceID];   
                    }
                }
            }
        }
        
        scalarField hGlobal = generateGlobalList(hLocal);
        vectorField hUGlobal = generateGlobalList(hULocal);
        
        
        forAll(SWData.VOFDonorProcCells_, SWi)//global loop
        {
            List<labelList> VOFDonorProcCells = SWData.VOFDonorProcCells_[SWi];
            
            scalarField alphaDonor(VOFDonorProcCells.size(), 0);
            vectorField UDonor(VOFDonorProcCells.size(), vector::zero);
            
            scalarField meshVDonor(VOFDonorProcCells.size(), 0);
            
            scalar minZ = GREAT;
            scalar maxZ = -GREAT;
            
            forAll(VOFDonorProcCells, i)
            {
                const label proci = VOFDonorProcCells[i][0];
                const label fci = VOFDonorProcCells[i][1];
                
                if(proci == Pstream::myProcNo())
                {
                    alphaDonor[i] = (*VOFData.alpha_)[fci];
                    UDonor[i] = (*VOFData.U_)[fci];
                    meshVDonor[i] = VOFMesh.V()[fci];
                    
                    const cell& c = VOFMesh.cells()[fci];
                    const pointField& pts = VOFMesh.points();

                    forAll(c, fi)
                    {
                        const face& f = VOFMesh.faces()[c[fi]];
                        forAll(f, pi)
                        {
                            scalar z = pts[f[pi]].z();
                            minZ = min(minZ, z);
                            maxZ = max(maxZ, z);
                        }
                    }
                }
                reduce(alphaDonor[i], sumOp<scalar>());
                reduce(UDonor[i], sumOp<vector>());
                reduce(meshVDonor[i], sumOp<scalar>());
                
                reduce(minZ, minOp<scalar>());
                reduce(maxZ, maxOp<scalar>());
            }
            scalar donorHeight = maxZ - minZ;
            
            scalar hFrom3D = 0.0;
            scalar totalV = 0.0;
            
            forAll(VOFDonorProcCells, i)
            {
                scalar alpha = alphaDonor[i];
                scalar V = meshVDonor[i];

                hFrom3D += alpha * V;
                totalV += V;
            }
            scalar bottomArea = totalV / donorHeight;
            
            hFrom3D /= bottomArea;
                        
            hGlobal[SWi] = hFrom3D;
        }
        
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
                 
                if(cellType[own] == 1 && cellType[nei] == 2)
                {
                    ghostCelli = nei;
                }
                else if(cellType[own] == 2 && cellType[nei] == 1)
                {
                    ghostCelli = own;
                }
                
                if(ghostCelli >= 0)
                {
                    SWData.h_->internalFieldRef()[ghostCelli] = hGlobal[globalID];
                    SWData.hU_->internalFieldRef()[ghostCelli] = hUGlobal[globalID];   
                }
                
                SWData.hBound_->internalFieldRef()[facei] = hGlobal[globalID];
                SWData.hUBound_->internalFieldRef()[facei] = hUGlobal[globalID];
            }
            else
            {
                label patchID = SWMesh.boundaryMesh().whichPatch(facei);
                const polyPatch& patch = SWMesh.boundaryMesh()[patchID];
                
                label localFaceID = facei - SWMesh.boundaryMesh()[patchID].start();

                if (patch.coupled())
                {
                    const labelUList& faceCells = patch.faceCells();
                    const label fCelli = faceCells[localFaceID]; 
                    
                    if(cellType[fCelli] == 1)
                    {
                        hBound.boundaryFieldRef()[patchID][localFaceID] = hGlobal[globalID];
                        hUBound.boundaryFieldRef()[patchID][localFaceID] = hUGlobal[globalID];
                        
                        hTrans[fCelli] = hGlobal[globalID];
                        hUTrans[fCelli] = hUGlobal[globalID];
                    }
                }
            }
        }

        labelList nbrType;
        scalarField nbrH;
        vectorField nbrHU;

        syncTools::swapBoundaryCellList
        (
            SWMesh,
            cellType,
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
                if(cellType[fCelli] == 2 && nbrType[faceI-SWMesh.nInternalFaces()] == 1)
                {
                    SWData.h_->internalFieldRef()[fCelli] = nbrH[faceI-SWMesh.nInternalFaces()];
                    SWData.hU_->internalFieldRef()[fCelli] = nbrHU[faceI-SWMesh.nInternalFaces()];
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::SWVOFArbitraryConnector::SWVOFArbitraryConnector
(
    const vector g,
    const dictionary& regionConnectorEntry,
    const fvMesh& SWRegion,
    const fvMesh& VOFRegion
)
:
    g_(g),
    updateConnection_(false)
{
    storeBasicData
    (
        regionConnectorEntry,
        SWRegion,
        VOFRegion
    );
    
    VOFRegionData& VOFData = VOFRegionData_; 
    
    vectorField VOFCBoundGlobal = generateGlobalList(VOFData.CBound_);
    VOFC_ = average(VOFCBoundGlobal);
    
    
    createConnection();
    
    updateConnectBoundary();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::SWVOFArbitraryConnector::update()
{
    const fvMesh* VOFMesh = VOFRegionData_.VOFMesh_;
    
    VOFRegionData& VOFData = VOFRegionData_; 
    
    VOFData.SfBound_ = VOFMesh->boundaryMesh()[VOFData.patchID_].faceAreas();
    VOFData.CBound_ = VOFMesh->boundaryMesh()[VOFData.patchID_].faceCentres();
    
    vectorField VOFCBoundGlobal = generateGlobalList(VOFData.CBound_);
    vector VOFCNew = average(VOFCBoundGlobal);
    
    Info << "VOFC-old: " << VOFC_ << ", VOFCNew: " << VOFCNew << endl;
    
    updateConnection_ = false;
    
    if (mag(VOFCNew - VOFC_) > SMALL)
    {
        VOFC_ = VOFCNew;
        
        Info << "update connection due to mesh motion" << endl;
        
        createConnection();
        
        updateConnection_ = true;
    }
    
    updateConnectBoundary();
}


// ************************************************************************* //
