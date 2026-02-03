#include "VOFLaplacianConnector.H"
#include "IOdictionary.H"
#include "Time.H"

// * * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::VOFLaplacianConnector::findNearestFace
(
	const vector position,
	const vectorList& positionList
)
{
	label patchFaceID = -1; //set initial ID
	scalar distance = GREAT;
		
	forAll(positionList,facei)
	{
		scalar distance2 = mag(positionList[facei] - position);
		
		if(distance2 < distance)
		{
			patchFaceID = facei;
			distance = distance2;
		}
	}
	
	if (patchFaceID == -1)
    {
        FatalErrorInFunction << "No suitable face found with acceptable angle" << exit(FatalError);
    }
	
	return patchFaceID;
}

template<class templateList>
templateList Foam::VOFLaplacianConnector::generateGlobalList
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
    
    shareList[proci] = list;
         
    Pstream::gatherList(shareList);
    Pstream::scatterList(shareList);
    
    forAll(shareList, listI)
    {
        globalList.append(shareList[listI]);
    }
    
    return globalList;
}

Foam::labelList Foam::VOFLaplacianConnector::generateGlobalIDs
(
    label localSize
)
{
    label nProcs = Pstream::nProcs();
    label myRank = Pstream::myProcNo();

    List<labelList> allLocalSizes(nProcs);
    forAll(allLocalSizes, i)
    {
        allLocalSizes[i] = labelList(0);
    }
    allLocalSizes[myRank] = labelList(1, localSize);

    Pstream::gatherList(allLocalSizes);
    Pstream::scatterList(allLocalSizes);

    List<label> offset(nProcs, 0);
    for (label i = 1; i < nProcs; ++i)
    {
        offset[i] = offset[i-1] + allLocalSizes[i-1][0];
    }

    label start = offset[myRank];

    labelList globalIDs(localSize);
    for (label i = 0; i < localSize; ++i)
    {
        globalIDs[i] = start + i;
    }

    return globalIDs;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::VOFLaplacianConnector::VOFLaplacianConnector
(
    const vector g,
    const dictionary& connectorDict,
    const fvMesh& VOFRegion,
    const fvMesh& laplacianRegion
)
:
    g_(g)
{
    VOFData_.VOFMesh_ = &VOFRegion;
    laplacianData_.laplacianMesh_ = &laplacianRegion;
    
    const fvMesh* VOFMesh = VOFData_.VOFMesh_;
    const fvMesh* laplacianMesh = laplacianData_.laplacianMesh_;
    
    VOFRegionPatchData& VOFData = VOFData_;
    laplacianRegionData& laplacianData = laplacianData_;    
    
    VOFPatchName_ = connectorDict.get<word>("VOFPatch");
    laplacianPatchName_ = connectorDict.get<word>("laplacianPatch");
    
    
    label VOFPatchID = VOFMesh->boundaryMesh().findPatchID(VOFPatchName_);
    if (VOFPatchID < 0)
    {
        FatalErrorInFunction
            << "Patch " << VOFPatchName_ << " not found in VOFRegion " << exit(FatalError);
    }
    else
    {
        VOFData.patchID_ = VOFPatchID;
    }
    
    label laplacianPatchID = laplacianMesh->boundaryMesh().findPatchID(laplacianPatchName_);
    if (laplacianPatchID < 0)
    {
        FatalErrorInFunction
            << "Patch " << laplacianPatchName_ << " not found in laplacianRegion " << exit(FatalError);
    }
    else
    {
        laplacianData.patchID_ = laplacianPatchID;
    }
    
    
    //VOFRegion
    {
        VOFData.p_ = VOFMesh->getObjectPtr<volScalarField>("p");
        
        VOFData.globalIDs_ = generateGlobalIDs(VOFData.p_->size());
    }
    
    //laplacianRegion
    {
        volScalarField* pptr = laplacianMesh->getObjectPtr<volScalarField>("p");
        laplacianData.p_ = &refCast<fixedValueFvPatchField<scalar>>(pptr->boundaryFieldRef()[laplacianPatchID]);
        
        laplacianData.globalIDs_ = generateGlobalIDs(laplacianData.p_->size());
    }
    
    const polyPatch& VOFPatch = VOFMesh->boundaryMesh()[VOFPatchID];
    const polyPatch& laplacianPatch = laplacianMesh->boundaryMesh()[laplacianPatchID];
    
    vectorList localVOFPositions = VOFPatch.faceCentres();
    vectorList globalVOFPositions = generateGlobalList(localVOFPositions);   
   	
   	vectorList localLaplacianPositions = laplacianPatch.faceCentres();
	vectorList globalLaplacianPositions = generateGlobalList(localLaplacianPositions);
    
    VOFIDsForLaplacian_.setSize(globalLaplacianPositions.size()); 
    forAll(globalLaplacianPositions,facei)
	{
		label fliudFacei = findNearestFace(globalLaplacianPositions[facei],  globalVOFPositions);
		
		VOFIDsForLaplacian_[facei] = fliudFacei;
	}

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::VOFLaplacianConnector::update()
{
    VOFRegionPatchData& VOFData = VOFData_;
    laplacianRegionData& laplacianData = laplacianData_;
    
    scalarField& localPVOF = VOFData.p_->boundaryFieldRef()[VOFData.patchID_];
    
    scalarField globalPVOF = generateGlobalList(localPVOF);
    
    forAll(laplacianData_.globalIDs_, laplacianI) // loop in local
    {
        label globalID =  laplacianData_.globalIDs_[laplacianI];
        
        label fID = VOFIDsForLaplacian_[globalID];
        
        (*laplacianData.p_)[laplacianI] = globalPVOF[fID];
    }
}


// ************************************************************************* //
