#include "regionConnector.H"
#include "IOdictionary.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * regionConnectors  * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Private Fuctions  * * * * * * * * * * * * * * //

void Foam::regionConnectors::setSWRegionData
(
    SWRegionData& SWData
)
{
    const fvMesh& SWMesh = *SWData.SWMesh_;
    
    const volScalarField& h = SWMesh.lookupObject<volScalarField>("h");
    const volVectorField& hU = SWMesh.lookupObject<volVectorField>("hU");
    
    SWData.cellType_ = new volScalarField
    (
        IOobject
        (
            "cellType",
            SWMesh.time().timeName(),
            SWMesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        SWMesh,
        dimensionedScalar("cellType", dimless, 0.0)
    );
    
    SWData.faceType_ = new surfaceScalarField
    (
        IOobject
        (
            "faceType",
            SWMesh.time().timeName(),
            SWMesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        SWMesh,
        dimensionedScalar("faceType", dimless, 0.0)
    );
    
    SWData.hBound_ = new surfaceScalarField
    (
        IOobject
        (
            "hBound",
            SWMesh.time().timeName(),
            SWMesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        linearInterpolate(h)
    );
    
    SWData.hUBound_ = new surfaceVectorField
    (
        IOobject
        (
            "hUBound",
            SWMesh.time().timeName(),
            SWMesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        linearInterpolate(hU)
    );
    
    
    const fvBoundaryMesh& fvp = SWMesh.boundary();

    forAll(fvp, patchI)
    { 
        if (fvp[patchI].coupled()) 
        {
            SWData.faceType_->boundaryFieldRef().set
            (
                patchI,
                fvsPatchField<scalar>::New
                (
                    "calculated",               
                    SWMesh.boundary()[patchI],    
                    *(SWData.faceType_)             
                )
            );
            
            SWData.hBound_->boundaryFieldRef().set
            (
                patchI,
                fvsPatchField<scalar>::New
                (
                    "calculated",               
                    SWMesh.boundary()[patchI],    
                    *(SWData.hBound_)             
                )
            );
            
            SWData.hUBound_->boundaryFieldRef().set
            (
                patchI,
                fvsPatchField<vector>::New
                (
                    "calculated",               
                    SWMesh.boundary()[patchI],    
                    *(SWData.hUBound_)             
                )
            );
        }
    }
}


void Foam::regionConnectors::updateSWRegionData
(
    SWRegionData& SWData
)
{
    const fvMesh& SWMesh = *SWData.SWMesh_;
    
    volScalarField& cellType = *SWData.cellType_;
    surfaceScalarField& faceType = *SWData.faceType_;
    
    cellType.internalFieldRef() = 0.0;
    faceType.internalFieldRef() = 0.0;
    forAll(faceType.boundaryField(), patchI)
    {
        faceType.boundaryFieldRef()[patchI] = 0.0;
    }    
    
    if(SWVOFFixedConnectorList_.size())
    {
        forAll(SWVOFFixedConnectorList_,i)
        {
            const labelList& SWCellType = SWVOFFixedConnectorList_[i].SWCellType();
            const labelList& SWFaceType = SWVOFFixedConnectorList_[i].SWFaceType();
            
            forAll(SWCellType, celli)
            {
                if(SWCellType[celli] > 0)
                {
                    cellType.internalFieldRef()[celli] = SWCellType[celli];
                }
            }
            
            forAll(SWFaceType, facei)
            {
                if(SWFaceType[facei] > 0)
                {
                    if(facei < SWMesh.nInternalFaces())
                    {
                        faceType.internalFieldRef()[facei] = SWFaceType[facei];
                    }
                    else
                    {
                        label patchID = SWMesh.boundaryMesh().whichPatch(facei);                    
                        label localFaceID = facei - SWMesh.boundaryMesh()[patchID].start();

                        faceType.boundaryFieldRef()[patchID][localFaceID] = SWFaceType[facei];                                
                    }
                }
            }
        }
    }
    
    if(SWVOFArbitraryConnectorList_.size())
    {
        forAll(SWVOFArbitraryConnectorList_,i)
        {
            const labelList& SWCellType = SWVOFArbitraryConnectorList_[i].SWCellType();
            const labelList& SWFaceType = SWVOFArbitraryConnectorList_[i].SWFaceType();
            
            forAll(SWCellType, celli)
            {
                if(SWCellType[celli] > 0)
                {
                    cellType.internalFieldRef()[celli] = SWCellType[celli];
                }
            }
            
            forAll(SWFaceType, facei)
            {
                if(SWFaceType[facei] > 0)
                {
                    if(facei < SWMesh.nInternalFaces())
                    {
                        faceType.internalFieldRef()[facei] = SWFaceType[facei];
                    }
                    else
                    {
                        label patchID = SWMesh.boundaryMesh().whichPatch(facei);                    
                        label localFaceID = facei - SWMesh.boundaryMesh()[patchID].start();

                        faceType.boundaryFieldRef()[patchID][localFaceID] = SWFaceType[facei];                                
                    }
                }
            }
        }
    }
    
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionConnectors::regionConnectors
(
    const Time& runTime,
    const vector g,
    const PtrList<fvMesh>& SWRegions,
    const PtrList<fvMesh>& VOFRegions,
    const PtrList<fvMesh>& laplacianRegions,
    const wordList& SWRegionNames,
    const wordList& VOFRegionNames,
    const wordList& laplacianRegionNames
)
:
    g_(g),
    SWVOFFixedConnectorList_(0),
    SWVOFArbitraryConnectorList_(0),
    VOFLaplacianConnectorList_(0)
{
    regionConnectProperties_ = IOdictionary
    (
        IOobject
        (
            "regionConnectProperties",
            runTime.constant(),
            runTime,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    );
    
    // SWVOFConnectors
    if (regionConnectProperties_.found("SWVOFConnectors"))
    {
        const dictionary& SWVOFConnectDict =
            regionConnectProperties_.subDict("SWVOFConnectors");

        for (const entry& e : SWVOFConnectDict)
        {
            const dictionary& d = e.dict();
            
            if (!d.found("SWRegion") || !d.found("VOFRegion"))
            {
                FatalErrorInFunction << "choose Regions for SWVOFConnectors" << abort(FatalError);
            }
            
            const word SWRegionName = d.get<word>("SWRegion");
            const word VOFRegionName = d.get<word>("VOFRegion");
            
            label SWI = -1;
            label VOFI = -1;
            
            forAll(SWRegionNames, j)
            {  
                if(SWRegionName == SWRegionNames[j])
                {
                    SWI = j;
                }
            }
            
            forAll(VOFRegionNames, j)
            {  
                if(VOFRegionName == VOFRegionNames[j])
                {
                    VOFI = j;
                }
            }
            
            if(SWI == -1 || VOFI == -1)
            {   
                FatalErrorInFunction << "choose correct region" << abort(FatalError);
            }
            else
            {
                Info << "connect regions, SWRegion: " << SWRegionName << ", VOFRegion: " << VOFRegionName << endl;
            }
            
            label SWDataIndex = -1;
                    
            forAll(SWRegionData_, i)
            {
                if (SWRegionData_[i].SWMesh_ == &SWRegions[SWI])
                {
                    SWDataIndex = i;
                    break;
                }
            }
            
            if(SWDataIndex == -1)
            {
                SWRegionData_.append(new SWRegionData);
                SWDataIndex = SWRegionData_.size() - 1;
            
                SWRegionData& SWData = SWRegionData_.last();
                SWData.SWMesh_ = &SWRegions[SWI];
                
                setSWRegionData(SWData);
            }
            
            
            const word mode = d.get<word>("mode");
            if (mode == "fixed")
            {
                SWVOFFixedConnectorList_.append
                (
                    new SWVOFFixedConnector
                    (
                        g_,
                        d,
                        SWRegions[SWI],
                        VOFRegions[VOFI]
                    )
                );
                
                SWRegionData_[SWDataIndex].SWVOFFixedConnectorList_.append
                (
                    &SWVOFFixedConnectorList_.last()
                );
            }
            else if(mode == "arbitrary")
            {
                SWVOFArbitraryConnectorList_.append
                (
                    new SWVOFArbitraryConnector
                    (
                        g_,
                        d,
                        SWRegions[SWI],
                        VOFRegions[VOFI]
                    )
                );
                
                SWRegionData_[SWDataIndex].SWVOFArbitraryConnectorList_.append
                (
                    &SWVOFArbitraryConnectorList_.last()
                );
            }
            else
            {
                FatalErrorInFunction << "choose correct connect mode" << abort(FatalError);
            }
        }
        
        // Update SW label
        forAll(SWRegionData_, i)
        {
            updateSWRegionData(SWRegionData_[i]);
        }
    }
    else
    {
        Info << "No SWVOFConnectors defined." << endl;
    }
    
    //VOFLaplacianConnectors
    if (regionConnectProperties_.found("VOFLaplacianConnectors"))
    {    
        const dictionary& VOFLaplacianConnectDict =
            regionConnectProperties_.subDict("VOFLaplacianConnectors");
        
        for (const entry& e : VOFLaplacianConnectDict)
        {
            const dictionary& d = e.dict();
            
            if (!d.found("VOFRegion") || !d.found("laplacianRegion"))
            {
                FatalErrorInFunction << "choose Regions for VOFLaplacianConnectors" << abort(FatalError);
            }
                    
            const word VOFRegionName = d.get<word>("VOFRegion");
            const word laplacianRegionName = d.get<word>("laplacianRegion");
            
            label VOFI = -1;
            label laplacianI = -1;
            
            forAll(VOFRegionNames, j)
            {  
                if(VOFRegionName == VOFRegionNames[j])
                {
                    VOFI = j;
                }
            }
            
            forAll(laplacianRegionNames, j)
            {  
                if(laplacianRegionName == laplacianRegionNames[j])
                {
                    laplacianI = j;
                }
            }
            
            if(VOFI == -1 || laplacianI == -1)
            {   
                FatalErrorInFunction << "choose correct region" << abort(FatalError);
            }
            else
            {
                Info << "connect regions, VOFRegion: " << VOFRegionName << ", laplacianRegion: " << laplacianRegionName << endl;
            }

            VOFLaplacianConnectorList_.append
            (
                new VOFLaplacianConnector
                (
                    g_,
                    d,
                    VOFRegions[VOFI],
                    laplacianRegions[laplacianI]
                )
            );
        }
    }
    else
    {
        Info << "No VOFLaplacianConnectors defined." << endl;
    }
}




// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::regionConnectors::updateSWVOF()
{
    if(SWVOFFixedConnectorList_.size())
    {
        forAll(SWVOFFixedConnectorList_,i)
        {
            SWVOFFixedConnectorList_[i].update();
        }
    }
    
    if(SWVOFArbitraryConnectorList_.size())
    {
        forAll(SWVOFArbitraryConnectorList_,i)
        {
            SWVOFArbitraryConnectorList_[i].update();
        }
    }
    
    //uodate SW label
    forAll(SWRegionData_, i)
    {
        SWRegionData& SWData = SWRegionData_[i];
        
        bool updateConnection = false;
        
        forAll(SWVOFArbitraryConnectorList_,i)
        {
            if(SWVOFArbitraryConnectorList_[i].updateConnection())
            {
                updateConnection = true;
            }
        }
        
        if(updateConnection)
        {
            updateSWRegionData(SWData);
        }
    }
}

void Foam::regionConnectors::updateVOFLaplacian()
{    
    if(VOFLaplacianConnectorList_.size())
    {
        forAll(VOFLaplacianConnectorList_,i)
        {
            VOFLaplacianConnectorList_[i].update();
        }
    }
}


