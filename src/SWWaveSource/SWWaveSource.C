#include "SWWaveSource.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::SWWaveSource::SWWaveSource
(
    const Time& runTime,
    volScalarField& h,
    volVectorField& hU
)
:
runTime_(runTime),
mesh_(h.mesh()),
h_(h),
hU_(hU)
{
    IOdictionary waveSourceDict
    (
        IOobject
        (
            "waveSourceDict",
            runTime.constant(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    );

    label nSource = waveSourceDict.size();

    waveSources_.setSize(nSource);

    label i = 0;
    for (const entry& e : waveSourceDict)
    {
        const dictionary& d = e.dict();
        
        word csvName = d.get<word>("waveSource");

        fileName csvDir = waveSourceDict.globalObjectPath().path() / csvName;

        csvArray csv(csvDir, true);

        List<scalarList> waveData = csv.array();

        if(!waveData.size())
        {
            FatalErrorInFunction
                << "Fail to read waveSource: " << csvDir
                << abort(FatalError);
        }
        
        scalar coeff = d.lookupOrDefault<scalar>("coeff",1);
        
        //Info << h << endl;
        
        word sourceType = d.get<word>("sourceType");
        
        //labelList hList;
        labelHashSet cellSet;
        
        if(sourceType == "point")
        {
            vector position = d.get<vector>("pos");
            
            cellSet.insert(mesh_.findCell(position));
        }
        else if(sourceType == "line")
        {
            vector position1 = d.get<vector>("pos1");
            vector position2 = d.get<vector>("pos2");
            label nSegment = d.get<label>("nSeg");
            
            for (label i = 0; i <= nSegment; ++i)
            {
                scalar alpha = scalar(i) / scalar(nSegment);
                vector pos = (1.0 - alpha)*position1 + alpha*position2;
                label cellI = mesh_.findCell(pos);
                
                if (cellI != -1)
                {
                    cellSet.insert(cellI);
                }
                else
                {
                    Info << "Warning: Point along line is outside the mesh!" << endl;
                }
            }
        }
        
        labelList hList = cellSet.toc();
        scalarList hInitList(hList.size(), 0);
        
        forAll(hInitList, celli)
        {
            if((hList[celli] >= 0))
            {
                hInitList[celli] = h.internalField()[hList[celli]];
            }
        }
        
        waveSources_.set
        (
            i,
            new waveSourceData
            (
                hList,
                waveData,
                coeff,
                hInitList
            )
        );
        
        i++;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::SWWaveSource::update()
{
    if (!waveSources_.size()) return;
    
    scalar tNow = runTime_.value();
        
    forAll(waveSources_,i)
    {        
        waveSourceData& ws = waveSources_[i];

        const List<scalarList>& wv = ws.waveData();
        const labelList& hCells = ws.hList();
        const scalar& coeff = ws.coeff();
        const scalarList& hInitList = ws.hInitList();

        if (wv.size() < 2) continue;

        label k = 0;
        while (k < wv.size() - 2 && tNow >= wv[k+1][0])
        {
            k++;
        }

        scalar t0 = wv[k][0];
        scalar deltah0csv = wv[k][1];

        scalar t1 = wv[k+1][0];
        scalar deltah1csv = wv[k+1][1];

        scalar theta = 0.0;
        if (t1 > t0 + SMALL)
        {
            theta = (tNow - t0) / (t1 - t0);
        }
        theta = max(min(theta, scalar(1)), scalar(0));

        scalar deltahCSV = (1.0 - theta)*deltah0csv + theta*deltah1csv;

        forAll(hCells, ci)
        {
            label cellI = hCells[ci];
            if (cellI >= 0)
            {
                scalar hSet = coeff*deltahCSV + hInitList[ci];
                
                h_[cellI] = hSet;
                
                //hU_[cellI] *= 0;
            }
        }

    }
}

// ************************************************************************* //
