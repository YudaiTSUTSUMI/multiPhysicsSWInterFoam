#include "fvCFD.H"
#include "triSurface.H"
#include "triSurfaceSearch.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createNamedMesh.H"
    #include "createFields.H"
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
    {
        h = dimensionedScalar("", dimensionSet(0,1,0,0,0,0,0), 0.0);

        const word mode = setSWDict.get<word>("mode");
        
        if(mode == "dict")
        {
            scalar initialh0 = setSWDict.get<scalar>("initialh0");
            h0 = dimensionedScalar("", dimensionSet(0,1,0,0,0,0,0), initialh0);
            
            const volVectorField& C = mesh.C();
            
            const dictionary& setSlopeDict = setSWDict.subDict("slope");        

            for (const entry& e : setSlopeDict)
            {
                const dictionary& dict = e.dict();
                
                vector startPos = dict.get<vector>("startPosition");
                vector direction = dict.get<vector>("direction");
                direction /= mag(direction) + SMALL;
                
                scalar slope = dict.get<scalar>("slope");
                
                scalar baseh0 = dict.get<scalar>("baseh0");
                
                forAll(h0, celli)
                {
                    forAll(h0, celli)
                    {
                        const vector& cellPos = C[celli];

                        scalar dist = (cellPos - startPos) & direction;

                        if (dist > 0)
                        {
                            h0[celli] = baseh0 + slope * dist;
                        }
                    }
                }
            }
        }
        else if(mode == "triSurface")
        {
            const word triSurfaceName(setSWDict.get<word>("triSurface"));
            
            const fileName triPath = "constant/triSurface" / triSurfaceName;
            
            triSurface triMesh = triSurface(triPath);
            triSurfaceSearch querySurf(triMesh);

            const volVectorField& C = mesh.C();

            forAll(h0, celli)
            {
                const point& cellPos = C[celli];
                const scalar cl = setSWDict.getOrDefault<scalar>("checkLength", 100.0);

                point startDown = cellPos;
                point endDown   = cellPos - vector(0, 0, cl);

                point startUp   = cellPos;
                point endUp     = cellPos + vector(0, 0, cl);

                pointField startPts(1);
                pointField endPts(1);

                List<pointIndexHit> hitResults;

                scalar terrainZ = -GREAT;

                startPts[0] = startDown;
                endPts[0] = endDown;
                querySurf.findLine(startPts, endPts, hitResults);

                bool hitD = !hitResults.empty() && hitResults[0].hit();
                point hitPtD;
                if (hitD)
                {
                    hitPtD = hitResults[0].hitPoint();
                } 

                startPts[0] = startUp;
                endPts[0] = endUp;
                querySurf.findLine(startPts, endPts, hitResults);

                bool hitU = !hitResults.empty() && hitResults[0].hit();
                point hitPtU;
                if (hitU) 
                {
                    hitPtU = hitResults[0].hitPoint();
                }

                if (hitD && hitU)
                {
                    scalar dzD = mag(hitPtD.z() - cellPos.z());
                    scalar dzU = mag(hitPtU.z() - cellPos.z());
                    terrainZ = (dzD < dzU) ? hitPtD.z() : hitPtU.z();
                }
                else if (hitD)
                {
                    terrainZ = hitPtD.z();
                }
                else if (hitU)
                {
                    terrainZ = hitPtU.z();
                }

                if (terrainZ > -GREAT)
                {
                    h0[celli] = terrainZ;
                }
                else
                {
                    h0[celli] = 0.0;
                }
            }
        }
        else
        {
            FatalErrorInFunction << "Choose correct mode" << abort(FatalError);
        }
        
        const scalar depth = setSWDict.get<scalar>("depth");
        
        forAll(h, celli)
        {         
            if(depth > h0[celli])
            {
                h[celli] = depth - h0[celli];
            }
        }
    }
    
    h.write();
    h0.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
