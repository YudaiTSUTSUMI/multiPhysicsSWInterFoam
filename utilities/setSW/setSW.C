#include "fvCFD.H"
#include "pimpleControl.H"
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
