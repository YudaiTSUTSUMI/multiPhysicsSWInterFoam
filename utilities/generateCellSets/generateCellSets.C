#include "fvCFD.H"
#include "cellSet.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createNamedMesh.H"
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
    const word path = mesh.time().constant() + "/polyMesh/sets";
    mkDir(path);
    
    const label nCells = mesh.nCells();
    Info << "totalCells: " << nCells << endl;
    
    label startCelli = 0;
    List<bool> marked(nCells, false);
    
    label regionI = 0;
    
    while(startCelli < nCells)
    {
        DynamicList<label> connectedCells;
        DynamicList<label> targetCells;
        targetCells.append(startCelli);

        while (!targetCells.empty())
        {
            DynamicList<label> nextTargetCells;
            
            forAll(targetCells, i)
            {
                const label cellI = targetCells[i];
                
                const labelList& nbrs = mesh.cellCells()[cellI];
                
                forAll(nbrs, j)
                {
                    const label nbrI = nbrs[j];
                    if (!marked[nbrI])
                    {
                        connectedCells.append(nbrI);
                        nextTargetCells.append(nbrI);
                        
                        marked[nbrI] = true;
                    }
                }
            }
            
            targetCells = nextTargetCells;
        }
        
        Info << "Connected region found with " << connectedCells.size() << " cells." << endl;
        
        //generateSets
        {
            word setName = "c" + Foam::name(regionI);
            
            cellSet newSet
            (
                IOobject
                (
                    setName, 
                    mesh.facesInstance(),
                    polyMesh::meshSubDir/"sets",
                    mesh,             
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                )
            );
            
            forAll(connectedCells, i)
            {
                newSet.insert(startCelli + i);
            }

            newSet.write();
        }
        
        regionI++;
        startCelli = startCelli + connectedCells.size();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
