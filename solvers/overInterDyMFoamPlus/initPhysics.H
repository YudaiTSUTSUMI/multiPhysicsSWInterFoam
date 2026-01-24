physicsManager PhysicsManager(runTime, gValue);

List<bool> updateVOFOIDF(OIDFRegions.size(), true);

Info << nl << "Register mesh and solver to physicsManager" << endl;

forAll(OIDFNames, regionI)
{
    Info << nl << "Region: " << OIDFRegions[regionI].name() << nl << endl;
    bulletOversetFvMesh& mesh = OIDFRegions[regionI];
    
    forAll(mesh.bulletMotionSolvers(), solverI)
    {
        dictionary solverDict = mesh.getSolverDict(solverI);
        
        word objName = solverDict.get<word>("objName");
        
        word mode = solverDict.get<word>("mode");
        
        bulletMotionSolver& solver = mesh.bulletMotionSolvers()[solverI];
        
        Info << "SolverI: " << solverI << ", mode: " << mode << endl;
        
        if(mode == "outerForce")
        {
            PhysicsManager.registerOuterForce(objName, mesh, solver);
        }
        else if(mode == "innerForce")   
        {
            PhysicsManager.registerInnerForce(objName, mesh, solver);
        }
        else if(mode == "innerMass")   
        {
            PhysicsManager.registerInnerMass(objName, mesh, solver, rhoOIDF[regionI]);
            
            updateVOFOIDF[regionI] = solverDict.lookupOrDefault<bool>("updateVOF", true);
        }
        else if(mode == "tracking")
        {
            const vector normal = solverDict.lookupOrDefault<vector>("normal", vector(1, 1, 0));
            const vector axis = solverDict.lookupOrDefault<vector>("axis", vector(0, 0, 1));
            
            PhysicsManager.registerTracking(objName, mesh, solver, normal, axis);
        }
        else
        {
            FatalErrorInFunction << "Select correct mode: " << mode << abort(FatalError);
        }
    }
}   

Info << nl << "Update object properties" << endl;
PhysicsManager.reregisterObjects();

