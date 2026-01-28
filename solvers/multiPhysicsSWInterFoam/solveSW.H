forAll(SWRegions, i)
{
	#include "setSWFields.H"     
	
	waveSource.update();
	
	Info << "compute " << SWNames[i] << " SWRiemannSolver flux" << endl;
    riemannSolver.computeFlux();
    		        
    Info << "compute " << SWNames[i] << " SWRiemannSolver residue and update" << endl;
    riemannSolver.computeResidue();
    
    h = h - runTime.deltaT() * riemannSolver.hResidue();
	hU = hU - runTime.deltaT() * riemannSolver.hUResidue();
	
    wetdry = pos(h-hDry);
    h = wetdry*h;
    hU = wetdry*hU;   
    hTotal = h + h0;

    volScalarField h2 = h;
    h2.max(hDry);
    U = hU/h2;
    
    if(limitFrontVelocity)
    {
        volScalarField front = 0*wetdry;

        forAll(mesh.owner(), facei)
        {
            const label own = mesh.owner()[facei];
            const label nei = mesh.neighbour()[facei];

            if(wetdry[own] == 1 && wetdry[nei] == 0)
            {
                front[own] = 1;
            }
            else if(wetdry[own] == 0 && wetdry[nei] == 1)
            {
                front[nei] = 1;
            }
        }

        volScalarField UWet = wetdry*(1-front)*mag(U);
        scalar UMax = gMax(UWet);

        Info << "\n limited velocity: " << UMax << nl << endl; 

        forAll(U, celli)
        { 
            scalar magU = mag(U[celli]);
            if (front[celli] == 1.0 &&  magU > UMax)
            { 
                U[celli] *= UMax / magU;
                hU[celli] = h[celli]*U[celli];
            }
        }
    }
    
    h.correctBoundaryConditions();
    hU.correctBoundaryConditions();
}

