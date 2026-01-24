dimensionedScalar dt("dt", dimensionSet(0,0,1,0,0,0,0), runTime.deltaTValue());

forAll(SWRegions, i)
{
	#include "setSWFields.H"     
	
	waveSource.update();
	
	Info << "compute " << SWNames[i] << " SWRiemannSolver flux" << endl;
    SWRiemannSolver.computeFlux();
    		        
    Info << "compute " << SWNames[i] << " SWRiemannSolver residue and update" << endl;
    SWRiemannSolver.computeResidue();
    
    h = h - dt * SWRiemannSolver.hResidue();
	hU = hU - dt * SWRiemannSolver.hUResidue();
	
	wetdry = pos(h-hSMALL);
    h *= wetdry;
	    
    hTotal == h + h0;
	    
    hU = wetdry*hU;
    volScalarField h2 = h;
    h2.max(hSMALL);
    U = hU/h2;

    // Limiting velocity
    forAll(U,celli)
    {
        if(mag(U[celli]) > Umax.value()) // need to be option
        {
            h[celli] *= 0;
            hU[celli] *= 0;
            U[celli] *= 0;
            wetdry[celli] = 0;
        }
    }

    h.correctBoundaryConditions();
    hU.correctBoundaryConditions();

    wetdryf = 0*wetdryf;
    forAll(wetdryf, facei)
    {
        const int own = mesh.owner()[facei];
        const int nei = mesh.neighbour()[facei];
                 
        if((wetdry[own] == 1 && wetdry[nei] == 0) || (wetdry[own] == 0 && wetdry[nei] == 1))
        {
            wetdryf[facei] = 1; 
        }
        
        if(lsLabel[facei] == 1 && (max(hTotal[own],hTotal[nei]) > max(h0[own],h0[nei])+lsHeight[facei]))
	    {
		    lsLabel[facei] = 0;
	    }
    }
}

