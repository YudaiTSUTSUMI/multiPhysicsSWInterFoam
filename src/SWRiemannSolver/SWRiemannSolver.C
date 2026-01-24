#include "SWRiemannSolver.H"
#include "processorPolyPatch.H"
// * * * * * * * * * * * * * * * * Private Function  * * * * * * * * * * * * * * //

void Foam::SWRiemannSolver::evaluateFlux
(
    scalar& hFluxL,
    scalar& hFluxR,
    vector& hUFluxL,
    vector& hUFluxR,
    scalar& SLeft,
    scalar& SRight,
    const scalar& hLeft,
    const scalar& hRight,
    const scalar& h0Left,
    const scalar& h0Right,
    const vector& hULeft,
    const vector& hURight,
    const scalar& Nf,
    const scalar& wetdryLeft,
    const scalar& wetdryRight,    
    const vector& Sf,
    const scalar& magSf,
    const scalar& dX,
    const scalar& magg
)
{
    // normal vector
    const vector n = Sf/magSf;
    const vector nt = vector(-n[1], n[0], 0.0);


    // Compute qLeft and qRight (scalar of Velocity)
    const vector ULeft = hULeft / stabilise(hLeft, SMALL);
    const vector URight = hURight / stabilise(hRight, SMALL);
    
    scalar qLeft = ULeft & n;
    scalar qRight = URight & n;
    scalar qtLeft = ULeft & nt;
    scalar qtRight = URight & nt;
    
    label front = 1.0;
    
    // Compute total height and consider front
    const scalar hTotalLeft = hLeft + h0Left;
    const scalar hTotalRight = hRight + h0Right;
    
    if((wetdryLeft == 0.0 && wetdryRight == 1.0 && h0Left > hTotalRight) || (wetdryLeft == 1.0 && wetdryRight == 0.0 && h0Right > hTotalLeft)) //    
    {
        qLeft = 0.0; 
        qRight = 0.0;
        front = 0.0;
    }
    
    // Compute Roe averages
    const scalar hroe = (hLeft + hRight)/2;
        
    const scalar qroe = (Foam::sqrt(hLeft)*qLeft + Foam::sqrt(hRight)*qRight)/stabilise(Foam::sqrt(hLeft) + Foam::sqrt(hRight), SMALL);
    
    const scalar CLeft = Foam::sqrt(magg*hLeft);
    const scalar CRight = Foam::sqrt(magg*hRight);
    const scalar Croe = Foam::sqrt(magg*hroe);
    
    
    // compute mid-value
    const scalar deltah0 = h0Right - h0Left;
    const scalar deltahTotal = hTotalRight - hTotalLeft;
    const scalar hj = pos(deltah0)*hLeft + (1-pos(deltah0))*hRight;
    
    scalar deltah0d = 0.0;
    if(deltah0 >= 0 && h0Right > hTotalLeft)
    {
        deltah0d = hLeft;
    }
    else if(deltah0 < 0 && h0Left > hTotalRight)
    {
        deltah0d = -hRight;
    }
    else
    {
        deltah0d = deltah0; 
    }
    
    
		const scalar T = -magg*(hj-mag(deltah0d)/2)*deltah0d;
		scalar Tmax = 0.0;
		if(mag(magg*hroe*deltah0) > mag(T))
		{
		    Tmax = -magg*hroe*deltah0;
		}
		else
		{
		    Tmax = T;
		}
		
		scalar STopo = 0.0;    
		
		if(deltahTotal*deltah0 >= 0 && qroe*deltah0 > 0)
		{
		    STopo = Tmax;
		}
		else
		{
		    STopo = T;
		}
    
    
    //calculate friction source term
    
		const scalar qmin = min(mag(qLeft), mag(qRight));
		
		// Friction coefficient
		const scalar Cf = magg*Foam::sqr(Nf)/stabilise(Foam::pow(hroe, 0.3333), SMALL);
		scalar SFriction = -front*dX*Cf*qroe*mag(qroe);
		        
		scalar frictionE = mag(SFriction/stabilise(magg*hroe,SMALL));
		scalar kineticE = mag(qroe*qmin/(2*magg)); 
		        
		if(frictionE > kineticE)
		{
		    SFriction = SFriction/stabilise(frictionE, SMALL)*kineticE;
		}
        
    scalar S = STopo + SFriction;
    
    
    // Compute signal speeds for face:
	const scalar SLroe = qroe-Croe;
	const scalar SRroe = qroe+Croe;
	
	SLeft = min(SLroe,qLeft-CLeft);
	SLeft = min(SLeft,qRight-CRight);
	
	SRight = max(SRroe,qLeft+CLeft);
	SRight = max(SRight,qRight+CRight);
	
	if(S != 0) 
	{
	    SLeft = SLroe;
	    SRight = SRroe;
	}
    

    scalar Hh = -S/stabilise(SLroe*SRroe, SMALL);
    vector HhURight = -S/stabilise(SLroe*SRroe, SMALL)*qtLeft*nt;
    vector HhULeft = -S/stabilise(SLroe*SRroe, SMALL)*qtRight*nt;
    
    scalar denomR = hRight*(qRight - SRight) - hLeft*(qLeft - SLeft) + SLeft*Hh;
    scalar denomL = hRight*(qRight - SRight) - hLeft*(qLeft - SLeft) + SRight*Hh;
    
    scalar SmidR = 0.0;
    scalar SmidL = 0.0;
    const scalar eps = 1e-12;
    
    if(mag(denomR) > eps)
    {
        SmidR = (SLeft*hRight*(qRight - SRight) - SRight*hLeft*(qLeft - SLeft) + SRight*SLeft*Hh)/denomR;
    }
    
    if(mag(denomL) > eps)
    {
        SmidL = (SLeft*hRight*(qRight - SRight) - SRight*hLeft*(qLeft - SLeft) + SRight*SLeft*Hh)/denomL;
    }
    
    SmidR = 0.9*SmidR + 0.1*qroe;
    SmidL = 0.9*SmidL + 0.1*qroe;
    
    SmidR = (SmidR + mag(SmidR))/2;
    SmidL = (SmidL - mag(SmidL))/2;
    
    SmidR = neg(-SmidR*qroe)*SmidR;
    SmidL = neg(-SmidL*qroe)*SmidL;

    scalar hmidR = (hLeft*(qLeft - SLeft) - SLeft*Hh)/stabilise(SmidR-SLeft, SMALL);
    scalar hmidL = (hRight*(qRight - SRight) + SRight*Hh)/stabilise(SmidL-SRight, SMALL);
    
    vector hUmidR = hmidR*(SmidR*n + qtLeft*nt);
    vector hUmidL = hmidL*(SmidL*n + qtRight*nt);

	if(SmidR > 0 && SmidL == 0 && hmidR < 0)
	{
		S = SLroe*SRroe*(hRight*(qRight - SRight) + hLeft*(SLeft-qLeft))/stabilise(SLeft, SMALL);
		
		Hh = -S/stabilise(SLroe*SRroe, SMALL);
		HhURight = -S/stabilise(SLroe*SRroe, SMALL)*qtLeft*nt;
		
		scalar denomR = hRight*(qRight - SRight) - hLeft*(qLeft - SLeft) + SLeft*Hh;
		
		if(mag(denomR) > eps)
        {
            SmidR = (SLeft*hRight*(qRight - SRight) - SRight*hLeft*(qLeft - SLeft) + SRight*SLeft*Hh)/denomR;
        }
		else
		{
		    SmidR = 0;
		}
		
		hmidR = (hLeft*(qLeft - SLeft) - SLeft*Hh)/stabilise(SmidR-SLeft, SMALL);
		hUmidR = hmidR*(SmidR*n + qtLeft*nt);
	}
	else if(SmidL < 0 && SmidR == 0 && hmidL < 0)
	{
		S = SLroe*SRroe*(hRight*(qRight - SRight) + hLeft*(SLeft-qLeft))/stabilise(SRight, SMALL);
		
		Hh = -S/stabilise(SLroe*SRroe, SMALL);
		HhULeft = -S/stabilise(SLroe*SRroe, SMALL)*qtRight*nt;
		
		scalar denomL = hRight*(qRight - SRight) - hLeft*(qLeft - SLeft) + SRight*Hh;
		
		if(mag(denomL) > eps)
        {
            SmidL = (SLeft*hRight*(qRight - SRight) - SRight*hLeft*(qLeft - SLeft) + SRight*SLeft*Hh)/denomR;
        }
        else
        {
            SmidL = 0;
        }
		hmidL = (hRight*(qRight - SRight) + SRight*Hh)/stabilise(SmidL-SRight, SMALL);
		hUmidL = hmidL*(SmidL*n + qtRight*nt);
	}
	
	const vector ThU = S*Sf;
    
    
    //Calculate Flux
    const scalar hFluxLBase = (hLeft*qLeft)*magSf;
    const scalar hFluxRBase = (hRight*qRight)*magSf;
    const vector hUFluxLBase = (qLeft*hULeft + magg*sqr(hLeft)/2*n)*magSf;
    const vector hUFluxRBase = (qRight*hURight + magg*sqr(hRight)/2*n)*magSf;
    
    if(SLeft >= 0)
    {
        hFluxL = hFluxLBase;
        hFluxR = hFluxLBase;
        hUFluxL = hUFluxLBase;
        hUFluxR = hUFluxLBase + ThU;
    }
    else if(SRight <= 0)
    {
        hFluxL = hFluxRBase;
        hFluxR = hFluxRBase;
        hUFluxL = hUFluxRBase - ThU;
        hUFluxR = hUFluxRBase;
    }
    else if(SLeft < 0 && 0 < SRight && SmidR >= 0 && SmidL == 0)
    {
        hFluxL = hFluxLBase + SLeft*(hmidR - Hh - hLeft)*magSf;
        hFluxR = hFluxL;
        hUFluxL = hUFluxLBase + SLeft*(hUmidR - HhURight - hULeft)*magSf;
        hUFluxR = hUFluxL + ThU;
    }
    else if(SLeft < 0 && 0 < SRight && SmidL < 0 && SmidR == 0)
    {
        hFluxR = hFluxRBase + SRight*(hmidL + Hh - hRight)*magSf;
        hFluxL = hFluxR;
        hUFluxR = hUFluxRBase + SRight*(hUmidL + HhULeft - hURight)*magSf;
        hUFluxL = hUFluxR - ThU;
    }
    else
    {
        FatalErrorInFunction <<  abort(FatalError);
    }
}




// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::SWRiemannSolver::SWRiemannSolver
(
    const volScalarField& h,
    const volScalarField& h0,
    const volVectorField& hU,
    const volScalarField& N,
    const volScalarField& wetdry,
    const surfaceScalarField& irSfLab,
    const vector g
)
:
    mesh_(h.mesh()),
    h_(h),
    h0_(h0),
    hU_(hU),
    N_(N),
    Nf_(linearInterpolate(N)),
    wetdry_(wetdry),
    irSfLab_(irSfLab),
    g_(g),
    hFluxL_
    (
        IOobject
        (
            "hFluxL",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimensionSet(0,3,-1,0,0,0,0),0.0)
    ),
    hFluxR_
    (
        IOobject
        (
            "hFluxR",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimensionSet(0,3,-1,0,0,0,0),0.0)
    ),
    hUFluxL_
    (
        IOobject
        (
            "hUFluxL",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("", dimensionSet(0,4,-2,0,0,0,0), vector(0.0,0.0,0.0))
    ),
    hUFluxR_
    (
        IOobject
        (
            "hUFluxR",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("", dimensionSet(0,4,-2,0,0,0,0), vector(0.0,0.0,0.0))
    ),
    hResidue_
    (
        IOobject
        (
            "hResidue",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimensionSet(0,1,-1,0,0,0,0),0.0)
    ),
    hUResidue_
    (
        IOobject
        (
            "hUResidue",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("", dimensionSet(0,2,-2,0,0,0,0), vector(0.0,0.0,0.0))
    ),
    SLeft_
    (
        IOobject
        (
            "SLeft",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimensionSet(0,1,-1,0,0,0,0), 0.0)
    ),
    SRight_
    (
        IOobject
        (
            "SRight",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimensionSet(0,1,-1,0,0,0,0), 0.0)
    )   
{
    vector dimension = vector(1, 1, 1);

    forAll(g_, i)
    {
        if (mag(g_[i]) > SMALL)
        {
            dimension[i] = 0.0;
        }
    }

    dimension_ = dimension;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::SWRiemannSolver::computeFlux()
{
    // Get the face area vector
    const surfaceVectorField& Sf = mesh_.Sf();
    const surfaceScalarField& magSf = mesh_.magSf();
    
    const_cast<volScalarField&>(h_).correctBoundaryConditions();
    const_cast<volVectorField&>(hU_).correctBoundaryConditions();
    
    // To treat cellCenter in boundary
    volVectorField zeroField
    (
        IOobject
        (
            "",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("", dimensionSet(0,1,0,0,0,0,0), vector(0.0,0.0,0.0))
    );
    const volVectorField cellCenter = zeroField + mesh_.C();
    
    const surfaceVectorField& faceCenter = mesh_.Cf();
    
    // Calculate fluxes at internal faces
    Info << "    Calculate internal flux" << nl;
    forAll(mesh_.owner(), faceI)
    {
        const label own = mesh_.owner()[faceI];
        const label nei = mesh_.neighbour()[faceI];
        
        const scalar dX = mag(cellCenter[own]-cellCenter[nei]);
        
        scalar hLeft = h_[own];
        scalar hRight = h_[nei];
        scalar h0Left = h0_[own];
        scalar h0Right = h0_[nei];
        vector hULeft = hU_[own];
        vector hURight = hU_[nei];
        
        if(mesh_.foundObject<volScalarField>("cellType"))
        {   
            const volScalarField& cellType = mesh_.lookupObject<volScalarField>("cellType");
            const surfaceScalarField& hBound = mesh_.lookupObject<surfaceScalarField>("hBound");
            const surfaceVectorField& hUBound = mesh_.lookupObject<surfaceVectorField>("hUBound");
            
            if(cellType[own] == 1 && cellType[nei] == 2)
            {
                hRight = hBound[faceI];
                hURight = hUBound[faceI];
            }
            if(cellType[own] == 2 && cellType[nei] == 1)
            {
                hLeft = hBound[faceI];
                hULeft = hUBound[faceI];
            }
        }
        
        // calculate fluxes with reconstructed primitive variables at faces
        if(wetdry_[own] == 1.0 || wetdry_[nei] == 1.0)
        {
            if(irSfLab_[faceI] == 0.0)
            {
                evaluateFlux
                (
                    hFluxL_[faceI],  //changed by this function
                    hFluxR_[faceI],  //changed by this function
                    hUFluxL_[faceI], //changed by this function
                    hUFluxR_[faceI], //changed by this function
                    SLeft_[faceI],   //changed by this function
                    SRight_[faceI],  //changed by this function
                    hLeft,
                    hRight,
                    h0Left,
                    h0Right,
                    hULeft,
                    hURight,
                    Nf_[faceI],
                    wetdry_[own],
                    wetdry_[nei],
                    Sf[faceI],
                    magSf[faceI],
                    dX,
                    mag(g_)
                );
            }
            else
            {
                hFluxL_[faceI] *= 0.0;
                hFluxR_[faceI] *= 0.0;
                hUFluxL_[faceI] = (mag(g_)*sqr(h_[own])/2)*Sf[faceI];
                hUFluxR_[faceI] = (mag(g_)*sqr(h_[nei])/2)*Sf[faceI];
                SLeft_[faceI] *= 0.0;
                SRight_[faceI] *= 0.0;
            }
        }
        else
        {
            hFluxL_[faceI] *= 0.0;
            hFluxR_[faceI] *= 0.0;
            hUFluxL_[faceI] *= 0.0;
            hUFluxR_[faceI] *= 0.0;
            SLeft_[faceI] *= 0.0;
            SRight_[faceI] *= 0.0;
        }
    }

    // Update boundary field and values
    Info << "    Calculate boundary flux" << nl;
    forAll (mesh_.boundaryMesh(), patchi)
    {        
        // Fluxes
        fvsPatchField<scalar>& phFluxL  = hFluxL_.boundaryFieldRef()[patchi];
        fvsPatchField<scalar>& phFluxR  = hFluxR_.boundaryFieldRef()[patchi];
        fvsPatchField<vector>& phUFluxL = hUFluxL_.boundaryFieldRef()[patchi];
        fvsPatchField<vector>& phUFluxR = hUFluxR_.boundaryFieldRef()[patchi];
        fvsPatchField<scalar>& pSLeft = SLeft_.boundaryFieldRef()[patchi];
        fvsPatchField<scalar>& pSRight = SRight_.boundaryFieldRef()[patchi];

        // Patch fields
        const fvPatchScalarField& ph = h_.boundaryField()[patchi];
        const fvPatchScalarField& ph0 = h0_.boundaryField()[patchi];
        const vectorField& phU = hU_.boundaryField()[patchi];
        const fvsPatchScalarField& pNf = Nf_.boundaryField()[patchi];
        const fvsPatchScalarField& pirSfLab = irSfLab_.boundaryField()[patchi];

        // Face areas
        const fvsPatchVectorField& pSf = Sf.boundaryField()[patchi];
        const fvsPatchScalarField& pmagSf = magSf.boundaryField()[patchi];
        
        {
        	if (mesh_.boundaryMesh()[patchi].type() == "empty") continue;

            scalarField phLeft  =
                h_.boundaryField()[patchi].patchInternalField();

            scalarField phRight =
                h_.boundaryField()[patchi];//.patchNeighbourField();
                
            const scalarField& ph0Left  =
                h0_.boundaryField()[patchi].patchInternalField();

            const scalarField& ph0Right =
                h0_.boundaryField()[patchi];//.patchNeighbourField();

            vectorField phULeft  =
                hU_.boundaryField()[patchi].patchInternalField();

            vectorField phURight =
                hU_.boundaryField()[patchi];//.patchNeighbourField();

            const scalarField& pwetdryLeft  =
                wetdry_.boundaryField()[patchi].patchInternalField();

            const scalarField& pwetdryRight =
                wetdry_.boundaryField()[patchi];//.patchNeighbourField();
                
            const scalarField pdX =
                2*mag(faceCenter.boundaryField()[patchi] - cellCenter.boundaryField()[patchi].patchInternalField());
                
            if(mesh_.foundObject<volScalarField>("cellType"))
            {
                const scalarField& faceType = mesh_.lookupObject<surfaceScalarField>("faceType").boundaryField()[patchi];
                const scalarField& hBound = mesh_.lookupObject<surfaceScalarField>("hBound").boundaryField()[patchi];
                const vectorField& hUBound = mesh_.lookupObject<surfaceVectorField>("hUBound").boundaryField()[patchi];
                
                forAll(ph, facei)
                {
                    if(faceType[facei] == 1)
                    {   
                        phRight[facei] = hBound[facei];
                        phURight[facei] = hUBound[facei];
                        
                        //phLeft[facei] = hBound[facei]; //optional?
                        //phULeft[facei] = hUBound[facei]; //optional?
                    }
                }
            }
                
            forAll(ph, facei)
            {
                if(pwetdryLeft[facei] == 1.0 || pwetdryRight[facei] == 1.0)
                {
                    if(pirSfLab[facei] == 0.0)
                    {                
                        evaluateFlux
                        (
                            phFluxL[facei],  
                            phFluxR[facei],  
                            phUFluxL[facei],
                            phUFluxR[facei], 
                            pSLeft[facei],   
                            pSRight[facei],  
                            phLeft[facei],
                            phRight[facei],
                            ph0Left[facei],
                            ph0Right[facei],
                            phULeft[facei],
                            phURight[facei],
                            pNf[facei],
                            pwetdryLeft[facei],
                            pwetdryRight[facei],
                            pSf[facei],
                            pmagSf[facei],
                            pdX[facei],
                            mag(g_)
                    	);
                    }
                    else
                    {
                        phFluxL[facei] *= 0.0;
                        phFluxR[facei] *= 0.0;  
                        phUFluxL[facei] = (mag(g_)*sqr(phLeft[facei])/2)*pSf[facei];
                        phUFluxR[facei] = (mag(g_)*sqr(phRight[facei])/2)*pSf[facei];
                        pSLeft[facei] *= 0.0;
                        pSRight[facei] *= 0.0;
                    }
                }   
                else
                {
                    phFluxL[facei] *= 0.0;
                    phFluxR[facei] *= 0.0;
                    phUFluxL[facei] *= 0.0;
                    phUFluxR[facei] *= 0.0;
                    pSLeft[facei] *= 0.0;
                    pSRight[facei] *= 0.0;
                }                          
            }
            
        }
    }
    
    
    computeDeltaT();
}

void Foam::SWRiemannSolver::computeDeltaT()
{
	const scalarField& volume = mesh_.V();
	
	scalarField dXown(mesh_.magSf().size(),VGREAT);
	scalarField dXnei(mesh_.magSf().size(),VGREAT);
	
	
	forAll(mesh_.owner(), facei)
	{
		const int own = mesh_.owner()[facei];
		const int nei = mesh_.neighbour()[facei];

		dXown[facei]=volume[own]/mesh_.magSf()[facei];
		dXnei[facei]=volume[nei]/mesh_.magSf()[facei];
	}
	
	scalarField maxS = max(mag(SLeft()),mag(SRight()));
	forAll(maxS,i)
	{
	    if(maxS[i] < SMALL)
	    {
	    	maxS[i] = SMALL;
	    }
	}

	scalarField deltaTField = min(dXown,dXnei)/maxS;
	
	deltaT_ = gMin(deltaTField);
}

void Foam::SWRiemannSolver::computeResidue()
{
    // Get the cell volume
    const scalarField& cellVolume = mesh_.V();
    
    // Initialize residues
    hResidue_ = 0.0*hResidue_;
    hUResidue_ = 0.0*hUResidue_;

    // Calculate residues from internal fluxes
    forAll(mesh_.owner(), faceI)
    {
        const label own = mesh_.owner()[faceI];
        const label nei = mesh_.neighbour()[faceI];
        
        hResidue_[own] += hFluxL()[faceI]/cellVolume[own];
        hResidue_[nei] -= hFluxR()[faceI]/cellVolume[nei];
        
        vector hUFluxL = applyDimension(hUFluxL_[faceI], dimension_);
        vector hUFluxR = applyDimension(hUFluxR_[faceI], dimension_);
        
        hUResidue_[own] += hUFluxL/cellVolume[own];
        hUResidue_[nei] -= hUFluxR/cellVolume[nei];
        
    }
    
    // Calculate residues from boundary fluxes
    forAll(mesh_.boundaryMesh(), patchi)
    {
        if (mesh_.boundaryMesh()[patchi].type() == "empty") continue;
        
        const fvsPatchScalarField& phFluxL  = hFluxL().boundaryField()[patchi];
        const fvsPatchVectorField& phUFluxL = hUFluxL().boundaryField()[patchi];
        
        const UList<label> &bfaceCells=mesh_.boundaryMesh()[patchi].faceCells();
        
        forAll(phFluxL, facei)
        {
            hResidue_[bfaceCells[facei]] += phFluxL[facei]/cellVolume[bfaceCells[facei]];
            
            vector applyFlux = applyDimension(phUFluxL[facei], dimension_);
        
            hUResidue_[bfaceCells[facei]] += applyFlux/cellVolume[bfaceCells[facei]];
        }
    }
    
    if(mesh_.foundObject<volScalarField>("cellType"))
    {
        const volScalarField& cellType = mesh_.lookupObject<volScalarField>("cellType");
        
        forAll(cellType, celli)
        {
            if(cellType[celli] > 1)
            {
                hResidue_[celli] *= 0;
                hUResidue_[celli] *= 0;
            }
        }
    }   
}

// ************************************************************************* //
