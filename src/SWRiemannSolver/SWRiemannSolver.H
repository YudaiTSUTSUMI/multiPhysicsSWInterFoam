#ifndef SWRiemannSolver_H
#define SWRiemannSolver_H

#include "polyMesh.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class SWRiemannSolver Declaration
\*---------------------------------------------------------------------------*/

class SWRiemannSolver
{
    // Private data

        //- Reference to mesh
        const fvMesh& mesh_;


        // Reference to primitive fields

            //- Inundation depth
            const volScalarField& h_;
            
            //- Topology
            const volScalarField& h0_;

            //- Velocity
            const volVectorField& hU_;
            
            //- Manning coefficient
            const volScalarField& N_;
            
            //- Interpolated manning coefficient
            const surfaceScalarField Nf_;

            //- bool of wet or dry
            const volScalarField& wetdry_;
            
            //- bool of dike
            const surfaceScalarField& irSfLab_;
            
            //- magunitude of gravitational acceleration
            vector g_;
            
            vector dimension_;


        // Fluxes

            //- Inundation depth flux
            surfaceScalarField hFluxL_;
            surfaceScalarField hFluxR_;

            //- Velocity flux
            surfaceVectorField hUFluxL_;
            surfaceVectorField hUFluxR_;
            
        // Residues

            //- Residue of h
            volScalarField hResidue_;
            
            //- Residue of hU
            volVectorField hUResidue_;
            
        
        // Wave speed
            surfaceScalarField SLeft_;
            surfaceScalarField SRight_;
            
        // delta T
            scalar deltaT_;
         
    // Private function  
        void evaluateFlux
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
        );
        
        inline vector applyDimension(const vector& v, const vector& dimension) const
        {
            return vector
            (
                v.x() * dimension.x(),
                v.y() * dimension.y(),
                v.z() * dimension.z()
            );
        }

public:

    // Constructors

        //- Construct from components
        SWRiemannSolver
        (
            const volScalarField& h,
            const volScalarField& h0,
            const volVectorField& hU,
            const volScalarField& N,
            const volScalarField& wetdry,
            const surfaceScalarField& irSfLab,
            const vector g
        );


    //- Destructor
    ~SWRiemannSolver()
    {}


    // Member Functions

        //- Return mesh reference
        const fvMesh& mesh() const
        {
            return mesh_;
        }


        // Return fluxes

            //- Return height flux
            const surfaceScalarField hFluxL() const
            {
                return hFluxL_;
            }
            
            surfaceScalarField& hFluxL()
            {
                return hFluxL_;
            }
            
            const surfaceScalarField hFluxR() const
            {
                return hFluxR_;
            }
            
            surfaceScalarField& hFluxR() 
            {
                return hFluxR_;
            }

            //- Return momentum flux
            const surfaceVectorField hUFluxL() const
            {
                return hUFluxL_;
            }
            
            surfaceVectorField& hUFluxL()
            {
                return hUFluxL_;
            }
            
            const surfaceVectorField hUFluxR() const
            {
                return hUFluxR_;
            }
            
            surfaceVectorField& hUFluxR()
            {
                return hUFluxR_;
            }
               
        // Return residues

            //- Return h residue
            const volScalarField& hResidue() const
            {
                return hResidue_;
            }

            //- Return hU residue
            const volVectorField& hUResidue() const
            {
                return hUResidue_;
            }
            
            
        // Return wave speed

            //- Return height residue
            const surfaceScalarField& SLeft() const
            {
                return SLeft_;
            }

            //- Return momentum residue
            const surfaceScalarField& SRight() const
            {
                return SRight_;
            }
            
        // Return delta T 
        const scalar& deltaT() const
        {
            return deltaT_;
        }

        // Update fluxes and friction source turm based on current state

            //- Compute flux
            void computeFlux();
            
            void computeDeltaT();
        
        
        // Update residues based on current state    
        
            //- Compute residue
            void computeResidue();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

