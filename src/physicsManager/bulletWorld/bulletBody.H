#ifndef bulletBody_H
#define bulletBody_H

#include "fvCFD.H"
#include "triSurface.H"
#include "btBulletDynamicsCommon.h"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class HLLCS Declaration
\*---------------------------------------------------------------------------*/

class bulletBody
{
    // Private data
    const dictionary& bodyDict_;
    
    bool dynamic_;
    
    word shape_; //box, sphere, cylinder...
    
    word name_;
    
    bool writeVTK_;

    bool writeArbitraryVTK_;

    PtrList<Foam::triSurface> triSurfaces_;
    wordList triSurfaceNames_;
    
    //- Data for Motion//- Current position of the centre of mass of the body
        vector initialP_;
        
        vector P_;
        vector P0_;
        
        tensor Q_;
        tensor Q0_;
        
        btQuaternion q_;
        btQuaternion q0_;
        
        vector localCoM_;
        vector localCoM0_;

        //- Linear velocity of body
        vector v_;
        vector v0_;
        
        vector ev_;
        vector ev0_;

        //- Total linear acceleration of body
        vector a_;

        //- Angular momentum of body, in body local reference frame
        vector pi_;
        vector pi0_;

        //- Total torque on body, in body local reference frame
        vector tau_;
        
        scalar rhoSolid_;
        
        scalar massSolid_;
        
        tensor inertiaSolid_;
        
        vector CoMSolid_; // value at local origin
        
        scalar massFluid_;
        
        tensor inertiaFluid_;
        
        tensor inertiaTotal_;
        tensor inertiaTotal0_;
        
        vector CoMFluid_;
        
        vector planeConstraint_ = vector(1, 1, 1);
        
        vector angularConstraint_ = vector(1, 1, 1);
        
        bool firstIter_ = true;
            
public:

    // Constructors
    bulletBody
    (
        const dictionary& bodyDict
    );


    //- Destructor
    ~bulletBody()
    {}
    
    // Member Functions
    // Access
    inline const dictionary& bodyDict() const { return bodyDict_; }

    inline const bool& dynamic() const { return dynamic_; }

    inline const word& shape() const { return shape_; }

    inline const word& name() const { return name_; }
    
    inline const bool& writeVTK() const { return writeVTK_; }

    inline const bool& writeArbitraryVTK() const { return writeArbitraryVTK_; }

    inline const PtrList<Foam::triSurface>& triSurfaces() const { return triSurfaces_; }
    inline const wordList& triSurfaceNames() const { return triSurfaceNames_; }

    inline const vector& initialP() const { return initialP_; }

    inline const vector& P() const { return P_; }
    inline const vector& P0() const { return P0_; }

    inline const tensor& Q() const { return Q_; }
    inline const tensor& Q0() const { return Q0_; }
    
    inline const btQuaternion& q() const { return q_; }
    inline const btQuaternion& q0() const { return q0_; }
    
    inline const vector& localCoM() const { return localCoM_; }
    inline const vector& localCoM0() const { return localCoM0_; }

    inline const vector& v() const { return v_; }
    inline const vector& v0() const { return v0_; }
    
    inline const vector& ev() const { return ev_; }
    inline const vector& ev0() const { return ev0_; }

    inline const vector& a() const { return a_; }

    inline const vector& pi() const { return pi_; }
    inline const vector& pi0() const { return pi0_; }

    inline const vector& tau() const { return tau_; }

    inline const scalar& rhoSolid() const { return rhoSolid_; }

    inline const scalar& massSolid() const { return massSolid_; }

    inline const tensor& inertiaSolid() const { return inertiaSolid_; }

    inline const vector& CoMSolid() const { return CoMSolid_; }

    inline const scalar& massFluid() const { return massFluid_; }

    inline const tensor& inertiaFluid() const { return inertiaFluid_; }
    
    inline const tensor& inertiaTotal() const { return inertiaTotal_; }
    inline const tensor& inertiaTotal0() const { return inertiaTotal0_; }

    inline const vector& CoMFluid() const { return CoMFluid_; }

    inline const vector& planeConstraint() const { return planeConstraint_; }

    inline const vector& angularConstraint() const { return angularConstraint_; }

    inline const bool& firstIter() const { return firstIter_; }


    // Edit

    inline word& shape() { return shape_; }

    inline word& name() { return name_; }
    
    inline bool& writeVTK() { return writeVTK_; }

    inline bool& writeArbitraryVTK() { return writeArbitraryVTK_; }

    inline PtrList<Foam::triSurface>& triSurfaces() { return triSurfaces_; }
    inline wordList& triSurfaceNames() { return triSurfaceNames_; }

    inline bool& dynamic() { return dynamic_; }

    inline vector& initialP() { return initialP_; }

    inline vector& P() { return P_; }
    inline vector& P0() { return P0_; }

    inline tensor& Q() { return Q_; }
    inline tensor& Q0() { return Q0_; }
    
    inline btQuaternion& q() { return q_; }
    inline btQuaternion& q0() { return q0_; }

    inline vector& localCoM() { return localCoM_; }
    inline vector& localCoM0() { return localCoM0_; }
    
    inline vector& v() { return v_; }
    inline vector& v0() { return v0_; }
    
    inline vector& ev() { return ev_; }
    inline vector& ev0() { return ev0_; }

    inline vector& a() { return a_; }

    inline vector& pi() { return pi_; }
    inline vector& pi0() { return pi0_; }

    inline vector& tau() { return tau_; }

    inline scalar& rhoSolid() { return rhoSolid_; }

    inline scalar& massSolid() { return massSolid_; }

    inline tensor& inertiaSolid() { return inertiaSolid_; }

    inline vector& CoMSolid() { return CoMSolid_; }

    inline scalar& massFluid() { return massFluid_; }

    inline tensor& inertiaFluid() { return inertiaFluid_; }
    
    inline tensor& inertiaTotal() { return inertiaTotal_; }
    inline tensor& inertiaTotal0() { return inertiaTotal0_; }

    inline vector& CoMFluid() { return CoMFluid_; }

    inline vector& planeConstraint() { return planeConstraint_; }

    inline vector& angularConstraint() { return angularConstraint_; }

    inline bool& firstIter() { return firstIter_; }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

