#ifndef physicsManager_H
#define physicsManager_H

#include "fvCFD.H"
#include "bulletOversetFvMesh.H"
#include "bulletMotionSolver.H"

#include "bulletWorld.H"
#include "moorDynWorld.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

struct primitiveFoamData
{
    bulletOversetFvMesh* mesh_;
    bulletMotionSolver* solver_;
};

struct forceData
{
    vector force_ = Zero;
    vector torque_ = Zero;
};

struct solidData
{
    scalar solidMass_ = 0;
    tensor solidInertia_ = Zero;
    vector solidCoMLocal_ = Zero;
};

struct fluidData
{
    scalar fluidMass_ = 0;
    tensor fluidInertia_ = Zero;
    vector fluidCoMLocal_ = Zero;
    vectorField CLocal_;
    const scalarField* V_;
    const scalarField* fluidRho_;
};

// outerForceData
struct outerForceData
{
    primitiveFoamData base_;
    forceData forceData_;
};

// innerForceData
struct innerForceData
{
    primitiveFoamData base_;
    forceData forceData_;
    solidData solidData_;
};

// innerMassData
struct innerMassData
{
    primitiveFoamData base_;
    solidData solidData_;
    fluidData fluidData_;
};

// trackingData
struct trackingData
{
    primitiveFoamData base_;
    vector normal_;
    vector axis_;
    vector initialP_;
    tensor initialQ_;
};


struct mooringData
{
    label mooringID_;
    vector mooringLocalPosition_; // local coordinate
    
    vector mooringForce_ = Zero;
    vector mooringTorque_ = Zero;
};


struct objectData
{
    word objName_;
    
    word objShape_;
    
    scalar fluidForceRelaxation_ = 1.0;
    scalar fluidForceDamping_ = 1.0;
    
    vector totalFluidForce_;
    vector totalFluidTorque_;

    vector totalMooringForce_;
    vector totalMooringTorque_;
    
    dictionary objDict_;
    
    const bulletBody* bulletBody_;
    
    List<outerForceData> outerForceDataList_;
    List<innerForceData> innerForceDataList_;
    List<innerMassData> innerMassDataList_;
    List<trackingData> trackingDataList_;
    
    List<mooringData> mooringDataList_;

    autoPtr<OFstream> stateCSV_;
    autoPtr<OFstream> forceCSV_;
    
    //private menber function
    template<class Type>
    Type getFromMaster(std::function<Type()> valueGetter)
    {
        Type result = pTraits<Type>::zero;

        if (Pstream::master())
        {
            result = valueGetter();
        }

        reduce(result, sumOp<Type>());
        return result;
    }
    
    vector getPosition()
    {
        return getFromMaster<vector>([&]() { return bulletBody_->P(); });
    }
    
    vector getPosition0()
    {
        return getFromMaster<vector>([&]() { return bulletBody_->P0(); });
    }

    tensor getRotation()
    {
        return getFromMaster<tensor>([&]() { return bulletBody_->Q(); });
    }
    
    tensor getRotation0()
    {
        return getFromMaster<tensor>([&]() { return bulletBody_->Q0(); });
    }
    
    vector getVelocity()
    {
        return getFromMaster<vector>([&]() { return bulletBody_->ev(); });
    }
    
    vector getAngularVelocity()
    {
        return getFromMaster<vector>([&]() { return bulletBody_->pi(); });
    }

    scalar getSolidRho()
    {
        return getFromMaster<scalar>([&]() { return bulletBody_->rhoSolid(); });
    }
    
    scalar getSolidMass()
    {
        return getFromMaster<scalar>([&]() { return bulletBody_->massSolid(); });
    }
    
    tensor getSolidInertia()
    {
        return getFromMaster<tensor>([&]() { return bulletBody_->inertiaSolid(); });
    }
    
    vector getSolidCoM()
    {
        return getFromMaster<vector>([&]() { return bulletBody_->CoMSolid(); });
    }
    
    vector getCoMLocal()
    {
        return getFromMaster<vector>([&]() { return bulletBody_->localCoM(); });
    }
    
    vector getCoMWorld()
    {
        const vector P = getPosition();
        const tensor Q = getRotation();
        const vector CoMLocal = getCoMLocal();
        
        return P + (Q & CoMLocal);
    }
    
    bool getDynamic()
    {
        return getFromMaster<bool>([&]() { return bulletBody_->dynamic(); });
    }
};


class physicsManager
{
    // Private data
    const Time& runTime_;
        
    dictionary physicsDict_;
        
    autoPtr<bulletWorld> bulletWorld_;
    
    autoPtr<moorDynWorld> moorDynWorld_;
    
    PtrList<objectData> objectDataList_;
    
    label step_ = 0;
    
    label outputCounter_ = 0;
    
    label nOuterCorrectors_ = 1;
    
    label nIterBullet_ = 10;
    
    bool implicitEccentric_ = false;
    
    bool predictMooringForces_ = false;
    
    //private menber function
    //initialize
    void initBulletWorld(vector gValue);
    
    void initMoorDynWorld();
    
    void setJoints();
    
    //update
    void updateFluidProperties(const scalar& deltaT);
    
    void calculateFluidForces();
    
    void calculateMooringForces();
    
    void applyForces();
    
    void updateSolvers();
    
    void predictMoorDynWorld(const scalar& t, const scalar& deltaT);
    
    void updateMoorDynWorld(const scalar& t, const scalar& deltaT);
    
    //calculate
    void calculateMassProps
    (
        const scalarField& rho,
        const scalarField& V,
        const vectorField& CLocal,
        
        scalar& mass,
        tensor& inertia,
        vector& CoMLocal
    );
    
    void calculateFluidProps
    (
        const scalarField& rho,
        const scalarField& V,
        const vectorField& CLocal,
        const scalar& mass,
        
        tensor& inertia,
        vector& CoMLocal
    );
    
    //Access
    int getObjectID(word objName);

    //Write Data
    void prepareCSV();

    void writeCSV();
    
public:

    // Constructors
    physicsManager(const Time& runTime, vector gValue = vector(0, 0, -9.81));

    //- Destructor
    ~physicsManager()
    {}
    
    // Menber functions
    void update();
    
    void restoreStates();
    
    void writeAndStoreStates();
        
    void registerOuterForce
    (
        const word& objName,
        bulletOversetFvMesh& mesh,
        bulletMotionSolver& solver
    );
    
    void registerInnerForce
    (
        const word& objName,
        bulletOversetFvMesh& mesh,
        bulletMotionSolver& solver
    );
    
    void registerInnerMass
    (
        const word& objName,
        bulletOversetFvMesh& mesh,
        bulletMotionSolver& solver,
        scalarField& rhoFluid
    );
    
    void registerTracking
    (
        const word& objName,
        bulletOversetFvMesh& mesh,
        bulletMotionSolver& solver,
        vector normal,
        vector axis
    );
    
    void reregisterObjects();
    
    const label& nOuterCorrectors() const
    {
        return nOuterCorrectors_;
    }
    
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

