#ifndef bulletWorld_H
#define bulletWorld_H

#include "fvCFD.H"
#include "bulletBody.H"
#include "btBulletDynamicsCommon.h"
#include "btBulletCollisionCommon.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class HLLCS Declaration
\*---------------------------------------------------------------------------*/

class bulletWorld
{
    // Private data
    const dictionary dict_;
    PtrList<bulletBody> bulletBodies_;
    btDiscreteDynamicsWorld* dynamicsWorld_;
    btAlignedObjectArray<btCollisionShape*> collisionShapes_;
    int numOfBulletBodies_;
    int outputCounter_;
    
    // Private function
    void createEmptyBulletWorld(const vector& gValue = vector(0, 0, -9.81), const word& frictionMode = "");
    
    void setRigidBodyConditions
    (
        const dictionary& dict,
        btRigidBody* body,
        bulletBody& btBody
    );
    
    void setInitVelocity
    (
        const dictionary& dict,
        btRigidBody* body
    );
    
    void setConstraint
    (
        const dictionary& dict,
        btRigidBody* body,
        bulletBody& btBody
    );
    
    void addBoxBody
    (
        const dictionary& dict,
        btAlignedObjectArray<btCollisionShape*> collisionShapes,
        bulletBody& btBody
    );
    
    void addSphereBody
    (
        const dictionary& dict,
        btAlignedObjectArray<btCollisionShape*> collisionShapes,
        bulletBody& btBody
    );
    
    void addCylinderZBody
    (
        const dictionary& dict,
        btAlignedObjectArray<btCollisionShape*> collisionShapes,
        bulletBody& btBody
    );
    
    void addArbitraryFloatingBody
    (
        const dictionary& dict,
        btAlignedObjectArray<btCollisionShape*> collisionShapes,
        bulletBody& btBody
    );
    
    void addCompoundBody
    (
        const dictionary& dict,
        btAlignedObjectArray<btCollisionShape*> collisionShapes,
        bulletBody& btBody
    );
    
    void exportCompoundShapeAsVTK
    (
        const btRigidBody* body, 
        const word& filename
    ); 
    
        // parameters
        btVector3 btVector(const vector& foamVector)
        {
            const scalar foamVectorX = foamVector[0];
            const scalar foamVectorY = foamVector[1];
            const scalar foamVectorZ = foamVector[2];
            
            return btVector3(foamVectorX, foamVectorY, foamVectorZ); 
        }
        
        vector foamVector(const btVector3& btVector)
        {
            const scalar btVectorX = btVector[0];
            const scalar btVectorY = btVector[1];
            const scalar btVectorZ = btVector[2];
            
            return vector(btVectorX, btVectorY, btVectorZ); 
        }
        
        btMatrix3x3 btMatrix(const tensor& foamTensor)
        {
            const scalar xx = foamTensor[0];
            const scalar xy = foamTensor[1];
            const scalar xz = foamTensor[2];
            const scalar yx = foamTensor[3];
            const scalar yy = foamTensor[4];
            const scalar yz = foamTensor[5];
            const scalar zx = foamTensor[6];
            const scalar zy = foamTensor[7];
            const scalar zz = foamTensor[8];
            
            return btMatrix3x3(xx,xy,xz,yx,yy,yz,zx,zy,zz); 
        }
        
        tensor foamTensor(const btMatrix3x3& btMatrix)
        {
            const scalar xx = btMatrix[0][0];
            const scalar xy = btMatrix[0][1];
            const scalar xz = btMatrix[0][2];
            const scalar yx = btMatrix[1][0];
            const scalar yy = btMatrix[1][1];
            const scalar yz = btMatrix[1][2];
            const scalar zx = btMatrix[2][0];
            const scalar zy = btMatrix[2][1];
            const scalar zz = btMatrix[2][2];
            
            return tensor(xx,xy,xz,yx,yy,yz,zx,zy,zz); 
        }
    
            
public:
    
    bulletWorld(const vector& gValue = vector(0, 0, -9.81), const word& frictionMode = "");

    //- Destructor
    ~bulletWorld()
    {}
    
    const PtrList<bulletBody>& bulletBodies() const { return bulletBodies_; }
        
    void addBulletBodies
    (
        const dictionary& dict
    );
    
    void setHingeSelf
    (
        const dictionary& jointDict
    );
    
    void setHingePair
    (
        const dictionary& jointDict
    );
    
    void setGeneric6DoFJoint
    (
        const dictionary& jointDict
    );
    
    void setSolidProperties
    (
        const int i,
        const scalar& massSolid,
        const tensor& inertiaSolid,
        const vector& CoMSolid
    );
    
    void setFluidProperties
    (
        const int i,
        const scalar& massFluid,
        const tensor& inertiaFluid,
        const vector& CoMFluid
    );
    
    void reinitializeMassProperties
    (
        const int i,
        const bool implicitEccentric
    );    
    
    void updateFluidPropertiesAndAngularVelocity
    (
        const int i,
        const tensor& inertiaFluid,
        const vector& CoMFluid,
        const scalar& deltaT        
    );
    
    void applyForceTorque
    (
        const int i,
        const vector& force,
        const vector& torque
    );
    
    void stepSimulation
    (
        const scalar& deltaT,
        const int& maxSubSteps,
        const scalar& fixedTimeStep
    );
    
    void storeStates();
    
    void restoreStates();
    
    void writeVTK
    (
        const int outputCounter
    );
    
    btVector3 rotatePoint(const btVector3& point, const btTransform& transform)
    {
        return transform * point;
    }
    
    void writeBoxVTK
    (
        const word& objName,
        const int& objI,
        const int& outputCounter
    );
    
    void writeSphereVTK
    (
        const word& objName,
        const int& objI,
        const int& outputCounter
    );
    
    void writeCylinderZVTK
    (
        const word& objName,
        const int& objI,
        const int& outputCounter
    );

    void writeArbitraryVTK
    (
        const word& objName,
        const int& objI,
        const int& outputCounter
    );
    
    // Access        
    const PtrList<bulletBody>& bulletBodies()
    {
        return bulletBodies_;
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

