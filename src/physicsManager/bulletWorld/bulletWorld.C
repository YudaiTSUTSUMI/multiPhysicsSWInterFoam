#include "bulletWorld.H"

// * * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //


void Foam::bulletWorld::createEmptyBulletWorld(const vector& gValue, const word& frictionMode)
{
    collisionConfiguration_ = new btDefaultCollisionConfiguration();

    dispatcher_ = new btCollisionDispatcher(collisionConfiguration_);

    overlappingPairCache_ = new btDbvtBroadphase();

    solver_ = new btSequentialImpulseConstraintSolver;
    
    dynamicsWorld_ = new btDiscreteDynamicsWorld(dispatcher_, overlappingPairCache_, solver_, collisionConfiguration_);
    
    dynamicsWorld_->setGravity(btVector(gValue));
    
    FrictionSetting::Mode = COMBINE_MULTIPLY;
    FrictionSetting::SetModeFromString(frictionMode);
    
    rmDir("Bullet");
    mkDir("Bullet");
}


void Foam::bulletWorld::setRigidBodyConditions
(
    const dictionary& dict,
    btRigidBody* body,
    bulletBody& btBody
)
{
    //rigid body opptions
    const scalar restitution = dict.lookupOrDefault<scalar>("restitution", 0);
    body->setRestitution(restitution);
    
    const scalar friction = dict.lookupOrDefault<scalar>("friction", 0);
    body->setFriction(friction);
    
    if (dict.found("rollingFriction")) 
    {
        const scalar rollingFriction = dict.get<scalar>("rollingFriction");
        body->setRollingFriction(rollingFriction);
    }

    if (dict.found("spinningFriction")) 
    {
        const scalar spinningFriction = dict.get<scalar>("spinningFriction");
        body->setSpinningFriction(spinningFriction);
    }
    
    body->setActivationState(DISABLE_DEACTIVATION);
    //body->setSleepingThresholds(0.01f, 0.01f);
            
    setInitVelocity(dict, body);
        
    // constraint options 
    setConstraint(dict, body, btBody);
}


void Foam::bulletWorld::setInitVelocity
(
    const dictionary& dict,
    btRigidBody* body
)
{
        const vector initVel = dict.lookupOrDefault<vector>("initialVelocity", vector(Zero));
        if(mag(initVel) > 0)
        {
            body->setEccentricLinearVelocity(btVector(initVel));
            body->setLinearVelocity(btVector(initVel));
        }
        
        const vector initAngVel = dict.lookupOrDefault<vector>("initialAngularVelocity", vector(Zero));
        if(mag(initAngVel) > 0)
        {
            body->setAngularVelocity(btVector(initAngVel));
        }
}


void Foam::bulletWorld::setConstraint
(
    const dictionary& dict,
    btRigidBody* body,
    bulletBody& btBody
)
{
    if (dict.found("constraints"))
    {
        const dictionary& constraintDict = dict.subDict("constraints");
      
        for (const entry& dEntry : constraintDict)
        {
            if (dEntry.isDict())
            {
                const dictionary constraintSubDict = dEntry.dict();
                
                word constraintType;
                constraintSubDict.lookup("bulletMotionConstraint") >> constraintType;
                
                if(constraintType == "plane")
                {
                    const vector normal = constraintSubDict.get<vector>("normal");
                    
                    body -> setLinearFactor ( btVector(normal)); 
                    btBody.planeConstraint() = normal;
                }
                else if(constraintType == "axis")
                {
                    const vector axis = constraintSubDict.get<vector>("axis");
                    
                    body -> setAngularFactor ( btVector(axis)); 
                    btBody.angularConstraint() = axis;
                }
            }
        }
    }
}


void Foam::bulletWorld::addBoxBody
(
    const dictionary& dict,
    btAlignedObjectArray<btCollisionShape*> collisionShapes,
    bulletBody& btBody
)
{
    const scalar Lx = dict.get<scalar>("Lx");
    const scalar Ly = dict.get<scalar>("Ly");
    const scalar Lz = dict.get<scalar>("Lz");
    
    const vector position = dict.get<vector>("position");
    
    btBody.P() = position;
    btBody.initialP() = position;
    
    const scalar eLength = dict.lookupOrDefault<scalar>("expansionLength", 0);
    const scalar rhoSolid = dict.get<scalar>("rhoSolid");
    btBody.rhoSolid() = rhoSolid;
    
    const scalar initialMass = rhoSolid*Lx*Ly*Lz;
            
    // add body to bullet world
    btCollisionShape* boxShapeCollision = new btBoxShape(btVector3(Lx/2+eLength, Ly/2+eLength, Lz/2+eLength));
    
    collisionShapes.push_back(boxShapeCollision);
    btCollisionShape* boxShape = new btBoxShape(btVector3(Lx/2, Ly/2, Lz/2));

    btTransform startTransform;
    startTransform.setIdentity();

    btVector3 initialLocalInertia(0, 0, 0);
    
    if(rhoSolid>0)
    {
        boxShape->calculateLocalInertia(initialMass, initialLocalInertia);
    }

    startTransform.setOrigin(btVector(position));

    btDefaultMotionState* motionState = new btDefaultMotionState(startTransform);
    btRigidBody::btRigidBodyConstructionInfo rbInfo(initialMass, motionState, boxShapeCollision, initialLocalInertia);
    btRigidBody* body = new btRigidBody(rbInfo);
    
    setRigidBodyConditions(dict, body, btBody);
    
    btBody.massSolid() = initialMass;
    btBody.dynamic() = btBody.massSolid() > 0 ? true : false;
    btBody.inertiaSolid()[0] = initialLocalInertia[0];
    btBody.inertiaSolid()[4] = initialLocalInertia[1];
    btBody.inertiaSolid()[8] = initialLocalInertia[2];

	btBody.inertiaTotal()[0] = initialLocalInertia[0];
    btBody.inertiaTotal()[4] = initialLocalInertia[1];
    btBody.inertiaTotal()[8] = initialLocalInertia[2];
    
    //add body to world
    dynamicsWorld_->addRigidBody(body);
    
    // get rotation
    btQuaternion q = body->getOrientation();   
    btMatrix3x3 m; 
    m.setRotation(q);
        
    btBody.Q() = foamTensor(m);
}


void Foam::bulletWorld::addSphereBody
(
    const dictionary& dict,
    btAlignedObjectArray<btCollisionShape*> collisionShapes,
    bulletBody& btBody
)
{
    const scalar R = dict.get<scalar>("R");
    
    const vector position = dict.get<vector>("position");
    
    btBody.P() = position;
    btBody.initialP() = position;
            
    const scalar eLength = dict.lookupOrDefault<scalar>("expansionLength", 0);
    const scalar rhoSolid = dict.get<scalar>("rhoSolid");
    btBody.rhoSolid() = rhoSolid;
    
    const scalar initialMass = rhoSolid * (4.0/3) * M_PI * pow(R,3);
            
    // add body to bullet world	    
    btCollisionShape* sphereShapeCollision = new btSphereShape(R+eLength);
    
    collisionShapes.push_back(sphereShapeCollision);
    btCollisionShape* sphereShape = new btSphereShape(R);

    /// Create Dynamic Objects
    btTransform startTransform;
    startTransform.setIdentity();

    btVector3 initialLocalInertia(0, 0, 0);
    sphereShape->calculateLocalInertia(initialMass, initialLocalInertia);

    startTransform.setOrigin(btVector(position));

    btDefaultMotionState* motionState = new btDefaultMotionState(startTransform);
    btRigidBody::btRigidBodyConstructionInfo rbInfo(initialMass, motionState, sphereShapeCollision, initialLocalInertia);
    btRigidBody* body = new btRigidBody(rbInfo);
    
    setRigidBodyConditions(dict, body, btBody);
        
    btBody.massSolid() = initialMass;
    btBody.dynamic() = btBody.massSolid() > 0 ? true : false;
    btBody.inertiaSolid()[0] = initialLocalInertia[0];
    btBody.inertiaSolid()[4] = initialLocalInertia[1];
    btBody.inertiaSolid()[8] = initialLocalInertia[2];

	btBody.inertiaTotal()[0] = initialLocalInertia[0];
    btBody.inertiaTotal()[4] = initialLocalInertia[1];
    btBody.inertiaTotal()[8] = initialLocalInertia[2];
    
    //add body to world
    dynamicsWorld_->addRigidBody(body);
    
    // get rotation
    btQuaternion q = body->getOrientation();   
    btMatrix3x3 m; 
    m.setRotation(q);
        
    btBody.Q() = foamTensor(m);
}


void Foam::bulletWorld::addCylinderBody
(
    const dictionary& dict,
    btAlignedObjectArray<btCollisionShape*> collisionShapes,
    bulletBody& btBody
)
{
    const word axis = dict.get<word>("axis");
    const scalar radius = dict.get<scalar>("radius");
    const scalar height = dict.get<scalar>("height");    
    
    const vector position = dict.get<vector>("position");
    
    btBody.P() = position;
    btBody.initialP() = position;    
    
    const scalar eLength = dict.lookupOrDefault<scalar>("expansionLength", 0);
    const scalar rhoSolid = dict.get<scalar>("rhoSolid");
    btBody.rhoSolid() = rhoSolid;
    
    const scalar initialMass = rhoSolid * height * M_PI * pow(radius,2);
    
    btCollisionShape* cylinderShape = nullptr;
    btCollisionShape* cylinderShapeCollision = nullptr;

    if(axis == "x" || axis == "X")
    {
        cylinderShapeCollision = new btCylinderShapeX(btVector3(height/2 + eLength, radius + eLength, radius + eLength));
        cylinderShape = new btCylinderShapeX(btVector3(height/2, radius, radius));
    }
    else if(axis == "y" || axis == "Y")
    {
        cylinderShapeCollision = new btCylinderShape(btVector3(radius + eLength, height/2 + eLength, radius + eLength));
        cylinderShape = new btCylinderShape(btVector3(radius, height/2, radius));
    }
    else if(axis == "z" || axis == "Z")
    {
        cylinderShapeCollision = new btCylinderShapeZ(btVector3(radius + eLength, radius + eLength, height/2 + eLength));
        cylinderShape = new btCylinderShapeZ(btVector3(radius, radius, height/2));
    }
    else
    {
        FatalErrorInFunction
            << "Invalid axis value: " << axis
            << ". Must be 'x', 'y', or 'z'."
            << exit(FatalError);
    }
    
    /// Create Dynamic Objects
    btTransform startTransform;
    startTransform.setIdentity();

    btVector3 initialLocalInertia(0, 0, 0);
    cylinderShape->calculateLocalInertia(initialMass, initialLocalInertia);
    
    startTransform.setOrigin(btVector(position));

    btDefaultMotionState* motionState = new btDefaultMotionState(startTransform);
    btRigidBody::btRigidBodyConstructionInfo rbInfo(initialMass, motionState, cylinderShapeCollision, initialLocalInertia);
    btRigidBody* body = new btRigidBody(rbInfo);
    
    setRigidBodyConditions(dict, body, btBody);
    
    btBody.massSolid() = initialMass;
    btBody.dynamic() = btBody.massSolid() > 0 ? true : false;
    btBody.inertiaSolid()[0] = initialLocalInertia[0];
    btBody.inertiaSolid()[4] = initialLocalInertia[1];
    btBody.inertiaSolid()[8] = initialLocalInertia[2];

	btBody.inertiaTotal()[0] = initialLocalInertia[0];
    btBody.inertiaTotal()[4] = initialLocalInertia[1];
    btBody.inertiaTotal()[8] = initialLocalInertia[2];
    
    //add body to world
    dynamicsWorld_->addRigidBody(body);
    
    // get rotation
    btQuaternion q = body->getOrientation();   
    btMatrix3x3 m; 
    m.setRotation(q);
        
    btBody.Q() = foamTensor(m);
}


void Foam::bulletWorld::addArbitraryFloatingBody
(
    const dictionary& dict,
    btAlignedObjectArray<btCollisionShape*> collisionShapes,
    bulletBody& btBody
)
{
    const scalar L = 0.0001;;
    
    const vector position = dict.get<vector>("position");
    
    btBody.P() = position;
    btBody.initialP() = position;
    
    const scalar rhoSolid = dict.get<scalar>("rhoSolid");
    btBody.rhoSolid() = rhoSolid;
    
    const scalar initialMass = dict.get<scalar>("initialMass");
            
    // add body to bullet world            
    btCollisionShape* boxShapeCollision = new btBoxShape(btVector3(L/2, L/2, L/2));
    
    collisionShapes.push_back(boxShapeCollision);

    /// Create Dynamic Objects
    btTransform startTransform;
    startTransform.setIdentity();

    btVector3 initialLocalInertia(0, 0, 0);
    boxShapeCollision->calculateLocalInertia(initialMass, initialLocalInertia);

    startTransform.setOrigin(btVector(position));

    btDefaultMotionState* motionState = new btDefaultMotionState(startTransform);
    btRigidBody::btRigidBodyConstructionInfo rbInfo(initialMass, motionState, boxShapeCollision, initialLocalInertia);
    btRigidBody* body = new btRigidBody(rbInfo);
    
    setRigidBodyConditions(dict, body, btBody);
    
    tensor setInertiaFoam = dict.get<tensor>("initialInertia");
    
    btMatrix3x3 setInertia = btMatrix(setInertiaFoam);
    
    body->setMassProps(initialMass, setInertia, btVector3(0.0, 0.0, 0.0));
    
    btBody.massSolid() = initialMass;
    btBody.dynamic() = btBody.massSolid() > 0 ? true : false;
    btBody.inertiaSolid() = setInertiaFoam;
	btBody.inertiaTotal() = setInertiaFoam;
    
    //add body to world
    dynamicsWorld_->addRigidBody(body);
    
    // get rotation
    btQuaternion q = body->getOrientation();   
    btMatrix3x3 m; 
    m.setRotation(q);
        
    btBody.Q() = foamTensor(m);
}


void Foam::bulletWorld::addCompoundBody
(
    const dictionary& dict,
    btAlignedObjectArray<btCollisionShape*> collisionShapes,
    bulletBody& btBody
)
{
    btCompoundShape* compoundShape = new btCompoundShape();
    collisionShapes.push_back(compoundShape);
    
    scalar mass = 0;
    
    const dictionary& components = dict.subDict("compoundComponents");
    
    Info << "   number of components : " << components.size() << endl;
    
    btScalar masses[components.size()];
    
    btTransform transform;
    
    label componentNum = 0;
    
    for (const entry& dEntry : components)
    {
        const dictionary& componentDict = dEntry.dict();
        
        Info << "       Entry Num : " << componentNum << nl;
        
        word shape;
        componentDict.lookup("shape") >> shape;
        
        if(shape == "box")
        {
            const scalar Lx = componentDict.get<scalar>("Lx");
            const scalar Ly = componentDict.get<scalar>("Ly");
            const scalar Lz = componentDict.get<scalar>("Lz");
            Info << "       Box, Lx : " << Lx << ", Ly : " << Ly << ", Lz : " << Lz << nl;
            
            btBoxShape* box = new btBoxShape(btVector3(Lx/2, Ly/2, Lz/2));
            scalar margin = dict.lookupOrDefault<scalar>("margin", 0.000001);
            box->setMargin(margin);
            collisionShapes.push_back(box);
            
            const vector position = componentDict.get<vector>("position");
            Info << "       position : " << position << endl;
            
            transform.setIdentity();
            transform.setOrigin(btVector(position));
            compoundShape->addChildShape(transform, box);
            
            const scalar rhoSolid = componentDict.get<scalar>("rhoSolid");
            Info << "       rhoSolid : " << rhoSolid << nl << endl;
            
            masses[componentNum] = rhoSolid*Lx*Ly*Lz;
            mass += rhoSolid*Lx*Ly*Lz;
            
            componentNum++;
        }
        
    }
    
    Info << "   total mass : " << mass << endl;
    
    btTransform principal;
    btVector3 inertia;
    compoundShape->calculatePrincipalAxisTransform(masses, principal, inertia);
    compoundShape->recalculateLocalAabb();
    
    scalar margin = dict.lookupOrDefault<scalar>("margin", 0.000001);
    compoundShape->setMargin(margin);
    
    // new compound shape
    btCompoundShape* compound2 = new btCompoundShape();
    collisionShapes.push_back(compound2);
    
    Info << "   numChildShape : " << compoundShape->getNumChildShapes() << endl;
    
    for(int i = 0; i < compoundShape->getNumChildShapes(); i++)
    {
        compound2->addChildShape(principal.inverse()*compoundShape->getChildTransform(i), compoundShape->getChildShape(i));
    }
    
    delete compoundShape;

    /// Create Dynamic Objects
    btTransform startTransform;
    startTransform.setIdentity();

    btDefaultMotionState* motionState = new btDefaultMotionState(principal);
    btRigidBody::btRigidBodyConstructionInfo rbInfo(mass, motionState, compound2, inertia);
    btRigidBody* body = new btRigidBody(rbInfo);
            
    btBody.massSolid() = mass;
    btBody.dynamic() = btBody.massSolid() > 0 ? true : false;
    btBody.inertiaSolid()[0] = inertia[0];
    btBody.inertiaSolid()[4] = inertia[1];
    btBody.inertiaSolid()[8] = inertia[2];

	btBody.inertiaTotal()[0] = inertia[0];
    btBody.inertiaTotal()[4] = inertia[1];
    btBody.inertiaTotal()[8] = inertia[2];
    
    setRigidBodyConditions(dict, body, btBody);
    
    //add body to world
    dynamicsWorld_->addRigidBody(body);

    // position
    btVector3 btPosition = body->getCenterOfMassPosition();
    point foamPosition = foamVector(btPosition);
    btBody.P() = foamPosition;
    btBody.initialP() = foamPosition;
    
    // rotation
    btQuaternion q = body->getOrientation();   
    btMatrix3x3 m; 
    m.setRotation(q);
        
    btBody.Q() = foamTensor(m);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::bulletWorld::bulletWorld(const vector& gValue, const word& frictionMode)
:
    numOfBulletBodies_(0),
    outputCounter_(0)
{
    if(Pstream::master())
    {
        createEmptyBulletWorld(gValue, frictionMode);
    }
}


// * * * * * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * * * //


void Foam::bulletWorld::addBulletBodies
(
    const dictionary& dict
)
{
    if(Pstream::master())
    {
        bulletBodies_.append
        (
            new bulletBody(dict)
        );
                    
        numOfBulletBodies_ = bulletBodies_.size();
        label bodyI = numOfBulletBodies_ - 1;
        
        bulletBody& btBody = bulletBodies_[bodyI];
        
        word shape;
        dict.lookup("objShape") >> shape;
        
        word objName;
        dict.lookup("objName") >> objName;
        
        if(shape == "box")
        {
            addBoxBody(dict, collisionShapes_, btBody);
            btBody.shape() = "box";
            btBody.name() = objName;
        }
        else if(shape == "sphere")
        {
            addSphereBody(dict, collisionShapes_, btBody);
            btBody.shape() = "sphere";
            btBody.name() = objName;
        }
        else if(shape == "cylinder")
        {
            addCylinderBody(dict, collisionShapes_, btBody);
            btBody.shape() = "cylinder";
            btBody.name() = objName;
        }
        else if(shape == "compound")
        {
            addCompoundBody(dict, collisionShapes_, btBody);
            btBody.shape() = "compound";
            btBody.name() = objName;
        }
        else if(shape == "arbitraryFloating")
        {
            addArbitraryFloatingBody(dict, collisionShapes_, btBody);
            btBody.shape() = "arbitraryFloating";
            btBody.name() = objName;
        }
        else
        {
            FatalErrorInFunction<< "set valid shape" << abort(FatalError);
        }
               
        btBody.writeVTK() = dict.lookupOrDefault<bool>("writeVTK", true);
        btBody.writeArbitraryVTK() = false;
        if(btBody.writeVTK())
        {
            mkDir("Bullet/" + objName);

            if(dict.found("triSurface"))
            {
                btBody.writeArbitraryVTK() = true;

                const wordList triSurfaceNames(dict.lookup("triSurface"));
    
                btBody.triSurfaces().setSize(triSurfaceNames.size());
                
                btBody.triSurfaceNames().setSize(triSurfaceNames.size());
                forAll(triSurfaceNames, i)
                {
                    btBody.triSurfaceNames()[i] = triSurfaceNames[i].lessExt();
                }

                forAll(triSurfaceNames, i)
                {
                    const fileName triPath = "constant/triSurface" / triSurfaceNames[i];
                    
                    Info << "Loading surface[" << i << "]: " << triPath << endl;
                    
                    btBody.triSurfaces().set
                    (
                        i,
                        new triSurface(triPath)
                    );
                }
            }
        }
    }
}


void Foam::bulletWorld::setHingeSelf
(
    const dictionary& jointDict
)
{
    if(Pstream::master())
    {        
        word objName = jointDict.get<word>("obj");

        int bodyIndex = -1;
        for (label i = 0; i < bulletBodies_.size(); i++)
        {
            if (objName == bulletBodies_[i].name())
            {
                bodyIndex = i;
                break;
            }
        }

        if(bodyIndex == -1)
        {
            FatalErrorInFunction << "Body not found: " << objName << abort(FatalError);
        }

        btVector3 pivotInBody = btVector(jointDict.get<vector>("pivot"));
        btVector3 axisInBody = btVector(jointDict.get<vector>("axis"));

        btRigidBody* body = btRigidBody::upcast(dynamicsWorld_->getCollisionObjectArray()[bodyIndex]);

        btHingeConstraint* selfHinge = new btHingeConstraint(*body, pivotInBody, axisInBody);

        dynamicsWorld_->addConstraint(selfHinge, true);
    }
}


void Foam::bulletWorld::setHingePair
(
    const dictionary& jointDict
)
{
    if(Pstream::master())
    {        
        
        word objAName =  jointDict.get<word>("objA");
        word objBName =  jointDict.get<word>("objB");
        
        int bodyAI = -1;
        int bodyBI = -1;
        
        for (label i = 0; i < bulletBodies_.size(); i++)
        {            
            if(objAName == bulletBodies_[i].name())
            {
                bodyAI = i;
            }
            if(objBName == bulletBodies_[i].name())
            {
                bodyBI = i;
            }
        }
        
        if(bodyAI == -1 || bodyBI == -1)
        {
            FatalErrorInFunction<< "Choose correct body to joint" << abort(FatalError);
        }
        
        btVector3 axisA = btVector(jointDict.get<vector>("axisA")); 
        btVector3 axisB = btVector(jointDict.get<vector>("axisA")); 
        btVector3 pivotA = btVector(jointDict.get<vector>("pivotA"));
        btVector3 pivotB = btVector(jointDict.get<vector>("pivotB"));
                                
        btCollisionObject* bodyA = dynamicsWorld_->getCollisionObjectArray()[bodyAI];
        btRigidBody* BodyA = btRigidBody::upcast(bodyA);
        
        btCollisionObject* bodyB = dynamicsWorld_->getCollisionObjectArray()[bodyBI];
        btRigidBody* BodyB = btRigidBody::upcast(bodyB);
        
        btHingeConstraint* spHingeDynAB = new btHingeConstraint(*BodyA, *BodyB, pivotA, pivotB, axisA, axisB);
        
        dynamicsWorld_->addConstraint(spHingeDynAB, true);
    }
}


void Foam::bulletWorld::setGeneric6DoFJoint
(
    const dictionary& jointDict
)
{
    if(Pstream::master())
    {                
        word objAName = jointDict.get<word>("objA");
        word objBName = jointDict.get<word>("objB");
                
        int bodyAI = -1;
        int bodyBI = -1;
        
        for (label i = 0; i < bulletBodies_.size(); i++)
        {            
            if(objAName == bulletBodies_[i].name())
            {
                bodyAI = i;
            }
            if(objBName == bulletBodies_[i].name())
            {
                bodyBI = i;
            }
        }
        
        if(bodyAI == -1 || bodyBI == -1)
        {
            FatalErrorInFunction<< "Choose correct body to joint" << abort(FatalError);
        }
        
        btTransform frameInA, frameInB;
        frameInA = btTransform::getIdentity();
        frameInB = btTransform::getIdentity();
        frameInA.setOrigin(btVector(jointDict.get<vector>("frameInA")));
        frameInB.setOrigin(btVector(jointDict.get<vector>("frameInB")));
                            
        btVector3 linearLow = btVector(jointDict.lookupOrDefault<vector>("linearLimitLow", vector::zero));
        btVector3 linearHigh = btVector(jointDict.lookupOrDefault<vector>("linearLimitHigh", vector::zero));
        btVector3 angularLow = btVector(jointDict.lookupOrDefault<vector>("angularLimitLow", vector::zero));
        btVector3 angularHigh = btVector(jointDict.lookupOrDefault<vector>("angularLimitHigh", vector::zero));
                                
        btCollisionObject* bodyA = dynamicsWorld_->getCollisionObjectArray()[bodyAI];
        btRigidBody* BodyA = btRigidBody::upcast(bodyA);
        
        btCollisionObject* bodyB = dynamicsWorld_->getCollisionObjectArray()[bodyBI];
        btRigidBody* BodyB = btRigidBody::upcast(bodyB);
                           
        bool useLinearReferenceFrameA = true;
        
        btGeneric6DofConstraint* constraint = new btGeneric6DofConstraint(
            *BodyA, *BodyB, frameInA, frameInB, useLinearReferenceFrameA);
            
        dynamicsWorld_->addConstraint(constraint, true);

        // constraint
        constraint->setLinearLowerLimit(linearLow);
        constraint->setLinearUpperLimit(linearHigh);
        constraint->setAngularLowerLimit(angularLow);
        constraint->setAngularUpperLimit(angularHigh);
    }
}


void Foam::bulletWorld::setSolidProperties
(
    const int i,
    const scalar& massSolid,
    const tensor& inertiaSolid,
    const vector& CoMSolid
)
{
    if(Pstream::master())
    {
        bulletBodies_[i].massSolid() = massSolid;
        bulletBodies_[i].inertiaSolid() = inertiaSolid;
		bulletBodies_[i].inertiaTotal() = inertiaSolid;
        bulletBodies_[i].CoMSolid() = CoMSolid;
    }
}


void Foam::bulletWorld::setFluidProperties
(
    const int i,
    const scalar& massFluid,
    const tensor& inertiaFluid,
    const vector& CoMFluid
)
{
    if(Pstream::master())
    {
        bulletBodies_[i].massFluid() = massFluid;
        bulletBodies_[i].inertiaFluid() = inertiaFluid;
        bulletBodies_[i].CoMFluid() = CoMFluid;
        
        bulletBodies_[i].dynamic() = true;
    }
}


void Foam::bulletWorld::reinitializeMassProperties
(
    const int i,
    const bool implicitEccentric
)
{
    if(Pstream::master())
    {
        bulletBody& props = bulletBodies_[i];
        
        const scalar massSolid = props.massSolid();
        const tensor inertiaSolid = props.inertiaSolid();
        const vector CoMSolid = props.CoMSolid();
        
        const scalar massFluid = props.massFluid();
        const tensor inertiaFluid = props.inertiaFluid();
        const vector CoMFluid = props.CoMFluid();
        
        const scalar massTotal = massSolid + massFluid;
        const vector localCoM = (massSolid*CoMSolid + massFluid*CoMFluid)/massTotal;
        props.localCoM() = localCoM;
                
        btCollisionObject* obj = dynamicsWorld_->getCollisionObjectArray()[i];
        btRigidBody* btBody = btRigidBody::upcast(obj);
        
        if(mag(CoMSolid) > 0 || massFluid > 0) // eccentric
        {            
            btVector3 setCoM = btVector(localCoM); // eccentricity
                scalar x = setCoM[0];
                scalar y = setCoM[1];
                scalar z = setCoM[2];   
            
            tensor component(Zero);
                
                component[0] = sqr(y) + sqr(z);
                component[1] = -y*x;
                component[2] = -z*x;
                component[3] = -x*y;
                component[4] = sqr(x) + sqr(z);
                component[5] = -z*y;
                component[6] = -x*z;
                component[7] = -y*z;
                component[8] = sqr(x) + sqr(y);
            
            tensor setInertiaFoam = inertiaSolid + inertiaFluid - massTotal * component;
            
            props.inertiaTotal() = setInertiaFoam;
            props.inertiaTotal0() = setInertiaFoam;
            
            btMatrix3x3 setInertia = btMatrix(setInertiaFoam);
            
            btBody->setMassProps(massTotal, setInertia, setCoM);
            
            if(implicitEccentric)
            {
                btBody->setImplicitEccentric();
            }
        }
        else
        {
            btScalar setMass = massSolid;
            btVector3 setInertia = btVector3(inertiaSolid[0], inertiaSolid[4], inertiaSolid[8]);
            
            btBody->setMassProps(setMass, setInertia);
        }
    }
}


void Foam::bulletWorld::updateFluidPropertiesAndAngularVelocity
(
    const int i,
    const tensor& inertiaFluid,
    const vector& CoMFluid,
    const scalar& deltaT
)
{
    if(Pstream::master())
    {
        bulletBody& props = bulletBodies_[i];
        
        const scalar massSolid = props.massSolid();
        const tensor inertiaSolid = props.inertiaSolid();
        const vector CoMSolid = props.CoMSolid();
        
        const scalar massFluid = props.massFluid();
        props.inertiaFluid() = inertiaFluid;
        props.CoMFluid() = CoMFluid;
        
        const scalar massTotal = massSolid + massFluid;
        const vector localCoM = (massSolid*CoMSolid + massFluid*CoMFluid)/massTotal;
        props.localCoM() = localCoM;
        
        btVector3 setCoM = btVector(localCoM); // eccentricity
            scalar x = setCoM[0];
            scalar y = setCoM[1];
            scalar z = setCoM[2];   
        
        tensor component(Zero);
            
            component[0] = sqr(y) + sqr(z);
            component[1] = -y*x;
            component[2] = -z*x;
            component[3] = -x*y;
            component[4] = sqr(x) + sqr(z);
            component[5] = -z*y;
            component[6] = -x*z;
            component[7] = -y*z;
            component[8] = sqr(x) + sqr(y);
        
        tensor setInertiaFoam = inertiaSolid + inertiaFluid - massTotal * component;
        props.inertiaTotal() = setInertiaFoam;
        
        btMatrix3x3 setInertia = btMatrix(setInertiaFoam);
        
        btCollisionObject* obj = dynamicsWorld_->getCollisionObjectArray()[i];
        btRigidBody* btBody = btRigidBody::upcast(obj);
        
        btBody->updateEccentricity(setInertia, setCoM, deltaT);
    }
}


void Foam::bulletWorld::applyForceTorque
(
    const int i,
    const vector& force,
    const vector& torque
)
{
    if(Pstream::master())
    {
        btCollisionObject* obj = dynamicsWorld_->getCollisionObjectArray()[i];
        btRigidBody* body = btRigidBody::upcast(obj);
        
        body->applyCentralForce(btVector(force));
        body->applyTorque(btVector(torque));
    }
}


void Foam::bulletWorld::stepSimulation
(
    const scalar& deltaT,
    const int& maxSubSteps,
    const scalar& fixedTimeStep
)
{
    if(Pstream::master())
    {	    
        dynamicsWorld_->stepSimulation(deltaT, maxSubSteps, fixedTimeStep); //call clearForces
        
        for (label i = 0; i < bulletBodies_.size(); i++)
        {
            btCollisionObject* obj = dynamicsWorld_->getCollisionObjectArray()[i];
	        btRigidBody* body = btRigidBody::upcast(obj);
	        
	        bulletBodies_[i].P() = foamVector(body->getCenterOfMassPosition());
	        
	        // rotation
	        btQuaternion q = body->getOrientation();   
	        bulletBodies_[i].q() = q;
	        
            btMatrix3x3 m; 
            m.setRotation(q);
            
            bulletBodies_[i].Q() = foamTensor(m);
            
            // velocity
            btVector3 velocity = body->getLinearVelocity();
            btVector3 eccentricVelocity = body->getEccentric() ? body->getEccentricLinearVelocity() : velocity;
            bulletBodies_[i].a() = (foamVector(velocity) - bulletBodies_[i].v())/deltaT;
            bulletBodies_[i].v() = foamVector(velocity);
            bulletBodies_[i].ev() = foamVector(eccentricVelocity);
            
            // angularVelocity
            btVector3 angularVelocity = body->getAngularVelocity();
            bulletBodies_[i].tau() = (foamVector(angularVelocity) - bulletBodies_[i].pi())/deltaT;
            bulletBodies_[i].pi() = foamVector(angularVelocity);
        }
    }
}


void Foam::bulletWorld::storeStates()
{
    if(Pstream::master())
    {	    
        for (label i = 0; i < bulletBodies_.size(); i++)
        {
            btCollisionObject* obj = dynamicsWorld_->getCollisionObjectArray()[i];
	        btRigidBody* body = btRigidBody::upcast(obj);
            
            bulletBodies_[i].P0() = foamVector(body->getCenterOfMassPosition());
            bulletBodies_[i].q0() = body->getOrientation();
            btMatrix3x3 m; 
            m.setRotation(body->getOrientation());
            bulletBodies_[i].Q0() = foamTensor(m);            
            bulletBodies_[i].v0() = foamVector(body->getLinearVelocity());
            bulletBodies_[i].pi0() = foamVector(body->getAngularVelocity());
            
            if(mag(bulletBodies_[i].CoMSolid()) > 0 || bulletBodies_[i].massFluid() > 0) // eccentric
            {
                //setInertia, setCoM,
                bulletBodies_[i].ev0() = foamVector(body->getEccentricLinearVelocity());
                bulletBodies_[i].localCoM0() = bulletBodies_[i].localCoM();
                bulletBodies_[i].inertiaTotal0() = bulletBodies_[i].inertiaTotal();
            }
        }
    }
}


void Foam::bulletWorld::restoreStates()
{
    if(Pstream::master())
    {	    
        for (label i = 0; i < bulletBodies_.size(); i++)
        {
            btCollisionObject* obj = dynamicsWorld_->getCollisionObjectArray()[i];
	        btRigidBody* body = btRigidBody::upcast(obj);
            
            btVector3 p0 = btVector(bulletBodies_[i].P0());
            btQuaternion q0 = bulletBodies_[i].q0();
            btVector3 v0 = btVector(bulletBodies_[i].v0());
            btVector3 pi0 = btVector(bulletBodies_[i].pi0());
            
            btTransform trans(q0, p0);
            body->setCenterOfMassTransform(trans);
            body->setLinearVelocity(v0);
            body->setAngularVelocity(pi0);
            
            if(mag(bulletBodies_[i].CoMSolid()) > 0 || bulletBodies_[i].massFluid() > 0) // eccentric
            {
                btVector3 ev0 = btVector(bulletBodies_[i].ev0());
                btVector3 localCoM0 = btVector(bulletBodies_[i].localCoM0());
                btMatrix3x3 inertiaTotal0 = btMatrix(bulletBodies_[i].inertiaTotal0());
                
                body->setEccentricProperties(inertiaTotal0, localCoM0);
                
                body->setEccentricLinearVelocity(ev0);
            }
        }
    }
}


void Foam::bulletWorld::writeVTK
(
    const int outputCounter
)
{
    if(Pstream::master())
    {
        for (label i = 0; i < bulletBodies_.size(); i++)
        {
            const word objShape = bulletBodies_[i].shape();
            const word objName = bulletBodies_[i].name();
            const bool writeVTK = bulletBodies_[i].writeVTK();
            const bool arbitrary = bulletBodies_[i].writeArbitraryVTK();
            const bool dynamic = bulletBodies_[i].dynamic();
            
            if((outputCounter == 0 && writeVTK) || (outputCounter > 0 && writeVTK && dynamic))
            {
                if(arbitrary)
                {
                    writeArbitraryVTK(objName, i, outputCounter);
                }
                else
                {
                    if(objShape == "box")
                    {
                        writeBoxVTK(objName, i, outputCounter);
                    }
                    else if(objShape == "sphere")
                    {
                        writeSphereVTK(objName, i, outputCounter);
                    }
                    else if(objShape == "cylinder")
                    {
                        writeCylinderVTK(objName, i, outputCounter);
                    }
                    else if(objShape == "compound")
                    {
                    }
                    else if(objShape == "arbitraryFloating")
                    {                    
                    }
                    else if(objShape == "arbitrary")
                    {
                    }
                }
                
            }
        }
    }
}


void Foam::bulletWorld::writeBoxVTK
(
    const word& objName,
    const int& objI,
    const int& outputCounter
)
{
    btCollisionObject* obj = dynamicsWorld_->getCollisionObjectArray()[objI];
    
    word saveName = "Bullet/" + objName + "/" + objName + "_" +  std::__cxx11::to_string(outputCounter);
    OFstream vtkFile(saveName + ".vtk");
    
    // obtain Box shape
    const btBoxShape* box = dynamic_cast<const btBoxShape*>(obj->getCollisionShape());

    btVector3 halfExtents = box->getHalfExtentsWithMargin();
    btTransform transform = obj->getWorldTransform();

    // calculate 8 points
    std::vector<btVector3> corners;
    for (int x = -1; x <= 1; x += 2) {
        for (int y = -1; y <= 1; y += 2) {
            for (int z = -1; z <= 1; z += 2) {
                btVector3 corner(x * halfExtents.x(), y * halfExtents.y(), z * halfExtents.z());
                corners.push_back(rotatePoint(corner, transform));
            }
        }
    }

    // VTK header
    vtkFile << "# vtk DataFile Version 3.0\n";
    vtkFile << "Rigid body output\n";
    vtkFile << "ASCII\n";
    vtkFile << "DATASET POLYDATA\n";

    // output point
    vtkFile << "POINTS 8 float\n";
    for (const auto& pt : corners) {
        vtkFile << pt.x() << " " << pt.y() << " " << pt.z() << "\n";
    }

    // define face
    vtkFile << "POLYGONS 6 30\n";
    const int faces[6][4] = {
        {0, 1, 3, 2}, {4, 5, 7, 6},
        {0, 1, 5, 4}, {2, 3, 7, 6},
        {0, 2, 6, 4}, {1, 3, 7, 5}
    };
    for (const auto& face : faces) {
        vtkFile << "4 " << face[0] << " " << face[1] << " " << face[2] << " " << face[3] << "\n";
    }
    
    vtkFile << "\n";
}


void Foam::bulletWorld::writeSphereVTK
(
    const word& objName,
    const int& objI,
    const int& outputCounter
)
{
    btCollisionObject* obj = dynamicsWorld_->getCollisionObjectArray()[objI];
    const btSphereShape* sphere = dynamic_cast<const btSphereShape*>(obj->getCollisionShape());

    word saveName = "Bullet/" + objName + "/" + objName + "_" + std::to_string(outputCounter);
    OFstream vtkFile(saveName + ".vtk");

    btScalar radius = sphere->getRadius();
    btTransform transform = obj->getWorldTransform();

    const int stacks = 20;
    const int slices = 40;
    std::vector<btVector3> points;

    for (int i = 0; i <= stacks; ++i)
    {
        btScalar phi = M_PI * btScalar(i) / stacks;
        for (int j = 0; j <= slices; ++j)
        {
            btScalar theta = 2.0 * M_PI * btScalar(j) / slices;
            btVector3 p(
                radius * sin(phi) * cos(theta),
                radius * sin(phi) * sin(theta),
                radius * cos(phi)
            );
            points.push_back(rotatePoint(p, transform));
        }
    }

    vtkFile << "# vtk DataFile Version 3.0\n";
    vtkFile << "Sphere output\n";
    vtkFile << "ASCII\n";
    vtkFile << "DATASET POLYDATA\n";
    vtkFile << "POINTS " << points.size() << " float\n";
    for (const auto& pt : points)
    {
        vtkFile << pt.x() << " " << pt.y() << " " << pt.z() << "\n";
    }
    
    vtkFile << "POLYGONS " << stacks * slices * 2 << " " << stacks * slices * 2 * 4 << "\n";
    for (int i = 0; i < stacks; ++i)
    {
        for (int j = 0; j < slices; ++j)
        {
            int p1 = i * (slices + 1) + j;
            int p2 = p1 + slices + 1;

            vtkFile << "3 " << p1 << " " << p1 + 1 << " " << p2 << "\n";
            vtkFile << "3 " << p1 + 1 << " " << p2 + 1 << " " << p2 << "\n";
        }
    }
}


void Foam::bulletWorld::writeCylinderVTK
(
    const word& objName,
    const int& objI,
    const int& outputCounter
)
{
    btCollisionObject* obj = dynamicsWorld_->getCollisionObjectArray()[objI];
    btCollisionShape* shape = obj->getCollisionShape();

    btVector3 halfExtents;
    word axis = "z";

    if (const btCylinderShapeZ* cylinderZ = dynamic_cast<const btCylinderShapeZ*>(shape))
    {
        halfExtents = cylinderZ->getHalfExtentsWithMargin();
        axis = "z";
    }
    else if (const btCylinderShapeX* cylinderX = dynamic_cast<const btCylinderShapeX*>(shape))
    {
        halfExtents = cylinderX->getHalfExtentsWithMargin();
        axis = "x";
    }
    else if (const btCylinderShape* cylinderY = dynamic_cast<const btCylinderShape*>(shape))
    {
        halfExtents = cylinderY->getHalfExtentsWithMargin();
        axis = "y";
    }
    else
    {
        FatalErrorInFunction
            << "Unsupported shape type for cylinder VTK output." << exit(FatalError);
    }

    btScalar radius, height;
    int ax = 0, r1 = 1, r2 = 2; // default z-axis → height: z, radius: x/y

    if (axis == "x")
    {
        height = halfExtents.x() * 2.0;
        radius = halfExtents.y();
        ax = 0; r1 = 1; r2 = 2;
    }
    else if (axis == "y")
    {
        height = halfExtents.y() * 2.0;
        radius = halfExtents.x();
        ax = 1; r1 = 0; r2 = 2;
    }
    else // z
    {
        height = halfExtents.z() * 2.0;
        radius = halfExtents.x();
        ax = 2; r1 = 0; r2 = 1;
    }

    word saveName = "Bullet/" + objName + "/" + objName + "_" + std::to_string(outputCounter);
    OFstream vtkFile(saveName + ".vtk");

    btTransform transform = obj->getWorldTransform();

    const int slices = 40;
    std::vector<btVector3> points;

    for (int i = 0; i < 2; ++i)
    {
        btScalar h = (i == 0 ? -0.5 : 0.5) * height;

        for (int j = 0; j < slices; ++j)
        {
            btScalar theta = 2.0 * M_PI * j / slices;
            btVector3 p(0, 0, 0);
            p[r1] = radius * cos(theta);
            p[r2] = radius * sin(theta);
            p[ax] = h;

            points.push_back(rotatePoint(p, transform));
        }
    }

    // Add center points of bottom and top faces
    btVector3 bottomCenter(0, 0, 0), topCenter(0, 0, 0);
    bottomCenter[ax] = -0.5 * height;
    topCenter[ax] = 0.5 * height;
    points.push_back(rotatePoint(bottomCenter, transform));
    points.push_back(rotatePoint(topCenter, transform));

    vtkFile << "# vtk DataFile Version 3.0\n";
    vtkFile << "Cylinder output\n";
    vtkFile << "ASCII\n";
    vtkFile << "DATASET POLYDATA\n";

    vtkFile << "POINTS " << points.size() << " float\n";
    for (const auto& pt : points)
    {
        vtkFile << pt.x() << " " << pt.y() << " " << pt.z() << "\n";
    }

    // POLYGONS
    vtkFile << "POLYGONS " << slices * 4 << " " << slices * 4 * 4 << "\n";
    for (int j = 0; j < slices; ++j)
    {
        int next = (j + 1) % slices;
        int bottom0 = j;
        int bottom1 = next;
        int top0 = j + slices;
        int top1 = next + slices;

        vtkFile << "3 " << bottom0 << " " << top0 << " " << bottom1 << "\n";
        vtkFile << "3 " << bottom1 << " " << top0 << " " << top1 << "\n";

        vtkFile << "3 " << points.size() - 2 << " " << bottom1 << " " << bottom0 << "\n"; // bottom face
        vtkFile << "3 " << points.size() - 1 << " " << top0 << " " << top1 << "\n";       // top face
    }
}


void Foam::bulletWorld::writeArbitraryVTK
(
    const word& objName,
    const int& objI,
    const int& outputCounter
)
{
    bulletBody& btBody = bulletBodies_[objI];

    forAll(btBody.triSurfaces(), i)
    {
        fileName saveDir = "Bullet"/objName;
        word saveName = btBody.triSurfaceNames()[i] + "_" + Foam::name(outputCounter);
        fileName fullPath = saveDir / saveName + ".vtk";

        Foam::triSurface tri = btBody.triSurfaces()[i]; 
        pointField points = tri.points();

        const vector P = btBody.P();              
        const vector initialP = btBody.initialP();
        const tensor Q = btBody.Q();

        forAll(points, j)
        {
            vector relativePos = points[j] - initialP;
            points[j] = P + (Q & relativePos);
        }

        tri.movePoints(points);

        tri.write(fullPath);
    }
}

// ************************************************************************* //

