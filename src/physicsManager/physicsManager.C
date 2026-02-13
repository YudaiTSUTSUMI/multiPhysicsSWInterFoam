#include "physicsManager.H"
#include "IOdictionary.H"
#include "Time.H"

// * * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //
 

// * * * * * * * * * * * * * * * Initialize * * * * * * * * * * * * * //


void Foam::physicsManager::initBulletWorld(vector gValue)
{
    word frictionMode = physicsDict_.lookupOrDefault<word>("frictionMode", "");
    nIterBullet_ = physicsDict_.lookupOrDefault<int>("nIterBullet", 10);
    implicitEccentric_ = physicsDict_.lookupOrDefault<bool>("implicitEccentric", false);

    bulletWorld_.reset(new bulletWorld(gValue, frictionMode));
    
    Info << "Initialize Bullet World, friction: " << frictionMode << ", nunber of iterations: " << nIterBullet_ << ", implicitEccentric: " << implicitEccentric_ << endl;
        
    Info << nl << "Register objects: " << endl;
    
    const dictionary& objectsDict = physicsDict_.subDict("objects");
    
    int objectNum = 0;
        
    for (const entry& dEntry : objectsDict)
    {
        const dictionary& objectDict = dEntry.dict();
        
        objectDataList_.setSize(objectNum+1);
        objectDataList_.set(objectNum, new objectData);
        
        objectData& objectData = objectDataList_[objectNum];
        
        objectData.objName_ = objectDict.get<word>("objName");
        objectData.objShape_ = objectDict.get<word>("objShape");
        
        objectData.fluidForceRelaxation_ = objectDict.lookupOrDefault<scalar>("fluidForceRelaxation", 1);
        objectData.fluidForceDamping_ = objectDict.lookupOrDefault<scalar>("fluidForceDamping", 1);
        
        objectData.objDict_ = objectDict;
        
        Info << nl << "objName: " << objectData.objName_ << ", objShape: " << objectData.objShape_ << endl;
        
        bulletWorld_->addBulletBodies(objectDict);
        
        if(Pstream::master())
        {
            objectData.bulletBody_ = &bulletWorld_->bulletBodies()[objectNum];
        }
        
        Info << "   initialPosition: " <<  objectData.getPosition()  << endl;
        Info << "   initialRotation: " <<  objectData.getRotation()  << endl;
        
        objectNum++;
    }
    
    bulletWorld_->storeStates();
    
    bulletWorld_->writeVTK(outputCounter_);
}

void Foam::physicsManager::setJoints()
{
    if(physicsDict_.found("joints"))
    {
        const dictionary& jointsDict = physicsDict_.subDict("joints");
    
        int jointNum = 0;
        
        for (const entry& dEntry : jointsDict)
        {
            const dictionary& jointDict = dEntry.dict();
            
            const word jointType = jointDict.get<word>("type");
            
            if(jointType == "hingeSelf")
            {
                bulletWorld_->setHingeSelf(jointDict);
                Info << "   set hingeSelf joint" << endl;
            }
            else if(jointType == "hingePair")
            {
                bulletWorld_->setHingePair(jointDict);
                Info << "   set hingePair joint" << endl;
            }
            else if(jointType == "generic6DoF")
            {
                bulletWorld_->setGeneric6DoFJoint(jointDict);
                Info << "   set generic6DoF joint" << endl;
            }
            else
            {
                FatalErrorInFunction<< "Choose correct joint type: " << jointType << abort(FatalError);
            }
            
            jointNum++;
        }
    }
}

void Foam::physicsManager::initMoorDynWorld()
{
    if (!physicsDict_.found("MoorDyn"))
    {
        Info << "No 'MoorDyn' section in physicsDict. Skipping mooring setup." << endl;
        return;
    }

    const dictionary& MoorDynDict = physicsDict_.subDict("MoorDyn");

    if (!MoorDynDict.found("file") || !MoorDynDict.found("moorings"))
    {
        Info << "Empty or incomplete 'MoorDyn' dictionary found. Skipping mooring setup." << endl;
        return;
    }

    predictMooringForces_ = physicsDict_.lookupOrDefault<bool>("predictMooringForces", false);

    fileName dirPath("Mooring");
    if (!isDir(dirPath))
    {
        FatalErrorInFunction << "Directory 'Mooring' does not exist" << abort(FatalError);
    }

    const word file = MoorDynDict.get<word>("file");
    
    moorDynWorld_.reset(new moorDynWorld(file));
    
    const dictionary& mooringsDict = MoorDynDict.subDict("moorings");            
    for (const entry& dEntry : mooringsDict)
    {
        const dictionary& mooringDict = dEntry.dict();
        const word obj = mooringDict.get<word>("obj");
        
        int objI = getObjectID(obj);
        Info << "make mooring: " << obj << endl;
        
        vector P = objectDataList_[objI].getPosition();
        tensor Q = objectDataList_[objI].getRotation();
        tensor invQ = Q.T();
        
        const labelList pointIDs = mooringDict.lookup("pointIDs");
        
        objectDataList_[objI].mooringDataList_.setSize(pointIDs.size());
        
        forAll(pointIDs, i)
        {
            mooringData& mooringData = objectDataList_[objI].mooringDataList_[i];
        
            mooringData.mooringID_ = pointIDs[i];
            
            const vector mooringWorldPosition =  moorDynWorld_->getFairPosition(pointIDs[i]);                  
            
            mooringData.mooringLocalPosition_ = invQ & (mooringWorldPosition - P);
                
            Info << "   ID: " << pointIDs[i] << ", local position: " << mooringData.mooringLocalPosition_ << endl;
            
            moorDynWorld_->setPosVel(mooringData.mooringID_, mooringWorldPosition, vector::zero);
        }
        
    }

    moorDynWorld_->init();
    moorDynWorld_->writeVTK(outputCounter_);
    moorDynWorld_->save_mooring(runTime_.value());
}


// * * * * * * * * * * * * * * * * Update * * * * * * * * * * * * * * //


void Foam::physicsManager::updateFluidProperties(const scalar& deltaT)
{
    forAll(objectDataList_, objI)
    {
        objectData& objectData = objectDataList_[objI];
        
        List<innerMassData>& innerMassDataList = objectData.innerMassDataList_;
        
        if(innerMassDataList.size())
        {
            Info << nl << "Update fluid properties: " <<objectData.objName_ << endl;
            
            scalar totalFluidMass = 0;
            tensor setFluidInertia = Zero;
            vector setFluidCoM = Zero;
            
            forAll(innerMassDataList, i)
            {
                innerMassData& innerMassData = objectData.innerMassDataList_[i];
                fluidData& fluidData = innerMassData.fluidData_;
                
                if (!fluidData.fluidRho_)
                {
                    FatalErrorInFunction << "fluidRho_ is nullptr for object: "
                                         << objectData.objName_ << ", index: " << i << abort(FatalError);
                }
                if (!fluidData.V_)
                {
                    FatalErrorInFunction << "V_ is nullptr for object: "
                                         << objectData.objName_ << ", index: " << i << abort(FatalError);
                }
                
                calculateFluidProps
                (
                    *fluidData.fluidRho_,
                    innerMassData.base_.mesh_->V(),
                    fluidData.CLocal_,
                    fluidData.fluidMass_,
                    
                    fluidData.fluidInertia_,
                    fluidData.fluidCoMLocal_
                );
                
                totalFluidMass += fluidData.fluidMass_;
                setFluidInertia += fluidData.fluidInertia_;
                setFluidCoM += fluidData.fluidMass_*fluidData.fluidCoMLocal_;
            }
            
            setFluidCoM /= totalFluidMass;
            
            Info << nl << "Correct fluid properties: " <<objectData.objName_ << endl;
            Info << "   corrected fluidInertia: " << setFluidInertia << endl;
            Info << "   corrected fluidCoMLocal: " << setFluidCoM << endl;
            
            bulletWorld_->updateFluidPropertiesAndAngularVelocity(objI, setFluidInertia, setFluidCoM, deltaT);
        }
    }
}


void Foam::physicsManager::calculateFluidForces()
{
    forAll(objectDataList_, objI)
    {
        objectData& objectData = objectDataList_[objI];
        const vector CoMWorld = objectData.getCoMWorld();
        
        List<outerForceData>& outerForceDataList = objectData.outerForceDataList_;
        List<innerForceData>& innerForceDataList = objectData.innerForceDataList_;
    
        forAll(outerForceDataList, i)
        {
            forceData& forceData = outerForceDataList[i].forceData_;
            
            bulletMotionSolver& BMSolver = *outerForceDataList[i].base_.solver_;
                        
            BMSolver.calcForceMoment(CoMWorld);
            
            forceData.force_ = BMSolver.force();
            forceData.torque_ = BMSolver.moment();
            
            Info << "  outerForce: " << forceData.force_ << ", outerTorque: " << forceData.torque_ <<  nl;
        }
        
        forAll(innerForceDataList, i)
        {
            forceData& forceData = innerForceDataList[i].forceData_;
            
            bulletMotionSolver& BMSolver = *innerForceDataList[i].base_.solver_;
                        
            BMSolver.calcForceMoment(CoMWorld);
            
            forceData.force_ = BMSolver.force();
            forceData.torque_ = BMSolver.moment();
            
            Info << "  innerForce: " << forceData.force_ << ", innerTorque: " << forceData.torque_ <<  nl;
        }
    }
}

void Foam::physicsManager::calculateMooringForces()
{
    forAll(objectDataList_, objI)
    {
        objectData& objectData = objectDataList_[objI];
        const vector CoMWorld = objectData.getCoMWorld();
        
        List<mooringData>& mooringDataList = objectData.mooringDataList_;
        
        forAll(mooringDataList, i)
        { 
            mooringData& mooringData = mooringDataList[i];
            
            mooringData.mooringForce_ *= 0;
            mooringData.mooringTorque_ *= 0;
            
            moorDynWorld_->getForce(mooringData.mooringID_, mooringData.mooringForce_);
            
            reduce(mooringData.mooringForce_, sumOp<vector>());
            
            const vector P = objectData.getPosition();
            const tensor Q = objectData.getRotation();
            
            vector mooringWorldPosition = P + (Q & (mooringData.mooringLocalPosition_));
            vector torqueArm = mooringWorldPosition - CoMWorld;
            
            mooringData.mooringTorque_ = torqueArm^mooringData.mooringForce_;
            
            Info << "  mooringForce: " << mooringData.mooringForce_ << ", mooringTorque: " << mooringData.mooringTorque_  <<  nl;
        }
        
    }
}


void Foam::physicsManager::applyForces()
{
    forAll(objectDataList_, objI)
    {
        objectData& objectData = objectDataList_[objI];
                
        List<outerForceData>& outerForceDataList = objectData.outerForceDataList_;
        List<innerForceData>& innerForceDataList = objectData.innerForceDataList_;
        List<mooringData>& mooringDataList = objectData.mooringDataList_;
        
        if(outerForceDataList.size() || innerForceDataList.size() || mooringDataList.size())
        {        
            const scalar fRelax = objectData.fluidForceRelaxation_;
            const scalar fDamp = objectData.fluidForceDamping_;
            
            const vector oldTotalFluidForce = objectData.totalFluidForce_;
            const vector oldTotalFluidTorque = objectData.totalFluidTorque_;
            
            vector applyForce(Zero);
            vector applyTorque(Zero);
            
            vector totalFluidForce(Zero);
            vector totalFluidTorque(Zero);
            
            vector totalMooringForce(Zero);
            vector totalMooringTorque(Zero);
            
            //Fluid forces
            {
                forAll(outerForceDataList, i)
                {
                    forceData& forceData = outerForceDataList[i].forceData_;
                    
                    totalFluidForce += forceData.force_;
                    totalFluidTorque += forceData.torque_;
                }
                
                forAll(innerForceDataList, i)
                {
                    forceData& forceData = innerForceDataList[i].forceData_;
                    
                    totalFluidForce += forceData.force_;
                    totalFluidTorque +=forceData.torque_;
                }
                
                if(step_ > 1)
                {
                    totalFluidForce = fDamp*(fRelax*totalFluidForce + (1-fRelax)*oldTotalFluidForce);
                    totalFluidTorque = fDamp*(fRelax*totalFluidTorque + (1-fRelax)*oldTotalFluidTorque);
                }
                else
                {
                    totalFluidForce = fDamp*totalFluidForce;
                    totalFluidTorque = fDamp*totalFluidTorque;   
                }
                
                objectData.totalFluidForce_ = totalFluidForce;
                objectData.totalFluidTorque_ = totalFluidTorque;
            }
            
            //Mooring forces
            {        
                forAll(mooringDataList, i)
                {
                    totalMooringForce += mooringDataList[i].mooringForce_;
                    totalMooringTorque += mooringDataList[i].mooringTorque_;
                }

                objectData.totalMooringForce_ = totalMooringForce;
                objectData.totalMooringTorque_ = totalMooringTorque;
            }
            
            applyForce = totalFluidForce + totalMooringForce;
            applyTorque = totalFluidTorque + totalMooringTorque;

            Info << "  appliedForce: " << applyForce << ", appliedTorque: " << applyTorque  <<  nl;

            bulletWorld_->applyForceTorque(objI, applyForce, applyTorque);
        }
    }
}


void Foam::physicsManager::updateSolvers()
{
    forAll(objectDataList_, objI)
    {
        objectData& objectData = objectDataList_[objI];
        
        const vector P = objectData.getPosition();
        const tensor Q = objectData.getRotation();
                
        List<outerForceData>& outerForceDataList = objectData.outerForceDataList_;
        List<innerForceData>& innerForceDataList = objectData.innerForceDataList_;
        List<innerMassData>& innerMassDataList = objectData.innerMassDataList_;
        List<trackingData>& trackingDataList = objectData.trackingDataList_;
        
        if(outerForceDataList.size() || innerForceDataList.size() || innerMassDataList.size() || trackingDataList.size())
        { 
            forAll(outerForceDataList, i)
            {
                bulletMotionSolver& BMSolver = *outerForceDataList[i].base_.solver_;
                BMSolver.updateMotionState(P, Q);
            }
            
            forAll(innerForceDataList, i)
            {
                bulletMotionSolver& BMSolver = *innerForceDataList[i].base_.solver_;
                BMSolver.updateMotionState(P, Q);
            }
            
            forAll(innerMassDataList, i)
            {
                bulletMotionSolver& BMSolver = *innerMassDataList[i].base_.solver_;
                BMSolver.updateMotionState(P, Q);
            }
            
            forAll(trackingDataList, i)
            {
                bulletMotionSolver& BMSolver = *trackingDataList[i].base_.solver_;
                
                const vector& normal = trackingDataList[i].normal_;
                const vector& axis = trackingDataList[i].axis_;
                
                const vector& initialP = trackingDataList[i].initialP_;
                const tensor& initialQ = trackingDataList[i].initialQ_;
                
                vector setP = vector::zero;
                setP.x() = normal.x()*P.x() + (1-normal.x())*initialP.x();
                setP.y() = normal.y()*P.y() + (1-normal.y())*initialP.y();
                setP.z() = normal.z()*P.z() + (1-normal.z())*initialP.z();
                
                tensor R(I);
                
                if (axis.x() == 1.0)
                {
                    scalar theta = Foam::atan2(Q.zy(), Q.yy());
                    scalar cosT = Foam::cos(theta);
                    scalar sinT = Foam::sin(theta);

                    tensor Rx(
                        1,    0,     0,
                        0,  cosT, -sinT,
                        0,  sinT,  cosT
                    );

                    R = Rx & R;
                }

                if (axis.y() == 1.0)
                {
                    scalar theta = Foam::atan2(Q.xz(), Q.zz());
                    scalar cosT = Foam::cos(theta);
                    scalar sinT = Foam::sin(theta);

                    tensor Ry(
                        cosT, 0, sinT,
                           0, 1,    0,
                       -sinT, 0, cosT
                    );

                    R = Ry & R;
                }

                if (axis.z() == 1.0)
                {
                    scalar theta = Foam::atan2(Q.yx(), Q.xx());
                    scalar cosT = Foam::cos(theta);
                    scalar sinT = Foam::sin(theta);

                    tensor Rz(
                        cosT, -sinT, 0,
                        sinT,  cosT, 0,
                           0,     0, 1
                    );

                    R = Rz & R;
                }
                
                tensor setQ = R & initialQ;
                
                BMSolver.updateMotionState(setP, setQ);
            }
        }
    }       
}


void Foam::physicsManager::predictMoorDynWorld(const scalar& t, const scalar& deltaT)
{
    Info << nl << "predict mooring forces" << endl;
    
    forAll(objectDataList_, objI)
    {
        objectData& objectData = objectDataList_[objI];
        
        const scalar dt = runTime_.deltaTValue();
        
        const vector P = objectData.getPosition();
        const tensor Q = objectData.getRotation();
        
        const vector v = objectData.getVelocity();
        const vector pi = objectData.getAngularVelocity();
        
        const vector CoMLocal = objectData.getCoMLocal();
        
        vector piDt = pi * deltaT;
        
        tensor skewPiDt = tensor(
             0,     -piDt.z(),  piDt.y(),
             piDt.z(),  0,     -piDt.x(),
            -piDt.y(), piDt.x(),   0
        );
        const tensor Qp = Q + (skewPiDt & Q);
        
        const vector Pp = P + v * dt - ((Qp - Q) & CoMLocal);
                
        List<mooringData>& mooringDataList = objectData.mooringDataList_;
        
        if(mooringDataList.size())
        {
            Info << nl << "    objName: " << objectData.objName_ << endl;
        
            Info << "    velocity: " << v << endl;
            Info << "    predictedPosition: " <<  Pp  << endl;
            Info << "    predictedRotation: " <<  Qp  << endl;
        }
        
        forAll(mooringDataList, i)
        { 
            mooringData& mooringData = mooringDataList[i];
            const vector& mPosLocal = mooringData.mooringLocalPosition_;
            const int& mID = mooringData.mooringID_;
            
            vector mPos =  P + (Q & mPosLocal);
            vector mPosp =  Pp + (Qp & mPosLocal);
            
            vector mVelp = (mPosp - mPos)/deltaT;
            
            moorDynWorld_->setPosVel(mID, mPosp, mVelp);
        }
    }
    
    moorDynWorld_->step(t, deltaT);
}


void Foam::physicsManager::updateMoorDynWorld(const scalar& t, const scalar& deltaT)
{
    
    forAll(objectDataList_, objI)
    {
        objectData& objectData = objectDataList_[objI];
        
        const vector P0 = objectData.getPosition0();
        const tensor Q0 = objectData.getRotation0();
        const vector P = objectData.getPosition();
        const tensor Q = objectData.getRotation();
        
        List<mooringData>& mooringDataList = objectData.mooringDataList_;
        
        forAll(mooringDataList, i)
        { 
            mooringData& mooringData = mooringDataList[i];
            const vector& mPosLocal = mooringData.mooringLocalPosition_;
            const int& mID = mooringData.mooringID_;
            
            vector mPos =  P + (Q & mPosLocal);
            vector mPos0 =  P0 + (Q0 & mPosLocal);
            
            vector mVel = (mPos - mPos0)/deltaT;
            
            moorDynWorld_->setPosVel(mID, mPos, mVel);
        }
    }
    
    moorDynWorld_->step(t, deltaT);
}


// * * * * * * * * * * * * * * * * Calculater * * * * * * * * * * * * * * //

void Foam::physicsManager::calculateMassProps
(
    const scalarField& rho,
    const scalarField& V,
    const vectorField& CLocal,
    
    scalar& mass,
    tensor& inertia,
    vector& CoMLocal
)
{
    mass = 0;
    inertia = Zero;
    CoMLocal = Zero; 
    
    forAll(CLocal, celli)
    {
        scalar massElement = rho[celli]*V[celli];
        
        scalar x = CLocal[celli][0];
        scalar y = CLocal[celli][1];
        scalar z = CLocal[celli][2];  
        
        mass += massElement;
        CoMLocal += massElement*CLocal[celli];
        
        inertia[0] += massElement*(sqr(y) + sqr(z));
        inertia[1] += -massElement*y*x;
        inertia[2] += -massElement*z*x;
        inertia[3] += -massElement*x*y;
        inertia[4] += massElement*(sqr(x) + sqr(z));
        inertia[5] += -massElement*z*y;
        inertia[6] += -massElement*x*z;
        inertia[7] += -massElement*y*z;
        inertia[8] += massElement*(sqr(x) + sqr(y));
    }
    reduce(mass, sumOp<scalar>());
    reduce(inertia, sumOp<tensor>());
    reduce(CoMLocal, sumOp<vector>());
    
    CoMLocal /= mass;
}


void Foam::physicsManager::calculateFluidProps
(
    const scalarField& rho,
    const scalarField& V,
    const vectorField& CLocal,
    const scalar& mass,
    
    tensor& inertia,
    vector& CoMLocal
)
{
    inertia = Zero;
    CoMLocal = Zero; 
    
    forAll(CLocal, celli)
    {
        scalar massElement = rho[celli]*V[celli];
        
        scalar x = CLocal[celli][0];
        scalar y = CLocal[celli][1];
        scalar z = CLocal[celli][2];  
        
        CoMLocal += massElement*CLocal[celli];
        
        inertia[0] += massElement*(sqr(y) + sqr(z));
        inertia[1] += -massElement*y*x;
        inertia[2] += -massElement*z*x;
        inertia[3] += -massElement*x*y;
        inertia[4] += massElement*(sqr(x) + sqr(z));
        inertia[5] += -massElement*z*y;
        inertia[6] += -massElement*x*z;
        inertia[7] += -massElement*y*z;
        inertia[8] += massElement*(sqr(x) + sqr(y));
    }
    reduce(inertia, sumOp<tensor>());
    reduce(CoMLocal, sumOp<vector>());
    
    CoMLocal /= mass;
}


// * * * * * * * * * * * * * * * * Access * * * * * * * * * * * * * * //


int Foam::physicsManager::getObjectID(word objName)
{
    forAll(objectDataList_, i)
    {
        if(objectDataList_[i].objName_ == objName)
        {
            return i;
        }
    }
    
    FatalErrorInFunction<< "Choose correct object name: " << objName << abort(FatalError);
    return -1;
}


// * * * * * * * * * * * * * * * * Write Data * * * * * * * * * * * * * * //


void Foam::physicsManager::prepareCSV()
{
    forAll(objectDataList_, objI)
    {
        objectData& objectData = objectDataList_[objI];
        
        const word objName = objectData.objName_;
    
        objectData.stateCSV_.reset(new OFstream("Bullet/" + objName + "State.csv"));
        OFstream& stateCSV = *objectData.stateCSV_;
        
        stateCSV << "Time[s],posX[m],posY[m],posZ[m],CoMX[m],CoMY[m],CoMZ[m],rolX[deg],rolY[deg],rolZ[deg],velX[m/s],velY[m/s],velZ[m/s],omegaX[rad/s],omegaY[rad/s],omegaZ[rad/s],Ixx[kg*m2],Ixy[kg*m2],Ixz[kg*m2],Iyy[kg*m2],Iyz[kg*m2],Izz[kg*m2]" << endl;
    
        List<outerForceData>& outerForceDataList = objectData.outerForceDataList_;
        List<innerForceData>& innerForceDataList = objectData.innerForceDataList_;
        List<mooringData>& mooringDataList = objectData.mooringDataList_;
    
        if(outerForceDataList.size() || innerForceDataList.size() || mooringDataList.size())
        {     
            objectData.forceCSV_.reset(new OFstream("Bullet/" + objName + "Force.csv"));
            OFstream& forceCSV = *objectData.forceCSV_;
    
            forceCSV << "Time[s],";
            
            forAll(outerForceDataList, i)
            {
                forceCSV << "outer" << i << "ForceX[N],outer" << i << "ForceY[N],outer" << i << "ForceZ[N],outer" << i << "TorqueX[N*m],outer" << i << "TorqueY[N*m],outer" << i << "TorqueZ[N*m],";
            }
            
            forAll(innerForceDataList, i)
            {
                forceCSV << "inner" << i << "ForceX[N],inner" << i << "ForceY[N],inner" << i << "ForceZ[N],inner" << i << "TorqueX[N*m],inner" << i << "TorqueY[N*m],inner" << i << "TorqueZ[N*m],";
            }
    
            forAll(mooringDataList, i)
            {
                forceCSV << "mooring" << i << "ForceX[N],mooring" << i << "ForceY[N],mooring" << i << "ForceZ[N],mooring" << i << "TorqueX[N*m],mooring" << i << "TorqueY[N*m],mooring" << i << "TorqueZ[N*m],";
            }
    
            forceCSV << "totalForceX[N],totalForceY[N],totalForceZ[N],totalTorqueX[N*m],totalTorqueY[N*m],totalTorqueZ[N*m]" << endl;
        }
    }
}


void Foam::physicsManager::writeCSV()
{
    forAll(objectDataList_, objI)
    {
        objectData& objectData = objectDataList_[objI];

        OFstream& stateCSV = *objectData.stateCSV_;
        
        const scalar time = runTime_.value();
        const vector pos = objectData.getPosition();
        const vector CoM = objectData.getCoMWorld();
        const vector rol = objectData.getRotationVector();  
        
        const vector vel = objectData.getVelocity();
        const vector omega = objectData.getAngularVelocity();

        const tensor I = objectData.getTotalInertia();

        stateCSV << time << "," << pos[0] << "," << pos[1] << "," << pos[2] << "," << CoM[0] << "," << CoM[1] << "," << CoM[2] << "," << rol[0] << "," << rol[1] << "," << rol[2] << ","
            << vel[0] << "," << vel[1] << "," << vel[2] << "," << omega[0] << "," << omega[1] << "," << omega[2] << ","
            << I[0] << "," << I[1] << "," << I[2] << "," << I[4] << "," << I[5] << "," << I[8] << endl;

        List<outerForceData>& outerForceDataList = objectData.outerForceDataList_;
        List<innerForceData>& innerForceDataList = objectData.innerForceDataList_;
        List<mooringData>& mooringDataList = objectData.mooringDataList_;

        if(outerForceDataList.size() || innerForceDataList.size() || mooringDataList.size())
        {    
            OFstream& forceCSV = *objectData.forceCSV_;

            forceCSV << runTime_.value() << ",";
            
            forAll(outerForceDataList, i)
            {
                forceData& forceData = outerForceDataList[i].forceData_;
                const vector outerForce = forceData.force_;
                const vector outerTorque = forceData.torque_;

                forceCSV << outerForce[0] << "," << outerForce[1] << "," << outerForce[2] << "," << outerTorque[0] << "," << outerTorque[1] << "," << outerTorque[2] << ","; 
            }
            
            forAll(innerForceDataList, i)
            {
                forceData& forceData = innerForceDataList[i].forceData_;
                const vector innerForce = forceData.force_;
                const vector innerTorque = forceData.torque_;

                forceCSV << innerForce[0] << "," << innerForce[1] << "," << innerForce[2] << "," << innerTorque[0] << "," << innerTorque[1] << "," << innerTorque[2] << ","; 
            }

            forAll(mooringDataList, i)
            {
                const vector mooringForce =  mooringDataList[i].mooringForce_;
                const vector mooringTorque =  mooringDataList[i].mooringTorque_;

                forceCSV << mooringForce[0] << "," << mooringForce[1] << "," << mooringForce[2] << "," << mooringTorque[0] << "," << mooringTorque[1] << "," << mooringTorque[2] << ","; 
            }

            const vector totalForce = objectData.totalFluidForce_ + objectData.totalMooringForce_;
            const vector totalTorque = objectData.totalFluidTorque_ + objectData.totalMooringTorque_;
            
            forceCSV << totalForce[0] << "," << totalForce[1] << "," << totalForce[2] << "," << totalTorque[0] << "," << totalTorque[1] << "," << totalTorque[2] << endl; 
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::physicsManager::physicsManager(const Time& runTime, vector gValue)
:
    runTime_(runTime)
{    
    physicsDict_ = IOdictionary
    (
        IOobject
        (
            "physicsDict",
            runTime.time().constant(),
            runTime.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    
    nOuterCorrectors_ = physicsDict_.lookupOrDefault<int>("nOuterCorrectors", 1);
    
    initBulletWorld(gValue);
    
    setJoints();
    
    initMoorDynWorld();
}



// * * * * * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * * * //

void Foam::physicsManager::update()
{    
    updateFluidProperties(runTime_.deltaTValue());
    
    calculateFluidForces();
    
    if(moorDynWorld_)
    {
        if(step_ > 1)
        {
            predictMoorDynWorld(runTime_.value(), runTime_.deltaTValue());
        }
        
        calculateMooringForces();
    }
    
    applyForces();
    
    Info << nl << "deltaT: "  << runTime_.deltaTValue()  << ", deltaT in Bullet: " << runTime_.deltaTValue()/nIterBullet_ << endl;
    
    bulletWorld_->stepSimulation(runTime_.deltaTValue(), nIterBullet_, runTime_.deltaTValue()/nIterBullet_);
    
    forAll(objectDataList_, objI)
    {
        objectData& objectData = objectDataList_[objI];
        
        Info << "object: " << objectData.objName_ << endl;
        Info << "   position before step: " << objectData.getPosition0() << endl;
        Info << "   position after step: " << objectData.getPosition() << endl;
        Info << "   linear velocity: " << objectData.getVelocity() << endl;
        Info << "   angular velocity: " << objectData.getAngularVelocity() << endl;
    }
        
    updateSolvers();

    step_++;
}


void Foam::physicsManager::restoreStates()
{
    bulletWorld_->restoreStates();
    if(moorDynWorld_)
    {  
        moorDynWorld_->load_mooring();
    }
}


void Foam::physicsManager::writeAndStoreStates() 
{
    //store states
    bulletWorld_->storeStates();
    
    if(moorDynWorld_)
    {
        moorDynWorld_->load_mooring();
        updateMoorDynWorld(runTime_.value(), runTime_.deltaTValue());
        moorDynWorld_->save_mooring(runTime_.value());
    }
    
    writeCSV();
    
    if(runTime_.outputTime())
    {
        outputCounter_++;
        
        bulletWorld_->writeVTK(outputCounter_);
                
        if(moorDynWorld_)
        {
            moorDynWorld_->writeVTK(outputCounter_);
        }
    }
}


void Foam::physicsManager::registerOuterForce
(
    const word& objName,
    bulletOversetFvMesh& mesh,
    bulletMotionSolver& solver
)
{
    int objI = getObjectID(objName);
    
    objectData& objectData = objectDataList_[objI];
    
    int registerI = objectData.outerForceDataList_.size();
    objectData.outerForceDataList_.setSize(registerI+1);
    
    // initialize motion state
    const vector P = objectDataList_[objI].getPosition();
    const tensor Q = objectDataList_[objI].getRotation();
    solver.initializeMotionState(P, Q);
    
    // register base data
    outerForceData& data = objectData.outerForceDataList_[registerI];
    
    data.base_.mesh_ = &mesh;
    data.base_.solver_ = &solver;
}


void Foam::physicsManager::registerInnerForce
(
    const word& objName,
    bulletOversetFvMesh& mesh,
    bulletMotionSolver& solver
)
{
    int objI = getObjectID(objName);
    
    objectData& objectData = objectDataList_[objI];
    
    int registerI = objectData.innerForceDataList_.size();
    objectData.innerForceDataList_.setSize(registerI+1);
    
    // initialize motion state
    const vector P = objectDataList_[objI].getPosition();
    const tensor Q = objectDataList_[objI].getRotation();
    solver.initializeMotionState(P, Q);
    
    // register base data
    innerForceData& data = objectData.innerForceDataList_[registerI];
    data.base_.mesh_ = &mesh;
    data.base_.solver_ = &solver;
    
    //register inner solid data
    const tensor invQ = Q.T();
    const scalar solidRho = objectDataList_[objI].getSolidRho();
    scalarField rho(mesh.nCells(), solidRho);
    vectorField CLocal = invQ & (mesh.C().primitiveField() - P);
    
    if(solidRho > 0)
    {
        calculateMassProps(rho, mesh.V(), CLocal, data.solidData_.solidMass_, data.solidData_.solidInertia_, data.solidData_.solidCoMLocal_);
    }
    
    Info << "   object: " << objectData.objName_ << endl;
    Info << "   registered solidMass: " << data.solidData_.solidMass_ << endl;
    Info << "   registered solidInertia: " << data.solidData_.solidMass_ << endl;
    Info << "   registered solidCoMLocal: " << data.solidData_.solidMass_ << endl;
}


void Foam::physicsManager::registerInnerMass
(
    const word& objName,
    bulletOversetFvMesh& mesh,
    bulletMotionSolver& solver,
    scalarField& fluidRho
)
{    
    int objI = getObjectID(objName);
    
    objectData& objectData = objectDataList_[objI];
    
    int registerI = objectData.innerMassDataList_.size();
    objectData.innerMassDataList_.setSize(registerI+1);
    
    // initialize motion state
    const vector P = objectDataList_[objI].getPosition();
    const tensor Q = objectDataList_[objI].getRotation();
    solver.initializeMotionState(P, Q);
    
    // register base data
    innerMassData& data = objectData.innerMassDataList_[registerI];
    data.base_.mesh_ = &mesh;
    data.base_.solver_ = &solver;
    
    //register inner solid data
    const tensor invQ = Q.T();
    const scalar solidRho = objectDataList_[objI].getSolidRho();
    scalarField rho(mesh.nCells(), solidRho);
    vectorField CLocal = invQ & (mesh.C().primitiveField() - P);
    
    if(solidRho > 0)
    {
        calculateMassProps(rho, mesh.V(), CLocal, data.solidData_.solidMass_, data.solidData_.solidInertia_, data.solidData_.solidCoMLocal_);
    }
    
    Info << "   object: " << objectData.objName_ << endl;
    Info << "   registered solidMass: " << data.solidData_.solidMass_ << endl;
    Info << "   registered solidInertia: " << data.solidData_.solidInertia_ << endl;
    Info << "   registered solidCoMLocal: " << data.solidData_.solidCoMLocal_ << endl;
    
    //register inner solid data
    data.fluidData_.CLocal_ = CLocal;
    data.fluidData_.V_ = &mesh.V();
    data.fluidData_.fluidRho_ = &fluidRho;
    
    calculateMassProps(fluidRho, mesh.V(), CLocal, data.fluidData_.fluidMass_, data.fluidData_.fluidInertia_, data.fluidData_.fluidCoMLocal_);
    Info << "   registered fluidMass: " << data.fluidData_.fluidMass_ << endl;
    Info << "   registered fluidInertia: " << data.fluidData_.fluidInertia_ << endl;
    Info << "   registered fluidCoMLocal: " << data.fluidData_.fluidCoMLocal_ << endl;
}


void Foam::physicsManager::registerTracking
(
    const word& objName,
    bulletOversetFvMesh& mesh,
    bulletMotionSolver& solver,
    vector normal,
    vector axis
)
{
    int objI = getObjectID(objName);
    
    objectData& objectData = objectDataList_[objI];
    
    int registerI = objectData.trackingDataList_.size();
    objectData.trackingDataList_.setSize(registerI+1);
    
    // initialize motion state
    const vector P = objectDataList_[objI].getPosition();
    const tensor Q = objectDataList_[objI].getRotation();
    solver.initializeMotionState(P, Q);
    
    // register base data
    trackingData& data = objectData.trackingDataList_[registerI];
    data.base_.mesh_ = &mesh;
    data.base_.solver_ = &solver;
    
    // register tracking data
    data.normal_ = normal;
    data.axis_ = axis;
    
    data.initialP_ = P;
    data.initialQ_ = Q;
    
    Info << "   object: " << objectData.objName_ << endl;
    Info << "   registered tracking, normal: " << data.normal_ << ", axis: " << data.axis_ << endl;
}



void Foam::physicsManager::reregisterObjects()
{
    forAll(objectDataList_, objI)
    {
        objectData& objectData = objectDataList_[objI];
        const scalar solidRho = objectData.getSolidRho();
        
        List<innerForceData>& innerForceDataList = objectData.innerForceDataList_;
        List<innerMassData>& innerMassDataList = objectData.innerMassDataList_;
        
        if(innerForceDataList.size() || innerMassDataList.size())
        {
            Info << nl << "Update properties in Bullet: " <<objectData.objName_ << endl;
        }
        
        if(solidRho > 0)
        {
            if(innerForceDataList.size() || innerMassDataList.size())
            {
                scalar totalInnerMass = 0;
                tensor totalInnerInertia = Zero;
                vector totalInnerCoM = Zero;
                
                forAll(innerForceDataList, i)
                {
                    solidData& solidData = innerForceDataList[i].solidData_;
                    
                    totalInnerMass += solidData.solidMass_;
                    totalInnerInertia += solidData.solidInertia_;
                    totalInnerCoM += solidData.solidMass_*solidData.solidCoMLocal_;
                }
                
                forAll(innerMassDataList, i)
                {
                    solidData& solidData = innerMassDataList[i].solidData_;
                    
                    totalInnerMass += solidData.solidMass_;
                    totalInnerInertia += solidData.solidInertia_;
                    totalInnerCoM += solidData.solidMass_*solidData.solidCoMLocal_;
                }
                
                totalInnerCoM /= totalInnerMass;
                
                const scalar solidMass = objectData.getSolidMass();
                const tensor solidInertia = objectData.getSolidInertia();
                
                scalar setSolidMass = solidMass - totalInnerMass;
                tensor setSolidInertia = solidInertia - totalInnerInertia;
                vector setSolidCoM = -totalInnerMass*totalInnerCoM/setSolidMass;
                
                bulletWorld_->setSolidProperties(objI, setSolidMass, setSolidInertia, setSolidCoM);
                
                Info << nl << "Correct solid properties: " <<objectData.objName_ << endl;
                Info << "   corrected solidMass: " << setSolidMass << endl;
                Info << "   corrected solidInertia: " << setSolidInertia << endl;
                Info << "   corrected solidCoMLocal: " << setSolidCoM << endl;
            }
        }
        
        if(innerMassDataList.size())
        {
            scalar totalFluidMass = 0;
            tensor totalFluidInertia = Zero;
            vector totalFluidCoM = Zero;
                            
            forAll(innerMassDataList, i)
            {
                fluidData& fluidData = innerMassDataList[i].fluidData_;
                
                totalFluidMass += fluidData.fluidMass_;
                totalFluidInertia += fluidData.fluidInertia_;
                totalFluidCoM += fluidData.fluidMass_*fluidData.fluidCoMLocal_;
            }
            
            totalFluidCoM /= totalFluidMass;
            
            bulletWorld_->setFluidProperties(objI, totalFluidMass, totalFluidInertia, totalFluidCoM);
            
            Info << nl << "Add fluid properties: " <<objectData.objName_ << endl;
            Info << "   added fluidMass: " << totalFluidMass << endl;
            Info << "   added fluidInertia: " << totalFluidInertia << endl;
            Info << "   added fluidCoMLocal: " << totalFluidCoM << endl;
        }
        
        if(innerForceDataList.size() || innerMassDataList.size())
        {
            bulletWorld_->reinitializeMassProperties(objI, implicitEccentric_);
        }
    }
    
    bulletWorld_->storeStates();

    prepareCSV();
}



// ************************************************************************* //
