#include "moorDynWorld.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::moorDynWorld::moorDynWorld(const word fileName)
{
    if(Pstream::master())
    {        
        std::string dirName = "Mooring/" + fileName;
        moordyn_ = MoorDyn_Create(dirName.c_str());
       
        if (!moordyn_) 
        {
            FatalError << "MoorDyn v2 cannot be created!" << exit(FatalError);
        }
        
        unsigned int n = 0;
        MoorDyn_NCoupledDOF(moordyn_, &n);
        nDOF_ = n;
        
        Info << "nDOF: " << nDOF_ << endl;
        
        fairPos_ = vectorField(int(nDOF_/3), vector::zero);
        fairVel_ = vectorField(int(nDOF_/3), vector::zero);
        fairForce_ = vectorField(int(nDOF_/3), vector::zero);
        
        rmDir("Mooring/VTK");
        mkDir("Mooring/VTK");
        
        moordyn_backup_.t = 0.0;
        moordyn_backup_.data = nullptr;
    }
}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::moorDynWorld::init()
{
    if(Pstream::master())
    {
        MoorDyn_Init(moordyn_, &fairPos_[0][0], &fairVel_[0][0]);
    }
}


void Foam::moorDynWorld::step(const scalar& t, const scalar& dt)
{
    if(Pstream::master())
    {
        scalar tDummy = t;
        scalar dtDummy = dt;
        
        MoorDyn_Step(moordyn_, &fairPos_[0][0], &fairVel_[0][0], &fairForce_[0][0], &tDummy, &dtDummy);
    }
}


void Foam::moorDynWorld::setPosVel(const int& i, const vector& pos, const vector& vel)
{
    if(Pstream::master())
    {
        const int j = i-1;
        if(3*i > nDOF_)
        {
            FatalError << "Out of range!" << exit(FatalError);
        }
                
        fairPos_[j] = pos;
        fairVel_[j] = vel;
    }
}


vector Foam::moorDynWorld::getFairPosition(const int& i)
{
    vector fairPosition(Zero);
    
    if(Pstream::master())
    {
        MoorDynPoint pt = MoorDyn_GetPoint(moordyn_, i);
        double pos_d[3];
        
        MoorDyn_GetPointPos(pt, pos_d);
        
        fairPosition = vector(pos_d[0], pos_d[1], pos_d[2]);
    }
    reduce(fairPosition, sumOp<vector>());
    
    return fairPosition;
}


void Foam::moorDynWorld::getForce(const int& i, vector& force)
{
    if(Pstream::master())
    {
        const int j = i-1;
        if(3*i > nDOF_)
        {
            FatalError << "Out of range!" << exit(FatalError);
        }
        
        force = fairForce_[j];
    }
}


void Foam::moorDynWorld::writeVTK(const int& outputCounter)
{
    if(Pstream::master())
    {
        OFstream mps("Mooring/VTK/line_" +  std::__cxx11::to_string(outputCounter) + ".vtk");
        mps.precision(4);

        unsigned int nLines, nSeg;
        MoorDyn_GetNumberLines(moordyn_, &nLines);

        labelList nodesPerLine(nLines, -1);
        for(int i=0; i<int(nLines); i++)
        {
            MoorDynLine line = MoorDyn_GetLine(moordyn_, i+1);
            MoorDyn_GetLineN(line, &nSeg);
            nodesPerLine[i] = nSeg+1;
        }

        double coord[max(nodesPerLine)][3];

        // Writing header
        mps << "# vtk DataFile Version 3.0" << nl
            //<< "MoorDyn v2 vtk output time=" << runTime.timeName()
            << nl << "ASCII" << nl << "DATASET POLYDATA" << endl;
     
        // Writing points
        mps << "\nPOINTS " << sum(nodesPerLine) << " float" << endl;

        for(int i=0; i<int(nLines); i++)
        {   
            //map_.getNodeCoordinates(i, nodesPerLine_[i], &coord[0][0]);
            MoorDynLine line = MoorDyn_GetLine(moordyn_, i+1);
            
            for(int p=0; p<nodesPerLine[i]; p++)
            {
                MoorDyn_GetLineNodePos(line, p, &coord[p][0]);
                
                mps << coord[p][0] << " " << coord[p][1] << " " << coord[p][2] << endl;
            }
        }
        
        // Writing lines
        mps << "\nLINES " << nLines << " " << sum(nodesPerLine+1) << endl;

        label start_node(0);
        for(int i=0; i<int(nLines); i++)
        {       
            mps << nodesPerLine[i];
            
            for(int j=0; j<nodesPerLine[i]; j++)
            {
                mps << " " << start_node+j;
            }
            mps << endl;
            
            start_node += nodesPerLine[i];
        }

        // Writing tension (POINT_DATA)
        mps << "\nPOINT_DATA " << sum(nodesPerLine) << "\n";
        mps << "SCALARS tension float 1\n";
        mps << "LOOKUP_TABLE default\n";

        for (int i = 0; i < int(nLines); i++)
        {
            MoorDynLine line = MoorDyn_GetLine(moordyn_, i+1);

            for (int p = 0; p < nodesPerLine[i]; p++)
            {
                double tension = 0.0;
                MoorDyn_GetLineNodeTen(line, p, &tension);
                mps << tension << "\n";
            }
        }
    }
}



// ************************************************************************* //
