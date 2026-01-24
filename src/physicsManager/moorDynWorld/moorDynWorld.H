#ifndef moorDynWorld_H
#define moorDynWorld_H

#include "fvCFD.H"

#include "MoorDyn2.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class moorDynWorld
{
    // Private data    
    MoorDyn moordyn_;
    
    int nDOF_ = 0;
    
    vectorField fairPos_;
    
    vectorField fairVel_;
    
    vectorField fairForce_;
    
    
    // A backup of the Mooring system to can restore in between implicit time
    // steps
    typedef struct _moordyn_backup {
        double t;
        uint64_t *data;
    } moordyn_backup;
    mutable moordyn_backup moordyn_backup_;
    
            
public:

    // Constructors
    moorDynWorld(const word fileName);

    //- Destructor
    ~moorDynWorld()
    {}
    
    // Menber functions
    void init();
    
    void step(const scalar& t, const scalar& dt);
    
    void setPosVel(const int& i, const vector& pos, const vector& vel);
    
    vector getFairPosition(const int& i);
    
    void getForce(const int& i, vector& force);

    void writeVTK(const int& outputCounter);
    
    //- Create a backup of the mooring system
    inline void save_mooring(const double t) const
    {
        if(Pstream::master())
        {    
            int err;
            size_t data_size;
            err = MoorDyn_Serialize(moordyn_, &data_size, NULL);
            if (err != MOORDYN_SUCCESS) {
                FatalError << "MoorDyn system size cannot be retrieved"
                           << exit(FatalError);
            }
            if (moordyn_backup_.data)
                free(moordyn_backup_.data);
            //moordyn_backup_.data = (uint64_t*)malloc(data_size);
            moordyn_backup_.data = static_cast<uint64_t*>(malloc(data_size));
            if (!moordyn_backup_.data) {
                FatalError << "Failure allocating " << data_size << " bytes"
                           << exit(FatalError);
            }
            err = MoorDyn_Serialize(moordyn_, NULL, moordyn_backup_.data);
            if (err != MOORDYN_SUCCESS) {
                FatalError << "MoorDyn system cannot be saved"
                           << exit(FatalError);
            }
            moordyn_backup_.t = t;
        }
    }

    //- Restore the backup of the mooring system
    inline void load_mooring() const
    {
        if(Pstream::master())
        {   
            if (!moordyn_backup_.data) {
                FatalError << "No available backup"
                           << exit(FatalError);
            }
            const int err = MoorDyn_Deserialize(moordyn_, moordyn_backup_.data);
            if (err != MOORDYN_SUCCESS) {
                FatalError << "MoorDyn system cannot be restored"
                           << exit(FatalError);
            }
        }
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

