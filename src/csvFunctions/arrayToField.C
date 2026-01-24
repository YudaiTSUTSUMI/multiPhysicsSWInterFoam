#include "arrayToField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::arrayToField::arrayToField
(
    volScalarField& field,
    const List<scalarList>& array,
    const scalar& Xmin,
    const scalar& Xmax,
    const scalar& Ymin,
    const scalar& Ymax,
    const scalar& delta
):
    field_(field),
    array_(array),
    Xmin_(Xmin),
    Xmax_(Xmax),
    Ymin_(Ymin),
    Ymax_(Ymax),
    delta_(delta)
{
    const label nX = array_[0].size();
    const label nY = array_.size();
    
    Info << "nX = " << nX << ", nY = " << nY << endl;
    
    scalarList arrayX(nX,0);
    for(int i = 0; i < nX; i++)
    {
        if(i == 0)
        {
            arrayX[i] = Xmin_ + delta_/2;
        }
        else
        {
            arrayX[i] = arrayX[i-1] + delta;
        }
    }
    
    scalarList arrayY(nY,0);
    for(int i = 0; i < nY; i++)
    {
        if(i == 0)
        {
            arrayY[i] = Ymax_ - delta_/2;
        }
        else
        {
            arrayY[i] = arrayY[i-1] - delta;
        }
    }
    
    forAll(field,celli)
    {
        const scalar fieldX =  field_.mesh().C()[celli][0];
        const scalar fieldY =  field_.mesh().C()[celli][1];

        int iX = 0;
        int iY = 0;
	      
        if(Xmin <= fieldX && fieldX <= Xmax && Ymin <= fieldY && fieldY <= Ymax)
        {
            scalar Xdistance = mag(arrayX[0]-fieldX);
            scalar Ydistance = mag(arrayY[0]-fieldY);
	      
            for(int i = 0; i < nX; i++)
            {
                scalar Xdistance2 = mag(arrayX[i] - fieldX);
                if(Xdistance2 <= Xdistance)
                {
                    Xdistance = Xdistance2;
		     iX = i;
		 }
            }
	     
	    for(int i = 0; i < nY; i++)
            {
                scalar Ydistance2 = mag(arrayY[i] - fieldY);
                if(Ydistance2 <= Ydistance)
                {
                    Ydistance = Ydistance2;
		     iY = i;
		 }
            }
            
        field[celli] = array_[iY][iX];
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// ************************************************************************* //
