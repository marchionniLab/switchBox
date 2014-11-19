#include <R.h>
#include "math.h"

#define max(x,y)  ((x)<(y))? (y):(x)
// comp = vect2comp(vect);
// comp is comp(i,j)=\delta{vect(i)<vect(j)}
//
// if vect is a matrix, comp will be a three dimensional array
//

 extern "C"
{
   
void vect2compC( double *vect, int *np,int *mp, double *comp)
{
    
    //I/O variables
    //n number of elements in vect and m the number of vectors
    //vect to be processed
    // comp    should be of size dimN^2*dimM
    int n = *np;
    int m = *mp;
    
    if( n == 1)
    {
        n = m;
        m = 1;
    }
    
    
    //Internal variables
    int i,j,k,n2;//indices
    n2 = n*n;

    //Loading input    
    
    
    
    for(k=0;k<m;k++)
    {
        for(i=0;i<n;i++)
            for(j=i+1;j<n;j++)
            {
                if( vect[i+k*n] < vect[j+k*n])
                {
                    comp[i+j*n+k*n2] = 1;
                    comp[j+i*n+k*n2] = 0;
                }
                else
                {
                    comp[i+j*n+k*n2] = 0;
                    comp[j+i*n+k*n2] = 1;
                }

        }   
    }


    
    
}


}
