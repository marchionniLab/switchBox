#include <R.h>
#include "math.h"


#define max(x,y)  ((x)<(y))? (y):(x)
// comp = vect2comp(vect);
// comp is comp(i,j)=\delta{vect(i)<vect(j)}
extern "C"
{
    
    
    
/*// [score,a,M,k,N, P, Q] = CalculateSignedScoreCore(situation,data(indexup,:),data(indexdown,:));*/
    // score is the pairwise score between
    
    
void CalculateSignedScoreCore(int *situation, int *rowLenPtr, double *data1,int *colLen1Ptr, double *data2, int *colLen2Ptr,//inputs
    double *output1, double *output3, double *output4, double *output5, double *output6, double *output7,double *output8)//outputs
{
    //void CalculateSignedScoreCore(int *situation, int rowLen, double *data1,int colNo1, double *data2, int colNo2,//inputs
    //double *score, double *a, double *M, double *k, double *N, double *P,double *Q)//outputs
    //Calculates pairwise score between every gene in data1 and every gene in data2. It is signed because
    // score(i,j)=P(X_i<X_j|1)-P(X_i<X_j|0) - P(X_i>X_j|1) +P(X_i>X_j|0) + C
    // The third and the forth terms are for avoiding the equlities. The C term is proportion to the secondary score
    // to avoid the ties.  C = (E(X_j-X_i|1)-E(X_j-X_i|0))/K where K is big enough to make sure that the secondary score
    // does not intervene with the primary score.
    // a,M,k,N are usable for fisher exact test. M is the size of the class situation == 0 and N is the total number
    // of samples. a_ij = #(X_i<X_j|0) k_ij = #(X_i<X_j)
    // P_ij = P(X_i<X_j|0) and Q(X_i<X_j|1)
    
    // Note that the colLen1, colLen2 and rowLen represent row and column length respectively. It is due to the early changes in the program that 
    // they maintain their cabapility.
    
    int Numless[2][2];
    int Numbigger[2][2];
    int i,j,l,countplace,pairsNo;
    
    double n[2]={0,0};
    
    int rowLen, colLen1, colLen2;
    
    //double *output1,*output2,*output3,*output4,*output5, *output6, *output7, *output8;
    double x1,x2;
    double score,partialsecondaryscore[2],secondaryscore;
        
    rowLen = *rowLenPtr;
    colLen1 = *colLen1Ptr;
    colLen2 = *colLen2Ptr;
    
    
    pairsNo  = colLen1*colLen2;
    
    
    //Rprintf("situation[2]=%d,situation[3]=%d",situation[2],situation[3]);
    
    //return;
    //Rprintf("n[0]=%f,n[1]=%f\n",n[0],n[1]);
    
    for( l=0; l< rowLen; l++)
        n[situation[l]]++;
    
    //Rprintf("rowLen=%d,colLen1=%d,colLen2=%d",rowLen,colLen1,colLen2);
    //Rprintf("n[0]=%f,n[1]=%f",n[0],n[1]);
    //return;
    
    
    
    countplace =0;
    for(i =0 ; i< colLen1;i++)
    {
        
        for( j=0; j< colLen2;j++)
        {
            
            partialsecondaryscore[0] = 0;
            partialsecondaryscore[1] = 0;
            
            
            for( l = 0; l< 4; l++)
            {
                ((int*)Numless)[l]=0;
                ((int*)Numbigger)[l]=0;
            }
            
            for( l =0; l<rowLen;l++)
            {
                x1 = data1[(l*colLen1)+i];
                x2 = data2[(l*colLen2)+j];
                Numless[ situation[l]][x1<x2]++;
                Numbigger[ situation[l]][x1>x2]++;
                partialsecondaryscore[ situation[l]] += x2 - x1;
                
            }
            
            secondaryscore = (partialsecondaryscore[0]-partialsecondaryscore[1])/n[0]/n[1]/1e6;
            
            //Calculating score*/
            score = -Numless[0][0]/n[0]+Numless[1][0]/n[1]
            +Numbigger[0][0]/n[0]-Numbigger[1][0]/n[1];
            
            output1[ j*colLen1+i ] = score+secondaryscore;
            output7[ j*colLen1+i ] = (Numbigger[0][0])/n[0];
            output8[ j*colLen1+i ] = (Numbigger[1][0])/n[1];
            
            
            
            output3[ countplace ] = Numless[0][0];
            output5[ countplace ] = Numless[0][0] + Numless[1][0];
            
            countplace++;
            
            
        }
    }
    
    *output4 = n[0]+n[1];
    *output6 = n[0];
    
    
    
    }
    
}
