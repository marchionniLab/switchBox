#include <R.h>
#include "math.h"


inline double max(double x, double y){ return ((x<y)? y:x);}
inline double myabs(double x){ return ((x>0)?x:(-x));}
#define myepsilon 1e-9
// comp = vect2comp(vect);
// comp is comp(i,j)=\delta{vect(i)<vect(j)}
extern "C"
{
    
    
    
/*// [score, P, Q, type] = CalculateSignedScoreCoreTieHandler(situation,data(indexup,:),data(indexdown,:));*/
    
    
    
void CalculateSignedScoreCoreTieHandler(int *situation, int *rowLenPtr, double *data1,int *colLen1Ptr, double *data2, int *colLen2Ptr,//inputs
    double *output1, double *output2, double *output3, double *output4)//outputs
{
    //void CalculateSignedScoreCoreTieHandler(int *situation, int rowLen, double *data1,int colNo1, double *data2, int colNo2,//inputs
    //double *score, double *P, double *Q, double *type)//outputs
    //Calculates pairwise score between every gene in data1 and every gene in data2. It is signed because
    // score(i,j)= \max( P(X_i<X_j|1)-P(X_i<=X_j|0),  P(X_i<X_j|1)-P(X_i<X_j|0),P(X_i<=X_j|1)-P(X_i<=X_j|0)) + C
	// and type(i,j) \in \{0,1,2\}indicates which one of these probabilities was the maximum, i.e.,
	// type(i,j)= \arg\max( P(X_i<X_j|1)-P(X_i<=X_j|0),  P(X_i<X_j|1)-P(X_i<X_j|0),P(X_i<=X_j|1)-P(X_i<=X_j|0)) 
	// In case of a tie (or ties), we will choose the first one. In other words, 0 is preferred to 1 and 1 to 2.  
	// If P(X_i=X_j)=0 (which is a common case unless we have discrete data), 0 will be chosen.
    // The C term is proportional to the secondary score
    // to avoid the ties among pairs.  C = (E(X_j-X_i|1)-E(X_j-X_i|0))/K where K is big enough to make sure that the secondary score
    // does not intervene with the primary score.
    // P_ij = P(X_i<X_j|0) and Q(X_i<X_j|1)
    
    // Note that the colLen1, colLen2 and rowLen represent row and column length respectively. It is due to the early changes in the program that 
    // they maintain their cabapility.
    
    int Numless[2][2];
    int Numeq[2];
    int i,j,l;
    
    double n[2]={0,0};
    
    int rowLen, colLen1, colLen2;
    
    //double *output1,*output2,*output3,*output4,*output5, *output6, *output7, *output8;
    double diffx;
	double score, scorebase, scoretie,  partialsecondaryscore[2],secondaryscore;
    double score1, score2, type;
	
    rowLen = *rowLenPtr;
    colLen1 = *colLen1Ptr;
    colLen2 = *colLen2Ptr;
    
    
     
    
    //Rprintf("situation[2]=%d,situation[3]=%d",situation[2],situation[3]);
    
    //return;
    //Rprintf("n[0]=%f,n[1]=%f\n",n[0],n[1]);
    
    for( l=0; l< rowLen; l++)
        n[situation[l]]++;
    
    //Rprintf("rowLen=%d,colLen1=%d,colLen2=%d",rowLen,colLen1,colLen2);
    //Rprintf("n[0]=%f,n[1]=%f",n[0],n[1]);
    //return;
    
    
    
    for(i =0 ; i< colLen1;i++)
    {
        
        for( j=0; j< colLen2;j++)
        {
            
            partialsecondaryscore[0] = 0;
            partialsecondaryscore[1] = 0;
            
            
			// Initializing Numless and Numeq
			Numless[0][0] = 0;
			Numless[0][1] = 0;
			Numless[1][0] = 0;
			Numless[1][1] = 0;
			Numeq[0] = 0;
			Numeq[1] = 0;
            
            for( l =0; l<rowLen;l++)
            {
                diffx = data2[(l*colLen2)+j]-data1[(l*colLen1)+i]; // Calculating the difference
				if( myabs(diffx ) < myepsilon)//if there is a tie
					Numeq[ situation[l]]++;
				else Numless[ situation[l]][diffx>0]++;//if there is not update Numless
                
                partialsecondaryscore[ situation[l]] += diffx;
                
            }
            
            secondaryscore = (partialsecondaryscore[0]-partialsecondaryscore[1])/n[0]/n[1]/1e6;
            
            //Calculating score*/
            			
			scorebase = Numless[1][0]/n[1] - Numless[0][0]/n[0];
			scoretie = Numeq[1]/n[1] - Numeq[0]/n[0];
			score2 = scorebase + scoretie;// tie in favor of class 0, we use score to reduce the computation time
			
			if( myabs(scoretie)<myepsilon)//This is the only case where randomization is better
			{
				score = scorebase + 0.5 * scoretie; //randomized tie
				type = 0;
			}else if(myabs(scorebase) < myabs(score2) ){
				
				score = score2;
				type = 1;
			}else{
				score = scorebase;
				type = 2;
			}

			//Rprintf("(%d,%d) %3g (%3g)[%3g] type %f\n",i+1,j+1,score,scorebase,scoretie,type);

			
            output1[ j*colLen1+i ] = score+secondaryscore;
            output2[ j*colLen1+i ] = 1-(Numless[0][0])/n[0];
            output3[ j*colLen1+i ] = 1-(Numless[1][0])/n[1];
            output4[ j*colLen1+i ] = type;

            
            
                 
            
            
            
        }
    } 
    
    
    
    }
    
}
