#include <R.h>
#include "math.h"

inline double max(double x, double y){ return ((x<y)? y:x);}
inline double myabs(double x){ return ((x>0)?x:(-x));}
#define myepsilon 1e-9

/*// [score, P, Q, type] = CalculateSignedScoreRestrictedPairsCoreTieHandler(situation,data(indexup,:),data(indexdown,:),...
 int32(edges(:,1)-1),int32(edges(:,2))-1);*/
// It calculates the restricted scores for pairs with edges(:,1) and edges(:,2)

extern "C"
{


	void CalculateSignedScoreRestrictedPairsCoreTieHandler(int *situation, int *rowLenPtr, double *data1,int *colLen1Ptr, double *data2, 
		int *colLen2Ptr,int *edges1, int *edges2, int *nopairsPtr,//inputs
	    double *output1, double *output2, double *output3, double *output4)//outputs

	{

		int Numless[2][2];
        int Numeq[2];

		int i,l,countplace,pairsNo;
		double n[2]={0,0};
	    	//double *data1, *data2;
		int rowLen, colLen1, colLen2;
		//int  *situation;
		//int  *edges1, *edges2;
		int edges1i, edges2i;
		//double *output1,*output2,*output3;
		double diffx;
		double score, scorebase, scoretie,  partialsecondaryscore[2],secondaryscore;
		double score1, score2, type;
    
            
		rowLen = *rowLenPtr;
		colLen1 = *colLen1Ptr;
		colLen2 = *colLen2Ptr;
    
    
		pairsNo  = *nopairsPtr;
	
		    
		//Rprintf("n=%d",rowLen);
    
		for( l=0; l< rowLen; l++)
		{
			n[situation[l]]++;
			//Rprintf("%d,",situation[l]);
        }
    
		//Rprintf("\n");	
    
    
		for(i =0 ; i< pairsNo;i++)
		{
        
        
			partialsecondaryscore[0] = 0;
		    partialsecondaryscore[1] = 0;
        
        	
			Numless[0][0] = 0;
			Numless[0][1] = 0;
			Numless[1][0] = 0;
			Numless[1][1] = 0;
			Numeq[0] = 0;
			Numeq[1] = 0;
                
			edges1i = edges1[i];
			edges2i = edges2[i];
            
    

        
       
			for( l =0; l<rowLen;l++)
			{
				diffx = data2[(l*colLen2)+edges2i]-data1[(l*colLen1)+edges1i]; // Calculating the difference
				if( myabs(diffx ) < myepsilon)//if there is a tie
					Numeq[ situation[l]]++;
				else Numless[ situation[l]][diffx>0]++;//if there is not update Numless
                partialsecondaryscore[ situation[l]] += diffx;
                
                
                //Rprintf("sample %d: index gene 1:%d, index gene 2:%d \n",l,(l*colLen1)+edges1i,(l*colLen2)+edges2i);
            
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
			}
			else if(myabs(scorebase) < myabs(score2) ){
				score = score2;
				type = 1;
			}else{
				score = scorebase;
				type = 2;				
			}

			//Rprintf("gene %d gene %d Numless[0][0]=%d Numless[1][0]= %d Numbigger[0][0]=%d Numbigger[1][0]=%d,n[0]=%2.2g,n[1]=%2.2g\n", 
			//			Numless[0][0],Numless[1][0],Numbigger[0][0],Numbigger[1][0],n[0],n[1]);
             
			output1[ i ] = score+secondaryscore;
			output2[ i ] = 1-(Numless[0][0])/n[0];
			output3[ i ] = 1-(Numless[1][0])/n[1];
   			output4[ i ] = type;
    
        
		}
    
    
    
   
    
    
}


}
