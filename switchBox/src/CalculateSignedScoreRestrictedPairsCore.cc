#include <R.h>
#include "math.h"

#define max(x,y) ((x)>(y))?(x):y

/*// [score, P, Q] = CalculateSignedScoreRestrictedPairsCore(situation,data(indexup,:),data(indexdown,:),...
 int32(edges(:,1)-1),int32(edges(:,2))-1);*/
// It calculates the restricted scores for pairs with edges(:,1) and edges(:,2)

extern "C"
{


	void CalculateSignedScoreRestrictedPairsCore(int *situation, int *rowLenPtr, double *data1,int *colLen1Ptr, double *data2, 
		int *colLen2Ptr,int *edges1, int *edges2, int *nopairsPtr,//inputs
	    double *output1, double *output2, double *output3)//outputs

	{

		int Numless[2][2];
	    	int Numbigger[2][2];
		int i,l,countplace,pairsNo;
		double n[2]={0,0};
	    	//double *data1, *data2;
		int rowLen, colLen1, colLen2;
		//int  *situation;
		//int  *edges1, *edges2;
		int edges1i, edges2i;
		//double *output1,*output2,*output3;
		double x1,x2;
		double score,partialsecondaryscore[2],secondaryscore;
    
    
    
        
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
    
    
		countplace =0;
		for(i =0 ; i< pairsNo;i++)
		{
        
        
			partialsecondaryscore[0] = 0;
		    partialsecondaryscore[1] = 0;
        
        	for( l = 0; l< 4; l++)
			{
				((int*)Numless)[l]=0;
				((int*)Numbigger)[l]=0;
			}
        
			edges1i = edges1[i];
			edges2i = edges2[i];
            
    

        
       
			for( l =0; l<rowLen;l++)
			{
				x1 = data1[(l*colLen1)+edges1i];
				x2 = data2[(l*colLen2)+edges2i];
            
				Numless[ situation[l]][x1<x2]++;
				Numbigger[ situation[l]][x1>x2]++;
				partialsecondaryscore[ situation[l]] += x2 - x1;
                
                //Rprintf("sample %d: index gene 1:%d, index gene 2:%d \n",l,(l*colLen1)+edges1i,(l*colLen2)+edges2i);
            
			}
        
			secondaryscore = (partialsecondaryscore[0]-partialsecondaryscore[1])/n[0]/n[1]/1e6;
             
			//Calculating score
			score = -Numless[0][0]/n[0]+Numless[1][0]/n[1]
				  +Numbigger[0][0]/n[0]-Numbigger[1][0]/n[1];
        
             
			//Rprintf("gene %d gene %d Numless[0][0]=%d Numless[1][0]= %d Numbigger[0][0]=%d Numbigger[1][0]=%d,n[0]=%2.2g,n[1]=%2.2g\n", 
			//			Numless[0][0],Numless[1][0],Numbigger[0][0],Numbigger[1][0],n[0],n[1]);
             
			output1[ i ] = score+secondaryscore;
			output2[ i ] = (Numbigger[0][0])/n[0];
			output3[ i ] = (Numbigger[1][0])/n[1];
                
        
		}
    
    
    
   
    
    
}


}
