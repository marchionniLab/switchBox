/* Registration of C routines */

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

void CalculateSignedScoreCore(int *situation, int *rowLenPtr, double *data1, int *colLen1Ptr, double *data2, int *colLen2Ptr,
                              double *output1, double *output3, double *output4, double *output5, double *output6, double *output7, double *output8);

void CalculateSignedScoreRestrictedPairsCore(int *situation, int *rowLenPtr, double *data1, int *colLen1Ptr, double *data2,
                                             int *colLen2Ptr, int *edges1, int *edges2, int *nopairsPtr,
                                             double *output1, double *output2, double *output3);

void vect2compC( double *vect, int *np, int *mp, double *comp);


#if _MSC_VER >= 1000
__declspec(dllexport)
#endif

static const R_CMethodDef cMethods[] = {
    {"CalculateSignedScoreCore", (DL_FUNC)&CalculateSignedScoreCore,13},
    {"CalculateSignedScoreRestrictedPairsCore", (DL_FUNC)&CalculateSignedScoreRestrictedPairsCore,12},
    {"vect2compC", (DL_FUNC)&vect2compC,4},
    {NULL, NULL, 0},
};

void R_init_Test(DllInfo *info)
{
  R_registerRoutines(info, cMethods, NULL, NULL, NULL);
}
