/***************************************************************************/
/*      Name:       DCMF_Mex.c                               
                                                
        Performing the continuous max-flow algorithm to solve the
        continuous min-cut problem in 2D
 
        Usage: [u, erriter, i, timet, ps, pt, px, py] = 
 *                          DCMF_Mex(penalty, C_s, C_t, para, ps_in, 
 *                          pt_in, px_in, py_in, u_in);
 
        Inputs (penalty, C_s, C_t, para, ps_in, pt_in, px_in, py_in, u_in): 
 
               - penalty: point to the edge-weight penalty parameters to
                          total-variation function.
 
                 For the case without incorporating image-edge weights, 
                 penalty is given by the constant everywhere. For the case 
                 with image-edge weights, penalty is given by the pixelwise 
                 weight function:
 
                 for example, penalty(x) = b/(1 + a*| grad f(x)|) where b,a > 0 .
   
               - C_s: point to the capacities of source flows ps
 
               - C_t: point to the capacities of sink flows pt
 
               - para: a sequence of parameters for the algorithm
                    para[0,1]: rows, cols of the given image
                    para[2]: the maximum iteration number
                    para[3]: the error bound for convergence
                    para[4]: cc for the step-size of augmented Lagrangian method
                    para[5]: the step-size for the graident-projection step to the
                           total-variation function. Its optimal range is [0.1, 0.17].
 
        Outputs (u, erriter, i, timet):
 
              - u: the final results u(x) in [0,1]. As the following paper,
                   the global binary result can be available by threshholding u
                  by any constant alpha in (0,1):

                  Nikolova, M.; Esedoglu, S.; Chan, T. F. 
                  Algorithms for Finding Global Minimizers of Image Segmentation and Denoising Models 
                  SIAM J. App. Math., 2006, 66, 1632-1648

               - erriter: it returns the error evaluation of each iteration,
                  i.e. it shows the convergence rate. One can check the algorithm
                  performance.

               - i: gives the total number of iterations, when the algorithm converges.

               - timet: gives the total computation time.

       This software can be used only for research purposes, you should cite ALL of
       the aforementioned papers in any resulting publication.

       Please email Jing Yuan (cn.yuanjing@gmail.com) for any questions, suggestions and bug reports

       The Software is provided "as is", without warranty of any kind.


                           Version 1.0
               https://sites.google.com/site/wwwjingyuan/       

               Copyright 2011 Jing Yuan (cn.yuanjing@gmail.com)      

*/
/***************************************************************************/




/*  compilation command (under matlab): mex DCMF_mex02.c  */

#include <stdio.h>
#include <stdlib.h>
#include <mex.h>
#include <math.h>
#include <time.h>


#define YES 0
#define NO 1

#define PI 3.1415926

#define MAX(a,b) ( a > b ? a : b )
#define MIN(a,b) ( a <= b ? a : b )
#define SIGN(x) ( x >= 0.0 ? 1.0 : -1.0 )
#define ABS(x) ( (x) > 0.0 ? x : -(x) )
#define SQR(x) (x)*(x)

#ifndef HAVE_RINT 
#define rint(A) floor((A)+(((A) < 0)? -0.5 : 0.5)) 
#endif

float SQRT(float number) {
    long i;
    float x, y;
    const float f = 1.5F;
    
    x = number * 0.5F;
    y  = number;
    i  = * ( long * ) &y;
    i  = 0x5f3759df - ( i >> 1 );
    y  = * ( float * ) &i;
    y  = y * ( f - ( x * y * y ) );
    y  = y * ( f - ( x * y * y ) );
    return number * y;
}


/**********************************************/
/************** MAIN FUNCTION *****************/
/**********************************************/

/****************************************/
extern void mexFunction(int iNbOut, mxArray *pmxOut[],
int iNbIn, const mxArray *pmxIn[])
{
    
  /* iNbOut: number of outputs
     pmxOut: array of pointers to output arguments */
    
  /* iNbIn: number of inputs
     pmxIn: array of pointers to input arguments */
    
    
    float   *pfpenalty, *pfu, *pfCs, *pfCt, *pfcvg, *pfVecParameters;
    float   *pfbx, *pfby, *pfps, *pfpt, *pfgk, *tt, *pfdv;
    float   fLambda, fError, cc, steps;
    float   fpt, fps;
    int     *punum, iNy, iNx, iNdim, iDim[3], ix, iy, i, iNI;
    int     iNbIters, szImg, idx, index;
    time_t  start_time, end_time;
    
    /* Inputs */
    pfpenalty = mxGetData(pmxIn[0]); /* Given penalty */
    pfCs = mxGetData(pmxIn[1]); /* bound of source flows */
    pfCt = mxGetData(pmxIn[2]); /* bound of sink flows */
    pfVecParameters = mxGetData(pmxIn[3]); /* Vector of parameters */
    
    
    /* 
     *pfVecParameters Setting
     * [0] : number of columns 
     * [1] : number of rows
     * [2] : the maximum iteration number
     * [3] : error criterion
     * [4] : cc for the step-size of ALM
     * [5] : steps for the step-size of projected-gradient of p
     */
    
    /* Size */
    iNy = (int) pfVecParameters[0];
    iNx = (int) pfVecParameters[1];
    szImg = iNy * iNx;
    
    iNbIters = (int) pfVecParameters[2]; 
    fError = (float) pfVecParameters[3]; 
    cc = (float) pfVecParameters[4]; 
    steps = (float) pfVecParameters[5]; 
    
   
    /* Outputs */
    /* outputs the computed u(x)  */
    iNdim = 2;
    iDim[0] = iNy;
    iDim[1] = iNx;
    
    pmxOut[0] = mxDuplicateArray(pmxIn[8]);
    pfu = mxGetData(pmxOut[0]);
               
    /* outputs the convergence rate  */
    iNdim = 2;
    iDim[0] = 1;
    iDim[1] = iNbIters;
    pmxOut[1] = mxCreateNumericArray(iNdim,(const int*)iDim,mxSINGLE_CLASS,mxREAL);
    pfcvg = mxGetData(pmxOut[1]);
    
    /* outputs the iteration number  */
    iNdim = 2;
    iDim[0] = 1;
    iDim[1] = 1;
    pmxOut[2] = mxCreateNumericArray(iNdim,(const int*)iDim,mxUINT16_CLASS,mxREAL);
    punum = mxGetData(pmxOut[2]);
    
    /* outputs the computation time  */
    iNdim = 2;
    iDim[0] = 1;
    iDim[1] = 1;
    pmxOut[3] = mxCreateNumericArray(iNdim,(const int*)iDim,mxSINGLE_CLASS,mxREAL);
    tt = mxGetData(pmxOut[3]);
    
    /* outputs the four flows: ps, pt, px, py */
    iNdim = 2;
    iDim[0] = iNy;
    iDim[1] = iNx;
    
    pmxOut[4] = mxDuplicateArray(pmxIn[4]);
    pfps = mxGetData(pmxOut[4]);
    
    iNdim = 2;
    iDim[0] = iNy;
    iDim[1] = iNx;
    
    pmxOut[5] = mxDuplicateArray(pmxIn[5]);
    pfpt = mxGetData(pmxOut[5]);
    
    iNdim = 2;
    iDim[0] = iNy;
    iDim[1] = iNx+1;
    
    pmxOut[6] = mxDuplicateArray(pmxIn[6]);
    pfbx = mxGetData(pmxOut[6]);
    
    iNdim = 2;
    iDim[0] = iNy+1;
    iDim[1] = iNx;
    
    pmxOut[7] = mxDuplicateArray(pmxIn[7]);
    pfby = mxGetData(pmxOut[7]);
    
    /* Memory allocation */
        
    /* allocate the memory for gk */
    pfgk = (float *) calloc( (unsigned)(iNy*iNx), sizeof(float) );
    if (!pfgk)
        mexPrintf("Memory allocation failure\n");
    
    /* allocate the memory for div */
    pfdv = (float *) calloc( (unsigned)(iNy*iNx), sizeof(float) );
    if (!pfdv)
        mexPrintf("Memory allocation failure\n");
    
    /* Preparing initial values for flows and u */
    for (ix=0; ix< iNx; ix++){
        idx = ix*iNy;
        for (iy=0; iy< iNy; iy++){
           index = idx + iy; 

                pfdv[index] = pfbx[index+iNy] - pfbx[index] 
                         + pfby[index+1] - pfby[index];
        }
    }
    
    /*  Main iterations */
    
    iNI = 0;
     
    start_time = clock();
    
    while( iNI<iNbIters ) 
    { 
        
        /* update p */
        
        for (ix=0; ix< iNx; ix++){
            idx = ix*iNy;
            for (iy=0; iy< iNy; iy++){
                index = idx + iy;
                pfgk[index] = pfdv[index] - (pfps[index] - pfpt[index] 
                                 + pfu[index]/cc);
                
            }
        }
        
        /* update px */
        
        for (ix=0; ix < iNx-1; ix++){
            idx = (ix+1)*iNy;
            for (iy=0; iy < iNy; iy++){
                index = idx + iy;
                pfbx[index] = steps*(pfgk[index] - pfgk[index-iNy]) + pfbx[index];
            }
        }
    
        /* update py */

        for(ix = 0; ix < iNx; ix ++){
            idx = ix*iNy;
            for(iy = 0; iy < iNy-1; iy ++){
                index = idx + iy + 1;
                pfby[index] = steps*(pfgk[index] - pfgk[index-1]) + pfby[index];
            }
        }
        
        /* projection step */
        
        for (ix = 0; ix < iNx; ix++){
            idx = ix*iNy;
            for (iy = 0; iy < iNy; iy++){
                index = idx + iy;
                fpt = SQRT((SQR(pfbx[index]) + SQR(pfbx[index+iNy]) 
                      + SQR(pfby[index]) + SQR(pfby[index+1]))*0.5);
                
                if (fpt > pfpenalty[index])
                    fpt = fpt / pfpenalty[index];
                else
                    fpt = 1;
                
                pfgk[index] = 1/fpt;
            }
        }
        
        for (ix =0; ix < iNx-1; ix++){
            idx = (ix+1)*iNy;
            for (iy = 0; iy < iNy; iy++){
                index = idx + iy;
                pfbx[index] = (pfgk[index] + pfgk[index-iNy])*0.5*pfbx[index];
            }
        }
        
        for (ix = 0; ix < iNx; ix++){
            idx = ix*iNy;
            for (iy = 0; iy < iNy-1; iy++){
                index = idx + iy + 1;
                pfby[index] = (pfgk[index-1] + pfgk[index])*0.5*pfby[index];
            }
        }      
        
        /* compute the divergence  */
        
        for (ix = 0; ix < iNx; ix++){
            idx = ix*iNy;
            for (iy = 0; iy < iNy; iy++){
                index = idx + iy;
                pfdv[index] = pfbx[index+iNy] - pfbx[index] 
                         + pfby[index+1] - pfby[index];
            }
        }
        
        /* update ps  */
        
        for (ix = 0; ix < iNx; ix++){
            idx = ix*iNy;
            for (iy = 0; iy < iNy; iy++){
                index = idx + iy;
                fpt = pfpt[index] - pfu[index]/cc + pfdv[index] + 1/cc;
                fpt = MIN(fpt , pfCs[index]);
                pfps[index] = fpt;
            }
        }
        
        /* update pt  */
        
        for (ix = 0; ix < iNx; ix++){
            idx = ix*iNy;
            for (iy = 0; iy < iNy; iy++){
                index = idx + iy;
                fpt = pfps[index] + pfu[index]/cc - pfdv[index];
                fpt = MIN(fpt , pfCt[index]);
                pfpt[index] = fpt;
            }
        }
   
        /* update multipliers */
        
        fps = 0;
        for (ix = 0; ix < iNx; ix++){
            idx = ix*iNy;
            for (iy = 0; iy < iNy; iy++){
                index = idx + iy;
                fpt = cc*(pfpt[index] + pfdv[index] - pfps[index]);
                fps += ABS(fpt);
                
                pfu[index] -= fpt;
            }
        }
        

        pfcvg[iNI] = fps / szImg;
        
        if ((pfcvg[iNI] <= fError) && (iNI >= 2))
            break;
                
        iNI ++;
     }   
    
    mexPrintf("number of iterations = %i\n",iNI);
    
    /* Outputs (see above) */
    
    punum[0] = iNI;
    
    /* Free memory */
    
    free( (float *) pfgk );
    free( (float *) pfdv );
    
    end_time = clock();
    tt[0] = difftime(end_time, start_time)/1000;
    
    mexPrintf("\nComputing Time for max-flow = %.4f sec\n \n", tt[0]);
    
    
}
/****************************************/


