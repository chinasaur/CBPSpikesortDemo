#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double *p0, *p1, w0, w1, *pout0, *pout1;
    size_t n0, n1;
    int i, j;
    
    if (nrhs < 4)
        mexErrMsgTxt("You must provide inputs t0, t1, w0, w1.");
    if (nlhs < 2)
        mexErrMsgTxt("You must provide two output variables.");
    if (!mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1]) || !mxIsDouble(prhs[2]) || !mxIsDouble(prhs[3]))
        mexErrMsgTxt("Inputs must all be doubles.");
    

    p0  = mxGetPr(prhs[0]);
    p1  = mxGetPr(prhs[1]);
    n0 = mxGetNumberOfElements(prhs[0]);
    n1 = mxGetNumberOfElements(prhs[1]);
    w0 = mxGetScalar(prhs[2]);
    w1 = mxGetScalar(prhs[3]);
    
    plhs[0] = mxCreateDoubleMatrix(n0, 1, mxREAL);
    pout0 = mxGetPr(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(n1, 1, mxREAL);
    pout1 = mxGetPr(plhs[1]);
    for (i = 0, j = 0; i < n1; ++i) {
        double t = p1[i];
        
        /* Advance j0 past beginning of window; anything before we get to
         * window is a miss in p0, so mark as 0 */
        for (; j < n0; ++j) {
            if (p0[j] - t >= w0) break;
            pout0[j] = 0;
        }
        
        /* If j has fallen off the end of p0, set the rest of p1 to 0 and exit */
        if (j >= n0) {
            for (; i < n1; ++i) pout1[i] = 0;
            break;
        }
        
        /* If we are still before end of window, record match */
        if (p0[j] - t <= w1) {
            /* Match; mark p0 and p1; 1-index for Matlab... */
            pout0[j] = i+1;
            pout1[i] = j+1;
            
            /* Don't use spike more than once, so advance j0 */
            ++j;
            
        } else {
            /* No match, so mark a miss in p1 */
            pout1[i] = 0;
        }
    }

}