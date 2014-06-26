#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double *pE, *pT, t0, t1, *poutE, *poutT;
    size_t nE, nT;
    int i, j0, j1, sz, k;
    
    if (nrhs < 4)
        mexErrMsgTxt("You must provide inputs events, triggers, t0, t1.");
    if (nlhs < 2)
        mexErrMsgTxt("You must provide two output variables.");
    if (!mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1]) || !mxIsDouble(prhs[2]) || !mxIsDouble(prhs[3]))
        mexErrMsgTxt("Inputs must all be doubles.");
    

    pE  = mxGetPr(prhs[0]);
    pT  = mxGetPr(prhs[1]);
    nE = mxGetNumberOfElements(prhs[0]);
    nT = mxGetNumberOfElements(prhs[1]);
    t0 = mxGetScalar(prhs[2]);
    t1 = mxGetScalar(prhs[3]);
    
    /* Initial pass to see how big an output we need */
    for (i = 0, j0 = 0, sz = 0; i < nT; ++i) {
        double trig = pT[i];
        
        /* Advance j0 to beginning of window */
        for (j0 = 0; j0 < nE; ++j0)
            if (pE[j0] - trig >= t0) break;

        /* Advance j1 to end of window */
        for (j1 = j0; j1 < nE; ++j1)
            if (pE[j1] - trig > t1) break;

        /* Note how many events in window */
        sz += j1 - j0;
    }
    
    /* Fill output, similar logic to above */
    plhs[0] = mxCreateDoubleMatrix(sz, 1, mxREAL);
    poutE = mxGetPr(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(sz, 1, mxREAL);
    poutT = mxGetPr(plhs[1]);
    for (i = 0, j0 = 0, k = 0; i < nT; ++i) {
        double trig = pT[i];
        
        for (j0 = 0; j0 < nE; ++j0)
            if (pE[j0] - trig >= t0) break;

        for (j1 = j0; j1 < nE; ++j1) {
            double diff = pE[j1] - trig;
            if (diff > t1) break;
            poutE[k]   = diff;
            poutT[k++] = i+1;
        }
    }

}