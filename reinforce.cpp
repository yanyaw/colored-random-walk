#include "mex.h"

void mexFunction (int nlhs, mxArray *plhs[],
             int nrhs, const mxArray *prhs[]) {
    mwSize m, n, nz, m2, n2;
    mwIndex *jc, *ir, *jc2, *ir2;
    double *pr, *pr2, *b;
    mxArray *v;
        
    m = mxGetM (prhs[0]);
    n = mxGetN (prhs[0]);
    nz = mxGetNzmax (prhs[0]);
    m2 = mxGetM (prhs[1]);
    n2 = mxGetN (prhs[1]);
    
    if (nrhs != 2) mexErrMsgTxt ("Must be two ARGs");
    if (!mxIsSparse (prhs[0])) mexErrMsgTxt ("ARG1 must be sparse");
    if (m != n) mexErrMsgTxt ("ARG1 must be a square matrix");
    if (mxIsSparse (prhs[1])) mexErrMsgTxt ("ARG2 must be full");
    if (m2 != m) mexErrMsgTxt ("ARG1 and ARG2 must be matched");
    if (n2 != 1) mexErrMsgTxt ("ARG2 must be a column vector");
    
    jc = mxGetJc(prhs [0]); // column pointers
    ir = mxGetIr(prhs [0]); // row indices
    pr = mxGetPr(prhs [0]); // values
    
    b = mxGetPr(prhs[1]);
    
    v = mxCreateSparse(m, n, nz, mxREAL);
    ir2 = mxGetIr(v);
    jc2 = mxGetJc(v);
    pr2 = mxGetPr(v);
    
    int i, j, r;
    double s, w;
//     for (i=0; i<nz; i++) printf("%f ", pr2[i]);
    for (i=0; i<n; i++) {
        s = 0;
        for (j=jc[i]; j<jc[i+1]; j++) {
            r = ir[j];
            ir2[j] = r;
            w = pr[j] + pr[j] * b[r];
            if (w > 0) {
                s += w;
                pr2[j] = w;
            }
        }
        if (s > 0) s = 1.0 / s;
        for (j=jc[i]; j<jc[i+1]; j++) {
            pr2[j] = pr2[j] * s;
        }
    }
    for (i=0; i<n+1; i++) jc2[i] = jc[i];
    if (nlhs > 0) plhs[0] = v;
}