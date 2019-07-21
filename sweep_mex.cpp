#include "mex.h"
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <limits>

struct greater2nd {
  template <typename P> bool operator() (const P& p1, const P& p2) {
    return p1.second > p2.second;
  }
};

void mexFunction (int nlhs, mxArray *plhs[],
             int nrhs, const mxArray *prhs[]) {

    mwSize m1, n1, m2, n2, nz2, K;
    mwIndex *jc1, *ir1, *jc2, *ir2;
    double *pr1, *pr2;
    size_t i=0;
    
        
    m1 = mxGetM (prhs[0]);
    n1 = mxGetN (prhs[0]);
    ir1 = mxGetIr (prhs[0]);
    jc1 = mxGetJc (prhs[0]);
    pr1 = mxGetPr (prhs[0]);
    m2 = mxGetM (prhs[1]);
    n2 = mxGetN (prhs[1]);
    nz2 = mxGetNzmax (prhs[1]);
    ir2 = mxGetIr(prhs[1]);
    jc2 = mxGetJc(prhs[1]);
    pr2 = mxGetPr(prhs[1]);
    
    if (nrhs < 2) mexErrMsgTxt ("Must be at least two ARGs");
    if (!mxIsSparse (prhs[0])) mexErrMsgTxt ("ARG1 must be sparse");
    if (m1 != n1) mexErrMsgTxt ("ARG1 must be a square matrix");
    if (!mxIsSparse (prhs[1])) mexErrMsgTxt ("ARG2 must be sparse");
    if (m2 != m1) mexErrMsgTxt ("ARG1 and ARG2 must be matched");
    if (n2 != 1) mexErrMsgTxt ("ARG2 must be a column vector");
    if (nrhs == 3) K = mxGetScalar(prhs[2]); else K = nz2;
    if (K > nz2) K = nz2;
//     printf("%d\n", K);
    
    // Copy vector to map
    std::unordered_map<mwIndex,double> map;
    for (i=0; i<nz2; i++) {
        map[ir2[i]] = pr2[i] / (jc1[ir2[i]+1] - jc1[ir2[i]]);
//         map[ir2[i]] = pr2[i];
    }

    std::vector< std::pair<int, double> > vec(map.begin(), map.end());
    std::sort(vec.begin(), vec.end(), greater2nd());
    
//     for (auto it=vec.begin(); it!=vec.end(); it++)
//         printf("%d: %f\n", it->first, it->second);
    
    //compute cutsize, volume, and cond
    std::vector<double> cond(K);

    std::unordered_map<int,size_t> rank;
    i = 0;
    for (auto it=vec.begin(); i<K; ++it, ++i) {
        rank[it->first] = i;
    }
    mwIndex total_degree = jc1[m1];
    mwIndex cut = 0;
    mwIndex vol = 0;
    i = 0;
    for (auto it=vec.begin(); i<K; ++it, ++i) {
        mwIndex v = it->first;
        mwIndex deg = jc1[v+1] - jc1[v];
        mwIndex dif = deg;
        for (mwIndex j=jc1[v]; j<jc1[v+1]; ++j) {
            mwIndex nbr = ir1[j];
            if (rank.count(nbr) > 0) {
                if (rank[nbr] < rank[v]) {
                    dif -= 2;
                }
            }
        }
        cut += dif;
//         printf("%d, %d\n", v, cut);
        vol += deg;
        if (vol == 0 || total_degree-vol==0) {
            cond[i] = 1;
        } else {
            cond[i] = (double)cut/
                    (double)std::min(vol,total_degree-vol);
        }
        
    }
//     for (i=0; i<cond.size(); i++)
//             printf("%d: %f\n", vec[i].first, cond[i]);
    size_t lastind = i;
    double mincond = std::numeric_limits<double>::max();
    size_t mincondind = 0; // set to zero so that we only add one vertex 
    for (i=0; i<K; i++) {
            if (cond[i] < mincond) {
            mincond = cond[i];
            mincondind = i;
        }
    }
//     printf("mincondid = %d, mincond = %f\n", mincondind, mincond);
    
    // Output community
    plhs[0] = mxCreateNumericMatrix(mincondind+1, 1, mxINT32_CLASS, mxREAL);
    int *community  = (int*) mxGetData(plhs[0]);
    for (int i=0; i<=mincondind; i++) {
        community[i] = vec[i].first + 1;
    }
}

