#include "mex.h"
#include <algorithm>  // for std::max
#include <cmath> // for std::abs


// Normal force
inline double NormalForces(double xn, double kn, double xn0) {
    double u = xn - xn0;
    return std::max(0.0, kn * u);
}

inline void TangentialForces(double xt, double &w, double kt, double mu, double FN, double &T) {
    if (FN > 0) {
        T = kt * (xt - w);
        if (std::abs(T) < mu * FN) {
            w = w; // sticking case
        } else {
            T = (T > 0 ? 1 : -1) * mu * FN; // sliding case
            w = xt - T / kt;
        }
    } else {
        T = 0;
        w = xt;
    }
}

// Helper: get field/property from struct or object
const mxArray* getFieldOrProperty(const mxArray *obj, const char *name) {
    if (mxIsStruct(obj)) {
        return mxGetField(obj, 0, name);
    } else if (mxIsClass(obj, "handle") || mxIsObject(obj)) {
        return mxGetProperty(obj, 0, name);
    }
    mexErrMsgIdAndTxt("gf_mex:inputType", "fc must be a struct or object.");
    return nullptr;
}

// Entry point
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    if (nrhs != 2) {
        mexErrMsgIdAndTxt("gf_mex:nrhs", "Two inputs required: x, fc.");
    }
    if (nlhs != 2) {
        mexErrMsgIdAndTxt("gf_mex:nlhs", "Two outputs required: F, w.");
    }

    // === Input x ===
    const mxArray *x_mx = prhs[0];
    if (!mxIsDouble(x_mx) || mxIsComplex(x_mx)) {
        mexErrMsgIdAndTxt("gf_mex:xType", "x must be a real double vector.");
    }
    double *x = mxGetPr(x_mx);
    mwSize Nx = mxGetNumberOfElements(x_mx);

    // === Input fc struct/object ===
    const mxArray *fc = prhs[1];

    const mxArray *kn_mx  = getFieldOrProperty(fc, "kn");
    const mxArray *xn0_mx = getFieldOrProperty(fc, "xn0");
    const mxArray *mu_mx  = getFieldOrProperty(fc, "mu");
    const mxArray *kt_mx  = getFieldOrProperty(fc, "kt");
    const mxArray *w_mx   = getFieldOrProperty(fc, "w");

    if (!kn_mx || !xn0_mx || !mu_mx || !kt_mx || !w_mx) {
        mexErrMsgIdAndTxt("gf_mex:missingField",
                          "fc must contain fields/properties kn, xn0, mu, kt, w.");
    }

    double *kn  = mxGetPr(kn_mx);
    double *xn0 = mxGetPr(xn0_mx);
    double *mu  = mxGetPr(mu_mx);
    double *kt  = mxGetPr(kt_mx);
    double *w   = mxGetPr(w_mx);

    mwSize Np = mxGetNumberOfElements(kn_mx); // contact points
    if (Nx != 3*Np) {
        mexErrMsgIdAndTxt("gf_mex:dimMismatch", "x length must be 3*Np.");
    }

    // === Output F ===
    plhs[0] = mxCreateDoubleMatrix(3*Np, 1, mxREAL);
    double *F = mxGetPr(plhs[0]);

    // === Output w ===
    plhs[1] = mxDuplicateArray(w_mx);
    double *w_out = mxGetPr(plhs[1]);

    // === Loop over contact points ===
    for (mwSize i = 0; i < Np; i++) {
        double xn  = x[3*i + 2];
        double xt1 = x[3*i + 0];
        double xt2 = x[3*i + 1];

        // Normal force
        double FN = NormalForces(xn, kn[i], xn0[i]);

        // mexPrintf("i = %d, w1 = %f, w2 = %f\n", i, w_out[0 + i*2], w_out[1 + i*2]);
        // Tangential forces
        double T1, T2;
        TangentialForces(xt1, w_out[0 + i*2], kt[0 + i*2], mu[0 + i*2], FN, T1);
        TangentialForces(xt2, w_out[1 + i*2], kt[1 + i*2], mu[1 + i*2], FN, T2);
        // std::cout <<  "i: " << i << ", w1: " << w_out[0 + i*2] << ", w2: " << w_out[1 + i*2] << std::endl;
        
       
        // Save
        F[3*i + 0] = T1;
        F[3*i + 1] = T2;
        F[3*i + 2] = FN;
    }
}