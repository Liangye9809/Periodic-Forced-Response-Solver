#include "mex.h"
#include <algorithm>
#include <cmath>
#include <vector>

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
    mexErrMsgIdAndTxt("g_mex:inputType", "fc must be a struct or object.");
    return nullptr;
}

// Main gf calculation (for one row x)
void gf_core(const double *x, mwSize Np,
             double *kn, double *xn0, double *mu, double *kt,
             double *w_inout, double *Fout)
{
    // for debug
    // mexPrintf("In gf_core, xt:\n");
    for (mwSize i = 0; i < Np; i++) {
        double xn  = x[3*i + 2];
        double xt1 = x[3*i + 0];
        double xt2 = x[3*i + 1];
        // mexPrintf("%f, %f, %f, ", i, xt1, xt2, xn);
        double FN = NormalForces(xn, kn[i], xn0[i]);

        double T1, T2;
        TangentialForces(xt1, w_inout[0 + i*2], kt[0 + i*2], mu[0 + i*2], FN, T1);
        TangentialForces(xt2, w_inout[1 + i*2], kt[1 + i*2], mu[1 + i*2], FN, T2);

        Fout[3*i + 0] = T1;
        Fout[3*i + 1] = T2;
        Fout[3*i + 2] = FN;
    }
    // mexPrintf("\n");
}

// Entry: F = g_mex(xt, pc)
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    if (nrhs != 2) {
        mexErrMsgIdAndTxt("g_mex:nrhs", "Two inputs required: xt, fc.");
    }
    if (nlhs != 2) {
        mexErrMsgIdAndTxt("g_mex:nlhs", "Two output required: F, w.");
    }

    // === Input xt ===
    const mxArray *xt_mx = prhs[0];
    if (!mxIsDouble(xt_mx) || mxIsComplex(xt_mx)) {
        mexErrMsgIdAndTxt("g_mex:xtType", "xt must be a real double matrix.");
    }
    double *xt = mxGetPr(xt_mx);
    mwSize N  = mxGetM(xt_mx); // rows
    mwSize Nx = mxGetN(xt_mx); // cols

    /* // for debug
    mexPrintf("before reorder, xt:\n");
    for (mwSize i = 0; i < N; i++) {
        for (mwSize j = 0; j < Nx; j++) {
            mexPrintf("%f ", xt[j + i*Nx]);
        }
        mexPrintf("\n");
    }
    // redoer xt to row-major
    std::vector<double> xt_row_major(N*Nx);
    for (mwSize i = 0; i < N; i++) {
        for (mwSize j = 0; j < Nx; j++) {
            xt_row_major[i*Nx + j] = xt[j*N + i];
        }
    }
    mexPrintf("after reorder, xt:\n");
    for (mwSize i = 0; i < N; i++) {
        for (mwSize j = 0; j < Nx; j++) {
            mexPrintf("%f ", xt_row_major[j + i*Nx]);
        }
        mexPrintf("\n");
    } */

    // reorder to row-major
    // === Input fc ===
    const mxArray *fc = prhs[1];

    const mxArray *kn_mx  = getFieldOrProperty(fc, "kn");
    const mxArray *xn0_mx = getFieldOrProperty(fc, "xn0");
    const mxArray *mu_mx  = getFieldOrProperty(fc, "mu");
    const mxArray *kt_mx  = getFieldOrProperty(fc, "kt");
    const mxArray *w_mx   = getFieldOrProperty(fc, "w");

    if (!kn_mx || !xn0_mx || !mu_mx || !kt_mx || !w_mx) {
        mexErrMsgIdAndTxt("g_mex:missingField",
                          "fc must contain fields/properties kn, xn0, mu, kt, w.");
    }

    double *kn  = mxGetPr(kn_mx);
    double *xn0 = mxGetPr(xn0_mx);
    double *mu  = mxGetPr(mu_mx);
    double *kt  = mxGetPr(kt_mx);
    double *w   = mxGetPr(w_mx);

    mwSize Np = mxGetNumberOfElements(kn_mx);
    if (Nx != 3*Np) {
        mexErrMsgIdAndTxt("g_mex:dimMismatch", "xt must have 3*Np columns.");
    }

    /* // === Outputs: create struct with fields F, w ===
    const char *field_names[] = {"F","w"};
    plhs[0] = mxCreateStructMatrix(1,1,2,field_names);

    // Allocate F.F (N x Nx)
    mxArray *F_mx = mxCreateDoubleMatrix(N, Nx, mxREAL);
    double *F_data = mxGetPr(F_mx);

    // Copy w (will be updated)
    mxArray *w_mx_out = mxDuplicateArray(w_mx);
    double *w_data = mxGetPr(w_mx_out); */

    // === Output F ===
    plhs[0] = mxCreateDoubleMatrix(N, Nx, mxREAL);
    double *F = mxGetPr(plhs[0]);

    // === Output w ===
    plhs[1] = mxDuplicateArray(w_mx);
    double *w_out = mxGetPr(plhs[1]);


    // === Loop: j = 1:2 ===
    std::vector<double> Frow(Nx);
    for (int j = 0; j < 2; j++) {
        for (mwSize i = 0; i < N; i++) {
            // 1. copy row i of xt
            std::vector<double> xrow(Nx);
            for (mwSize j = 0; j < Nx; j++) {
                xrow[j] = xt[i + j * N];   // MATLAB (i+1, j+1)
            }
            // 2. call gf_core
            gf_core(xrow.data(), Np, kn, xn0, mu, kt, w_out, Frow.data());
            // 3. Store row into F
            for (mwSize k = 0; k < Nx; k++) {
                F[i + k*N] = Frow[k]; // column-major
            }
        }
    }

    /* // Assign outputs
    mxSetField(plhs[0], 0, "F", F_mx);
    mxSetField(plhs[0], 0, "w", w_mx_out);
     */
}
