#include "mex.h"
#include "functionMEX.h"


using namespace std;
double NormalForcesMEX(double xn, double kn, double xn0)
{
    double fn, u;
    u = xn - xn0;
    fn = max(0.0, kn * u);
    
    return fn;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Check for proper number of arguments.
    if (nrhs != 3) {
        mexErrMsgIdAndTxt("MATLAB:NormalForcesMEX:invalidNumInputs", "Three inputs required.");
    }
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("MATLAB:NormalForcesMEX:invalidNumOutputs", "One output required.");
    }
    
    // Get the input values
    double xn = mxGetScalar(prhs[0]);
    double kn = mxGetScalar(prhs[1]);
    double xn0 = mxGetScalar(prhs[2]);
    
    // Call the NormalForcesMEX function
    double fn = NormalForcesMEX(xn, kn, xn0);
    
    // Set the output value
    plhs[0] = mxCreateDoubleScalar(fn);
}