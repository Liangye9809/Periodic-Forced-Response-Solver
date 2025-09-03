
#include "functionMEX.h"
#include <cmath>
#include "mex.h"
using namespace std;

TengentialForceAndDisplacement TangentialForcesMEX(double xt, double wt, double kt, double mu, double fn)
{
    TengentialForceAndDisplacement result;
    double w, ft;
    if (fn > 0)
    {
        ft = kt * (xt - wt);
        if (std::abs(ft) < mu * fn)
        {
            w = wt;
        }
        else
        {
            double sign = (ft > 0) - (ft < 0); // sign function
            ft = mu * fn * sign;
            w = xt - sign * mu * (fn / kt);
        }
    }
    else
    {
        ft = 0;
        w = wt;
    }
    
    result.ft = ft;
    result.w = w;
    
    return result;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Check for proper number of arguments.
    if (nrhs != 5) {
        mexErrMsgIdAndTxt("MATLAB:TangentialForcesMEX:invalidNumInputs", "Five inputs required.");
    }
    if (nlhs != 2) {
        mexErrMsgIdAndTxt("MATLAB:TangentialForcesMEX:invalidNumOutputs", "Two outputs required.");
    }
    
    // Get the input values
    double xt = mxGetScalar(prhs[0]);
    double wt = mxGetScalar(prhs[1]);
    double kt = mxGetScalar(prhs[2]);
    double mu = mxGetScalar(prhs[3]);
    double fn = mxGetScalar(prhs[4]);
    
    // Call the TangentialForcesMEX function
    TengentialForceAndDisplacement result = TangentialForcesMEX(xt, wt, kt, mu, fn);
    
    // Set the output values
    plhs[0] = mxCreateDoubleScalar(result.ft);
    plhs[1] = mxCreateDoubleScalar(result.w);
}