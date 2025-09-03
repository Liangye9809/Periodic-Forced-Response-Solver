#pragma once
#include <iostream>
#include <cmath>

double NormalForcesMEX(double xn, double kn, double xn0);

struct TengentialForceAndDisplacement
{
    double ft;
    double w;
};
TengentialForceAndDisplacement TangentialForcesMEX(double xt, double wt, double kt, double mu, double fn);

