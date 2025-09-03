#include<iostream>
#include<cmath>
#include"functionMEX.h"

using namespace std;





int main()
{
    double fn;
    fn = NormalForcesMEX(-0.1, 1000, 0.05);
    cout << "Normal Force: " << fn << endl;
    
    return 0;
}