#include <iostream>

#include <lawa/lawa.h>
#include <lawa/lawa.tcc>

using namespace lawa;
using namespace std;

typedef Basis<double,Primal,R,CDF> RealBasis;

int
main()
{
    const unsigned int d  = 4;
    const unsigned int d_ = 8;
    const unsigned int deriv = 1;
    RealBasis basis(d,d_);

    Integral<Gauss,RealBasis,RealBasis> integral(basis,basis);

    int first = -basis.mra.phi.l2+1,
        last  = -basis.mra.phi.l1;
    for (int k=first; k<=last; ++k) { 
        cout << integral(0,k,XBSpline,0,
                         0,0,XBSpline,0) << endl;
    }

    return 0;
}
