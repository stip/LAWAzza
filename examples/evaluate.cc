#include <fstream>
#include <iostream>

#include <lawa/lawa.h>
#include <lawa/lawa.tcc>

using namespace lawa;
using namespace std;

#ifndef SET_D
#define SET_D 2
#endif

#ifndef SET_D_
#define SET_D_ 4
#endif

int
main()
{
    const int          j     = 0;
    const unsigned int d     = SET_D;
    const unsigned int d_    = SET_D_;


    typedef Basis<double,Primal,R,CDF> PrimalBasis;
    typedef Basis<double,Dual,R,CDF>   DualBasis;

    PrimalBasis primal(d,d_);
    DualBasis   dual(d,d_);


    const double a = -10;
    const double b = 10;
    const int    N = 10000*(b-a);


    double sumY1 = 0;
    double sumY2 = 0;
    double sumY3 = 0;
    double sumY4 = 0;

    for (int i=0; i<N; ++i) {
        double x = a + i*(b-a)/N;
        double y1 = primal.mra.phi(x, j, 0, 0);
        double y2 = primal.psi(x, j, 0, 0);
        double y3 = dual.mra_.phi_(x, j, 0);
        double y4 = dual.psi_(x, j, 0);
        cout << x << " " << y1 << " " << y2 << " " << y3 << " " << y4 << endl;

        sumY1 += y1;
        sumY2 += y2;
        sumY3 += y3;
        sumY4 += y4;
    }

    cerr << "primal.mra.phi.mask().range() = " << primal.mra.phi.mask().range() << endl;
    cerr << "primal.mra.phi.mask() = " << primal.mra.phi.mask() << endl;

    cerr << "dual.mra.phi_.mask().range() = " << dual.mra_.phi_.mask().range() << endl;
    cerr << "dual.mra.phi_.mask() = " << dual.mra_.phi_.mask() << endl;

    cerr << "dual.psi_.mask() = " << dual.psi_.mask() << endl;

    cerr << "sumY1 = " << sumY1 << endl;
    cerr << "sumY2 = " << sumY2 << endl;
    cerr << "sumY3 = " << sumY3 << endl;
    cerr << "sumY4 = " << sumY4 << endl;

    return 0;
}
