#include <fstream>
#include <iostream>

#include <lawa/lawa.h>
#include <lawa/lawa.tcc>

using namespace lawa;
using namespace std;

int
main()
{
    const int          j     = 0;
    const unsigned int d     = 2;
    const unsigned int d_    = 8;


    typedef Basis<double,Primal,R,CDF> PrimalBasis;
    typedef Basis<double,Dual,R,CDF>   DualBasis;

    PrimalBasis primal(d,d_);
    DualBasis   dual(d,d_);


    DenseVector<double>  coef(3, -1);
    coef(-1) = 0;
    coef( 0) = 1;
    coef( 1) = 0;

    const int    N = 200;
    const double a = -4;
    const double b = 4;

    /*
    for (int i=0; i<N; ++i) {
        double x = a + i*(b-a)/N;
        double y = primal.mra.phi(x, j, 0, 0);
        cout << x << " " << y << endl;
    }
    cout << endl;
    */

    /*
    for (int i=0; i<N; ++i) {
        double x = a + i*(b-a)/N;
        double y = primal.psi(x, j, 0, 0);
        cout << x << " " << y << endl;
    }
    cout << endl;
    */

    /*
    for (int i=0; i<N; ++i) {
        double x = a + i*(b-a)/N;
        double y = dual.mra_.phi_(x, j, 0);
        cout << x << " " << y << endl;
    }
    cout << endl;
    */

    for (int i=0; i<N; ++i) {
        double x = a + i*(b-a)/N;
        double y1 = primal.mra.phi(x, j, 0, 0);
        double y2 = primal.psi(x, j, 0, 0);
        double y3 = dual.mra_.phi_(x, j, 0);
        double y4 = dual.psi_(x, j, 0);
        cout << x << " " << y1 << " " << y2 << " " << y3 << " " << y4 << endl;
    }


    return 0;
}
