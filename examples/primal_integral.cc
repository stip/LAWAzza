#include <fstream>
#include <iostream>

#include <lawa/lawa.h>
#include <lawa/lawa.tcc>

using namespace lawa;
using namespace std;

int
main()
{
    const int          j0    = 0;
    const int          J     = 1;
    const unsigned int d     = 2;
    const unsigned int d_    = 4;
    const unsigned int deriv = 0;

///
/// Setup a B-Spline of order `d`
///
    MRA<double,Primal,R,CDF> mra(d);

    DenseVector<double>  c(3, -1);
    c(-1) = 0;
    c( 0) = 0;
    c( 1) = 1;

    const int    N = 400;
    const double a = -4;
    const double b = 4;

    for (int i=0; i<N; ++i) {
        double x = a + i*(b-a)/N;
        double y = evaluate(mra, x, J, c, deriv);
        cout << x << " " << y << endl;
    }
    cout << endl;


    typedef Basis<double,Primal,R,CDF> PrimalBasis;
    typedef Basis<double,Dual,R,CDF>   DualBasis;

    PrimalBasis primal(d,d_);
    DualBasis   dual(d,d_);

    Integral<Trapezoidal,PrimalBasis,DualBasis> integral(primal, dual);
    
    integral.quadrature.setN(1500);


    DenseVector<int>     w(4, j0-1);
    DenseVector<double>  coef(21*3);
    DenseVector<int>     k_(21*3);


    w(j0-1) = 1;         // start of scaling coefficients level j0
    w(  j0) = 1 +   21;  // start of wavlet coefficients level j0
    w(j0+1) = 1 + 2*21;  // start of wavlet coefficients level j0+1
    w(j0+2) = 1 + 2*21;  // start of wavlet coefficients level j0+2




    const int   K = 1;

    int count = 1;

    for (int k=-10; k<=10; ++k) {
        coef(count) = integral(J, K, XBSpline, deriv,
                               0, k, XBSpline, deriv);
        k_(count) = k;
        ++count;
    }

    for (int j=0; j<=J; ++j) {
        for (int k=-10; k<=10; ++k) {
            coef(count) = integral(J, K, XBSpline, deriv,
                                   j, k, XWavelet, deriv);
            k_(count) = k;
            ++count;
        }
    }

    for (int i=0; i<N; ++i) {
        double x = a + i*(b-a)/N;
        double y = evaluate(primal, x, j0, j0+2, w, k_, coef, deriv);

        cout << x << " " << y << endl;
    }

/*
    for (int k1=-10; k1<=10; ++k1) {
        for (int k2=-10; k2<=10; ++k2) {
            cout << integral(0, k1, XBSpline, deriv,
                             0, k2, XWavelet, deriv)
                 << ",  ";
        }
        cout << endl;
    }
*/

    return 0;
}
