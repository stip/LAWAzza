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
    const unsigned int d     = 4;
    const unsigned int d_    = 8;
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


    WaveletCoeffCoord<double>       CoeffCoord;

    const int   K = 1;

    for (int k=-20; k<=20; ++k) {
        CoeffCoord(j0-1,k) += integral(J, K, XBSpline, deriv,
                                       0, k, XBSpline, deriv);
    }

    for (int j=j0; j<=J+5; ++j) {
        for (int k=-20; k<=20; ++k) {
            CoeffCoord(j,k) += integral(J, K, XBSpline, deriv,
                                        j, k, XWavelet, deriv);
        }
    }

    WaveletCoeff<double>  Coeff = CoeffCoord;

    for (int i=0; i<N; ++i) {
        double x = a + i*(b-a)/N;
        double y = evaluate(primal, x, Coeff, deriv);

        cout << x << " " << y << endl;
    }
    
    cerr << "firstIndex.range(): " << Coeff._FirstIndex.range() << endl;
    cerr << "firstIndex: " << Coeff._FirstIndex << endl;
    cerr << "K: "          << Coeff._K << endl;
    cerr << "values: "     << Coeff._values << endl;
    

    return 0;
}
