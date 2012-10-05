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


typedef GeMatrix<double>        RealGeMatrix;
typedef DenseVector<double>     RealDenseVector;


int
main()
{
    const unsigned int d     = SET_D;
    const unsigned int d_    = SET_D_;

    const int J     = 3;
    const int j0    = 0;


    // coeff on level J
    RealDenseVector   c(_(-2,4));
    c    = 0;
    c(0) = 1;

    typedef Basis<double,Primal,R,CDF> PrimalBasis;
    typedef Basis<double,Dual,R,CDF>   DualBasis;

    PrimalBasis primal(d,d_);
    DualBasis   dual(d,d_);

    Integral<Trapezoidal,PrimalBasis,DualBasis>  integral(primal, dual);
    integral.quadrature.setN(1000);

    WaveletCoeffCoord<double>       CoeffCoord, CoeffCoordCheck;

//
//  Compute the wavelet coefficient with brute force: Integration!
//
    const int k1 = -10;
    const int k2 =  10;
    for (int j=j0-1; j<=J; ++j) {
        for (int k=k1; k<=k2; ++k) {
            double value = 0;
            if (j>=j0) {
                for (int m=c.firstIndex(); m<=c.lastIndex(); ++m) {
                    value += c(m)*integral(J, m, XBSpline, 0,
                                           j, k, XWavelet, 0);
                }
            } else {
                for (int m=c.firstIndex(); m<=c.lastIndex(); ++m) {
                    value += c(m)*integral(J, m, XBSpline, 0,
                                           0, k, XBSpline, 0);
                }
            }
            CoeffCoordCheck(j,k) += value;
        }
    }

//
//  Compute the wavelet coefficient the right way: FWT!!!
//
    for (int j=J-1; j>=j0; --j) {
        RealDenseVector d0, c0;
        c0 = transpose(dual.mra_.M0_) * c;
        d0 = transpose(dual.M1_) * c;

        RealDenseVector  cCheck;
        compose(c0, primal.mra.phi.a, d0, primal.psi.b, cCheck);
        cCheck /= sqrt(2);

        cerr << "c = " << c << endl;
        cerr << "cCheck = " << cCheck << endl << endl;

        cerr << "mask: a = " << primal.mra.phi.a << endl;
        cerr << "mask: b = " << primal.psi.b << endl;

        RealDenseVector  diff = c - cCheck(_(c.firstIndex(), c.lastIndex()));
        cerr << "diff = " << diff << endl;

        for (int i=d0.firstIndex(); i<=d0.lastIndex(); ++i) {
            CoeffCoord(j,i) += d0(i);
        }

        if (j==j0) {
            for (int i=c0.firstIndex(); i<=c0.lastIndex(); ++i) {
                CoeffCoord(j-1,i) += c0(i);
            }
        }

        c.resize(c0);
        c = c0;
    }

    WaveletCoeff<double>  CoeffCheck = CoeffCoordCheck;
    WaveletCoeff<double>  Coeff      = CoeffCoord;

    const int    N = 400;
    const double a = -4;
    const double b = 4;

    for (int i=0; i<N; ++i) {
        double x = a + i*(b-a)/N;
        double y1 = evaluate(primal, x, Coeff, 0);
        double y2 = evaluate(primal, x, CoeffCheck, 0);
        double y3 = primal.mra.phi(x, J, 0, 0);

        cout << x << " " << y1 << " " << y2 << " " << y3 << endl;
    }



    return 0;
}
