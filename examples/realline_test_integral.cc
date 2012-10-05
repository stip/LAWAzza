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


typedef lawa::GeMatrix<double>  RealGeMatrix;


int
main()
{
    const int J              = 4;
    const int K              = 10;
    const unsigned int d     = SET_D;
    const unsigned int d_    = SET_D_;
    const unsigned int deriv = 0;

    typedef Basis<double,Primal,R,CDF> PrimalBasis;
    typedef Basis<double,Dual,R,CDF>   DualBasis;

    PrimalBasis primal(d,d_);
    DualBasis   dual(d,d_);

    Integral<Trapezoidal,PrimalBasis,DualBasis>  integral(primal, dual);
    integral.quadrature.setN(1000);

    // phi * phi_
    {
        RealGeMatrix   A(_(-K,K),_(-K,K));

        for (int k1=-K; k1<=K; ++k1) {
            for (int k2=-K; k2<=K; ++k2) {
                double value = integral(J, k1, XBSpline, deriv,
                                        J, k2, XBSpline, 0);
                A(k1,k2) = value;
            }
        }
        cout << "A1 = [" << A << "]" << endl;
    }

    // phi * psi_
    {
        RealGeMatrix   A(_(-K,K),_(-K,K));

        for (int k1=-K; k1<=K; ++k1) {
            for (int k2=-K; k2<=K; ++k2) {
                double value = integral(J, k1, XBSpline, deriv,
                                        J, k2, XWavelet, 0);
                A(k1,k2) = value;
            }
        }
        cout << "A2 = [" << A << "]" << endl;
    }

    // psi * phi_
    {
        RealGeMatrix   A(_(-K,K),_(-K,K));

        for (int k1=-K; k1<=K; ++k1) {
            for (int k2=-K; k2<=K; ++k2) {
                double value = integral(J, k1, XWavelet, deriv,
                                        J, k2, XBSpline, 0);
                A(k1,k2) = value;
            }
        }
        cout << "A3 = [" << A << "]" << endl;
    }

    // psi * psi_
    {
        RealGeMatrix   A(_(-K,K),_(-K,K));

        for (int k1=-K; k1<=K; ++k1) {
            for (int k2=-K; k2<=K; ++k2) {
                double value = integral(J, k1, XWavelet, deriv,
                                        J, k2, XWavelet, 0);
                A(k1,k2) = value;
            }
        }
        cout << "A4 = [" << A << "]" << endl;
    }

    return 0;
}
