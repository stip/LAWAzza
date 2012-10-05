#include <fstream>
#include <iostream>

#include <lawa/lawa.h>
#include <lawa/lawa.tcc>

using namespace lawa;
using namespace std;

typedef lawa::GeMatrix<double>  RealGeMatrix;

int
main()
{
    
    //const int          J     = 1;
    const unsigned int d     = 2;
    const unsigned int d_    = 4;
    //const unsigned int deriv = 0;

    typedef Basis<double,Primal,Interval,Dijkema>  PrimalBasis;
    typedef Basis<double,Dual,Interval,Dijkema>    DualBasis;

    PrimalBasis primal(d, d_);
    DualBasis   dual(d, d_);

    primal.enforceBoundaryConditions<NoBC, DirichletBC>();
    dual.enforceBoundaryConditions<NoBC, DirichletBC>();

    // primal.enforceBoundaryConditions<DirichletBC, DirichletBC>();
    // dual.enforceBoundaryConditions<DirichletBC, DirichletBC>();


    const int j0 = std::max(primal.j0, dual.j0);

    //dual.mra_.setLevel(j0+1);


    const int k0 = dual.rangeJ_(j0+1).firstIndex();
    const int k1 = dual.rangeJ_(j0+1).lastIndex();
    
    for (int k=k0; k<=k1; ++k) {
        cout << "M1_(" << j0+1 << ", _, " << k << ") = "
             << dual.M1_(j0+1,_,k)
             << endl;
    }

    RealGeMatrix M0, M1, M0_, M1_;

    cout << "primal.mra.M0.rows() = " << primal.mra.M0.rows() << endl;
    cout << "primal.mra.M0.cols() = " << primal.mra.M0.cols() << endl;

    densify(NoTrans, primal.mra.M0, M0);
    cout << "M0 = " << M0 << endl;


    cout << "primal.M1.rows() = " << primal.M1.rows() << endl;
    cout << "primal.M1.cols() = " << primal.M1.cols() << endl;

    densify(NoTrans, primal.M1, M1);
    cout << "M1 = " << M1 << endl;


    cout << "dual.mra_.M0_.rows() = " << dual.mra_.M0_.rows() << endl;
    cout << "dual.mra_.M0_.cols() = " << dual.mra_.M0_.cols() << endl;

    densify(NoTrans, dual.mra_.M0_, M0_);
    cout << "M0_ = " << M0_ << endl;



    cout << "dual.M1_.rows() = " << dual.M1_.rows() << endl;
    cout << "dual.M1_.cols() = " << dual.M1_.cols() << endl;

    densify(NoTrans, dual.M1_, M1_);
    cout << "M1_ = " << M1_ << endl;


    RealGeMatrix  A, B, C, D;
    
    
    A = transpose(M0)*M0_;
    cerr << "1  A = " << A << endl;

    B = transpose(M0)*M1_;
    cerr << "0  B = " << B << endl;

    C = transpose(M1)*M0_;
    cerr << "0  C = " << C << endl;

    D = transpose(M1)*M1_;
    cerr << "1  D = " << D << endl;


/*
    cout << "dual.M1_.rows() = " << dual.M1_.rows() << endl;
    cout << "dual.M1_.cols() = " << dual.M1_.cols() << endl;

    M1_.resize(0,0);
    densify(NoTrans, dual.M1_, M1_);
    cout << "M1_ = " << M1_ << endl;
*/

    Integral<Trapezoidal,PrimalBasis,DualBasis>  integral(primal, dual);
    integral.quadrature.setN(15000);

    {
        cerr << "case1" << endl;
        const int k0 = primal.mra.rangeI(j0).firstIndex();
        const int k1 = primal.mra.rangeI(j0).lastIndex();

        const int K0 = dual.mra_.rangeI_(j0).firstIndex();
        const int K1 = dual.mra_.rangeI_(j0).lastIndex();


        RealGeMatrix   A(_(k0, k1), _(K0, K1));

        for (int k=k0; k<=k1; ++k) {
            for (int K=K0; K<=K1; ++K) {
                double value = integral(j0, k, XBSpline, 0,
                                        j0, K, XBSpline, 0);
                A(k, K) = value;
            }
        }
        cout << "A(XBSpline, XBSpline) = " << A << endl;
    }

/*
    {
        cerr << "case2" << endl;
        const int k0 = primal.mra.rangeI(j0).firstIndex();
        const int k1 = primal.mra.rangeI(j0).lastIndex();

        const int K0 = dual.rangeJ_(j0).firstIndex();
        const int K1 = dual.rangeJ_(j0).lastIndex();


        RealGeMatrix   A(_(k0, k1), _(K0, K1));

        for (int k=k0; k<=k1; ++k) {
            for (int K=K0; K<=K1; ++K) {
                double value = integral(j0, k, XBSpline, 0,
                                        j0, K, XWavelet, 0);
                A(k, K) = value;
            }
        }
        cout << "A(XBSpline, XWavelet) = " << A << endl;
    }
*/

    {
        cerr << "case3" << endl;
        const int k0 = primal.rangeJ(j0).firstIndex();
        const int k1 = primal.rangeJ(j0).lastIndex();

        const int K0 = dual.mra_.rangeI_(j0).firstIndex();
        const int K1 = dual.mra_.rangeI_(j0).lastIndex();


        RealGeMatrix   A(_(k0, k1), _(K0, K1));

        for (int k=k0; k<=k1; ++k) {
            for (int K=K0; K<=K1; ++K) {
                double value = integral(j0, k, XWavelet, 0,
                                        j0, K, XBSpline, 0);
                A(k, K) = value;
            }
        }
        cout << "A(XWavelet, XBSpline) = " << A << endl;
    }

/*
    {
        cerr << "case4" << endl;
        const int k0 = primal.rangeJ(j0).firstIndex();
        const int k1 = primal.rangeJ(j0).lastIndex();

        const int K0 = dual.rangeJ_(j0).firstIndex();
        const int K1 = dual.rangeJ_(j0).lastIndex();


        RealGeMatrix   A(_(k0, k1), _(K0, K1));

        for (int k=k0; k<=k1; ++k) {
            for (int K=K0; K<=K1; ++K) {
                double value = integral(j0, k, XWavelet, 0,
                                        j0, K, XWavelet, 0);
                A(k, K) = value;
            }
        }
        cout << "A(XWavelet, XWavelet) = " << A << endl;
    }
*/
    return 0;
}
