#include <iostream>
#include <fstream>

#include <lawa/lawa.h>
#include <lawa/lawa.tcc>

using namespace lawa;
using namespace std;


int
main()
{
    const int j = 1;
    const int k = 0;
    const unsigned int d = 4;
    const unsigned int d_ = 6;

    BSpline<double,Primal,R,CDF> phi(d);
    BSpline<double,Dual,R,CDF>   phi_(d,d_);

    Wavelet<double,Primal,R,CDF> psi(d, d_);
    Wavelet<double,Dual,R,CDF>   psi_(d,d_);

    cout << "phi.mask() = " << phi.mask() << endl;
    cout << "phi_.mask() = " << phi_.mask() << endl;

    lawa::DenseVector<double> c1(_(-25,45));
    c1 = 0;
    c1(0) = 1;

    lawa::DenseVector<double> c0, d0;
    

    flens::RefinementMatrix<double,R,CDF>  M0_(phi_);
    flens::RefinementMatrix<double,R,CDF>  M1_(psi_);

    c0 = transpose(M0_) * c1;
    
    c0(_(-5,2)) = 0;
    
    
    cout << "c1.range() = " << c1.range() << endl;
    cout << "c0.range() = " << c0.range() << endl;
    cout << "c1 = " << c1 << endl;
    cout << "c0 = " << c0 << endl;


    d0 = transpose(M1_) * c1;
    cout << "d0.range() = " << d0.range() << endl;
    d0(_(-5,8)) = 0;

    cout << "c0 = " << c0 << endl;
    cout << "d0 = " << d0 << endl;
    lawa::DenseVector<double> C1, D1;

    flens::RefinementMatrix<double,R,CDF>  M0(phi);
    flens::RefinementMatrix<double,R,CDF>  M1(psi);

    C1 = M0 * c0;
    D1 = M1 * d0;

    cout << endl;

    cout << "C1.range() = " << C1.range() << endl;
    cout << "D1.range() = " << D1.range() << endl;
    cout << "C1 = " << C1(_(-20,-2)) << endl;
    cout << "D1 = " << D1(_(-4,14)) << endl;

    C1(_(-20,-2)) += D1(_(-4,14));

    cout << "C1 = " << C1 << endl;

    return 0;
}
