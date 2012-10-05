#include <fstream>
#include <iostream>

#include <lawa/lawa.h>
#include <lawa/lawa.tcc>

using namespace lawa;
using namespace std;

int
main()
{
    const int d  = 2;
    const int d_ = 6;

    MRA<double,Primal,Interval,Dijkema>  mra(d);
    //mra.enforceBoundaryConditions<DirichletBC, NoBC>();

    MRA<double,Dual,Interval,Dijkema>  mra_(d, d_);
    
#   ifdef ENFORCE_BC
    std::cerr << "enforceBoundaryConditions:" << std::endl;
    //mra_.enforceBoundaryConditions<DirichletBC, NoBC>();
    mra_.enforceBoundaryConditions<NoBC, DirichletBC>();
#   endif // ENFORCE_BC


    {
        const int m = mra.M0.numRows();
        const int n = mra.M0.numCols();
        lawa::GeMatrix<double> M0(m, n);
        lawa::DenseVector<double> e(n);
        for (int i=1; i<=n; ++i) {
            lawa::DenseVector<double> y(m);
            e(i) = 1.;
            M0(_,i) = mra.M0*e;
            e(i) = 0.;
        }
        cout << "M0.rows() = " << M0.rows() << endl;
        cout << "M0.cols() = " << M0.cols() << endl;
        cout << "M0 = " << M0 << endl;
    }


    {
        const int m = mra_.M0_.numRows();
        const int n = mra_.M0_.numCols();
        lawa::GeMatrix<double> M0_(m, n);
        lawa::DenseVector<double> e(n);
        for (int i=1; i<=n; ++i) {
            lawa::DenseVector<double> y(m);
            e(i) = 1.;
            cerr << "i = " << i << endl;
            M0_(_,i) = mra_.M0_*e;
            e(i) = 0.;
        }
        cout << "M0_.rows() = " << M0_.rows() << endl;
        cout << "M0_.cols() = " << M0_.cols() << endl;
        cout << "M0_ = " << M0_ << endl;
    }

#   ifdef ENFORCE_BC
    std::cerr << "compiled with ENFORCE_BC" << std::endl;
#   else
    std::cerr << "compiled without ENFORCE_BC" << std::endl;
#   endif // ENFORCE_BC

}
