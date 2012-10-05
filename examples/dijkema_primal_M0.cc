#include <fstream>
#include <iostream>

#include <lawa/lawa.h>
#include <lawa/lawa.tcc>

using namespace lawa;
using namespace std;

int
main()
{
    const int d = 2;

    MRA<double,Primal,Interval,Dijkema>  mra(d);

    mra.enforceBoundaryConditions<DirichletBC, DirichletBC>();
    
    cout << mra.cardIL() << " " << mra.cardIR() << endl;
    cout << mra.M0.left.range() << " " << mra.M0.right.range() << endl;
//    mra.setLevel(4);
    
    const int m = mra.M0.numRows();
    const int n = mra.M0.numCols();
    lawa::GeMatrix<double> D(m,n);
    lawa::DenseVector<double> e(n);
    for (int i=1; i<=n; ++i) {
        lawa::DenseVector<double> y(m);
        e(i) = 1.;
        D(_,i) = mra.M0*e;
        e(i) = 0.;
    }
    cout << D << endl;
}
