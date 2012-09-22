#include <fstream>
#include <iostream>

#include <lawa/lawa.h>
#include <lawa/lawa.tcc>

using namespace lawa;
using namespace std;

int
main()
{
    const int d = 5;

    MRA<double,Primal,Interval,Dijkema>  mra(d);
    mra.enforceBoundaryConditions<DirichletBC, DirichletBC>();

    const int N = 1000;
    const int j = mra.min_j0;

    const int k0 = mra.rangeI(j).firstIndex();
    const int k1 = mra.rangeI(j).lastIndex();

    for (int k=k0; k<=k1; ++k) {
        for (int i=0; i<=N; ++i) {
            const double x = double(i)/N;
            const double y = mra.phi(x, j, k, 0);
            cout << x << " " << y << endl;
        }
        cout << endl;
    }
    return 0;
}
