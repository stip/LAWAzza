#include <iostream>

#include <lawa/lawa.h>
#include <lawa/lawa.tcc>

using namespace lawa;
using namespace std;

int
main()
{
    BSpline<double,Primal,R,CDF> phi(3);

    for (double x=phi.support(0,0).l1; x<=phi.support(0,0).l2; x+=0.0125) {
        cout << x << " " << phi(x,0,0,0) << endl;
    }
    cout << endl;

    return 0;
}
