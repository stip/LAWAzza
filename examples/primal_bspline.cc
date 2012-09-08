#include <fstream>
#include <iostream>

#include <lawa/lawa.h>
#include <lawa/lawa.tcc>

using namespace lawa;
using namespace std;

int
main()
{
    const int j = 0;
    const unsigned int d = 2;
    const unsigned int deriv = 0;

///
/// Setup a B-Spline of order `d`
///
    BSpline<double,Primal,R,CDF> phi(d);
    ofstream data("primal_bspline.data");

///
/// For `j=0` the support of `phi` is the interval `[phi.l1, phi.l2]`.  We
/// want to plot all splines that overlap with `[0,1]`.  As the support of
/// $\varphi\left(x - k\right)$ is `[phi.l1+k, phi.l2+k]` we select only
/// indices of `k` that satify:
///
/// - `phi.l1+k < 1`, i.e. `k < 1-phi.l1` or
/// - `phi.l2+k > 0`, i.e. `k > -phi.l2`
///
/// So the smallest relevant value of `k` is `-phi.l2+1` and the largest is
/// `-phi.l1`.
///
    int first = -phi.l2+1,
        last  = -phi.l1;

///
/// Next we evaluate $\varphi\left(x - k\right)$ for all these values of $k$
///
    int k = 0;
    //for (int k=first; k<=last; ++k) {
        Support<double> supp = phi.support(j,k);
        for (double x=supp.l1; x<=supp.l2; x+=1.0/256) {
            data << x << " " << phi(x,j,k,deriv) << endl;
        }
        data << endl << endl;
    //}
    data << endl;
    data.close();

///
///  For plotting we create a gnuplot script.
///
    ofstream gps("primal_bspline.gps");
    gps << "set terminal pngcairo;" << endl;
    gps << "set output 'primal_bspline.png'" << endl;
    gps << "set size ratio -1" << endl;
    gps << "plot 'bspline.data' i 0 w l t''";
    for (int i=first+1; i<=last; ++i) {
        gps << ", 'primal_bspline.data' i " << i-first << " w l t''";
    }
    gps << endl;
    gps.close();


    return 0;
}
