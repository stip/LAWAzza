#include <fstream>
#include <iostream>

#include <lawa/lawa.h>
#include <lawa/lawa.tcc>

using namespace lawa;
using namespace std;

int
main()
{
    const int j = 1;
    const unsigned int d = 2;
    const unsigned int d_ = 4;

///
/// Setup a B-Spline of order `d`
///
    BSpline<double,Dual,R,CDF> phi_(d,d_);

    ofstream data("dual_bspline.data");

///
/// For `j=0` the support of `phi` is the interval `[phi.l1, phi.l2]`.  We
/// want to plot all splines that overlap with `[0,1]`.  As the support of
/// $\varphi\left(x - k\right)$ is `[phi.l1+k, phi.l2+k]` we select only
/// indices of `k` that satify:
///
///   - `phi.l1+k < 1`, i.e. `k < 1-phi.l1` or
///   - `phi.l2+k > 0`, i.e. `k > -phi.l2`
///
/// So the smallest relevant value of `k` is `-phi.l2+1` and the largest is
/// `-phi.l1`.
///
    int first = -phi_.l2_+1,
        last  = -phi_.l1_;

///
/// Next we evaluate $\varphi\left(x - k\right)$ for all these values of $k$
///
    int k = 0;
    //for (int k=first; k<=last; ++k) {
        Support<double> supp = phi_.support(j,k);
        for (double x=supp.l1; x<=supp.l2; x+=1.0/256) {
            data << x << " " << phi_(x,j,k) << endl;
        }
        data << endl << endl;
    //}
    data << endl;
    data.close();

///
/// For plotting we create a gnuplot script.
///
    ofstream gps("dual_bspline.gps");
    gps << "set terminal pngcairo;" << endl;
    gps << "set output 'dual_bspline.png'" << endl;
    gps << "set size ratio -1" << endl;
    gps << "plot 'dual_bspline.data' i 0 w l t''";
    for (int i=first+1; i<=last; ++i) {
        gps << ", 'dual_bspline.data' i " << i-first << " w l t''";
    }
    gps << endl;
    gps.close();

    return 0;
}
