#include <iostream>
#include <fstream>

#include <lawa/lawa.h>
#include <lawa/lawa.tcc>

using namespace lawa;
using namespace std;

int
main()
{
    const int j = 0;
    const int k = 0;
    const unsigned int d = 4;
    const unsigned int deriv = 0;

    ofstream data("test_bspline_refinement.data");

    BSpline<double,Primal,R,CDF> phi(d);
    Support<double> supp = phi.support(j,k);

    cout << "supp.l1 = " << supp.l1 << endl;
    cout << "supp.l2 = " << supp.l2 << endl;


    cout << "phi.mask().range() = " << phi.mask().range() << endl;
    cout << "phi.mask() = " << phi.mask() << endl;

    int l1 = supp.l1;
    int l2 = supp.l2;

    cout << "l1 = " << l1 << endl;
    cout << "l2 = " << l2 << endl;

    lawa::DenseVector<double>  values(_(l1,l2));

/*
    for (double x=supp.l1; x<=supp.l2; x+=1) {
        double y = phi(x,j,k,deriv);
        
        cout << "x = " << x << ", y = " << y << endl;
        values(int(x)) = y;
        
        data << x << " " << y << endl;
    }
    data << endl << endl;

    cout << endl;
    supp = phi.support(j-1,k);
    for (double x=supp.l1; x<=supp.l2; x+=1) {
        double y = phi(x,j-1,k,deriv);
        cout << "x = " << x << ", y = " << y << endl;

        data << x << " " << y << endl;
    }
    data << endl << endl;
    data.close();
*/

    MRA<double,Primal,R,CDF>    mra(d, j);
    lawa::DenseVector<double>   x(_(-4,4));

    for (int i=-3; i<=3; ++i) {
        lawa::DenseVector<double> y;
        x = 0;
        x(i) = 1;

        mv(NoTrans, double(1), mra.M0, x, double(0), y);

        y *= sqrt(2);
        cout << "y = " << y << endl;
    }


    for (int i=-3; i<=3; ++i) {
        lawa::DenseVector<double> y;
        x = 0;
        x(i) = 1;

        mv(Trans, double(1), mra.M0, x, double(0), y);

        y *= sqrt(2);
        cout << "y = " << y << endl;
    }



//
//  For plotting we create a gnuplot script.
//
    ofstream gps("test_bspline_refinement.gps");
    gps << "set terminal pngcairo;" << endl;
    gps << "set output 'test_bspline_refinement.png'" << endl;
    gps << "set size ratio -1" << endl;
    gps << "plot 'test_bspline_refinement.data' i 0 w l t''";
    for (int i=2; i<=2; ++i) {
        gps << ", 'test_bspline_refinement.data' i " << i-1 << " w l t''";
    }
    gps << endl;
    gps.close();


    return 0;
}
