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
    const unsigned int d_ = 2;
    BSpline<double,Dual,R,CDF> phi_(d,d_);

    DenseVector<double> phi_cascade;
    evalAtDyadicGrid_Cascade(phi_.a_, 10, phi_cascade);
    phi_cascade.changeIndexBase(1);

    ofstream data("dual_bspline.data");
    int first = -phi_.l2_+1,
        last  = -phi_.l1_;
    last = first;
    for (int k=first; k<=last; ++k) { 
        Support<double> supp = phi_.support(j,k);
        DenseVector<double> x = linspace(supp.l1, supp.l2, phi_cascade.length());
        x.changeIndexBase(1);
        for (int i=1; i<=x.length(); ++i) {
            data << x(i) << " " << phi_(x(i),j,k) << " "
                 << /*phi_(x(i),j,k)-*/phi_cascade(i) << endl;
        }
        data << endl << endl; 
    }
    data << endl;
    data.close();

    ofstream gps("dual_bspline.gps");
    gps << "set terminal png;" << endl;
    gps << "set output 'dual_bspline" << d << "_" << d_ << ".png'" << endl;
    gps << "plot 'dual_bspline.data' i 0 w l t''";
    for (int i=first+1; i<=last; ++i) {
        gps << ", 'dual_bspline.data' i " << i-first << " w l t''";
    }
    gps << endl;
    gps.close();

    return 0;
}
