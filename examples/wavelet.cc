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
    const int k = 0;
    const unsigned int d = 3;
    const unsigned int d_ = 5;
    const unsigned int deriv = 1;

///
/// Setup a wavelet of order `d,d_`
///
    Wavelet<double,Primal,R,CDF> psi(d,d_);

    ofstream data("wavelet.data");

    Support<double> supp = psi.support(j,k);
    
    for (double x=supp.l1; x<=supp.l2; x+=0.0125) {
        data << x << " " << psi(x,j,k,deriv) << endl;
    }
    data << endl;
    data.close();

///
/// For plotting we create a gnuplot script.
///
    ofstream gps("wavelet.gps");
    gps << "set terminal pngcairo;" << endl;
    gps << "set output 'wavelet.png'" << endl;
    gps << "set size ratio -1" << endl;
    gps << "plot 'wavelet.data' w l t''";
    gps << endl;
    gps.close();

    return 0;
}
