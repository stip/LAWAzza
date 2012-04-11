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
    const unsigned int d = 4;
    const unsigned int deriv = 1;
    BSpline<double,Primal,R,CDF> phi(d);
   
    ofstream data("schalk.data");
    int first = -phi.l2+1,
        last  = -phi.l1;
    for (int k=first; k<=last; ++k) { 
        Support<double> supp = phi.support(j,k);
        for (double x=supp.l1; x<=supp.l2; x+=0.0125) {
            data << x << " " << phi(x,j,k,deriv) << endl;
        }
        data << endl << endl; 
    }
    data << endl;
    data.close();

    ofstream gps("schalk.gps");
    gps << "set terminal png;" << endl;
    gps << "set output 'N" << d << "_" << deriv << ".png'" << endl;
    gps << "plot 'schalk.data' i 0 w l t''";
    for (int i=first+1; i<=last; ++i) {
        gps << ", 'schalk.data' i " << i-first << " w l t''";
    }
    gps << endl;
    gps.close();

    return 0;
}
