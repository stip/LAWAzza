#include <fstream>
#include <iostream>

#include <lawa/lawa.h>
#include <lawa/lawa.tcc>

using namespace lawa;
using namespace std;

int
main()
{
    const int          j     = 3;
    const unsigned int d     = 2;
    const unsigned int d_    = 4;
    const unsigned int deriv = 0;


///
/// Setup a B-Spline of order `d`
///
    MRA<double,Primal,R,CDF> mra(d);

    DenseVector<double>  coef(3, -1);
    coef(-1) = 0;
    coef( 0) = 1;
    coef( 1) = 0;

    const int    N = 2000;
    const double a = -4;
    const double b = 4;

    for (int i=0; i<N; ++i) {
        double x = a + i*(b-a)/N;
        double y = evaluate(mra, x, j, coef, deriv);
        cout << x << " " << y << endl;
    }


/*

    Basis<double,Primal,R,CDF> basis(d, d_);

    const int j0 = 0;

    DenseVector<int>     w(4, j0-1);
    DenseVector<double>  coef(7);
    DenseVector<int>     k(7);
    
    w(j0-1) = 1;         // start of scaling coefficients level j0
    w(  j0) = 4;         // start of wavlet coefficients level j0
    w(j0+1) = 6;         // start of wavlet coefficients level j0+1
    w(j0+2) = 8;         // start of wavlet coefficients level j0+2

    coef =  1, 2, 1,     // scaling function on level j0
            4, 4,        // wavelets on level j0
            3, 3;        // wavelets on level j0+1

    k     = -1, 0, 1,     // scaling function on level j0
            -2, 0,        // wavelets on level j0
            -1, 0;        // wavelets on level j0+1

    const int    N = 1;
    const double a = -4;
    const double b = 4;

    for (int i=0; i<N; ++i) {
        double x = a + i*(b-a)/N;
        double y = evaluate(basis, x, j0, j0+2, w, k, coef, deriv);

//        cout << x << " " << y << endl;
    }

*/

    return 0;
}
