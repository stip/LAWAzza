#include <fstream>
#include <iostream>

#include <lawa/lawa.h>
#include <lawa/lawa.tcc>

using namespace lawa;
using namespace std;

int
main()
{

    WaveletCoeffCoord<double>   CoeffCoord;
    
    CoeffCoord(-1, 2) += 1;
    CoeffCoord(-1, 2) += 2;
    CoeffCoord( 1, 2) += 3;
    CoeffCoord( 0, 1) += 4;
    CoeffCoord( 1, 2) += 5;
    CoeffCoord( 0, 1) += 6;
    CoeffCoord(-1, 2) += 1;
    CoeffCoord(-1, 2) += 2;
    CoeffCoord( 1, 2) += 3;
    CoeffCoord( 0, 1) += 4;
    CoeffCoord( 1, 2) += 5;
    CoeffCoord( 0, 1) += 6;
    CoeffCoord(-1, 2) += 1;
    CoeffCoord(-1, 2) += 2;
    CoeffCoord( 1, 2) += 3;
    CoeffCoord( 0, 1) += 4;
    CoeffCoord( 1, 2) += 5;
    CoeffCoord( 0, 1) += 6;
    CoeffCoord(-1, 2) += 1;
    CoeffCoord(-1, 2) += 2;
    CoeffCoord( 1, 2) += 3;
    CoeffCoord( 0, 1) += 4;
    CoeffCoord( 0, 2) += 4;
    CoeffCoord( 0, 4) += 4;
    CoeffCoord( 1, 3) += 5;
    CoeffCoord( 0, 1) += 6;

    CoeffCoord.accumulate();

    auto coordVector = CoeffCoord.coordVector();

    for (int i=0; i<CoeffCoord.numNonZeros(); ++i) {
        auto coord = coordVector[i];
        cout << "j = " << coord.j
             << ", k= " << coord.k
             << ", value = " << coord.value
             << endl;
    }
    

    WaveletCoeff<double>   Coeff = CoeffCoord;

    cout << "firstIndex.range(): " << Coeff._FirstIndex.range() << endl;
    cout << "firstIndex: " << Coeff._FirstIndex << endl;
    cout << "K: "          << Coeff._K << endl;
    cout << "values: "     << Coeff._values << endl;

    return 0;
}
