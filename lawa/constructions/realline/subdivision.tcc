/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008-2012 Schalk, Alexander Stippler.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */

#include <cassert>
#include <numeric>
#include <lawa/math/math.h>

namespace lawa {

// calculate values at integer nodes by solving an eigenvalue problem (EVP).
template <typename T>
void
_evalAtIntegersByEVP(const DenseVector<T> &a, DenseVector<T> &valuesAtIntegers)
{
    using std::abs;
    
    int l1 = a.firstIndex(),
        l2 = a.lastIndex();

    GeMatrix<T> A(a.length(), a.length(), l1, l1);
    
    // fill matrix A of eigenvalue problem.
    for (int m=l1; m<=l2; ++m) {
        int from = std::max(l1, 2*m-l2),
              to = std::min(l2, 2*m-l1);
        for (int k=from; k<=to; ++k) {
            A(m,k) = a(2*m-k);
        }
    }

    // calculate eigenvalues of A.
    DenseVector<T> wi, wr, work; // eigenvalues
    GeMatrix<T> VL, VR; // left and right eigenvectors
    A.changeIndexBase(1,1);
    // only calculate right eigenvectors.
    flens::lapack::ev(false, true, A, wr, wi, VL, VR, work);

    // we choose the eigenvector correspondig to the eigenvalue closest to 1.

    T minDiff = abs(wr(1) - 1);
    int pos = 1;
    for (int i=2; i<=wr.length(); ++i) {
        if (abs(wr(i)-1)<minDiff) {
            minDiff = abs(wr(i)-1);
            pos = i;
        }
    }
    valuesAtIntegers = VR(_,pos); // choosing the corresponding eigenvector.
    
    // the elements of the eigenvector have to sum up to 1.
    T sum = std::accumulate(valuesAtIntegers.engine().data(),
                            valuesAtIntegers.engine().data()
                          + valuesAtIntegers.length(), 0.0);
    
    assert(sum);
    valuesAtIntegers /= sum;
    valuesAtIntegers.engine().changeIndexBase(l1);
}

template <typename T>
void
_evalAtIntegers(const BSpline<T,Primal,R,CDF> phi,
                DenseVector<T> &valuesAtIntegers)
{
    _evalAtIntegersByEVP(phi.a, valuesAtIntegers);
}

template <typename T>
void
_evalAtIntegers(const BSpline<T,Dual,R,CDF> phi_,
                DenseVector<T> &valuesAtIntegers)
{
    if ((phi_.d==2) && (phi_.d_==2)) {
        valuesAtIntegers.engine().resize(5,-2);
        valuesAtIntegers = 0., -2., 5., -2., 0.;
    } else {
        _evalAtIntegersByEVP(phi_.a_, valuesAtIntegers);
    }
}

template <typename T>
void
subdivide(const BSpline<T,Primal,R,CDF> &phi,
          const int                     j,
          DenseVector<T>                &dyadicValues)
{
    DenseVector<T> valuesAtIntegers;
    _evalAtIntegers(phi, valuesAtIntegers);
    // set values at integer positions of result.
    int twoJ = pow2i<T>(j);
    int from = twoJ*phi.l1,
          to = twoJ*phi.l2;
    dyadicValues.engine().resize(to-from+1,from);
    dyadicValues(_(from,twoJ,to)) = valuesAtIntegers;

    // calculate values inbetween on dyadic grid.
    for (int l=1; l<=j; ++l) {
        for (int k=from+pow2i<T>(j-l); k<=to; k+=pow2i<T>(j-l+1)) {
            int mFrom = std::max(phi.l1, ((2*k-to)>>j)+1);
            int   mTo = std::min(phi.l2,  (2*k-from)>>j);
            for (int m=mFrom; m<=mTo; ++m) {
                dyadicValues(k) += phi.a(m)*dyadicValues(2*k-twoJ*m);
            }
        }
    }
}

template <typename T>
void
subdivide(const BSpline<T,Dual,R,CDF> &phi_, int j,
          DenseVector<T> &dyadicValues)
{
    DenseVector<T> valuesAtIntegers;
    _evalAtIntegers(phi_, valuesAtIntegers);

    // set values at integer positions of result.
    int twoJ = pow2i<T>(j);
    int from = twoJ*phi_.l1_,
          to = twoJ*phi_.l2_;

    dyadicValues.engine().resize(to-from+1,from);
    
    dyadicValues(_(from,twoJ,to)) = valuesAtIntegers;

    // calculate values inbetween on dyadic grid.
    for (int l=1; l<=j; ++l) {
        for (int k=from+pow2i<T>(j-l); k<=to; k+=pow2i<T>(j-l+1)) {
            int mFrom = std::max(phi_.l1_, ((2*k-to)>>j)+1);
            int   mTo = std::min(phi_.l2_,  (2*k-from)>>j);
            for (int m=mFrom; m<=mTo; ++m) {
                dyadicValues(k) += phi_.a_(m)*dyadicValues(2*k-twoJ*m);
            }
        }
    }
}

} // namespace lawa
