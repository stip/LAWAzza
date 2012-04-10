/*
 LAWA - Library for Adaptive Wavelet Applications.
 Copyright (C) 2008-2011 Sebastian Kestler, Kristina Steih,
                         Alexander Stippler, Schalk.

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

#ifndef LAWA_CONSTRUCTIONS_REALLINE_PRIMAL_BSPLINE_TCC
#define LAWA_CONSTRUCTIONS_REALLINE_PRIMAL_BSPLINE_TCC 1

#include <cassert>
#include <cmath>
#include <lawa/flensforlawa.h>
#include <lawa/math/math.h>

namespace lawa {

template <typename T>
    const DenseVector<T>
    _bspline_mask(int d);

template <typename T>
BSpline<T,Primal,R,CDF>::BSpline(int _d)
    : d(_d), mu(d&1), l1(.5*(-d+mu)), l2(.5*(d+mu)),
      a(_bspline_mask<T>(d))
{
    assert(_d>0);
}

/*TODO
template <typename T>
BSpline<T,Primal,R,CDF>::BSpline(const MRA<T,Primal,R,CDF> &mra)
    : d(mra.d), mu(mra.d&1), l1(.5*(-d+mu)), l2(.5*(d+mu)), 
      a(_bspline_mask<T>(d))
{
    assert(mra.d>0);
}
*/

template <typename T>
BSpline<T,Primal,R,CDF>::~BSpline()
{
}

template <typename T>
const T
BSpline<T,Primal,R,CDF>::operator()(T x, int j, Integer k, unsigned short deriv) const
{
    using std::abs;
    using std::pow;

    if (inner(x,support(j,k))) {
        T ret = T(0);
        x = pow2i<T>(j)*x-k - mu/2.;
        if (deriv==0) {
            x = abs(x);
            for (int p=0; p<=floor<int>(d/2.-x); ++p) {
                int sign = (p&1) ? -1 : 1;
                ret +=   sign * binomial(d, p) * pow(T(d/2.-x-p), d-1);
            }
            ret /= factorial(d - 1);
        } else {
            for (int p=0; p<=floor<int>(d/2.-abs(x)); ++p) {
                int sign = ( (p&1)==( (x>0)&&(deriv&1) ) ) ? 1 : -1;
                ret += sign * binomial(d, p)
                            * pow(T(d/2.-abs(x)-p), d-1-deriv);
            }
            ret /= factorial(d-1-deriv);
        }
        // 2^(j*deriv) * 2^(j/2) * ret
        return pow2ih<T>(2*j*deriv+j)*ret;
    }
    return T(0);
}

template <typename T>
const Support<T>
BSpline<T,Primal,R,CDF>::support(int j, Integer k) const
{
    return pow2i<T>(-j) * Support<T>(l1 + k, l2 + k);
}

/*TODO
template <typename T>
const DenseVector<Array<T> >
BSpline<T,Primal,R,CDF>::singularSupport(int j, Integer k) const
{
    return linspace(support(j,k).l1, support(j,k).l2, d+1);
}

template <typename T>
const T
BSpline<T,Primal,R,CDF>::tic(int j) const
{
    return pow2i<T>(-j);
}
*/

template <typename T>
const DenseVector<T> &
BSpline<T,Primal,R,CDF>::mask() const
{
    return a;
}

//------------------------------------------------------------------------------

template <typename T>
const DenseVector<T>
_bspline_mask(int d)
{
    assert(d>=0);

    if (d==0) {
        DenseVector<T> res(1);
        res = T(1.);
        return res;
    }
    T factor = 1 << (d-1);
    int kappa = d & 1;

    int from = -(d+kappa)/2;
    int to   =  (d+kappa)/2;
    DenseVector<T> res(_(from,to));
    for (int i=from; i<=to; ++i) {
        res(i) = binomial(d,i-from) / factor;
    }
    return res;
}

} // namespace lawa

#endif // LAWA_CONSTRUCTIONS_REALLINE_PRIMAL_BSPLINE_TCC
