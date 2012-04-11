/*
  This file is part of LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008-2012  Schalk, Alexander Stippler.

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

#include <lawa/math/math.h>

namespace lawa {

using namespace flens;

template <typename T>
Wavelet<T,Dual,R,CDF>::Wavelet(int _d, int _d_)
    : d(_d), d_(_d_), mu(d&1), l1_((2-(d+d_))/2), l2_((d+d_)/2),
      b_(mask(d,d_)), phi(d), phi_(d,d_)
{
    assert(d<=d_);
    assert(((d+d_)&1)==0);
}

template <typename T>
Wavelet<T,Dual,R,CDF>::Wavelet(const BSpline<T,Primal,R,CDF> &_phi,
                               const BSpline<T,Dual,R,CDF> &_phi_)
    : d(_phi.d), d_(_phi_.d_), mu(d&1), l1_((2-(d+d_))/2), l2_((d+d_)/2),
      b_(mask(d,d_)), phi(_phi), phi_(_phi_)
{
}

template <typename T>
T
Wavelet<T,Dual,R,CDF>::operator()(T x, int j, long k, unsigned short deriv) const
{
    assert(deriv==0);
    
    T ret = T(0);
    x = pow2i<T>(j)*x-k;
    for (int i=b_.firstIndex(); i<=b_.lastIndex(); ++i) {
        ret += b_(i)*phi_(2*x-i, 0, 0);
    }
    return pow2ih<T>(j) * ret;
}

template <typename T>
Support<T>
Wavelet<T,Dual,R,CDF>::support(int j, long k) const
{
    return pow2i<T>(-j) * Support<T>(l1_+k, l2_+k);
}

template <typename T>
const DenseVector<T> &
Wavelet<T,Dual,R,CDF>::mask() const
{
    return b_;
}

template <typename T>
DenseVector<T>
Wavelet<T,Dual,R,CDF>::mask(int d, int d_)
{
    assert(d<=d_);
    assert(((d+d_)&1)==0);

    int mu = d & 1;
    BSpline<T,Primal,R,CDF>  phi(d);
    DenseVector<T> b_(_(1-(d+mu)/2, 1+(d-mu)/2));
    for (int k=b_.firstIndex(); k<=b_.lastIndex(); ++k) {
        int sign = (k&1) ? -1 : 1;
        b_(k) = sign * phi.a(1-k);
    }
    return b_;
}

} // namespace lawa

