/*
  This file is part of LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008-2011  Mario Rometsch, Alexander Stippler.

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

#ifndef LAWA_CONSTRUCTIONS_INTERVAL_DIJKEMA_DUAL_WAVELET_TCC
#define LAWA_CONSTRUCTIONS_INTERVAL_DIJKEMA_DUAL_WAVELET_TCC 1

#include <cassert>
#include <lawa/constructions/interval/dijkema/primal/wavelet.h>

namespace lawa {

template <typename T>
Wavelet<T,Dual,Interval,Dijkema>::Wavelet(const Basis<T,Dual,
                                                      Interval,
                                                      Dijkema> &_basis_)
    : basis_(_basis_)
{
}

template <typename T>
const T
Wavelet<T,Dual,Interval,Dijkema>::operator()(T x, int j, Integer k,
                                             unsigned short deriv) const
{
    assert(deriv==0);
    assert(x>=0.);
    assert(x<=1.);
    assert(j>=basis_.min_j0);
    assert(k>=basis_.rangeJ_(j).firstIndex());
    assert(k<=basis_.rangeJ_(j).lastIndex());
    
    const typename DenseVector<T>::ConstView coeffs = basis_.M1_(j,_,k);
    T ret = 0;
    for (int r=coeffs.firstIndex(); r<=coeffs.lastIndex(); ++r) {
        ret += coeffs(r) * basis_.mra_.phi_(x,j+1,r);
    }
    return ret;
}


template <typename T>
const Support<T>
Wavelet<T,Dual,Interval,Dijkema>::support(int j, int k) const
{
    /*
    assert(j>=basis_.min_j0);
    assert(k>=basis_.rangeJ_(j).firstIndex());
    assert(k<=basis_.rangeJ_(j).lastIndex());
    
    if (k<=basis_.M1_.left.lastIndex()) {
        return Support<T>(0.,pow2i<T>(-j-1)*basis_.M1_.lengths(k));
    }
    if (k>pow2i<T>(j)-basis_.M1_.right.length()) {
        return Support<T>(1-pow2i<T>(-j-1)*(basis_.M1_.lengths(k-1-pow2i<T>(j))), 1.);
    }
    return pow2i<T>(-j-1)*Support<T>(basis_.M1_.lengths(0)+1-basis_.d+2*(k-basis_.M1_.left.lastIndex()-1), 
                                     basis_.M1_.lengths(0)+basis_.M1_.leftband.length()+2*(k-basis_.M1_.left.lastIndex()-1));
    */
    return Support<T>(0, 1);
}


} // namespace lawa

#endif // LAWA_CONSTRUCTIONS_INTERVAL_DIJKEMA_DUAL_WAVELET_TCC
