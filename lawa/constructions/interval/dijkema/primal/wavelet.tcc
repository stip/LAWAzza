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

#ifndef LAWA_CONSTRUCTIONS_INTERVAL_DIJKEMA_PRIMAL_WAVELET_TCC
#define LAWA_CONSTRUCTIONS_INTERVAL_DIJKEMA_PRIMAL_WAVELET_TCC 1

#include <cassert>
#include <lawa/constructions/interval/dijkema/primal/wavelet.h>

namespace lawa {

template <typename T, Construction Cons>
Wavelet<T,Primal,Interval,Cons>::Wavelet(const Basis<T,Primal,
                                                     Interval,Cons> &_basis)
    : basis(_basis)
{
}

template <typename T, Construction Cons>
const T
Wavelet<T,Primal,Interval,Cons>::operator()(T x, int j, Integer k,
                                            unsigned short deriv) const
{
    assert(x>=0.);
    assert(x<=1.);
    assert(j>=basis.min_j0);
    assert(k>=basis.rangeJ(j).firstIndex());
    assert(k<=basis.rangeJ(j).lastIndex());

    const typename DenseVector<T>::ConstView coeffs = basis.M1(j,_,k);
    T ret = 0;
    for (int r=coeffs.firstIndex(); r<=coeffs.lastIndex(); ++r) {
        ret += coeffs(r) * basis.mra.phi(x,j+1,r,deriv);
    }
    return ret;
}

template <typename T, Construction Cons>
const Support<T>
Wavelet<T,Primal,Interval,Cons>::support(int j, int k) const
{
    /*
    assert(j>=basis.min_j0);
    assert(k>=basis.rangeJ(j).firstIndex());
    assert(k<=basis.rangeJ(j).lastIndex());
    if (k<=basis.M1.left.lastIndex()) {
        return Support<T>(0.,pow2i<T>(-j-1)*basis.M1.lengths(k));
    }
    if (k>pow2i<T>(j)-basis.M1.right.length()) {
        return Support<T>(1-pow2i<T>(-j-1)
                        *(basis.M1.lengths(k-1-pow2i<T>(j))), 1.);
    }
    // FIXME: remove std::max: left support end cannot be less than 0. Check for error (Primbs!!!)
    return pow2i<T>(-j-1)*Support<T>(std::max(0L,basis.M1.lengths(0)+1-basis.d+2*(k-basis.M1.left.lastIndex()-1)),
                                     basis.M1.lengths(0)+basis.M1.leftband.length()+2*(k-basis.M1.left.lastIndex()-1));
    */
    return Support<T>(0, 1);
}

template <typename T, Construction Cons>
const DenseVector<T>
Wavelet<T,Primal,Interval,Cons>::singularSupport(int j, long k) const
{
    assert(j>=basis.min_j0);
    assert(k>=basis.rangeJ(j).firstIndex());
    assert(k<=basis.rangeJ(j).lastIndex());

    if (k<=basis.M1.left.lastIndex()) {
        return linspace(0.,
                        pow2i<T>(-j-1)*basis.M1.lengths(k),
                        2*basis.M1.lengths(k)+1.);
    }
    if (k>pow2i<T>(j)-basis.M1.right.length()) {
        return linspace(1-pow2i<T>(-j-1)*(basis.M1.lengths(k-1-pow2i<T>(j))), 
                        1.,
                        2*basis.M1.lengths(k-1-pow2i<T>(j))+1.);
    }
    // FIXME: remove std::max: left support end cannot be less than 0. Check for error (Primbs!!!)
    return pow2i<T>(-j-1)*linspace(std::max(0.,basis.M1.lengths(0)+1-basis.d+2*(k-basis.M1.left.lastIndex()-1.)),
                                   basis.M1.lengths(0)+basis.M1.leftband.length()+2*(k-basis.M1.left.lastIndex()-1.),
                                   // FIXME: understand why value in last line is too large
                                   2*(basis.d+basis.d_)-1);
                                   // FIXME: understand why 2*n instead of 2*n-1  ... (+d+1)
                                   //2*(basis.M1.leftband.length())+basis.d+1.);
}

template <typename T, Construction Cons>
const int
Wavelet<T,Primal,Interval,Cons>::vanishingMoments(int j, long k) const
{
    assert(j>=basis.min_j0);
    assert(k>=basis.rangeJ(j).firstIndex());
    assert(k<=basis.rangeJ(j).lastIndex());

    assert(0);
    return 0;
}

template <typename T, Construction Cons>
const T
Wavelet<T,Primal,Interval,Cons>::tic(int j) const
{
    return pow2i<T>(-j-1);
}

} // namespace lawa

#endif // LAWA_CONSTRUCTIONS_INTERVAL_DIJKEMA_PRIMAL_WAVELET_TCC
