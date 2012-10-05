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

#ifndef LAWA_CONSTRUCTIONS_INTERVAL_DIJKEMA_PRIMAL_WAVELET_H
#define LAWA_CONSTRUCTIONS_INTERVAL_DIJKEMA_PRIMAL_WAVELET_H

#include <lawa/constructions/support.h>
#include <lawa/constructions/basis.h>
#include <lawa/constructions/basisfunction.h>
#include <lawa/constructions/wavelet.h>

namespace lawa {

template <typename _T, Construction _Cons>
struct Wavelet<_T,Primal,Interval,_Cons>
    : public BasisFunction<_T,Primal,Interval,_Cons>
{
    typedef _T T;
    static const FunctionSide   Side   = Primal;
    static const DomainType     Domain = Interval;
    static const Construction   Cons   = _Cons;

    Wavelet(const Basis<T,Primal,Interval,Cons> &_basis);

    const T
    operator()(T x, int j, Integer k, unsigned short deriv) const;

    const Support<T>
    support(int j, int k) const;

    const DenseVector<T>
    singularSupport(int j, long k) const;

    const int
    vanishingMoments(int j, long k) const;

    const T
    tic(int j) const;

    const Basis<T,Primal,Interval,Cons> &basis;
};

} // namespace lawa

#endif // LAWA_CONSTRUCTIONS_INTERVAL_DIJKEMA_PRIMAL_WAVELET_H

