/*
  This file is part of LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008-201  Schalk, Alexander Stippler.

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

#ifndef LAWA_CONSTRUCTIONS_REALLINE_DUAL_WAVELET_H
#define LAWA_CONSTRUCTIONS_REALLINE_DUAL_WAVELET_H 1

#include <lawa/enum.h>
#include <lawa/flensforlawa.h>
#include <lawa/constructions/basisfunction.h>
#include <lawa/constructions/realline/primal/bspline.h>
#include <lawa/constructions/realline/dual/bspline.h>
#include <lawa/constructions/wavelet.h>

namespace lawa {

using namespace flens;

template <typename _T>
struct Wavelet<_T,Dual,R,CDF>
    : public BasisFunction<_T,Dual,R,CDF>
{
    typedef _T T;
    static const FunctionSide Side = Dual;
    static const DomainType Domain = R;
    static const Construction Cons = CDF;

    Wavelet(int _d, int _d_);

    Wavelet(const BSpline<T,Primal,R,CDF> &_phi,
            const BSpline<T,Dual,R,CDF> &_phi_);

    T
    operator()(T x, int j, long k, unsigned short deriv=0) const;

    Support<T>
    support(int j, long k) const;

    const DenseVector<T> &
    mask() const;

    static DenseVector<T>
    mask(int d, int d_);

    const int d, d_, mu;
    const int l1_, l2_;
    const DenseVector<T> b_;
    const BSpline<T,Primal,R,CDF> phi;
    const BSpline<T,Dual,R,CDF> phi_;        
};

} // namespace lawa

#endif // LAWA_CONSTRUCTIONS_REALLINE_DUAL_WAVELET_H

