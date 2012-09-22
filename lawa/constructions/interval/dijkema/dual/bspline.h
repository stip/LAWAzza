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

#ifndef LAWA_CONSTRUCTIONS_INTERVAL_DIJKEMA_DUAL_BSPLINE_H
#define LAWA_CONSTRUCTIONS_INTERVAL_DIJKEMA_DUAL_BSPLINE_H 1

#include <lawa/flensforlawa.h>

#include <lawa/constructions/basisfunction.h>
#include <lawa/constructions/bspline.h>
#include <lawa/constructions/support.h>
#include <lawa/constructions/mra.h>
#include <lawa/enum.h>

namespace lawa {

template <typename _T>
struct BSpline<_T,Dual,Interval,Dijkema>
    : public BasisFunction<_T,Dual,Interval,Dijkema>
{
    typedef _T T;
    static const FunctionSide Side   = Dual;
    static const DomainType   Domain = Interval;
    static const Construction Cons   = Dijkema;

    BSpline(const MRA<T,Dual,Interval,Dijkema> &mra_);

    const T
    operator()(T x, int j, int k) const;

    const Support<T>
    support(int j, int k) const;

    const MRA<T,Dual,Interval,Dijkema> &mra_;

    private:
        const BSpline<T,Dual,R,CDF>  _phi_R;
};

} // namespace lawa

#endif // LAWA_CONSTRUCTIONS_INTERVAL_DIJKEMA_DUAL_BSPLINE_H

