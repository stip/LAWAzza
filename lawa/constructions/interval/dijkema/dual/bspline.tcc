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

#ifndef LAWA_CONSTRUCTIONS_INTERVAL_DIJKEMA_DUAL_BSPLINE_TCC
#define LAWA_CONSTRUCTIONS_INTERVAL_DIJKEMA_DUAL_BSPLINE_TCC 1

#include <cassert>

namespace lawa {

template <typename T>
BSpline<T,Dual,Interval,Dijkema>::BSpline(const MRA<T,Dual,Interval,Dijkema>
                                                                         &_mra_)
    : mra_(_mra_), _phi_R(mra_.d, mra_.d_)
{
}

template <typename T>
const T
BSpline<T,Dual,Interval,Dijkema>::operator()(T x, int j, int k) const
{
    assert(j>=mra_.j0);
    assert(k>=mra_.rangeI_(j).firstIndex());
    assert(k<=mra_.rangeI_(j).lastIndex());

    if (k<=mra_.rangeI_L().lastIndex()) {
        T value = 0.0;

        int l = -_phi_R.a_.lastIndex()+1;

        const int r0 = mra_.R_Right.firstRow();
        const int r1 = mra_.R_Right.lastRow();

        for (int r=r0; r<=r1; ++r, ++l) {
            value += mra_.R_Left(r,k) * _phi_R(x,j,l);
        }
//        assert(l==-_phi_R.a_.firstIndex());
        return value;
    }

    if (k>=mra_.rangeI_R(j).firstIndex()) {
        k -= mra_.rangeI_R(j).firstIndex()-1;
        T value = 0.0;

        int l = pow2i<int>(j)-_phi_R.a_.lastIndex()+1;

        const int r0 = mra_.R_Right.firstRow();
        const int r1 = mra_.R_Right.lastRow();

        for (int r=r0; r<=r1; ++r, ++l) {
            value += mra_.R_Right(r,k) * _phi_R(x,j,l);
        }
        return value;
    }

    return _phi_R(x,j,k-(mra_.d + mra_.d_ - 1)-_phi_R.l1_);
}

template <typename T>
const Support<T>
BSpline<T,Dual,Interval,Dijkema>::support(int j, int k) const
{
    if (k<=mra_.rangeI_L().lastIndex()) {
        return Support<T>(0.,pow2i<T>(-j)*(_phi_R.a_.length()-2));
    }
    if (k>=mra_.rangeI_R(j).firstIndex()) {
        return Support<T>(1.-pow2i<T>(-j)*(_phi_R.a_.length()-2), 1.);
    }
    return pow2i<T>(-j)*Support<T>(k-(mra_.d+mra_.d_-1),
                                   k-(mra_.d+mra_.d_-1)+_phi_R.l2_-_phi_R.l1_);
}

} // namespace lawa

#endif // LAWA_CONSTRUCTIONS_INTERVAL_DIJKEMA_DUAL_BSPLINE_TCC
