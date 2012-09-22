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

#ifndef LAWA_CONSTRUCTIONS_REALLINE_PRIMAL_BASIS_TCC
#define LAWA_CONSTRUCTIONS_REALLINE_PRIMAL_BASIS_TCC 1

#include <lawa/constructions/realline/primal/wavelet.h>
#include <lawa/constructions/realline/refinementmatrix.h>

namespace lawa {

template <typename T>
Basis<T,Primal,R,CDF>::Basis(int _d, int _d_, int j)
    : d(_d), d_(_d_), j0(j), mra(d), mra_(d,d_), psi(d,d_), M1(psi)
{
}

template <typename T>
const int
Basis<T,Primal,R,CDF>::level() const
{
    return _j;
}

template <typename T>
void
Basis<T,Primal,R,CDF>::setLevel(int j) const
{
    _j = j;
}

template <typename T>
const BasisFunction<T,Primal,R,CDF> &
Basis<T,Primal,R,CDF>::generator(XType xtype) const
{
    if (xtype==XBSpline) {
        return mra.phi;
    } else {
        return psi;
    }
}

} // namespace lawa

#endif // LAWA_CONSTRUCTIONS_REALLINE_PRIMAL_BASIS_TCC
