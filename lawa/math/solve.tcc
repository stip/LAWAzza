/*
 LAWA - Library for Adaptive Wavelet Applications.
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

#ifndef LAWA_MATH_SOLVE_TCC
#define LAWA_MATH_SOLVE_TCC 1

#include <cmath>
#include <lawa/math/solve.h>

namespace lawa {

template <typename T>
int
solve(GeMatrix<T> &A, DenseVector<T> &b)
{
    const int i0 = A.firstRow();
    const int j0 = A.firstCol();
    const int I0 = b.firstIndex();

    A.changeIndexBase(1, 1);
    b.changeIndexBase(1);

    DenseVector<int> piv;
    int info = flens::lapack::sv(A, piv, b);

    A.changeIndexBase(i0, j0);
    b.changeIndexBase(I0);
    return info;
}

} // namespace lawa

#endif // LAWA_MATH_SOLVE_TCC
