/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008,2009  Mario Rometsch, Alexander Stippler.

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

#ifndef LAWA_AUXILIARY_DENSIFY_TCC
#define LAWA_AUXILIARY_DENSIFY_TCC 1

#include <lawa/auxiliary/densify.h>

namespace lawa {

template <typename I>
void
densify(cxxblas::Transpose trans, const Matrix<I> &A,
        GeMatrix<typename I::ElementType> &D,
        int firstRow, int firstCol)
{
    typedef typename I::ElementType  ElementType;

    if (trans==cxxblas::NoTrans) {
        int m = A.impl().numRows();
        int n = A.impl().numCols();
        D.resize(m,n,firstRow,firstCol);
        DenseVector<ElementType> e(n);
        for (int i=1; i<=n; ++i) {
            DenseVector<ElementType> y(m);
            e(i) = 1.;
            D(_,i+firstCol-1) = A.impl()*e;
            e(i) = 0.;
        }
    } else {
        assert(trans==cxxblas::Trans);

        int m = A.impl().numRows();
        int n = A.impl().numCols();
        D.resize(n,m,firstCol,firstRow);
        DenseVector<ElementType> e(m);
        for (int i=1; i<=m; ++i) {
            DenseVector<ElementType> y(n);
            e(i) = 1.;
            D(_,i+firstCol-1) = transpose(A.impl())*e;
            e(i) = 0.;
        }
    }
}

} // namespace lawa

#endif // LAWA_AUXILIARY_DENSIFY_TCC
