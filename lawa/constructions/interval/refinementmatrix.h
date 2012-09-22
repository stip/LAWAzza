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

#ifndef LAWA_CONSTRUCTIONS_INTERVAL_REFINEMENTMATRIX_H
#define LAWA_CONSTRUCTIONS_INTERVAL_REFINEMENTMATRIX_H 1

#include <lawa/constructions/refinementmatrix.h>
#include <lawa/flensforlawa.h>
#include <lawa/typedefs.h>

namespace flens {

using namespace lawa;
using namespace cxxblas;

template <typename T, Construction Cons>
class RefinementMatrix<T, Interval, Cons>
    : public Matrix<RefinementMatrix<T, Interval, Cons> >
{
    public:
        typedef T ElementType;

        RefinementMatrix();

        RefinementMatrix(int nLeft, int nRight,
                         const lawa::GeMatrix<T> &A,
                         int _min_j0, int cons_j);

        void
        operator=(const RefinementMatrix &rhs);

        const typename lawa::DenseVector<T>::ConstView
        operator()(int j, const Underscore<int> &u, int col) const;
        
        const Range<int>
        rows() const;

        const Range<int>
        cols() const;

        const int
        numRows() const;

        const int
        numCols() const;

        const int
        firstRow() const;

        const int
        lastRow() const;

        const int
        firstCol() const;

        const int
        lastCol() const;

        const int
        level() const;

        void
        setLevel(int j) const;

        lawa::DenseVector<lawa::DenseVector<T> >  left, right;
        lawa::DenseVector<T>                      leftBand, rightBand;
        int                                       min_j0;

    private:

        void
        _extractMasks(const lawa::GeMatrix<T> &A);

        int          _cons_j;
        mutable int  _j;
        int          _firstRow, _firstCol, _lastRow, _lastCol;
        mutable int  _additionalRows, _additionalCols;
};


template <typename X, Construction Cons, typename Y>
void
mv(Transpose transA, typename X::ElementType alpha,
   const RefinementMatrix<typename X::ElementType,Interval,Cons> &A,
   const DenseVector<X> &x, typename X::ElementType beta, DenseVector<Y> &y);

} // namespace flens

#endif // LAWA_CONSTRUCTIONS_REFINEMENTMATRIX_H
