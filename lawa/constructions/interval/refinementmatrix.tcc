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

#ifndef LAWA_CONSTRUCTIONS_INTERVAL_REFINEMENTMATRIX_TCC
#define LAWA_CONSTRUCTIONS_INTERVAL_REFINEMENTMATRIX_TCC 1

#include <cassert>
#include <cmath>
#include <lawa/auxiliary/auxiliary.h>
#include <lawa/flensforlawa.h>

namespace flens {

using namespace lawa;

template <typename T, Construction Cons>
RefinementMatrix<T,Interval,Cons>::RefinementMatrix()
{
}

template <typename T, Construction Cons>
RefinementMatrix<T,Interval,Cons>::RefinementMatrix(int nLeft, int nRight,
                                                    int leftBandOffset,
                                                    int rightBandOffset,
                                                    const lawa::GeMatrix<T> &A,
                                                    int _min_j0, int cons_j)
    : left(nLeft, nLeft>0 ? A.firstCol() : 1),
      right(nRight, nRight>0 ? A.lastCol()-nRight+1 : A.lastCol()+2),
      min_j0(_min_j0), _cons_j(cons_j), _j(cons_j),
      _firstRow(A.firstRow()), _firstCol(A.firstCol()),
      _lastRow(A.lastRow()), _lastCol(A.lastCol()),
      _additionalRows(0), _additionalCols(0)
{
    assert(nLeft>=0);
    assert(nRight>=0);
    assert(_cons_j>=min_j0);
    assert(_firstCol>=1);
    _extractMasks(A, leftBandOffset, rightBandOffset);
}

template <typename T, Construction Cons>
void
RefinementMatrix<T,Interval,Cons>::operator=(const RefinementMatrix &rhs)
{
    if (this!=&rhs) {
        leftBand.resize(rhs.leftBand);
        rightBand.resize(rhs.rightBand);
        left.resize(0);
        right.resize(0);
    
        leftBand  = rhs.leftBand;
        rightBand = rhs.rightBand;
        left      = rhs.left;
        right     = rhs.right;
    
        min_j0          = rhs.min_j0;
        _cons_j         = rhs._cons_j;
        _j              = rhs._j;
        _firstRow       = rhs._firstRow;
        _firstCol       = rhs._firstCol;
        _lastRow        = rhs._lastRow;
        _lastCol        = rhs._lastCol;
        _additionalRows = rhs._additionalRows;
        _additionalCols = rhs._additionalCols;
    }
}

template <typename T, Construction Cons>
const typename lawa::DenseVector<T>::ConstView
RefinementMatrix<T,Interval,Cons>::operator()(int j,
                                              const Underscore<int> &/*_*/,
                                              int col) const
{
    assert(j>=min_j0);

    int additionalCols = 0, additionalRows = 0;
    if (j>_cons_j) {
        for (int l=_cons_j; l<j; ++l) {
            additionalCols += pow2i<int>(l);
        }
        additionalRows = 2*additionalCols;
    } else if (j<_cons_j) {
        for (int l=_cons_j-1; l>=j; --l) {
            additionalCols -= pow2i<int>(l);
        }
        additionalRows = 2*additionalCols;
    }

    assert(col>=_firstCol);
    assert(col<=_lastCol+additionalCols);

    if (col<=left.lastIndex()) {
        return left(col);
    }

    if (col>=right.firstIndex()+additionalCols) {
        const auto &rightCol = right(col-additionalCols);
        const int firstIndex = rightCol.firstIndex()+additionalRows;
        return rightCol( _ , firstIndex);
    }

    const int middle     = (_firstCol+_lastCol+additionalCols)/2;
    const int firstIndex = leftBand.firstIndex() + 2*(col-left.lastIndex()-1);

    return (col>middle) ? rightBand( _ , firstIndex)
                        : leftBand( _ , firstIndex);
}

template <typename T, Construction Cons>
const Range<int>
RefinementMatrix<T,Interval,Cons>::rows() const
{
    return _(firstRow(), lastRow());
}

template <typename T, Construction Cons>
const Range<int>
RefinementMatrix<T,Interval,Cons>::cols() const
{
    return _(firstCol(), lastCol());
}

template <typename T, Construction Cons>
const int
RefinementMatrix<T,Interval,Cons>::numRows() const
{
    return lastRow()-firstRow()+1;
}

template <typename T, Construction Cons>
const int
RefinementMatrix<T,Interval,Cons>::numCols() const
{
    return lastCol()-firstCol()+1;
}

template <typename T, Construction Cons>
const int
RefinementMatrix<T,Interval,Cons>::firstRow() const
{
    return _firstRow;
}

template <typename T, Construction Cons>
const int
RefinementMatrix<T,Interval,Cons>::lastRow() const
{
    return _lastRow + _additionalRows;
}

template <typename T, Construction Cons>
const int
RefinementMatrix<T,Interval,Cons>::firstCol() const
{
    return _firstCol;
}

template <typename T, Construction Cons>
const int
RefinementMatrix<T,Interval,Cons>::lastCol() const
{
    return _lastCol + _additionalCols;
}

template <typename T, Construction Cons>
const int
RefinementMatrix<T,Interval,Cons>::level() const
{
    return _j;
}

template <typename T, Construction Cons>
void
RefinementMatrix<T,Interval,Cons>::setLevel(int j) const
{
    if (j<_j) {
        assert(j>=min_j0);
        for (int l=_j-1; l>=j; --l) {
            _additionalCols -= pow2i<int>(l);
        }
        _additionalRows = 2*_additionalCols;
        _j = j;
        return;
    }
    if (j>_j) {
        for (int l=_j; l<j; ++l) {
            _additionalCols += pow2i<int>(l);
        }
        _additionalRows = 2*_additionalCols;
        _j = j;
        return;
    }
}

template <typename T, Construction Cons>
void
RefinementMatrix<T,Interval,Cons>::_extractMasks(const lawa::GeMatrix<T> &A,
                                                 int leftBandOffset,
                                                 int rightBandOffset)
{
    using std::abs;

    assert(leftBandOffset>=0);
    assert(rightBandOffset>=0);

/*
    std::cerr << "_extractMasks:" << std::endl;
    std::cerr << "A = " << A << std::endl;

    std::cerr << "left =  " << left << std::endl;
    std::cerr << "right = " << right << std::endl;
*/

//
//  extract left block
//
    for (int c=A.firstCol(); c<A.firstCol()+left.length(); ++c) {
        int r = A.lastRow();
        while (abs(A(r,c))<=1e-12) {
            --r;
            assert(r>=A.firstRow());
        }
        left(c) = A(_(A.firstRow(),r),c);
        left(c).changeIndexBase(A.firstRow());
    }

//
//  extract right block
//
    for (int c=A.lastCol()-right.length()+1; c<=A.lastCol(); ++c) {
        int r = A.firstRow();
        while (abs(A(r,c))<=1e-12) {
            ++r;
            assert(r<=A.lastRow());
        }
        right(c) = A(_(r,A.lastRow()), c);
        right(c).changeIndexBase(r);
    }

//
//  extract band (left, middle]
//
    int c = A.firstCol()+left.length()+leftBandOffset;
    
    int first = A.firstRow();
    while (abs(A(first,c))<=1e-12) {
        ++first;
        assert(first<=A.lastRow());
    }    
    int last = A.lastRow();
    while (abs(A(last,c))<=1e-12) {
        --last;
        assert(last>=A.firstRow());
    }
    leftBand = A(_(first,last), c);
    leftBand.changeIndexBase(first);
#   ifdef CHECK_INTERVAL_CONSTRUCTION
    for (++c; c<=(A.firstCol()+A.lastCol())/2; ++c) {
        int i=leftBand.firstIndex();
        first += 2; last += 2;
        for (int r=first; r<=last; ++r, ++i) {
            assert(abs(leftBand(i)-A(r,c))<=1e-12);
        }
    }
#   endif

//
//  extract band (middle, right)
//
    c = A.lastCol()-right.length()-rightBandOffset;
    first = A.firstRow();
    while (abs(A(first,c))<=1e-12) {
        ++first;
        assert(first<=A.lastRow());
    }    
    last = A.lastRow();
    while (abs(A(last,c))<=1e-12) {
        --last;
        assert(last>=A.firstRow());
    }
    rightBand = A(_(first,last), c);
    rightBand.changeIndexBase(first);
#   ifdef CHECK_INTERVAL_CONSTRUCTION
    for (--c; c>(A.firstCol()+A.lastCol())/2; --c) {
        int i=rightBand.firstIndex();
        first -= 2; last -= 2;
        for (int r=first; r<=last; ++r, ++i) {
            assert(abs(rightBand(i)-A(r,c))<=1e-12);
        }
    }

//
//  Consistency check for inner bands
//
    assert(leftBand.length()==rightBand.length());
#   endif

/*
    std::cerr << "leftBand = " << leftBand << std::endl;
    std::cerr << "rightBand = " << rightBand << std::endl;
*/
}

//------------------------------------------------------------------------------

template <typename X, Construction Cons, typename Y>
void
mv(Transpose                                                       transA,
   typename X::ElementType                                         alpha,
   const RefinementMatrix<typename X::ElementType,Interval,Cons>   &A,
   const DenseVector<X>                                            &x,
   typename X::ElementType beta, DenseVector<Y>                    &y)
{
    typedef typename X::ElementType T;
    assert(alpha==T(1));
    assert(x.stride()==1);
    assert(y.stride()==1);

/*
    cerr << "x = " << x << endl;
    
    cerr << "A.left = " << A.left << endl;
    cerr << "A.right = " << A.right << endl;

    cerr << "A.leftBand = " << A.leftBand << endl;
    cerr << "A.rightBand = " << A.rightBand << endl;
*/
    if (transA==NoTrans) {
        assert(A.numCols()==x.length());

        if (beta==T(0)) {
            y.engine().resize(A.rows()) || y.fill(T(0));
        } else {
            assert(y.length()==A.numRows());
            y.changeIndexBase(A.firstRow());
        }

        // left upper block
        int ix = x.firstIndex();

/*
        cerr << "left upper block:" << endl;
        cerr << "   from: " << A.left.firstIndex() << endl;
        cerr << "   to:   " << A.left.lastIndex() << endl;
*/
        for (int c=A.left.firstIndex(); c<=A.left.lastIndex(); ++c, ++ix) {
            int n = A.left(c).length();
            cxxblas::axpy(n, 
                          x(ix),
                          A.left(c).data(), 1,
                          y.data(), 1);
        }

        // central band (up to middle)
        int iy = A.leftBand.firstIndex()-A.firstRow();
        int n = A.leftBand.length();

        //int middle = iceil<int>(x.length()/2.)+A.firstCol();
        //int middle = iceil<int>((A.firstCol()+A.lastCol())/2.0);
        //const int middle = iceil<int>((A.firstCol()+A.lastCol()+A.left.length()-A.right.length())/2.0);
        const int middle = A.firstCol()+A.left.length()-1
                         + iceil<int>((x.length()-A.left.length()-A.right.length())/2.0);

/*
        cerr << "central band (up to middle):" << endl;
        cerr << "   from: " << A.left.lastIndex()+1 << endl;
        cerr << "   to:   " << middle << endl;
*/
        for (int c=A.left.lastIndex()+1; c<=middle; ++c, iy+=2, ++ix) {
            cxxblas::axpy(n,
                          x(ix),
                          A.leftBand.data(), 1,
                          y.data()+iy, 1);
        }
        
        // central band (right of middle)
        int end = A.left.firstIndex() + x.length() - A.right.length();
/*
        cerr << "central band (right of middle):" << endl;
        cerr << "   from: " << middle+1 << endl;
        cerr << "   to:   " << end-1 << endl;
*/
        for (int c=middle+1; c<end; ++c, iy+=2, ++ix) {
            cxxblas::axpy(n,
                          x(ix),
                          A.rightBand.data(), 1,
                          y.data()+iy, 1);
        }

        // right lower block
/*
        cerr << "right lower block:" << endl;
        cerr << "   from: " << A.right.firstIndex() << endl;
        cerr << "   to:   " << A.right.lastIndex() << endl;
*/
        for (int c=A.right.firstIndex(); c<=A.right.lastIndex(); ++c, ++ix) {
            int n = A.right(c).length();
            cxxblas::axpy(n, 
                          x(ix),
                          A.right(c).data(), 1,
                          y.data()+y.length()-1-n+1, 1);
        }
    } else { // transA==Trans
        assert(A.numRows()==x.length());

        if (beta==T(0)) {
            y.engine().resize(A.cols());
        } else {
            assert(y.length()==A.numCols());
            y.changeIndexBase(A.firstCol());
        }
        int iy = y.firstIndex();
        // left upper block
        for (int c=A.left.firstIndex(); c<=A.left.lastIndex(); ++c, ++iy) {
            int n = A.left(c).length();
            cxxblas::dot(n,
                         A.left(c).data(), 1,
                         x.data(), 1, 
                         y(iy));
        }

        // central band (up to middle)
        int middle = y.length()/2;
        int ix = A.leftBand.firstIndex() - A.firstRow();
        for (int i=A.left.length()+1; i<=middle; ++i, ix+=2, ++iy) {
            cxxblas::dot(A.leftBand.length(),
                         A.leftBand.data(), 1,
                         x.data()+ix, 1,
                         y(iy));
        }
        // central band (right of middle)
        int end = y.length() - A.right.length();
        for (int i=middle+1; i<=end; ++i, ix+=2, ++iy) {
            cxxblas::dot(A.rightBand.length(),
                         A.rightBand.data(), 1,
                         x.data()+ix, 1,
                         y(iy));
        }
        // right lower block
        for (int c=A.right.firstIndex(); c<=A.right.lastIndex(); ++c, ++iy) {
            int n = A.right(c).length();
            cxxblas::dot(n, 
                         A.right(c).data(), 1,
                         x.data() + A.numRows() - n, 1, 
                         y(iy));
        }
    }
}

} // namespace flens

#endif // LAWA_CONSTRUCTIONS_REFINEMENTMATRIX_TCC
