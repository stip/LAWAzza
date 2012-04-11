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

#ifndef LAWA_CONSTRUCTIONS_REALLINE_REFINEMENTMATRIX_H
#define LAWA_CONSTRUCTIONS_REALLINE_REFINEMENTMATRIX_H 1

#include <lawa/flensforlawa.h>
#include <lawa/constructions/refinementmatrix.h>
#include <lawa/constructions/bspline.h>
#include <lawa/constructions/wavelet.h>
#include <lawa/enum.h>

namespace flens {

template <typename T>
class RefinementMatrix<T,R,CDF>
    : public Matrix<RefinementMatrix<T,R,CDF> >
{
    public:
        typedef T ElementType;

        template <FunctionSide Side>
            RefinementMatrix(const lawa::BSpline<T,Side,R,CDF> &bw);

        template <FunctionSide Side>
            RefinementMatrix(const lawa::Wavelet<T,Side,R,CDF> &bw);

        DenseVector<Array<T> > band;
};

template <typename T>
struct TypeInfo<RefinementMatrix<T,R,CDF> >
{
    typedef RefinementMatrix<T,R,CDF> Impl;
};

//------------------------------------------------------------------------------

template <typename X, typename Y>
void
mv(cxxblas::Transpose transA, typename X::ElementType alpha,
   const RefinementMatrix<typename X::ElementType,R,CDF> &A,
   const DenseVector<X> &x, typename X::ElementType beta, DenseVector<Y> &y);

} // namespace flens

#endif // LAWA_CONSTRUCTIONS_REALLINE_REFINEMENTMATRIX_H

