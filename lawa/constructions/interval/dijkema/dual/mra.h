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

#ifndef LAWA_CONSTRUCTIONS_INTERVAL_DIJKEMA_DUAL_MRA_H
#define LAWA_CONSTRUCTIONS_INTERVAL_DIJKEMA_DUAL_MRA_H 1

#include <lawa/constructions/basisfunction.h>
#include <lawa/constructions/mra.h>
#include <lawa/constructions/interval/refinementmatrix.h>
#include <lawa/enum.h>

namespace lawa {

using namespace flens;

template <typename _T>
class MRA<_T,Dual,Interval,Dijkema>
{
    public:
        typedef _T T;
        static const FunctionSide  Side   = Dual;
        static const DomainType    Domain = Interval;
        static const Construction  Cons   = Dijkema;

        typedef BasisFunction<T,Dual,Interval,Dijkema>  BasisFunctionType;
        typedef BSpline<T,Dual,Interval,Dijkema>        BSplineType;

        MRA(int d, int d_, int j=-1);

        ~MRA();

        // cardinalities of whole, left, inner, right index sets.
        const int
        cardI_(int j) const;

        const int
        cardI_L() const;

        const int
        cardI_I(int j) const;

        const int
        cardI_R() const;

        // ranges of whole left, inner, right index sets.
        const Range<int>
        rangeI_(int j) const;

        const Range<int>
        rangeI_L() const;

        const Range<int>
        rangeI_I(int j) const;

        const Range<int>
        rangeI_R(int j) const;

        const int
        level() const;

        void
        setLevel(int j) const;

        template <BoundaryCondition LeftBC, BoundaryCondition RightBC>
            void
            enforceBoundaryConditions();

        BoundaryCondition
        boundaryCondition(BoundarySide side) const;

        const int d, d_, mu;   // mu = mu(d) = d&1.

    public:
        const int min_j0;      // minimal allowed(!) level;
        const int j0;          // minimal used(!) level.

        GeMatrix<T>                                  R_Left, R_Right;
        RefinementMatrix<T,Interval,Dijkema>         M0_;

    private:

        void
        _calcC_LR(BoundarySide side, GeMatrix<T> &C_LR);

        void
        _calcC_L(GeMatrix<T> &C_L);

        void
        _calcC_R(GeMatrix<T> &C_R);

        void
        _calcM0_();

        DenseVector<BoundaryCondition>              _bc;
        mutable int                                 _j;
    
    public:
        BSpline<T,Dual,Interval,Dijkema>            phi_;
};

} // namespace lawa

#endif // LAWA_CONSTRUCTIONS_INTERVAL_DIJKEMA_DUAL_MRA_H

