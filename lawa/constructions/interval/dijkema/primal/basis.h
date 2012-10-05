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

#ifndef LAWA_CONSTRUCTIONS_INTERVAL_DIJKEMA_PRIMAL_BASIS_H
#define LAWA_CONSTRUCTIONS_INTERVAL_DIJKEMA_PRIMAL_BASIS_H 1

#include <lawa/constructions/basis.h>
#include <lawa/constructions/basisfunction.h>
#include <lawa/constructions/interval/dijkema/dual/mra.h>
#include <lawa/constructions/interval/dijkema/primal/mra.h>
#include <lawa/enum.h>

namespace lawa {
    
template <typename _T>
class Basis<_T,Primal,Interval,Dijkema>
{
    public:
        typedef _T T;
        static const FunctionSide   Side   = Primal;
        static const DomainType     Domain = Interval;
        static const Construction   Cons   = Dijkema;

        typedef BasisFunction<T,Primal,Interval,Dijkema> BasisFunctionType;
        typedef BSpline<T,Primal,Interval,Dijkema>       BSplineType;
        typedef Wavelet<T,Primal,Interval,Dijkema>       WaveletType;

        Basis(int _d, int _d_, int j=-1);
        
        const int
        level() const;

        void
        setLevel(int j) const;

        template <BoundaryCondition BCLeft, BoundaryCondition BCRight>
            void
            enforceBoundaryConditions();

        const BasisFunctionType &
        generator(XType xtype) const;

        // cardinalities of whole, left, inner, right index sets (primal).
        const int
        cardJ(int j) const;

        const int
        cardJL() const;

        const int
        cardJI(int j) const;

        const int
        cardJR() const;

        // ranges of whole, left, inner, right index sets (primal).
        const Range<int>
        rangeJ(int j) const;

        const Range<int>
        rangeJL() const;

        const Range<int>
        rangeJI(int j) const;

        const Range<int>
        rangeJR(int j) const;

        MRA<T,Primal,Interval,Dijkema>          mra;
        MRA<T,Dual,Interval,Dijkema>            mra_;

        RefinementMatrix<T,Interval,Dijkema>    M1;

        const int d, d_, mu;   // mu = mu(d) = d&1.
        const int min_j0;      // minimal allowed(!) level;
        const int j0;          // minimal used(!) level.

    private:
        DenseVector<BoundaryCondition>      _bc;
        mutable int                         _j;    // the current level.

    public:
        Wavelet<T,Primal,Interval,Dijkema>  psi;
};

} // namespace lawa

#endif // LAWA_CONSTRUCTIONS_INTERVAL_DIJKEMA_PRIMAL_BASIS_H

