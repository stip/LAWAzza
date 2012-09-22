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

#ifndef LAWA_CONSTRUCTIONS_INTERVAL_DIJKEMA_PRIMAL_MRA_H
#define LAWA_CONSTRUCTIONS_INTERVAL_DIJKEMA_PRIMAL_MRA_H 1

#include <lawa/flensforlawa.h>
#include <lawa/constructions/basisfunction.h>
#include <lawa/constructions/mra.h>
#include <lawa/constructions/support.h>
#include <lawa/constructions/interval/refinementmatrix.h>


/*
 *  Disclaimer:
 *
 *  Historically this MRA was developed by Miriam Primbs first.  Note that
 *  the Primbs and Dijkema bases have the same primal MRA.
 *
 */

namespace lawa {

using namespace flens;

template <typename _T>
class MRA<_T,Primal,Interval,Dijkema>
{
    public:
        typedef _T T;
        static const FunctionSide Side = Primal;
        static const DomainType Domain = Interval;
        static const Construction Cons = Dijkema;

        typedef BasisFunction<T,Primal,Interval,Dijkema>  BasisFunctionType;
        typedef BSpline<T,Primal,Interval,Dijkema>        BSplineType;

        MRA(int d, int j=-1);

        ~MRA();

        // cardinalities of whole, left, inner, right index sets.
        const int
        cardI(int j) const;

        const int
        cardIL() const;

        const int
        cardII(int j) const;

        const int
        cardIR() const;

        // ranges of whole left, inner, right index sets.
        const Range<int>
        rangeI(int j) const;

        const Range<int>
        rangeIL() const;

        const Range<int>
        rangeII(int j) const;

        const Range<int>
        rangeIR(int j) const;

        const int
        level() const;

        void
        setLevel(int j) const;

        template <BoundaryCondition BCLeft, BoundaryCondition BCRight>
            void
            enforceBoundaryConditions();

        BoundaryCondition
        boundaryCondition(BoundarySide side) const;

        const int d, mu;       // mu = mu(d) = d&1.
        const int min_j0;      // minimal allowed(!) level;
        const int j0;          // minimal used(!) level.


        RefinementMatrix<T,Interval,Dijkema>    M0;

    private:

        void
        _calcM0();

        DenseVector<BoundaryCondition>  _bc;
        mutable int                     _j;     // the current level.

    public:
        BSpline<T,Primal,Interval,Dijkema> phi;

};

} // namespace lawa

#endif // LAWA_CONSTRUCTIONS_INTERVAL_DIJKEMA_PRIMAL_MRA_H

