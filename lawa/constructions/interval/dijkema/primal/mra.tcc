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

#ifndef LAWA_CONSTRUCTIONS_INTERVAL_DIJKEMA_PRIMAL_MRA_TCC
#define LAWA_CONSTRUCTIONS_INTERVAL_DIJKEMA_PRIMAL_MRA_TCC 1

#include <cassert>
#include <cmath>
#include <list>

#include <lawa/auxiliary/arrow.h>
#include <lawa/constructions/interval/dijkema/primal/splinehelper.h>
#include <lawa/constructions/realline/primal/bspline.h>
#include <lawa/enum.h>
#include <lawa/math/math.h>

namespace lawa {

template <typename T>
MRA<T,Primal,Interval,Dijkema>::MRA(int _d, int j)
    : d(_d),
      mu(d&1),
      min_j0(iceil<int>(log2(2*d))),     // TODO: Check with Primbs thesis
      j0(j==-1 ? min_j0 : j),
      _bc(2,0),
      _j(j0),
      phi(*this)
{
    assert(d>1);
    assert(_j>=min_j0);

    _calcM0();
    setLevel(_j);
}

template <typename T>
MRA<T,Primal,Interval,Dijkema>::~MRA()
{
}

//--- cardinalities of whole, left, inner, right index sets. -------------------

template <typename T>
const int
MRA<T,Primal,Interval,Dijkema>::cardI(int j) const
{
    return cardIL() + cardII(j) + cardIR();


    assert(j>=min_j0);
    return pow2i<int>(j) + d - 1 - (_bc(0)+_bc(1));
}

template <typename T>
const int
MRA<T,Primal,Interval,Dijkema>::cardIL() const
{
    return d - 1 -_bc(0);
}

template <typename T>
const int
MRA<T,Primal,Interval,Dijkema>::cardII(int j) const
{
    assert(j>=min_j0);
    return pow2i<int>(j) - d + 1;
}

template <typename T>
const int
MRA<T,Primal,Interval,Dijkema>::cardIR() const
{
    return d - 1 - _bc(1);
}

//--- ranges of whole, left, inner, right index sets. --------------------------

template <typename T>
const Range<int>
MRA<T,Primal,Interval,Dijkema>::rangeI(int j) const
{
    assert(j>=min_j0);
    Range<int>  range(1 + _bc(0), pow2i<int>(j) + d - 1 - _bc(1));
    assert(range.length()==cardI(j));
    return range;
}

template <typename T>
const Range<int>
MRA<T,Primal,Interval,Dijkema>::rangeIL() const
{
    Range<int>  range(1 + _bc(0), d - 1);
    assert(range.length()==cardIL());
    return range;
}

template <typename T>
const Range<int>
MRA<T,Primal,Interval,Dijkema>::rangeII(int j) const
{
    assert(j>=min_j0);
    Range<int>  range(d, pow2i<int>(j));
    assert(range.length()==cardII(j));
    return range;
}

template <typename T>
const Range<int>
MRA<T,Primal,Interval,Dijkema>::rangeIR(int j) const
{
    assert(j>=min_j0);
    Range<int>  range(pow2i<int>(j) + 1, pow2i<int>(j) + d - 1 - _bc(1));
    assert(range.length()==cardIR());
    return range;
}

template <typename T>
const int
MRA<T,Primal,Interval,Dijkema>::level() const
{
    return _j;
}

template <typename T>
void
MRA<T,Primal,Interval,Dijkema>::setLevel(int j) const
{
    assert(j>=min_j0);
    _j = j;
    M0.setLevel(_j);
}

template <typename T>
template <BoundaryCondition BCLeft, BoundaryCondition BCRight>
void
MRA<T,Primal,Interval,Dijkema>::enforceBoundaryConditions()
{
    M0.left.resize(0);
    M0.right.resize(0);
    M0.leftBand.resize(0);
    M0.rightBand.resize(0);

    _bc(0) = BCLeft;
    _bc(1) = BCRight;
    _calcM0();
}

template <typename T>
BoundaryCondition
MRA<T,Primal,Interval,Dijkema>::boundaryCondition(BoundarySide side) const
{
    assert(side==Left || side==Right);

    if (side==Left) {
        return _bc(0);
    }
    return _bc(1);
}

//
// See Dijkema page 19ff
//
template <typename T>
void
MRA<T,Primal,Interval,Dijkema>::_calcM0()
{
    using std::max;
    using std::min;

    const T Zero(0), One(1);

    const BSpline<T,Primal,R,CDF>   phiR(d);
    const int                       l1 = phiR.a.firstIndex();
    const int                       l2 = phiR.a.lastIndex();
    const int                       j  = min_j0;

//
//  Apply the Oslo algorithm
//
    DenseVector<T> knots = linspace(-d+1., d-1., 2*d-1);
    knots.engine().changeIndexBase(1);
    GeMatrix<T>  Transformation(knots.length()-d, knots.length()-d);
    Transformation.diag(0) = 1.;
    for (int i=1; i<d; ++i) {
        GeMatrix<T> Tmp = insertKnot(d-1,knots,0.), Tmp2;
        blas::mm(cxxblas::NoTrans,cxxblas::NoTrans,
                 1.,Tmp,Transformation,0.,Tmp2);
        Transformation.resize(Tmp2);
        Transformation = Tmp2;
    }

//
//  R is square with dimensions [1-l2,...,2^j-l1-1] x [1-l2,...,2^j-l1-1]
//
    GeMatrix<T>  R(_(1-l2, pow2i<int>(j)-l1-1),
                   _(1-l2, pow2i<int>(j)-l1-1));

    R.diag(0) = 1.;

    const int         m            = Transformation.lastRow();
    GeMatrix<T>       TransInvLeft = Transformation(_(m-(d-1)+1,m), _);

    inv(TransInvLeft);

    TransInvLeft.changeIndexBase(R.firstRow(),R.firstCol());
    R(TransInvLeft) = TransInvLeft;

    GeMatrix<T>  TransInvRight;
    arrow(TransInvLeft, TransInvRight);
    TransInvRight.changeIndexBase(R.lastRow()-TransInvRight.numRows()+1,
                                  R.lastCol()-TransInvRight.numCols()+1);
    R(TransInvRight) = TransInvRight;

//
//  R on level j+1
//  (square with dimensions [1-l2,...,2^(j+1)-l1-1] x [1-l2,...,2^(j+1)-l1-1])
//
    GeMatrix<T>  RjPlus1(_(-l2+1,pow2i<int>(j+1)-l1-1),
                         _(-l2+1,pow2i<int>(j+1)-l1-1));
    RjPlus1.diag(0) = 1.;
    TransInvLeft.changeIndexBase(RjPlus1.firstRow(),RjPlus1.firstCol());
    RjPlus1(TransInvLeft) = TransInvLeft;
    TransInvRight.changeIndexBase(RjPlus1.lastRow()-TransInvRight.numRows()+1,
                                  RjPlus1.lastCol()-TransInvRight.numCols()+1);
    RjPlus1(TransInvRight) = TransInvRight;

//
//  Setup inner band structure: page 25
//
    GeMatrix<T> ExtM0(_(-l2+1,pow2i<int>(j+1)-l1-1),
                      _(-l2+1,pow2i<int>(  j)-l1-1));
    for (int q=ExtM0.firstCol(); q<=ExtM0.lastCol(); ++q) {
        const int p0 = max(l1+2*q, ExtM0.firstRow());
        const int p1 = min(l2+2*q, ExtM0.lastRow());
        for (int p = p0; p<=p1; ++p) {
            ExtM0(p,q) = phiR.a(p-2*q);
        }
    }

//
//  Finalize M0
//
    GeMatrix<T>  Mj0;
    GeMatrix<T>  M0Tmp;

//  compute RjPlus1Inv = RjPlus1^{-1}
    GeMatrix<T>  RjPlus1Inv = RjPlus1;
    inv(RjPlus1Inv);

    blas::mm(NoTrans, NoTrans, One, RjPlus1Inv, ExtM0, Zero, M0Tmp);
    blas::mm(NoTrans, NoTrans, One, M0Tmp, R, Zero, Mj0);
    blas::scal(Const<T>::R_SQRT2, Mj0);

    const int I0 = (_bc(0)==NoBC)
                 ? Mj0.firstRow()
                 : Mj0.firstRow() + 1;
    const int J0 = (_bc(0)==NoBC)
                 ? Mj0.firstCol()
                 : Mj0.firstCol() + 1;

    const int I1 = (_bc(1)==NoBC)
                 ? Mj0.lastRow()
                 : Mj0.lastRow() - 1;
    const int J1 = (_bc(1)==NoBC)
                 ? Mj0.lastCol()
                 : Mj0.lastCol() - 1;

    GeMatrix<T> _Mj0 = Mj0(_(I0,I1), _(J0,J1));
    M0 = RefinementMatrix<T,Interval,Dijkema>(d-1-_bc(0), d-1-_bc(1),
                                              _Mj0, j, j);

//
//  Set current level
//
    M0.setLevel(_j);
}

} // namespace lawa

#endif // LAWA_CONSTRUCTIONS_INTERVAL_DIJKEMA_PRIMAL_MRA_TCC
