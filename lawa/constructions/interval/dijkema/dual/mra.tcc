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

#ifndef LAWA_CONSTRUCTIONS_INTERVAL_DIJKEMA_DUAL_MRA_TCC
#define LAWA_CONSTRUCTIONS_INTERVAL_DIJKEMA_DUAL_MRA_TCC

#include <algorithm>
#include <cassert>
#include <cmath>

#include <lawa/auxiliary/auxiliary.h>
#include <lawa/constructions/interval/dijkema/dual/selectboundaryfunctions.h>
#include <lawa/constructions/interval/dijkema/primal/splinehelper.h>
#include <lawa/constructions/realline/dual/bspline.h>
#include <lawa/flensforlawa.h>
#include <lawa/math/math.h>

namespace lawa {

using namespace std;

template <typename T>
MRA<T,Dual,Interval,Dijkema>::MRA(int _d, int _d_, int j)
    : d(_d), d_(_d_), mu(d&1),
      min_j0(iceil<int>(log2(d+2*d_-3))),
      j0(j==-1 ? min_j0 : j),
      _bc(2,0), _j(j0),
      phi_(*this)
{
    assert(d>1);
    assert((d+d_)%2==0);
    assert(_j>=min_j0);

    _calcM0_();
}

template <typename T>
MRA<T,Dual,Interval,Dijkema>::~MRA()
{
}

//--- cardinalities of whole, left, inner, right index sets. -------------------

template <typename T>
const int
MRA<T,Dual,Interval,Dijkema>::cardI_(int j) const
{
    assert(j>=min_j0);
    return pow2i<T>(j) + d - 1 - (_bc(0)+_bc(1));
}

template <typename T>
const int
MRA<T,Dual,Interval,Dijkema>::cardI_L() const
{
    return d + d_ - 2 -_bc(0);
}

template <typename T>
const int
MRA<T,Dual,Interval,Dijkema>::cardI_I(int j) const
{
    assert(j>=min_j0);
    return pow2i<T>(j) -d -2*d_ + 3;
}

template <typename T>
const int
MRA<T,Dual,Interval,Dijkema>::cardI_R() const
{
    return d + d_ - 2 - _bc(1);
}

//--- ranges of whole, left, inner, right index sets. --------------------------

template <typename T>
const Range<int>
MRA<T,Dual,Interval,Dijkema>::rangeI_(int j) const
{
    assert(j>=min_j0);
    return Range<int>(1 + _bc(0), pow2i<T>(j) + d - 1 - _bc(1));
}

template <typename T>
const Range<int>
MRA<T,Dual,Interval,Dijkema>::rangeI_L() const
{
    return Range<int>(1 + _bc(0), d + d_ - 2);
}

template <typename T>
const Range<int>
MRA<T,Dual,Interval,Dijkema>::rangeI_I(int j) const
{
    assert(j>=min_j0);
    return Range<int>(d + d_ - 1, pow2i<T>(j) - d_ + 1);
}

template <typename T>
const Range<int>
MRA<T,Dual,Interval,Dijkema>::rangeI_R(int j) const
{
    assert(j>=min_j0);
    return Range<int>(pow2i<T>(j) - d_ + 2, pow2i<T>(j) + d - 1 - _bc(1));
}

template <typename T>
const int
MRA<T,Dual,Interval,Dijkema>::level() const
{
    return _j;
}

template <typename T>
void
MRA<T,Dual,Interval,Dijkema>::setLevel(int j) const
{
    assert(j>=min_j0);
    _j = j;
    M0_.setLevel(_j);
}

template <typename T>
template <BoundaryCondition LeftBC, BoundaryCondition RightBC>
void
MRA<T,Dual,Interval,Dijkema>::enforceBoundaryConditions()
{
    _bc(Left)  = LeftBC;
    _bc(Right) = RightBC;

    _calcM0_();
}

template <typename T>
BoundaryCondition
MRA<T,Dual,Interval,Dijkema>::boundaryCondition(BoundarySide side) const
{
    assert(side==Left || side==Right);

    return _bc(side);
}

template <typename T>
void
MRA<T,Dual,Interval,Dijkema>::_calcC_LR(BoundarySide side, GeMatrix<T> &C_LR)
{
    using std::abs;
    using std::swap;

    const int  l1  = (mu-d)/2;
    const int  l2  = (mu+d)/2;
    const int  l1_ = l1-d_+1;
    const int  l2_ = l2+d_-1;

    const BoundaryCondition bc = boundaryCondition(side);

    BSpline<T,Dual,R,CDF>   phi_R(d, d_);

    const int extra = std::max(2-d+bc,0);

    GeMatrix<T>  A_(_(-l2_+1,-l1_-1+extra),
                    _(-l2_+1,-l1_-1+extra));

    for (int k=A_.firstRow(); k<=A_.lastRow(); ++k) {
        for (int m=A_.firstCol(); m<=A_.lastCol(); ++m) {
            if (l1_<=k-2*m && k-2*m<=l2_) {
                A_(k,m) = phi_R.a_(k-2*m);
            }
        }
    }

    GeMatrix<T>     VL;
    GeMatrix<T>     VR(A_.numRows(), A_.numRows());
    DenseVector<T>  r(VR.numRows()), i(VR.numRows());
    DenseVector<T>  work(1000);

    A_.changeIndexBase(1,1);
    flens::lapack::ev(false, true, A_, r, i, VL, VR, work);

    DenseVector<bool>  r_skip(r.length(), r.firstIndex());
    DenseVector<bool>  r_used(r.length(), r.firstIndex());
    DenseVector<int>   r_card(d_, r.firstIndex());
    int                numSkip = 0;
    int                numUsed = 0;
    bool               swapped;
    
    
    cerr << "r =  " << r << endl;
    cerr << "VR = " << VR << endl;


    r_skip = false;
    r_used = false;

    const T tol = 1e-7;

//
//  Skip eigenvalue if corresponding eigenvector is zero
//
    for (int i=1; i<=r.lastIndex(); ++i) {
        T sum = flens::blas::asum(VR(_,i));
        if (sum<tol) {
            r_skip(i) = true;
            ++numSkip;
        }
    }

//
//  Skip eigenvalue if its multitude is more than one
//
    for (int i=r.firstIndex(); i<=r.lastIndex(); ++i) {
        if (r_skip(i)) {
            continue;
        }
        for (int j=i+1; j<=r.lastIndex(); ++j) {
            if (r_skip(j)) {
                continue;
            }
            if (abs(r(i)-r(j))<tol) {
                r_skip(j) = true;
                ++numSkip;
            }
        }
    }

//
//  Eigenvalues of the form 2^{-j} always will be used
//
    for (int i=r.firstIndex(), I=i; i<=r.lastIndex(); ++i) {
        if (r_skip(i)) {
            continue;
        }
        if (r(i)>=0) {
            const T   fExp = log2(r(i));
            const int iExp = round(fExp);
            const T   diff = abs(fExp - iExp);
            
            if ((iExp<=0) && (-d_ < iExp) && (diff<tol)) {
                r_used(i) = true;
                r_card(I) = i;
                ++numUsed;
                ++I;
            }
        }
    }

//
//  Bubble sort the cardinal eigenvalues
//
    do {
        swapped = false;
        for (int i=r_card.firstIndex(); i<=r_card.lastIndex()-1; ++i) {
            if (r(r_card(i))<r(r_card(i+1))) {
                swap(r_card(i), r_card(i+1));
                swapped = true;
            }
        }
    } while (swapped);

//
//  Get indices of free (i.e. not used and not skipped) eigenvalues
//
    DenseVector<int>  r_free(r.length()-numUsed-numSkip, r.firstIndex());
    
    for (int i=r.firstIndex(), I=i; i<=r.lastIndex(); ++i) {
        if (r_skip(i) || r_used(i)) {
            continue;
        }
        r_free(I) = i;
        ++I;
    }

//
//  Bubble sort the free eigenvalues
//
    do {
        swapped = false;
        for (int i=r_free.firstIndex(); i<=r_free.lastIndex()-1; ++i) {
            if (r(r_free(i))>r(r_free(i+1))) {
                swap(r_free(i), r_free(i+1));
                swapped = true;
            }
        }
    } while (swapped);

    DenseVector<T>  indices;
    defaultSelectBoundaryFunctions(*this, side, indices);
    assert(d_+indices.length()==d+d_-2);


    C_LR.resize(d+2*d_-3+extra, d+d_-2+extra-bc);

//
//  Copy eigenvectors for cardinal eigenvalues
//
    const int i0 = r_card.firstIndex()+bc;
    const int i1 = r_card.lastIndex();
    
    cerr << "r_card = "  << r_card << endl;
    cerr << "r = "  << r << endl;
    cerr << "extra = "  << extra << endl;
    
    for (int i=i0, I=C_LR.firstCol(); i<=i1; ++i, ++I) {
        C_LR(_,I) = VR(_,r_card(i));
        cerr << "r_card(" << i << ") = " << r_card(i) << endl;
    }

//
//  Copy selection of the remaining eigenvectors
//
    const int I0 = C_LR.firstCol()+(i1-i0+1);

    for (int i=r_free.firstIndex(), I=I0; I<=C_LR.lastCol(); ++i, ++I) {
        C_LR(_,I) = VR(_,r_free(i));
    }

    if (side==Right) {
        GeMatrix<T> tmp = C_LR;
        arrow(tmp, C_LR);
    }

}

template <typename T>
void
MRA<T,Dual,Interval,Dijkema>::_calcC_L(GeMatrix<T> &C_L)
{
    _calcC_LR(Left, C_L);
}

template <typename T>
void
MRA<T,Dual,Interval,Dijkema>::_calcC_R(GeMatrix<T> &C_R)
{
    _calcC_LR(Right, C_R);
}


template <typename T>
void
MRA<T,Dual,Interval,Dijkema>::_calcM0_()
{
    using std::max;
    using std::min;
    using cxxblas::NoTrans;
    using cxxblas::Trans;

    const int cons_j = ((d==2) && ((d_==2) || (d_==4))) ? min_j0+1 : min_j0;

    const int  l1  = (mu-d)/2;
    const int  l2  = (mu+d)/2;
    const int  l1_ = l1-d_+1;
    const int  l2_ = l2+d_-1;

    BSpline<T,Dual,R,CDF>  phi_R(d,d_);

    GeMatrix<T>  C_L, C_R;
    
    _calcC_L(C_L);
    _calcC_R(C_R);
    
    cerr << "C_L = " << C_L << endl;
    cerr << "C_R = " << C_R << endl;

//
//  Page 25 (top)
//
    C_L.changeIndexBase(-l2_+1, -l2+1+_bc(0));

    GeMatrix<T>  R_init(cardI_I(cons_j)+2*(d+2*d_-3),
                        cardI_(cons_j),
                        C_L.firstRow(),
                        C_L.firstCol());

    C_R.changeIndexBase(R_init.lastRow()-C_R.numRows()+1,
                        R_init.lastCol()-C_R.numCols()+1);
    
    assert(R_init.numRows()==pow2i<int>(cons_j)+l2_-l1_-1);
    assert(R_init.numCols()==pow2i<int>(cons_j)+l2-l1-1-_bc(0)-_bc(1));
    assert(R_init.firstRow()==-l2_+1);
    assert(R_init.firstCol()==-l2+1+_bc(0));

    for (int c=R_init.firstCol(); c<=R_init.lastCol(); ++c) {
        R_init(c,c) = 1.;
    }

    R_init(C_L) = C_L;
    R_init(C_R) = C_R;


//    // cout << "R_init = " << R_init << endl;

    BSpline<T,Primal,R,CDF> phi(d);
    int kmin = -l2+1, mmin = -l2_+1, pmin = -l2+1, qmin = -l2_+1,
        kmax = -l1-1, mmax = -l1_-1, pmax = -l1-1, qmax = -l1_-1;
    GeMatrix<T>  C(_(pmin,pmax), _(qmin,qmax));
    int ZLength = C.numRows()*C.numCols();
    GeMatrix<T>  Z(_(pmin*(qmax-qmin+1)+qmin,
                     pmin*(qmax-qmin+1)+qmin+ZLength-1),
                   _(pmin*(qmax-qmin+1)+qmin, 
                     pmin*(qmax-qmin+1)+qmin+ZLength-1));
    for (int p=pmin; p<=pmax; ++p) {
        for (int q=qmin; q<=qmax; ++q) {
            for (int k=kmin; k<=kmax; ++k) {
                for (int m=mmin; m<=mmax; ++m) {
                    if (l1<=k-2*p && k-2*p<=l2 &&
                        l1_<=m-2*q && m-2*q<=l2_) {
                            Z(p*(qmax-qmin+1)+q, k*(mmax-mmin+1)+m) =
                                           0.5 * phi.a(k-2*p) * phi_R.a_(m-2*q);
                    }
                }
            }
        }
    }
    for (int i=Z.firstRow(); i<=Z.lastRow(); ++i) {
        Z(i,i) -= 1.;
    }
    Z *= -1;

    DenseVector<T> f(Z.numCols(), Z.firstCol());
    for (int p=pmin; p<=pmax; ++p) {
        for (int q=qmin; q<=qmax; ++q) {
            T sum = 0.;
            for (int k=max(-l2+1,l1+2*p); k<=l2+2*p; ++k) {
                for (int m=max(max(-l2_+1,l1_+2*q),l1-l2_+1+k);
                    m<=min(l2_+2*q, l2-l1_-1+k); ++m) {
                    if (!((kmin<=k && k<=kmax) && (mmin<=m && m<=mmax))) {
                        sum += 0.5 * phi.a(k-2*p) * phi_R.a_(m-2*q) * (k==m);
                    }
                }
            }
            f(p*(qmax-qmin+1)+q) = sum;
        }
    }
    solve(Z, f);

    for (int p=pmin; p<=pmax; ++p) {
        for (int q=qmin; q<=qmax; ++q) {
            C(p,q) = f(p*(qmax-qmin+1)+q);
        }
    }

    GeMatrix<T> ExtMass(pow2i<int>(cons_j)+d-1,
                        pow2i<int>(cons_j)+l2_-l1_-1,
                        -l2+1, -l2_+1);
    for (int r=ExtMass.firstRow(); r<=ExtMass.lastRow(); ++r) {
        ExtMass(r,r) = 1.;
    }
    ExtMass(C) = C;
    
    GeMatrix<T>  CR;
    arrow(C,CR);
    CR.changeIndexBase(ExtMass.lastRow()-C.numRows()+1,
                       ExtMass.lastCol()-C.numCols()+1);
    ExtMass(CR) = CR;

    GeMatrix<T>  ExtMassjPlus1(pow2i<T>(cons_j+1)+d-1, 
                               pow2i<T>(cons_j+1)+l2_-l1_-1,
                               -l2+1, -l2_+1);
    for (int r=ExtMassjPlus1.firstRow(); r<=ExtMassjPlus1.lastRow(); ++r) {
        ExtMassjPlus1(r,r) = 1.;
    }
    ExtMassjPlus1(C) = C;
    CR.changeIndexBase(ExtMassjPlus1.lastRow()-C.numRows()+1,
                       ExtMassjPlus1.lastCol()-C.numCols()+1);
    ExtMassjPlus1(CR) = CR;


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
    GeMatrix<T>  R(_(1-l2, pow2i<int>(cons_j)-l1-1),
                   _(1-l2+_bc(0), pow2i<int>(cons_j)-l1-1-_bc(1)));

    //  cerr << R.rows() << "x" << R.cols() << endl;
    int rr = R.firstRow()+_bc(0);
    for (int c=R.firstCol(); c<=R.lastCol(); ++c, ++rr) {
        R(rr,c) = 1.;
    }

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
    GeMatrix<T>  RjPlus1(_(-l2+1,pow2i<int>(cons_j+1)-l1-1),
                         _(-l2+1+_bc(0),pow2i<int>(cons_j+1)-l1-1-_bc(1)));
    RjPlus1.diag(0) = 1.;
    TransInvLeft.changeIndexBase(RjPlus1.firstRow(),RjPlus1.firstCol());
    RjPlus1(TransInvLeft) = TransInvLeft;
    TransInvRight.changeIndexBase(RjPlus1.lastRow()-TransInvRight.numRows()+1,
                                  RjPlus1.lastCol()-TransInvRight.numCols()+1);
    RjPlus1(TransInvRight) = TransInvRight;

    GeMatrix<T>  Mass, MassTmp;

    blas::mm(NoTrans, NoTrans, 1., ExtMass, R_init, 0., MassTmp);
    blas::mm(Trans, NoTrans, 1., R, MassTmp, 0., Mass);

    GeMatrix<T>  InvMass = Mass;

    inv(InvMass);

    // std::cerr << "R_init = " << R_init << std::endl;
    // std::cerr << "InvMass = " << InvMass << std::endl;


    GeMatrix<T>   R_(R_init.numRows(), R_init.numCols(),
                     R_init.firstRow(), R_init.firstCol());
    blas::mm(NoTrans, NoTrans, 1., R_init, InvMass, 0., R_);


    R_Left.resize(0,0);
    R_Right.resize(0,0);

    R_Left = R_(_(C_L.firstRow(),C_L.lastRow()),
                _(C_L.firstCol(),C_L.lastCol()+_bc(0)));
    R_Right = R_(_(C_R.firstRow(),C_R.lastRow()),
                 _(C_R.firstCol()-_bc(1),C_R.lastCol()));


    // std::cerr << "R_ = " << R_ << std::endl;

    GeMatrix<T>   ExtM0_(_(-l2_+1,pow2i<int>(cons_j+1)-l1_-1),
                         _(-l2_+1,pow2i<int>(cons_j)-l1_-1));
    for (int q=ExtM0_.firstCol(); q<=ExtM0_.lastCol(); ++q) {
        for (int p=max(l1_+2*q, ExtM0_.firstRow());
             p<=min(l2_+2*q, ExtM0_.lastRow()); ++p) {
            ExtM0_(p,q) = phi_R.a_(p-2*q);
        }
    }

    GeMatrix<T>  Mj0_, M0_Tmp;
    blas::mm(Trans, NoTrans, 1., RjPlus1, ExtMassjPlus1, 0., Mj0_);
    blas::mm(NoTrans, NoTrans, 1., Mj0_, ExtM0_, 0., M0_Tmp);

    Mj0_.resize(0,0);
    blas::mm(NoTrans, NoTrans, 1., M0_Tmp, R_, 0., Mj0_);

    Mj0_.changeIndexBase(1+_bc(0),1+_bc(1));
    blas::scal(Const<T>::R_SQRT2, Mj0_);
    R_Left.changeIndexBase(1,1+_bc(0));
    R_Right.changeIndexBase(1,1);
    
    //cerr << Mj0_ << endl;
    Mj0_.changeIndexBase(1,1);
    
    M0_ = RefinementMatrix<T,Interval,Dijkema>(d+d_-2-_bc(0), d+d_-2-_bc(1),
                                               Mj0_, min_j0, cons_j);
    setLevel(_j);
}

} // namespace lawa

#endif // LAWA_CONSTRUCTIONS_INTERVAL_DIJKEMA_DUAL_MRA_TCC