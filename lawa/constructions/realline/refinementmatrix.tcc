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

#include <cassert>

#include <lawa/flensforlawa.h>
#include <lawa/constructions/realline/refinementmatrix.h>
#include <lawa/math/const.h>

namespace flens {

template <typename T>
template <FunctionSide Side>
RefinementMatrix<T,R,CDF>::RefinementMatrix(const BSpline<T,Side,R,CDF> &spline)
      : band(Const<T>::R_SQRT2 * spline.mask())
{
}

template <typename T>
template <FunctionSide Side>
RefinementMatrix<T,R,CDF>::RefinementMatrix(const Wavelet<T,Side,R,CDF> &wavelet)
      : band(Const<T>::R_SQRT2 * wavelet.mask())
{
}


//------------------------------------------------------------------------------

/*
template <typename X, typename Y>
void
mv(cxxblas::Transpose                                       transA,
   typename X::ElementType                                  alpha,
   const RefinementMatrix<typename X::ElementType,R,CDF>    &A,
   const flens::DenseVector<X>                              &x, 
   typename X::ElementType                                  beta,
   flens::DenseVector<Y>                                    &y)
{
    using namespace cxxblas;
    using std::min;

    typedef typename X::ElementType T;
    typedef typename X::IndexType   IndexType;

    const lawa::DenseVector<T> &a = A.band;
    IndexType l1 = x.firstIndex();
    IndexType l2 = x.lastIndex();

    if (transA==cxxblas::NoTrans) {
        if (beta==0) {
            y.resize(2*(l2-l1)+1, 2*l1) || y.fill();
        } else {
            assert(y.length()==2*(l2-l1));
            assert(y.firstIndex()==2*l1);
        }

        const Underscore<IndexType>  _;

        for (IndexType l=l1, L=2*l1; l<=l2; ++l, ++L) {
            const IndexType Ly = min(2*l2, L+a.length()-1);
            const IndexType La = min(l2, 2*(l2-l1)-l-1);

            auto        _y = y(_(L,Ly));
            const auto  _a = a(_(l1,La));

            _y += alpha*x(l) * _a;

        }

    } else { // (transA==cxxblas::Trans)
        assert(0);
    }
}
*/


/*
//
// Urban: p25-26
//
template <typename X, typename Y>
void
mv(cxxblas::Transpose                                       transA,
   typename X::ElementType                                  alpha,
   const RefinementMatrix<typename X::ElementType,R,CDF>    &A,
   const flens::DenseVector<X>                              &x, 
   typename X::ElementType                                  beta,
   flens::DenseVector<Y>                                    &y)
{
    using namespace cxxblas;
    
    typedef typename X::ElementType T;



    assert(alpha==1.);
    assert(x.engine().stride()==1);

    const lawa::DenseVector<T> &a = A.band;
    int lx = x.length();
    int la = a.length();
    
    if (transA==cxxblas::NoTrans) {
        if (beta==0) {
            y.resize(2*lx,x.firstIndex()-a.firstIndex()) || y.fill();
        } else {
            assert(y.length()==2*lx);
            y.changeIndexBase(x.firstIndex()-a.firstIndex());
        }
        const T *xp = x.engine().data();
        for (int k=0; k<lx; ++k, ++xp) {
            int mMin = a.firstIndex() + 2*k;
            const T *abegin = a.engine().data();
                  T *ybegin = y.engine().data();
            cxxblas::axpy(la, *xp, abegin, 1, ybegin+mMin, 1);
        }
    } else { // (transA==cxxblas::Trans)
        if (beta==0) {
            y.engine().resize(lx/2,x.firstIndex()+a.firstIndex()) || y.engine().fill();
        } else {
            assert(y.length()==lx/2);
            y.engine().changeIndexBase(x.firstIndex()+a.firstIndex());
        }
        T *iter = y.engine().data();
        for (int m=0; m<lx/2; ++m, ++iter) {
            int kMin = 2*m;
            const T *abegin = a.engine().data();
            T dotValue;
            cxxblas::dot(la, abegin, 1, x.engine().data()+kMin, 1, dotValue);
            *iter += dotValue;
        }
    }
}
*/

template <typename VY0, typename VMASK, typename BETA, typename VY1>
void
compose(const DenseVector<VY0> &y0, const DenseVector<VMASK> &mask, 
        BETA beta, DenseVector<VY1> &y1)
{
    using std::max;
    using std::min;
    using std::sqrt;

    typedef typename DenseVector<VY1>::ElementType  T;
    
    const int l1 = mask.firstIndex();
    const int l2 = mask.lastIndex();

    const int m1 = y0.firstIndex();
    const int m2 = y0.lastIndex();

    const int M1 = 2*m1+l1;
    const int M2 = 2*m2+l2;

#   ifndef NDEBUG
    assert (beta==BETA(0) || y1.length()==M2-M1+1);
#   endif

    std::cerr << "l1 = " << l1 << ", l2 = " << l2 << std::endl;
    std::cerr << "M1 = " << M1 << ", M2 = " << M2 << std::endl;


    if (beta==BETA(0)) {
        y1.resize(M2-M1+1, M1) || y1.fill(T(0));
    }

    for (int l=M1; l<=M2; ++l) {
        const int m1_ = max(m1, iceil<int>((l-l2)/2.0));
        const int m2_ = min(m2, ifloor<int>((l-l1)/2.0));

        for (int m=m1_; m<=m2_; ++m) {
            y1(l) += y0(m) * mask(l-2*m);
        }
    }
}

template <typename X, typename Y>
void
mv(cxxblas::Transpose                                       transA,
   typename X::ElementType                                  alpha,
   const RefinementMatrix<typename X::ElementType,R,CDF>    &A,
   const flens::DenseVector<X>                              &x, 
   typename X::ElementType                                  beta,
   flens::DenseVector<Y>                                    &y)
{
    using namespace cxxblas;
    using std::max;
    using std::min;
    using std::sqrt;

    typedef typename X::ElementType T;


    assert(alpha==1.);
    assert(x.engine().stride()==1);

    const lawa::DenseVector<T> &a = A.band;
    
    const int l1 = a.firstIndex();
    const int l2 = a.lastIndex();
    
    const int L1 = x.firstIndex();
    const int L2 = x.lastIndex();
    
    const int m1 = iceil<int>((L1 - l2)/2.);
    const int m2 = ifloor<int>((L2 - l1)/2.);
    
    if (transA==cxxblas::NoTrans) {
        assert(0);
    } else { // (transA==cxxblas::Trans)
        if (beta==0) {
            y.resize(m2-m1+1, m1) || y.engine().fill();
        } else {
            assert(y.length()==m2-m1+1);
            y.engine().changeIndexBase(m1);
        }
        for (int m=m1; m<=m2; ++m) {
            y(m) = 0;
            for (int l=max(l1+2*m,L1); l<=min(l2+2*m,L2); ++l) {
                y(m) += a(l-2*m) * x(l);
            }
            //y(m) /= sqrt(2);
        }
    }
}


} // namespace flens

