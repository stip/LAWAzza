/*
  This file is part of LAWA - Library for Adaptive Wavelet Applications.
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

#ifndef LAWA_CONSTRUCTIONS_REALLINE_DUAL_BSPLINE_H
#define LAWA_CONSTRUCTIONS_REALLINE_DUAL_BSPLINE_H 1

#include <lawa/flensforlawa.h>

#include <lawa/constructions/basisfunction.h>
#include <lawa/constructions/bspline.h>
#include <lawa/constructions/support.h>

namespace lawa {

template <typename _T>
struct BSpline<_T,Dual,R,CDF>
    : public BasisFunction<_T,Dual,R,CDF>
{
    typedef _T T;
    static const FunctionSide Side = Dual;
    static const DomainType Domain = R;
    static const Construction Cons = CDF;

    BSpline(int _d, int _d_);

    virtual
    ~BSpline();

    const T
    operator()(T x, int j, Integer k, unsigned short deriv=0) const;

    const Support<T>
    support(int j, Integer k) const;
    
    const DenseVector<T> &
    mask() const;

    const int d, d_, mu;
    const int l1_, l2_;
    const DenseVector<T> a_;

    static unsigned int resolution;
};

template <typename T>
BSpline<T,Dual,R,CDF>
N_(int d, int d_);

template <typename T>
DenseVector<T>
N_1_1();

template <typename T>
DenseVector<T>
N_1_3();

template <typename T>
DenseVector<T>
N_1_5();

template <typename T>
DenseVector<T>
N_1_7();

template <typename T>
DenseVector<T>
N_1_9();

template <typename T>
DenseVector<T>
N_2_2();

template <typename T>
DenseVector<T>
N_2_4();

template <typename T>
DenseVector<T>
N_2_6();

template <typename T>
DenseVector<T>
N_2_8();

template <typename T>
DenseVector<T>
N_2_10();

template <typename T>
DenseVector<T>
N_3_3();

template <typename T>
DenseVector<T>
N_3_5();

template <typename T>
DenseVector<T>
N_3_7();

template <typename T>
DenseVector<T>
N_3_9();

template <typename T>
DenseVector<T>
N_4_4();

template <typename T>
DenseVector<T>
N_4_6();

template <typename T>
DenseVector<T>
N_4_8();

template <typename T>
DenseVector<T>
N_4_10();

template <typename T>
DenseVector<T>
N_5_5();

template <typename T>
DenseVector<T>
N_5_7();

template <typename T>
DenseVector<T>
N_5_9();

template <typename T>
DenseVector<T>
N_6_6();

template <typename T>
DenseVector<T>
N_6_8();

template <typename T>
DenseVector<T>
N_6_10();

template <typename T>
DenseVector<T>
N_7_7();

template <typename T>
DenseVector<T>
N_7_9();

template <typename T>
DenseVector<T>
N_8_8();

template <typename T>
DenseVector<T>
N_8_10();

template <typename T>
DenseVector<T>
N_9_9();

template <typename T>
DenseVector<T>
N_10_10();

} // namespace lawa

#endif // LAWA_CONSTRUCTIONS_REALLINE_DUAL_BSPLINE_H

