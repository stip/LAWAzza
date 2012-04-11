/*
 LAWA - Library for Adaptive Wavelet Applications.
 Copyright (C) 2008-2012 Sebastian Kestler, Kristina Steih,
                         Alexander Stippler, Schalk.

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

#ifndef LAWA_CONSTRUCTIONS_REALLINE_PRIMAL_BSPLINE_H
#define LAWA_CONSTRUCTIONS_REALLINE_PRIMAL_BSPLINE_H 1

#include <lawa/flensforlawa.h>

#include <lawa/constructions/basisfunction.h>
#include <lawa/constructions/bspline.h>
//TODO #include <lawa/constructions/support.h>

namespace lawa {

template <typename _T>
struct BSpline<_T,Primal,R,CDF>
    : public BasisFunction<_T,Primal,R,CDF>
{
    typedef _T T;
    static const FunctionSide Side = Primal;
    static const DomainType Domain = R;
    static const Construction Cons = CDF;

    BSpline(int _d);

//TODO    BSpline(const MRA<T,Primal,R,CDF> &mra);

    virtual
    ~BSpline();

    const T
    operator()(T x, int j, Integer k, unsigned short deriv) const;

    const Support<T>
    support(int j, Integer k) const;

//TODO    const DenseVector<T>
//TODO    singularSupport(int j, Integer k) const;

    const T
    tic(int j) const;

    const DenseVector<T> &
    mask() const;

    const int d, mu;
    const Integer l1, l2;
    const DenseVector<T> a;
};

} // namespace lawa

#endif // LAWA_CONSTRUCTIONS_REALLINE_PRIMAL_BSPLINE_H
