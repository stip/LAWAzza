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

#ifndef LAWA_INTEGRALS_INTEGRAL_TCC
#define LAWA_INTEGRALS_INTEGRAL_TCC 1

#include <cassert>
#include <lawa/constructions/basisfunction.h>
#include <lawa/constructions/bspline.h>
#include <lawa/constructions/support.h>
#include <lawa/math/math.h>
#include <lawa/typetraits.h>

namespace lawa {
   
//--- primal * primal or orthogonal * orthogonal
template <typename First, typename Second>
typename RestrictTo<BothPrimal<First,Second>::value
                    or BothOrthogonal<First,Second>::value, typename First::T>::Type
_integrate(const Integral<Gauss,First,Second> &integral)
{
    typedef typename First::T T;

    const typename First::BasisFunctionType &first = integral.first.generator(integral.e1);
    const typename Second::BasisFunctionType &second = integral.second.generator(integral.e2);

    // the (minimal) width of the polynomial pieces.
    T unit = std::min(first.tic(integral.j1), second.tic(integral.j2));
    integral.quadrature.setOrder((integral.first.d-integral.deriv1+integral.second.d-integral.deriv2)/2+1);

    T ret = 0.;
    Support<T> common;
    if (overlap(first.support(integral.j1,integral.k1),
               second.support(integral.j2,integral.k2),common)) {
        T a = common.l1;
        for (T b=a+unit; b<=common.l2; b+=unit) {
            ret += integral.quadrature(a,b);
            a = b;
        }
    }

    return ret;
}

/*
//--- (primal * dual) or (dual * primal) or (dual * dual) or (dual * orthogonal) or (orthogonal * dual)
template <QuadratureType Quad, typename First, typename Second>
typename RestrictTo<IsDual<First>::value or IsDual<Second>::value, typename First::T>::Type
_integrate(const Integral<Quad,First,Second> &integral)
{
    typedef typename First::T T;
    const typename First::BasisFunctionType &first = integral.first.generator(integral.e1);
    const typename Second::BasisFunctionType &second = integral.second.generator(integral.e2);

    Support<T> common;
    if (overlap(first.support(integral.j1,integral.k1),
                second.support(integral.j2,integral.k2),
                common)) {
        return integral.quadrature(common.l1, common.l2);
    } else {
        return 0;
    }
}
    
//--- function * primal/dual/orthogonal
template <QuadratureType Quad, typename First, typename Second>
typename First::T
_integrate_f1(const IntegralF<Quad,First,Second> &integral)
{
    typedef typename First::T T;
    const typename First::BasisFunctionType &first = integral.first.generator(integral.e1);
    Integer p = integral.function.singularPoints.length();
    DenseVector<T> singularPoints;
    if (IsPrimal<First>::value or IsOrthogonal<First>::value) {
        if (p>0) {
            const DenseVector<T> & firstSingularPoints = first.singularSupport(integral.j1,integral.k1);
            Integer m = firstSingularPoints.length();
            singularPoints.engine().resize(m+p);

            std::merge(firstSingularPoints.engine().data(),
                       firstSingularPoints.engine().data() + m,
                       integral.function.singularPoints.engine().data(),
                       integral.function.singularPoints.engine().data() + p,
                       singularPoints.engine().data());
        } else {
            singularPoints = first.singularSupport(integral.j1,integral.k1);
        }
    } else {
        singularPoints.engine().resize(2L);
        Support<T> supp = first.support(integral.j1,integral.k1);
        singularPoints = supp.l1, supp.l2;
    }

    T ret = 0.0;
    for (Integer i=singularPoints.firstIndex(); i<singularPoints.lastIndex(); ++i) {
        ret += integral.quadrature(singularPoints(i),singularPoints(i+1));
    }
    return ret;
}

//--- function * primal/dual/orthogonal * primal/dual/orthogonal)
template <QuadratureType Quad, typename First, typename Second>
typename First::T
_integrate_f2(const IntegralF<Quad,First,Second> &integral)
{
    typedef typename First::T T;
    const typename First::BasisFunctionType &first = integral.first.generator(integral.e1);
    const typename Second::BasisFunctionType &second = integral.second.generator(integral.e2);
    const Function<T> & function = integral.function;

    DenseVector<T> singularPoints;
    Integer p = function.singularPoints.length();

    if (!BothDual<First,Second>::value) { // we do have (additional) singular points
        // merge singular points of bsplines/wavelets and function to one list.
        // -> implicit assumption: singular points are sorted!
        if (BothPrimal<First,Second>::value or BothOrthogonal<First,Second>::value) {
            const DenseVector<T> & firstSingularPoints = first.singularSupport(integral.j1,integral.k1);
            const DenseVector<T> & secondSingularPoints = second.singularSupport(integral.j2,integral.k2);
            int m = firstSingularPoints.length();
            int n = secondSingularPoints.length();
            singularPoints.engine().resize(m+n+p);

            std::merge(firstSingularPoints.engine().data(),
                       firstSingularPoints.engine().data() + m,
                       secondSingularPoints.engine().data(),
                       secondSingularPoints.engine().data() + n,
                       singularPoints.engine().data());
            if (p>0) {
                singularPoints(_(m+n+1,m+n+p)) = function.singularPoints;
                std::inplace_merge(singularPoints.engine().data(),
                                   singularPoints.engine().data() + m + n,
                                   singularPoints.engine().data() + m + n + p);
            }
        } else {
            if (IsPrimal<First>::value or IsOrthogonal<First>::value) {
                if (p>0) {
                    const DenseVector<T> & firstSingularPoints = first.singularSupport(integral.j1,integral.k1);
                    int m = firstSingularPoints.length();
                    singularPoints.engine().resize(m+p);

                    std::merge(firstSingularPoints.engine().data(),
                               firstSingularPoints.engine().data() + m,
                               function.singularPoints.engine().data(),
                               function.singularPoints.engine().data() + p,
                               singularPoints.engine().data());
                } else {
                    singularPoints = first.singularSupport(integral.j1,integral.k1);
                }
            } else {
                if (p>0) {
                    const DenseVector<T> & secondSingularPoints = second.singularSupport(integral.j2,integral.k2);
                    int n = secondSingularPoints.length();
                    singularPoints.engine().resize(n+p);

                    std::merge(secondSingularPoints.engine().data(),
                               secondSingularPoints.engine().data() + n,
                               function.singularPoints.engine().data(),
                               function.singularPoints.engine().data() + p,
                               singularPoints.engine().data());
                } else {
                    singularPoints = second.singularSupport(integral.j2,integral.k2);
                }
            }
        }

    }

    T ret = 0.0;
    for (int i=singularPoints.firstIndex(); i<singularPoints.lastIndex(); ++i) {
        ret += integral.quadrature(singularPoints(i),singularPoints(i+1));
    }

    return ret;
}
*/
//-----------------------------------------------------------------------------

//--- primal/dual/orthogonal * primal/dual/orthogonal
template <QuadratureType Quad, typename First, typename Second>
typename First::T
_integrand(const Integral<Quad,First,Second> &integral, typename First::T x)
{
    const typename First::BasisFunctionType &first = integral.first.generator(integral.e1);
    const typename Second::BasisFunctionType &second = integral.second.generator(integral.e2);
    return first(x,integral.j1,integral.k1,integral.deriv1) * second(x,integral.j2,integral.k2,integral.deriv2);
}
/*
//--- function * primal/dual/orthogonal
template <QuadratureType Quad, typename First, typename Second>
typename RestrictTo<PrimalOrDualOrOrthogonal<First>::value, typename First::T>::Type
_integrand_f1(const IntegralF<Quad,First,Second> &integral, typename First::T x)
{
    const typename First::BasisFunctionType &first = integral.first.generator(integral.e1);
    return integral.function(x) * first(x,integral.j1,integral.k1,integral.deriv1);
}

//--- function * primal/dual/orthogonal * primal/dual/orthogonal
template <QuadratureType Quad, typename First, typename Second>
typename RestrictTo<PrimalOrDualOrOrthogonal<First>::value 
                    and PrimalOrDualOrOrthogonal<Second>::value, typename First::T>::Type
_integrand_f2(const IntegralF<Quad,First,Second> &integral, typename First::T x)
{
    const typename First::BasisFunctionType &first = integral.first.generator(integral.e1);
    const typename Second::BasisFunctionType &second = integral.second.generator(integral.e2);
    return integral.function(x) * first(x,integral.j1,integral.k1,integral.deriv1)
                                * second(x,integral.j2,integral.k2,integral.deriv2);
}
*/
//------------------------------------------------------------------------------

template <QuadratureType Quad, typename First, typename Second>
Integral<Quad,First,Second>::Integral(const First &_first, 
                                      const Second &_second)
    : first(_first), second(_second), quadrature(*this)
{
}

template <QuadratureType Quad, typename First, typename Second>
typename First::T
Integral<Quad,First,Second>::operator()(int _j1, Integer _k1, XType _e1, int _deriv1, 
                                        int _j2, Integer _k2, XType _e2, int _deriv2) const
{
    j1 = _j1; k1 = _k1; e1 = _e1; deriv1 = _deriv1;
    j2 = _j2; k2 = _k2; e2 = _e2; deriv2 = _deriv2;
    return _integrate(*this);
}

template <QuadratureType Quad, typename First, typename Second>
typename First::T
Integral<Quad,First,Second>::integrand(typename First::T x) const
{
    return _integrand(*this, x);
}

//------------------------------------------------------------------------------

template <QuadratureType Quad, typename First, typename Second>
IntegralF<Quad,First,Second>::IntegralF(const Function<T> &_function,
                                        const First &_first)
    : function(_function), first(_first), second(_first), quadrature(*this),
      _f2(false)
{
}

template <QuadratureType Quad, typename First, typename Second>
IntegralF<Quad,First,Second>::IntegralF(const Function<T> &_function,
                                        const First &_first, 
                                        const Second &_second)
    : function(_function), first(_first), second(_second), quadrature(*this),
      _f2(true)

{
}

template <QuadratureType Quad, typename First, typename Second>
typename First::T
IntegralF<Quad,First,Second>::operator()(int _j, Integer _k, XType _e, int _deriv) const
{
    j1 = _j; k1 = _k; e1 = _e; deriv1 = _deriv;
    return _integrate_f1(*this);
}

template <QuadratureType Quad, typename First, typename Second>
typename First::T
IntegralF<Quad,First,Second>::operator()(int _j1, Integer _k1, XType _e1, int _deriv1,
                                         int _j2, Integer _k2, XType _e2, int _deriv2) const
{
    j1 = _j1; k1 = _k1; e1 = _e1; deriv1 = _deriv1;
    j2 = _j2; k2 = _k2; e2 = _e2; deriv2 = _deriv2;
    return _integrate_f2(*this);
}

template <QuadratureType Quad, typename First, typename Second>
typename First::T
IntegralF<Quad,First,Second>::integrand(typename First::T x) const
{
    return _f2 ? _integrand_f2(*this,x) : _integrand_f1(*this,x);
}

} // namespace lawa

#endif // LAWA_INTEGRALS_INTEGRAL_TCC
