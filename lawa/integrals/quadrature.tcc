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

#ifndef LAWA_INTEGRALS_QUADRATURE_TCC
#define LAWA_INTEGRALS_QUADRATURE_TCC 1

#include <cassert>
#include <cmath>
#include <lawa/math/math.h>

namespace lawa {

//--- Gauss-Legendre Quadrature-------------------------------------------------

template <typename Integral>
Quadrature<Gauss,Integral>::Quadrature(const Integral &_integral)
    : integral(_integral), _order(-1)
{
    _legendre(_precalculatedOrder);
    setOrder(4);
}

template <typename Integral>
const typename Integral::T
Quadrature<Gauss,Integral>::operator()(T a, T b) const
{
    T ret = 0.0;
    for (int i=1; i<=_order; ++i) {
        ret +=   _weights(_order,i)
               * integral.integrand(0.5*(b-a)*_knots(_order,i)+0.5*(b+a));
    }
    ret *= 0.5*(b-a);

    return ret;
}

template <typename Integral>
void
Quadrature<Gauss,Integral>::setOrder(Integer order)
{
    assert(order>0);

    if (order>=_precalculatedOrder) {
        _legendre(order);
        _precalculatedOrder = order;
    }
    _order = order;
}

template <typename Integral>
const Integer
Quadrature<Gauss,Integral>::order() const
{
    return _order;
}

template <typename Integral>
void
Quadrature<Gauss,Integral>::_legendre(Integer order)
{
    using std::abs;
    using std::cos;
    
    T eps = 1e-15;
    _knots.engine().resize(order, order);
    _weights.engine().resize(order, order);

    T x1 = -1,
      x2 =  1;
    for (int k=1; k<=order; ++k) {
        int     m = (k+1)/2;
        T xm = 0.5 * (x2+x1),
          xl = 0.5 * (x2-x1);
        for (int i=1; i<=m; ++i) {
            T z = cos(M_PI*(i-0.25)/(k+0.5)),
              z1, pp;
            do {
                T p1 = 1.0,
                  p2 = 2.0;
                for (int j=1; j<=k; ++j) {
                    T p3 = p2;
                    p2 = p1;
                    p1 = ((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
                }
                pp = k * (z*p1-p2)/(z*z-1.0);
                z1 = z;
                z = z1-p1/pp;
            } while (abs(z-z1) > eps);
            _knots(k,i)     = xm - xl*z;
            _knots(k,k+1-i) = xm + xl*z;
            _weights(k,i)     = 2.0*xl/((1.0-z*z)*pp*pp);
            _weights(k,k+1-i) = _weights(k,i);
        }
    }
}

template <typename Integral>
GeMatrix<typename Integral::T>
Quadrature<Gauss,Integral>::_weights;

template <typename Integral>
GeMatrix<typename Integral::T>
Quadrature<Gauss,Integral>::_knots;

template <typename Integral>
Integer
Quadrature<Gauss,Integral>::_precalculatedOrder = 10;

//---  Trapezoidal rule -------------------------------------------------------
template <typename Integral>
Quadrature<Trapezoidal,Integral>::Quadrature(const Integral &_integral)
    : integral(_integral), _n(-1)
{
    setN(100);
}

template <typename Integral>
const typename Integral::T
Quadrature<Trapezoidal,Integral>::operator()(T a, T b) const
{
    T h = (b-a) / _n;
    T ret = .5 * h * integral.integrand(a);
    a += h;
    for (int i=1; i<_n; ++i, a+=h) {
        ret += h * integral.integrand(a);
    }
    ret += .5 * h * integral.integrand(b);

    return ret;
}

template <typename Integral>
const Integer
Quadrature<Trapezoidal,Integral>::n() const
{
    return _n;
}


template <typename Integral>
void
Quadrature<Trapezoidal,Integral>::setN(const Integer n)
{
    assert(n>0);

    _n = n;
}

} // namespace lawa

#endif // LAWA_INTEGRALS_QUADRATURE_TCC
