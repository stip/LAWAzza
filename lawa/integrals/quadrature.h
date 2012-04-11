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

#ifndef LAWA_INTEGRALS_QUADRATURE_H
#define LAWA_INTEGRALS_QUADRATURE_H 1

#include <lawa/enum.h>
#include <lawa/flensforlawa.h>

namespace lawa {

template <QuadratureType Quad, typename Integral>
class Quadrature
{
};

template <typename Integral>
class Quadrature<Gauss,Integral>
{
    public:
        typedef typename Integral::T T;

        Quadrature(const Integral &_integral);

        const T
        operator()(T a, T b) const;

        const Integer
        order() const;

        void
        setOrder(Integer order);

        const Integral &integral;
    private:
        static void
        _legendre(Integer order);

        Integer _order;
        static Integer _precalculatedOrder;
        static GeMatrix<T> _knots;
        static GeMatrix<T> _weights;
};

//------------------------------------------------------------------------------

template <typename Integral>
class Quadrature<Trapezoidal,Integral>
{
    public:
        typedef typename Integral::T T;
        
        Quadrature(const Integral &_integral);

        const T
        operator()(T a, T b) const;

        const Integer
        n() const;

        void
        setN(Integer n);

        const Integral &integral;
    private:
        Integer _n;
};

} // namespace lawa

#endif // LAWA_QUADRATURE_H
