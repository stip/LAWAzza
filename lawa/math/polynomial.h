/*
 LAWA - Library for Adaptive Wavelet Applications.
 Copyright (C) 2008,2009  Mario Rometsch, Alexander Stippler.

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

#ifndef LAWA_MATH_POLYNOMIAL_H
#define LAWA_MATH_POLYNOMIAL_H 1

#include <lawa/flensforlawa.h>

namespace lawa {

template <typename T>
class Polynomial
{
    public:
        Polynomial(int n=0);

        Polynomial(const DenseVector<T> &coefficients);

        const T &
        operator()(int n) const;

        T &
        operator()(int n);

        Polynomial<T> &
        operator+=(const Polynomial<T> &rhs);

        int
        degree() const;

    private:
        DenseVector<T> _coefficients;
};

template <typename T>
Polynomial<T>
operator*(const Polynomial<T> &lhs, const Polynomial<T> &rhs);

template <typename S, typename T>
Polynomial<T>
operator*(const S &lhs, const Polynomial<T> &rhs);

template <typename T>
Polynomial<T>
pow(const Polynomial<T> &p, int n);

} // namespace lawa

#endif // LAWA_MATH_POLYNOMIAL_H

