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

#ifndef LAWA_MATH_CONST_H
#define LAWA_MATH_CONST_H 1

#include <cmath>

namespace lawa {

template <typename T>
struct Const
{
};

template <>
struct Const<double>
{
    static double EQUALITY_EPS;
    static double SQRT2;
    static double R_SQRT2;
};

} // namespace lawa

#endif // LAWA_MATH_CONST_H

