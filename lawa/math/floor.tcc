/*
 LAWA - Library for Adaptive Wavelet Applications.
 Copyright (C) 2008-2011 Sebastian Kestler, Kristina Steih,
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

#ifndef LAWA_MATH_FLOOR_TCC
#define LAWA_MATH_FLOOR_TCC 1

#include <cmath>

namespace lawa {

template <typename T>
const T
floor(double x)
{
    return static_cast<T>(std::floor(x));

/*  to be generalized i.e. more than 64 bit.
    register double twoTo52 = 4503599627370496.0;
    double c = (x>=0.0 ? -twoTo52 : twoTo52);
    double result = (x-c)+c;
    if (x<result) {
        result -= 1.0;
    }
    return static_cast<int>(result);
*/
}

} // namespace lawa

#endif // LAWA_MATH_FLOOR_TCC
