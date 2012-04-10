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

#ifndef LAWA_TYPEDEFS_H
#define LAWA_TYPEDEFS_H 1

//TODO
#include <flens/flens.h>

namespace lawa {

typedef flens::IndexOptions<>::IndexType Integer;

template <typename T>
    using DenseVector = flens::DenseVector<flens::Array<T> >;

template <typename T>
    using GeMatrix = flens::GeMatrix<flens::FullStorage<T,flens::ColMajor> >;

} // namespace lawa

#endif // LAWA_TYPEDEFS_H

