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

#ifndef LAWA_CONSTRUCTIONS_SUPPORT_H
#define LAWA_CONSTRUCTIONS_SUPPORT_H 1

#include <lawa/flensforlawa.h>

namespace lawa {

template <typename T>
struct Support
{
    Support();

    Support(const T &a, const T &b);

    ~Support();

    T
    length() const;

    T l1, l2;
};

template <typename T>
    bool
    inner(T x, const Support<T> &supp);

template <typename T>
    T
    overlap(const Support<T> &supp1, const Support<T> &supp2);

template <typename T>
    T
    overlap(const Support<T> &supp1, const Support<T> &supp2, 
            Support<T> &common);

template <typename T>
    T
    distance(const Support<T> &supp1, const Support<T> &supp2);

template <typename T>
    T
    distance(const Support<T> &supp1, const Support<T> &supp2, Support<T> &common);

template <typename T>
    T
    distance(const DenseVector<T> &singsupp1, const Support<T> &supp2);

template <typename T>
    T
    distance(const Support<T> &supp1, const DenseVector<T> &singsupp2);

template <typename T, typename S>
    Support<T>
    operator+(const Support<T> &supp, S shift);

template <typename S, typename T>
    Support<T>
    operator*(S factor, const Support<T> &supp);

template <typename T>
    std::ostream &
    operator<<(std::ostream &out, const Support<T> &supp);

} // namespace lawa

#endif // LAWA_CONSTRUCTIONS_SUPPORT_H
