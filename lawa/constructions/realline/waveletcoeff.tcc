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

#ifndef LAWA_CONSTRUCTIONS_REALLINE_WAVELETCOORD_TCC
#define LAWA_CONSTRUCTIONS_REALLINE_WAVELETCOORD_TCC 1

#include <vector>

namespace lawa {

//== WaveletCoord ==============================================================

template <typename T, typename IndexType>
WaveletCoord<T,IndexType>::WaveletCoord(IndexType _j, IndexType _k, T _value)
    : j(_j), k(_k), value(_value)
{
}

//== WaveletCoordProxy =========================================================

template <typename T>
WaveletCoordProxy<T>::WaveletCoordProxy(T *_value)
    : value(_value)
{
}

template <typename T>
void
WaveletCoordProxy<T>::operator+=(const T &x)
{
    *value = x;
}

template <typename T>
void
WaveletCoordProxy<T>::operator-=(const T &x)
{
    *value = x;
}

//== WaveletCoordCmp ===========================================================

// return true if a < b
template <typename T, typename IndexType>
bool
WaveletCoordCmp::operator()(const WaveletCoord<T, IndexType> &a,
                            const WaveletCoord<T, IndexType> &b) const
{
    if (a.j<b.j) {
        return true;
    }
    if (a.j==b.j && a.k<b.k) {
        return true;
    }
    return false;

}

//== WaveletCoeffCoord =========================================================

template <typename T>
WaveletCoeffCoord<T>::WaveletCoeffCoord()
    : _lastSortedCoord(0), _isSorted(true), _isAccumulated(true)
{
    _coord.reserve(10);
}

template <typename T>
WaveletCoeffCoord<T>::~WaveletCoeffCoord()
{
}

//-- operators -----------------------------------------------------------------

template <typename T>
typename WaveletCoeffCoord<T>::ElementProxy
WaveletCoeffCoord<T>::operator()(IndexType j, IndexType k)
{
    if (_coord.size()>=_coord.capacity()) {
        accumulate();
        _coord.reserve(_coord.capacity()*2);
    }

    _coord.push_back(CoordType(j, k, T(0)));
    _isAccumulated = false;

    size_t lastIndex = _coord.size()-1;
    if ((lastIndex>0) && _isSorted) {
        if (_less(_coord[lastIndex-1], _coord[lastIndex])) {
            _lastSortedCoord = lastIndex;
        } else {
            _isSorted = false;
        }
    }

    return &_coord[lastIndex].value;
}

//-- methods -------------------------------------------------------------------

template <typename T>
void
WaveletCoeffCoord<T>::accumulate() const
{
//
//  Quick return if possible
//
    if (_isAccumulated) {
        ASSERT(_isSorted);
        return;
    }

//
//  sort
//
    if (!_isSorted) {
        std::sort(_coord.begin()+_lastSortedCoord+1, _coord.end(), _less);
    }
    if (_lastSortedCoord<_coord.size()) {
        std::inplace_merge(_coord.begin(),
                           _coord.begin() + _lastSortedCoord+1,
                           _coord.end(),
                           _less);
    }
    _isSorted = true;

//
//  accumulate values
//
    size_t k, K;
    for (k=0, K=1; K<_coord.size(); ++k, ++K) {
        while ((K<_coord.size()) && (!_less(_coord[k], _coord[K]))) {
            _coord[k].value += _coord[K].value;
            _coord[K].value = T(0);
            ++K;
        }
        if (K<_coord.size()) {
            _coord[k+1] = _coord[K];
        }
    }
    if ((k<_coord.size()) && (K-k-1>0)) {
#       ifndef NDEBUG
        size_t oldCapacity = _coord.capacity();
#       endif

        _coord.erase(_coord.end()-(K-k-1), _coord.end());

#       ifndef NDEBUG
        if (oldCapacity!=_coord.capacity()) {
            std::cerr << "[WARNING] Possible performance bottleneck in "
                      << "CoordStorage<T,Cmp,I>::accumulate()"
                      << std::endl;
        }
#       endif
    }
    _lastSortedCoord = _coord.size()-1;
    _isAccumulated = true;
}

template <typename T>
const typename WaveletCoeffCoord<T>::CoordVector &
WaveletCoeffCoord<T>::coordVector() const
{
    return _coord;
}

template <typename T>
const typename WaveletCoeffCoord<T>::IndexType
WaveletCoeffCoord<T>::minLevel() const
{
    accumulate();
    if (numNonZeros()>0) {
        return _coord[0].j;
    }
    return 0;
}

template <typename T>
const typename WaveletCoeffCoord<T>::IndexType
WaveletCoeffCoord<T>::maxLevel() const
{
    accumulate();
    if (numNonZeros()>0) {
        return _coord[numNonZeros()-1].j;
    }
    return 0;
}

template <typename T>
const typename WaveletCoeffCoord<T>::IndexType
WaveletCoeffCoord<T>::numNonZeros() const
{
    accumulate();
    return _coord.size();
}

//== WaveletCoeff ==============================================================

template <typename T>
WaveletCoeff<T>::WaveletCoeff()
{
}

template <typename T>
WaveletCoeff<T>::WaveletCoeff(const WaveletCoeffCoord<T> &coeffCoords)
{
    _compress(coeffCoords);
}

template <typename T>
WaveletCoeff<T>::~WaveletCoeff()
{
}

//-- operators -----------------------------------------------------------------

template <typename T>
void
WaveletCoeff<T>::operator=(const WaveletCoeffCoord<T> &coeffCoords)
{
    _compress(coeffCoords);
}

//-- methods -------------------------------------------------------------------

template <typename T>
const typename WaveletCoeff<T>::IndexType
WaveletCoeff<T>::numNonZeros() const
{
    return _values.length();
}

template <typename T>
const typename WaveletCoeff<T>::IndexType
WaveletCoeff<T>::minLevel() const
{
    return _FirstIndex.firstIndex();
}

template <typename T>
const typename WaveletCoeff<T>::IndexType
WaveletCoeff<T>::maxLevel() const
{
    return _FirstIndex.lastIndex()-1;
}

template <typename T>
const typename WaveletCoeff<T>::IndexType &
WaveletCoeff<T>::firstIndex(const IndexType j) const
{
    return _FirstIndex(j);
}

template <typename T>
const typename WaveletCoeff<T>::IndexType &
WaveletCoeff<T>::k(const IndexType i) const
{
    return _K(i);
}

template <typename T>
const typename WaveletCoeff<T>::ElementType &
WaveletCoeff<T>::value(const IndexType i) const
{
    return _values(i);
}

template <typename T>
void
WaveletCoeff<T>::_compress(const WaveletCoeffCoord<T> &coeffCoords)
{
    coeffCoords.accumulate();
    IndexType nnz = coeffCoords.numNonZeros();

    IndexType j0 = coeffCoords.minLevel()+1;
    IndexType J  = coeffCoords.maxLevel();

    _FirstIndex.resize(J-j0+3, j0-1);
    _K.resize(nnz);
    _values.resize(nnz);

    const auto &coord = coeffCoords.coordVector();

    IndexType j = j0-1;
    _FirstIndex(j) = 1;

    for (size_t p=0; p<coord.size(); ++p) {
        while (coord[p].j>j) {
            _FirstIndex(j+1) = p+1;
            ++j;
        }
        _K(p+1)      = coord[p].k;
        _values(p+1) = coord[p].value;
    }
    while (j<=J) {
        _FirstIndex(j+1) = coord.size() + 1;
        ++j;
    }
}

} // namespace lawa

#endif // LAWA_CONSTRUCTIONS_REALLINE_WAVELETCOORD_TCC
