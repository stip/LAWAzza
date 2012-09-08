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

#ifndef LAWA_CONSTRUCTIONS_REALLINE_WAVELETCOORD_H
#define LAWA_CONSTRUCTIONS_REALLINE_WAVELETCOORD_H 1

#include <vector>

namespace lawa {

//== WaveletCoord ==============================================================

template <typename T, typename IndexType>
struct WaveletCoord
{
    WaveletCoord(IndexType j, IndexType k, T value);
    
    IndexType   j;
    IndexType   k;
    T           value;
};

//== WaveletCoordProxy =========================================================

template <typename T>
struct WaveletCoordProxy
{
    WaveletCoordProxy(T *value);

    void
    operator+=(const T &x);

    void
    operator-=(const T &x);

    private:
        T *value;

    //
    //  We only allow accumulation of data.  Assignment would be too expensive.
    //
        void
        operator=(const T &x);
};

//== WaveletCoordCmp ===========================================================

struct WaveletCoordCmp
{
        // return true if a < b
        template <typename T, typename IndexType>
            bool
            operator()(const WaveletCoord<T, IndexType> &a,
                       const WaveletCoord<T, IndexType> &b) const;
};

//== WaveletCoeffCoord =========================================================

template <typename T>
class WaveletCoeffCoord
{
    public:
        typedef T                           ElementType;
        typedef int                         IndexType;

        typedef WaveletCoordProxy<T>        ElementProxy;

        typedef WaveletCoord<T, IndexType>  CoordType;
        typedef std::vector<CoordType>      CoordVector;

        WaveletCoeffCoord();

        ~WaveletCoeffCoord();

        //-- operators ---------------------------------------------------------

        ElementProxy
        operator()(IndexType j, IndexType k);

        //-- methods -----------------------------------------------------------

        void
        accumulate() const;

        const CoordVector &
        coordVector() const;

        const IndexType
        minLevel() const;

        const IndexType
        maxLevel() const;

        const int
        numNonZeros() const;

        WaveletCoeffCoord(const WaveletCoeffCoord &rhs);

    private:

        mutable CoordVector  _coord;
        mutable size_t       _lastSortedCoord;
        mutable bool         _isSorted;
        mutable bool         _isAccumulated;
        WaveletCoordCmp      _less;
};

//== WaveletCoeff ==============================================================

template <typename T>
class WaveletCoeff
{
    public:
        typedef T                           ElementType;
        typedef int                         IndexType;

        typedef DenseVector<ElementType>    ElementTypeVector;
        typedef DenseVector<IndexType>      IndexTypeVector;

        WaveletCoeff();

        WaveletCoeff(const WaveletCoeffCoord<T> &coeffCoords);

        ~WaveletCoeff();

        //-- operators ---------------------------------------------------------

        void
        operator=(const WaveletCoeffCoord<T> &coeffCoords);

        //-- methods -----------------------------------------------------------

        const IndexType
        numNonZeros() const;

        const IndexType
        minLevel() const;
        
        const IndexType
        maxLevel() const;

        const IndexType &
        firstIndex(const IndexType j) const;

        const IndexType &
        k(const IndexType i) const;

        const ElementType &
        value(const IndexType i) const;

        void
        _compress(const WaveletCoeffCoord<T> &coeffCoords);

        // not implemented on purpose!
        WaveletCoeff(const WaveletCoeff &rhs);


        DenseVector<IndexType>  _FirstIndex;
        DenseVector<IndexType>  _K;
        DenseVector<T>          _values;
};


} // namespace lawa

#endif // LAWA_CONSTRUCTIONS_REALLINE_WAVELETCOORD_H
