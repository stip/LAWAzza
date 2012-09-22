/*
  This file is part of LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008-2011  Mario Rometsch, Alexander Stippler.

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

#ifndef LAWA_CONSTRUCTIONS_INTERVAL_DIJKEMA_PRIMAL_SPLINEHELPER_TCC
#define LAWA_CONSTRUCTIONS_INTERVAL_DIJKEMA_PRIMAL_SPLINEHELPER_TCC 1

#include <list>

namespace lawa {
    
template <typename T>
T
w(int i, int d, const DenseVector<T> &knots, T x)
{
    assert(1<=i);
    assert(i<=knots.length()-d+1);
    
    if (x<=knots(i)) {
        return 0.0;
    } else if (x>=knots(i+d-1)) {
        return 1.0;
    } else {
        return (x-knots(i)) / (knots(i+d-1)-knots(i));
    }
}

template <typename T>
GeMatrix<T>
insertKnot(int d, DenseVector<T> &knots, T x)
{
    assert(knots.length()-d-1>=1);

    GeMatrix<T>  ret(knots.length()-d, knots.length()-d-1);

    for (int i=ret.firstCol(); i<=ret.lastCol(); ++i) {
        ret(i,i) = w(i,d+1,knots,x);
        ret(i+1,i) = 1-w(i+1,d+1,knots,x);
    }

    std::list<T> temp;
    for (int i=knots.firstIndex(); i<=knots.lastIndex(); ++i) {
        temp.push_back(knots(i));
    }
    temp.push_back(x);
    temp.sort();

    knots.resize(knots.length()+1);
    typename std::list<T>::const_iterator it=temp.begin();
    for (int i=1; it!=temp.end(); ++it, ++i) {
        knots(i) = *it;
    }
    return ret;
}

} // namespace lawa

#endif // LAWA_CONSTRUCTIONS_INTERVAL_DIJKEMA_PRIMAL_SPLINEHELPER_TCC
