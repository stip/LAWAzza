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

#ifndef LAWA_AUXILIRAY_ARROW_TCC
#define LAWA_AUXILIRAY_ARROW_TCC 1

#include <lawa/auxiliary/arrow.h>

namespace lawa {

template <typename T>
void
arrow(const GeMatrix<T> &In, GeMatrix<T> &Out)
{
    const int i0 = In.firstRow();
    const int i1 = In.lastRow();

    const int j0 = In.firstCol();
    const int j1 = In.lastCol();
    
    Out.resize(In);
    for (int j=j0; j<=j1; ++j) {
        for (int i=i0; i<=i1; ++i) {
            Out(i,j) = In(i1+i0-i, j1+j0-j);
        }
    }
}

} // namespace lawa

#endif // LAWA_AUXILIRAY_ARROW_TCC
