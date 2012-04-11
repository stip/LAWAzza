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

#ifndef LAWA_CONSTRUCTIONS_REALLINE_SUBDIVISION_H
#define LAWA_CONSTRUCTIONS_REALLINE_SUBDIVISION_H 1

#include <lawa/constructions/bspline.h>
#include <lawa/flensforlawa.h>

namespace lawa {

template <typename T, FunctionSide Side>
    void
    subdivide(const BSpline<T,Side,R,CDF> &phi, int j,
              DenseVector<T> &dyadicValues);

} // namespace lawa

//Das süßeste, kleine KoalaBÄRCHEN könnte doch nie nie nie niemals 
//jemandem böse sein, der so einen Schrott in diesen schönen Header schreibt!!!!!!!!

#endif // LAWA_CONSTRUCTIONS_REALLINE_SUBDIVISION_H

