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

#ifndef LAWA_CONSTRUCTIONS_INTERVAL_DIJKEMA_DUAL_SELECTBOUNDARYFUNCTIONS_H
#define LAWA_CONSTRUCTIONS_INTERVAL_DIJKEMA_DUAL_SELECTBOUNDARYFUNCTIONS_H 1

#include <cassert>
#include <lawa/enum.h>
#include <lawa/constructions/interval/dijkema/dual/mra.h>

namespace lawa {

template <typename T>
    void
    defaultSelectBoundaryFunctions(const MRA<T,Dual,Interval,Dijkema> &mra_,
                                   BoundarySide                       side,
                                   DenseVector<T>                     &indices);

} // namespace lawa

#endif // LAWA_CONSTRUCTIONS_INTERVAL_DIJKEMA_DUAL_SELECTBOUNDARYFUNCTIONS_H

