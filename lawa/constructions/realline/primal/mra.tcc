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

namespace lawa {

template <typename T>
MRA<T,Primal,R,CDF>::MRA(int _d, int j)
    : d(_d), j0(j), phi(d), M0(phi)
{
}

template <typename T>
const int
MRA<T,Primal,R,CDF>::level() const
{
    return _j;
}

template <typename T>
void
MRA<T,Primal,R,CDF>::setLevel(int j) const
{
    _j = j;
}

} // namespace lawa

