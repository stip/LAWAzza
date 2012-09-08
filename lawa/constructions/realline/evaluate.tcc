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

#ifndef LAWA_CONSTRUCTIONS_REALLINE_EVALUATE_TCC
#define LAWA_CONSTRUCTIONS_REALLINE_EVALUATE_TCC 1

#include <lawa/constructions/realline/evaluate.h>
#include <lawa/typedefs.h>

namespace lawa {

template <typename T, typename VCOEF>
T
evaluate(const MRA<T,Primal,R,CDF>          &mra,
         const T                            &x,
         const int                          j,
         const flens::DenseVector<VCOEF>    &coef,
         const int                          deriv)
{
    T result = T(0);
    for (int k=coef.firstIndex(); k<=coef.lastIndex(); ++k) {
        result += coef(k)*mra.phi(x, j, k, deriv);
    }
    return result;
}

template <typename T>
T
evaluate(const Basis<T,Primal,R,CDF>        &basis,
         const T                            &x,
         const WaveletCoeff<T>              &W,
         const int                          deriv)
{
    T result = T(0);

    const int j0 = W.minLevel()+1;
    const int J  = W.maxLevel();

    for (int i=W.firstIndex(j0-1); i<W.firstIndex(j0); ++i) {
        result += W.value(i)*basis.mra.phi(x, j0, W.k(i), deriv);
    }
    for (int j=j0; j<=J; ++j) {
        for (int i=W.firstIndex(j); i<W.firstIndex(j+1); ++i) {
            result += W.value(i)*basis.psi(x, j, W.k(i), deriv);
        }
    }
    return result;
}




template <typename T, typename VW, typename VK, typename VCOEF>
T
evaluate(const Basis<T,Primal,R,CDF>        &basis,
         const T                            &x,
         const int                          j0,
         const int                          J,
         const flens::DenseVector<VW>       &w,
         const flens::DenseVector<VK>       &k,
         const flens::DenseVector<VCOEF>    &coef,
         const int                          deriv)
{
    T result = T(0);
    for (int i=w(j0-1); i<w(j0); ++i) {
        result += coef(i)*basis.mra.phi(x, j0, k(i), deriv);
    }
    for (int j=j0; j<J; ++j) {
        for (int i=w(j); i<w(j+1); ++i) {
            result += coef(i)*basis.psi(x, j, k(i), deriv);
        }
    }
    return result;
}

} // namespace lawa

#endif // LAWA_CONSTRUCTIONS_REALLINE_EVALUATE_TCC
