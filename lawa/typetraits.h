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

#ifndef LAWA_TYPETRAITS_H
#define LAWA_TYPETRAITS_H 1

#include <lawa/flensforlawa.h>
#include <lawa/constructions/bspline.h>
#include <lawa/constructions/wavelet.h>

namespace lawa {

//--- IsPrimal

template <typename X>
struct IsPrimal
{
    static const bool value = false;
};

template <typename T, DomainType Domain, Construction Cons,
          template <typename, FunctionSide, DomainType, Construction> class Some>
struct IsPrimal<Some<T, Primal, Domain, Cons> >
{
    static const bool value = true;
};


//--- IsOrthogonal

template <typename X>
struct IsOrthogonal
{
    static const bool value = false;
};

template <typename T, DomainType Domain, Construction Cons,
    template <typename, FunctionSide, DomainType, Construction> class Some>
struct IsOrthogonal<Some<T, Orthogonal, Domain, Cons> >
{
    static const bool value = true;
};

//--- IsDual

template <typename X>
struct IsDual
{
    static const bool value = false;
};

template <typename T, DomainType Domain, Construction Cons,
          template <typename, FunctionSide, DomainType, Construction> class Some>
struct IsDual<Some<T, Dual, Domain, Cons> >
{
    static const bool value = true;
};

//--- PrimalOrDual
template <typename X>
struct PrimalOrDual
{
    static const bool value = (IsPrimal<X>::value || IsDual<X>::value);
};

//--- PrimalOrDualOrOrthogonal
template <typename X>
struct PrimalOrDualOrOrthogonal
{
    static const bool value = (IsPrimal<X>::value || IsDual<X>::value || IsOrthogonal<X>::value);
};

//--- BothPrimal
template <typename X, typename Y>
struct BothPrimal
{
    static const bool value = IsPrimal<X>::value && IsPrimal<Y>::value;
};


//--- BothDual
template <typename X, typename Y>
struct BothDual
{
    static const bool value = IsDual<X>::value && IsDual<Y>::value;
};

//--- BothOrthogonal
template <typename X, typename Y>
struct BothOrthogonal
{
    static const bool value = IsOrthogonal<X>::value && IsOrthogonal<Y>::value;
};

//--- PrimalOrDual
template <typename X, typename Y>
struct PrimalAndDual
{
    static const bool value = (IsPrimal<X>::value && IsDual<Y>::value);
};

//--- IsWavelet
template <typename X>
struct IsWavelet
{
    static const bool value = false;
};

template <typename T, FunctionSide Side, DomainType Domain, Construction Cons>
struct IsWavelet<Wavelet<T,Side,Domain,Cons> >
{
    static const bool value = true;
};

//--- IsBSpline
template <typename X>
struct IsBSpline
{
    static const bool value = false;
};

template <typename T, FunctionSide Side, DomainType Domain, Construction Cons>
struct IsBSpline<BSpline<T,Side,Domain,Cons> >
{
    static const bool value = true;
};

//--- IsPeriodic
template <typename X>
struct IsPeriodic
{
    static const bool value = false;
};

template <typename T, FunctionSide Side, Construction Cons>
struct IsPeriodic<Wavelet<T,Side,Periodic,Cons> >
{
    static const bool value = true;
};

template <typename T, FunctionSide Side, Construction Cons>
struct IsPeriodic<BSpline<T,Side,Periodic,Cons> >
{
    static const bool value = true;
};

} // namespace lawa

#endif // LAWA_TYPETRAITS_H
