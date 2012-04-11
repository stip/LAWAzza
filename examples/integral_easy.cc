#include <iostream>

#include <lawa/lawa.h>
#include <lawa/lawa.tcc>

using namespace lawa;
using namespace std;
namespace lawa {

template <QuadratureType Quad, typename First, typename Second>
struct Integral
{
    typedef typename First::T T;

    Integral(const First &first, const Second &second);

    T
    operator()(int _j1, Integer _k1, int _deriv1,
               int _j2, Integer _k2, int _deriv2) const;

    T
    integrand(T x) const;

    const First &first;
    const Second &second;
    mutable Quadrature<Quad,Integral<Quad,First,Second> > quadrature;
    mutable int j1, deriv1,
                j2, deriv2;
    mutable Integer k1, k2;
};

template <QuadratureType Quad, typename First, typename Second>
Integral<Quad,First,Second>::Integral(const First &_first, 
                                      const Second &_second)
    : first(_first), second(_second), quadrature(*this)
{
}

template <QuadratureType Quad, typename First, typename Second>
typename First::T
Integral<Quad,First,Second>::operator()(int _j1, Integer _k1, int _deriv1, 
                                        int _j2, Integer _k2, int _deriv2) const
{
    j1 = _j1; k1 = _k1; deriv1 = _deriv1;
    j2 = _j2; k2 = _k2; deriv2 = _deriv2;
    return _integrate(*this);
}

template <QuadratureType Quad, typename First, typename Second>
typename First::T
Integral<Quad,First,Second>::integrand(typename First::T x) const
{
    return _integrand(*this, x);
}

template <typename First, typename Second>
typename First::T
_integrate(const Integral<Gauss,First,Second> &integral)
{
    typedef typename First::T T;

    const First &first = integral.first;
    const Second &second = integral.second;
    // the (minimal) width of the polynomial pieces.
    T unit = std::min(first.tic(integral.j1), second.tic(integral.j2));
    integral.quadrature.setOrder((integral.first.d - integral.deriv1 + 
                                  integral.second.d - integral.deriv2)/2 + 1);

    T ret = 0.;
    Support<T> common;
    if (overlap(first.support(integral.j1,integral.k1),
               second.support(integral.j2,integral.k2),common)) {
        T a = common.l1;
        for (T b=a+unit; b<=common.l2; b+=unit) {
            ret += integral.quadrature(a,b);
            a = b;
        }
    }

    return ret;
}

template <QuadratureType Quad, typename First, typename Second>
typename First::T
_integrand(const Integral<Quad,First,Second> &integral, typename First::T x)
{
    const First &first = integral.first;
    const Second &second = integral.second;
    return first(x,integral.j1,integral.k1,integral.deriv1)
         * second(x,integral.j2,integral.k2,integral.deriv2);
}

} // namespace lawa

int
main()
{
    const unsigned int d  = 3;
    typedef BSpline<double,Primal,R,CDF> RealBSpline;

    RealBSpline phi(d);
    Integral<Gauss,RealBSpline,RealBSpline> integral(phi,phi);

    int first = -phi.l2+1,
        last  = -phi.l1;
    for (int k=first; k<=last; ++k) { 
        cout << integral(0,k,0,
                         0,0,0) << endl;
    }

    return 0;
}
