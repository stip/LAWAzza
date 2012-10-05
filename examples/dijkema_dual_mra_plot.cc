#include <cmath>
#include <fstream>
#include <iostream>

#include <lawa/lawa.h>
#include <lawa/lawa.tcc>

namespace lawa {

template <QuadratureType Quad, typename First, typename Second>
struct SimpleIntegral
{
    typedef typename First::T T;

    SimpleIntegral(const First &first, const Second &second);

    T
    operator()(int _j1, Integer _k1, int _deriv1,
               int _j2, Integer _k2, int _deriv2) const;

    T
    integrand(T x) const;

    const First &first;
    const Second &second;
    mutable Quadrature<Quad, SimpleIntegral<Quad, First, Second> > quadrature;
    mutable int j1, deriv1,
                j2, deriv2;
    mutable Integer k1, k2;
};

template <QuadratureType Quad, typename First, typename Second>
SimpleIntegral<Quad,First,Second>::SimpleIntegral(const First &_first, 
                                                  const Second &_second)
    : first(_first), second(_second), quadrature(*this)
{
}

template <QuadratureType Quad, typename First, typename Second>
typename First::T
SimpleIntegral<Quad,First,Second>::operator()(int _j1, Integer _k1, int _deriv1, 
                                              int _j2, Integer _k2, int _deriv2) const
{
    j1 = _j1; k1 = _k1; deriv1 = _deriv1;
    j2 = _j2; k2 = _k2; deriv2 = _deriv2;
    return _integrate(*this);
}

template <QuadratureType Quad, typename First, typename Second>
typename First::T
SimpleIntegral<Quad,First,Second>::integrand(typename First::T x) const
{
    return _integrand(*this, x);
}

template <typename First, typename Second>
typename First::T
_integrate(const SimpleIntegral<Trapezoidal,First,Second> &integral)
{
    typedef typename First::T T;
    
    const First &first = integral.first;
    const Second &second = integral.second;

    Support<T> common;
    if (overlap(first.support(integral.j1,integral.k1),
                second.support(integral.j2,integral.k2),
                common)) 
    {
        return integral.quadrature(common.l1, common.l2);
    } else {
        return 0;
    }
}


template <QuadratureType Quad, typename First, typename Second>
typename First::T
_integrand(const SimpleIntegral<Quad,First,Second> &integral, typename First::T x)
{
    const First &first = integral.first;
    const Second &second = integral.second;
    return first(x,integral.j1, integral.k1, integral.deriv1)
         * second(x,integral.j2, integral.k2);
}


} // namespace lawa




using namespace lawa;
using namespace std;

int
main()
{
    const int d  = 2;
    const int d_ = 4;

    MRA<double,Primal,Interval,Dijkema>  mra(d);
    MRA<double,Dual,Interval,Dijkema>    mra_(d, d_);

#   ifdef ENFORCE_BC_LEFT
    std::cerr << "compiled with ENFORCE_BC_LEFT" << std::endl;
    std::cerr << "enforceBoundaryConditions:" << std::endl;
    mra.enforceBoundaryConditions<DirichletBC, NoBC>();
    mra_.enforceBoundaryConditions<DirichletBC, NoBC>();
#   elif defined ENFORCE_BC_RIGHT
    std::cerr << "compiled with ENFORCE_BC_RIGHT" << std::endl;
    std::cerr << "enforceBoundaryConditions:" << std::endl;
    mra.enforceBoundaryConditions<NoBC, DirichletBC>();
    mra_.enforceBoundaryConditions<NoBC, DirichletBC>();
#   else
    std::cerr << "compiled without ENFORCE_BC_LEFT or ENFORCE_BC_RIGHT"
              << std::endl;
#   endif // ENFORCE_BC

    typedef BSpline<double,Primal,Interval,Dijkema>  RealBSpline;
    typedef BSpline<double,Dual,Interval,Dijkema>    RealBSpline_;

    const int j0 = std::max(mra.j0,  mra_.j0);


    const int k0 = mra.rangeI(j0).firstIndex();
    const int k1 = mra.rangeI(j0).lastIndex();

    RealBSpline bSpline(mra);

    cerr << "j0 = " << j0 << endl;
    cerr << "primal range = " << 0 << ":" << k1-k0 << endl;

    const int N = 1000;
    for (int k=k0; k<=k1; ++k) {
        for (int i=0; i<=N; ++i) {
            const double x = i/double(N);
            cout << x << " " << bSpline(x, j0, k, 0) << endl;
        }
        cout << endl << endl;
    }


    const int k0_ = mra_.rangeI_(j0).firstIndex();
    const int k1_ = mra_.rangeI_(j0).lastIndex();

    RealBSpline_ bSpline_(mra_);

    cerr << "dual range = " << k1-k0+1 << ":" << k1_-k0_ + (k1-k0+1) << endl;

    for (int k=k0_; k<=k1_; ++k) {
        for (int i=0; i<=N; ++i) {
            const double x = i/double(N);
            cout << x << " " << bSpline_(x, j0, k) << endl;
        }
        cout << endl << endl;
    }



#   ifdef ENFORCE_BC_LEFT
    std::cerr << "compiled with ENFORCE_BC_LEFT" << std::endl;
#   elif defined ENFORCE_BC_RIGHT
    std::cerr << "compiled with ENFORCE_BC_RIGHT" << std::endl;
#   else
    std::cerr << "compiled without ENFORCE_BC_LEFT or ENFORCE_BC_RIGHT"
              << std::endl;
#   endif // ENFORCE_BC


    return 0;
}
