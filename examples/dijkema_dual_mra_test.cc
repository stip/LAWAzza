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
    const int d  = 3;
    const int d_ = 3;

    MRA<double,Primal,Interval,Dijkema>  mra(d);
    MRA<double,Dual,Interval,Dijkema>    mra_(d, d_);
    mra_.enforceBoundaryConditions<DirichletBC, DirichletBC>();


    const int m = mra_.M0_.numRows();
    const int n = mra_.M0_.numCols();
    cerr << m << " x " << n << endl;
    lawa::GeMatrix<double> D(m,n);
    lawa::DenseVector<double> e(n);
    for (int i=1; i<=n; ++i) {
        e(i) = 1.;
        D(_,i) = mra_.M0_*e;
        e(i) = 0.;
    }
    cout << D << endl;


/*
    const int N = 2000;

    for (int i=0; i<N; ++i) {
        double x = double(i)/N;
        double y = mra_.phi_(x, mra_.min_j0, 1);
        cout << x << " " << y << endl;
    }
*/

    typedef BSpline<double,Primal,Interval,Dijkema>  RealBSpline;
    typedef BSpline<double,Dual,Interval,Dijkema>    RealBSpline_;

    SimpleIntegral<Trapezoidal,RealBSpline,RealBSpline_>  integral(mra.phi, mra_.phi_);

    const int j  = std::max(mra.min_j0, mra_.min_j0);

    const int k0 = mra.rangeI(j).firstIndex();
    const int k1 = mra.rangeI(j).lastIndex();

    const int k_0 = mra_.rangeI_(j).firstIndex();
    const int k_1 = mra_.rangeI_(j).lastIndex();

    cout << "j = " << j << endl;

    lawa::GeMatrix<double>  C(_(k0,k1), _(k_0,k_1));

    integral.quadrature.setN(6000);
    
    for (int k=k0; k<=k1; ++k) {
        for (int k_=k_0; k_<=k_1; ++k_) {
            C(k,k_) = integral(j, k, 0, j, k_, 0);
        }
    }

    cout << "C = " << C << endl;

    const double tol = 1e-4;
    for (int k=k0; k<=k1; ++k) {
        for (int k_=k_0; k_<=k_1; ++k_) {
            if (std::abs(C(k,k_))<tol) {
                C(k,k_) = 0;
            }
        }
    }

    cout << "C = " << C << endl;

    return 0;
}
