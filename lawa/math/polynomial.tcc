#include <iostream>

namespace lawa {

template <typename T>
Polynomial<T>::Polynomial(int n)
    : _coefficients(DenseVector<T>(_(0,n)))
{
}

template <typename T>
Polynomial<T>::Polynomial(const DenseVector<T> &coefficients)
    : _coefficients(coefficients)
{
}


template <typename T>
const T &
Polynomial<T>::operator()(int n) const
{
    assert(0<=n);
    assert(n<=this->degree());

    return _coefficients(n);
}

template <typename T>
T &
Polynomial<T>::operator()(int n)
{
    assert(0<=n);
    assert(n<=this->degree());

    return _coefficients(n);
}


template <typename T>
Polynomial<T> &
Polynomial<T>::operator+=(const Polynomial<T> &rhs)
{
    if (this->degree()>=rhs.degree()) {
        for (int i=0; i<=rhs.degree(); ++i) {
            _coefficients(i) += rhs._coefficients(i);
        }
    } else {
        DenseVector<T> tmp = rhs._coefficients;
        for (int i=0; i<=this->degree(); ++i) {
            tmp(i) += _coefficients(i);
        }
        this->_coefficients = tmp;
    }
    return *this;
}

template <typename T>
int
Polynomial<T>::degree() const
{
    return _coefficients.lastIndex();
}

//------------------------------------------------------------------------------

template <typename T>
Polynomial<T>
operator*(const Polynomial<T> &lhs, const Polynomial<T> &rhs)
{
    int degree = lhs.degree() + rhs.degree();

    Polynomial<T> res(degree);
    for (int i=0; i<=lhs.degree(); i++) {
        for (int j=0; j<=rhs.degree(); j++) {
            res(i+j) += lhs(i)*rhs(j);
        }
    }
    return res;
}

template <typename S, typename T>
Polynomial<T>
operator*(const S &lhs, const Polynomial<T> &rhs)
{
    Polynomial<T> res(rhs);
    for (int i=0; i<=rhs.degree(); ++i) {
        res(i) *= lhs;
    }
    return res;
}

template <typename T>
Polynomial<T>
pow(const Polynomial<T> &p, int n)
{
    if (!n) {
        Polynomial<T> res(0);
        res(0) = 1.;
        return res;
    }
    Polynomial<T> res(p);
    for (int k=1; k<n; ++k) {
        res = res*p;
    }
    return res;
}

template <typename T>
std::ostream &
operator<<(std::ostream &out, const Polynomial<T> &rhs)
{
    out << "(";
    for (int i=0; i<=rhs.degree(); ++i) {
        out << rhs(i) << " ";
    }
    out << ")";
    return out;
}

} // namespace lawa

