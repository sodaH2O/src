#ifndef CLASSIC_FMatrix_HH
#define CLASSIC_FMatrix_HH

// ------------------------------------------------------------------------
// $RCSfile: FMatrix.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.2.2.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Template class: FMatrix<T,R,C>
//
// ------------------------------------------------------------------------
// Class category: FixedAlgebra
// ------------------------------------------------------------------------
//
// $Date: 2004/11/18 22:18:06 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "FixedAlgebra/FArray2D.h"
#include "FixedAlgebra/FVector.h"
#include "Utilities/SizeError.h"
#include <algorithm>
#include <numeric>


// Template class FMatrix<T,R,C>
// ------------------------------------------------------------------------
/// A templated representation for matrices.
//  This class implements the basic algebraic operations on general matrices,
//  which need not be square.
//  The destructor generated by the compiler does the correct thing,
//  and is not defined for speed.

template<class T, int R, int C>
class FMatrix: public FArray2D<T, R, C> {

public:

    /// Default constructor.
    //  Construct zero matrix.
    FMatrix();

    /// Copy constructor.
    FMatrix(const FMatrix &);

    /// Conversion from two-dimensional array.
    FMatrix(const FArray2D<T, R, C> &);

    /// Constructor.
    //  Set all matrix elements to [b]t[/b].
    explicit FMatrix(const T &t);

    /// Assignment.
    FMatrix &operator=(const FMatrix &);

    /// Constructor.
    //  Fill all matrix elements from the C-array [b]t[/b].
    explicit FMatrix(const T *t);

    /// Convert and assign.
    FMatrix &operator=(const FArray2D<T, R, C> &);

    /// Multiply by scalar and assign.
    FMatrix &operator*=(const T &);

    /// Divide by scalar and assign.
    FMatrix &operator/=(const T &);

    /// Add matrix and assign.
    FMatrix &operator+=(const FMatrix &);

    /// Subtract ,atrix and assign.
    FMatrix &operator-=(const FMatrix &);

    /// FMatrix transpose.
    FMatrix<T, C, R> transpose() const;
};


// ------------------------------------------------------------------------

/// Unary minus.
template <class T, int R, int C>
FMatrix<T, R, C> operator-(const FMatrix<T, R, C> &);

/// Add.
template <class T, int R, int C>
FMatrix<T, R, C> operator+(const FMatrix<T, R, C> &, const FMatrix<T, R, C> &);

/// Subtract.
template <class T, int R, int C>
FMatrix<T, R, C> operator-(const FMatrix<T, R, C> &, const FMatrix<T, R, C> &);

/// Multiply by column vector.
template <class T, int R, int C>
FVector<T, R> operator*(const FMatrix<T, R, C> &, const FVector<T, C> &);

/// Row vector times Matrix.
template <class T, int R, int C>
FVector<T, C> operator*(const FVector<T, R> &, const FMatrix<T, R, C> &);

/// Matrix times scalar.
template <class T, int R, int C>
FMatrix<T, R, C> operator*(const FMatrix<T, R, C> &, const T &);

/// Scalar times Matrix.
template <class T, int R, int C>
FMatrix<T, R, C> operator*(const T &, const FMatrix<T, R, C> &);

/// Matrix divided by scalar.
template <class T, int R, int C>
FMatrix<T, R, C> operator/(const FMatrix<T, R, C> &, const T &);

/// Add unit matrix times [b]t[/b] to matrix.
template <class T, int C>
FMatrix<T, C, C> operator+(const FMatrix<T, C, C> &, const T &t);

/// Subtract unit matrix times [b]t[/b] from matrix.
template <class T, int C>
FMatrix<T, C, C> operator-(const FMatrix<T, C, C> &, const T &t);


// Public methods of FMatrix<T,R,C>.
// ------------------------------------------------------------------------

template<class T, int R, int C>
FMatrix<T, R, C>::FMatrix():
    FArray2D<T, R, C>()
{}


template<class T, int R, int C>
FMatrix<T, R, C>::FMatrix(const FMatrix<T, R, C> &rhs):
    FArray2D<T, R, C>(rhs)
{}


template<class T, int R, int C>
FMatrix<T, R, C>::FMatrix(const FArray2D<T, R, C> &rhs):
    FArray2D<T, R, C>(rhs)
{}


template<class T, int R, int C>
FMatrix<T, R, C>::FMatrix(const T &rhs):
    FArray2D<T, R, C>(rhs)
{}


template<class T, int R, int C>
FMatrix<T, R, C>::FMatrix(const T *rhs):
    FArray2D<T, R, C>(rhs)
{}


template<class T, int R, int C>
FMatrix<T, R, C> &FMatrix<T, R, C>::operator=(const FMatrix<T, R, C> &rhs) {
    FArray2D<T, R, C>::operator=(rhs);
    return *this;
}


template<class T, int R, int C>
FMatrix<T, R, C> &FMatrix<T, R, C>::operator=(const FArray2D<T, R, C> &rhs) {
    FArray2D<T, R, C>::operator=(rhs);
    return *this;
}


template<class T, int R, int C>
FMatrix<T, R, C> &FMatrix<T, R, C>::operator*=(const T &rhs) {
    std::transform(this->begin(), this->end(), this->begin(),
                   std::bind2nd(std::multiplies<T>(), rhs));
    return *this;
}


template<class T, int R, int C>
FMatrix<T, R, C> &FMatrix<T, R, C>::operator/=(const T &rhs) {
    std::transform(this->begin(), this->end(), this->begin(),
                   std::bind2nd(std::divides<T>(), rhs));
    return *this;
}


template<class T, int R, int C>
FMatrix<T, R, C> &FMatrix<T, R, C>::operator+=(const FMatrix<T, R, C> &array) {
    std::transform(this->begin(), this->end(), array.begin(), this->begin(),
                   std::plus<T>());
    return *this;
}


template<class T, int R, int C>
FMatrix<T, R, C> &FMatrix<T, R, C>::operator-=(const FMatrix<T, R, C> &array) {
    std::transform(this->begin(), this->end(), array.begin(), this->begin(),
                   std::minus<T>());
    return *this;
}


template <class T, int R, int C>
FMatrix<T, C, R> FMatrix<T, R, C>::transpose() const {
    FMatrix<T, C, R> result;

    for(int i = 0; i < R; ++i) {
        std::copy(this->row_begin(i), this->row_end(i), result.col_begin(i));
    }

    return result;
}


// Global functions acting on FMatrix<T,R,C>.
// ------------------------------------------------------------------------

template<class T, int R, int C>
FMatrix<T, R, C> operator-(const FMatrix<T, R, C> &rhs) {
    FMatrix<T, R, C> result;
    std::transform(rhs.begin(), rhs.end(), result.begin(), std::negate<T>());
    return result;
}


template<class T, int R, int C>
FMatrix<T, R, C> operator+(const FMatrix<T, R, C> &lhs, const FMatrix<T, R, C> &rhs) {
    FMatrix<T, R, C> result(lhs);
    return result += rhs;
}


template<class T, int R, int C>
FMatrix<T, R, C> operator-(const FMatrix<T, R, C> &lhs, const FMatrix<T, R, C> &rhs) {
    FMatrix<T, R, C> result(lhs);
    return result -= rhs;
}


template <class T, int R, int I, int C>
FMatrix<T, R, C> operator*(const FMatrix<T, R, I> &lhs, const FMatrix<T, I, C> &rhs) {
    FMatrix<T, R, C> result;

    for(int i = 0; i < R; i++) {
        typename FMatrix<T, R, I>::const_row_iterator i1 = lhs.row_begin(i);
        typename FMatrix<T, R, I>::const_row_iterator i2 = lhs.row_end(i);

        for(int j = 0; j < C; j++) {
            result(i, j) = std::inner_product(i1, i2, rhs.col_begin(j), T(0));
        }
    }

    return result;
}


template<class T, int R, int C>
FVector<T, R> operator*(const FMatrix<T, R, C> &lhs, const FVector<T, C> &rhs) {
    FVector<T, R> result;

    for(int i = 0; i < R; ++i) {
        result[i] = std::inner_product(lhs.row_begin(i), lhs.row_end(i),
                                       rhs.begin(), T(0));
    }

    return result;
}


template<class T, int R, int C>
FVector<T, C> operator*(const FVector<T, R> &lhs, const FMatrix<T, R, C> &rhs) {
    FVector<T, C> result;

    for(int j = 0; j < C; j++) {
        result[j] = std::inner_product(lhs.begin(), lhs.end(),
                                       rhs.col_begin(j), T(0));
    }

    return result;
}


template<class T, int R, int C>
FMatrix<T, R, C> operator*(const FMatrix<T, R, C> &lhs, const T &rhs) {
    FMatrix<T, R, C> result(lhs);
    return result *= rhs;
}


template<class T, int R, int C>
FMatrix<T, R, C> operator*(const T &lhs, const FMatrix<T, R, C> &rhs) {
    FMatrix<T, R, C> result(rhs);
    return result *= lhs;
}


template<class T, int R, int C>
FMatrix<T, R, C> operator/(const FMatrix<T, R, C> &lhs, const T &rhs) {
    FMatrix<T, R, C> result(lhs);
    return result /= rhs;
}


// Special functions on square matrices.
// ------------------------------------------------------------------------

template<class T, int C>
FMatrix<T, C, C> operator+(const FMatrix<T, C, C> &lhs, const T &rhs) {
    FMatrix<T, C, C> result(lhs);
    for(int i = 0; i < C; i++) result[i][i] += rhs;
    return result;
}


template<class T, int C>
FMatrix<T, C, C> operator-(const FMatrix<T, C, C> &lhs, const T &rhs) {
    FMatrix<T, C, C> result(lhs);
    for(int i = 0; i < C; i++) result[i][i] -= rhs;
    return result;
}

#endif // CLASSIC_FMatrix_HH