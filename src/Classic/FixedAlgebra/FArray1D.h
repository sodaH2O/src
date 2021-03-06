#ifndef CLASSIC_FArray1D_HH
#define CLASSIC_FArray1D_HH

// ------------------------------------------------------------------------
// $RCSfile: FArray1D.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1.2.3 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Template class: FArray1D<T,N>
//
// ------------------------------------------------------------------------
// Class category: FixedAlgebra
// ------------------------------------------------------------------------
//
// $Date: 2004/11/18 22:18:06 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "Utilities/CLRangeError.h"
#include "Utilities/SizeError.h"
#include <algorithm>
#include <initializer_list>
#include <ostream>


// Template class FArray1D<T,N>
// ------------------------------------------------------------------------
/// A templated representation for one-dimensional arrays.
//  This version has fixed dimensions.
//  It implements array access, but contains no arithmetic operations.
//  The destructor generated by the compiler performs the correct operation.
//  For speed reasons it is not implemented.

template<class T, int N>
class FArray1D {

public:

    /// The value type of this array.
    typedef T value_type;

    /// Iterator for the array.
    typedef T *iterator;

    /// Iterator for constant array.
    typedef const T *const_iterator;

    /// Default constructor.
    //  Constructs a zero array.
    FArray1D();

    /// Constructor.
    //  Set all array elements to [b]t[/b].
    explicit FArray1D(const T &t);

    /// Copy constructor.
    FArray1D(const FArray1D &);

    /// Consructor with initializer list (needs C++11) (see http://en.cppreference.com/w/cpp/utility/initializer_list)
    FArray1D(const std::initializer_list<T>&);
    
    /// Assignment.
    const FArray1D &operator=(const FArray1D &);

    /// Get element.
    //  Return a reference to element [b]n[/b].
    //  Throw RangeError, if [b]n[/b] is out of range.
    T &operator()(int n);

    /// Get element.
    //  Return a constant reference to element [b]n[/b].
    //  Throw RangeError, if [b]n[/b] is out of range.
    const T &operator()(int n) const;

    /// Get element.
    //  Return a reference to element [b]n[/b].
    //  Result is undefined, if [b]n[/b] is out of range.
    T &operator[](int n);

    /// Get element.
    //  Return a reference to element [b]n[/b].
    //  Result is undefined, if [b]n[/b] is out of range.
    const T &operator[](int n) const;

    /// Get iterator pointing to beginning of array.
    //  Version for non-constant array.
    iterator begin();

    /// Get iterator pointing past end of array.
    //  Version for non-constant array.
    iterator end();

    /// Get iterator pointing to beginning of array.
    //  Version for constant array.
    const_iterator begin() const;

    /// Get iterator pointing past end of array.
    //  Version for constant array.
    const_iterator end() const;

    /// Get array size.
    int size() const;

protected:

    // The array data.
    T data[N];
};


// Template function implementation.
// ------------------------------------------------------------------------

template<class T, int N>
FArray1D<T, N>::FArray1D() {
    std::fill(begin(), end(), T(0));
}


template<class T, int N>
FArray1D<T, N>::FArray1D(const FArray1D &rhs) {
    std::copy(rhs.begin(), rhs.end(), begin());
}


template<class T, int N>
FArray1D<T, N>::FArray1D(const T &val) {
    std::fill(begin(), end(), val);
}

#include <iostream>

template<class T, int N>
FArray1D<T, N>::FArray1D(const std::initializer_list<T>& init) {
    
    
    if (init.size() > N) {
        throw SizeError("FArray1D<T,N>(const std::initializer_list<T>&)","Size exceeds array dimension.");
    }
    
    // copy list to array
    std::copy(init.begin(), init.end(), begin());
    
    // if initializer list smaller than array dimension, fill rest with zero
    if (init.size() < N) {
        std::fill(begin()+init.size(), end(), T(0));
    }
}

template<class T, int N>
const FArray1D<T, N> &FArray1D<T, N>::operator=(const FArray1D &rhs) {
    std::copy(rhs.begin(), rhs.end(), begin());
    return *this;
}


template<class T, int N> inline
T &FArray1D<T, N>::operator[](int i) {
    return data[i];
}


template<class T, int N> inline
const T &FArray1D<T, N>::operator[](int i) const {
    return data[i];
}


template<class T, int N>
T &FArray1D<T, N>::operator()(int i) {
    if(i >= N) {
        throw CLRangeError("FArray1D::operator()", "Index out of range.");
    }
    return data[i];
}


template<class T, int N>
const T &FArray1D<T, N>::operator()(int i) const {
    if(i >= N) {
        throw CLRangeError("FArray1D::operator()", "Index out of range.");
    }
    return data[i];
}


template<class T, int N> inline
typename FArray1D<T, N>::iterator FArray1D<T, N>::begin() {
    return data;
}


template<class T, int N> inline
typename FArray1D<T, N>::iterator FArray1D<T, N>::end() {
    return data + N;
}


template<class T, int N> inline
typename FArray1D<T, N>::const_iterator FArray1D<T, N>::begin() const {
    return data;
}


template<class T, int N> inline
typename FArray1D<T, N>::const_iterator FArray1D<T, N>::end() const {
    return data + N;
}


template<class T, int N> inline
int FArray1D<T, N>::size() const {
    return N;
}


template <class T, int N>
std::ostream &operator<<(std::ostream &os, const FArray1D<T, N> &v) {
    for(int i = 0; i < N; ++i) {
        os << v[i] << " ";
    }

    os << std::endl;
    return os;
}

#endif // CLASSIC_FArray1D_HH
