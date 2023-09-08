/// \file src/DenseVector.cpp
/// Implementation of the DenseVector class -- a dense, variable size, vector object.

#include <vector>
#include <complex>
#include <algorithm>
#include <functional>
#include <numeric>

#include <DenseVector.h>
#include <Exceptions.h>
#include <Functors.h>

namespace CppNoddy {

  template <typename _Type>
  DenseVector<_Type>::DenseVector()
  {}

  template <typename _Type>
  DenseVector<_Type>::DenseVector(const std::size_t& size, const _Type* p) {
    m_vec.reserve(size);
    // assign the array contents to the vector
    m_vec.assign(p, p + size);
  }

  template <typename _Type>
  DenseVector<_Type>::DenseVector(const std::size_t& size,
                                  const _Type& fill) : m_vec(size, fill)
  {}

  template <typename _Type>
  void DenseVector<_Type>::scale(const _Type& m) {
    operator*=(m);
  }

  template <typename _Type>
  void DenseVector<_Type>::add(const DenseVector<_Type>& x) {
    operator+=(x);
  }

  template <typename _Type>
  void DenseVector<_Type>::sub(const DenseVector<_Type>& x) {
    operator-=(x);
  }

  template <typename _Type>
  double DenseVector<_Type>::one_norm() const {
    return accumulate(begin(), end(), 0.0, absAdd_functor<_Type>());
  }

  template <typename _Type>
  double DenseVector<_Type>::two_norm() const {
    return std::sqrt(accumulate(begin(), end(), 0.0, absSquareAdd_functor<_Type>()));
  }

  template <typename _Type>
  double DenseVector<_Type>::inf_norm() const {
    return std::abs(*max_element(begin(), end(), abs_predicate<_Type>()));
  }

  template <typename _Type>
  void DenseVector<_Type>::dump() const {
    std::cout << "size = " << size() << "\n";
    for(std::size_t i = 0; i < size(); ++i) {
      std::cout << m_vec[ i ] << ", ";
    }
    std::cout << "\n";
  }


  template <>
  void DenseVector<double>::dump_file(std::string filename, int precision) const {
    std::ofstream dump;
    dump.open(filename.c_str());
    dump.precision(precision);
    dump.setf(std::ios::showpoint);
    dump.setf(std::ios::showpos);
    dump.setf(std::ios::scientific);
    for(std::size_t i = 0; i < m_vec.size(); ++i) {
      dump << i << " " << m_vec[ i ] << "\n";
    }
  }

  template <>
  void DenseVector<std::complex<double> >::dump_file(std::string filename, int precision) const {
    std::ofstream dump;
    dump.open(filename.c_str());
    dump.precision(precision);
    dump.setf(std::ios::showpoint);
    dump.setf(std::ios::showpos);
    dump.setf(std::ios::scientific);
    for(std::size_t i = 0; i < m_vec.size(); ++i) {
      dump << i << " " << m_vec[ i ].real() << " " << m_vec[ i ].imag()  << "\n";
    }
  }

  template <typename _Type>
  void DenseVector<_Type>::swap(const std::size_t& i,
                                const std::size_t& j) {
#ifdef PARANOID
    if((i >= size()) || (j >= size())) {
      std::string problem;
      problem = " The DenseVector.swap method is trying to access \n";
      problem += " outside the max/min number of elements. \n";
      if(i > size())
        throw ExceptionRange(problem, size(), i);
      if(j > size())
        throw ExceptionRange(problem, size(), j);
    }
#endif
    std::swap<_Type>(m_vec[ i ], m_vec[ j ]);
  }

  // the templated versions we require are:
  template class DenseVector<double>
  ;
  template class DenseVector<std::complex<double> >
  ;
  template class DenseVector<int>
  ;

} // end namespace
