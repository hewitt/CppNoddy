/// \file Sequential_Matrix_base.h
/// A base matrix class to ensure a consistent interface between
/// the inheriting dense/banded matrix classes.

#ifndef SEQUENTIALMATRIXBASE_H
#define SEQUENTIALMATRIXBASE_H

#include <DenseVector.h>

namespace CppNoddy {

  /// A base matrix class for sequential matrices.
  /// The sparse/banded/dense matrix classes inherit from
  /// this which allows for

  template < typename _Type >
  class Sequential_Matrix_base {

   public:

    /// An empty constructor
    Sequential_Matrix_base();

    virtual ~Sequential_Matrix_base();

    // ACCESS TO ELEMENTS -- virtual access operators.
    // These are slow if used via the base class
    virtual const _Type& operator()(const std::size_t& row, const std::size_t& col) const = 0;
    virtual _Type& operator()(const std::size_t& row, const std::size_t& col) = 0;
    virtual const _Type& get(const std::size_t& row, const std::size_t& col) const = 0;
    virtual _Type& set(const std::size_t& row, const std::size_t& col) = 0;

    // ENQUIRIES
    virtual std::size_t nrows() const = 0;
    virtual std::size_t ncols() const = 0;
    virtual std::size_t nelts() const = 0;

    // scaling
    virtual void scale( const _Type& mult ) = 0;

    // OUTPUT
    virtual void dump() const = 0;

   private:

  }
  ; // END CLASS

  template <typename _Type>
  Sequential_Matrix_base<_Type>::Sequential_Matrix_base()
  {}

  template <typename _Type>
  Sequential_Matrix_base<_Type>::~Sequential_Matrix_base()
  {}

} // end namespace


#endif
