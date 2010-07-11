/// \file Matrix_base.h
/// A base matrix class to ensure a consistent interface between
/// the inheriting dense/banded matrix classes.

#ifndef MATRIXBASE_H
#define MATRIXBASE_H

#include <DenseVector.h>

namespace CppNoddy
{

  /// A base matrix class. The sparse/banded/dense matrix classes inherit from
  /// this to ensure a consitent sub-interface across the three classes.

  template < typename _Type >
  class Matrix_base
  {

  public:

    /// An empty constructor
    Matrix_base();

    virtual ~Matrix_base();

    // ACCESS TO ELEMENTS -- virtual access operators.
    // These are slow if used via the base class
    virtual const _Type& operator()( const std::size_t& row, const std::size_t& col ) const = 0;
    virtual _Type& operator()( const std::size_t& row, const std::size_t& col ) = 0;
    virtual const _Type& get( const std::size_t& row, const std::size_t& col ) const = 0;
    virtual _Type& set( const std::size_t& row, const std::size_t& col ) = 0;

    // ENQUIRIES
    virtual std::size_t nrows() const = 0;
    virtual std::size_t ncols() const = 0;
    virtual std::size_t nelts() const = 0;

    // MATRIX OPERATIONS ON FULL MATIRX
    virtual void scale( const _Type& mult ) = 0;
    virtual void transpose() = 0;
    virtual double one_norm() const = 0;
    virtual double two_norm() const = 0;
    virtual double inf_norm() const = 0;
    virtual double frob_norm() const = 0;
    virtual DenseVector<_Type> multiply( const DenseVector<_Type>& X ) const = 0;

    // OUTPUT
    virtual void dump() const = 0;

  private:

  }
  ; // END CLASS

  template <typename _Type>
  Matrix_base<_Type>::Matrix_base()
  {}

  template <typename _Type>
  Matrix_base<_Type>::~Matrix_base()
  {}

} // end namespace


#endif
