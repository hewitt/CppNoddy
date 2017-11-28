/// \file DenseVector.h
/// Specification for a templated DenseVector class -- a dense, dynamic, vector object.

#ifndef DENSEVECTOR_H
#define DENSEVECTOR_H

#include <vector>
#include <complex>
#include <algorithm>
#include <functional>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>

#include <Exceptions.h>
#include <Functors.h>

namespace CppNoddy
{

  // forward declare the matrix classes that will be friends
  // conversion from banded to vector (for LAPACK routines)
  // is best done if we allow BandedMatrix direct access to the
  // storage container vec.
  template <typename _Type>
  class BandedMatrix;

  /// An DenseVector class -- a dense vector object.
  /// This is templated but intended ONLY for double or
  /// std::complex<double>. We just encapsulate the STL vector
  /// container and pass through a few simple iterators whilst
  /// adding appropriate operator overloading and norms.
  template <typename _Type>
  class DenseVector
  {
  public:

    typedef typename std::vector<_Type>::iterator elt_iter;
    typedef typename std::vector<_Type>::const_iterator elt_citer;
    typedef typename std::vector<_Type>::reverse_iterator elt_riter;
    typedef typename std::vector<_Type>::const_reverse_iterator elt_criter;

    /// Constructor for a non-filled vector, to be filled by the user.
    DenseVector();

    /// Constructor with specified fill-in initialization
    /// \param fill Data to be initialised to each entry
    /// \param size The size of the vector to be instantiated
    DenseVector( const std::size_t& size, const _Type& fill );

    /// Construct a Noddy vector from a contiguous set of real data.
    /// This will be nasty if you pass the wrong pointer, but is
    /// useful in interfacing with external libraries
    /// \param size The number of elements in the vector.
    /// \param p A pointer to the start of the data.
    DenseVector( const std::size_t& size, const _Type* p );

    /// A templated implicitly converting copy constructor.
    /// Note this is specialised for double to double and
    /// std::complex<double> to std::complex<double> copy construction.
    /// \param source The DenseVector to be used in the initialising.
    template <typename _sourceType>
    DenseVector( const DenseVector<_sourceType>& source )
    {
      // size the current vector
      VEC.resize( source.size() );
      elt_iter p_local( VEC.begin() );
      for ( typename DenseVector<_sourceType>::elt_citer
            p_from = source.begin();
            p_from != source.end();
            ++p_from, ++p_local )
      {
        *p_local = *p_from;
      }
    }

    /// Copy assignment
    /// \param source Object to copy
    /// \return The new object
    DenseVector& operator=( const DenseVector& source );

    ~DenseVector();

    /// Pass through to the storage container.
    elt_iter begin()
    {
      return VEC.begin();
    }

    /// Pass through to the storage container.
    elt_riter rbegin()
    {
      return VEC.rbegin();
    }

    /// Pass through to the storage container.
    elt_citer begin() const
    {
      return VEC.begin();
    }

    /// Pass through to the storage container.
    elt_criter rbegin() const
    {
      return VEC.rbegin();
    }

    /// Pass through to the storage container.
    elt_iter end()
    {
      return VEC.end();
    }

    /// Pass through to the storage container.
    elt_citer end() const
    {
      return VEC.end();
    }

    /// Pass through to the storage container.
    elt_riter rend()
    {
      return VEC.rend();
    }

    /// Pass through to the storage container.
    elt_criter rend() const
    {
      return VEC.rend();
    }

    /// Operator overloading for addition
    /// \param x The dense vector to be added
    /// \return The sum of 'this' and x
    DenseVector<_Type> operator+( const DenseVector<_Type>& x ) const;

    /// Overloading for +
    /// \return +this
    DenseVector<_Type> operator+() const;

    /// Operator overloading for subtraction
    /// \param x The dense vector to be subtracted from 'this'
    /// \return The subtraction of x from 'this'
    DenseVector<_Type> operator-( const DenseVector<_Type>& x ) const;

    /// Overloading for -
    /// \return -this
    DenseVector<_Type> operator-() const;

    /// Operator overloading for scalar multiplication
    /// \param m The scalar multiplier
    /// \return The 'this' vector scaled by the constant 'm'
    DenseVector<_Type> operator*( const _Type& m ) const;

    /// Operaotr overloading for scalar division
    /// \param m The scalar divisor
    /// \return The 'this' vector divided by the constant 'm'
    DenseVector<_Type> operator/( const _Type& m ) const;

    /// Overloading of the [] operator
    /// \param i The index of the element to be accessed
    /// \return The element stored at index i (read only)
    const _Type& operator[] ( const std::size_t& i ) const;

    /// Overloading of the [] operator
    /// \param i The index of the element to be accessed
    /// \return The element  stored at index i (read/write)
    _Type& operator[] ( const std::size_t& i );

    /// Overloading *= for scalar multiplication
    /// \param m The scalar multiplier
    /// \return A reference to the 'this' vector multiplied by the scalar
    DenseVector<_Type>& operator*=( const _Type& m );

    /// Overloading /= for scalar multiplication
    /// \param m The scalar divisor
    /// \return A reference to the 'this' vector divided by the scalar
    DenseVector<_Type>& operator/=( const _Type& m );

    /// Overloading the -= operator
    /// \param x A dense vector that will be subtracted from 'this'
    /// \return A reference to the 'this' vector after subtraction by x
    DenseVector<_Type>& operator-=( const DenseVector<_Type>& x );

    /// Overloading the += operator
    /// \param x A dense vector that will be added to 'this'
    /// \return A reference to the 'this' vector after addition by x
    DenseVector<_Type>& operator+=( const DenseVector<_Type>& x );

    /// A pass-thru definition of push_back
    /// \param fill The data element to be pushed into the vector
    void push_back( const _Type& fill );

    /// A pass-thru definition of resize
    /// \param length The target length after resizing
    void resize( const std::size_t& length );

    /// A pass-thru definition of assign
    /// \param n The number of elements to assign
    /// \param elem The element copy to be used in the assign
    void assign( const std::size_t n, const _Type elem )
    {
      VEC.assign( n, elem );
    }

    /// A pass-thru definition of clear
    void clear();

    /// l1-norm.
    /// \return The square-root of the sum of the squares.
    double one_norm() const;

    /// l2-norm.
    /// No attention paid to possible overflow for large vectors.
    /// \return The square-root of the sum of the squares.
    double two_norm() const;

    /// Infinity norm.
    /// \return The maximum (abs) element in the vector.
    double inf_norm() const;

    /// Scale each element of the vector, equivalent to *=
    /// \param scale The value to scale each element by.
    void scale( const _Type& scale );

    /// Add a vector, element wise, equivalent to +=
    /// \param x The vector to be added to this object.
    void add( const DenseVector<_Type>& x );

    /// Subtract a vector, element wise, equivalent to -=
    /// \param x The vector to be subtracted from this object.
    void sub( const DenseVector<_Type>& x );

    /// A pass-thru definition to get the size of the vector.
    /// Since the vector is dense, the number of elements is the
    /// size.
    /// \return The number of elements in the vector.
    std::size_t size() const;

    /// Get the number of elements in the vector
    /// Since the vector is dense, the number of elements is the
    /// size.
    /// \return The number of elements in the vector.
    std::size_t nelts() const;

    /// Reserve space for the vector.
    /// \param n The number of elements to reserve space for
    void reserve( const std::size_t& n );

    /// Swap elements i and j.
    /// \param i The index of the element to swap.
    /// \param j The index of the other element to swap.
    void swap( const std::size_t& i, const std::size_t& j );

    /// Dump to std::cout
    void dump() const;

    /// Dump the contents to a file, each element on a separate line
    /// \param filename The name of the file to write to
    /// \param precision Precision of the output strings
    void dump_file( std::string filename, int precision = 10 ) const;
    // {
    //   std::ofstream dump;
    //   dump.open( filename.c_str() );
    //   dump.precision( precision );
    //   dump.setf( std::ios::showpoint );
    //   dump.setf( std::ios::showpos );
    //   dump.setf( std::ios::scientific );
    //   for ( std::size_t i = 0; i < VEC.size(); ++i )
    //   {
    //     dump << VEC[ i ] << "\n";
    //   }
    // }

  private:

    // friend classes that may access VEC directly
    friend class BandedMatrix<_Type>;

    // private storage of the data - encapsulated in std::vector
    std::vector<_Type> VEC;

  }
  ; // end class


  // INLINED METHODS BELOW - the dense vector class is used heavily, so
  // inlining simple access methods and operator overloads will give us
  // a speed improvement throughout.

  template <typename _Type>
  inline DenseVector<_Type>::~DenseVector()
  {}

  // specialise some copy constructors
  template <>
  template <>
  inline DenseVector<double>::DenseVector( const DenseVector<double>& source )
  {
    VEC = source.VEC;
  }

  template <>
  template <>
  inline DenseVector<std::complex<double> >::DenseVector( const DenseVector<std::complex<double> >& source )
  {
    VEC = source.VEC;
  }

  template <typename _Type>
  inline _Type& DenseVector<_Type>::operator[] ( const std::size_t& i )
  {
#ifdef PARANOID
    if ( ( i < 0 ) || ( i >= size() ) )
    {
      std::string problem;
      problem = " The DenseVector.operator[] method has a range error. \n";
      throw ExceptionRange( problem, size(), i );
    }
#endif
    return VEC[ i ];
  }

  template <typename _Type>
  inline const _Type& DenseVector<_Type>::operator[] ( const std::size_t& i ) const
  {
#ifdef PARANOID
    if ( ( i < 0 ) || ( i >= size() ) )
    {
      std::string problem;
      problem = " The DenseVector.operator[] method has a range error. \n";
      throw ExceptionRange( problem, size(), i );
    }
#endif
    return VEC[ i ];
  }


  template <typename _Type>
  inline void DenseVector<_Type>::push_back( const _Type& fill )
  {
    VEC.push_back( fill );
  }

  template <typename _Type>
  inline void DenseVector<_Type>::resize( const std::size_t& length )
  {
    VEC.resize( length );
  }

  template <typename _Type>
  inline void DenseVector<_Type>::clear()
  {
    VEC.clear();
  }

  template <typename _Type>
  inline DenseVector<_Type> DenseVector<_Type>::operator+() const
  {
    return * this;
  }

  template <typename _Type>
  inline std::size_t DenseVector<_Type>::size() const
  {
    return VEC.size();
  }

  template <typename _Type>
  inline std::size_t DenseVector<_Type>::nelts() const
  {
    return VEC.size();
  }

  template <typename _Type>
  inline void DenseVector<_Type>::reserve( const std::size_t& n )
  {
    VEC.reserve( n );
  }

  template <typename _Type>
  inline DenseVector<_Type> DenseVector<_Type>::operator*( const _Type& m ) const
  {
    DenseVector<_Type> temp( *this );
    temp *= m;
    return temp;
  }

  template <typename _Type>
  inline DenseVector<_Type>& DenseVector<_Type>::operator*=( const _Type& m )
  {
    // run thru VEC from begin to end, transforming into
    // VEC container starting from begin using functor
    std::transform( VEC.begin(), VEC.end(), VEC.begin(), scale_functor<_Type, _Type>( m ) );
    return *this;
  }

  template <typename _Type>
  inline DenseVector<_Type> DenseVector<_Type>::operator/( const _Type& m ) const
  {
    DenseVector<_Type> temp( *this );
    temp *= ( 1. / m );
    return temp;
  }

  template <typename _Type>
  inline DenseVector<_Type>& DenseVector<_Type>::operator/=( const _Type& m )
  {
    // run thru VECfrom begin to end, transforming into
    // VEC container starting from begin using functor
    std::transform( VEC.begin(), VEC.end(), VEC.begin(), scale_functor<_Type, _Type>( 1.0 / m ) );
    return *this;
  }

  template <typename _Type>
  inline DenseVector<_Type>& DenseVector<_Type>::operator-=( const DenseVector<_Type>& x )
  {
#ifdef PARANOID
    if ( x.size() != size() )
    {
      std::string problem;
      problem = " The DenseVector.operator-= method is trying to use \n";
      problem += " two vectors of unequal length \n";
      throw ExceptionGeom( problem, size(), x.size() );
    }
#endif
    std::transform( VEC.begin(), VEC.end(), x.begin(), VEC.begin(), std::minus<_Type>() );
    return *this;
  }

  template <typename _Type>
  inline DenseVector<_Type>& DenseVector<_Type>::operator+=( const DenseVector<_Type>& x )
  {
#ifdef PARANOID
    if ( x.size() != size() )
    {
      std::string problem;
      problem = " The DenseVector.operator+= method is trying to use \n";
      problem += " two vectors of unequal length \n";
      throw ExceptionGeom( problem, size(), x.size() );
    }
#endif
    std::transform( VEC.begin(), VEC.end(), x.begin(), VEC.begin(), std::plus<_Type>() );
    return *this;
  }

  template <typename _Type>
  inline DenseVector<_Type> DenseVector<_Type>::operator-( const DenseVector<_Type>& x ) const
  {
#ifdef PARANOID
    if ( x.size() != size() )
    {
      std::string problem;
      problem = " The DenseVector.operator- method is trying to use \n";
      problem += " two vectors of unequal length \n";
      throw ExceptionGeom( problem, size(), x.size() );
    }
#endif
    DenseVector<_Type> temp( *this );
    temp -= x;
    return temp;
  }

  template <typename _Type>
  inline DenseVector<_Type> DenseVector<_Type>::operator-() const
  {
    DenseVector<_Type> temp( *this );
    temp *= -1;
    return temp;
  }

  template <typename _Type>
  inline DenseVector<_Type> DenseVector<_Type>::operator+( const DenseVector<_Type>& x ) const
  {
#ifdef PARANOID
    if ( x.size() != size() )
    {
      std::string problem;
      problem = " The DenseVector.operator+ method is trying to use \n";
      problem += " two vectors of unequal length \n";
      throw ExceptionGeom( problem, size(), x.size() );
    }
#endif
    DenseVector<_Type> temp( *this );
    temp += x;
    return temp;
  }

  template <typename _Type>
  inline DenseVector<_Type>& DenseVector<_Type>::operator=( const DenseVector<_Type>& source )
  {
    if ( this == &source )
      return *this;
    VEC = source.VEC;
    return *this;
  }

}
#endif
