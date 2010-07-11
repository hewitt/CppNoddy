/// \file Functors.h
/// Some Function Objects that CppNoddy makes use of
/// in algorithms applied to STL containers.

#ifndef FUNCTORS_H
#define FUNCTORS_H

#include <cmath>

namespace CppNoddy
{

  /// A function object predicate that compares the absolute value
  /// of two elements and returns a true of elt1 < elt2.
  template <typename _Type>
  class abs_predicate
  {
  public:
    /// Return true if argument 1 < argument 2.
    /// \param elt1 First element to compare
    /// \param elt2 Second element to compare
    bool operator() ( _Type elt1, _Type elt2 )
    {
      return std::abs( elt2 ) > std::abs( elt1 );
    }
  };

  /// A function object predicate that first computes the absolute
  /// difference between two elements and a specified value, then
  /// compares these two resulting values. Used for computing which
  /// of two components is nearest to the specified value.
  template <typename _Type>
  class absDiff_predicate : public std::binary_function< _Type, _Type, bool>
  {
  public:
    /// \param value The target value that is to used.
    absDiff_predicate( _Type value ) : std::binary_function< _Type, _Type, bool>()
    {
      this -> target_elt = value;
    }

    /// Returns true of the first argument is closer to the specified
    /// (in constructor) value than the second argument.
    /// \param elt1 The first value to be used in the comparison
    /// \param elt2 The second value to be used in the comparison
    bool operator() ( _Type elt1, _Type elt2 )
    {
      return std::abs( elt2 - target_elt ) > std::abs( elt1 - target_elt );
    }

  private:
    _Type target_elt;
  };


  /// A function object used to accumulate the absolute value of a
  /// container.
  template <typename _Type>
  class absAdd_functor
  {
  public:
    double operator() ( double x, _Type y )
    {
      return x + std::abs( y );
    }
  };

  /// A function object used to accumulate the square of the
  /// absolute value of a container.
  template <typename _Type>
  class absSquareAdd_functor
  {
  public:
    double operator() ( double x, _Type y )
    {
      return x + std::pow( std::abs( y ), 2. );
    }
  };

  /// A unary pure function object that scales through multiplication.
  /// The object is type templated for any object that defines
  /// operator* and the MULTIPLIER is set in the functor constructor.
  template <class _containerType, typename _valueType>
  class scale_functor : public std::unary_function< _valueType, _containerType >
  {
  public:
    /// \param value The value to be used in the scale operation
    scale_functor( _valueType value ) : std::unary_function< _valueType, _containerType>()
    {
      MULTIPLIER = value;
    }

    /// \return The operator*( elt, value) operation where value
    /// is specified in the constructor
    _containerType operator() ( _containerType elt )
    {
      return elt * MULTIPLIER;
    }

  private:
    _valueType MULTIPLIER;
  };


} // end namespace

#endif // FUNCOBJS_H
