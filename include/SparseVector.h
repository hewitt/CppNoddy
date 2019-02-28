/// \file SparseVector.h
/// A templated SparseVector class -- a sparse, variable size, vector object.

#ifndef SPARSEVECTOR_H
#define SPARSEVECTOR_H

#include <map>
#include <Exceptions.h>

#if defined(PETSC_D) || defined(PETSC_Z)
  #include "petsc.h"
#endif

namespace CppNoddy {

  /// An SparseVector class -- a sparse vector object.
  /// This is templated but intended ONLY for double or
  /// std::complex<double>.
  template <typename _Type>
  class SparseVector {
    typedef typename std::map< std::size_t, _Type >::const_iterator citer;
    typedef typename std::map< std::size_t, _Type >::iterator iter;

   public:

    /// Constructor for a non-filled vector, to be filled by the user.
    /// \param max A maximum size for the sparse vector, used in geometry
    /// checking.
    explicit SparseVector(const std::size_t& max);

    /// Copy constructor.
    /// \param source The source object to be copied
    SparseVector(const SparseVector& source);

    /// Assignment operator.
    /// \param source The source object for the assignment
    /// \return The newly assigned object
    SparseVector& operator=(const SparseVector& source);

    /// Not a full iterator implementation -- just a pass
    /// through to the storage container.
    /// Return an iterator pointing to the beginning of
    /// the encpsulated STL map storage
    iter begin();

    /// Not a full iterator implementation -- just a pass
    /// through to the storage container.
    /// Return a constant iterator pointing to the beginning of
    /// the encapsulated STL map storage
    citer begin() const;

    /// Not a full iterator implementation -- just a pass
    /// through to the storage container.
    /// Return an iterator pointing to the end of the
    /// encapsulated STL map storage
    iter end();

    /// Not a full iterator implementation -- just a pass
    /// through to the storage container.
    /// Return a constant iterator pointing to the end of
    /// the encapsulated STL map storage
    citer end() const;

    /// Resize the maximum length of the sparse vector.
    /// The maximum length is used in geometry checking.
    /// \param length The new maximum length of the vector
    void resize(const std::size_t& length);

    /// Remove all elements from the sparse vector
    void clear();

    /// Erase an element from the vector
    /// \param pos An std::map< std::size_t, _Type >::iterator to the map element
    void erase(const iter& pos);

    /// Erase an element from th vector
    /// \param index The index of the element to be erased
    void erase(const std::size_t& index);

    /// Swap elements i and j.
    /// \param i The index of the element to swap.
    /// \param j The index of the other element to swap.
    void swap(const std::size_t& i, const std::size_t& j);

    /// Swap ALL elements with those of another SparseVector
    /// \param X The sparse vector to exchange with
    void swap(SparseVector<_Type>& X);

    /// Find the (max) size of the vector.
    std::size_t size() const;

    /// Find the number of non-zero elements in the vector.
    std::size_t nelts() const;

    /// Set an element of the vector
    /// \param i The index to be accessed
    /// \return A reference to the element
    _Type& set(const std::size_t& i);

    /// Get an element of the vector
    /// \param i The index to be accessed
    /// \return The contents of the element, including
    /// a zero if the element has not been set
    const _Type& get(const std::size_t& i) const;

    /// Equivalent to the 'get' method
    const _Type& operator[](const std::size_t& i) const;

    /// Equivalent to the 'set' method
    _Type& operator[](const std::size_t& i);

    /// Operator overloading for sparse vector addition
    /// \param X A sparse vector to be added
    /// \return The sum of this and the passed vector X
    SparseVector<_Type> operator+(const SparseVector<_Type>& X) const;

    /// Overloading for +
    /// \return +this
    SparseVector<_Type> operator+() const;

    /// Operator overloading for sparse vector subtraction
    /// \param X A sparse vector to be subtracted
    /// \return The subtraction of this and the passed vector X
    SparseVector<_Type> operator-(const SparseVector<_Type>& X) const;

    /// Overloading for -
    /// \return The negative of the vector
    SparseVector<_Type> operator-() const;

    /// Overloading multiplication for a scalar
    /// \param m The scalar to multiply by
    /// \result The scalar multiplication of the vector
    SparseVector<_Type> operator*(const _Type& m) const;

    /// Overloading *= for scalar multiplication
    /// \param m The scalar to multiply by
    SparseVector<_Type>& operator*=(const _Type& m);

    /// Overloading -= for sparse vectors
    /// \param X The sparse vector to be subtracted
    /// \return Subtracts X from the 'this' vector
    SparseVector<_Type>& operator-=(const SparseVector<_Type>& X);

    /// Overloading += for sparse vectors
    /// \param X The sparse vector to be added
    /// \return Adds X to the 'this' vector
    SparseVector<_Type>& operator+=(const SparseVector<_Type>& X);

    /// Look for an element index
    /// \param i The index of the element to search for
    citer find(std::size_t i) const {
      citer pos;
      pos = MAP_VEC.find(i);
      return pos;
    }

    /// A dot product.
    /// \param x The vector to be "dotted" with.
    const _Type dot(const SparseVector<_Type>& x) const;

    /// l1-norm.
    /// \return The square-root of the sum of the squares divided by number of elts.
    double one_norm() const;

    /// l2-norm.
    /// No attention paid to possible overflow for large vectors.
    /// \return The square-root of the sum of the squares divided by number of elts.
    double two_norm() const;

    /// Infinity norm.
    /// \return The maximum (abs) element in the vector.
    double inf_norm() const;

    /// Scale each element of the vector.
    /// \param scale The value to scale each element by.
    void scale(const _Type& scale);

    /// Add a vector, element wise.
    /// \param X The vector to be added to this object.
    void add(const SparseVector<_Type>& X);

    /// Subtract a vector, element wise.
    /// \param X The vector to be subtracted from this object.
    void sub(const SparseVector<_Type>& X);

    /// Find the index of the element NEAREST in value to that specified
    /// \param value The value to search for.
    /// \return The index of the element with value that minimises the difference.
    std::size_t nearest_index(const _Type& value) const;

    /// Find the index of the maximum element in the vector
    /// \return The index of the maximum element.
    std::size_t maxabs_index() const;

    /// Find the index of the maximum element in the vector
    /// \return The index of the maximum element.
    std::size_t minabs_index() const;

    /// Output the sparse vector's contents.
    void dump() const;

#if defined(PETSC_D) || defined(PETSC_Z)

    void get_petsc( PetscScalar* storage, PetscInt* cols );

    /*{
      citer pos;
      std::size_t i(0);
      //
      // start at the begining of this row
      pos = MAP_VEC.begin();
      do {
        // store the value and column
        storage[ i ] = pos -> second;
        cols[ i ] = pos -> first;
        // increment the iterator
        ++pos;
        // increment the array index
        ++i;
        } while(pos != MAP_VEC.end());
        }*/
    
    
#endif

    
   private:
    // private storage of the data - encapsulated in std::vector
    std::map< std::size_t, _Type > MAP_VEC;
    _Type ZERO;
    std::size_t MAX_SIZE;

  }
  ; // end class


  // INLINE METHODS BELOW

  template <typename _Type>
  inline typename std::map< std::size_t, _Type >::iterator SparseVector<_Type>::begin() {
    return MAP_VEC.begin();
  }

  template <typename _Type>
  inline typename std::map< std::size_t, _Type >::const_iterator SparseVector<_Type>::begin() const {
    return MAP_VEC.begin();
  }

  template <typename _Type>
  inline typename std::map< std::size_t, _Type >::iterator SparseVector<_Type>::end() {
    return MAP_VEC.end();
  }

  template <typename _Type>
  inline typename std::map< std::size_t, _Type >::const_iterator SparseVector<_Type>::end() const {
    return MAP_VEC.end();
  }

  template <typename _Type>
  inline _Type& SparseVector<_Type>::operator[](const std::size_t& i) {
#ifdef PARANOID
    if(i >= size()) {
      std::string problem;
      problem = " The SparseVector.operator[] method is trying to access \n";
      problem += " outside the container. \n";
      throw ExceptionRange(problem, size(), i);
    }
#endif
    return MAP_VEC[ i ];
  }

  template <typename _Type>
  inline const _Type& SparseVector<_Type>::operator[](const std::size_t& i) const {
#ifdef PARANOID
    if(i >= size()) {
      std::string problem;
      problem = " The SparseVector.operator[] method is trying to access \n";
      problem += " outside the container. \n";
      throw ExceptionRange(problem, size(), i);
    }
#endif
    return get(i);
  }

  template <typename _Type>
  inline _Type& SparseVector<_Type>::set(const std::size_t& i) {
#ifdef PARANOID
    if(i >= size()) {
      std::string problem;
      problem = " The SparseVector.set method is trying to access \n";
      problem += " outside the container. \n";
      throw ExceptionRange(problem, size(), i);
    }
#endif
    return MAP_VEC[ i ];
  }

  template <typename _Type>
  inline const _Type& SparseVector<_Type>::get(const std::size_t& i) const {
#ifdef PARANOID
    if(i >= size()) {
      std::string problem;
      problem = " The SparseVector.get method is trying to access \n";
      problem += " outside the container. \n";
      throw ExceptionRange(problem, size(), i);
    }
#endif
    citer pos;
    pos = MAP_VEC.find(i);
    if(pos != MAP_VEC.end()) {
      return pos -> second;
    } else {
      return ZERO;
    }
  }

  template <typename _Type>
  inline std::size_t SparseVector<_Type>::size() const {
    return MAX_SIZE;
  }

  template <typename _Type>
  inline std::size_t SparseVector<_Type>::nelts() const {
    return MAP_VEC.size();
  }

  template <typename _Type>
  inline void SparseVector<_Type>::erase(const std::size_t& index) {
    MAP_VEC.erase(index);
  }

  template <typename _Type>
  inline void SparseVector<_Type>::erase(const iter& pos) {
    MAP_VEC.erase(pos);
  }

  template <typename _Type>
  inline void SparseVector<_Type>::swap(SparseVector<_Type>& X) {
#ifdef PARANOID
    if(X.size() != size()) {
      std::string problem;
      problem = " The SparseVector.swap method is trying to use \n";
      problem += " two vectors of unequal length \n";
      throw ExceptionGeom(problem, size(), X.size());
    }
#endif
    MAP_VEC.swap(X.MAP_VEC);
  }

  template <typename _Type>
  inline SparseVector<_Type> SparseVector<_Type>::operator+() const {
    return * this;
  }
  

} // end namespace

#endif

