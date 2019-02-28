/// \file Residual_with_coords.h
/// A specification of a (double/complex) residual class that
/// not only defines a vector residual of a vector of state variables
/// but may also depend upon N additional (double) variables. This is
/// useful for the specification of boundary conditions in the PDE_IBVP
/// and PDE_double_IBVP classes where the boundary conditions are
/// dependent on the time/spatial location.

#ifndef RESIDUAL_WITH_COORDS_H
#define RESIDUAL_WITH_COORDS_H

#include <Residual.h>
#include <DenseVector.h>
#include <DenseMatrix.h>

namespace CppNoddy {
  /// A base class to be inherited by objects that define residuals
  template < typename _Type, typename _Xtype = double >
  class Residual_with_coords : public Residual<_Type> {
   public:
    /// Constructor for a 'square' residual object
    /// that is, N residuals for N unknowns.
    /// \param order The order of the residual vector
    /// \param ncoords The number of coordinates to store
    Residual_with_coords(const unsigned& order, const unsigned& ncoords);

    /// Constructor for a 'non-square' residual object
    /// that is, there are less residual constraints than unknowns.
    /// \param order The number of residuals
    /// \param nvars The number of unknowns/variables
    /// \param ncoords The number of coordinates to store
    Residual_with_coords(const unsigned& order, const unsigned& nvars, const unsigned& ncoords);

    /// An empty destructor
    virtual ~Residual_with_coords();

    /// General handle access to the coordinates
    /// \return A handle to the i-th coordinate
    _Xtype& coord(const unsigned& i);

    /// General handle access to the coordinates
    /// \return A handle to the i-th coordinate
    const _Xtype& coord(const unsigned& i) const;

   protected:

    /// The coordinates stored for this residual
    std::vector<_Xtype> coords;

  }
  ;  // end class

  template <typename _Type, typename _Xtype>
  inline _Xtype& Residual_with_coords<_Type, _Xtype>::coord(const unsigned& i) {
#ifdef PARANOID
    /// \todo check range on coord
#endif
    return coords[ i ];
  }

  template <typename _Type, typename _Xtype>
  inline const _Xtype& Residual_with_coords<_Type, _Xtype>::coord(const unsigned& i) const {
#ifdef PARANOID
#endif
    return coords[ i ];
  }

}   // end namespace

#endif
