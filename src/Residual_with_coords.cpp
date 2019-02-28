/// \file Residual_with_coords.cpp
/// A specification of a (double/complex) residual class that
/// not only defines a vector residual of a vector of state variables
/// but may also depend upon N additional (double) variables. This is
/// useful for the specification of boundary conditions in the PDE_IBVP
/// and PDE_double_IBVP classes where the boundary conditions are
/// dependent on the time/spatial location.

#include <Residual_with_coords.h>

namespace CppNoddy {
  template <typename _Type, typename _Xtype>
  Residual_with_coords<_Type, _Xtype>::Residual_with_coords(const unsigned& order, const unsigned& ncoords) : Residual<_Type>(order) {
    coords = std::vector<_Xtype>(ncoords, 0.0);
  }

  template <typename _Type, typename _Xtype>
  Residual_with_coords<_Type, _Xtype>::Residual_with_coords(const unsigned& order, const unsigned& nvars, const unsigned& ncoords) : Residual<_Type>(order, nvars) {
    coords = std::vector<_Xtype>(ncoords, 0.0);
  }

  template <typename _Type, typename _Xtype>
  Residual_with_coords<_Type, _Xtype>::~Residual_with_coords() {
  }

  // the templated versions we require are:
  template class Residual_with_coords<double>
  ;
  template class Residual_with_coords<std::complex<double> >
  ;
  template class Residual_with_coords<std::complex<double>, std::complex<double> >
  ;

}   // end namespace
