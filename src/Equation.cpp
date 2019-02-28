/// \file Equation.cpp
/// Implementation for an equations class that can be inherited from
/// to allow instantiation of ODE objects using the resulting class.

#include <Equation.h>
#include <Residual_with_coords.h>

namespace CppNoddy {

  template <typename _Type, typename _Xtype>
  Equation<_Type, _Xtype>::Equation(const unsigned& order) :
    Residual_with_coords<_Type, _Xtype>(order, 1) {
  }

  template <typename _Type, typename _Xtype>
  Equation<_Type, _Xtype>::~Equation() {
  }

  // the required templated versions are:
  template class Equation<double>
  ;
  template class Equation<std::complex<double> >
  ;
  template class Equation<std::complex<double>, std::complex<double> >
  ;

} // end namespace
