/// \file Equation.h
/// A templated class for equations that can be inherited from
/// to allow instantiation of ODE/PDE objects using the resulting class.
/// An equation class is simply a square residual class with one
/// extra coordinate for the ODE_BVP solver.

#ifndef EQUATION_H
#define EQUATION_H

#include <Residual_with_coords.h>
#include <DenseVector.h>
#include <DenseMatrix.h>

namespace CppNoddy {

  /// An equation object base class used in the BVP/IVP classes.
  /// An equation object is essentially a ('square') residual object
  /// with an independent variable
  /// data member and access methods. By 'square' we mean that it defines
  /// N residuals and N state variables.
  template < typename _Type, typename _Xtype = double >
  class Equation : public Residual_with_coords<_Type, _Xtype> {
   public:
    /// Constructor for equation class.
    /// \param order The order of the system
    Equation(const unsigned& order);

    /// An empty destructor, virtual since we have virtual methods.
    virtual ~Equation();

  }
  ; // end class

} // end namespace

#endif
