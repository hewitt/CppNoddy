/// \file Equation_with_mass.h
/// A templated class for equations that can be inherited from
/// to allow instantiation of PDE_IBVP objects.

#ifndef EQUATION_WITH_MASS_H
#define EQUATION_WITH_MASS_H

#include <Residual_with_coords.h>

namespace CppNoddy
{

  /// An equation object base class used in the IBVP classes.
  /// An equation object is essentially a ('square') residual object (although
  /// it doesn't currently inherit) with an independent variable
  /// data member and access methods. By 'square' we mean that it defines
  /// N residuals and N state variables. In this case the equation also defines
  /// a mass matrix (amongst other data), which is used in the Crank-Nicolson
  /// time stepping of the PDE_IBVP class.
  template < typename _Type, typename _Xtype = double >
  class Equation_with_mass : public Residual_with_coords<_Type, _Xtype>
  {
  public:
    /// Constructor for equation class.
    /// \param order The order of the system
    explicit Equation_with_mass( const unsigned &order );

    /// An empty destructor, virtual since we have virtual methods.
    virtual ~Equation_with_mass();

    /// Update the Equation object for the current set of state variables
    /// \param state The state vector at which to set the equation object
    void update( const DenseVector<_Type> &state );

    /// Return a handle to the mass matrix
    const DenseMatrix<_Type>& mass() const;

    /// Return the product of the Jacobian-of-the-mass-matrix and a vector 'vec'
    /// when the equation has a given 'state'. The user should overload this
    /// if concerned about performance of the solver. If not overloaded, the
    /// default is to finite difference the Jacobian-of-the-mass-matrix.
    /// \param state The current state variables -- used for clarity when
    /// overloaded by the user instead of expecting the user to access the member data.
    /// \param vec The vector that will be multiplied by the Jacobian-of-the-mass-matrix
    /// \param h The resulting 2D matrix
    virtual void get_jacobian_of_mass_mult_vector( const DenseVector<_Type> &state, const DenseVector<_Type> &vec, DenseMatrix<_Type> &h ) const;

  protected:

    /// Define the mass matrix for parabolic problems in terms of the
    /// current state vector.
    /// \param x The current state vector.
    /// \param m The mass matrix.
    virtual void mass( const DenseVector<_Type> &x, DenseMatrix<_Type> &m ) const
    {
      std::string problem;
      problem = "The equation::mass method has not been implemented.\n";
      problem += "You have to implement this method to define the equation.\n";
      throw ExceptionRuntime( problem );
    }

  private:
    /// Mass matrix for the last state vector
    DenseMatrix<_Type> MASS_AT_LAST_STATE;

  }
  ; // end class


  template <typename _Type, typename _Xtype>
  inline const DenseMatrix<_Type>& Equation_with_mass<_Type, _Xtype>::mass() const
  {
    return MASS_AT_LAST_STATE;
  }

} // end namespace

#endif
