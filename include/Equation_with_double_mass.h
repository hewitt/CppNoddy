/// \file Equation_with_double_mass.h
/// A templated class for equations that can be inherited from
/// to allow instantiation of PDE_double_IBVP objects using the resulting class.

#ifndef EQUATION_WITH_DOUBLE_MASS_H
#define EQUATION_WITH_DOUBLE_MASS_H

#include <Residual_with_coords.h>

namespace CppNoddy
{

  /// An equation object base class used in the PDE_double_IBVP class.
  /// An equation object is essentially a ('square') residual object (although
  /// it doesn't currently inherit) with some independent variable
  /// data members and access methods. By 'square' we mean that it defines
  /// N residuals and N state variables. In this case the equation also defines
  /// 2 mass matrices (amongst other data), which are used in the Crank-Nicolson
  /// time stepping of the PDE_double_IBVP class.
  template <typename _Type>
  class Equation_with_double_mass : public Residual_with_coords<_Type>
  {
  public:

    /// Constructor for equation class.
    /// \param order The order of the system
    explicit Equation_with_double_mass( const unsigned &order );

    /// An empty destructor, virtual since we have virtual methods.
    virtual ~Equation_with_double_mass();

    /// Update the Equation object for the current set of state variables
    /// \param state The state vector at which to set the equation object
    void update( const DenseVector<_Type> &state );

    /// Return a handle to the mass1 matrix member data
    const DenseMatrix<_Type>& mass1() const;

    /// Return a handle to the mass1 matrix member data
    const DenseMatrix<_Type>& mass2() const;

    /// Return the product of the Jacobian-of-the-mass1-matrix and a vector 'vec'
    /// when the equation has a given 'state'. The user should overload this
    /// if concerned about performance of the solver. If not overloaded, the
    /// default is to finite difference the Jacobian-of-the-mass-matrix.
    /// \param state The current state variables -- used for clarity when
    /// overloaded by the user instead of expecting the user to access the member data.
    /// \param vec The vector that will be multiplied by the Jacobian-of-the-mass-matrix
    /// \param h The resulting 2D matrix
    virtual void get_jacobian_of_mass1_mult_vector( const DenseVector<_Type> &state, const DenseVector<_Type> &vec, DenseMatrix<_Type> &h ) const;

    /// Return the product of the Jacobian-of-the-mass2-matrix and a vector 'vec'
    /// when the equation has a given 'state'. The user should overload this
    /// if concerned about performance of the solver. If not overloaded, the
    /// default is to finite difference the Jacobian-of-the-mass-matrix.
    /// \param state The current state variables -- used for clarity when
    /// overloaded by the user instead of expecting the user to access the member data.
    /// \param vec The vector that will be multiplied by the Jacobian-of-the-mass-matrix
    /// \param h The resulting 2D matrix
    virtual void get_jacobian_of_mass2_mult_vector( const DenseVector<_Type> &state, const DenseVector<_Type> &vec, DenseMatrix<_Type> &h ) const;

  protected:

    /// Define the mass matrix for parabolic problems in terms of the
    /// current state vector.
    /// \param state The current state vector.
    /// \param m The mass matrix.
    virtual void mass1( const DenseVector<_Type> &state, DenseMatrix<_Type> &m ) const
    {
      std::string problem;
      problem = "The equation::mass1 method has not been implemented.\n";
      problem += "You have to implement this method to define the equation.\n";
      throw ExceptionRuntime( problem );
    }

    /// Define the mass matrix for parabolic problems in terms of the
    /// current state vector.
    /// \param state The current state vector.
    /// \param m The mass matrix.
    virtual void mass2( const DenseVector<_Type> &state, DenseMatrix<_Type> &m ) const
    {
      std::string problem;
      problem = "The equation::mass2 method has not been implemented.\n";
      problem += "You have to implement this method to define the equation.\n";
      throw ExceptionRuntime( problem );
    }

  private:

    /// Mass matrices for the last state vector
    DenseMatrix<_Type> MASS1_AT_LAST_STATE;
    DenseMatrix<_Type> MASS2_AT_LAST_STATE;

  }
  ; // end class

  template <typename _Type>
  inline const DenseMatrix<_Type>& Equation_with_double_mass<_Type >::mass1() const
  {
    return MASS1_AT_LAST_STATE;
  }

  template <typename _Type>
  inline const DenseMatrix<_Type>& Equation_with_double_mass<_Type >::mass2() const
  {
    return MASS2_AT_LAST_STATE;
  }

} // end namespace

#endif
