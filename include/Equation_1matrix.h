/// \file Equation_1matrix.h
/// A templated class for equations that can be inherited from
/// to allow instantiation of PDE_IBVP objects (amongst others).

#ifndef EQUATION_1MATRIX_H
#define EQUATION_1MATRIX_H

#include <Residual_with_coords.h>

namespace CppNoddy
{

  /// An equation object base class used in the IBVP classes (and others).
  /// An equation object is essentially a ('square') residual object with an independent variable
  /// data member and access methods. By 'square' we mean that it defines
  /// N residuals and N state variables. The equation is defined
  /// using an NxN matrix that multiplies the derivative of unknowns and a residual RHS.
  template < typename _Type, typename _Xtype = double >
  class Equation_1matrix : public Residual_with_coords<_Type, _Xtype>
  {
  public:
    /// Constructor for equation class.
    /// \param order The order of the system
    explicit Equation_1matrix( const unsigned &order );

    /// An empty destructor, virtual since we have virtual methods.
    virtual ~Equation_1matrix();

    /// Update the Equation object for the current set of state variables
    /// \param state The state vector at which to set the equation object
    void update( const DenseVector<_Type> &state );

    /// Return a handle to the matrix
    const DenseMatrix<_Type>& matrix0() const;

    /// Return the product of the Jacobian-of-the-matrix and a vector 'vec'
    /// when the equation has a given 'state'. The user should overload this
    /// if concerned about performance of the solver. If not overloaded, the
    /// default is to finite difference the Jacobian-of-the-matrix.
    /// \param state The current state variables -- used for clarity when
    /// overloaded by the user instead of expecting the user to access the member data.
    /// \param vec The vector that will be multiplied by the Jacobian-of-the-matrix
    /// \param h The resulting 2D matrix
    virtual void get_jacobian_of_matrix0_mult_vector( const DenseVector<_Type> &state, 
       const DenseVector<_Type> &vec, DenseMatrix<_Type> &h ) const;


  protected:

    /// Define the matrix in terms of the current state vector.
    /// \param x The current state vector.
    /// \param m The matrix.
    virtual void matrix0( const DenseVector<_Type> &x, DenseMatrix<_Type> &m ) const
    {
      std::string problem;
      problem = "The equation::matrix0 method has not been implemented!\n";
      problem += "You have to implement this method to define the equation.\n";
      throw ExceptionRuntime( problem );
    }

  private:
    /// Matrix0 evaluated for the last state vector
    DenseMatrix<_Type> MATRIX0_AT_LAST_STATE;

  }
  ; // end class

  template <typename _Type, typename _Xtype>
  inline const DenseMatrix<_Type>& Equation_1matrix<_Type, _Xtype>::matrix0() const
  {
    return MATRIX0_AT_LAST_STATE;
  }

} // end namespace

#endif
