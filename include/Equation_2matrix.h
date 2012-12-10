/// \file Equation_2matrix.h
/// A templated class for equations that can be inherited from
/// to allow instantiation of PDE_double_IBVP objects (amongst others).

#ifndef EQUATION_2MATRIX_H
#define EQUATION_2MATRIX_H

#include <Equation_1matrix.h>

namespace CppNoddy
{

  /// An equation object base class used in the PDE_double_IBVP class.
  /// An equation object is essentially a ('square') residual object (although
  /// it doesn't currently inherit) with some independent variable
  /// data members and access methods. By 'square' we mean that it defines
  /// N residuals and N state variables. In this case the equation also defines
  /// 2 matrices (amongst other data). This inherits from the Equation_1matrix
  /// and adds the functionality for the additional matrix.

  template < typename _Type, typename _Xtype = double >
  class Equation_2matrix : public Equation_1matrix<_Type, _Xtype>
  {    
  public:

    /// Constructor for equation class.
    /// \param order The order of the system
    explicit Equation_2matrix( const unsigned &order );

    /// An empty destructor, virtual since we have virtual methods.
    virtual ~Equation_2matrix();

    /// Update the Equation object for the current set of state variables
    /// \param state The state vector at which to set the equation object
    void update( const DenseVector<_Type> &state );

    /// Return a handle to the matrix member data
    const DenseMatrix<_Type>& matrix1() const;

    /// Return the product of the Jacobian-of-the-matrix and a vector 'vec'
    /// when the equation has a given 'state'. The user should overload this
    /// if concerned about performance of the solver. If not overloaded, the
    /// default is to finite difference the Jacobian-of-the-matrix.
    /// \param state The current state variables -- used for clarity when
    /// overloaded by the user instead of expecting the user to access the member data.
    /// \param vec The vector that will be multiplied by the Jacobian-of-the-matrix
    /// \param h The resulting 2D matrix
    virtual void get_jacobian_of_matrix1_mult_vector( const DenseVector<_Type> &state, const DenseVector<_Type> &vec, DenseMatrix<_Type> &h ) const;

  protected:

    /// Define the matrix in terms of the current state vector.
    /// \param state The current state vector.
    /// \param m The matrix.
    virtual void matrix1( const DenseVector<_Type> &state, DenseMatrix<_Type> &m ) const
    {
      std::string problem;
      problem = "The equation::matrix1 method has not been implemented.\n";
      problem += "You have to implement this method to define the equation.\n";
      throw ExceptionRuntime( problem );
    }

  private:

    /// Matrix for the last state vector
    DenseMatrix<_Type> MATRIX1_AT_LAST_STATE;

  }
  ; // end class

  template <typename _Type, typename _Xtype>
  inline const DenseMatrix<_Type>& Equation_2matrix<_Type, _Xtype>::matrix1() const
  {
    return MATRIX1_AT_LAST_STATE;
  }

} // end namespace

#endif
