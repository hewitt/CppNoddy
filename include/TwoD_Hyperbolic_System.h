/// \file TwoD_Hyperbolic_System.h

#ifndef TWOD_HYPERBOLIC_SYSTEM_H
#define TWOD_HYPERBOLIC_SYSTEM_H

#include <Types.h>
#include <Uncopyable.h>
#include <Timer.h>

namespace CppNoddy
{
  /// A class to represent a two-dimensional hyperbolic system of equations.
  class TwoD_Hyperbolic_System : private Uncopyable
  {

  public:

    /// \param order The order of the hyperbolic system
    explicit TwoD_Hyperbolic_System( const unsigned& order ) : ORDER_OF_SYSTEM( order )
    {}

    /// An empty destructor, virtual since we have virtual methods.
    virtual ~TwoD_Hyperbolic_System()
    {}

    /// A virtual flux function for the x-derivative
    /// \param x The vector position
    /// \param q The unknowns
    /// \param f The flux function
    virtual void flux_fn_x( const DenseVector<double> &x, const DenseVector<double> &q, DenseVector<double> &f ) const
    {
      std::string problem;
      problem = "The Hyperbolic_Conservative_System::flux_fn_x method has not been implemented.\n";
      problem += "You have to implement this method to define the system.\n";
      throw ExceptionRuntime( problem );
    }

    /// A virtual flux function for the y-derivative
    /// \param x The vector position
    /// \param q The unknowns
    /// \param f The flux function
    virtual void flux_fn_y( const DenseVector<double> &x, const DenseVector<double> &q, DenseVector<double> &f ) const
    {
      std::string problem;
      problem = "The Hyperbolic_Conservative_System::flux_fn_y method has not been implemented.\n";
      problem += "You have to implement this method to define the system.\n";
      throw ExceptionRuntime( problem );
    }

    /// A virtual function function to define the Jacobian of the
    /// x-flux function. The default method uses first-order finite
    /// differencing to compute the Jacobian if not otherwise specified
    /// by the user.
    /// \param x The vector position
    /// \param q The unknowns
    /// \param J The Jacobian of the flux function
    virtual void Jac_flux_fn_x( const DenseVector<double> &x, const DenseVector<double> &q, DenseMatrix<double> &J ) const
    {
      /// cheap & nasty differencing
      const double delta( 1.e-8 );
      DenseVector<double> state( q );
      DenseVector<double> temp1( ORDER_OF_SYSTEM, 0.0 );
      DenseVector<double> temp2( ORDER_OF_SYSTEM, 0.0 );
      flux_fn_x( x, state, temp1 );
      // default is to FD the Jacobian
      for ( std::size_t i = 0; i < ORDER_OF_SYSTEM; ++i )
      {
        state[ i ] += delta;
        flux_fn_x( x, state, temp2 );
        state[ i ] -= delta;
        J.set_col( i, ( temp2 - temp1 ) / delta );
      }
    }

    /// A virtual function function to define the Jacobian of the
    /// y-flux function. The default method uses first-order finite
    /// differencing to compute the Jacobian if not otherwise specified
    /// by the user.
    /// \param x The vector position
    /// \param q The unknowns
    /// \param J The Jacobian of the flux function
    virtual void Jac_flux_fn_y( const DenseVector<double> &x, const DenseVector<double> &q, DenseMatrix<double> &J ) const
    {
      /// cheap & nasty differencing
      const double delta( 1.e-8 );
      DenseVector<double> state( q );
      DenseVector<double> temp1( ORDER_OF_SYSTEM, 0.0 );
      DenseVector<double> temp2( ORDER_OF_SYSTEM, 0.0 );
      flux_fn_y( x, state, temp1 );
      // default is to FD the Jacobian
      for ( std::size_t i = 0; i < ORDER_OF_SYSTEM; ++i )
      {
        state[ i ] += delta;
        flux_fn_y( x, state, temp2 );
        state[ i ] -= delta;
        J.set_col( i, ( temp2 - temp1 ) / delta );
      }
    }

    /// A virtual method that is used to bound the characteristic speed in both directions.
    /// It determines the time step in order that the CFL constraint holds. It is sometimes
    /// required to specify a speed that is a function of position and to even specify the
    /// two-directional speeds. eg. your mesh may actually be w.r.t. a polar coordinate
    /// system in which the physical mesh size grows linearly with radius, so one of the
    /// speeds should be scaled accordingly.
    /// \param x The global coordinate
    /// \param q The unknowns
    /// \param c The speed
    virtual void max_charac_speed( const DenseVector<double> &x, const DenseVector<double> &q, DenseVector<double> &c ) const
    {
      std::string problem;
      problem = "The Hyperbolic_Conservative_System::max_charac_speed method has not\n";
      problem += "been implemented. You have to implement this method to define the system.\n";
      throw ExceptionRuntime( problem );
    }

    /// Define the edge boundary conditions.
    /// \param face_index An index for the face:
    /// 0,1,2,3 for S,E,N,W on the TwoD_TVDLF_Mesh
    /// \param x The global position vector
    /// \param q Specify the unknowns specified along the face
    /// \return A vector that defines which components have been set as inflow
    virtual std::vector<bool> edge_values( const int& face_index, const DenseVector<double>& x, DenseVector<double>& q ) const
    {
      return std::vector<bool>( ORDER_OF_SYSTEM, false );
    }

    /// Define the edge boundary condition slopes.
    /// These default to zero unless otherwise over-ridden.
    /// \param face_index An index for the face:
    /// 0,1,2,3 for S,E,N,W on the TwoD_TVDLF_Mesh
    /// \param x The global position vector
    /// \param sigma_n Specify the slope normal to the face
    virtual void edge_slopes( const int& face_index, const DenseVector<double>& x, DenseVector<double>& sigma_n ) const
    {
    }

    virtual void source_fn( const DenseVector<double> &x, const DenseVector<double> &q, DenseVector<double>& r ) const
    {
    }

    unsigned get_order()
    {
      return ORDER_OF_SYSTEM;
    }

  protected:

    /// The order of the system of equations
    const unsigned ORDER_OF_SYSTEM;
  }
  ; // end class

} // end namespace

#endif
