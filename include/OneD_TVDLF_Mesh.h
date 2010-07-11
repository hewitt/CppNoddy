/// \file OneD_TVDLF_Mesh.h
/// Specification of an object that represents a one dimensional
/// mesh for TVD LX methods.

#ifndef ONED_TVDLF_MESH_H
#define ONED_TVDLF_MESH_H

#include <vector>
#include <fstream>
#include <limits>

#include <Types.h>
#include <OneD_TVDLF_Elt.h>
#include <OneD_Node_Mesh.h>
#include <OneD_Hyperbolic_System.h>

namespace CppNoddy
{

  class OneD_TVDLF_Mesh
  {
    /// iterators for the vector of elements
    typedef std::vector<OneD_TVDLF_Elt> vector_of_elts;
    typedef vector_of_elts::const_iterator celt_iter;
    typedef vector_of_elts::iterator elt_iter;
    /// function pointer used in the initial conditions
    typedef void ( *fn_ptr ) ( const double&, DenseVector<double>& );

  public:

    /// Constructor for the Finite Volume Mesh using linear elements
    /// \param X A vector of nodal locations at which the element
    ///   FACES will positioned
    /// \param ptr A pointer to the hyperbolic system applied to this mesh
    /// \param init_ptr A pointer to a function that defines the initial conditions
    OneD_TVDLF_Mesh( const DenseVector<double>& X, OneD_Hyperbolic_System* ptr, fn_ptr init_ptr );

    /// Empty desctructor
    ~OneD_TVDLF_Mesh();

    /// Get the nodal positions in the middle of each element
    /// \return An NVector<double> of the nodal points
    DenseVector<double> get_mid_node_vector() const;

    /// Get the positions of the element faces
    /// \return An NVector<double> of the spatial points
    DenseVector<double> get_face_pos_vector() const;

    /// Set the limiter type to be applied in the slope
    /// values. 0 is no limiter, 1 is Lax-Wendroff, 2 is
    /// Beam-Warming, 3 is MC and 4 is Superbee.
    /// \param id The identifier of the limiter.
    void set_limiter( unsigned id );

    /// Update all the elements in this mesh to
    /// a new time step.
    /// \param CFL The CFL value to be used to determine the time step
    /// \param max_dt Do not take a time step larger than this irrespective
    /// of the CFL value
    double update( const double& CFL, const double& max_dt = std::numeric_limits<long double>::max() );

    /// Update all the elements in this mesh to a USER SPECIFIED time step.
    /// \param CFL The CFL value to be used to determine the time step
    /// \param t_end The time level to compute to
    void update_to( const double& CFL, const double& t_end );

    double update_to_red( const double& CFL, const double& max_dt )
    {
      // the black mesh slopes are set in the constructor
      // and at the end of every 'update' call

      // integrate the black mesh data onto the red mesh
      {
        elt_iter er = RED_ELTS.begin();
        while ( er != RED_ELTS.end() )
        {
          er -> set_Q_mid( er -> contributed_Q() / er -> get_dx() );
          ++er;
        }
      }

      // step thru the black elements to find a max time step
      double first_dt;
      {
        elt_iter eb = BLACK_ELTS.begin();
        first_dt = eb -> get_max_dt();
        ++eb;
        while ( eb != BLACK_ELTS.end() )
        {
          first_dt = std::min( first_dt, eb -> get_max_dt() );
          ++eb;
        }
        first_dt *= CFL;
      }
      if ( first_dt > max_dt )
      {
        std::cout << "Wanting to take a step of size " << first_dt << " but max_dt = " << max_dt << "\n";
        first_dt = max_dt;
      }

      calc_slopes( &RED_ELTS );
      // compute the fluxes through the element boundaries
      // in the red mesh, using the nodal values of the black mesh
      // then update the concentrations in the red mesh elements
      {
        elt_iter er = RED_ELTS.begin();
        DenseVector<double> flux_in_left( ORDER_OF_SYSTEM, 0.0 );
        DenseVector<double> flux_out_right( ORDER_OF_SYSTEM, 0.0 );
        // start with the left most elt & work out the flux in from the left
        er -> contributed_flux_in_left( first_dt, flux_in_left );
        while ( er != RED_ELTS.end() )
        {
          // work out the flux out to the right
          er -> contributed_flux_out_right( first_dt, flux_out_right );
          // we now have the flux difference
          DenseVector<double> deltaQ = ( flux_in_left - flux_out_right ) * first_dt / er -> get_dx();
          // contribution from the source integral
          {
            double x_mid( er -> get_x( 0.0 ) );
            DenseVector<double> slope( er -> get_slope() );
            DenseVector<double> q_mid( er -> get_Q( 0.0 ) );
            q_mid += ( er -> get_source_fn( 0.0 )
                       - er -> get_Jac_flux_fn( 0.0 ).multiply( slope ) ) * 0.5 * first_dt;
            DenseVector<double> r_mid( ORDER_OF_SYSTEM, 0.0 );
            er -> system_ptr -> source_fn( x_mid, q_mid, slope, r_mid );
            deltaQ += r_mid * first_dt;
          }
          er -> set_Q_mid( er -> get_Q( 0.0 ) + deltaQ );
          // the flux out right in this elt is the flux in left of the next one
          flux_in_left = flux_out_right;
          ++er;
        }
      }

      // compute the slopes using the specified limiter for the red elements
      calc_slopes( &RED_ELTS );
      MESH_TIME += first_dt;
      return first_dt;
    }

    double update_to_black( const double& CFL, const double& max_dt )
    {
      // integrate the red mesh data back onto the black mesh
      {
        elt_iter eb = BLACK_ELTS.begin();
        while ( eb != BLACK_ELTS.end() )
        {
          eb -> set_Q_mid( eb -> contributed_Q() / eb -> get_dx() );
          ++eb;
        }
      }
      // step thru the red elements to find a max time step
      double second_dt;
      {
        elt_iter er = RED_ELTS.begin();
        second_dt = er -> get_max_dt();
        ++er;
        while ( er != RED_ELTS.end() )
        {
          second_dt = std::min( second_dt, er -> get_max_dt() );
          ++er;
        }
        second_dt *= CFL;
      }
      if ( second_dt > max_dt )
      {
        second_dt = max_dt;
      }

      calc_slopes( &BLACK_ELTS );
      // compute the fluxes through the element boundaries
      // in the black mesh, using the nodal values of the red mesh
      // then update the concentrations in the black mesh elements
      {
        elt_iter eb = BLACK_ELTS.begin();
        DenseVector<double> flux_in_left( ORDER_OF_SYSTEM, 0.0 );
        DenseVector<double> flux_out_right( ORDER_OF_SYSTEM, 0.0 );
        // start with the left most elt & work out the flux in from the left
        eb -> contributed_flux_in_left( second_dt, flux_in_left );
        while ( eb != BLACK_ELTS.end() )
        {
          // work out the flux out to the right
          eb -> contributed_flux_out_right( second_dt, flux_out_right );
          // we now have the flux difference
          DenseVector<double> deltaQ = ( flux_in_left - flux_out_right ) * second_dt / eb -> get_dx();
          // contribution from the source integral
          {
            double x_mid( eb -> get_x( 0.0 ) );
            DenseVector<double> q_mid( eb -> get_Q( 0.0 ) );
            DenseVector<double> slope( eb -> get_slope() );
            q_mid += ( eb -> get_source_fn( 0.0 )
                       - eb -> get_Jac_flux_fn( 0.0 ).multiply( slope ) ) * 0.5 * second_dt;
            DenseVector<double> r_mid( ORDER_OF_SYSTEM, 0.0 );
            eb -> system_ptr -> source_fn( x_mid, q_mid, slope, r_mid );
            deltaQ += r_mid * second_dt;
          }
          eb -> set_Q_mid( eb -> get_Q( 0.0 ) + deltaQ );
          // the flux out right in this elt is the flux in left of the next one
          flux_in_left = flux_out_right;
          ++eb;
        }
      }

      // compute the slopes using the specified limiter
      // for the black elements
      calc_slopes( &BLACK_ELTS );
      MESH_TIME += second_dt;
      return second_dt;

    }

    /// Get a const reference to the time value for the current mesh.
    /// \return time The time level at which the data in the mesh applies.
    const double& get_time() const;

    /// Integrate the concentration values across the entire mesh.
    /// \return The values of the integral of each component.
    DenseVector<double> integrate() const;

    /// Get a OneD_Mesh<double> object containing the one dimensional
    /// data in the usual format.
    /// \param location Use "centre" for mid-elt values and "face_average" for
    /// the average values at the (discontinuous) element boundaries
    /// \param colour Which mesh to output, unless debugging, this should be "black"
    /// otherwise time values will be slightly out
    /// \return The required mesh object
    OneD_Node_Mesh<double> get_soln( std::string location = "centre", std::string colour = "black" );

    OneD_Node_Mesh<double> get_slope()
    {
      vector_of_elts* elts( get_elts_from_colour( "black" ) );
      DenseVector<double> X( get_mid_node_vector() );
      // the variables are the slope for each variable
      OneD_Node_Mesh<double> slope_mesh( X, ORDER_OF_SYSTEM );
      for ( std::size_t i = 0; i < X.size(); ++i )
      {
        slope_mesh.set_nodes_vars( i, ( *elts )[i].get_slope( ) );
      }
      return slope_mesh;
    }

  private:

    /// Loops through all the elements. For any elements on the
    /// boundary, we set the mid-nodal concentration to be the value
    /// on the wall, as specified by the edge_values method
    /// of the Hyperbolic sytem object.
    /// \param elts A vector of elts to be looped through
    void set_boundary_Q( vector_of_elts* elts );

    /// Given a choice of black or red mesh, return a pointer
    /// to the appropriate vector of elements.
    /// \param mesh_colour A mesh identifier (black or red)
    /// \return A pointer to a std::vector of elements
    vector_of_elts* get_elts_from_colour( std::string mesh_colour );

    /// Use the appropriate limiter to approximate the slope in each
    /// element in the mesh. The slopes will be sent down to the element
    /// objects.
    /// \param elt_vector The vector of elements to set the slope for
    void calc_slopes( vector_of_elts* elt_vector );

    /// Compute the left-face difference.
    /// We dont worry about which way is up/down wind because
    /// we will only use limiters that treat up/down symmetrically.
    /// \param e An iterator reference to an element
    /// \return An approximation to the leftward differece
    DenseVector<double> left_diff( elt_iter e );

    /// Compute the right-face difference.
    /// We dont worry about which way is up/down wind because
    /// we will only use limiters that treat up/down symmetrically.
    /// \param k An iterator reference to an element
    /// \return An approximation to the rightward differece
    DenseVector<double> right_diff( elt_iter e );

    /// Sign of a double.
    /// \param a The value to return the sign of
    /// \return The sign of the value
    int sgn( double a );

    /// Returns the minimum of the absolute value of the two arguments
    /// if both arguments are of the same sign.
    /// \param a First element to compare
    /// \param b Second element to compare
    /// \return a if |a|<|b| and a*b > 0, b if |b|<|a| and a*b >0,
    ///   otherwise return zero
    double minmod( double a, double b );

    /// Returns the maximum of the absolute value of the two arguments
    /// if both arguments are of the same sign.
    /// \param a First element to compare
    /// \param b Second element to compare
    /// \return a if |a|>|b| and a*b > 0, b if |b|>|a| and a*b >0,
    ///   otherwise return zero
    double maxmod( double a, double b );

    /// A vector version of the minmod operator
    /// \param A A vector to compare
    /// \param B A vector to compare
    /// \return A component-wise minmod vector, i.e., each component
    /// of the returned vector is the minmod of the components of A & B.
    DenseVector<double> minmod( DenseVector<double> A, DenseVector<double> B );

    /// A vector version of the maxmod operator
    /// \param A A vector to compare
    /// \param B A vector to compare
    /// \return A component-wise maxmod vector, i.e., each component
    /// of the returned vector is the maxmod of the components of A & B.
    DenseVector<double> maxmod( DenseVector<double> A, DenseVector<double> B );

    /// The time level of this mesh solution
    double MESH_TIME;
    /// Slope limiter method
    unsigned LIMITER;
    /// order of the conservative system
    std::size_t ORDER_OF_SYSTEM;

    /// An STL vector of linear elements -- the black mesh
    vector_of_elts BLACK_ELTS;
    /// An STL vector of linear elements -- the red mesh
    vector_of_elts RED_ELTS;
    /// function pointer to a function that defines the initial distribution
    fn_ptr p_Q_INIT;

  };

}

#endif // end OneD_TVDLF_Mesh_H
