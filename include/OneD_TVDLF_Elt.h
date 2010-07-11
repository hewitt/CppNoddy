/// \file OneD_TVDLF_Elt.h
/// Specification of a one dimensional linear element for
/// use in a TVD Lax-Friedrichs scheme.

#ifndef ONED_TVDLF_ELT_H
#define ONED_TVDLF_ELT_H

#include <list>

#include <OneD_Hyperbolic_System.h>

namespace CppNoddy
{

  /// A Linear Element class.
  class OneD_TVDLF_Elt
  {

    struct contribution
    {
      OneD_TVDLF_Elt* elt_ptr;
      DenseVector<double> s;
      int face_index;
    };

  public:

    /// Construct a linear element
    /// \param a The left hand edge
    /// \param b The right hand edge
    /// \param ptr A pointer to the hyperbolic system for this element
    /// \param flag A boolean indicator to specify an elt with an external face
    /// \param index Index of the external face if flag=true
    OneD_TVDLF_Elt( double a, double b, OneD_Hyperbolic_System* ptr,
                    bool flag = false, int index = 0 );

    /// An empty destructor
    ~OneD_TVDLF_Elt();

    /// Allow the connection of this element to another element
    /// \param ptr A pointer to another element of the same type
    /// \param s_range A vector range of local coordinates (s_min,s_max) in
    /// this 1D case.
    /// \param index An index to indicate which face of the current elt
    /// this contribution is associated with for flux computations -1=left, +1=right, 0=none
    void add_contribution( OneD_TVDLF_Elt* ptr, const DenseVector<double>& s_range, int index );

    /// Clear the record of contributions to this element.
    void clear_contributions()
    {
      CONT_LIST.clear();
    }

    /// Find the integral contributions to this black/red element from the
    /// contributed red/black element.
    /// \return The integrated value of all components
    DenseVector<double> contributed_Q();

    /// Compute the flux into this element through the left face over a given
    /// time step.
    /// \param dt The time step to compute the flux over
    /// \param flux_in_left A vector to return the flux components in
    void contributed_flux_in_left( const double &dt, DenseVector<double> &flux_in_left );

    /// Compute the flux out of this element through the right face over a given
    /// time step.
    /// \param dt The time step to compute the flux over
    /// \param flux_out_right A vector to return the flux components in
    void contributed_flux_out_right( const double &dt, DenseVector<double> &flux_out_right );

    /// Convert a global coordinate into a local coordinate within this element
    /// If the global coordinate is outside the element, then we throw an exception.
    /// \param x The global coordinate to convert
    double get_s( double x ) const
    {
#ifdef PARANOID
      if ( ( x < LEFT ) || ( x > RIGHT ) )
      {
        std::string problem;
        problem = " The OneD_TVDLF_Elt::get_s method has been called for an element \n";
        problem += " whose range does not bracket the given global coordinate.\n";
        throw ExceptionRuntime( problem );
      }
#endif
      return -1.0 + 2 * ( x - LEFT ) / ( RIGHT - LEFT );
    }

    /// Get the global position of the local coordinate
    /// \param s The local coordinate
    /// \return The global position
    double get_x( double s ) const;

    /// Get the size of the element
    /// \return The separation of the left and right faces
    double get_dx() const;

    /// Get the external flag.
    /// \return The value of the flag (true if this is an elt with an external face)
    bool get_external_flag() const;

    /// Set the external flag.
    /// \param flag is the value to set
    void set_external_flag( bool flag );

    /// Get the exterior boundary number
    /// \return An integer boundary indicator +1 for right, -1 for left.
    int get_external_face_i() const;

    /// Set the local coordinate of the external face.
    /// Only one face can be external.
    /// \param i The index of the external face -1=left, +1=right
    void set_external_face_i( int i );

    /// Get the value of the 'concentration' stored in this
    /// element.
    /// \param s The local coordinate in the elt at which Q is requested
    /// \return The concentration value of the elt for each component
    DenseVector<double> get_Q( double s ) const;

    /// Get the integral of Q over a sub-element
    /// \param s0 The left local coordinate
    /// \param s1 The right local coordinate
    /// \return The integral of the concentration over the local range
    DenseVector<double> get_int_Q( const double s0, const double s1 );

    /// Set the value of the 'concentration' stored in this
    /// element. To second order, this is the value of the
    /// concentration at the mid-point of the element.
    /// \param value The value to be assigned to the
    ///  concentration member data Q
    void set_Q_mid( const DenseVector<double> value );

    /// The concentration is approximated by a linear function
    /// in each element. The slope of this function is set
    /// through this method.
    /// \param value The value to be assigned to the member data 'slope'
    void set_slope( const DenseVector<double> value );

    /// The concentration is approximated by a linear function
    /// in each element. The slope of this function is got
    /// through this method.
    /// \return The value of the 'slope' member data
    DenseVector<double> get_slope() const;

    /// Get the flux function evaluated for the concentration
    /// value stored in this elt.
    /// \param s The local coordinate at which to evaluate
    /// \return The flux function vector
    DenseVector<double> get_flux_fn( const double s ) const;

    /// Get the Jacobian of the flux function evaluated for the
    /// concentration value stored in this elt.
    /// \param s The local coordinate at which to evaluate
    /// \return The Jacobian matrix
    DenseMatrix<double> get_Jac_flux_fn( const double s ) const;

    /// Get the value of the source function term in this element at
    /// a given local coordinate.
    /// \param s The local coordinate to evaluate at
    DenseVector<double> get_source_fn( const double s ) const;

    /// Get the maximum allowable time step for this element by using
    /// information about the size of the element and the maximum wave
    /// speed set in the conservative system. The time step is guranteed
    /// to satisfy the CFL < 0.5 constraint.
    /// \return The maximum time step.
    double get_max_dt() const;

    OneD_Hyperbolic_System* system_ptr;

  private:

    std::list< contribution > CONT_LIST;

    DenseVector<double> Q;
    DenseVector<double> SLOPE;
    double LEFT, RIGHT;
    int EXTERNAL_FACE_I;

    bool EXTERNAL;

  };


  inline double OneD_TVDLF_Elt::get_dx() const
  {
    return RIGHT - LEFT;
  }

  inline double OneD_TVDLF_Elt::get_x( double s ) const
  {
    return ( LEFT + RIGHT + s * get_dx() ) / 2;
  }

  inline bool OneD_TVDLF_Elt::get_external_flag() const
  {
    return EXTERNAL;
  }

  inline int OneD_TVDLF_Elt::get_external_face_i() const
  {
    return EXTERNAL_FACE_I;
  }

  inline DenseVector<double> OneD_TVDLF_Elt::get_Q( double s ) const
  {
    return Q + SLOPE * s * ( RIGHT - LEFT ) / 2;
  }

  inline DenseVector<double> OneD_TVDLF_Elt::get_int_Q( const double s0, const double s1 )
  {
    // linear elt so integral is 0.5 * dx * ( f(s0) + f(s1) )
    return ( get_Q( s0 ) + get_Q( s1 ) ) * 0.5 * ( ( s1 - s0 ) * get_dx() / 2 );
  }

  inline void OneD_TVDLF_Elt::set_external_flag( bool flag )
  {
    EXTERNAL = flag;
  }

  inline void OneD_TVDLF_Elt::set_external_face_i( int i )
  {
    EXTERNAL_FACE_I = i;
  }

  inline void OneD_TVDLF_Elt::set_Q_mid( const DenseVector<double> value )
  {
    Q = value;
  }

  inline void OneD_TVDLF_Elt::set_slope( const DenseVector<double> value )
  {
    SLOPE = value;
  }

  inline DenseVector<double> OneD_TVDLF_Elt::get_slope() const
  {
    return SLOPE;
  }

} // CppNoddy namespace

#endif // OneD_TVDLF_Elt_H
