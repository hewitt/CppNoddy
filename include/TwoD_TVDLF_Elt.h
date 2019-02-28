/// \file TwoD_TVDLF_Elt.h
/// Specification of a two dimensional linear element for
/// use in a TVD Lax-Friedrichs scheme.

#ifndef TWOD_TVDLF_ELT_H
#define TWOD_TVDLF_ELT_H

#include <list>
#include <set>
#include <algorithm>

#include <TwoD_Hyperbolic_System.h>

//#define DEBUG

namespace CppNoddy {

  /// A Linear Element class.
  class TwoD_TVDLF_Elt {

    typedef std::vector<TwoD_TVDLF_Elt> vector_of_elts;
    typedef vector_of_elts::const_iterator celt_iter;
    typedef vector_of_elts::iterator elt_iter;

    /// A contributory element's data in the projection scheme.
    struct contribution {
      TwoD_TVDLF_Elt* elt_ptr;
      DenseVector<double> sw, ne;
      std::set< int > face_indices;
    };

   public:

    /// Construct a linear rectangular element
    /// \param west The western face location
    /// \param east The eastern face location
    /// \param south The western face location
    /// \param north The eastern face location
    /// \param ptr A pointer to the hyperbolic system
    /// \param flag A boolean indicator to specify an elt with an external face
    /// \param faces A set of indices of the external faces
    TwoD_TVDLF_Elt(double west, double east,
                   double south, double north,
                   TwoD_Hyperbolic_System* ptr,
                   bool flag = false, std::set<int> faces = std::set<int>());

    /// An empty destructor
    ~TwoD_TVDLF_Elt();

    /// Test to see if a face index of an element is external to the mesh.
    /// \param index The face index to test (0,1,2,3) for (S,E,N,W)
    /// \return True if the face is external to the mesh.
    bool face_is_external(const int& index) const;

    /// Test to see if a face index of an element is internal to the mesh.
    /// \param index The face index to test (0,1,2,3) for (S,E,N,W)
    /// \return True if the face is internal to the mesh.
    bool face_is_internal(const int& index) const;

    /// Each element stores pointers to any adjacent elements in the
    /// four compass directions 0,1,2,3 for S,E,N,W
    /// \param index The index of the direction to set the pointer for
    /// \param ptr A pointer to the element in the chosen direction
    void set_ptrs(const int& index, TwoD_TVDLF_Elt* ptr);

    /// Return a pointer to an element in the given compass direction
    /// \param index The appropriate index of the direction 0,1,2,3 for
    ///   S,E,N,W
    /// \return The pointer to the element
    TwoD_TVDLF_Elt* get_ptrs(const int& index) const;

    /// Add a contribution to this element from another element.
    /// \param ptr A pointer to the element that contributes to this one
    /// \param sw The SW local coordinate in the contributory element
    /// \param ne The NE local coordinate in the contributory element
    /// \param indices A set of indices of faces that this element contributes to
    void add_contribution(TwoD_TVDLF_Elt* ptr,
                          const DenseVector<double>& sw,
                          const DenseVector<double>& ne,
                          std::set< int > indices);

    /// Dump the details of this element.
    void dump() const;

    /// Reset the contributions to this element.
    void clear_contributions();

    /// Find the integral contributions to this black/red element from the
    /// corresponding contributing red/black element.
    /// \return The integrated value of all components
    DenseVector<double> contributed_Q() const;

    /// Adds the contribution of each face's flux to the correction
    /// for this element unless it has already been added by a flux
    /// computation across the same face in an adjacent elt.
    /// \param dt The time step over which the flux is to be computed
    void add_flux_contributions(const double& dt);

    /// Compute the flux into this element across the western face.
    /// \param dt The time step over which the flux is to be computed
    /// \param flux_in_west A vector flux that is overwritten by this method
    ///    and returned containing the appropriate flux values.
    void contributed_flux_in_west(const double &dt, DenseVector<double> &flux_in_west);

    /// Compute the flux into this element across the southern face.
    /// \param dt The time step over which the flux is to be computed
    /// \param flux_in_south A vector flux that is overwritten by this method
    ///    and returned containing the appropriate flux values.
    void contributed_flux_in_south(const double &dt, DenseVector<double> &flux_in_south);

    /// Compute the flux out of this element across the eastern face.
    /// \param dt The time step over which the flux is to be computed
    /// \param flux_out_east A vector flux that is overwritten by this method
    ///    and returned containing the appropriate flux values.
    void contributed_flux_out_east(const double &dt, DenseVector<double> &flux_out_east);

    /// Compute the flux out of this element across the northern face.
    /// \param dt The time step over which the flux is to be computed
    /// \param flux_out_north A vector flux that is overwritten by this method
    ///    and returned containing the appropriate flux values.
    void contributed_flux_out_north(const double &dt, DenseVector<double> &flux_out_north);

    /// Get the local coordinate that corresponds to a given global coordinate.
    /// If the global coordinate is outside this element, an exception is thrown.
    /// \param x The global coordinate
    /// \return The local coordinate
    DenseVector<double> get_s(const DenseVector<double> &x) const;

    /// Get the global position of a local coordinate in this element.
    /// \param s The local coordinate
    /// \return The global position
    DenseVector<double> get_x(DenseVector<double> s) const;

    /// Get the global position of a central node in this element.
    /// \return The global position of the central node
    DenseVector<double> get_x_mid() const;

    /// Get the size of the element as a vector (delta_x, delta_y).
    /// \return The separation of the rectangles faces (delta_x, delta_y)
    DenseVector<double> get_dx() const;

    /// Get the area of the element
    /// \return The area of the element
    double get_dA() const;

    /// Set the external flag.
    /// \param flag is the value to set
    void set_external_flag(bool flag);

    /// Get the external flag.
    /// \return The value of the flag (true if this is an elt with an external face)
    bool get_external_flag() const;

    /// Get a set of indices of faces that are external
    /// \return A set of integers (0,1,2,3 for south, east, north, west)
    ///   indicating which faces are external.
    std::set<int> get_external_faces() const;

    /// Get the value of the 'concentration' stored in this element.
    /// \param s The local coordinate in the elt at which Q is requested
    /// \return The concentration value of the elt for each component
    DenseVector<double> get_Q(const DenseVector<double>& s) const;

    /// Get the integral of Q over a sub-element
    /// \param sw A vector containing the SW local coordinates
    /// \param ne A vector containing the NE local coordinates
    /// \return The integral of the concentration over the local range
    DenseVector<double> get_int_Q(const DenseVector<double>& sw, const DenseVector<double>& ne) const;

    /// Add an external face to the stored stl::set.
    /// \param i The index of the external face 0,1,2,3 for S,E,N,W.
    void add_external_face(const int& i);

    /// Remove an external face to the stored stl::set.
    /// \param i The index of the external face 0,1,2,3 for S,E,N,W.
    void remove_external_face(const int& i);

    /// Set the value of the vector 'concentration' stored in this
    /// element. To second order, this is the value of the
    /// concentration at the mid-point of the element.
    /// \param value The vector value to be assigned to the
    ///  concentration member data Q
    void set_Q_mid(const DenseVector<double>& value);

    /// Get the nodal concentration
    /// \return The concentration at the central node
    DenseVector<double> get_Q_mid() const;

    /// Set the x-slope vector for this element.
    /// \param value The slope vector to be set
    void set_slope_x(const DenseVector<double>& value);

    /// Set the y-slope vector for this element.
    /// \param value The slope vector to be set
    void set_slope_y(const DenseVector<double>& value);

    /// The concentration is approximated by a linear function
    /// in each element. The slope of this function in the x-direction
    /// is found through this method.
    /// \return The value of the 'slope_x' member data
    DenseVector<double> get_slope_x() const;

    /// The concentration is approximated by a linear function
    /// in each element. The slope of this function in the y-direction
    /// is found through this method.
    /// \return The value of the 'slope_x' member data
    DenseVector<double> get_slope_y() const;

    /// Get the Taylor series expansion of the concentration
    /// for a given time step and local coordinate for this element
    /// \param s The local coordinate
    /// \param dt The time step
    DenseVector<double> get_Q_Taylor_series(const DenseVector<double> &s, const double &dt) const;

    /// Evaluate the contribution to this element by the hyperbolic
    /// system's source function over a given time step using a
    /// mid-point in time evaluation.
    /// \param dt The time step over which source function is to be integrated
    /// \return The total contribution
    DenseVector<double> get_source_contribution(const double &dt) const;

    /// Get the maximum allowable time step for this element by using
    /// information about the size of the element and the maximum wave
    /// speed set in the conservative system. The time step is guranteed
    /// to satisfy the CFL < 0.5 constraint in the two principle directions.
    /// \return The maximum time step for this element
    double get_max_dt() const;

    /// Get the flux function in the x direction evaluated for
    /// the concentration value stored in this elt.
    /// \param s The local vector coordinate at which to evaluate
    /// \return The flux function vector
    DenseVector<double> get_flux_fn_x(const DenseVector<double>& s) const {
      DenseVector<double> f(p_system -> get_order(), 0.0);
      p_system -> flux_fn_x(get_x(s), get_Q(s), f);
      return f;
    }

    /// Get the flux function in the y direction evaluated for
    /// the concentration value stored in this elt.
    /// \param s The local vector coordinate at which to evaluate
    /// \return The flux function vector
    DenseVector<double> get_flux_fn_y(const DenseVector<double>& s) const {
      DenseVector<double> g(p_system -> get_order(), 0.0);
      p_system -> flux_fn_y(get_x(s), get_Q(s), g);
      return g;
    }

    /// Get the Jacobian of the x flux function evaluated for the
    /// concentration value stored in this elt.
    /// \param s The local coordinate at which to evaluate
    /// \return The Jacobian matrix
    DenseMatrix<double> get_Jac_flux_fn_x(const DenseVector<double>& s) const {
      DenseMatrix<double> J(p_system -> get_order(), p_system -> get_order(), 0.0);
      p_system -> Jac_flux_fn_x(get_x(s), get_Q(s), J);
      return J;
    }

    /// Get the Jacobian of the y flux function evaluated for the
    /// concentration value stored in this elt.
    /// \param s The local coordinate at which to evaluate
    /// \return The Jacobian matrix
    DenseMatrix<double> get_Jac_flux_fn_y(const DenseVector<double>& s) const {
      DenseMatrix<double> J(p_system -> get_order(), p_system -> get_order(), 0.0);
      p_system -> Jac_flux_fn_y(get_x(s), get_Q(s), J);
      return J;
    }


    /// Get the source function evaluated for the
    /// concentration value stored in this elt
    /// \param s The local coordinate at which to evaluate
    /// \return The source function value
    DenseVector<double> get_source_fn(const DenseVector<double>& s) const {
      DenseVector<double> r(p_system -> get_order(), 0.0);
      p_system -> source_fn(get_x(s), get_Q(s), r);
      return r;
    }


    TwoD_Hyperbolic_System* p_system;

   private:

    /// Since any face is common to at least two elements this
    /// is a private method that allows any element to set a
    /// flux contribution in another (adjacent) elt.
    /// \param face_index The index (0,1,2,3) of the face that this
    /// contribution is for
    /// \param dq The contribution to be made
    void add_to_delta_Q(const int& face_index, const DenseVector<double>& dq);

    /// a list of elts that contribute to this one
    std::list< contribution > CONT_LIST;
    /// the concentration at the center of this elt
    DenseVector<double> Q;
    /// the 2 slopes for this elt
    DenseVector<double> SLOPE_X;
    DenseVector<double> SLOPE_Y;
    /// keep track of the flux into the elt so far
    DenseVector<double> DELTA_Q;
    /// the global coordinates of the 4 faces
    double SOUTH, WEST;
    /// elt sizes
    double DX, DY;
    /// a vector of iterators to neighbouring elts
    std::vector< TwoD_TVDLF_Elt* > p_ELTS;
    /// a set of indices of external faces
    std::set<int> EXTERNAL_FACES;
    /// a vector of booleans that are set when a face's flux contribution has
    /// been computed
    std::vector<bool> FLUX_FACE_DONE;
    /// a boolean to indicate that this elt has external faces
    bool EXTERNAL;

  };

  // used too many times for inlining apparently
  inline DenseVector<double> TwoD_TVDLF_Elt::get_Q(const DenseVector<double>& s) const {
    return Q + (SLOPE_X * s[ 0 ] * DX
                + SLOPE_Y * s[ 1 ] * DY) / 2;
  }

  inline void TwoD_TVDLF_Elt::set_Q_mid(const DenseVector<double>& value) {
    Q = value;
  }

  inline DenseVector<double> TwoD_TVDLF_Elt::get_Q_mid() const {
    return Q;
  }

  inline void TwoD_TVDLF_Elt::set_slope_x(const DenseVector<double>& value) {
    SLOPE_X = value;
  }

  inline void TwoD_TVDLF_Elt::set_slope_y(const DenseVector<double>& value) {
    SLOPE_Y = value;
  }

  inline DenseVector<double> TwoD_TVDLF_Elt::get_slope_x() const {
    return SLOPE_X;
  }

  inline DenseVector<double> TwoD_TVDLF_Elt::get_slope_y() const {
    return SLOPE_Y;
  }

  inline DenseVector<double> TwoD_TVDLF_Elt::get_int_Q(const DenseVector<double>& s_sw, const DenseVector<double>& s_ne) const {
    const double sub_dx(DX * (s_ne[ 0 ] - s_sw[ 0 ]) / 2);
    const double sub_dy(DY * (s_ne[ 1 ] - s_sw[ 1 ]) / 2);
    return get_Q((s_ne + s_sw) / 2) * sub_dx * sub_dy;
  }

  inline double TwoD_TVDLF_Elt::get_max_dt() const {
    DenseVector<double> c(2, 0.0);
    DenseVector<double> s(2, 0.0);
    p_system -> max_charac_speed(get_x(s), get_Q_mid(), c);
    return std::min(DX / c[ 0 ], DY / c[ 1 ]);
  }

} // CppNoddy namespace

#endif // TwoD_TVDLF_Elt_H
