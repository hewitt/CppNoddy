/// \file OneD_TVDLF_Elt.cpp
/// Implementation of a one dimensional linear element for
/// use in a TVD Lax-Friedrichs scheme.

#include <list>

#include <OneD_Hyperbolic_System.h>
#include <OneD_TVDLF_Elt.h>
#include <Types.h>

namespace CppNoddy {

  OneD_TVDLF_Elt::OneD_TVDLF_Elt(double a, double b,
                                 OneD_Hyperbolic_System* ptr,
                                 bool flag, int index) {
#ifdef DEBUG
    std::cout << "DEBUG: Constructing a OneD_TVDLF_Elt object between "
              << a << " and " << b << "\n";
#endif
    system_ptr = ptr;
    LEFT = a;
    RIGHT = b;
    EXTERNAL = flag;
    EXTERNAL_FACE_I = index;
    SLOPE = DenseVector<double>(system_ptr -> get_order(), 0.0);
    Q = DenseVector<double>(system_ptr -> get_order(), 0.0);
  }

  OneD_TVDLF_Elt::~OneD_TVDLF_Elt()
  {}

  void OneD_TVDLF_Elt::add_contribution(OneD_TVDLF_Elt* ptr, const DenseVector<double>& s_range, int index) {
    contribution cont;
    cont.elt_ptr = ptr;
    cont.s = s_range;
    cont.face_index = index;
    // push the contribution data into a list
    CONT_LIST.push_back(cont);
  }

  DenseVector<double> OneD_TVDLF_Elt::contributed_Q() {
    std::list< contribution >::iterator c = CONT_LIST.begin();
    // start with zero
    DenseVector<double> sum(system_ptr -> get_order(), 0.0);
    // loop over contributions
    while(c != CONT_LIST.end()) {
      // get integral of each
      sum += c -> elt_ptr -> get_int_Q(c -> s[0], c -> s[1]);
      ++c;
    }
    return sum;
  }

  void OneD_TVDLF_Elt::contributed_flux_in_left(const double &dt, DenseVector<double> &flux_in_left, const double &current_time) {
    std::list< contribution >::iterator c = CONT_LIST.begin();
    // loop over all contributing elements
    while(c != CONT_LIST.end()) {
      // find contributions to the left face
      if(c -> face_index == -1) {
        // which local coord in the contributing elt?
        const double s = c -> s[ 0 ];
        DenseVector<double> q_half_step(c -> elt_ptr -> get_Q(s));
        // add half step correction
        q_half_step += (c -> elt_ptr -> get_source_fn(s)
                        - c -> elt_ptr -> get_Jac_flux_fn(s).multiply(c -> elt_ptr -> SLOPE))
                       * 0.5 * dt;
        //DenseVector<double> x( 1, 0.0 ); x[ 0 ] = c -> elt_ptr -> get_x( s );
        double x(c -> elt_ptr -> get_x(s));
        // evaluate the flux at this Q value
        c -> elt_ptr -> system_ptr -> flux_fn(x, q_half_step, flux_in_left);
      }
      ++c;
    }
    // treat the user specified external faces separately
    if(EXTERNAL) {
      // we are in an external elt
      c = CONT_LIST.begin();
      while(c != CONT_LIST.end()) {
        // external faces are external for the current elt *and* its contributions
        if(c -> elt_ptr -> get_external_flag()) {
          if(EXTERNAL_FACE_I < 0) {
            // LH boundary
            DenseVector<double> q_left(c -> elt_ptr -> get_Q(-1.0));
            // convert x at boundary to a vector
            //DenseVector<double> x_left( 1, 0.0 ); x_left[ 0 ] = c -> elt_ptr -> get_x( -1.0 );
            double x_left(c -> elt_ptr -> get_x(-1.0));
            DenseVector<double> q_left_half_step(q_left);
            //compute flux at the half time step
            q_left_half_step += (c -> elt_ptr -> get_source_fn(-1.0)
                                 - c -> elt_ptr -> get_Jac_flux_fn(-1.0).multiply(c -> elt_ptr -> SLOPE))
                                * 0.5 * dt;
            // get the edge values
            system_ptr -> edge_values(-1, x_left, q_left_half_step, current_time + 0.5*dt);
            // get the flux
            system_ptr -> flux_fn(x_left, q_left_half_step, flux_in_left);
          }
        }
        ++c;
      }
    }
  }

  void OneD_TVDLF_Elt::contributed_flux_out_right(const double &dt, DenseVector<double> &flux_out_right, const double &current_time) {
    std::list< contribution >::iterator c = CONT_LIST.begin();
    // loop over all contributing elements
    while(c != CONT_LIST.end()) {
      // find contributions to the right face
      if(c -> face_index == 1) {
        // which local coord in the contributing elt?
        const double s = c -> s[ 1 ];
        DenseVector<double> q_half_step(c -> elt_ptr -> get_Q(s));
        // add half step correction
        q_half_step += (c -> elt_ptr -> get_source_fn(s)
                        - c -> elt_ptr -> get_Jac_flux_fn(s).multiply(c -> elt_ptr -> SLOPE))
                       * 0.5 * dt;
        //DenseVector<double> x( 1, 0.0 ); x[ 0 ] = c -> elt_ptr -> get_x( s );
        double x(c -> elt_ptr -> get_x(s));
        // evaluate the flux at this Q value
        c -> elt_ptr -> system_ptr -> flux_fn(x, q_half_step, flux_out_right);
      }
      ++c;
    }
    // treat the user specified external faces separately
    if(EXTERNAL) {
      // we are in an external elt
      c = CONT_LIST.begin();
      while(c != CONT_LIST.end()) {
        // external faces are external for the current elt *and* its contributions
        if(c -> elt_ptr -> get_external_flag()) {
          if(EXTERNAL_FACE_I > 0) {
            // RH boundary
            DenseVector<double> q_right(c -> elt_ptr -> get_Q(1.0));
            // convert x at boundary to a vector
            //DenseVector<double> x_right( 1, 0.0 ); x_right[ 0 ] = c -> elt_ptr -> get_x( 1.0 );
            double x_right(c -> elt_ptr -> get_x(1.0));
            DenseVector<double> q_right_half_step(q_right);
            // compute flux at the half time step
            q_right_half_step += (c -> elt_ptr -> get_source_fn(1.0)
                                  - c -> elt_ptr -> get_Jac_flux_fn(1.0).multiply(c -> elt_ptr -> SLOPE))
                                 * 0.5 * dt;
            // get the edge values
            system_ptr -> edge_values(1, x_right, q_right_half_step, current_time + 0.5*dt);
            // get the flux
            system_ptr -> flux_fn(x_right, q_right_half_step, flux_out_right);
          }
        }
        ++c;
      }
    }
  }

  DenseVector<double> OneD_TVDLF_Elt::get_flux_fn(const double s) const {
    DenseVector<double> f(system_ptr -> get_order(), 0.0);
    double x(get_x(s));
    system_ptr -> flux_fn(x, get_Q(s), f);
    return f;
  }

  DenseMatrix<double> OneD_TVDLF_Elt::get_Jac_flux_fn(const double s) const {
    DenseMatrix<double> J(system_ptr -> get_order(), system_ptr -> get_order(), 0.0);
    double x(get_x(s));
    system_ptr -> Jac_flux_fn(x, get_Q(s), J);
    return J;
  }

  DenseVector<double> OneD_TVDLF_Elt::get_source_fn(const double s) const {
    DenseVector<double> r(system_ptr -> get_order(), 0.0);
    double x(get_x(s));
    system_ptr -> source_fn(x, get_Q(s), get_slope(), r);
    return r;
  }

  double OneD_TVDLF_Elt::get_max_dt() const {
    // left face & right face speeds
    return get_dx() / std::max(std::abs(system_ptr -> max_charac_speed(get_Q(-1.0))),
                               std::abs(system_ptr -> max_charac_speed(get_Q(1.0))));
  }

} // CppNoddy namespace
