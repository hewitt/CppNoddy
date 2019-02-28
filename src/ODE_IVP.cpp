/// \file ODE_IVP.cpp
/// Implementation of an \f$n\f$-th order system
/// of ODEs that form the IVP:
/// \f[ \underline{ \dot f} (t) = \underline R ( \underline f (t), t) \,, \f]
/// where \f$ \underline f (0) \f$ is known.

#include <string>

#include <Types.h>
#include <Equation.h>
#include <OneD_Node_Mesh.h>
#include <ODE_IVP.h>
#include <Exceptions.h>
#include <Utility.h>

namespace CppNoddy {

  template <typename _Type>
  ODE_IVP<_Type>::ODE_IVP(Equation<_Type > *ptr,
                          const double &x1, const double &x2,
                          const std::size_t &num_of_points) :
    X_INIT(x1),
    X_FINAL(x2),
    H_INIT((x2 - x1) / num_of_points),
    N(num_of_points),
    p_EQUATION(ptr),
    STORE_EVERY(1) {
    p_EQUATION -> coord(0) = X_INIT;
  }

  template <typename _Type>
  ODE_IVP<_Type>::~ODE_IVP()
  {}

  template <typename _Type>
  DenseVector<_Type> ODE_IVP<_Type>::shoot4(DenseVector<_Type> u) {

    double x = X_INIT;
    const double h = H_INIT;
    const double hby2 = h / 2.;
    const double hby6 = h / 6.;
    const double hby3 = h / 3.;
    const int order = u.size();

    DenseVector<_Type> z(order, 0.0), k1(order, 0.0), k2(order, 0.0), k3(order, 0.0), k4(order, 0.0);

    DenseVector< double > coords;
    std::vector< DenseVector<_Type> > values;

    coords.push_back(x);
    values.push_back(u);

    for(unsigned i = 0; i < N; i++) {
      // k1 = F(u,x)
      p_EQUATION -> coord(0) = x;
      p_EQUATION -> residual_fn(u, k1);
      z = u + k1 * hby2;

      x += hby2;

      // k2 = F(z,xhh)
      p_EQUATION -> coord(0) = x;
      p_EQUATION -> residual_fn(z, k2);
      z = u + k2 * hby2;

      // k3 = F(z,xhh)
      p_EQUATION -> coord(0) = x;
      p_EQUATION -> residual_fn(z, k3);
      z = u + k3 * h;

      x += hby2;

      // k4 = F(z,xh)
      p_EQUATION -> coord(0) = x;
      p_EQUATION -> residual_fn(z, k4);
      u += k1 * hby6 + k2 * hby3 + k3 * hby3 + k4 * hby6;

      if(i % STORE_EVERY == 0) {
        coords.push_back(x);
        values.push_back(u);
      }

    } //for loop stepping across domain

    // construct the solution mesh stored in this object
    SOLN = OneD_Node_Mesh<_Type>(coords, p_EQUATION -> get_order());
    for(unsigned i = 0; i < coords.size(); ++i) {
      // fill mesh
      SOLN.set_nodes_vars(i, values[i]);
    }
    return u;
  }

  template <typename _Type>
  DenseVector<_Type> ODE_IVP<_Type>::shoot45(DenseVector<_Type> u, const double& tol, const double& h_init) {
    bool ok(false);
    unsigned step = 0;
    double x = X_INIT;
    double h = h_init;
    double c, diff;

    static const double X2 = 1. / 4.;
    static const double X3 = 3. / 8.;
    static const double X4 = 12. / 13.;
    static const double X5 = 1.;
    static const double X6 = 1. / 2.;

    static const double W21 = 1. / 4.;

    static const double W31 = 3. / 32.;
    static const double W32 = 9. / 32.;

    static const double W41 = 1932. / 2197.;
    static const double W42 = -7200. / 2197.;
    static const double W43 = 7296. / 2197.;

    static const double W51 = 439. / 216.;
    static const double W52 = -8.;
    static const double W53 = 3680. / 513.;
    static const double W54 = -845. / 4104;

    static const double W61 = -8. / 27.;
    static const double W62 = 2.;
    static const double W63 = -3544. / 2565.;
    static const double W64 = 1859. / 4104.;
    static const double W65 = -11. / 40.;

    static const double U1 = 25. / 216.;
    static const double U3 = 1408. / 2565.;
    static const double U4 = 2197. / 4104.;
    static const double U5 = -1. / 5.;

    static const double Z1 = 16. / 135.;
    static const double Z3 = 6656. / 12825.;
    static const double Z4 = 28561. / 56430.;
    static const double Z5 = -9. / 50.;
    static const double Z6 = 2. / 55.;

    const unsigned order = u.size();

    DenseVector<_Type> z(order, 0.0), e(order, 0.0), k1(order, 0.0),
                k2(order, 0.0), k3(order, 0.0), k4(order, 0.0), k5(order, 0.0), k6(order, 0.0);

    DenseVector< double > coords;
    std::vector< DenseVector<_Type> > values;

    coords.push_back(x);
    values.push_back(u);

    do {
      step += 1;
      // k1 = F(u,x)
      p_EQUATION -> coord(0) = x;
      p_EQUATION -> residual_fn(u, k1);
      k1 *= h;
      z = u + k1 * W21;

      // k2 = F(z,x+X2*h)
      p_EQUATION -> coord(0) = x + X2 * h;
      p_EQUATION -> residual_fn(z, k2);
      k2 *= h;
      z = u + k1 * W31 + k2 * W32;

      // k3 = F(z,x+X3*h)
      p_EQUATION -> coord(0) = x + X3 * h;
      p_EQUATION -> residual_fn(z, k3);
      k3 *= h;
      z = u + k1 * W41 + k2 * W42 + k3 * W43;

      // k4 = F(z,x+X4*h)
      p_EQUATION -> coord(0) = x + X4 * h;
      p_EQUATION -> residual_fn(z, k4);
      k4 *= h;
      z = u + k1 * W51 + k2 * W52 + k3 * W53 + k4 * W54;

      // k5 = F(z,x+X5*h)
      p_EQUATION -> coord(0) = x + X5 * h;
      p_EQUATION -> residual_fn(z, k5);
      k5 *= h;
      z = u + k1 * W61 + k2 * W62 + k3 * W63 + k4 * W64 + k5 * W65;

      // k6 = F(z,x+X6*h)
      p_EQUATION -> coord(0) = x + X6 * h;
      p_EQUATION -> residual_fn(z, k6);
      k6 *= h;

      e = k1 * U1 + k3 * U3 + k4 * U4 + k5 * U5;
      z = k1 * Z1 + k3 * Z3 + k4 * Z4 + k5 * Z5 + k6 * Z6;

      e -= z;

      // diff = ||e|| -- here use "abs" to deal with Complex systems
      diff = e.inf_norm();

      c = sqrt(sqrt(tol * h / (2 * diff)));
      ok = true;

      // is the first step ok? or does it need reducing?

      if((step == 1) && (c < 1.0)) {
        // step needs reducing so start from initial value again
        ok = false;
        step = 1;
      }

      if(ok) {
        x += h;
        u += z;

        if(step % STORE_EVERY == 0) {
          coords.push_back(x);
          values.push_back(u);
        }

      }

      h *= c;

      if(x + h > X_FINAL) {
        h = (X_FINAL - x);
      }

      if(step >= N) {
        std::string problem;
        problem = "The ODE.shoot45 method reached the maximum \n";
        problem += "number of steps specified by the user. \n";
        throw ExceptionRuntime(problem);
      }

    } while(std::abs(x - X_FINAL) > tol);     // end loop stepping across domain

    // construct the solution mesh stored in this object
    SOLN = OneD_Node_Mesh<_Type>(coords, p_EQUATION -> get_order());
    for(unsigned i = 0; i < coords.size(); ++i) {
      // fill mesh
      SOLN.set_nodes_vars(i, values[i]);
    }

    return u;
  }

  template <typename _Type>
  OneD_Node_Mesh<_Type>& ODE_IVP<_Type>::get_mesh() {
    return SOLN;
  }

  template <typename _Type>
  unsigned& ODE_IVP<_Type>::store_every() {
    return STORE_EVERY;
  }

  // the templated versions we require are:
  template class ODE_IVP<double>
  ;
  template class ODE_IVP<std::complex<double> >
  ;

} // end namespace

