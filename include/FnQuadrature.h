/// \file FnQuadrature.h
/// A specification for quadrature classes. This
/// class is only useful for explicitly known functions
/// and is instantiated with a function pointer.

#ifndef FNQUADRATURE_H
#define FNQUADRATURE_H

#include <Uncopyable.h>

namespace CppNoddy {

  /// A quadrature class that takes a function pointer.
  class FnQuadrature : private Uncopyable {

   public:

    /// The function pointer associated with this instance.
    typedef void (*fn_ptr)(const double&, double&);

    /// Constructor.
    /// \param ptr_to_fn the function that defines the integrand.
    /// \param x1 left hand boundary of the domain.
    /// \param x2 right hand boundary of the domain.
    /// \param num_of_regions initial number of sub-regions to divide the domain into.
    FnQuadrature(fn_ptr ptr_to_fn,
                 const double& x1, const double& x2,
                 const unsigned& num_of_regions);

    /// Constructor.
    /// \param ptr_to_fn the function that defines the integrand.
    /// \param x1 left hand boundary of the domain.
    /// \param x2 right hand boundary of the domain.
    /// \param nodes A vector of nodal positions.
    FnQuadrature(fn_ptr ptr_to_fn,
                 const double& x1, const double& x2,
                 const DenseVector<double>& nodes);

    /// A set method to define a UNIFORM number of sub intervals.
    /// \param n Number of sub intervals
    void set_subintervals(const unsigned& n);

    /// n-point Gauss rule inefficiently written!
    /// \param n The number of points.
    /// \return Approximation to the integral
    double Gauss(const int& n);

    /// Evaluate the integral by applying an
    /// n-point Gauss rule on each of N sub-intervals.
    /// <b> This is inefficient in terms of class instantiation
    /// for each sub-interval, but not currently an issue. </b>
    /// \param n  The order of the Gauss rule.
    /// \return Approximation to the integral.
    double sub_Gauss(const int& n);

    /// Quick trapezium summation again for sanity checking.
    /// Should be essentially equivalent to the 1-point
    /// Gauss rule.
    /// \return An approximation to the integral.
    double trapezium();

   private:
    double A, B;
    fn_ptr p_FN;
    DenseVector<double> NODES;

  };


}

#endif // FNQUADRATURE_H

