/// \file Types.h
/// Some standard typedefs. These are kept for legacy reasons. I've tried
/// to remove them from the code, which means things are a little more
/// verbose, but perhaps also a little less obfuscated.

/*! \mainpage Overview
 *
 * \section intro0 Introduction
 * \subsection intro1 What is it?
 *
 *
 * A simple (aka Noddy) collection of object orientated numerical routines written in C++,
 * aimed at undergraduate projects and starting graduate students.
 * In the current version, the test/example cases solve (amongst others):
 *
 *
 * - Two-dimensional parabolic problems
 *      ( eg., the unsteady boundary-layer equations ).
 * - Boundary-value ODE problems
 *      ( eg., the Karman rotating-disk, Blasius boundary-layer and other similarity solutions ).
 * - Arc-length continuation of problems involving limit points
 *      ( eg., the Karman rotating disk equations, Falkner-Skan equation,
 *              the plane Poiseuille flow linear neutral curve. ).
 * - One-dimensional eigenvalue problems ( eg., the (bi-)harmonic equation, Orr-Sommerfeld equation ) solved
 *      both directly and via local methods.
 * - One-dimensional hyperbolic problems
 *      ( eg., Sod's shocktube problem,
 *              linear acoustic waves with reflection in non-uniform medium,
 *              shallow water eqiations ).
 * - Two-dimensional hyperbolic problems
 *      ( eg., compressible Euler problems,
 *              linear acoustic waves,
 *              shallow water eqiations ).
 * - Initial-boundary-value problems
 *      ( eg., the heat diffusion equation in 1D, and the unsteady Karman rotating-disk equations ).
 * - Initial-value problems ( eg., the Lorenz equations ).
 * - Poisson problems in Cartesian and cylindrical geometries.
 *
 *
 *
 * A breakdown of examples into groups is found under the 'Modules' link above.
 * Alternatively, a complete list of examples can be found at this link \link Tests \endlink
 *
 * The library provides:

 * - Both \link #CppNoddy::DenseVector dense vector \endlink and \link #CppNoddy::SparseVector sparse vector \endlink classes (including the usual vector operations).
 * - \link #CppNoddy::DenseMatrix Dense \endlink, \link #CppNoddy::BandedMatrix banded \endlink and \link #CppNoddy::SparseMatrix sparse \endlink matrix classes.
 * - A class for \link #CppNoddy::ODE_IVP ODE IVPs \endlink with 4th-order Runge-Kutta(-Fehlberg) method(s).
 * - A class for \link #CppNoddy::ODE_BVP ODE BVPs \endlink with second-order finite-difference methods and adaptive refinement.
 * - A class for \link #CppNoddy::ODE_EVP ODE EVPs \endlink with second-order finite-difference methods.
 * - A class for \link #CppNoddy::PDE_IBVP IBVPs \endlink with second-order methods in both `space'
 and `time'.
 * - A class for \link #CppNoddy::PDE_double_IBVP Two dimensional parabolic problems \endlink with a second-order box scheme.
 * - Classes for both \link #CppNoddy::OneD_TVDLF_Mesh 1-D \endlink and \link #CppNoddy::TwoD_TVDLF_Mesh 2-D \endlink hyperbolic problems via central scheme algorithms.
 * - \link #CppNoddy::Newton vector \endlink Newton iteration classes.
 * - Arc-length continuation solvers exist for Residual objects and boundary value problems.
 * - An ability to link to selected BLAS, LAPACK and PETSc routines via a simplified API
 *      (these currently include the real/complex generalised eigenproblem solvers,
 *       dense/banded/sparse LU solvers.).
 *
 * \subsection intro2 What is it for?
 * It exists for two reasons:
 * - It's an introduction/framework for final-year undergraduate project students or graduate students.
 * - Just for fun.
 *
 * \subsection add0 I think it needs a CppNoddy::foo<bar> class
 *
 * Feel free to add something. If you're an undergraduate looking for a final-year project or an MSc. student and have  an idea of something to include (or wish to redesign something that I did in a stupid way), then let me know.
 *
 *
<br/>
<p>
Content created by R.E. Hewitt, 2007. MIT license.
</p>

 * \endhtmlonly
 */

#ifndef TYPES_H
#define TYPES_H

#include <complex>
#include <cmath>
#include <sys/stat.h>

#include <DenseVector.h>
#include <DenseMatrix.h>
#include <BandedMatrix.h>


/// A collection of OO numerical routines aimed at simple
/// (typical) applied problems in continuum mechanics.
namespace CppNoddy {

  /// A complex double precision number using std::complex
  typedef std::complex<double> D_complex;

}

#endif // NTYPES_H
