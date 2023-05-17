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
 * - 2-D Poisson objects (\link #CppNoddy::Poisson_Cartesian Cartesian \endlink and
 *        \link #CppNoddy::Poisson_meridional axisymmetric \endlink cylindrical polars).
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
 * \section get0 Getting and running it
 *
 * You need a machine with a recent C++ compiler and 'git' to clone the latest version.
 * The build system uses 'meson'. The source is hosted on Github
 * and can be obtained using 'git' via
 *
 *       git clone git://github.com/hewitt/CppNoddy.git
 *
 * A minimal install can be achieved via
 *
 *       meson --buildtype=debugoptimized build
 *
 *       cd build
 *
 *       ninja
 *
 * The self-tests can be executed via
 *
 *       ninja tests
 *
 *
 * Optionally you can link to PETSc (for sparse matrix problems) and SLEPc (for sparse matrix eigenvalue problems) via
 *
 *       meson configure -Dpetscz=true -Dslepc=true
 *
 * for complex arithmetic, or
 *
 *       meson configure -Dpetscd=true -Dslepc=true
 *
 * for real arithmetic.
 * If PETSc/SLEPc are not installed system-wide, you will need to set PKG_CONFIG_PATH to point to PETSc.pc and SLEPc.pc via
 *
 *       export PKG_CONFIG_PATH=$PETSC_DIR/$PETSC_ARCH/lib/pkgconfig:$SLEPC/$PETSC_ARCH/lib/pkgconfig
 *
 * assuming a bash shell. The environment variables $PETSC_DIR/$SLEPC_DIR point to the base locations of the respective source trees, whilst $PETSC_ARCH is the name of the build directory (assumed to be the same in both cases). Suggested configure options for PETSc are
 *
 *      ./configure --download-superlu_dist --download-mpich --download-metis --download-parmetis --download-cmake --download-scalapack --download-mumps --with-scalar-type=complex
 *
 * for complex arithmetic or use --with-scalar-type=real for real arithmetic.
 *
 * See the \link Tests \endlink for a starting point.
 *
 * \section best0 Is it fast/accurate?
 *
 * The matrix classes have native solvers that are naive unoptimised Gaussian elimination algorithms. These routines will only be practical (if at all!) for rather `small' matrix/band sizes and do not scale well. If the problem is of even moderate size, then you should link to your local LAPACK/PETSc routines. LAPACK/PETSc/SLEPc libraries are not shipped with CppNoddy, you have to install them separately yourself if they are not available by default.
 *
 * The code is not especially optimised, in fact in many places the code is deliberately  un-optimised for greater transparency; it is not intended for 'heavy duty' problems. The only sanity checks applied are those listed in the test/example codes \link Tests \endlink.
 *
 * \section flames0 It made my machine burst into flames
 *
 * I never said it wouldn't ;-) The code comes with no guarantees.
 *
 * \section add0 I think it needs a CppNoddy::foo<bar> class
 *
 * Feel free to add something. If you're an undergraduate looking for a final-year project or an MSc. student and have  an idea of something to include (or wish to redesign something that I did in a stupid way), then let me know.
 *
 *
 * \htmlonly
<p style="text-align: center;"
<!-- Creative Commons License -->
<a href="http://creativecommons.org/licenses/GPL/2.0/">
<img alt="CC-GNU GPL" border="0" src="http://creativecommons.org/images
/public/cc-GPL-a.png" /></a><br />
This software is licenced under the <a href="http://creativecommons.org/licenses/GPL/2.0/">CC-GNU GPL</a>.
<!-- /Creative Commons License -->
<!--
<rdf:RDF xmlns="http://web.resource.org/cc/"
    xmlns:dc="http://purl.org/dc/elements/1.1/"
    xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
<Work rdf:about="">
   <license rdf:resource="http://creativecommons.org/licenses/GPL/2.0/" />
   <dc:type rdf:resource="http://purl.org/dc/dcmitype/Software" />
</Work>
<License rdf:about="http://creativecommons.org/licenses/GPL/2.0/">
   <permits rdf:resource="http://web.resource.org/cc/Reproduction" />
   <permits rdf:resource="http://web.resource.org/cc/Distribution" />
   <requires rdf:resource="http://web.resource.org/cc/Notice" />
   <permits rdf:resource="http://web.resource.org/cc/DerivativeWorks" />
   <requires rdf:resource="http://web.resource.org/cc/ShareAlike" />
   <requires rdf:resource="http://web.resource.org/cc/SourceCode" />
</License>
</rdf:RDF>
-->
<br/>
&#169; Content created by R.E. Hewitt, 2007.
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
