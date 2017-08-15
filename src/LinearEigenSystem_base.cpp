/// \file LinearEigenSystem_base.cpp
/// Implementation for the LinearEigenSystem_base class.
/// The specific eigesolvers inherit from here.

#include <vector>
#include <set>

#include <LinearEigenSystem_base.h>
#include <Exceptions.h>
#include <Types.h>

namespace CppNoddy
{

  LinearEigenSystem_base::LinearEigenSystem_base( ) :
      SHIFT( D_complex( 0., 0. ) ),
      CALC_EIGENVECTORS( true )
  {
  }


  LinearEigenSystem_base::~LinearEigenSystem_base()
  {}


  void LinearEigenSystem_base::set_shift( const D_complex& z )
  {
    SHIFT = z;
  }


  std::complex<double> LinearEigenSystem_base::get_shift() const
  {
    return SHIFT;
  }


  void LinearEigenSystem_base::eigensolve()
  {
    std::string problem;
    problem = "The LinearEigenSystem_base::eigensolve method has been called\n";
    problem += "but the method has not been implemented ... this should be \n";
    problem += "implemented in the sub-class.";
    throw ExceptionExternal( problem );
  }


  void LinearEigenSystem_base::set_calc_eigenvectors( bool flag )
  {
    CALC_EIGENVECTORS = flag;
  }


  DenseVector<D_complex> LinearEigenSystem_base::get_tagged_eigenvalues() const
  {
    if ( TAGGED_INDICES.size() == 0 )
    {
      std::string problem;
      problem = "In LinearEigenSystem_base.get_tagged_eigenvalues() : there are\n";
      problem += "no eigenvalues that have been tagged. This set is empty.\n";
      throw ExceptionRuntime( problem );
    }
    // storage for the eigenvalues
    DenseVector<D_complex> evals;
    // loop through the tagged set
    for ( iter p = TAGGED_INDICES.begin(); p != TAGGED_INDICES.end(); ++p )
    {
      // get the index of the relevant eigenvalue from the set
      std::size_t j = *p;
      std::cout << " number " << j << " is a tagged ev.\n";
      // work out the complex eigenvalue associated with this index
      // and add it to the vector
      evals.push_back( ALL_EIGENVALUES[ j ] );
    }
    // return the complex vector of eigenvalues
    return evals;
  }


  DenseMatrix<D_complex> LinearEigenSystem_base::get_tagged_eigenvectors() const
  {
    if ( TAGGED_INDICES.size() == 0 )
    {
      std::string problem;
      problem = "In LinearEigenSystem_base.get_tagged_eigenvectors() : there are\n";
      problem += "no eigenvalues that have been tagged. This set is empty.\n";
      throw ExceptionRuntime( problem );
    }
    // number of degrees of freedom in each eigenvector
    std::size_t N = ALL_EIGENVECTORS[0].size();
    // eigenvector storage : size() eigenvectors each of length N
    DenseMatrix<D_complex> evecs( TAGGED_INDICES.size(), N, 0.0 );
    std::size_t row = 0;
    // loop through the tagged set
    for ( iter p = TAGGED_INDICES.begin(); p != TAGGED_INDICES.end(); ++p )
    {
      // get the index of the relevant eigenvalue from the set
      std::size_t j = *p;
      // put the eigenvector in the matrix
      evecs[ row ] = ALL_EIGENVECTORS[ j ];
      // next row/eigenvector
      ++row;
    }
    return evecs;
  }

  // EIGENVALUE/VECTOR TAGGING


  void LinearEigenSystem_base::tag_eigenvalues_all()
  {
    for ( unsigned i = 0; i < ALL_EIGENVALUES.size(); ++i )
    {
      TAGGED_INDICES.insert( TAGGED_INDICES.end(), i );
    }
    std::cout << "** I have tagged " << ALL_EIGENVALUES.size() << " eigenvalues.\n";
  }


  void LinearEigenSystem_base::tag_eigenvalues_disc( const int &val, const double& radius )
  {
    // loop through all the eigenvalues
    for ( std::size_t i = 0; i < ALL_EIGENVALUES.size(); ++i )
    {
      // if the eigenvalue is in the disc centred at shift then include it
      if ( std::abs( ALL_EIGENVALUES[ i ] - SHIFT ) < radius )
      {
        if ( val > 0 )
        {
          // add it to our set of tagged eigenvalues
          TAGGED_INDICES.insert( TAGGED_INDICES.end(), i );
        }
        else
        {
          // remove it from the set if it exists
          TAGGED_INDICES.erase( i );
        }
      }
    }
  }


  void LinearEigenSystem_base::tag_eigenvalues_right( const int &val )
  {
    // loop through all the eigenvalues
    for ( std::size_t i = 0; i < ALL_EIGENVALUES.size(); ++i )
    {
      // if the eigenvalue is in the disc centred at SHIFT then include it
      if ( ( ALL_EIGENVALUES[ i ] - SHIFT ).real() > 0.0 )
      {
        if ( val > 0 )
        {
          // add it to our set of tagged eigenvalues
          TAGGED_INDICES.insert( TAGGED_INDICES.end(), i );
        }
        else
        {
          // remove it from the set if it exists
          TAGGED_INDICES.erase( i );
        }
      }
    }
  }


  void LinearEigenSystem_base::tag_eigenvalues_left( const int &val )
  {
    // loop through all the eigenvalues
    for ( std::size_t i = 0; i < ALL_EIGENVALUES.size(); ++i )
    {
      // if the eigenvalue is in the disc centred at SHIFT then include it
      if ( ( ALL_EIGENVALUES[ i ] - SHIFT ).real() < 0.0 )
      {
        if ( val > 0 )
        {
          // add it to our set of tagged eigenvalues
          TAGGED_INDICES.insert( TAGGED_INDICES.end(), i );
        }
        else
        {
          // remove it from the set if it exists
          TAGGED_INDICES.erase( i );
        }
      }
    }
  }


  void LinearEigenSystem_base::tag_eigenvalues_upper( const int &val )
  {
    // loop through all the eigenvalues
    for ( std::size_t i = 0; i < ALL_EIGENVALUES.size(); ++i )
    {
      // if the eigenvalue is in the disc centred at SHIFT then include it
      if ( ( ALL_EIGENVALUES[ i ] - SHIFT ).imag() > 0.0 )
      {
        if ( val > 0 )
        {
          // add it to our set of tagged eigenvalues
          TAGGED_INDICES.insert( TAGGED_INDICES.end(), i );
        }
        else
        {
          // remove it from the set if it exists
          TAGGED_INDICES.erase( i );
        }
      }
    }
  }

  void LinearEigenSystem_base::tag_eigenvalues_lower( const int &val )
  {
    // loop through all the eigenvalues
    for ( std::size_t i = 0; i < ALL_EIGENVALUES.size(); ++i )
    {
      // if the eigenvalue is in the disc centred at SHIFT then include it
      if ( ( ALL_EIGENVALUES[ i ] - SHIFT ).imag() < 0.0 )
      {
        if ( val > 0 )
        {
          // add it to our set of tagged eigenvalues
          TAGGED_INDICES.insert( TAGGED_INDICES.end(), i );
        }
        else
        {
          // remove it from the set if it exists
          TAGGED_INDICES.erase( i );
        }
      }
    }
  }


  class LinearEigenSystem_base;

} // end namespace
