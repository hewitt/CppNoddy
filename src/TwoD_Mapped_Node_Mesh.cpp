/// \file TwoD_Mapped_Node_Mesh.cpp
/// Implementation of a two dimensional (mapped) mesh object. Data
/// is stored on a (mapped) nodal mesh.

#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cassert>

#include <Exceptions.h>
#include <DenseVector.h>
#include <DenseMatrix.h>
#include <TwoD_Mapped_Node_Mesh.h>
#include <Utility.h>

namespace CppNoddy
{

  template <typename _Type>
  void TwoD_Mapped_Node_Mesh<_Type>::set_nodes_vars( const std::size_t nodex, const std::size_t nodey, const DenseVector<_Type>& U )
  {
#ifdef PARANOID
    if ( U.size() > NV )
    {
      std::string problem;
      problem = " The TwoD_Mapped_Node_Mesh.set_nodes_vars method is trying to use a \n";
      problem += " vector that has more entries than variables stored in the mesh. \n";
      throw ExceptionRuntime( problem );
    }
#endif
    // assign contents of U to the member data
    std::size_t offset( ( nodex * NY + nodey ) * NV );
    for ( std::size_t var = 0; var < NV; ++var )
    {
      VARS[ offset++ ] = U[ var ];
    }
  }

  template <typename _Type>
  DenseVector<_Type> TwoD_Mapped_Node_Mesh<_Type>::get_nodes_vars( const std::size_t nodex, const std::size_t nodey ) const
  {
#ifdef PARANOID
    if ( nodex > NX - 1 || nodey > NY - 1 )
    {
      std::string problem;
      problem = " The TwoD_Mapped_Node_Mesh.get_nodes_vars method is trying to \n";
      problem += " access a nodal point that is not in the mesh. \n";
      throw ExceptionRange( problem, NX, nodex, NY, nodey );
    }
#endif
    // construct a vector with NV elements starting from a pointer
    DenseVector<_Type> nodes_vars( NV, &VARS[ ( nodex * NY + nodey ) * NV ] );
    return nodes_vars;
  }

  template <typename _Type>
  void TwoD_Mapped_Node_Mesh<_Type>::init_mapping()
  {
    std::cout << "[DEBUG] Physical domain is [" << LEFT << "," << RIGHT
	      << "] x [" << BOTTOM << "," << TOP << "]\n";
    std::cout << "[DEBUG] Computational domain is ["
	      << FnComp_X(LEFT) << "," << FnComp_X(RIGHT)
	      << "] x [" << FnComp_Y(BOTTOM) << "," << FnComp_Y(TOP) << "]\n";
    {
      // a uniform mesh in the computational coordinates                          
      double comp_left( FnComp_X(LEFT) );
      double comp_right( FnComp_X(RIGHT) );
      const double comp_delta_x = ( comp_right - comp_left ) / ( NX - 1 );
      for ( std::size_t i = 0; i < NX; ++i )
	{
	  COMP_X[i] = comp_left + comp_delta_x * i;
	}
    }
    {
      // a uniform mesh in the computational coordinates                          
      double comp_bottom( FnComp_Y(BOTTOM) );
      double comp_top( FnComp_Y(TOP) );
      const double comp_delta_y = ( comp_top - comp_bottom ) / ( NY - 1 );
      for ( std::size_t j = 0; j < NY; ++j )
	{
	  COMP_Y[j] = comp_bottom + comp_delta_y * j;
	}
    }

    // temporary fill to span the domain -- corrected below
    X = Utility::uniform_node_vector(LEFT,RIGHT,NX);
    Y = Utility::uniform_node_vector(BOTTOM,TOP,NY);
    
    // we now need the corresponding Y coordinates in the physical domain           
    {
      for ( unsigned j = 1; j < NY; ++j )
	{
	  // for each node in COMP_Y                                              
	  unsigned kmin(0); double min(99e9);
	  for ( unsigned k = 0; k < NY; ++k )
	    {
	      // find the y value that is closest to it
	      if ( std::abs( FnComp_Y( Y[k] ) - COMP_Y[j] ) < min )
		{
		  min = std::abs( FnComp_Y( Y[k] ) - COMP_Y[j] );
		  kmin = k;
		}
	    }
	  double y = Y[kmin];
	  double delta = 1.e-8;
	  double correction = 1.0;
	  do
	    {
	      double newY = FnComp_Y( y + delta ) - COMP_Y[j];
	      double oldY = FnComp_Y( y ) - COMP_Y[j];
	      double deriv = (newY-oldY)/delta;
	      correction = -oldY/deriv;
	      y += correction;
	    } while ( fabs(correction) > 1.e-8 );
	  Y[ j ] = y;
	}
    }
    // we now need the corresponding X coordinates in the physical domain           
    {
      for ( unsigned i = 1; i < NX; ++i )
	{
	  // for each node in COMP_Y                                              
	  unsigned kmin(0); double min(99e9);
	  for ( unsigned k = 0; k < NY; ++k )
	    {
	      // find the y value that is closest to it
	      if ( std::abs( FnComp_X( X[k] ) - COMP_X[i] ) < min )
		{
		  min = std::abs( FnComp_X( X[k] ) - COMP_X[i] );
		  kmin = k;
		}
	    }
	  double x = X[kmin];
	  double delta = 1.e-8;
	  double correction = 1.0;
	  do
	    {
	      double newX = FnComp_X( x + delta ) - COMP_X[i];
	      double oldX = FnComp_X( x ) - COMP_X[i];
	      double deriv = (newX-oldX)/delta;
	      correction = -oldX/deriv;
	      x += correction;
	    } while ( fabs(correction) > 1.e-8 );
	  X[ i ] = x;
	}
    }
  }

  template <typename _Type>
  std::pair< double, double> TwoD_Mapped_Node_Mesh<_Type>::get_comp_step_sizes() const
  {
    std::pair<double,double> steps;
    steps.first = COMP_X[1]-COMP_X[0];
    steps.second = COMP_Y[1]-COMP_Y[0];
    return steps;
  }
  
  template <typename _Type>
  std::pair< std::size_t, std::size_t > TwoD_Mapped_Node_Mesh<_Type>::get_nnodes() const
  {
    std::pair< std::size_t, std::size_t > nodes;
    nodes.first = NX;
    nodes.second = NY;
    return nodes;
  }

  template <typename _Type>
  std::size_t TwoD_Mapped_Node_Mesh<_Type>::get_nvars() const
  {
    return NV;
  }

  template <typename _Type>
  DenseVector<double>& TwoD_Mapped_Node_Mesh<_Type>::xnodes()
  {
    return X;
  }

  template <typename _Type>
  DenseVector<double>& TwoD_Mapped_Node_Mesh<_Type>::ynodes()
  {
    return Y;
  }

  template<>
  void TwoD_Mapped_Node_Mesh<double>::normalise( const std::size_t& var )
  {
    double maxval( max(var) );
    VARS.scale( 1./maxval );
  }

  template<>
  void TwoD_Mapped_Node_Mesh<D_complex>::normalise( const std::size_t& var )
  {
    //std::cout << "[DEBUG] asked to normalise a complex mesh\n";
    unsigned max_nx( 0 );
    unsigned max_ny( 0 );
    double max( 0.0 );
    // step through the nodes
    for ( unsigned nodex = 0; nodex < X.size(); ++nodex )
    {
      for ( unsigned nodey = 0; nodey < Y.size(); ++nodey )
      {
        if ( std::abs( VARS[ ( nodex * NY + nodey ) * NV + var ] ) > max )
        {
          max = std::abs( VARS[ ( nodex * NY + nodey ) * NV + var ] );
          max_nx = nodex;
          max_ny = nodey;
        }
      }
    }
    D_complex factor( VARS[ ( max_nx * NY + max_ny ) * NV + var ] );
    //std::cout << "[DEBUG] MAX |variable| had complex value of " << factor << "\n";
    VARS.scale( 1./factor );
  }

  template<>
  void TwoD_Mapped_Node_Mesh<double>::dump_gnu( std::string filename ) const
  {
    std::ofstream dump;
    dump.open( filename.c_str() );
    dump.precision( 15 );
    dump.setf( std::ios::showpoint );
    dump.setf( std::ios::showpos );
    dump.setf( std::ios::scientific );
    //
    for ( std::size_t i = 0; i < NX; ++i )
    {
      for ( std::size_t j = 0; j < NY; ++j )
      {
        dump << X[ i ] << " " << Y[ j ] << " ";
        for ( std::size_t var = 0; var < NV; ++var )
        {
          dump << VARS[ ( i * NY + j ) * NV + var ] << " ";
        }
        dump << "\n";
      }
      dump << "\n";
    }
    dump.close();
  }

 template <>
 void TwoD_Mapped_Node_Mesh<D_complex>::dump_gnu( std::string filename ) const
  {
    std::ofstream dump;
    dump.open( filename.c_str() );
    dump.precision( 15 );
    dump.setf( std::ios::showpoint );
    dump.setf( std::ios::showpos );
    dump.setf( std::ios::scientific );
    //
    for ( std::size_t i = 0; i < NX; ++i )
    {
      for ( std::size_t j = 0; j < NY; ++j )
      {
        dump << X[ i ] << " " << Y[ j ] << " ";
        for ( std::size_t var = 0; var < NV; ++var )
        {
          dump << real(VARS[ ( i * NY + j ) * NV + var ]) << " ";
          dump << imag(VARS[ ( i * NY + j ) * NV + var ]) << " ";
        }
        dump << "\n";
      }
      dump << "\n";
    }
    dump.close();
  }

  //the templated versions we require are:
  template class TwoD_Mapped_Node_Mesh<double>
  ;
  template class TwoD_Mapped_Node_Mesh<std::complex<double> >
  ;

}
