/// \file TwoD_MappedNode_Mesh.h
/// A specification for a two dimensional mesh object on a mapped mesh.
///
#ifndef TWOD_MAPPEDNODE_MESH_H
#define TWOD_MAPPEDNODE_MESH_H

#include <Types.h>
#include <Utility.h>
#include <TwoD_Node_Mesh.h>

namespace CppNoddy
{

  /// A two dimensional mesh utility object.
  template <typename _Type>
  class TwoD_MappedNode_Mesh : public TwoD_Node_Mesh<_Type>
  {

  public:

    TwoD_MappedNode_Mesh()
    {}

    /// ctor
    TwoD_MappedNode_Mesh( const double left, const double right,
        const double bottom, const double top,
        const std::size_t nx, const std::size_t ny, const std::size_t nvars ) :
        TwoD_Node_Mesh<_Type>( left, right, bottom, top, nx, ny, nvars )
    {
      init_mesh_mapping();
    }

    /// dtor
    virtual ~TwoD_MappedNode_Mesh()
    {}

    // mesh mapping defaults
    double CY = 4.0;
    double BY = 30.0;
    double CX = 4.0;
    double BX = 30.0;
    // non-uniform mesh functions -- can be overloaded
    // these take the UNIFORM mesh in the computational domain
    // to the non-uniform mesh in the physical domain
    virtual double FX( const double& comp_x ) const
    {
      return BX+comp_x-BX*exp(-comp_x/CX);
    }
    virtual double FXd( const double& comp_x ) const
    {
      return 1+BX*exp(-comp_x/CX)/CX;
    }
    virtual double FXdd( const double& comp_x ) const
    {
      return -BX*exp(-comp_x/CX)/(CX*CX);
    }
    virtual double FY( const double& comp_y ) const
    {
      return BY+comp_y-BY*exp(-comp_y/CY);
    }
    virtual double FYd( const double& comp_y ) const
    {
      return 1+BY*exp(-comp_y/CY)/CY;
    }
    virtual double FYdd( const double& comp_y ) const
    {
      return -BY*exp(-comp_y/CY)/(CY*CY);
    }


    void init_mesh_mapping()
    {
      double left = TwoD_Node_Mesh<_Type>::X[0];
      double right = TwoD_Node_Mesh<_Type>::X[TwoD_Node_Mesh<_Type>::NX-1];
      double bottom = TwoD_Node_Mesh<_Type>::Y[0];
      double top = TwoD_Node_Mesh<_Type>::Y[TwoD_Node_Mesh<_Type>::NY-1];
      //
      std::cout << "[DEBUG] initializing a MAPPED 2D mesh.\n";
      std::cout << "[DEBUG] physical domain is: (" << left << " to " << right << ") x (" << bottom << " to " << top << ")\n";
      std::cout << "[DEBUG] mapped computational domain is: (" << FX(left) << " to " << FX(right)
        << ") x (" << FY(bottom) << " to " << FY(top) << ")\n";
      std::cout << "[DEBUG] mesh resolution is: " << TwoD_Node_Mesh<_Type>::NX << " x " << TwoD_Node_Mesh<_Type>::NY << "\n";
      // these are the computational coords and are UNIFORMLY spaced
      MAPPED_X_NODES_ = Utility::uniform_node_vector( FX(left), FX(right), TwoD_Node_Mesh<_Type>::NX );
      MAPPED_Y_NODES_ = Utility::uniform_node_vector( FY(bottom), FY(top), TwoD_Node_Mesh<_Type>::NY );
      // The physical nodes are non-uniformly spaced and stored in the base class
      // but we don't know where they are since we only have the mapping from
      // physical to computational (not vice versa). So we make a uniform node meshes
      // then invert the mapping for each node using Newton iteration.
      // TwoD_Node_Mesh<_Type>::X = Utility::uniform_node_vector( left, right, TwoD_Node_Mesh<_Type>::NX );
      // TwoD_Node_Mesh<_Type>::Y = Utility::uniform_node_vector( bottom, top, TwoD_Node_Mesh<_Type>::NY );
      // evaluate what (X,Y) nodes map to the uniformly distributed (MAPPED_X,MAPPED_Y) nodes
      for ( unsigned j = 1; j < TwoD_Node_Mesh<_Type>::NY-1; ++j )
      {
        // for each node in Y
        unsigned kmin(0); double min(99e9);
        for ( unsigned k = 0; k < TwoD_Node_Mesh<_Type>::NY; ++k )
        {
          std::cout << " j = " << j << " k = "<< k << " comparing " << TwoD_Node_Mesh<_Type>::Y[k]
            << " with F = " << FY(TwoD_Node_Mesh<_Type>::Y[k]) << " to " << MAPPED_Y_NODES_[j] << "\n";
          // find the mapped value that is closest to it
          if ( std::abs( FY( TwoD_Node_Mesh<_Type>::Y[k] ) - MAPPED_Y_NODES_[j] ) < min )
          {
            min = std::abs( FY( TwoD_Node_Mesh<_Type>::Y[k] ) - MAPPED_Y_NODES_[j] );
            kmin = k;
            std::cout << " kmin = " << kmin << "\n";
          }
        }
        double y = TwoD_Node_Mesh<_Type>::Y[kmin];
        double delta = 1.e-8;
        double dy = 1.0;
        do
        {
          double newY = FY( y + delta ) - MAPPED_Y_NODES_[j];
          double oldY = FY( y ) - MAPPED_Y_NODES_[j];
          double deriv = (newY-oldY)/delta;
          dy = -oldY/deriv;
          y += dy;
        } while ( fabs(dy) > 1.e-8 );
        TwoD_Node_Mesh<_Type>::Y[ j ] = y;
      }
      // // do the same for the X nodes
      // for ( unsigned j = 1; j < TwoD_Node_Mesh<_Type>::NX-1; ++j )
      // {
      //   // for each node in X_nodes
      //   unsigned kmin(0); double min(99e9);
      //   for ( unsigned k = 0; k < TwoD_Node_Mesh<_Type>::NX; ++k )
      //   {
      //     // find the zeta value that is closest to it
      //     if ( std::abs( FX( TwoD_Node_Mesh<_Type>::X[k] ) - MAPPED_X_NODES_[j] ) < min )
      //     {
      //       min = std::abs( FX( TwoD_Node_Mesh<_Type>::X[k] ) - MAPPED_X_NODES_[j] );
      //       kmin = k;
      //     }
      //   }
      //   double x = TwoD_Node_Mesh<_Type>::X[kmin];
      //   double delta = 1.e-8;
      //   double dx = 1.0;
      //   do
      //   {
      //     double newX = FX( x + delta ) - TwoD_Node_Mesh<_Type>::X[j];
      //     double oldX = FX( x ) - TwoD_Node_Mesh<_Type>::X[j];
      //     double deriv = (newX-oldX)/delta;
      //     dx = -oldX/deriv;
      //     x += dx;
      //   } while ( fabs(dx) > 1.e-8 );
      //   TwoD_Node_Mesh<_Type>::X[ j ] = x;
      // }
    }

    /// Access the nodal position - as a pair.
    /// \param nodex The x nodal position to return
    /// \param nodey The y nodal position to return
    /// \return The spatial position of this node as a pair
    std::pair<double, double> coord( const std::size_t nodex, const std::size_t nodey ) const;

  protected:
    // store UNIFORM x nodal points for the COMPUTATIONAL domain
    DenseVector<double> MAPPED_X_NODES_;
    // store UNIFORM y nodal points for the COMPUTATIONAL domain
    DenseVector<double> MAPPED_Y_NODES_;
  };

}

#endif // TWOD_MAPPEDNODE_MESH_H
