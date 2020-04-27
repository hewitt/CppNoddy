/// \file EVP2DOrrSommerfeld_slepcz.cpp
/// \ingroup Tests
/// \ingroup EVP
/// Solves the linear eigenvalue problem for values \f$ c \f$ formulated
/// in the paper of Tatsumi and Yoshimura (J Fluid Mechanics, 1990).
/// The macros defined below determine the symmetry of the sought eigenmode
/// again using the classification introduced by Tatsumi and Yoshimura.
/// The test is for the existence of a near (within 1.e-4) neutral wave
/// for a duct flow of aspect ratio 6.


#include <cassert>

#include <Utility.h>
#include <SparseMatrix.h>
#include <SparseLinearEigenSystem.h>
#include <TwoD_Mapped_Node_Mesh.h>
#include <TwoD_Node_Mesh.h>
#include <SlepcSession.h>
#include <TrackerFile.h>

//#define DENSE
//#define UNIFORM
#define TYPEI

#if defined(TYPEI)
  #define V_EVEN_Y
  #define V_EVEN_Z
#endif
#if defined(TYPEII)
  #define V_EVEN_Y
  #define V_ODD_Z
#endif
#if defined(TYPEIII)
  #define V_ODD_Y
  #define V_EVEN_Z
#endif
#if defined(TYPEIV)
  #define V_ODD_Y
  #define V_ODD_Z
#endif

#if defined(V_ODD_Y)
  #define W_EVEN_Y
#elif defined(V_EVEN_Y)
  #define W_ODD_Y
#endif

#if defined(V_EVEN_Z)
  #define W_ODD_Z
#elif defined(V_ODD_Z)
  #define W_EVEN_Z
#endif

using namespace CppNoddy;
using namespace std;

namespace Example
{
  
  // instanitate a mapped mesh
  template <typename _Type>
  class My_Mesh : public TwoD_Mapped_Node_Mesh<_Type>
  {
  public:
    My_Mesh( double left, double right, double bottom, double top, std::size_t nz, std::size_t ny, std::size_t nv ) : TwoD_Mapped_Node_Mesh<_Type>( left, right, bottom, top, nz, ny, nv )
    {
      // mesh stretching parameters
       BX = 10.0;
       BY = 10.0;
       CX = 0.2;
       CY = 0.8;
       this -> init_mapping();
    }

#if defined(UNIFORM)
    // the computational x coord as a function of physical x
    double FnComp_X( const double& x ) const
    {
      return x;
    }
    double FnComp_Xd( const double& x ) const
    {
      return 1.0;
    }
    double FnComp_Xdd( const double& x ) const
    {
      return 0.0;
    }
    // the computational y coord as a function of physical y
    double FnComp_Y( const double& y ) const
    {
      return y;
    }
    double FnComp_Yd( const double& y ) const
    {
      return 1.0;
    }
    double FnComp_Ydd( const double& y ) const
    {
      return 0.0;
    }
#else
    // the computational x coord as a function of physical x
    double FnComp_X( const double& x ) const
    {
      return BX + x - BX*exp(-(x - this -> m_left)/CX );
    }
    double FnComp_Xd( const double& x ) const
    {
      return 1 + BX*exp(-(x - this -> m_left)/CX )/CX;
    }
    double FnComp_Xdd( const double& x ) const
    {
      return - BX*exp(-(x - this -> m_left)/CX )/(CX*CX);
    }
    // the computational y coord as a function of physical y
    double FnComp_Y( const double& y ) const
    {
      return BY + y - BY*exp(-(y - this -> m_bottom)/CY );
    }
    double FnComp_Yd( const double& y ) const
    {
      return 1 + BY*exp(-(y - this -> m_bottom)/CY )/CY;
    }
    double FnComp_Ydd( const double& y ) const
    {
      return - BY*exp(-(y - this -> m_bottom)/CY )/(CY*CY);
    }
#endif

  private:
    double BX,CX,BY,CY;
  };
  
}


int main( int argc, char *argv[] )
{
  SlepcSession::getInstance(argc,argv);
  
  cout << "\n";
  cout << "=== EVP: Temporal spectra of the 2D Orr-Sommerfeld eqn =\n";
  cout << "===  with a matrix problem assembled by hand using     =\n";
  cout << "===  the formulation of Tatsumi and Yoshimura (1990).  =\n";
  cout << "\n";

  enum {v,w,F,G};

  // duct aspect ratio
  const double A(6.0);
  // streamwise wavenumber and Reynolds number
  const double alpha ( 0.94 );
  double Re ( 8200.0 ); 
  const D_complex eye( 0.0, 1.0 );
  
  // discretise with these many nodal points
  const std::size_t NZ( 251 );
  const std::size_t NY( 251 );
  //
  //
  Example::My_Mesh<double> base( -A, 0.0, -1.0, 0.0, NZ, NY, 1 );
  //
  for ( std::size_t i = 0; i < NZ; ++i )
    {
      for ( std::size_t j = 0; j < NY; ++j )
        {
          double Z = base.coord(i,j).first;
          double Y = base.coord(i,j).second;
          double sum = 1.0 - Y*Y;
          for ( std::size_t n = 0; n < 1000; ++n )
            {
              double C = (2*n+1)*M_PI/2.0;
              double term = cos(C*Y)*(exp(C*(Z-A))+exp(-C*(Z+A)))/(1+exp(-2*C*A));
              sum -= 4*std::pow(2./M_PI,3)*std::pow(-1.0,n)*term/std::pow(2*n+1,3);
            }
          base(i,j,0) = sum;
        }
    }
  base.dump_gnu("./DATA/U.dat");
  

  // we'll solve as FOUR equations
  const std::size_t N( 4*NZ*NY  );

  cout << "There are " << N << " degrees of freedom in this problem.\n";
  
  // spatial step for the uniform computational-domain mesh
  const double DZ = base.get_comp_step_sizes().first;
  const double DY = base.get_comp_step_sizes().second;
  //std::cout << " DZ = " << DZ << " DY = " << DY << "\n";

  // matrices for the EVP, initialised with zeroes
#ifdef DENSE
  DenseMatrix<D_complex> a( N, N, 0.0 );
  DenseMatrix<D_complex> b( N, N, 0.0 );
#else
  SparseMatrix<D_complex> a( N, N );
  SparseMatrix<D_complex> b( N, N );
#endif
  
  std::size_t dof_per_i( 4*NY );
  //
  for ( std::size_t j = 0; j < NY; ++j )
    {
      const std::size_t i(0);
      // left boundary conditions for the duct
      const std::size_t O( i*4*NY );
      std::size_t k( j*4 );
      // coordinate mapping
      double Z( base.coord(i,j).first );
      double KZd( base.FnComp_Xd( Z ) );
      double KZdd( base.FnComp_Xdd( Z ) );
      // v = 0
      a( O+k+0, O+k + v ) = 1.0;
      // w = 0
      a( O+k+1, O+k + w ) = 1.0;
      // definition of F with w=0 for all y and 2nd order 1-sided diff.
      a( O+k+2, O+k + dof_per_i + v ) = -5.0*(eye/(alpha*Re))*KZd*KZd/(DZ*DZ)
        + 4.0*(eye/(alpha*Re))*KZdd/(2.0*DZ); 
      a( O+k+2, O+k + 2*dof_per_i + v ) = 4.0*(eye/(alpha*Re))*KZd*KZd/(DZ*DZ)
        - 1.0*(eye/(alpha*Re))*KZdd/(2.0*DZ); 
      a( O+k+2, O+k + 3*dof_per_i + v ) = -1.0*(eye/(alpha*Re))*KZd*KZd/(DZ*DZ);
      a( O+k+2, O+k + F ) = -1.0;
      // definition of G -- building in w_z = 0 with high order scheme
      a( O+k+3, O+k + dof_per_i + w ) = 10.0*(eye/(alpha*Re))*KZd*KZd/(3*DZ*DZ);
      a( O+k+3, O+k + 2*dof_per_i + w ) = -1.0*(eye/(alpha*Re))*KZd*KZd/(3*DZ*DZ);
      a( O+k+3, O+k + G ) = -1.0;
    }
  //
  for ( std::size_t i = 1; i <= NZ - 2; ++i )
  {
    {
      const std::size_t j(0);
      // bottom boundary conditions for the duct
      std::size_t O( i*4*NY );
      std::size_t k( j*4 );
      // coordinate mapping
      double Y( base.coord(i,j).second );
      double KYd( base.FnComp_Yd( Y ) );
      double KYdd( base.FnComp_Ydd( Y ) );
      // v=0 on y=0
      a( O+k+0, O+k + v ) = 1.0;
      // w=0 on y=0
      a( O+k+1, O+k + w ) = 1.0;
      // definition of F, building in v_y=0 with high order scheme
      a( O+k+2, O+k + 4 + v ) = 10.0*(eye/(alpha*Re))*KYd*KYd/(3*DY*DY);
      a( O+k+2, O+k + 8 + v ) = -1.0*(eye/(alpha*Re))*KYd*KYd/(3*DY*DY);
      a( O+k+2, O+k + F ) = -1.0;
      // definition of G with w=0 for all y and 2nd order 1-sided diff.
      a( O+k+3, O+k + 4 + w ) = -5.0*(eye/(alpha*Re))*KYd*KYd/(DY*DY)
        + 4.0*(eye/(alpha*Re))*KYdd/(2.0*DY);
      a( O+k+3, O+k + 8 + w ) = 4.0*(eye/(alpha*Re))*KYd*KYd/(DY*DY)
        - 1.0*(eye/(alpha*Re))*KYdd/(2.0*DY);
      a( O+k+3, O+k + 12 + w ) = -1.0*(eye/(alpha*Re))*KYd*KYd/(DY*DY);
      a( O+k+3, O+k + G ) = -1.0; 
    }
    for ( std::size_t j = 1; j < NY-1; ++j )
      {
        // coordinate mapping
        double Z( base.coord(i,j).first );
        double KZd( base.FnComp_Xd( Z ) );
        double KZdd( base.FnComp_Xdd( Z ) );
        double Y( base.coord(i,j).second );
        double KYd( base.FnComp_Yd( Y ) );
        double KYdd( base.FnComp_Ydd( Y ) );
        //
        double Uy  = (base(i,j+1,0)-base(i,j-1,0))*KYd/(2*DY);
        double Uz  = (base(i+1,j,0)-base(i-1,j,0))*KZd/(2*DZ);
        double Uyy = (base(i,j+1,0)-2*base(i,j,0)+base(i,j-1,0))*KYd*KYd/(DY*DY)
          + (base(i,j+1,0)-base(i,j-1,0))*KYdd/(2*DY);
        double Uzz = (base(i+1,j,0)-2*base(i,j,0)+base(i-1,j,0))*KZd*KZd/(DZ*DZ)
          + (base(i+1,j,0)-base(i-1,j,0))*KZdd/(2*DZ);
        double Uyz = (base(i+1,j+1,0)-base(i+1,j-1,0)
                      - base(i-1,j+1,0)+base(i-1,j-1,0))*KYd*KZd/(4*DY*DZ);
        // interior equations within the duct
        std::size_t O( i*4*NY );
        std::size_t k( j*4 );
        //
        // L{v}=F
        //
        // v_i,j
        a( O+k+0, O+k + v ) = (eye/(alpha*Re))*( - 2.0*KYd*KYd/(DY*DY)
                                                 - 2.0*KZd*KZd/(DZ*DZ)
                                                 - alpha*alpha ) + base(i,j,0);
        // v_i+1,j
        a( O+k+0, O+k + dof_per_i + v ) = (eye/(alpha*Re))*KZd*KZd/(DZ*DZ)
          + (eye/(alpha*Re))*KZdd/(2.0*DZ);
        // v_i-1,j
        a( O+k+0, O+k - dof_per_i + v ) = (eye/(alpha*Re))*KZd*KZd/(DZ*DZ)
          - (eye/(alpha*Re))*KZdd/(2.0*DZ);
        // v_i,j+1
        a( O+k+0, O+k + 4 + v ) = (eye/(alpha*Re))*KYd*KYd/(DY*DY)
          + (eye/(alpha*Re))*KYdd/(2.0*DY);
        // v_i,j-1
        a( O+k+0, O+k - 4 + v ) = (eye/(alpha*Re))*KYd*KYd/(DY*DY)
          - (eye/(alpha*Re))*KYdd/(2.0*DY);
        // F_i,j
        a( O+k+0, O+k + F ) = -1.0;
        // ev -- wave speed c
        b( O+k+0, O+k + v ) = 1.0;
        //
        // L{w}=G
        //
        // w_i,j
        a( O+k+1, O+k + w ) = (eye/(alpha*Re))*( - 2.0*KYd*KYd/(DY*DY)
                                                 - 2.0*KZd*KZd/(DZ*DZ)
                                                 - alpha*alpha ) + base(i,j,0);
        // w_i+1,j
        a( O+k+1, O+k + dof_per_i + w ) = (eye/(alpha*Re))*KZd*KZd/(DZ*DZ)
          + (eye/(alpha*Re))*KZdd/(2.0*DZ);
        // w_i-1,j
        a( O+k+1, O+k - dof_per_i + w ) = (eye/(alpha*Re))*KZd*KZd/(DZ*DZ)
          - (eye/(alpha*Re))*KZdd/(2.0*DZ);
        // w_i,j+1
        a( O+k+1, O+k + 4 + w ) = (eye/(alpha*Re))*KYd*KYd/(DY*DY)
          + (eye/(alpha*Re))*KYdd/(2.0*DY);
        // w_i,j-1
        a( O+k+1, O+k - 4 + w ) = (eye/(alpha*Re))*KYd*KYd/(DY*DY)
          - (eye/(alpha*Re))*KYdd/(2.0*DY);
        // G_i,j
        a( O+k+1, O+k + G ) = -1.0;
        // ev -- wave speed c
        b( O+k+1, O+k + w ) = 1.0;
        //
        // E(y,z){v} = O(y,z){w}
        //
        // F_i,j
        a( O+k+2, O+k + F ) = -2.0*KYd*KYd/(DY*DY) - alpha*alpha;
        // F_i,j+1
        a( O+k+2, O+k + 4 + F ) = 1.0*KYd*KYd/(DY*DY) + KYdd/(2.0*DY);
        // F_i,j-1
        a( O+k+2, O+k - 4 + F ) = 1.0*KYd*KYd/(DY*DY) - KYdd/(2.0*DY);
        //
        // G_i+1,j+1
        a( O+k+2, O+k + dof_per_i + 4 + G ) = 0.25*KYd*KZd/(DY*DZ);
        // G_i+1,j-1
        a( O+k+2, O+k + dof_per_i - 4 + G ) = -0.25*KYd*KZd/(DY*DZ);
        // G_i-1,j-1
        a( O+k+2, O+k - dof_per_i - 4 + G ) = 0.25*KYd*KZd/(DY*DZ);
        // G_i-1,j+1
        a( O+k+2, O+k - dof_per_i + 4 + G ) = -0.25*KYd*KZd/(DY*DZ);
        //
        // v_i,j
        a( O+k+2, O+k + v ) = -2.0*Uyy;
        // v_i,j+1
        a( O+k+2, O+k + 4 + v ) = -2.0*Uy*KYd/(2.0*DY);
        // v_i,j-1
        a( O+k+2, O+k - 4 + v ) = +2.0*Uy*KYd/(2.0*DY);
        //
        // w_i,j
        a( O+k+2, O+k + w ) = -2.0*Uyz;
        // w_i,j+1
        a( O+k+2, O+k + 4 + w ) = -2.0*Uz*KYd/(2.0*DY);
        // w_i,j-1
        a( O+k+2, O+k - 4 + w ) = +2.0*Uz*KYd/(2.0*DY);
        //
        //
        // E(z,y){w} = O(z,y){v}
        //
        // G_i,j
        a( O+k+3, O+k + G ) = -2.0*KZd*KZd/(DZ*DZ) - alpha*alpha; 
        // G_i+1,j
        a( O+k+3, O+k + dof_per_i + G ) = 1.0*KZd*KZd/(DZ*DZ) + KZdd/(2.0*DZ);
        // G_i-1,j
        a( O+k+3, O+k - dof_per_i + G ) = 1.0*KZd*KZd/(DZ*DZ) - KZdd/(2.0*DZ);
        //
        // F_i+1,j+1
        a( O+k+3, O+k + dof_per_i + 4 + F ) = 0.25*KYd*KZd/(DY*DZ);
        // F_i+1,j-1
        a( O+k+3, O+k + dof_per_i - 4 + F ) = -0.25*KYd*KZd/(DY*DZ);
        // F_i-1,j-1
        a( O+k+3, O+k - dof_per_i - 4 + F ) = 0.25*KYd*KZd/(DY*DZ);
        // F_i-1,j+1
        a( O+k+3, O+k - dof_per_i + 4 + F ) = -0.25*KYd*KZd/(DY*DZ);
        //
        // v_i,j
        a( O+k+3, O+k + v ) = -2.0*Uyz;
        // v_i+1,j
        a( O+k+3, O+k + dof_per_i + v ) = -2.0*Uy*KZd/(2.0*DZ);
        // v_i-1,j
        a( O+k+3, O+k - dof_per_i + v ) = +2.0*Uy*KZd/(2.0*DZ);
        // 
        // w_i,j
        a( O+k+3, O+k + w  ) = -2.0*Uzz;
        // w_i+1,j
        a( O+k+3, O+k + dof_per_i + w ) = -2.0*Uz*KZd/(2.0*DZ);
        // w_i-1,j
        a( O+k+3, O+k - dof_per_i + w ) = +2.0*Uz*KZd/(2.0*DZ);
      }
    {
      // top boundary
      //
      const std::size_t j(NY-1);
      std::size_t O( i*4*NY );
      const std::size_t k( j*4 );
      // coordinate mapping
      double Y( base.coord(i,j).second );
      double KYd( base.FnComp_Yd( Y ) );
      //double KYdd( base.FnComp_Ydd( Y ) );
      //
#ifdef V_EVEN_Y
      // v_y = 0
      a( O+k+0, O+k + v ) = 3.0*KYd/(2.0*DY);
      a( O+k+0, O+k - 4 + v ) = -4.0*KYd/(2.0*DY);
      a( O+k+0, O+k - 8 + v ) = 1.0*KYd/(2.0*DY);
      // v_yyy = 0  => F_y = 0
      a( O+k+2, O+k + F ) = 3.0*KYd/(2.0*DY);
      a( O+k+2, O+k - 4 + F ) = -4.0*KYd/(2.0*DY); 
      a( O+k+2, O+k - 8 + F ) = 1.0*KYd/(2.0*DY);
#elif defined(V_ODD_Y)
      // v=0
      a( O+k+0, O+k + v ) = 1.0;
      // v_yy =0 => F=0
      a( O+k+2, O+k + F ) = 1.0;
#endif
      
#ifdef W_EVEN_Y
      // w_y = 0
      a( O+k+1, O+k + w ) = 3.0*KYd/(2.0*DY);
      a( O+k+1, O+k - 4 + w ) = -4.0*KYd/(2.0*DY);
      a( O+k+1, O+k - 8 + w ) = 1.0*KYd/(2.0*DY);
      // w_yyy = 0  => G_y = 0
      a( O+k+3, O+k + G ) = 3.0*KYd/(2.0*DY);
      a( O+k+3, O+k - 4 + G ) = -4.0*KYd/(2.0*DY); 
      a( O+k+3, O+k - 8 + G ) = 1.0*KYd/(2.0*DY);
#elif defined(W_ODD_Y)
      // w = 0
      a( O+k+1, O+k + w ) = 1.0;
      // G = 0
      a( O+k+3, O+k + G ) = 1.0;
#endif
    }
  }
  {
    for ( std::size_t j = 0; j < NY; ++j )
      {
        // right hand boundary
        const std::size_t i(NZ-1);
        // left boundary conditions for the duct
        const std::size_t O( i*4*NY );
        std::size_t k( j*4 );
        // coordinate mapping
        double Z( base.coord(i,j).first );
        double KZd( base.FnComp_Xd( Z ) );
        //double KZdd( base.FnComp_Xdd( Z ) );
#ifdef V_EVEN_Z
        // v_z = 0
        a( O+k+0, O+k + v ) = 3.0*KZd/(2.0*DZ);
        a( O+k+0, O+k - dof_per_i + v ) = -4.0*KZd/(2.0*DZ);
        a( O+k+0, O+k - 2*dof_per_i + v ) = 1.0*KZd/(2.0*DZ);
        // v_zzz = 0  => F_z = 0
        a( O+k+2, O+k + F ) = 3.0*KZd/(2.0*DZ);
        a( O+k+2, O+k - dof_per_i + F ) = -4.0*KZd/(2.0*DZ); 
        a( O+k+2, O+k - 2*dof_per_i + F ) = 1.0*KZd/(2.0*DZ);
#elif defined(V_ODD_Z)
        // v=0
        a( O+k+0, O+k + v ) = 1.0;
        // v_zz =0 => F=0
        a( O+k+2, O+k + F ) = 1.0;
#endif
        
#ifdef W_EVEN_Z
        // w_z = 0
        a( O+k+1, O+k + w ) = 3.0*KZd/(2.0*DZ);
        a( O+k+1, O+k - dof_per_i + w ) = -4.0*KZd/(2.0*DZ);
        a( O+k+1, O+k - 2*dof_per_i + w ) = 1.0*KZd/(2.0*DZ);
        // w_zzz = 0  => G_z = 0
        a( O+k+3, O+k + G ) = 3.0*KZd/(2.0*DZ);
        a( O+k+3, O+k - dof_per_i + G ) = -4.0*KZd/(2.0*DZ); 
        a( O+k+3, O+k - 2*dof_per_i + G ) = 1.0*KZd/(2.0*DZ);
#elif defined(W_ODD_Z)
        // w = 0
        a( O+k+1, O+k + w ) = 1.0;
        // w_zz = 0  => G = 0
        a( O+k+3, O+k + G ) = 1.0;
#endif
      }
  }
  
  // a vector for the eigenvalues
  DenseVector<D_complex> ev_c;
#ifdef DENSE
  DenseLinearEigenSystem<D_complex> system( &a, &b );
#else
  SparseLinearEigenSystem<D_complex> system( &a, &b );
  system.set_nev(1);
  system.set_region(0.0, 1.0, -1.0, 1.0);
  system.set_target( D_complex( 0.244, 0.0 ) );
  system.set_order("EPS_TARGET_IMAGINARY");
#endif
  
  system.eigensolve();
  
  system.tag_eigenvalues_disc( + 1, 1.0 );
  ev_c = system.get_tagged_eigenvalues();
  TrackerFile spectrum( "./DATA/spectrum_disc_2d.dat" );
  spectrum.push_ptr( &ev_c, "evs" );
  spectrum.update();
  ev_c.dump();
  
  DenseMatrix<D_complex> vec_mtx = system.get_tagged_eigenvectors();  
  Example::My_Mesh<D_complex> eigfn( -A, 0.0, -1.0, 0.0, NZ, NY, 4 );

  for ( std::size_t n = 0; n < ev_c.size(); ++n )
    {
      for ( std::size_t i = 0; i < NZ; ++i )
        {
          for ( std::size_t j = 0; j < NY; ++j )
            {
              eigfn(i,j,v) = vec_mtx(n, i*4*NY + 4*j + v );
              eigfn(i,j,w) = vec_mtx(n, i*4*NY + 4*j + w );
              eigfn(i,j,F) = vec_mtx(n, i*4*NY + 4*j + F );
              eigfn(i,j,G) = vec_mtx(n, i*4*NY + 4*j + G );
            }
        }
      eigfn.dump_gnu("./DATA/eigfn_"+Utility::stringify(n,4)+".dat");
    }
  
  if ( std::abs( ev_c[0].imag() ) < 1.e-4 )
    {
      cout << "\033[1;32;48m  * PASSED \033[0m\n";
      return 0;
    }

  cout << "\033[1;31;48m  * FAILED \033[0m\n";
  cout << "    Final error = " << ev_c[0].imag() << "\n";
  return 1;
  
}
