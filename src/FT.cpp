/// \file FT.cpp
/// An implementation for a collection of Fourier methods

#include <DenseVector.h>
#include <TwoD_Node_Mesh.h>
#include <Types.h>
#include <OneD_Node_Mesh.h>
#include <TwoD_Node_Mesh.h>
#include <FT.h>


namespace CppNoddy {
  namespace FT {

    OneD_Node_Mesh<D_complex> shift( const OneD_Node_Mesh<D_complex>& ft ) {
      OneD_Node_Mesh<D_complex> ft_shifted( ft );
      std::size_t N = ft.get_nnodes();
      for ( std::size_t i = 0; i < N/2; ++i ) {
        for ( std::size_t var = 0; var < ft.get_nvars(); ++var ) {
          // first half of shifted output
          ft_shifted( i, var ) = ft( N/2+i, var );
          // second half of shifted output
          ft_shifted( N/2 + i, var ) = ft( i, var );
        }
        ft_shifted.coord( i ) = -ft.coord( N/2 - i );
        ft_shifted.coord( N/2 + i ) = ft.coord( i );
      }
      return ft_shifted;
    }    

    OneD_Node_Mesh<D_complex> ishift( const OneD_Node_Mesh<D_complex>& ft ) {
      OneD_Node_Mesh<D_complex> ft_ishifted( ft );
      double df = ft.coord(1)-ft.coord(0);
      std::size_t N = ft.get_nnodes();
      for ( std::size_t i = 0; i < N/2; ++i ) {
        for ( std::size_t var = 0; var < ft.get_nvars(); ++var ) {
          // first half of shifted output
          ft_ishifted( i, var ) = ft( N/2 + i, var );
          // second half of shifted output
          ft_ishifted( N/2 + i, var ) = ft( i, var );
        }
        ft_ishifted.coord( i ) = ft.coord( N/2 + i );
        ft_ishifted.coord( N/2 + i ) = ft.coord(N-1) + i*df;
      }
      return ft_ishifted;
    }        
    
    OneD_Node_Mesh<D_complex> dft_with_shift( const OneD_Node_Mesh<D_complex>& f ) {
      return shift(dft(f));
    }

    TwoD_Node_Mesh<D_complex> dft_with_shift( const TwoD_Node_Mesh<D_complex>&f ) {
      // assume that we idft each ft(.,node)
      unsigned ny = f.get_nnodes().second;
      unsigned nx = f.get_nnodes().first;      
      TwoD_Node_Mesh<D_complex> FTsoln( f.xnodes(), f.ynodes(), f.get_nvars());
      for ( unsigned jy = 0; jy < ny; ++jy ) {
        OneD_Node_Mesh<D_complex> mesh = f.get_xsection_at_ynode(jy);
        OneD_Node_Mesh<D_complex> ftmesh = FT::dft_with_shift(mesh);
        for ( unsigned jx = 0; jx < nx; ++jx ) {
          FTsoln.xcoord(jx)=ftmesh.coord(jx);
          FTsoln.set_nodes_vars( jx, jy, ftmesh.get_nodes_vars(jx) );
        }
      }
      return FTsoln;
    }
    
    OneD_Node_Mesh<D_complex> idft_with_ishift( const OneD_Node_Mesh<D_complex>& ft, double origin ) {
      return idft(ishift(ft),origin);
    }

    TwoD_Node_Mesh<D_complex> idft_with_ishift( const TwoD_Node_Mesh<D_complex>&ft, double origin ) {
      // assume that we idft each ft(.,node)
      unsigned ny = ft.get_nnodes().second;
      unsigned nx = ft.get_nnodes().first;      
      TwoD_Node_Mesh<D_complex> soln( ft.xnodes(), ft.ynodes(), ft.get_nvars());
      for ( unsigned jy = 0; jy < ny; ++jy ) {
        OneD_Node_Mesh<D_complex> ftmesh = ft.get_xsection_at_ynode(jy);
        OneD_Node_Mesh<D_complex> mesh = FT::idft_with_ishift(ftmesh,origin);
        for ( unsigned jx = 0; jx < nx; ++jx ) {
          soln.xcoord(jx)=mesh.coord(jx);
          soln.set_nodes_vars( jx, jy, mesh.get_nodes_vars(jx) );
        }
      }
      return soln;
    }
    
    OneD_Node_Mesh<D_complex> dft( const OneD_Node_Mesh<D_complex>& f ) {
      const D_complex I(0.0,1.0); // imaginary unit
      // number of samples in the signal
      std::size_t N = f.get_nnodes();
      // time step in the signal (has to be constant)
      double dt = f.coord(1)-f.coord(0);
      // sample frequency (rad) is related to inverse of "time" step
      double fs = 2*M_PI/dt;
      // frequency step between frequency nodes
      // - largest frequency is (N-1)*df=(N-1)*fs/N
      double df = fs/N;      
      // container to return the Fourier spectrum
      OneD_Node_Mesh<D_complex> ft( f.nodes(), f.get_nvars() );
      //
      const double Wn = 2*M_PI/N;
      // slow FT
      for ( std::size_t k = 0; k < N; ++k ) {
        // set the frequency
        ft.coord(k) = k*df;  // 0,df,...,(N-1)*df
        //
        //
        // t = n*dt
        // omega = k*df
        // df = (2*pi/dt)/N 
        //
        // variables are all initialised to zero
        DenseVector<D_complex> vars( f.get_nvars(), 0.0 );
        for ( std::size_t n = 0; n < N; ++n ) {
            vars += f.get_nodes_vars(n)*exp(-I*Wn*double(k*n));
        }
        // multiply by "dt" to make this correlate with analytical
        // results from the FT integral 
        ft.set_nodes_vars(k,vars*dt);
      }
      return ft;
    } 

    
    // OneD_Node_Mesh<D_complex> idft( const OneD_Node_Mesh<D_complex>& ft, double origin ) {
    //   const D_complex I(0.0,1.0); // imaginary unit
    //   // number of samples in the signal
    //   std::size_t N = ft.get_nnodes();
    //   // time step in the signal (has to be constant)
    //   double df = ft.coord(1)-ft.coord(0);
    //   // sample frequency
    //   double fs = N*df;
    //   // time step in output signal
    //   double dt = 2*M_PI/fs;
    //   // containter to return the "time series" in
    //   OneD_Node_Mesh<D_complex> f( ft );
    //   for ( std::size_t n = 0; n < N; ++n ) {
    //     f.coord(n) = origin + n*dt;
    //     for ( std::size_t var = 0; var < f.get_nvars(); ++var ) {
    //       f(n,var) = 0.0;
    //       for (  std::size_t k = 0; k < N; ++k ) {
    //         double angle = 2*M_PI*k*n/N;
    //         f(n,var) += ft(k,var)*exp(I*angle);
    //       }
    //       // the dt factor is for consistency with continuous FT integral
    //       f(n,var) /= (N*dt);
    //     }
    //   }
    //   return f;
    // }

    OneD_Node_Mesh<D_complex> idft( const OneD_Node_Mesh<D_complex>& ft, double origin ) {
      const D_complex I(0.0,1.0); // imaginary unit
      // number of "frequencies" in the data
      std::size_t N = ft.get_nnodes();
      // (uniform) frequency step in the data
      double df = ft.coord(1)-ft.coord(0);
      // sample frequency
      double fs = N*df;
      // time step in output signal
      double dt = 2*M_PI/fs;
      // NOTE: dt*df = 2*pi/N
      //
      // containter to return the "time series" in
      OneD_Node_Mesh<D_complex> f( ft );
      //
      const double Wn = 2*M_PI/N;
      // number of variables at each time/frequency
      const std::size_t nVars( ft.get_nvars() );
      // slow inverse FT
      for ( std::size_t n = 0; n < N; ++n ) {
        // set the coordinate/time using the offset passed to this method
        f.coord(n) = origin + n*dt;
        // invert for all variables in the mesh at the same time
        DenseVector<D_complex> vars( nVars, 0.0 );
        for (  std::size_t k = 0; k < N; ++k ) {
          vars += ft.get_nodes_vars(k)*exp(+I*Wn*double(k*n));
        }
        // the dt factor is for consistency with continuous FT integral
        f.set_nodes_vars(n,vars/(N*dt));
      }
      return f;
    }


    

  } // end FT namespace

} // end CppNoddy namespace
