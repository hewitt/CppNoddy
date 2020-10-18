/// \file FT.cpp
/// An implementation for a collection of Fourier methods

#include <Types.h>
#include <OneD_Node_Mesh.h>
#include <FT.h>
#include <bits/c++config.h>

namespace CppNoddy {
  namespace FT {

    OneD_Node_Mesh<D_complex> dft( const OneD_Node_Mesh<double>& f ) {
      OneD_Node_Mesh<D_complex> fc( f.nodes(), f.get_nvars() );
      fc.set_vars_from_vector( f.vars_as_vector() );
      return dft(fc);
    }
    
    OneD_Node_Mesh<D_complex> dft( const OneD_Node_Mesh<D_complex>& f ) {
      const D_complex I(0.0,1.0); // imaginary unit
      // number of samples in the signal
      std::size_t N = f.get_nnodes();
      // time step in the signal (has to be constant)
      double dt = f.coord(1)-f.coord(0);
      // sample frequency is inverse of "time" step
      double fs = 1/dt;
      // frequency step between frequency nodes
      double df = fs/N;
      // container to return the Fourier spectrum
      OneD_Node_Mesh<D_complex> ft( f.nodes(), f.get_nvars() );
      // slow FT
      for ( std::size_t k = 0; k < N; ++k ) {
        // set the frequency
        ft.coord(k) = k*df;  // 0,df,...,(N-1)*df
        for ( std::size_t var = 0; var < ft.get_nvars(); ++var ) {
          // initiate the entry
          ft(k,var) = 0.0;
          // compute the FT
          //
          // t = n*dt
          // omega = k*df
          // df = (2*pi/dt)/N 
          //
          for ( std::size_t n = 0; n < N; ++n ) {
            const double Wn = 2*M_PI/N;
            ft(k,var) += f(n,var)*exp(-I*Wn*double(k*n));
          }
        }
      }
      return ft;
    } 

    OneD_Node_Mesh<D_complex> ft_shift( const OneD_Node_Mesh<D_complex>& ft ) {
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

    OneD_Node_Mesh<D_complex> ft_ishift( const OneD_Node_Mesh<D_complex>& ft ) {
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
    
    
    OneD_Node_Mesh<D_complex> idft( const OneD_Node_Mesh<D_complex>& ft, double origin ) {
      const D_complex I(0.0,1.0); // imaginary unit
      // number of samples in the signal
      std::size_t N = ft.get_nnodes();
      // time step in the signal (has to be constant)
      double df = ft.coord(1)-ft.coord(0);
      // sample frequency
      double fs = N*df;
      // time step in output signal
      double dt = 1.0/fs;
      // containter to return the "time series" in
      OneD_Node_Mesh<D_complex> f( ft );
      for ( std::size_t n = 0; n < N; ++n ) {
        f.coord(n) = origin + n*dt;
        for ( std::size_t var = 0; var < f.get_nvars(); ++var ) {
          f(n,var) = 0.0;
          for (  std::size_t k = 0; k < N; ++k ) {
            double angle = 2*M_PI*k*n/N;
            f(n,var) += ft(k,var)*exp(I*angle);
          }
          f(n,var) /= N;
        }
      }
      return f;
    }

    

  } // end FT namespace

} // end CppNoddy namespace
