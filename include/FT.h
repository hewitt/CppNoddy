/// \file FT.h
/// A spec for a collection of Fourier methods that act on Noddy containers

#ifndef FT_H
#define FT_H

#include <OneD_Node_Mesh.h>
#include <Types.h>

namespace CppNoddy {
  /// Some algorithms associated with CppNoddy containers.
  namespace FT {

    /// A wrapper that converts (x_i,f_i) to complex then calls the
    /// DFT routine below.
    /// \param f The data to be (discrete) Fourier transformed.
    /// \return The Fourier transform in "frequency" space (omega_i,F_i).
    /// The frequencies are omega_i = 0, domega,..., (N-1)*domega.
    /// The SECOND half of the spectrum provides the -ve frequencies (as in MATLAB).    
    OneD_Node_Mesh<D_complex> dft( const OneD_Node_Mesh<double>& f );
    
    /// (Slow) DFT of the real data (x_i,f_i), i = 0, ... N-1; N must be EVEN.
    /// \param f The data to be (discrete) Fourier transformed.
    /// \return The Fourier transform in "frequency" space (omega_i,F_i).
    /// The frequencies are omega_i = 0, domega,..., (N-1)*domega.
    /// The SECOND half of the spectrum provides the -ve frequencies (as in MATLAB).
    OneD_Node_Mesh<D_complex> dft( const OneD_Node_Mesh<D_complex>& f );

    /// (Slow) Inverse DFT of the complex data (omega_i,F_i), i = 0,...,N-1; N must be EVEN.
    /// \param ft The frequency spectrum in the (unshifted) format.
    /// \param origin The position of the first data point returned.
    /// \return The inversion data (origin+x_i, f_i).
    OneD_Node_Mesh<D_complex> idft( const OneD_Node_Mesh<D_complex>& ft, double origin = 0 );

    /// Shift the frequency spectrum obtained from dft to give positive and negative freq.
    /// Note: input should have an even number of nodes.
    /// \param ft The frequency spectrum in the (unshifted) format.
    /// \return The frequency spectrum in the (shifted) format.
    /// -(N/2)domega,...,-domega,0,domega,...,(N/2-1)*domega
    OneD_Node_Mesh<D_complex> ft_shift( const OneD_Node_Mesh<D_complex>& ft );

    /// Invert the shif operation to recover a spectrum in the form that is expected
    /// by the idft routine.
    /// \param ft A previously shifted frequency spectrum.
    /// \return The frequency spectrum in the (unshifted) format.
    /// 0,domega,...,(N/2)*domega,  (N/2+1)*domega,...,(N-1)*domega.
    OneD_Node_Mesh<D_complex> ft_ishift( const OneD_Node_Mesh<D_complex>& ft );
    
  }

} // end namespace

#endif // FT_H
