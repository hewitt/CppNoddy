/// \file FT.cpp
/// \ingroup Tests
/// \ingroup FT
/// A simple 1D DFT test case. Takes a signal, evaluates the
/// (slow) DFT, then inverts the same DFT and subtracts the
/// result from the original signal data.

#include <FT.h>
#include <OneD_Node_Mesh.h>
#include <Types.h>
#include <Utility.h>

using namespace CppNoddy;
using namespace std;

int main()
{
  
  cout << "\n";
  cout << "=== FT: one-dimensional Fourier Transform (slow) ====\n";
  cout << "\n";

  // physical domain size
  const double xMax(20);
  // number of points/wavenumbers/frequencies
  std::size_t N(256);
  // nodes
  const DenseVector<double> xNodes = Utility::uniform_node_vector(-xMax,xMax,N);
  // the physical signal
  OneD_Node_Mesh<D_complex> h( xNodes, 2 );
  // two *complex* values stored across the mesh
  for ( std::size_t i = 0; i < N; ++i ) {
    // sum of 2 "humps" of the form 1/(1+x^2)
    h(i,0) = 1./(1+(xNodes[i]-0)*(xNodes[i]-0))
      + D_complex(0.,1.)/(1+(xNodes[i]-4)*(xNodes[i]-4));
    // a simple sinusoid
    h(i,1) = sin(xNodes[i]);
  }
  // original (physical) signal
  // h.dump_gnu("./DATA/h.dat");

  OneD_Node_Mesh<D_complex> testft = FT::dft_with_shift(h);
  testft.dump_gnu("./DATA/testft.dat");

  OneD_Node_Mesh<D_complex> testf = FT::idft_with_ishift(testft);
  testf.dump_gnu("./DATA/testf.dat");


  
  // Fourier transformed signal
  OneD_Node_Mesh<D_complex> hft = FT::dft(h);
  // hft.dump_gnu("./DATA/hft.dat");

  // shifted Fourier transformed signal
  // this puts the spectrum in a -k_max to +k_max format
  OneD_Node_Mesh<D_complex> hft_shifted = FT::shift(hft);
  // hft_shifted.dump_gnu("./DATA/hft_shifted.dat");
  
  // invert the DFT (must be done using the UNshifted spectrum)
  // the starting (physcial)_ x-coordinate is included as
  // an optional second argument here, otherwise the reconstruction
  // will start at x=0.
  OneD_Node_Mesh<D_complex> hReconstructed = FT::idft(hft,-xMax);
  // hReconstructed.dump_gnu("./DATA/h_reconstruct.dat");

  // subtract the initial and reconstructed data for both data elements
  // which should give zero
  for ( std::size_t i = 0; i < N; ++i ) {
    h(i,0) -= hReconstructed(i,0);
    h(i,1) -= hReconstructed(i,1);
  }

  // check the biggest error in the reconstruction
  if ( max(h.max_abs(0),h.max_abs(1)) < 1.e-12 ) {
    cout << "\033[1;32;48m  * PASSED \033[0m\n";
    return 0;
  }
  cout << "\033[1;31;48m  * FAILED \033[0m\n";
  cout << "    Final |error| = " << h.max_abs(0) << "\n";
  return 1;

}
