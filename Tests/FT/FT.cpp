// IN PROGRESS
// modified to have a non-oscillating bump

#include <FT.h>
#include <OneD_Node_Mesh.h>
#include <Types.h>
#include <Utility.h>

using namespace CppNoddy;
using namespace std;

int main()
{
  const double xMax(20);
  std::size_t N(256);
  // nodes
  const DenseVector<double> xNodes = Utility::uniform_node_vector(-xMax,xMax,N);
  // the signal
  OneD_Node_Mesh<D_complex> h( xNodes, 2 );
  // two complex values stored across the mesh
  for ( std::size_t i = 0; i < N; ++i ) {
    h(i,0) = 1./(1+(xNodes[i]-0)*(xNodes[i]-0))
      + D_complex(0.,1.)/(1+(xNodes[i]-4)*(xNodes[i]-4));
    h(i,1) = sin(xNodes[i]);
  }
  // original signal
  //h.dump_gnu("./DATA/h.dat");

  // Fourier transformed signal
  OneD_Node_Mesh<D_complex> hft = FT::dft(h);
  //hft.dump_gnu("./DATA/hft.dat");

  // shifted Fourier transformed signal
  OneD_Node_Mesh<D_complex> hft_shifted = FT::ft_shift(hft);
  //hft_shifted.dump_gnu("./DATA/hft_shifted.dat");
  
  // invert the DFT using the unshifted spectrum
  OneD_Node_Mesh<D_complex> hReconstructed = FT::idft(hft,-xMax);
  //hReconstructed.dump_gnu("./DATA/h_reconstruct.dat");

  // subtract the initial and reconstructed data
  for ( std::size_t i = 0; i < N; ++i ) {
    h(i,0) -= hReconstructed(i,0);
    h(i,1) -= hReconstructed(i,1);
  }

  if ( max(h.max(0),h.max(1)) < 1.e-12 ) {
    cout << "\033[1;32;48m  * PASSED \033[0m\n";
    return 0;
  }
  cout << "\033[1;31;48m  * FAILED \033[0m\n";
  cout << "    Final error = " << h.max(0) << "\n";
  return 1;

}
