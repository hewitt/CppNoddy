#
# nix-shell --pure 
#
{ pkgs ? import <nixpkgs> { overlays = [ (import ./nix/overlay.nix) ]; } }: # use MKL instead of default blas/lapack
#{ pkgs ? import <nixpkgs> {} }:  # to skip overlay
pkgs.mkShell {
  # nativeBuildInputs is usually what you want in the dev environment
  nativeBuildInputs = with pkgs.buildPackages; [ meson ninja pkg-config git #basic dev pkgs
                                                 gnumake binutils 
                                                 python39 gfortran13 gcc13  #langs
                                                 blas lapack                #libs -- overlay maps to MKL
                                                 flex bison                 #petsc dependencies
                                                 cmake hwloc mpich libfabric
                                                 curl cacert                #allow download-<> in petsc
                                               ];
  # buildInputs is for dependencies you'd need "at run time"
  # BuildInputs = with pkgs.buildPackages; [ ];

  # set env variables. Default to a REAL version of PETSc
  shellHook = ''
    #export PETSC_ARCH=nixos_complex
    export PETSC_ARCH=nixos_real  
    export MESON_TESTTHREADS=1
    export MKL_NUM_THREADS=4
    export OMP_NUM_THREADS=4    
    export CPPNODDY_DIR=$PWD
    export PETSC_DIR=$PWD/3rdParty/petsc
    export SLEPC_DIR=$PWD/3rdParty/slepc
    export LD_LIBRARY_PATH=~/.local/noddy/$PETSC_ARCH/lib:$LD_LIBRARY_PATH
    export PKG_CONFIG_PATH=~/.local/noddy/$PETSC_ARCH/lib/pkgconfig:$PKG_CONFIG_PATH
  '';

}
