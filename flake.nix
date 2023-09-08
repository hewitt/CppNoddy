{
  description = "A Nix-flake-based development environment with CppNoddy/PETSc/SLEPc.";

  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs/nixos-23.05";
  };

  outputs = { self , nixpkgs ,... }:
    let
      system = "x86_64-linux";
      pkgs = import nixpkgs {
        inherit system;
        overlays = [ (import ./nix/overlay.nix) ];
        config = {
          # Disable if you don't want unfree packages
          allowUnfree = true;
        };                             
      };
      packages = with pkgs; [
        git meson ninja pkg-config
        gnumake cmake
        python39 gfortran12 gcc12
        curl cacert
        vim more
      ];
      nativeBuildInputs = with pkgs.buildPackages; [
        blas lapack
        flex bison
        hwloc mpich libfabric
      ];

      envVars = ''
            export MESON_TESTTHREADS=1
            export MKL_NUM_THREADS=4
            export OMP_NUM_THREADS=4    
            export CPPNODDY_DIR=$PWD
            export PETSC_DIR=$PWD/3rdParty/petsc
            export SLEPC_DIR=$PWD/3rdParty/slepc
            export LD_LIBRARY_PATH=~/.local/noddy/$PETSC_ARCH/lib:$LD_LIBRARY_PATH
            export PKG_CONFIG_PATH=~/.local/noddy/$PETSC_ARCH/lib/pkgconfig:$PKG_CONFIG_PATH
      '';      
    in
      {
        # PETSc real arithmetic
        devShells.x86_64-linux.real = pkgs.mkShell {

          inherit packages;
          inherit nativeBuildInputs;
      
          shellHook = ''
            export PETSC_ARCH=nixos_real
            echo "Dev environment for type: REAL"
          '' + envVars;
        };

        # SLEPc real arithmetic
        devShells.x86_64-linux.complex = pkgs.mkShell {

          inherit packages;
          inherit nativeBuildInputs;
      
          shellHook = ''
            export PETSC_ARCH=nixos_complex
            echo "Dev environment for type: COMPLEX"
          '' + envVars;
        };

      };
}
