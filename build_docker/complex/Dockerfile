FROM ubuntu:22.04 AS base
LABEL maintainer="Rich Hewitt, richard.hewitt@manchester.ac.uk"

#########################
# base system configure #
#########################
ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=Europe/London
RUN apt-get update
# utils + build-tools + libs
RUN apt-get install -y wget tar \
     build-essential meson ninja-build git python3 gcc gfortran make pkg-config git \
     libmkl-dev && rm -rf /var/lib/apt/lists/*
#######################################################
# we want the MKL version of libs in this base system #
#######################################################
# Probably not needed because we specify MKL explicitly below now in the PETSc ./configure
RUN update-alternatives  --install /usr/lib/x86_64-linux-gnu/blas/libblas.so  libblas.so-x86_64-linux-gnu /usr/lib/x86_64-linux-gnu/libmkl_rt.so 99 &&\
  update-alternatives  --install /usr/lib/x86_64-linux-gnu/blas/libblas.so.3  libblas.so.3-x86_64-linux-gnu /usr/lib/x86_64-linux-gnu/libmkl_rt.so 99 &&\
  update-alternatives  --install /usr/lib/x86_64-linux-gnu/blas/liblapack.so  liblapack.so-x86_64-linux-gnu /usr/lib/x86_64-linux-gnu/libmkl_rt.so 99 &&\
  update-alternatives  --install /usr/lib/x86_64-linux-gnu/blas/libblapack.so.3  liblapack.so.3-x86_64-linux-gnu /usr/lib/x86_64-linux-gnu/libmkl_rt.so 99
# setup the user
RUN useradd -ms /bin/bash noddy && echo "noddy:noddy" | chpasswd
######################
# switch to the user #
######################
USER noddy
WORKDIR /home/noddy


#
# PETSC/SLEPC build layer
#


FROM base AS petscslepc
ENV PETSCTYPE complex
#################################
# grab PETSc/SLEPc as the user  #
#################################
#RUN wget --no-verbose https://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-$VPETSC.tar.gz && tar -xzf petsc-$VPETSC.tar.gz
#RUN wget --no-verbose https://slepc.upv.es/download/distrib/slepc-$VSLEPC.tar.gz && tar -xzf slepc-$VSLEPC.tar.gz
RUN git clone -b release https://gitlab.com/petsc/petsc.git petsc
RUN git clone -b release https://gitlab.com/slepc/slepc
# set up PETSc/SLEPc env vars
ENV PETSC_DIR /home/noddy/petsc
ENV SLEPC_DIR /home/noddy/slepc
ENV PETSC_ARCH x86_64-linux-gnu-$PETSCTYPE-opt
#################
# compile PETSc #
#################
WORKDIR /home/noddy/petsc
RUN ./configure --prefix=/home/noddy/local COPTFLAGS="-O3 -ffast-math -march=native" CXXOPTFLAGS="-O3 -ffast-math -march=native" FOPTFLAGS="-O3 -ffast-math -march=native" --download-superlu_dist \
    --download-mpich --download-scalapack --download-mumps --download-cmake --with-blaslapack-lib=/usr/lib/x86_64-linux-gnu/libmkl_rt.so --with-blaslapack-include=/usr/include/mkl --with-mkl_pardiso-dir=/usr/include/mkl \
    --with-scalar-type=$PETSCTYPE && make all && make install
#################
# compile SLEPc #
#################
WORKDIR /home/noddy/slepc
RUN ./configure --prefix=/home/noddy/local && make all && make install


#
# final CppNoddy layer cuts out the build layer's source files
#


FROM base AS cppnoddy
USER noddy
WORKDIR /home/noddy/local
COPY --from=petscslepc /home/noddy/local .
######################
# configure CppNoddy #
######################
WORKDIR /home/noddy
RUN git clone https://github.com/hewitt/CppNoddy.git
WORKDIR /home/noddy/CppNoddy
# we need to set the *_PATH variables to build CppNoddy with PETSc/SLEPc
ENV LD_LIBRARY_PATH /home/noddy/local/lib:$LD_LIBRARY_PATH
ENV PKG_CONFIG_PATH /home/noddy/local/lib/pkgconfig:$PKG_CONFIG_PATH
RUN meson build --prefix=/home/noddy/local
# CHANGE petscd/petscz for real/complex
RUN meson configure -Dslepc=true -Dpetscz=true build
WORKDIR /home/noddy/CppNoddy/build
RUN ninja reconfigure
RUN ninja
ENV MESON_TESTTHREADS 1
#RUN ninja test
RUN ninja install
# now UPDATE *_PATH to add CppNoddy
ENV LD_LIBRARY_PATH /home/noddy/local/lib/x86_64-linux-gnu:$LD_LIBRARY_PATH
ENV PKG_CONFIG_PATH /home/noddy/local/lib/x86_64-linux-gnu/pkgconfig:$PKG_CONFIG_PATH
# default THREADS, MKL takes priority I think
ENV MKL_NUM_THREADS 4
ENV OMP_NUM_THREADS 4
