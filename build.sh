#!/usr/bin/env bash

meson setup build --buildtype=plain --optimization=2 --prefix=~/.local/noddy/$PETSC_ARCH

if [ ${PETSC_ARCH} = "nixos_real" ]; then
  meson configure -Dslepc=true -Dpetscd=true build
else
  meson configure -Dslepc=true -Dpetscz=true build
fi
