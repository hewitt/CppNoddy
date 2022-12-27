#!/bin/sh

# Meson build process:
# meson --buildtype=debugoptimized build
# meson --buildtype=plain --prefix=/home/hewitt/CURRENT/Projects/CppNoddy/build build

meson build_complex --prefix=/home/hewitt/CURRENT/Projects/CppNoddy/build_complex
meson configure -Dslepc=true -Dpetscz=true build_complex
#meson configure build

cd build_complex
ninja reconfigure
ninja
#ninja test
#ninja install

#DESTDIR=/home/hewitt/CURRENT/Projects/CppNoddy/build ninja install

#See:
#
#http://hewitt.ddns.net/Dev/CppNoddy/doc/html/index.html

