#!/bin/sh

# Meson build process:
# meson --buildtype=debugoptimized build
meson --buildtype=plain --prefix=/home/hewitt/CURRENT/Projects/CppNoddy/build build
meson configure -Dslepc=true -Dpetscz=true build
cd build
ninja reconfigure
ninja
ninja test
ninja install

#See:
#
#http://hewitt.ddns.net/Dev/CppNoddy/doc/html/index.html

