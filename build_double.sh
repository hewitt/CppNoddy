#!/bin/sh

# Meson build process:
# meson --buildtype=debugoptimized build
# meson --buildtype=plain --prefix=/home/hewitt/CURRENT/Projects/CppNoddy/build build

meson build_double --prefix=/home/hewitt/local
meson configure -Dslepc=true -Dpetscd=true build_double
#meson configure build

cd build_double
ninja reconfigure
ninja
#ninja test
#ninja install

#DESTDIR=/home/hewitt/CURRENT/Projects/CppNoddy/build ninja install

#See:
#
#http://hewitt.ddns.net/Dev/CppNoddy/doc/html/index.html

