#!/bin/sh

#Meson build process:
#meson --buildtype=debugoptimized build
meson --buildtype=plain --prefix=/home/hewitt/CURRENT/Projects/CppNoddy/build build
meson configure -Dslepc=true -Dpetscz=true build
cd build
#meson configure -Dprefix=/home/hewitt/CURRENT/Projects/CppNoddy/build
ninja
#ninja test
#DESTDIR=/home/hewitt/CURRENT/Projects/CppNoddy/build ninja install
#DESTDIR=install ninja install
ninja install

#See:
#
#http://hewitt.gotdns.org/Dev/CppNoddy/doc/html/index.html

