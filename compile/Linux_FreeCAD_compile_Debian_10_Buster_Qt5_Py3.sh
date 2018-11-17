#!/bin/bash


base_dir=`pwd`

# building FreeCAD
cd $base_dir
cd freecad
mkdir build
cd build


# do not forget the backslash !!!
#cmake \
#-DCMAKE_INSTALL_PREFIX:PATH=/opt/local/FreeCAD-0.17  \
#-DFREECAD_USE_OCC_VARIANT="Official Version"  -DOCC_INCLUDE_DIR=/opt/local/FreeCAD-0.17/include/opencascade  -DOCC_LIBRARY=/opt/local/FreeCAD-0.17/lib/libTKernel.so \
#../freecad

sudo cmake \
-DCMAKE_INSTALL_PREFIX:PATH=/opt/local/FreeCAD-0.17  \
-DBUILD_QT5=ON \
-DFREECAD_BUILD_DEBIAN=ON \
-DLIB_SUFFIX="" \
-DOCC_INCLUDE_DIR="/usr/include/occt" \
-DOCCT_CMAKE_FALLBACK=True \
-DPYTHON_CONFIG_SUFFIX="-python2.7.x86_64-linux-gnu" \
../freecad


# https://forum.freecadweb.org/viewtopic.php?f=4&t=32328
# -DOCCT_CMAKE_FALLBACK=True \

# https://forum.freecadweb.org/viewtopic.php?f=4&p=269835&sid=75e815077edeffb01d7dfc19bc8f7f0c#p269806
# -DPYTHON_CONFIG_SUFFIX="-python2.7.x86_64-linux-gnu" \

#-DPYTHON_CONFIG_SUFFIX="-python2.7.${DEB_HOST_MULTIARCH}" \
#-DCMAKE_CXX_FLAGS="-Wall -DHAVE_SWIG=1 -fpermissive $(shell dpkg-buildflags --get CXXFLAGS) $(shell dpkg-buildflags --get CPPFLAGS)" \
#-DCMAKE_C_FLAGS="-Wall -fpermissive $(shell dpkg-buildflags --get CFLAGS) $(shell dpkg-buildflags --get CPPFLAGS)" \
#-DCMAKE_SHARED_LINKER_FLAGS="$(shell dpkg-buildflags --get LDFLAGS)" \ 


make -j 2
sudo make install
sudo ln -s /opt/local/FreeCAD-0.17/bin/FreeCAD /usr/local/bin/FreeCAD  # make link
sudo ln -s /opt/local/FreeCAD-0.17/bin/FreeCADCmd /usr/local/bin/FreeCADCmd  # make link


