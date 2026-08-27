#!/bin/bash 

export DEPENDENCIES_ROOT_DIR=$SALOME_DIST/dependencies

echo "=========================================================="
echo "Salome Distribution ${SALOME_DIST}"
echo "=========================================================="

export ABSOLUTE_APPLI_PATH=$SALOME_DIST/salome
export OMNIORB_USER_PATH="/tmp"
export PYTHONPATH=$SALOME_DIST/kernel/bin/salome:$PYTHONPATH 
export PYTHONPATH=$SALOME_DIST/kernel/lib/salome:$PYTHONPATH
export PYTHONPATH=$SALOME_DIST/kernel/lib/python3.14/site-packages/salome:$PYTHONPATH
export DYLD_LIBRARY_PATH=$DEPENDENCIES_ROOT_DIR/pv-5.13.3/lib

export DYLD_LIBRARY_PATH=$SALOME_DIST/kernel/lib/salome:$DYLD_LIBRARY_PATH
export PATH=$SALOME_DIST/kernel/lib/salome:$PATH
export DYLD_FALLBACK_LIBRARY_PATH=$DYLD_LIBRARY_PATH
export SALOME_VERBOSE=1

export KERNEL_ROOT_DIR=$SALOME_DIST/kernel
export PATH=$KERNEL_ROOT_DIR/bin/salome:$PATH
export DYLD_LIBRARY_PATH=$KERNEL_ROOT_DIR/lib/salome:$DYLD_LIBRARY_PATH
export PATH=$KERNEL_ROOT_DIR/lib/salome:$PATH
export GUI_ROOT_DIR=$SALOME_DIST/gui 
export DYLD_LIBRARY_PATH=$GUI_ROOT_DIR/lib/salome:$DYLD_LIBRARY_PATH
export PATH=$GUI_ROOT_DIR/lib/salome:$PATH
export PATH=$GUI_ROOT_DIR/bin/salome:$PATH

export CommonGeomLib_ROOT_DIR=$SALOME_DIST/CommonGeom 
export DYLD_LIBRARY_PATH=$CommonGeomLib_ROOT_DIR/lib/salome:$DYLD_LIBRARY_PATH
export PATH=$CommonGeomLib_ROOT_DIR/lib/salome:$PATH

export GEOM_ROOT_DIR=$SALOME_DIST/GEOM
export DYLD_LIBRARY_PATH=$GEOM_ROOT_DIR/lib/salome:$DYLD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=$GEOM_ROOT_DIR/lib/python3.14/site-packages:$DYLD_LIBRARY_PATH

export PATH=$GEOM_ROOT_DIR/lib/salome:$PATH
export PATH=$GEOM_ROOT_DIR/lib/python3.14/site-packages:$PATH

export PYTHONPATH=$GUI_ROOT_DIR/lib/salome:$PYTHONPATH
export PYTHONPATH=$GUI_ROOT_DIR/lib/python3.14/site-packages/salome:$PYTHONPATH
export PYTHONPATH=$GEOM_ROOT_DIR/lib/python3.14/site-packages/salome:$PYTHONPATH
export PYTHONPATH=$GEOM_ROOT_DIR/bin/salome:$PYTHONPATH

export LD_LIBRARY_PATH=$DYLD_LIBRARY_PATH

export SalomeAppConfig=$GUI_ROOT_DIR/share/salome/resources/gui
export SALOME_MODULES="GEOM"
export PYTHONPATH=$DEPENDENCIES_ROOT_DIR/PyQt-dist:$PYTHONPATH
export QT_QPA_PLATFORM_PLUGIN_PATH=/opt/homebrew/opt/qt@5/plugins/platforms
export QT_PLUGIN_PATH=/opt/homebrew/opt/qt@5/plugins
