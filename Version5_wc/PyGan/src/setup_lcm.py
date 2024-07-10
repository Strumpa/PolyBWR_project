#
# python3 setup_lcm.py install --home=.
#
from sys import version_info
if version_info[0] == 3 and version_info[1] >= 12:
  from setuptools import setup, Extension
elif version_info[0] > 3:
  from setuptools import setup, Extension
else:
  from distutils.core import setup, Extension
import sysconfig
import numpy

def main():
  import os
  incdir = os.path.join(os.path.split(os.path.abspath(numpy.__file__))[0], "core/include")
  mach = os.path.basename(os.getcwd())
  Compiler = os.environ.get("COMPILER", None) # Compiler selection
  if Compiler == "NVTOOLS":
    libdir="../../lib/"+mach+"_nvidia"
    libNv=os.environ.get("NVTOOLS", None)+"/../lib"
    extralink=["-lnvc","-lnvcpumath"]
  elif Compiler == "INTELTOOLS":
    libdir="../../lib/"+mach+"_intel"
    libNv=" "
    extralink=[ ]
  else:
    libdir="../../lib/"+mach
    libNv=" "
    extralink=[ ]
  setup(name="Lcm",
          version="5.0",
          description="Python interface for the lcm C library API",
          author="Alain Hebert",
          author_email="alain.hebert@polymtl.ca",
          license="LGPL",
          ext_modules=[Extension("lcm", ["lcmmodule.c"],
                      extra_link_args = extralink,
                      include_dirs=["../../../Ganlib/src",incdir],
                      library_dirs=[libdir,libNv],
                      runtime_library_dirs=[libNv],
                      libraries=["Ganlib"] ) ])

if __name__ == "__main__":
    main()
