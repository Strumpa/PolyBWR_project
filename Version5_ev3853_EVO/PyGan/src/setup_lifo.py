#
# python3 setup_lifo.py install --home=.
#
from sys import version_info
if version_info[0] == 3 and version_info[1] >= 12:
  from setuptools import setup, Extension
elif version_info[0] > 3:
  from setuptools import setup, Extension
else:
  from distutils.core import setup, Extension

def main():
  import os
  mach = os.path.basename(os.getcwd())
  Compiler = os.environ.get("COMPILER", None) # Compiler selection
  if Compiler == "NVTOOLS":
    libdir="../../lib/"+mach+"_nvidia"
  elif Compiler == "LLVMTOOLS":
    libdir="../../lib/"+mach+"_llvm"
  elif Compiler == "INTELTOOLS":
    libdir="../../lib/"+mach+"_intel"
  else:
    libdir="../../lib/"+mach
  setup(name="Lifo",
          version="5.0",
          description="Python interface for the lifo C library function",
          author="Alain Hebert",
          author_email="alain.hebert@polymtl.ca",
          license="LGPL",
          ext_modules=[Extension("lifo", ["lifomodule.c"],
                      include_dirs=["../../../Ganlib/src"],
                      library_dirs=[libdir],
                      libraries=["Ganlib"] ) ])

if __name__ == "__main__":
    main()
