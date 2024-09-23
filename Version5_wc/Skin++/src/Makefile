#---------------------------------------------------------------------------
#
#  Makefile for building the Skin++ library and load module
#  Author : A. Hebert (2018-5-10)
#
#---------------------------------------------------------------------------
#
ARCH = $(shell uname -m)
ifeq ($(ARCH),aarch64)
  nbit =
else
  ifeq (,$(filter $(ARCH),i386 i686))
    nbit = -m64
  else
    nbit = -m32
  endif
endif

DIRNAME = $(shell uname -sm | sed 's/[ ]/_/')
OS = $(shell uname -s | cut -d"_" -f1)
opt = -O -g
ifeq ($(openmp),1)
  COMP = -fopenmp
  FOMP = -fopenmp -D_OPENMP
else
  COMP =
  FOMP =
endif

ifeq ($(intel),1)
  fcompiler = ifort
else
  ifeq ($(nvidia),1)
    fcompiler = nvfortran
  else
    fcompiler = gfortran
  endif
endif

clib = -lstdc++
ifeq ($(OS),Darwin)
  ifeq ($(openmp),1)
    C = g++ -std=c++17
  else
    C = clang++ -std=c++17
  endif
  F90 = $(fcompiler)
  FLAGS = -DLinux -DUnix
  CFLAGS = -Wall $(nbit) -fPIC
  LFLAGS = $(nbit)
  clib = -lc++
else
ifeq ($(OS),Linux)
  ifeq ($(nvidia),1)
    C = nvc++ -std=c++17
  else
    C = g++ -std=c++17
  endif
  F90 = $(fcompiler)
  FLAGS = -DLinux -DUnix
  CFLAGS = -Wall $(nbit) -fPIC
  LFLAGS = $(nbit)
else
ifeq ($(OS),CYGWIN)
  C = g++ -std=c++17
  F90 = $(fcompiler)
  FLAGS = -DLinux -DUnix
  CFLAGS = -Wall $(nbit) -fPIC
  LFLAGS = $(nbit)
else
ifeq ($(OS),SunOS)
  MAKE = gmake
  F90 = f90
  C = cc
  FLAGS = -DSunOS -DUnix
  CFLAGS = $(nbit)
  LFLAGS = $(nbit)
else
ifeq ($(OS),AIX)
  opt = -O4
  MAKE = gmake
  DIRNAME = AIX
  F90 = xlf90
  C = xlc
  FLAGS = -DAIX -DUnix
  CFLAGS = -qstrict
  LFLAGS = -qstrict -bmaxdata:0x80000000 -qipa
else
  $(error $(OS) is not a valid OS)
endif
endif
endif
endif
endif
INCLUDE = -I../../Ganlib/src/

ifeq ($(intel),1)
  FFLAGS = -fPIC
  FFLAG77 = -fPIC
  lib = ../lib/$(DIRNAME)_intel
  bin = ../bin/$(DIRNAME)_intel
  libUtl = ../../Utilib/lib/$(DIRNAME)_intel
  libGan = ../../Ganlib/lib/$(DIRNAME)_intel
  libTri = ../../Trivac/lib/$(DIRNAME)_intel
  libDra = ../../Dragon/lib/$(DIRNAME)_intel
  libDon = ../../Donjon/lib/$(DIRNAME)_intel
else
  ifeq ($(nvidia),1)
    lib = ../lib/$(DIRNAME)_nvidia
    bin = ../bin/$(DIRNAME)_nvidia
    libUtl = ../../Utilib/lib/$(DIRNAME)_nvidia
    libGan = ../../Ganlib/lib/$(DIRNAME)_nvidia
    libTri = ../../Trivac/lib/$(DIRNAME)_nvidia
    libDra = ../../Dragon/lib/$(DIRNAME)_nvidia
    libDon = ../../Donjon/lib/$(DIRNAME)_nvidia
  else
    lib = ../lib/$(DIRNAME)
    bin = ../bin/$(DIRNAME)
    libUtl = ../../Utilib/lib/$(DIRNAME)
    libGan = ../../Ganlib/lib/$(DIRNAME)
    libTri = ../../Trivac/lib/$(DIRNAME)
    libDra = ../../Dragon/lib/$(DIRNAME)
    libDon = ../../Donjon/lib/$(DIRNAME)
  endif
endif

ifeq ($(hdf5),1)
  CFLAGS += -DHDF5_LIB -I${HDF5_INC}
  LFLAGS += -L${HDF5_API} -lhdf5
endif

SRCC = $(shell ls *.cxx)
OBJC = $(SRCC:.cxx=.o)
all : sub-make Skin++
ifeq ($(openmp),1)
	@echo 'Skin++: openmp is defined'
endif
ifeq ($(intel),1)
	@echo 'Skin++: intel is defined'
endif
ifeq ($(nvidia),1)
  @echo 'Skin++: nvidia is defined'
endif
ifeq ($(hdf5),1)
	@echo 'Skin++: hdf5 is defined'
endif
sub-make:
	$(MAKE) openmp=$(openmp) intel=$(intel) hdf5=$(hdf5) -C ../../Donjon/src
%.o : %.cxx
	$(C) $(CFLAGS) $(FLAGS) $(opt) $(COMP) $(INCLUDE) -c $< -o $@
$(lib)/:
	mkdir -p $(lib)/
libSkin++.a: $(OBJC) $(lib)/
	ar r $@ $(OBJC)
	cp $@ $(lib)/$@
$(bin)/:
	mkdir -p $(bin)/
Skin++: libSkin++.a Skin++.o $(bin)/ sub-make
	$(F90) $(opt) $(FOMP) Skin++.o $(lib)/libSkin++.a $(libDon)/libDonjon.a $(libDra)/libDragon.a $(libTri)/libTrivac.a \
	$(libDra)/libDragon.a $(libUtl)/libUtilib.a $(libGan)/libGanlib.a $(clib) $(LFLAGS) -o Skin++
	cp $@ $(bin)/$@
clean:
	$(MAKE) -C ../../Donjon/src clean
	/bin/rm -f *.o *.a sub-make Skin++
