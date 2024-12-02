#---------------------------------------------------------------------------
#
#  Makefile for building the Yacs++ library and load module
#  Author : A. Hebert (2018-5-10)
#
#---------------------------------------------------------------------------
#
nbit = -m32
ARCH = $(shell uname -m)
ifeq (,$(filter $(ARCH),i386 i686))
  nbit = -m64
endif

DIRNAME = $(shell uname -sm | sed 's/[ ]/_/')
OS = $(shell uname -s | cut -d"_" -f1)
opt = -O -g
ifeq ($(openmp),1)
  COMP = -fopenmp
else
  COMP =
endif

ifeq ($(OS),Darwin)
  C = g++ -std=c++11
  FLAGS = -DLinux -DUnix
  CFLAGS = -Wall $(nbit) -fPIC
else
ifeq ($(OS),Linux)
  ifeq ($(nvidia),1)
    C = nvc++ -std=c++11
  else
    C = g++ -std=c++11
  endif
  FLAGS = -DLinux -DUnix
  CFLAGS = -Wall $(nbit) -fPIC
else
ifeq ($(OS),CYGWIN)
  C = g++ -std=c++11
  FLAGS = -DLinux -DUnix
  CFLAGS = -Wall $(nbit) -fPIC
else
ifeq ($(OS),SunOS)
  MAKE = gmake
  C = cc
  FLAGS = -DSunOS -DUnix
  CFLAGS = $(nbit)
else
ifeq ($(OS),AIX)
  opt = -O4
  MAKE = gmake
  DIRNAME = AIX
  C = xlc
  FLAGS = -DAIX -DUnix
  CFLAGS = -qstrict
else
  $(error $(OS) is not a valid OS)
endif
endif
endif
endif
endif
lib = ../lib/$(DIRNAME)
INCLUDE = -I../../Ganlib/src/ -I../../Skin++/src/

ifeq ($(hdf5),1)
  CFLAGS += -L${HDF5_DIR}/lib -lhdf5
endif

SRCC = $(shell ls *.cxx)
OBJC = $(SRCC:.cxx=.o)
all : sub-make libYacs++.a
ifeq ($(openmp),1)
	@echo 'Yacs++: openmp is defined'
endif
ifeq ($(intel),1)
	@echo 'Yacs++: intel is defined'
endif
ifeq ($(nvidia),1)
	@echo 'Yacs++: nvidia is defined'
endif
ifeq ($(hdf5),1)
	@echo 'Yacs++: hdf5 is defined'
endif
sub-make:
	$(MAKE) openmp=$(openmp) intel=$(intel) nvidia=$(nvidia) hdf5=$(hdf5) -C ../../Skin++/src
%.o : %.cxx
	$(C) $(CFLAGS) $(FLAGS) $(opt) $(COMP) $(INCLUDE) -I$(BOOST_ROOT) -I$(SALOME_KERNEL)/include/salome/ -c $< -o $@
$(lib)/:
	mkdir -p $(lib)/
libYacs++.a: $(OBJC) $(lib)/
	ar r $@ $(OBJC)
	cp $@ $(lib)/$@
clean:
	$(MAKE) -C ../../Skin++/src clean
	/bin/rm -f *.o *.a sub-make temp.
