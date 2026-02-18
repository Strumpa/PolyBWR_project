
#---------------------------------------------------------------------------
#
#  Makefile for executing the Pygan non-regression tests
#  Author : A. Hebert (2020-7-3)
#
#---------------------------------------------------------------------------
#
OS = $(shell uname -s | cut -d"_" -f1)
ifneq (,$(filter $(OS),SunOS AIX))
  MAKE = gmake
endif
ifeq ($(openmp),1)
  nomp = 16
else
  nomp = 0
endif
ifeq ($(intel),1)
  fcompilerSuite = intel
else
  ifeq ($(nvidia),1)
    fcompilerSuite = nvidia
  else
    ifeq ($(llvm),1)
      fcompilerSuite = llvm
    else
      fcompilerSuite = custom
    endif
  endif
endif

ifeq ($(pip),1)
  MAKE_PY = $(MAKE) -f Makefile_pip
  options = -p $(nomp) -q -e
else
  MAKE_PY = $(MAKE)
  options = -p $(nomp) -q
endif

all :
	$(MAKE_PY) donjon -C src
ganlib :
	$(MAKE_PY) ganlib -C src
trivac :
	$(MAKE_PY) trivac -C src
dragon :
	$(MAKE_PY) dragon -C src
donjon :
	$(MAKE_PY) donjon -C src
clean :
	$(MAKE_PY) clean -C src
tests :
	./rpython $(options) -c $(fcompilerSuite) test_lcm.py
	./rpython $(options) -c $(fcompilerSuite) test_lifo.py
	./rpython $(options) -c $(fcompilerSuite) fact.py
	./rpython $(options) -c $(fcompilerSuite) Exemple_jdd.py
	./rpython $(options) -c $(fcompilerSuite) simplePOW.py
	./rpython $(options) -c $(fcompilerSuite) testVVER.py
ifeq ($(hdf5),1)
	./rpython $(options) -c $(fcompilerSuite) Equivalence_SPH_APEX.py
	./rpython $(options) -c $(fcompilerSuite) pincell_mpo_concat.py
endif
ifeq ($(vtk),1)
	./rpython $(options) -c $(fcompilerSuite) rDragView.py
endif
