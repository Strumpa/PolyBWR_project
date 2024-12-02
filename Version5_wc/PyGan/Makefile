
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
ifeq ($(intel),1)
  fcompiler = intel
else
  ifeq ($(nvidia),1)
    fcompiler = nvidia
  else
    ifeq ($(flang),1)
      fcompiler = flang
    else
      fcompiler = custom
    endif
  endif
endif

ifeq ($(pip),1)
  MAKE_PY = $(MAKE) -f Makefile_pip
  options = -q -e
else
  MAKE_PY = $(MAKE)
  options = -q
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
	./rpython $(options) -c $(fcompiler) test_lcm.py
	./rpython $(options) -c $(fcompiler) test_lifo.py
	./rpython $(options) -c $(fcompiler) fact.py
	./rpython $(options) -c $(fcompiler) Exemple_jdd.py
	./rpython $(options) -c $(fcompiler) simplePOW.py
	./rpython $(options) -c $(fcompiler) testVVER.py
ifeq ($(hdf5),1)
	./rpython $(options) -c $(fcompiler) Equivalence_SPH_APEX.py
	./rpython $(options) -c $(fcompiler) pincell_mpo_concat.py
endif
