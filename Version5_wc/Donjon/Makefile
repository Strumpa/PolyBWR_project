#---------------------------------------------------------------------------
#
#  Makefile for executing the Donjon non-regression tests
#  Author : A. Hebert (2018-5-10)
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
all :
	$(MAKE) -C src
clean :
	$(MAKE) clean -C src
tests :
	./rdonjon -c $(fcompilerSuite) -p $(nomp) -q Reflector_beavrs_DF-RT_Jeff3p1SHEM295_mco.x2m
ifeq ($(hdf5),1)
	./rdonjon -c $(fcompilerSuite) -p $(nomp) -q Reflector_beavrs_APEX_DF-RT.x2m
	./rdonjon -c $(fcompilerSuite) -p $(nomp) -q Reflector_beavrs_DF-NEM.x2m
	./rdonjon -c $(fcompilerSuite) -p $(nomp) -q Reflector_beavrs_DF-ANM.x2m
	./rdonjon -c $(fcompilerSuite) -p $(nomp) -q Reflector_beavrs_DF-RT.x2m
	./rdonjon -c $(fcompilerSuite) -p $(nomp) -q Reflector_beavrs_DF-RT-SPN.x2m
	./rdonjon -c $(fcompilerSuite) -p $(nomp) -q Reflector_beavrs_ERM-NEM.x2m
	./rdonjon -c $(fcompilerSuite) -p $(nomp) -q Reflector_beavrs_ERM-ANM.x2m
	./rdonjon -c $(fcompilerSuite) -p $(nomp) -q Reflector_beavrs_Koebke.x2m
endif
