#---------------------------------------------------------------------------
#
#  Makefile for executing the Trivac non-regression tests
#  Author : A. Hebert (2018-5-10)
#
#---------------------------------------------------------------------------
#
OS = $(shell uname -s | cut -d"_" -f1)
ifneq (,$(filter $(OS),SunOS AIX))
  MAKE = gmake
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
	./rtrivac -c $(fcompilerSuite) -q DIFtst.x2m
	./rtrivac -c $(fcompilerSuite) -q Ktests.x2m
	./rtrivac -c $(fcompilerSuite) -q SPNtst.x2m
	./rtrivac -c $(fcompilerSuite) -q iaea3d.x2m
	./rtrivac -c $(fcompilerSuite) -q monju3D_spn.x2m
	./rtrivac -c $(fcompilerSuite) -q multigroup_albedo_2d.x2m
	./rtrivac -c $(fcompilerSuite) -q multigroup_albedo_3d.x2m
	./rtrivac -c $(fcompilerSuite) -q takedaM4_spn.x2m
	./rtrivac -c $(fcompilerSuite) -q NodalTests.x2m
