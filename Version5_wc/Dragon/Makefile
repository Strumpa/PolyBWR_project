#---------------------------------------------------------------------------
#
#  Makefile for executing the Dragon non-regression tests
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
	./rdragon -c $(fcompilerSuite) -q iaea2d.x2m
	./rdragon -c $(fcompilerSuite) -q g2s_prestation.x2m
	./rdragon -c $(fcompilerSuite) -q salmacro.x2m
	./rdragon -c $(fcompilerSuite) -q tmacro.x2m
	./rdragon -c $(fcompilerSuite) -q VanDerGucht.x2m
	./rdragon -c $(fcompilerSuite) -q VanDerGucht-295.x2m
	./rdragon -c $(fcompilerSuite) -q uo2_295_kec1_openMP.x2m 
	./rdragon -c $(fcompilerSuite) -q tdraglib.x2m
	./rdragon -c $(fcompilerSuite) -q twimsE.x2m
	./rdragon -c $(fcompilerSuite) -q twlup.x2m
	./rdragon -c $(fcompilerSuite) -q tndas.x2m
	./rdragon -c $(fcompilerSuite) -q tmatxs2.x2m
	./rdragon -c $(fcompilerSuite) -q trowland_shem295_jeff3.1.x2m
	./rdragon -c $(fcompilerSuite) -q lumpSS.x2m
	./rdragon -c $(fcompilerSuite) -q OSC_CASEA_1level_rse.x2m
	./rdragon -c $(fcompilerSuite) -q OSC_CASEA_2level_rse.x2m
	./rdragon -c $(fcompilerSuite) -q OSC_openMP_tiso.x2m
	./rdragon -c $(fcompilerSuite) -q sens.x2m
	./rdragon -c $(fcompilerSuite) -q testVVER7.x2m
	./rdragon -c $(fcompilerSuite) -q fbr_colorset.x2m
	./rdragon -c $(fcompilerSuite) -q fbr_tone.x2m
	./rdragon -c $(fcompilerSuite) -q TEST_GEO_hex_sect_tspc.x2m
	./rdragon -c $(fcompilerSuite) -q testDuo.x2m
	./rdragon -c $(fcompilerSuite) -q testDuo_B1.x2m
	./rdragon -c $(fcompilerSuite) -q uo2_evo_xsm.x2m
	./rdragon -c $(fcompilerSuite) -q C2D20.x2m
	./rdragon -c $(fcompilerSuite) -q CFC-CELL.x2m
	./rdragon -c $(fcompilerSuite) -q ErmBeavrsPwrRefl.x2m
	./rdragon -c $(fcompilerSuite) -q rep900_het_gff_jef2p2.x2m
	./rdragon -c $(fcompilerSuite) -q RegtestCNG_mccg.x2m
	./rdragon -c $(fcompilerSuite) -q RegtestLZC_mccg.x2m
	./rdragon -c $(fcompilerSuite) -q pincell_mco.x2m
	./rdragon -c $(fcompilerSuite) -q pincell_sap.x2m
ifeq ($(hdf5),1)
	./rdragon -c $(fcompilerSuite) -q UOX_5x5_TG6_sym8_multiDom.x2m
	./rdragon -c $(fcompilerSuite) -q uo2_evo_hdf.x2m
	./rdragon -c $(fcompilerSuite) -q pincell_mpo.x2m
	./rdragon -c $(fcompilerSuite) -q ASSBLY_CASEA_1level_mpo.x2m
	./rdragon -c $(fcompilerSuite) -q ASSBLY_VVER_1level_mpo.x2m
endif
