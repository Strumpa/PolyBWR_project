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
	./rdragon -c $(fcompilerSuite) -p $(nomp) -q iaea2d.x2m
	./rdragon -c $(fcompilerSuite) -p $(nomp) -q g2s_prestation.x2m
	./rdragon -c $(fcompilerSuite) -p $(nomp) -q salmacro.x2m
	./rdragon -c $(fcompilerSuite) -p $(nomp) -q tmacro.x2m
	./rdragon -c $(fcompilerSuite) -p $(nomp) -q VanDerGucht.x2m
	./rdragon -c $(fcompilerSuite) -p $(nomp) -q VanDerGucht-295.x2m
	./rdragon -c $(fcompilerSuite) -p $(nomp) -q uo2_295_kec1_openMP.x2m 
	./rdragon -c $(fcompilerSuite) -p $(nomp) -q tdraglib.x2m
	./rdragon -c $(fcompilerSuite) -p $(nomp) -q twimsE.x2m
	./rdragon -c $(fcompilerSuite) -p $(nomp) -q twlup.x2m
	./rdragon -c $(fcompilerSuite) -p $(nomp) -q tndas.x2m
	./rdragon -c $(fcompilerSuite) -p $(nomp) -q tmatxs2.x2m
	./rdragon -c $(fcompilerSuite) -p $(nomp) -q trowland_shem295_jeff3.1.x2m
	./rdragon -c $(fcompilerSuite) -p $(nomp) -q uo2_kec1_ecco1962_light.x2m
	./rdragon -c $(fcompilerSuite) -p $(nomp) -q lumpSS.x2m
	./rdragon -c $(fcompilerSuite) -p $(nomp) -q OSC_CASEA_1level_rse.x2m
	./rdragon -c $(fcompilerSuite) -p $(nomp) -q OSC_CASEA_2level_rse.x2m
	./rdragon -c $(fcompilerSuite) -p $(nomp) -q OSC_openMP_tiso.x2m
	./rdragon -c $(fcompilerSuite) -p $(nomp) -q sens.x2m
	./rdragon -c $(fcompilerSuite) -p $(nomp) -q testVVER7.x2m
	./rdragon -c $(fcompilerSuite) -p $(nomp) -q fbr_colorset.x2m
	./rdragon -c $(fcompilerSuite) -p $(nomp) -q fbr_tone.x2m
	./rdragon -c $(fcompilerSuite) -p $(nomp) -q TEST_GEO_hex_sect_tspc.x2m
	./rdragon -c $(fcompilerSuite) -p $(nomp) -q TEST_GEO_latt_tspc_S30.x2m
	./rdragon -c $(fcompilerSuite) -p $(nomp) -q testDuo.x2m
	./rdragon -c $(fcompilerSuite) -p $(nomp) -q testDuo_B1.x2m
	./rdragon -c $(fcompilerSuite) -p $(nomp) -q uo2_evo_xsm.x2m
	./rdragon -c $(fcompilerSuite) -p $(nomp) -q C2D20.x2m
	./rdragon -c $(fcompilerSuite) -p $(nomp) -q CFC-CELL.x2m
	./rdragon -c $(fcompilerSuite) -p $(nomp) -q ErmBeavrsPwrRefl.x2m
	./rdragon -c $(fcompilerSuite) -p $(nomp) -q DF_RTBeavrsPwrRefl.x2m
	./rdragon -c $(fcompilerSuite) -p $(nomp) -q rep900_het_gff_jef2p2.x2m
	./rdragon -c $(fcompilerSuite) -p $(nomp) -q RegtestCNG_mccg.x2m
	./rdragon -c $(fcompilerSuite) -p $(nomp) -q RegtestLZC_mccg.x2m
	./rdragon -c $(fcompilerSuite) -p $(nomp) -q pincell_mco.x2m
	./rdragon -c $(fcompilerSuite) -p $(nomp) -q pincell_sap.x2m
ifeq ($(hdf5),1)
	./rdragon -c $(fcompilerSuite) -p $(nomp) -q UOX_5x5_TG6_sym8_multiDom.x2m
	./rdragon -c $(fcompilerSuite) -p $(nomp) -q uo2_evo_hdf.x2m
	./rdragon -c $(fcompilerSuite) -p $(nomp) -q pincell_mpo.x2m
	./rdragon -c $(fcompilerSuite) -p $(nomp) -q ASSBLY_CASEA_1level_mpo.x2m
	./rdragon -c $(fcompilerSuite) -p $(nomp) -q ASSBLY_VVER_1level_mpo.x2m
endif
