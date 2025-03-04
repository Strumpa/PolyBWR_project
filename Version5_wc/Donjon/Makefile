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
	./rdonjon -c $(fcompilerSuite) -q 1Sap_vs_3Sap.x2m
	./rdonjon -c $(fcompilerSuite) -q AFMtest.x2m
	./rdonjon -c $(fcompilerSuite) -q Candu6.x2m
	./rdonjon -c $(fcompilerSuite) -q candu6_iram.x2m
	./rdonjon -c $(fcompilerSuite) -q Equivalence_SPH_SN_MCO.x2m
	./rdonjon -c $(fcompilerSuite) -q Example1.x2m
	./rdonjon -c $(fcompilerSuite) -q Example2.x2m
	./rdonjon -c $(fcompilerSuite) -q Fessenheim.x2m
	./rdonjon -c $(fcompilerSuite) -q FourChannels.x2m
	./rdonjon -c $(fcompilerSuite) -q Hex_Model3D.x2m
	./rdonjon -c $(fcompilerSuite) -q NCRtest.x2m
	./rdonjon -c $(fcompilerSuite) -q SCR_RES_Interp.x2m
	./rdonjon -c $(fcompilerSuite) -q SCWR64Core.x2m
	./rdonjon -c $(fcompilerSuite) -q TC_EVO.x2m
	./rdonjon -c $(fcompilerSuite) -q TC_XENON.x2m
	./rdonjon -c $(fcompilerSuite) -q TINST.x2m
	./rdonjon -c $(fcompilerSuite) -q Test4x4coupl_nompi.x2m
	./rdonjon -c $(fcompilerSuite) -q Test_microdepl_NCR.x2m
	./rdonjon -c $(fcompilerSuite) -q Test_IDET.x2m
	./rdonjon -c $(fcompilerSuite) -q Test_PMAXS_Macro.x2m
	./rdonjon -c $(fcompilerSuite) -q Test_PMAXS_Map.x2m
	./rdonjon -c $(fcompilerSuite) -q Test_SPH.x2m
	./rdonjon -c $(fcompilerSuite) -q Test_SPH_void_diff_picard.x2m
	./rdonjon -c $(fcompilerSuite) -q Test_SPH_void.x2m
	./rdonjon -c $(fcompilerSuite) -q Test_SPH_void_bfgs.x2m
	./rdonjon -c $(fcompilerSuite) -q Test_SPH_void_lbfgs.x2m
	./rdonjon -c $(fcompilerSuite) -q Test_SPH_void_newton.x2m
	./rdonjon -c $(fcompilerSuite) -q Triga_Model3D.x2m
	./rdonjon -c $(fcompilerSuite) -q VHTR_SPH.x2m
	./rdonjon -c $(fcompilerSuite) -q VHTRc.x2m
	./rdonjon -c $(fcompilerSuite) -q channel_mphy.x2m
	./rdonjon -c $(fcompilerSuite) -q d2p_mco.x2m
	./rdonjon -c $(fcompilerSuite) -q d2p_rfl.x2m
	./rdonjon -c $(fcompilerSuite) -q d2p_sap.x2m
	./rdonjon -c $(fcompilerSuite) -q fbr_core.x2m
	./rdonjon -c $(fcompilerSuite) -q history_mphy.x2m
	./rdonjon -c $(fcompilerSuite) -q pincell_mphy_thm.x2m
	./rdonjon -c $(fcompilerSuite) -q pmox_sap.x2m
	./rdonjon -c $(fcompilerSuite) -q pulseTHM_0d.x2m
	./rdonjon -c $(fcompilerSuite) -q rep1300_reflector.x2m
	./rdonjon -c $(fcompilerSuite) -q rep1300_reflector_sp1.x2m
	./rdonjon -c $(fcompilerSuite) -q rep900.x2m
	./rdonjon -c $(fcompilerSuite) -q rep900_anm.x2m
	./rdonjon -c $(fcompilerSuite) -q rep900cluster.x2m
	./rdonjon -c $(fcompilerSuite) -q rep900EnrichCOMPO.x2m
	./rdonjon -c $(fcompilerSuite) -q rep900EnrichCOMPOhet.x2m
	./rdonjon -c $(fcompilerSuite) -q rep900_iram.x2m
	./rdonjon -c $(fcompilerSuite) -q rep900_mphy.x2m
	./rdonjon -c $(fcompilerSuite) -q rep900_mphy_alt.x2m
	./rdonjon -c $(fcompilerSuite) -q rep900_msap.x2m
	./rdonjon -c $(fcompilerSuite) -q rep900_sim_defrod.x2m
	./rdonjon -c $(fcompilerSuite) -q rep900_sim_recopy.x2m
	./rdonjon -c $(fcompilerSuite) -q rep900_sim.x2m
	./rdonjon -c $(fcompilerSuite) -q ssrcell_thm.x2m
	./rdonjon -c $(fcompilerSuite) -q testNAPhet.x2m
	./rdonjon -c $(fcompilerSuite) -q test_channBWR.x2m
	./rdonjon -c $(fcompilerSuite) -q thm_salt.x2m
	./rdonjon -c $(fcompilerSuite) -q thm_salt_lf.x2m
	./rdonjon -c $(fcompilerSuite) -q thm_salt_lf_texp.x2m
	./rdonjon -c $(fcompilerSuite) -q thm_salt_lf_trans.x2m
	./rdonjon -c $(fcompilerSuite) -q vver_test.x2m
ifeq ($(hdf5),1)
	./rdonjon -c $(fcompilerSuite) -q AFA_180_310_type1_dual.x2m
	./rdonjon -c $(fcompilerSuite) -q Equivalence_SPH_RT_APEX.x2m
	./rdonjon -c $(fcompilerSuite) -q Equivalence_SPH_SN_APEX.x2m
	./rdonjon -c $(fcompilerSuite) -q Equivalence_SPH_RT_MPO.x2m
	./rdonjon -c $(fcompilerSuite) -q SmallCore_rt_mpo.x2m
	./rdonjon -c $(fcompilerSuite) -q SmallCore_anm_mpo_mpo.x2m
	./rdonjon -c $(fcompilerSuite) -q SmallCore_anm_mpo_mpo_erm.x2m
	./rdonjon -c $(fcompilerSuite) -q UOX_5x5_TG6_sym8_multiDom_rt2.x2m
	./rdonjon -c $(fcompilerSuite) -q mpo_mco_refl.x2m
	./rdonjon -c $(fcompilerSuite) -q pincell_mpo_concat.x2m
	./rdonjon -c $(fcompilerSuite) -q Reflector_beavrs_DF-NEM.x2m
	./rdonjon -c $(fcompilerSuite) -q Reflector_beavrs_DF-ANM.x2m
	./rdonjon -c $(fcompilerSuite) -q Reflector_beavrs_DF-RT.x2m
	./rdonjon -c $(fcompilerSuite) -q Reflector_beavrs_DF-RT-SPN.x2m
	./rdonjon -c $(fcompilerSuite) -q Reflector_beavrs_ERM-NEM.x2m
	./rdonjon -c $(fcompilerSuite) -q Reflector_beavrs_ERM-ANM.x2m
	./rdonjon -c $(fcompilerSuite) -q Reflector_beavrs_Koebke.x2m
	./rdonjon -c $(fcompilerSuite) -q SmallCore_BaffRefl_Benchmark_mpo.x2m
endif
