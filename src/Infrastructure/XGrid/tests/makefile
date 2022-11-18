# $Id$ 

ALL:      build_unit_tests
run:      run_unit_tests
run_uni:  run_unit_tests_uni

LOCDIR	  = src/Infrastructure/XGrid/tests

.NOTPARALLEL:
TESTS_BUILD   = $(ESMF_TESTDIR)/ESMF_XGridUTest \
	        $(ESMF_TESTDIR)/ESMF_AltXGridUTest \
	        $(ESMF_TESTDIR)/ESMF_XGridMaskingUTest \
	        $(ESMF_TESTDIR)/ESMC_XGridUTest \
	        $(ESMF_TESTDIR)/ESMF_XGridOffUTest \
	        $(ESMF_TESTDIR)/ESMF_TimeStepUTest

TESTS_RUN     = RUN_ESMF_XGridUTest RUN_ESMF_AltXGridUTest RUN_ESMF_XGridMaskingUTest RUN_ESMC_XGridUTest RUN_ESMF_XGridOffUTest RUN_ESMF_TimeStepUTest

TESTS_RUN_UNI = RUN_ESMF_XGridUTestUNI RUN_ESMF_AltXGridUTestUNI RUN_ESMF_XGridMaskingUTestUNI RUN_ESMC_XGridUTestUNI RUN_ESMD_XGridOffUTestUNI RUN_ESMF_TimeStepUTestUNI

include ${ESMF_DIR}/makefile

CLEANDIRS   = 
CLEANFILES  = $(TESTS_BUILD) $(CPFILES)
CLOBBERDIRS =

DIRS        = 


#
# XGrid
#

RUN_ESMF_XGridUTest:
	$(MAKE) TNAME=XGrid NP=2 ftest

RUN_ESMF_XGridUTestUNI:
	$(MAKE) TNAME=XGrid NP=1 ftest

RUN_ESMF_AltXGridUTest:
	$(MAKE) TNAME=AltXGrid NP=2 ftest

RUN_ESMF_AltXGridUTestUNI:
	$(MAKE) TNAME=AltXGrid NP=1 ftest

RUN_ESMF_XGridMaskingUTest:
	$(MAKE) TNAME=XGridMasking NP=2 ftest

RUN_ESMF_XGridMaskingUTestUNI:
	$(MAKE) TNAME=XGridMasking NP=1 ftest

RUN_ESMC_XGridUTest:
	$(MAKE) TNAME=XGrid NP=4 ctest

RUN_ESMC_XGridUTestUNI:
	$(MAKE) TNAME=XGrid NP=1 ctest

RUN_ESMF_XGridOffUTest:
	$(MAKE) TNAME=XGridOff NP=2 ftest

RUN_ESMF_XGridOffUTestUNI:
	$(MAKE) TNAME=XGridOff NP=1 ftest

RUN_ESMF_TimeStepUTest:
	$(MAKE) TNAME=TimeStep NP=2 ftest

RUN_ESMF_TimeStepUTestUNI:
	$(MAKE) TNAME=TimeStep NP=1 ftest
