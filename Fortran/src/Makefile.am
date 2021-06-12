################################################################################
##
## Process this file with automake to produce Makefile.in
##
AM_LD = $MPIFC
AM_FCFLAGS = -Jmod

bin_PROGRAMS  = CoSMic

CoSMic_SOURCES = mod_precision.f90 mod_global_constants.f90 mod_global_types.f90 \
                 mod_list_variable.f90  quicksort_non_recur.f90 mod_qsort_c_module.f90 \
                 mod_support_fun.f90 mod_kernel.f90 mod_CoSMic_IO.f90 \
                 mod_data_preprocessing.f90 CoSMic.f90

clean-local:
	rm -f *.mod mod/*.mod

distclean-local:
	rm -f *~

