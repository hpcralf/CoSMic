################################################################################
##
## Process this file with automake to produce Makefile.in
##
CoSMic_FCFLAGS = -fopenmp
CoSMic_LDFLAGS = -fopenmp

bin_PROGRAMS  = CoSMic

CoSMic_SOURCES = mod_precision.f90 mod_globals.f90  \
                 mod_times.f90 \
                 mod_strings.f90 \
                 mod_urandom.f90 \
		 mod_pt_precision.f90 mod_pt_constants.f90 mod_pt_types.f90 mod_pt.f90 \
                 mod_list_variable.f90  quicksort_non_recur.f90 mod_qsort_c_module.f90 \
                 mod_support_fun.f90 mod_CoSMic_IO.f90 mod_CoSMic.f90 mod_GA.f90  \
                 CoSMic.f90

distclean-local:
	rm -f *~
	rm -f *.mod

clean-local:
	rm -f *.mod
