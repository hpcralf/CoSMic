AC_INIT([CoSMic], [1.0], [ralf.schneider@hlrs.de])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])

# If --with-mpi=auto is used, try to find MPI, but use standard FC compiler if it is not found.
# If --with-mpi=yes is used, try to find MPI and fail if it isn't found.
# If --with-mpi=no is used, use a standard FC compiler instead.
AC_ARG_WITH(mpi, [AS_HELP_STRING([--with-mpi],
    [compile with MPI (parallelization) support. If none is found,
    MPI is not used. Default: auto])
],,[with_mpi=yes])
  
AC_OPENMP

AX_PROG_FC_MPI([test x"$with_mpi" != xno],[use_mpi=yes],[
  use_mpi=no
  if test x"$with_mpi" = xyes; then
    AC_MSG_FAILURE([MPI compiler requested and mandatory, but couldn't use MPI.])
  else
    AC_MSG_FAILURE([No MPI compiler found or --with-mpi=no was given but MPI is mandatory.])
  fi
])

AC_CONFIG_HEADERS([config.h])

AC_CONFIG_FILES([
 Makefile
 src/Makefile
])

AC_OUTPUT
