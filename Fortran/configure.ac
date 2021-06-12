AC_INIT([CoSMic], [1.0], [ralf.schneider@hlrs.de])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])


AC_CONFIG_MACRO_DIRS([/usr/share/aclocal])
AX_PROG_FC_MPI(
   [test 0 == 0],
   [echo "MPI compiler found and set."],
   [echo "!! ERROR - No MPI compiler found !!";exit 1])

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([
 Makefile
 src/Makefile
])
AC_OUTPUT
