AC_INIT(sHijing, [1.00])
AC_CONFIG_SRCDIR([configure.ac])

AM_INIT_AUTOMAKE

AC_PROG_CXX(CC g++)
AC_PROG_F77(gfortran)
LT_INIT([disable-static])

dnl test for root 6
if test `root-config --version | gawk '{print $1>=6.?"1":"0"}'` = 1; then
CINTDEFS=" -noIncludePaths  -inlineInputHeader "
AC_SUBST(CINTDEFS)
fi

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
