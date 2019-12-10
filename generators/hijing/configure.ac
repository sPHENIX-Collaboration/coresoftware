AC_INIT(hijing, [1.383bs.2])
AC_CONFIG_SRCDIR(src/hijing.f)
AM_INIT_AUTOMAKE

AC_PROG_CXX
AC_PROG_F77(gfortran)
LT_INIT([disable-static])

dnl test for root 6
if test `root-config --version | awk '{print $1>=6.?"1":"0"}'` = 1; then
CINTDEFS=" -noIncludePaths  -inlineInputHeader "
AC_SUBST(CINTDEFS)
fi

AC_CONFIG_FILES([Makefile])
AC_OUTPUT

