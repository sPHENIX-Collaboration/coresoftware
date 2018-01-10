AC_INIT(configure.ac)

AM_INIT_AUTOMAKE(g4tpc, 1.00)
AC_PROG_CXX(CC g++)
AC_ENABLE_STATIC(no)
LT_INIT

if test $ac_cv_prog_gxx = yes; then
  CXXFLAGS="$CXXFLAGS -Wall -Werror"
dnl  AM_CONDITIONAL(GCC_GE_48, test `g++ -dumpversion | gawk '{print $1>=4.8?"1":"0"}'` = 1)
fi

dnl test for root 6
if test `root-config --version | gawk '{print $1>=6.?"1":"0"}'` = 1; then
CINTDEFS=" -noIncludePaths  -inlineInputHeader "
AC_SUBST(CINTDEFS)
fi

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
