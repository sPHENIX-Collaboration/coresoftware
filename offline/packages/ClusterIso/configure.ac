AC_INIT(configure.ac)

AM_INIT_AUTOMAKE(ClusterIso, 1.00)
AC_PROG_CXX(CC g++)
AC_ENABLE_STATIC(no)
LT_INIT

if test $ac_cv_prog_gxx = yes; then
  CXXFLAGS="$CXXFLAGS -Wall -Werror"
fi
CXXFLAGS="$CXXFLAGS"
ROOTLIBS=`root-config --libs`
AC_SUBST(ROOTLIBS)

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
