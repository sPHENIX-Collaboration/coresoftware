AC_INIT(hijing, [1.383bs.2])
AC_CONFIG_SRCDIR(src/hijing.f)
AM_INIT_AUTOMAKE

AC_PROG_CXX
AC_PROG_F77(gfortran)
LT_INIT([disable-static])

if test $ac_cv_prog_gxx = yes; then
  CXXFLAGS="$CXXFLAGS -Wno-deprecated-declarations"
fi

CINTDEFS=" -noIncludePaths  -inlineInputHeader -Wno-inconsistent-missing-override"
AC_SUBST(CINTDEFS)

AC_CONFIG_FILES([Makefile])
AC_OUTPUT

