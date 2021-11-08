AC_INIT(hijing, [1.383bs.2])
AC_CONFIG_SRCDIR(src/hijing.f)
AM_INIT_AUTOMAKE

AC_PROG_CXX
AC_PROG_F77(gfortran)
LT_INIT([disable-static])

dnl deprecated declaraion from fastjet header
if test $ac_cv_prog_gxx = yes; then
  CXXFLAGS="$CXXFLAGS -Wall -Werror -Wno-deprecated-declarations -Wextra"
fi

CINTDEFS=" -noIncludePaths  -inlineInputHeader -Wno-inconsistent-missing-override -Wextra "

AC_SUBST(CINTDEFS)

AC_CONFIG_FILES([Makefile])
AC_OUTPUT

