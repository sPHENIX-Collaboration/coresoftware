AC_INIT(phgeom,[1.00])
AC_CONFIG_SRCDIR([configure.ac])

AM_INIT_AUTOMAKE
AC_PROG_CXX(CC g++)
LT_INIT([disable-static])

if test $ac_cv_prog_gxx = yes; then
  CXXFLAGS="$CXXFLAGS -Wall -Werror -Wextra"
  AM_CONDITIONAL(GCC_GE_48, test `g++ -dumpversion | awk '{print $1>=4.8?"1":"0"}'` = 1)
fi

CINTDEFS=" -noIncludePaths  -inlineInputHeader "
AC_SUBST(CINTDEFS)


AC_CONFIG_FILES([Makefile])
AC_OUTPUT
