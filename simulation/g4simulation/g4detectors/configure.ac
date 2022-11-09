AC_INIT(g4detectors, [1.00])
AC_CONFIG_SRCDIR([configure.ac])

AM_INIT_AUTOMAKE
AC_PROG_CXX(CC g++)
LT_INIT([disable-static])

if test $ac_cv_prog_gxx = yes; then
  CXXFLAGS="-Wall -Wextra -Werror -Wshadow $CXXFLAGS"
fi

case $CXX in
  clang++)
    CXXFLAGS="$CXXFLAGS -Wno-c11-extensions -Wno-deprecated-copy -Wno-tautological-overlap-compare -Wno-null-pointer-subtraction -Wno-unknown-warning-option"
  ;;
esac

CINTDEFS=" -noIncludePaths  -inlineInputHeader "
AC_SUBST(CINTDEFS)

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
