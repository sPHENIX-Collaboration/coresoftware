AC_INIT(phoolraw, [2.0])
AC_CONFIG_SRCDIR([configure.ac])

AM_INIT_AUTOMAKE
AC_PROG_CXX(CC g++)
LT_INIT([disable-static])

dnl   unused param warnings come from newbasic include
if test $ac_cv_prog_gxx = yes; then
  CXXFLAGS="$CXXFLAGS -Wall -Wextra -Werror -Wno-unused-parameter"
fi

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
