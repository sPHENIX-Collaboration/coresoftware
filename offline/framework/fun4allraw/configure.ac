AC_INIT(fun4allraw,[2.00])
AC_CONFIG_SRCDIR([configure.ac])

AM_INIT_AUTOMAKE
AC_PROG_CXX(CC g++)

LT_INIT([disable-static])

dnl   unused param warnings come from newbasic include
if test $ac_cv_prog_gxx = yes; then
   CXXFLAGS="$CXXFLAGS -Wall -Werror -Wno-unused-parameter"
fi


AC_CONFIG_FILES([Makefile])
AC_OUTPUT
