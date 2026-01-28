AC_INIT( dEdxFitter,[1.00])
AC_CONFIG_SRCDIR([configure.ac])

AM_INIT_AUTOMAKE
AC_PROG_CXX(CC g++)

LT_INIT([disable-static])

dnl   enable more warnings and make them fatal
dnl this package needs openmp which requires -fopenmp for clang
if test $ac_cv_prog_gxx = yes; then
   CXXFLAGS="$CXXFLAGS -Wall -Wshadow -Wextra -Werror"
fi

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
