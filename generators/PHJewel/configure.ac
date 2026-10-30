AC_INIT(PHJewel,[1.00])
AC_CONFIG_SRCDIR([configure.ac])
AM_INIT_AUTOMAKE([subdir-objects])

AC_PROG_CXX(CC g++)
AC_PROG_F77([gfortran])
LT_INIT([disable-static])

if test $ac_cv_prog_gxx = yes; then
   CXXFLAGS="$CXXFLAGS -Wall -pedantic -Wextra -Werror -Wshadow"
fi

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
