AC_INIT(PHPythia8,[1.00])
AC_CONFIG_SRCDIR([configure.ac])
AM_INIT_AUTOMAKE

AC_PROG_CXX(CC g++)
LT_INIT([disable-static])

dnl   no point in suppressing warnings people should
dnl   at least see them, so here we go for g++: -Wall
dnl   make warnings fatal errors: -Werror
dnl leaving this here in case we want to play with different compiler 
dnl specific flags
if test $ac_cv_prog_gxx = yes; then
   CXXFLAGS="$CXXFLAGS -Wall -pedantic -Wextra -Werror -Wshadow"
fi

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
