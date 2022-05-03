AC_INIT(sHEPGen,[1.00])
AC_CONFIG_SRCDIR([configure.ac])
AM_INIT_AUTOMAKE

AC_PROG_CXX(g++)
LT_INIT([disable-static])

dnl   no point in suppressing warnings people should
dnl   at least see them, so here we go for g++: -Wall
dnl   make warnings fatal errors: -Werror
if test $ac_cv_prog_gxx = yes; then
   CXXFLAGS="$CXXFLAGS -Wall -Werror -Wextra"
fi

AC_CONFIG_FILES([Makefile])
AC_OUTPUT

