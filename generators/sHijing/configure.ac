AC_INIT(sHijing, [1.00])
AC_CONFIG_SRCDIR([configure.ac])

AM_INIT_AUTOMAKE

AC_PROG_CXX(CC g++)
AC_PROG_F77(gfortran)
LT_INIT([disable-static])

dnl   no point in suppressing warnings people should
dnl   at least see them, so here we go for g++: -Wall
dnl   make warnings fatal errors: -Werror
if test $ac_cv_prog_gxx = yes; then
dnl turn off warnings about unused functions in cfortran.h
   CXXFLAGS="$CXXFLAGS -Wall  -Werror -Wextra -Wno-unused-function -Wno-deprecated-declarations"
fi

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
