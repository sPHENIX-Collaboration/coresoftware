AC_INIT(calib_emc_pi0,[1.00])
AC_CONFIG_SRCDIR([configure.ac])

AM_INIT_AUTOMAKE
AC_PROG_CXX(CC g++)

LT_INIT([disable-static])

dnl   no point in suppressing warnings people should 
dnl   at least see them, so here we go for g++: -Wall
if test $ac_cv_prog_gxx = yes; then
   CXXFLAGS="-Wall -Werror -Wextra -Wshadow $CXXFLAGS"
fi

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
