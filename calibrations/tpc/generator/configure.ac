AC_INIT(tpcdistgenerator,[2.0])
AC_CONFIG_SRCDIR([configure.ac])

AM_INIT_AUTOMAKE
AC_PROG_F77(gfortran)
AC_PROG_CXX(CC g++)

LT_INIT([disable-static])

dnl   no point in suppressing warnings people should 
dnl   at least see them, so here we go for g++: -Wall
if test $ac_cv_prog_gxx = yes; then
  CXXFLAGS="$CXXFLAGS -Wall -Werror -Wextra"
fi

case $CXX in
 clang++)
  CXXFLAGS="$CXXFLAGS -Wno-unused-parameter -Wno-bitwise-instead-of-logical -Wno-unused-but-set-variable -Wno-unknown-warning-option"
 ;;
esac

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
