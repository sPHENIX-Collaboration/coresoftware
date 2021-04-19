AC_INIT(g4trackfastsim, [1.00])
AC_CONFIG_SRCDIR([configure.ac])

AM_INIT_AUTOMAKE

AC_PROG_CXX(CC g++)
LT_INIT([disable-static])

dnl leaving this here in case we want to play with different compiler 
dnl specific flags
case $CXX in
 clang++)
  CXXFLAGS="$CXXFLAGS -Wall -Werror"
 ;;
 g++)
  CXXFLAGS="$CXXFLAGS -Wall -Werror"
 ;;
esac


AC_CONFIG_FILES([Makefile])
AC_OUTPUT
