AC_INIT(g4vertex, [1.00])
AC_CONFIG_SRCDIR([configure.ac])

AM_INIT_AUTOMAKE

AC_PROG_CXX(CC g++)
LT_INIT([disable-static])

dnl leaving this here in case we want to play with different compiler 
dnl specific flags
case $CXX in
 clang++)
  CXXFLAGS="$CXXFLAGS -Wall -Werror -Wextra"
 ;;
 *g++)
  CXXFLAGS="$CXXFLAGS -Wall -Werror -Wextra"
 ;;
esac


CINTDEFS=" -noIncludePaths  -inlineInputHeader "
AC_SUBST(CINTDEFS)

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
