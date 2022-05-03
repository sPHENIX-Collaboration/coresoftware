AC_INIT(nodedump,[1.0])
AC_CONFIG_SRCDIR([configure.ac])

AM_INIT_AUTOMAKE

AC_PROG_CXX(CC g++)
LT_INIT([disable-static])

case $CXX in
 clang++)
  CXXFLAGS="$CXXFLAGS -Wall -Werror -Wextra"
 ;;
 *g++)
  CXXFLAGS="$CXXFLAGS -Wall -Werror -Wextra"
 ;;
esac

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
