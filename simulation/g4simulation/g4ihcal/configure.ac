AC_INIT(g4ihcal, [1.00])
AC_CONFIG_SRCDIR([configure.ac])

AM_INIT_AUTOMAKE

AC_PROG_CXX(CC g++)
LT_INIT([disable-static])

CXXFLAGS="$CXXFLAGS -Wall -Werror -Wextra -Wshadow"

case $CXX in
  clang++)
    CXXFLAGS="$CXXFLAGS -Wno-c11-extensions"
  ;;
esac

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
