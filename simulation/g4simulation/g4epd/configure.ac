AC_INIT(g4epd, [1.00])
AC_CONFIG_SRCDIR([configure.ac])

AM_INIT_AUTOMAKE
AC_PROG_CXX(CC g++)

LT_INIT([disable-static])

CXXFLAGS="$CXXFLAGS -Wall -Werror -Wextra -Wshadow"

case $CXX in
 clang++)
  CXXFLAGS="$CXXFLAGS -Wno-undefined-var-template"
 ;;
esac

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
