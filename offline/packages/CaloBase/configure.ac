AC_INIT(calobase,[1.00])
AC_CONFIG_SRCDIR([configure.ac])

AM_INIT_AUTOMAKE
AC_PROG_CXX(CC g++)
LT_INIT([disable-static])

case $CXX in
 clang++)
  CXXFLAGS="$CXXFLAGS -Wall  -Werror -Wextra -Wshadow"
 ;;
 *g++)
  CXXFLAGS="$CXXFLAGS -Wall -Werror -pedantic -Wextra -Wshadow"
 ;;
esac

CINTDEFS=" -noIncludePaths  -inlineInputHeader "
AC_SUBST(CINTDEFS)

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
