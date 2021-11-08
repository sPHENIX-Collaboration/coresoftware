AC_INIT(g4tpc,[1.00])
AC_CONFIG_SRCDIR([configure.ac])
AM_INIT_AUTOMAKE

AC_PROG_CXX(CC g++)
LT_INIT([disable-static])

case $CXX in
 clang++)
  CXXFLAGS="$CXXFLAGS -Wall -Werror -Wno-undefined-var-template -Wno-misleading-indentation -Wextra"
 ;;
 *g++)
  CXXFLAGS="$CXXFLAGS -Wall -Werror -Wextra"
 ;;
esac

CINTDEFS=" -noIncludePaths  -inlineInputHeader -Wno-inconsistent-missing-override "
AC_SUBST(CINTDEFS)

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
