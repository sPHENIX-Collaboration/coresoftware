AC_INIT(g4tpc,[1.00])
AC_CONFIG_SRCDIR([configure.ac])
AM_INIT_AUTOMAKE

AC_PROG_CXX(CC g++)
LT_INIT([disable-static])

CXXFLAGS="$CXXFLAGS -Wall -Werror -Wextra -Wshadow"

case $CXX in
 *analyzer)
   CXXFLAGS="$CXXFLAGS -Wno-vla"
 ;;
 clang++)
   CXXFLAGS="$CXXFLAGS -Wno-vla-extension -Wno-c11-extensions -Wno-deprecated-copy"
 ;;
 *g++)
   CXXFLAGS="$CXXFLAGS -Wno-vla"
 ;;
esac

CINTDEFS=" -noIncludePaths  -inlineInputHeader -Wno-inconsistent-missing-override "
AC_SUBST(CINTDEFS)

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
