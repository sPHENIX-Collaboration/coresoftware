AC_INIT(g4micromegas,[1.00])
AC_CONFIG_SRCDIR([configure.ac])

AM_INIT_AUTOMAKE
AC_PROG_CXX(CC g++)

LT_INIT([disable-static])

case $CXX in
 clang++)
  CXXFLAGS="$CXXFLAGS -Wall -Werror -Wno-undefined-var-template -Wno-return-stack-address"
 ;;
 g++)
  CXXFLAGS="$CXXFLAGS -Wall -Werror"
 ;;
esac

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
