AC_INIT(eastphysicslist, [1.00])
AC_CONFIG_SRCDIR([configure.ac])

AM_INIT_AUTOMAKE

AC_PROG_CXX(CC g++)
LT_INIT([disable-static])

CXXFLAGS="$CXXFLAGS -Wall -Werror -Wextra -Wshadow"

case $CXX in
  clang++)
    CXXFLAGS="$CXXFLAGS -Wno-final-dtor-non-final-class"
  ;;
esac
dnl   AM_CONDITIONAL(GCC_GE_48, test `g++ -dumpversion | awk '{print $1>=4.8?"1":"0"}'` = 1)

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
