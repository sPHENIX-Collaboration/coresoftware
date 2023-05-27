AC_INIT(sphenixnpc,[1.00])
AC_CONFIG_SRCDIR([configure.ac])

AM_INIT_AUTOMAKE
AC_PROG_CXX(CC g++)
LT_INIT([disable-static])

if test $ac_cv_prog_gxx = yes; then
  CXXFLAGS="$CXXFLAGS -Wall -Werror -Wextra -Wshadow"
fi

case $CXX in
  clang++)
    CXXFLAGS="$CXXFLAGS -Wno-overloaded-virtual"
  ;;
esac

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
