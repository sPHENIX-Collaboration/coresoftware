AC_INIT(PHTpcTracker,[1.00])
AC_CONFIG_SRCDIR([configure.ac])

AM_INIT_AUTOMAKE
AC_PROG_CXX(CC g++)
LT_INIT([disable-static])


if test $ac_cv_prog_gxx = yes; then
   CXXFLAGS="$CXXFLAGS -Wall -Wextra -Werror"
fi

case $CXX in
 clang++)
  CXXFLAGS="$CXXFLAGS -Wno-unused-parameter -Wno-deprecated-copy"
 ;;
esac

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
