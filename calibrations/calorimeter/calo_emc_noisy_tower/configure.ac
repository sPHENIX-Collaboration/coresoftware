AC_INIT(emcNoisyTowerFinder,[1.00])
AC_CONFIG_SRCDIR([configure.ac])

AM_INIT_AUTOMAKE
AC_PROG_CXX(CC g++)
AC_PROG_LIBTOOL

LT_INIT([disable-static])

ROOTLIBS=`root-config --libs`
AC_SUBST(ROOTLIBS)

dnl   no point in suppressing warnings people should 
dnl   at least see them, so here we go for g++: -Wall
if test $ac_cv_prog_gxx = yes; then
   CXXFLAGS="$CXXFLAGS -Wall -Wextra -Wshadow -Werror"
fi

CINTDEFS=" -noIncludePaths -inlineInputHeader "
AC_SUBST(CINTDEFS)

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
