AC_INIT([HelixHough], [1.00])
AC_CONFIG_SRCDIR([configure.ac])

AM_INIT_AUTOMAKE

AC_PROG_CXX(CC g++)
AC_ENABLE_STATIC(no)

LT_INIT([disable-static])

dnl   no point in suppressing warnings people should 
dnl   at least see them, so here we go for g++: -Wall
dnl if test $ac_cv_prog_gxx = yes; then
dnl    CXXFLAGS="$CXXFLAGS -Wall -Werror"
dnl fi

ROOTLIBS=`$ROOTSYS/bin/root-config --libs`
AC_SUBST(ROOTLIBS)

dnl test for root 6
if test `root-config --version | awk '{print $1>=6.?"1":"0"}'` = 1; then
CINTDEFS=" -noIncludePaths  -inlineInputHeader "
AC_SUBST(CINTDEFS)
fi

AC_CONFIG_FILES([
	Makefile
	FitNewton/Makefile
	Seamstress/Makefile
	helix_hough/Makefile
	helix_hough/circlegen/Makefile
	helix_hough/test/Makefile
])
AC_OUTPUT
