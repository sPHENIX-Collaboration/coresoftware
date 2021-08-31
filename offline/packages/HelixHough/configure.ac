AC_INIT([HelixHough], [1.00])
AC_CONFIG_SRCDIR([configure.ac])

AM_INIT_AUTOMAKE

AC_PROG_CXX(CC g++)
AC_ENABLE_STATIC(no)

LT_INIT([disable-static])

dnl   code produces warnings - helixhough is on its way out, not worth bothering
if test $ac_cv_prog_gxx = yes; then
    CXXFLAGS="$CXXFLAGS -Wno-error"
fi

ROOTLIBS=`$ROOTSYS/bin/root-config --libs`
AC_SUBST(ROOTLIBS)

dnl test for root 6
if test `root-config --version | awk '{print $1>=6.?"1":"0"}'` = 1; then
CINTDEFS=" -noIncludePaths  -inlineInputHeader -Wno-inconsistent-missing-override "
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
