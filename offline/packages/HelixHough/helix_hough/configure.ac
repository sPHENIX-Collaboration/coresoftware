AC_INIT(HelixHough, 1.0.0, alan.dion@gmail.com)

AC_CONFIG_AUX_DIR(config)
AC_CONFIG_SRCDIR(HelixHough_init.cpp)

AM_INIT_AUTOMAKE

AC_PROG_CXX
AC_PROG_LIBTOOL

PKG_CHECK_MODULES([eigen3], [eigen3])
PKG_CHECK_MODULES([fitnewton], [fitnewton])
PKG_CHECK_MODULES([Seamstress], [seamstress])

AC_SUBST(sPHENIX_INC)

AC_ARG_ENABLE(
  sPHENIX,
  AS_HELP_STRING(
    [--enable-sPHENIX],
    [build the sPHENIX tracker]
  ),
  [sPHENIX_INC="sPHENIX/sPHENIXTracker.h"
   ac_cv_use_sphenix=$enableval
  ],
  [sPHENIX_INC=''
   ac_cv_use_sphenix=no
  ]
)

AM_CONDITIONAL(SPHENIX_ENABLED, test x"$ac_cv_use_sphenix" = x"yes")


AC_SUBST(AVX_FLAG)

AC_ARG_ENABLE(
  avx,
  AS_HELP_STRING(
    [--enable-avx],
    [use avx instructions]
  ),
  [AVX_FLAG="-DAVXHOUGH"
  ac_cv_use_avx=$enableval
  ],
  [AVX_FLAG=''
  ac_cv_use_avx=no
  ]
)

CXXFLAGS="-Werror -O2 -DNDEBUG  ${CXXFLAGS} ${AVX_FLAG}"

dnl test for root 6
if test `root-config --version | gawk '{print $1>=6.?"1":"0"}'` = 1; then
CINTDEFS=" -noIncludePaths  -inlineInputHeader "
AC_SUBST(CINTDEFS)
fi

AC_CONFIG_FILES([Makefile helix_hough.pc])

AC_OUTPUT

