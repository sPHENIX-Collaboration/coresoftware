AC_INIT(kfparticle_sphenix,[0.01])
AC_CONFIG_SRCDIR([configure.ac])

AM_INIT_AUTOMAKE
AC_PROG_CXX(CC g++)
LT_INIT([disable-static])

if test $ac_cv_prog_gxx = yes; then
   CXXFLAGS="$CXXFLAGS -Wall -Wextra -Werror"
fi

dnl KFParticle KFParticleDatabase.h includes Vc/Vc which collides with
dnl clangs inbuild Vc by redefining a macro
case $CXX in
  clang++)
   CXXFLAGS="$CXXFLAGS -Wno-macro-redefined -Wno-unused-parameter -Wno-deprecated-copy -Wno-bitwise-instead-of-logical -Wno-null-pointer-subtraction -Wno-unknown-warning-option"
  ;;
esac

CINTDEFS=" -noIncludePaths  -inlineInputHeader "
AC_SUBST(CINTDEFS)

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
