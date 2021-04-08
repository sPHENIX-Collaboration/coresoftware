AC_INIT(kfparticle_sphenix,[0.01])
AC_CONFIG_SRCDIR([configure.ac])

AM_INIT_AUTOMAKE
AC_PROG_CXX(CC g++)
LT_INIT([disable-static])

dnl KFParticle KFParticleDatabase.h includes Vc/Vc which collides with
dnl clangs inbuild Vc by redefining a macro
case $CXX in
  clang++)
   CXXFLAGS="$CXXFLAGS -Wall -Werror -Wno-macro-redefined -Wno-return-stack-address"
  ;;
  g++)
   CXXFLAGS="$CXXFLAGS -Wall -Werror"
  ;;
esac

CINTDEFS=" -noIncludePaths  -inlineInputHeader "
AC_SUBST(CINTDEFS)

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
