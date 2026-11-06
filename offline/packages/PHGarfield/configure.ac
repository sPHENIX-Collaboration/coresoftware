AC_INIT(PHGarfield, [1.00])
AC_CONFIG_SRCDIR([configure.ac])

AM_INIT_AUTOMAKE

AC_PROG_CXX(CC g++)
LT_INIT([disable-static])

AC_LANG_PUSH([C++])
AC_CHECK_HEADERS([openssl/evp.h], [],
  [AC_MSG_ERROR([OpenSSL development headers are required for compact-map SHA256])])
AC_CHECK_LIB([crypto], [EVP_DigestInit_ex], [CRYPTO_LIBS="-lcrypto"],
  [AC_MSG_ERROR([OpenSSL libcrypto is required for compact-map SHA256])])
AC_LANG_POP([C++])
AC_SUBST([CRYPTO_LIBS])

dnl leaving this here in case we want to play with different compiler 
dnl specific flags
CXXFLAGS="$CXXFLAGS -Wall -Wextra -Wshadow -Werror"

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
