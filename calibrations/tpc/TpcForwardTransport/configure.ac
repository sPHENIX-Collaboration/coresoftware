AC_INIT(TpcForwardTransport, [1.00])
AC_CONFIG_SRCDIR([configure.ac])
AM_INIT_AUTOMAKE
AC_PROG_CXX(CC g++)
LT_INIT([disable-static])

CXXFLAGS="$CXXFLAGS -std=c++20 -Wall -Wextra -Wshadow -Werror"
AC_LANG_PUSH([C++])
AC_MSG_CHECKING([for C++20 standard library support])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[#include <bit>
#include <cstdint>
#include <numbers>]], [[static_assert(__cplusplus >= 202002L);
  const auto bits = std::bit_cast<std::uint64_t>(std::numbers::pi);
  return bits == 0;]])],
  [AC_MSG_RESULT([yes])],
  [AC_MSG_RESULT([no])
   AC_MSG_ERROR([TpcForwardTransport requires a compiler and standard library supporting C++20])])
AC_LANG_POP([C++])
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
