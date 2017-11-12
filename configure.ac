# Prelude
AC_INIT([waterwaves], [1.0], [michiel.eggermont@gmail.com])
AC_PREREQ([2.69])
AM_INIT_AUTOMAKE([foreign -Wall -Werror])

# Checks for programs
AC_PROG_F77

# Checks for libraries via pkg-config
PKG_CHECK_MODULES(blas, [blas])
PKG_CHECK_MODULES(lapack, [lapack])

AC_F77_FUNC(idamax)
AC_CHECK_LIB([blas], [$idamax])

AC_F77_FUNC(dptsv)
AC_CHECK_LIB([lapack], [$dptsv])

AC_F77_FUNC(nf_create)
AC_CHECK_LIB([netcdff], [$nf_create])

# Check for header files
AC_CONFIG_HEADERS([config.h])

# Output files
AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
