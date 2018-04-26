# This is where the version information is set
AC_INIT([CRAFT_HEIMDALL], [0.1.0], [andrew.jameson@gmail.com])

AC_CONFIG_MACRO_DIR(config)
AC_CONFIG_SRCDIR([src/craft_utils.c])
AC_PROG_LIBTOOL
AC_DISABLE_SHARED

# Automake initialization
AM_INIT_AUTOMAKE
AM_CONFIG_HEADER([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

# Checks for libraries.
SWIN_LIB_PSRDADA
SWIN_LIB_SOFA

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
AC_CONFIG_FILES([Makefile
  src/Makefile  
  sigproc/Makefile
])

AC_OUTPUT

