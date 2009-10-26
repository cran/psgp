dnl ------------------------------------------
dnl Configure script for R package PSGP (Projected 
dnl Sequential Gaussian Process)
dnl 
dnl 2009, Remi Barillec <r.barillec@aston.ac.uk>
dnl
dnl ------------------------------------------
AC_INIT([PSGP], 0.2.7)

dnl ------------------------------------------
dnl Option for a custom itpp-config path
dnl ------------------------------------------
ITPP_CONFIG=""
AC_ARG_WITH([itpp-config],
            AC_HELP_STRING([--with-itpp-config=DIR],
                           [specify the path to the itpp-config script]),
            [itpp_config_path=$withval])            

dnl If no custom path specify, assume itpp-config is in the path
if test -z "$itpp_config_path"; then
  ITPP_CONFIG="itpp-config"
else
  ITPP_CONFIG="$itpp_config_path/itpp-config"
fi

dnl ------------------------------------------
dnl Check that itpp-config is found
dnl ------------------------------------------
AC_CHECK_PROG(ITPP_CONFIG_FOUND, $ITPP_CONFIG, "found", "not found") 
if test "$ITPP_CONFIG_FOUND" = "not found"; then
    AC_MSG_ERROR("itpp-config script not found. Make sure you have installed the IT++ library. If the itpp-config script is not in your path, please use "R CMD INSTALL --configure-args='--with-itpp-config=DIR'" where DIR is the folder containing itpp-config. For further information, see installation instructions in the doc/ folder.)
    exit 1
fi

dnl ------------------------------------------
dnl Set itpp flags and libraries
dnl ------------------------------------------
ITPP_CFLAGS=`${ITPP_CONFIG} --cflags`
ITPP_LIBS=`${ITPP_CONFIG} --libs`

dnl ------------------------------------------
dnl Check that R is found
dnl ------------------------------------------
dnl Test for R
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "Could not determine R_HOME. Is R installed?"
  exit 1
fi

dnl ------------------------------------------
dnl C++ compiler
dnl ------------------------------------------
CPP=`"${R_HOME}/bin/R" CMD config CPP`

dnl ------------------------------------------
dnl C++ flags and libraries
dnl ------------------------------------------
PSGP_CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`
PSGP_CPPFLAGS="$PSGP_CPPFLAGS $ITPP_CFLAGS -I."

LIBS=$ITPP_LIBS

dnl ------------------------------------------
dnl Test C++ compiler
dnl ------------------------------------------
AC_PROG_CPP

dnl ------------------------------------------
dnl Substitude PSGP_CPPFLAGS and LIBS
dnl ------------------------------------------
AC_SUBST(PSGP_CPPFLAGS)
AC_SUBST(LIBS)
AC_CONFIG_HEADERS([src/config.h])

dnl ------------------------------------------
dnl Substitute into the Makevars.in and src/config.h
dnl ------------------------------------------
AC_CONFIG_FILES([src/Makevars])

dnl ------------------------------------------
dnl Create Makefile
dnl ------------------------------------------
AC_OUTPUT
