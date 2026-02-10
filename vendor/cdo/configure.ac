# Process this file with autoconf to produce a configure script.
#
#  automake 1.14
#  autoconf 2.69
#  libtool  2.4.2

AC_PREREQ([2.69])
AC_INIT([cdo],[2.5.3],[https://mpimet.mpg.de/cdo])

AC_DEFINE_UNQUOTED(CDO, ["$PACKAGE_VERSION"], [CDO version])

CONFIG_ABORT=yes
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_MACRO_DIRS([m4 libcdi/m4])
AC_CANONICAL_HOST
AC_CANONICAL_BUILD

# AM_INIT_AUTOMAKE([foreign -Wall -Werror])
AM_INIT_AUTOMAKE([subdir-objects foreign -Wall])
AC_CONFIG_SRCDIR([src/cdo.cc])
AC_CONFIG_HEADERS([src/config.h])
AM_MAINTAINER_MODE([disable])

# Set up libtool.
AM_PROG_AR
#
AC_MSG_NOTICE([setting up libtool])
LT_INIT

#python
AM_PATH_PYTHON([3.6], [:])
AM_CONDITIONAL([HAVE_PYTHON], [test "$PYTHON" != :])

# bash
AC_PATH_TOOL([BASH],[bash],[:])

# Check building environment
AC_CHECK_TOOL([CC],[gcc],[:])
AC_CHECK_TOOL([CXX],[g++],[:])
AC_PROG_CC_C99
AM_PROG_CC_C_O
AC_PROG_CXX
AC_PROG_F77
AC_C_RESTRICT
AX_CXX_COMPILE_STDCXX(20)
AC_OPENMP
AC_CHECK_TOOL([AR],[ar],[:])
AC_CHECK_TOOL([CPP],[cpp],[:])
AC_CHECK_TOOL([LD],[ld],[:])
AC_CHECK_TOOL([NM],[nm],[:])
AC_CHECK_TOOL([AR],[ar],[:])
AC_CHECK_TOOL([AS],[as],[:])
AC_CHECK_TOOL([DLLTOOL],[dlltool],[:])
AC_CHECK_TOOL([OBJDUMP],[objdump],[:])
AC_CHECK_TOOL([STRIP],[strip],[:])
AC_CHECK_TOOL([RANLIB],[ranlib],[:])
AC_PROG_INSTALL
#
FC=no
export FC
#
#  ----------------------------------------------------------------------
# Check large file support on 32 bit system
AC_SYS_LARGEFILE
#  ----------------------------------------------------------------------
# Checks for library functions.
AC_FUNC_MMAP
#  ----------------------------------------------------------------------
# Checks for structures.
AC_CHECK_MEMBERS([struct stat.st_blksize])
#  ----------------------------------------------------------------------
# Checks for the availability of functions
AC_CHECK_FUNCS(mallinfo)
AC_CHECK_FUNCS(backtrace)
#  ----------------------------------------------------------------------
# Checks for the availability of ANSI-C99 functions
AC_CHECK_FUNCS(getrlimit)

# Check compiler version
case "$CXX" in
  *pgc*)   CXX_VERSION=`$CXX -V | head -2 | tail -n 1`;;
  *sxc*)   CXX_VERSION=`$CXX -V 2>&1 | tail -n 1`;;
  *xlc*)   CXX_VERSION=`$CXX -qversion 2>&1 | head -n 1`;;
  *)       CXX_VERSION=`$CXX --version 2>&1 | head -n 1 | grep -v error`;;
esac

case "$CC" in
  *pgc*)   C_VERSION=`$CC -V | head -2 | tail -n 1`;;
  *sxc*)   C_VERSION=`$CC -V 2>&1 | tail -n 1`;;
  *xlc*)   C_VERSION=`$CC -qversion 2>&1 | head -n 1`;;
  *)       C_VERSION=`$CC --version 2>&1 | head -n 1 | grep -v error`;;
esac

if test -n "$F77" ; then
case "$F77" in
  *pgf*)      F77_VERSION=`$F77 -V | head -2 | tail -n 1`;;
  *)          F77_VERSION=`$F77 --version 2>&1 | head -n 1 | grep -v error`;;
esac
fi

if test -z "$CXX_VERSION" ; then CXX_VERSION="unknown"; fi;
AC_DEFINE_UNQUOTED(CXX_VERSION, ["$CXX_VERSION"], [CXX Compiler version])

if test -z "$C_VERSION" ; then C_VERSION="unknown"; fi;
AC_DEFINE_UNQUOTED(C_VERSION, ["$C_VERSION"], [C Compiler version])

if test -n "$F77" ; then
if test -z "$F77_VERSION" ; then F77_VERSION="unknown"; fi;
AC_DEFINE_UNQUOTED(F77_VERSION, ["$F77_VERSION"], [F77 Compiler version])
fi

# Check for system type
AC_DEFINE_UNQUOTED([SYSTEM_TYPE],["$ac_cv_build"], [System type])
AC_SUBST([SYSTEM_TYPE],["$ac_cv_build"])
#  ----------------------------------------------------------------------
# check for feenableexcept needs to define macro __USE_GNU and to include fenv.h
# AC_CHECK_FUNCS(feenableexcept)
#
AC_CHECK_MEMBERS([fenv_t.__control, fenv_t.__mxcsr],,,[[#include <fenv.h>]])

#  ----------------------------------------------------------------------
#  Create the Interface to Fortran77 routines via cfortran.h
AC_ARG_ENABLE([fortran],
              [AS_HELP_STRING([--disable-fortran],
                              [Omit building of Fortran routines])],
              [enable_fortran=${enableval}],[enable_fortran=yes])
AS_IF([test "x${enable_fortran}" = "xyes"], [
dnl ######################################################################
dnl Test whether cfortran.h works correctly
dnl ######################################################################
ACX_FIND_CFORTRAN_DEF
ACX_XLF_QEXTNAME_ADD_APPENDUS
AS_IF([test -n "$F77" -a X"$F77" != Xno],
  [ACX_CHECK_CFORTRAN([$srcdir/src],
    [AC_DEFINE([HAVE_CF_INTERFACE],[1],
       [Defined to 1 if C / Fortran interface cfortran.h works])],
    [AS_CASE([x"$acx_cv_cfortran_works"],
       [x"error"],
       [AC_MSG_NOTICE([Linking/Running with C EXTERNAL built with cfortran.h does not work!])],
       [x"compiling with cfortran.h failed"],
       [AC_MSG_NOTICE([Compilation with cfortran.h is not working!])],
       [x"error compiling Fortran subroutine"],
       [AC_MSG_NOTICE([compilation of simple Fortran source failed!])],
       [AC_MSG_NOTICE([Unexpected error when linking C and Fortran via cfortran.h!])])
     AC_MSG_NOTICE([Disabling cfortran.h bindings generation])
     acx_cv_cfortran_works=no])])
])
AS_IF([test x$enable_fortran = 'xno'],[AC_SUBST([ENABLE_FORTRAN],[false])],[AC_SUBST([ENABLE_FORTRAN],[true])])
#
AM_CONDITIONAL([USE_F77],[test -n "$F77" -a X"$F77" != Xno -a x"$acx_cv_cfortran_works" = xyes])
AS_IF([test x$acx_cv_cfortran_works = 'xno'],[AC_SUBST([FORTRAN_WORKS],[no])],[AC_SUBST([FORTRAN_WORKS],[yes])])
#  ----------------------------------------------------------------------
CFLAGS="$CFLAGS ${OPENMP_CFLAGS}"
CXXFLAGS="$CXXFLAGS ${OPENMP_CFLAGS}"
#  ----------------------------------------------------------------------
# Add configure options
ACX_CDO_OPTIONS
ACX_CDI_OPTIONS
#  ----------------------------------------------------------------------
AC_SUBST([CXXFLAGS])
AC_SUBST([CPPFLAGS])
AC_SUBST([FFLAGS])
AC_SUBST([CLDFLAGS])
AC_SUBST([CLIBS])
#AC_SUBST([INCLUDES])
AC_SUBST([AM_CPPFLAGS])
AC_SUBST([AR])

# Checks for compiler
C_COMPILER="$CC $CFLAGS"
AC_DEFINE_UNQUOTED(C_COMPILER, ["$C_COMPILER"], [C Compiler])
CXX_COMPILER="$CXX $CXXFLAGS"
AC_DEFINE_UNQUOTED(CXX_COMPILER, ["$CXX_COMPILER"], [CXX Compiler])
F77_COMPILER="$F77 $FFLAGS"
AC_DEFINE_UNQUOTED(F77_COMPILER, ["$F77_COMPILER"], [F77 Compiler])

AC_REQUIRE_AUX_FILE([tap-driver.sh])
AC_PROG_AWK
#test util
AC_CONFIG_FILES([test/cdoTestFunctions.test ],[chmod a+x "$ac_file"])

AC_CONFIG_FILES([test/pytest/cdoTest.py
                test/pytest/CDO_test.py.test
                test/pytest/Adisit.py.test
                test/pytest/Afterburner.py.test
                test/pytest/Arith.py.test
                test/pytest/Arithc.py.test
                test/pytest/Arith_extra.py.test
                test/pytest/Cat.py.test
                test/pytest/Change.py.test
                test/pytest/CMOR.py.test
                test/pytest/Collgrid.py.test
                test/pytest/Comp.py.test
                test/pytest/Compc.py.test
                test/pytest/Cond.py.test
                test/pytest/Cond2.py.test
                test/pytest/Condc.py.test
                test/pytest/Consecstat.py.test
                test/pytest/Copy_netcdf.py.test
                test/pytest/Dayarith.py.test
                test/pytest/Detrend.py.test
                test/pytest/Eca.py.test
                test/pytest/Enspctl.py.test
                test/pytest/Ensstat.py.test
                test/pytest/Ensstat2.py.test
                test/pytest/EOFcoeff.py.test
                test/pytest/EOF.py.test
                test/pytest/Etccdi.py.test
                test/pytest/Etccdi2.py.test
                test/pytest/Expr.py.test
                test/pytest/File.py.test
                test/pytest/Filter.py.test
                test/pytest/Fldpctl.py.test
                test/pytest/Fldstat2.py.test
                test/pytest/Fldstat.py.test
                test/pytest/Genweights.py.test
                test/pytest/Gradsdes.py.test
                test/pytest/Gridarea.py.test
                test/pytest/Gridboxstat.py.test
                test/pytest/Importcmsaf.py.test
                test/pytest/Intgrid.py.test
                test/pytest/Inttime.py.test
                test/pytest/Intyear.py.test
                test/pytest/Isosurface.py.test
                test/pytest/Maggraph.py.test
                test/pytest/Magplot.py.test
                test/pytest/Magvector.py.test
                test/pytest/MapReduce.py.test
                test/pytest/Maskregion.py.test
                test/pytest/Math.py.test
                test/pytest/Merge.py.test
                test/pytest/Mergetime.py.test
                test/pytest/Merstat.py.test
                test/pytest/Monarith.py.test
                test/pytest/Multiyearstat.py.test
                test/pytest/Ninfo.py.test
                test/pytest/Pack.py.test
                test/pytest/Percentile.py.test
                test/pytest/Read_grib.py.test
                test/pytest/Read_netcdf.py.test
                test/pytest/Remap.py.test
                test/pytest/Remap_noweights.py.test
                test/pytest/Remap2.py.test
                test/pytest/Remap3.py.test
                test/pytest/Remap4.py.test
                test/pytest/Remap_global_5_grid.py.test
                test/pytest/Remap_gme.py.test
                test/pytest/Remap_healpix.py.test
                test/pytest/Remap_small.py.test
                test/pytest/Remap_knn.py.test
                test/pytest/Remap_extra.py.test
                test/pytest/Remap_extra_file.py.test
                test/pytest/Remap_extra_file2.py.test
                test/pytest/Remapeta.py.test
                test/pytest/Remapstat.py.test
                test/pytest/Runpctl.py.test
                test/pytest/Seasstat.py.test
                test/pytest/Runstat.py.test
                test/pytest/Select.py.test
                test/pytest/Selregion.py.test
                test/pytest/Setmiss.py.test
                test/pytest/Smooth.py.test
                test/pytest/Spectral.py.test
                test/pytest/Split.py.test
                test/pytest/Splittime.py.test
                test/pytest/threads.py.test
                test/pytest/Timfillmiss.py.test
                test/pytest/Timpctl.py.test
                test/pytest/Timselpctl.py.test
                test/pytest/Timselstat.py.test
                test/pytest/Timstat2.py.test
                test/pytest/Timstat3.py.test
                test/pytest/Timstat.py.test
                test/pytest/tsformat.py.test
                test/pytest/userInput.py.test
                test/pytest/Varsstat.py.test
                test/pytest/Vertfillmiss.py.test
                test/pytest/Vertint.py.test
                test/pytest/Vertstat.py.test
                test/pytest/wildcard.py.test
                test/pytest/cdoReturnValues.py.test
                test/pytest/Wind.py.test
                test/pytest/Ydayarith.py.test
                test/pytest/Ydrunpctl.py.test
                test/pytest/Ydrunstat.py.test
                test/pytest/Yeararith.py.test
                test/pytest/Yearmonstat.py.test
                test/pytest/Ymonarith.py.test
                test/pytest/Zonstat.py.test],[chmod a+x "$ac_file"])

#internal tests
AC_CONFIG_FILES([test/pytest/Vargen.py.test ],[chmod a+x "$ac_file"])

#tests which only execute operators and dont check about correctness
AC_CONFIG_FILES([test/executeOnly/Splitsel.test], [chmod a+x "$ac_file"])

#Cdo Options
AC_CONFIG_FILES([contrib/Makefile  cdo.spec cdo.settings])
#Makefiles
AC_CONFIG_FILES([Makefile src/lib/ncl/Makefile src/lib/yac/src/Makefile src/lib/gradsdes/Makefile src/lib/healpix/Makefile src/Makefile])
AC_CONFIG_FILES([src/lib/yac/clapack/Makefile])
#Test Makefiles
AC_CONFIG_FILES([test/data/Makefile test/executeOnly/Makefile test/pytest/Makefile])
AC_OUTPUT

#  ----------------------------------------------------------------------
#  Show configuration
AC_MSG_NOTICE([CDO is configured with the following options:])
cat cdo.settings

AS_IF([test "$ac_cv_prog_cc_c99" = "no"],
      [AC_MSG_NOTICE([

  Warning: The C compiler does not accept ANSI C99 source code!
])])
AC_MSG_NOTICE([

  Configuration completed.

    You can now say 'make' to compile the CDO package and 'make install' to install it afterwards.
])
