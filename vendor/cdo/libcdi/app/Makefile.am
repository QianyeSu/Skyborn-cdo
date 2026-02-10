bin_PROGRAMS    =
noinst_PROGRAMS =
check_PROGRAMS  =

if ENABLE_CDI_APP
if ENABLE_CDI_LIB
bin_PROGRAMS += cdi
else !ENABLE_CDI_LIB
noinst_PROGRAMS += cdi
endif !ENABLE_CDI_LIB
noinst_PROGRAMS += createtable
else !ENABLE_CDI_APP
if with_on_demand_check_programs
check_PROGRAMS += cdi
else !with_on_demand_check_programs
noinst_PROGRAMS += cdi
endif !with_on_demand_check_programs
endif !ENABLE_CDI_APP

AM_CPPFLAGS = -I$(top_srcdir)/src

AM_LDFLAGS =
if ENABLE_ALL_STATIC
AM_LDFLAGS += -all-static
endif

LDADD = $(top_builddir)/src/libcdi.la

cdi_SOURCES = \
  cdi.c \
  printinfo.c \
  printinfo.h

createtable_SOURCES = createtable.c

clean-local:
	-rm -rf *.dSYM

include $(top_srcdir)/config/lt_overrides_Makefile.inc
