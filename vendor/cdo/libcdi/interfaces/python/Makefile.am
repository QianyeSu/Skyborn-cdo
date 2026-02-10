if ENABLE_PYTHON_INTERFACE

python_PYTHON = Cdi.py
pyexec_LTLIBRARIES = _Cdi.la
_Cdi_la_SOURCES = cdi_wrapper.cpp
_Cdi_la_CPPFLAGS = $(AM_CPPFLAGS) $(PYTHON_CPPFLAGS)
AM_CXXFLAGS = -shared
AM_LDFLAGS = -module -avoid-version
AM_CPPFLAGS = -I$(srcdir)/.. -I$(top_srcdir)/src
_Cdi_la_LIBADD = ../libcdipp.la $(PYTHON_LIBS)

../libcdipp.la:
	$(am__cd) .. && $(MAKE) libcdipp.la

endif ENABLE_PYTHON_INTERFACE

if MAINTAINER_MODE
Cdi.py: ../cdi.i ../cdi.hpp
	$(AM_V_GEN)$(SWIG) -python -c++ -o cdi_wrapper.cpp $<

cdi_wrapper.cpp: Cdi.py
	$(AM_V_at)test -f $@ || rm -f $<
	$(AM_V_at)test -f $@ || $(MAKE) $(AM_MAKEFLAGS) $<
endif
