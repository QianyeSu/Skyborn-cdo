
# Process this file with automake to produce Makefile.in
SUBDIRS = libcdi src/lib/ncl src/lib/yac/clapack src/lib/yac/src src/lib/gradsdes  src/lib/healpix src contrib test/data test/pytest
if MAINTAINER_MODE
SUBDIRS += test/executeOnly
endif
#
ACLOCAL_AMFLAGS = -I m4 -I libcdi/m4
#
EXTRA_DIST=config/default LICENSE OPERATORS doc/cdo.pdf doc/cdo_cmor.pdf doc/cdo_eca.pdf doc/cdo_magics.pdf doc/cdo_refcard.pdf cdo.spec README
#
CLEANFILES  = `ls *~ 2> /dev/null`
