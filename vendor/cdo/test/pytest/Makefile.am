export

# COLOR ALWAYS also if no color capable terminal detected
#     this is needed i.e. to run the tests via a oneline ssh cmd:
#     ssh host '... make check'
AM_COLOR_TESTS=always

TEST_LOG_DRIVER = env AM_TAP_AWK='$(AWK)' $(SHELL) \
                  $(top_srcdir)/config/tap-driver.sh

TESTS =\
	    cdoTest.py\
        CDO_test.py.test\
		Adisit.py.test\
		Afterburner.py.test\
		Arith.py.test\
		Arithc.py.test\
        Arith_extra.py.test\
        Cat.py.test\
	    Change.py.test\
		CMOR.py.test\
		Collgrid.py.test\
		Comp.py.test\
		Compc.py.test\
		Cond.py.test\
		Cond2.py.test\
		Condc.py.test\
		Consecstat.py.test\
		Copy_netcdf.py.test\
		Dayarith.py.test\
		Detrend.py.test\
		EOF.py.test\
		EOFcoeff.py.test\
		Eca.py.test\
		Enspctl.py.test\
		Ensstat.py.test\
		Ensstat2.py.test\
		Etccdi.py.test\
		Etccdi2.py.test\
		Expr.py.test\
		Extra.py.test\
		File.py.test\
		Filter.py.test\
		Fldpctl.py.test\
		Fldstat.py.test\
		Fldstat2.py.test\
		Genweights.py.test\
		Gradsdes.py.test\
		Gridarea.py.test\
		Gridboxstat.py.test\
        Importcmsaf.py.test\
		Intgrid.py.test\
		Inttime.py.test\
		Intyear.py.test\
		Isosurface.py.test\
		Maggraph.py.test\
		Magplot.py.test\
		Magvector.py.test\
		MapReduce.py.test\
		Maskregion.py.test\
		Math.py.test\
		Merge.py.test\
		Mergetime.py.test\
		Merstat.py.test\
		Merstat.py.test\
		Monarith.py.test\
		Multiyearstat.py.test\
		Ninfo.py.test\
		Pack.py.test\
		Percentile.py.test\
		Read_grib.py.test\
		Read_netcdf.py.test\
		Remap.py.test\
		Remap_noweights.py.test\
		Remap2.py.test\
		Remap3.py.test\
		Remap4.py.test\
		Remap_global_5_grid.py.test\
		Remap_gme.py.test\
		Remap_healpix.py.test\
		Remap_small.py.test\
		Remap_knn.py.test\
		Remap_extra.py.test\
		Remap_extra_file.py.test\
		Remap_extra_file2.py.test\
		Remapeta.py.test\
		Remapstat.py.test\
		Runpctl.py.test\
		Runstat.py.test\
		Seasstat.py.test\
		Select.py.test\
		Selregion.py.test\
		Setmiss.py.test\
		Smooth.py.test\
		Spectral.py.test\
		Split.py.test\
		Splittime.py.test\
		Timcumsum.py.test\
		Timfillmiss.py.test\
		Timpctl.py.test\
		Timselpctl.py.test\
		Timselstat.py.test\
		Timstat.py.test\
		Timstat2.py.test\
		Timstat3.py.test\
		Varsstat.py.test\
		Vertfillmiss.py.test\
		Vertint.py.test\
		Vertstat.py.test\
		Wind.py.test\
		Ydayarith.py.test\
		Ydrunpctl.py.test\
		Ydrunstat.py.test\
		Yeararith.py.test\
		Yearmonstat.py.test\
		Ymonarith.py.test\
		Zonstat.py.test\
		threads.py.test\
		tsformat.py.test\
		userInput.py.test\
		wildcard.py.test\
	    cdoReturnValues.py.test

INTERNAL_TESTS=Vargen.py.test
TESTS+=$(INTERNAL_TESTS)


CDO        = $(top_builddir)/src/cdo
DATAPATH   = $(top_srcdir)/test/data

CLEANFILES = `ls  cerr* cout* ifile* ofile* __pycache__/*`
