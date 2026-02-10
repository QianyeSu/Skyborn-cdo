EXTRA_DIST = makecompl.rb cdoCompletion.bash cdoCompletion.tcsh cdoCompletion.zsh
CLEANFILES =
DISTCLEANFILES = `ls cdoCompletion.*`

completions:
	if hash ruby >/dev/null 2>&1 ; then \
	  ruby -KN $(srcdir)/makecompl.rb -o cdoCompletion; \
	fi

cdoCompletion.bash cdoCompletion.zsh cdoCompletion.tcsh: completions

test:
	if hash ruby >/dev/null 2>&1 ; then \
	    cd ruby;ruby test/test_cdo.rb \
	fi
	if hash python >/dev/null 2>&1 ; then \
	    cd python; python test/test_cdo.py \
	fi


if MAINTAINER_MODE
all-local: completions

CLEANFILES += $(DISTCLEANFILES)
endif
