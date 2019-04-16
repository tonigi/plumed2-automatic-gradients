USE=core tools colvar

STANMATH=math
ADDCPPFLAGS += -I$(STANMATH) -I$(STANMATH)/lib/eigen_3.3.3/ \
	-I$(STANMATH)/lib/cvodes_2.9.0/include  \
	-I$(STANMATH)/lib/boost_1.62.0/

# This is mostly for macports, which extracts but does not do submodule
ifeq ($(shell test -e $(STANMATH)/stan && echo -n yes),yes)
    $(info STAN submodule already initialized)
else
    $(warning Initializing STAN git submodule)
    TMP := $(shell git submodule update --init --recursive)
endif


# generic makefile
include ../maketools/make.module
