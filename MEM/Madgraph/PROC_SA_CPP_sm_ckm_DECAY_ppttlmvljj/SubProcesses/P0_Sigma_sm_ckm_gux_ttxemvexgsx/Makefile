LIBDIR=../../lib
INCDIR=../../src
MODELLIB=model_sm_ckm_ppttlmvljj
CXXFLAGS= -O -I$(INCDIR) -I. -fPIC
LIBFLAGS= -L$(LIBDIR) -l$(MODELLIB)

main=check

all: $(main)

objects=CPPProcess_P0_Sigma_sm_ckm_gux_ttxemvexgsx.o check_sa.o

$(main): $(objects) $(LIBDIR)/lib$(MODELLIB).a
	$(CXX) -o $@ $(objects) $(LIBFLAGS)

$(LIBDIR)/lib$(MODELLIB).a:
	cd ../../src
	make

.PHONY: clean

clean:
	rm -f $(target)
	rm -f $(objects)

