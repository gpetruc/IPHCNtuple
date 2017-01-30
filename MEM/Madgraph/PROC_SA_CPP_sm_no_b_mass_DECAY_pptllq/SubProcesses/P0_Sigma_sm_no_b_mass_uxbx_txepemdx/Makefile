LIBDIR=../../lib
INCDIR=../../src
MODELLIB=model_sm_no_b_mass_pptllq
CXXFLAGS= -O -I$(INCDIR) -I. -fPIC
LIBFLAGS= -L$(LIBDIR) -l$(MODELLIB)

main=check

all: $(main)

objects=CPPProcess_P0_Sigma_sm_no_b_mass_uxbx_txepemdx.o check_sa.o

$(main): $(objects) $(LIBDIR)/lib$(MODELLIB).a
	$(CXX) -o $@ $(objects) $(LIBFLAGS)

$(LIBDIR)/lib$(MODELLIB).a:
	cd ../../src
	make

.PHONY: clean

clean:
	rm -f $(target)
	rm -f $(objects)

