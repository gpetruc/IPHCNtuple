LIBDIR=../../lib
INCDIR=../../src
MODELLIB=model_sm_hw2l2nu
CXXFLAGS= -O -I$(INCDIR) -I. -fPIC
LIBFLAGS= -L$(LIBDIR) -l$(MODELLIB)

main=check

all: $(main)

objects=CPPProcess_hw2l2nu.o check_sa.o

$(main): $(objects) $(LIBDIR)/lib$(MODELLIB).a
	$(CXX) -o $@ $(objects) $(LIBFLAGS)

$(LIBDIR)/lib$(MODELLIB).a:
	cd ../../src
	make

.PHONY: clean

clean:
	rm -f $(target)
	rm -f $(objects)

