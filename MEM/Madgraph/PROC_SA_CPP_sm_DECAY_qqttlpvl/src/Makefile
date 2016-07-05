LIBDIR=../lib
CXXFLAGS= -O -I. -fPIC

target=$(LIBDIR)/libmodel_sm_qqttlpvl.a

all: $(target)

objects=HelAmps_sm_qqttlpvl.o Parameters_sm_qqttlpvl.o rambo.o ../../PROC_SA_CPP_sm_4/src/read_slha.o

$(target): $(objects)
	$(AR) cru $@ $(objects)
	ranlib $(target)

.PHONY: clean

clean:
	rm -f $(target)
	rm -f $(objects)

