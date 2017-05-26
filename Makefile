ROOT_CFLAGS := -Wno-deprecated-declarations -pthread -m64
ROOT_CFLAGS += -I$(shell root-config --incdir)
ROOT_LIBS   := $(shell root-config --libs)

.PHONY: all clean

EXE := cor convert test_hist test_cov

BINNER_DEP := axis.hh binner.hh default_bin_filler.hh exception.hh \
	      sequence_traits.hh type_traits.hh utility.hh

all: cor convert
tests: test_hist test_cov

cor: include/mat.hh include/mapper.hh include/catstr.hh \
     include/timed_counter.hh $(BINNER_DEP:%=include/%)
convert: include/catstr.hh
test_cov: include/mat.hh include/timed_counter.hh
test_hist: include/timed_counter.hh

$(EXE): %: src/%.cc
	$(CXX) -std=c++14 -Wall -O3 -flto -Iinclude -fmax-errors=3 \
	  $(ROOT_CFLAGS) \
	  $(filter %.cc,$^) -o $@ \
	  $(ROOT_LIBS) -lTreePlayer

clean:
	@rm -fv $(EXE)

