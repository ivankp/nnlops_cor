ROOT_CFLAGS := -Wno-deprecated-declarations -pthread -m64
ROOT_CFLAGS += -I$(shell root-config --incdir)
ROOT_LIBS   := $(shell root-config --libs)

.PHONY: all clean

EXE := cov test_hist test_cov

all: cov
test: test_hist test_cov

cov test_cov: include/mat.hh

$(EXE): %: src/%.cc include/timed_counter.hh
	$(CXX) -std=c++14 -Wall -O3 -flto -Iinclude -fmax-errors=3 \
	  $(ROOT_CFLAGS) \
	  $(filter %.cc,$^) -o $@ \
	  $(ROOT_LIBS) -lTreePlayer

clean:
	@rm -fv $(EXE)

