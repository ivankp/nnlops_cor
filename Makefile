ROOT_CFLAGS := -Wno-deprecated-declarations -pthread -m64
ROOT_CFLAGS += -I$(shell root-config --incdir)
ROOT_LIBS   := $(shell root-config --libs)

.PHONY: all clean

EXE := cor test_hist test_cov test_with_data

all: cor test_with_data
test: test_hist test_cov test_with_data

cov test_cov test_with_data: include/mat.hh

$(EXE): %: src/%.cc include/timed_counter.hh
	$(CXX) -std=c++14 -Wall -g -Iinclude -fmax-errors=3 \
	  $(ROOT_CFLAGS) \
	  $(filter %.cc,$^) -o $@ \
	  $(ROOT_LIBS) -lTreePlayer

clean:
	@rm -fv $(EXE)

