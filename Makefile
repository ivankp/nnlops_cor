ROOT_CFLAGS := -Wno-deprecated-declarations -pthread -m64
ROOT_CFLAGS += -I$(shell root-config --incdir)
ROOT_LIBS   := $(shell root-config --libs)

.PHONY: all clean

EXE := cor cor2 test_hist test_cov test_with_data cor3

all: cor test_with_data cor3
test: test_hist test_cov test_with_data

cor test_cov test_with_data cor2: include/mat.hh

$(EXE): %: src/%.cc include/timed_counter.hh
	$(CXX) -std=c++14 -Wall -g -Iinclude -fmax-errors=3 \
	  $(ROOT_CFLAGS) \
	  $(filter %.cc,$^) -o $@ \
	  $(ROOT_LIBS) -lTreePlayer

clean:
	@rm -fv $(EXE)

