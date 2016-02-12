# Unit tests
TESTBINS := test_kde6 test_kde5 #test_kde1 test_kde2 test_kde3 test_kde4 #test_kdtree3 test_kdtree2 test_kdtree1 test_kernels #test_rectangle #test_kde_attributes #test_interval #test_float_utils test_point_weight_attributes test_decorated_point #test_point  #test_fft generate_bimodal_gauss prodkde2d_scan prodkde2d_fftcv prodadakde2d_scan

# Remainder of this file controls the compilation process

IDIR := ../include
LDIR := ../lib
DDIR := .d

CXX ?= g++
CXXFLAGS = -Wall -Werror -pedantic -std=c++11 -O2
INCFLAGS = -I$(IDIR)
LDFLAGS = -L$(LDIR) -lbbrcit_kde

$(shell mkdir -p $(DDIR) > /dev/null)
DEPFLAGS = -MT $@ -MMD -MP -MF $(DDIR)/$*.Td
POSTCOMPILE = @mv -f $(DDIR)/$*.Td $(DDIR)/$*.d
SRCS = $(addsuffix .cc, $(TESTBINS))

.PHONY : all clean cleanall

all : $(TESTBINS)

% : %.cc
	$(CXX) $(DEPFLAGS) $(CXXFLAGS) $(INCFLAGS) $(LDFLAGS) $< -o $@
	$(POSTCOMPILE)

$(DDIR)/%.d : ;

include $(patsubst %,$(DDIR)/%.d,$(basename $(SRCS)))

clean:
	@rm -rf *~ *.o $(TESTBINS)

cleanall:
	@rm -rf *~ *.o $(TESTBINS) *.dat *.csv *.out *.dSYM
