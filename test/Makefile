# Unit tests
TESTBINS := test_fft generate_bimodal_gauss

# Remainder of this file controls the compilation process

IDIR := ../include
LDIR := ../lib
DDIR := .d

CXX ?= g++
CXXFLAGS = -Wall -Werror -pedantic -std=c++11
INCFLAGS = -I$(IDIR)
LDFLAGS = -L$(LDIR) -lbbrcit_kde

$(shell mkdir -p $(DDIR) > /dev/null)
DEPFLAGS = -MT $@ -MMD -MP -MF $(DDIR)/$*.Td
POSTCOMPILE = @mv -f $(DDIR)/$*.Td $(DDIR)/$*.d
SRCS = $(addsuffix .cc, $(TESTBINS))

.PHONY : all clean

all : $(TESTBINS)

% : %.cc
	$(CXX) $(DEPFLAGS) $(CXXFLAGS) $(INCFLAGS) $(LDFLAGS) $< -o $@
	$(POSTCOMPILE)

$(DDIR)/%.d : ;

include $(patsubst %,$(DDIR)/%.d,$(basename $(SRCS)))

clean:
	@rm -rf *~ *.o $(TESTBINS)
