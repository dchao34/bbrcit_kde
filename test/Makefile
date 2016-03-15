# Unit tests
CPUBINS := test_kde20
#CPUBINS := $(shell cat cpu_binary_targets.txt)
GPUBINS := test_kde20_cuda
#GPUBINS := $(shell cat gpu_binary_targets.txt)

# additional include and library search paths
IDIR := ../include
LDIR := ../lib
DDIR := .d
INCFLAGS = $(addprefix -I, $(IDIR))
LDFLAGS = $(addprefix -L, $(LDIR)) -lbbrcit_kde

# host compiler
CXX ?= g++
CXXFLAGS = -Wall -Werror -pedantic -pthread -std=c++11 -O2

# device compiler
NVCC := nvcc -ccbin=$(CXX)
NVCCFLAGS = -std=c++11 -O2 -arch=sm_35 -x cu

$(shell mkdir -p $(DDIR) > /dev/null)
CXXDEPFLAGS = -MT $@ -MMD -MP -MF $(DDIR)/$*.Td
CXXPOSTCOMPILE = @mv -f $(DDIR)/$*.Td $(DDIR)/$*.d
NVCCPOSTCOMPILE = @$(NVCC) -M -MT $@ $(NVCCFLAGS) $(INCFLAGS) $< -o $(DDIR)/$@.d

.PHONY : all clean cleanall

all : cpu gpu

cpu : $(CPUBINS)

gpu : $(GPUBINS)

% : %.cc
	$(CXX) $(CXXDEPFLAGS) $(CXXFLAGS) $(INCFLAGS) $(LDFLAGS) $< -o $@
	$(CXXPOSTCOMPILE)

%_cuda : %.cc
	$(NVCC) $(NVCCFLAGS) $(INCFLAGS) $< -o $@
	$(NVCCPOSTCOMPILE)

$(DDIR)/%.d : ;

include $(patsubst %,$(DDIR)/%.d,$(CPUBINS))
include $(patsubst %,$(DDIR)/%.d,$(GPUBINS))

clean:
	@rm -rf *~ *.o $(CPUBINS) $(GPUBINS)

cleanall:
	@rm -rf *~ *.o $(CPUBINS) $(GPUBINS) \
		*.dat *.csv *.out *.dSYM $(DDIR)/*.d
