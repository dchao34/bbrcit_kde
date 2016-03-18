BINARIES = grid_search

OBJECTS = cv_utils.o custom_program_option_utils.o

SRCS = $(wildcard *.cc *.cpp)

CUDA_ROOT = /usr/local/cuda-7.5
CUDA_LIBS = $(CUDA_ROOT)/lib64

BOOST_ROOT = /home/dchao/workspace/boost_1_59_0
BOOST_LIBS = $(BOOST_ROOT)/stage/lib
LIBS = -lboost_program_options

BBRCITKDE_ROOT = /home/dchao/workspace/analyses/2016/March11/kde_cross_validation/bbrcit_kde
INCFLAGS = -I$(BBRCITKDE_ROOT)/include -I$(BOOST_ROOT)
LDFLAGS = -L$(CUDA_LIBS) -lcudart \
          -L$(BOOST_LIBS) \
					-Wl,-rpath,$(BOOST_LIBS) $(LIBS)

BUILDDIR = build

CXX := g++
CXXFLAGS = -Wall -Werror -pedantic -pthread -std=c++11 -O2

DEPDIR = .d
$(shell mkdir -p $(DEPDIR) > /dev/null)
DEPFLAGS = -MT $@ -MMD -MP -MF $(DEPDIR)/$*.Td
POSTCOMPILE = @mv -f $(DEPDIR)/$*.Td $(DEPDIR)/$*.d

NVCC := nvcc -ccbin=$(CXX)
NVCCFLAGS = -std=c++11 -O2 -arch=sm_35 -x cu
NVCCMAKEDEP = @$(NVCC) -M -MT $@ $(NVCCFLAGS) $(INCFLAGS) $< -o\
							$(patsubst $(BUILDDIR)/%.o, $(DEPDIR)/%.d, $@)


all : $(BINARIES)

grid_search : $(addprefix $(BUILDDIR)/, grid_search.o $(OBJECTS))
	$(CXX) $(LDFLAGS) $^ -o $@

$(BUILDDIR)/%.o : %.cpp
$(BUILDDIR)/%.o : %.cpp $(DEPDIR)/%.d
	$(CXX) $(DEPFLAGS) $(CXXFLAGS) -c $< -o $@
	$(POSTCOMPILE)

$(BUILDDIR)/%.o : %.cc
$(BUILDDIR)/%.o : %.cc $(DEPDIR)/%.d
	$(NVCCMAKEDEP)
	$(NVCC) $(NVCCFLAGS) $(INCFLAGS) -c $< -o $@

$(DEPDIR)/%.d: ;

.PRECIOUS: $(DEPDIR)/%.d

-include $(patsubst %,$(DEPDIR)/%.d,$(basename $(SRCS)))

clean : 
	@rm -f *~ $(BINARIES) $(BUILDDIR)/*

cleanall : clean
	@rm -f $(DEPDIR)/* *.csv
