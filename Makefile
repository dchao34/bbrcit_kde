# Object files to package into the library. 
LIBOBJS := Kde1d.o Point.o Kde2d.o ProdKde2d.o ProdAdaKde2d.o fft.o file_io_utils.o

# Library name
LIBNAME := bbrcit_kde

# Remainder of this file controls the compilation process

IDIR := ./include
LDIR := ./lib
SDIR := ./src
ODIR := ./build
DDIR := .d

CXX ?= g++
CXXFLAGS = -Wall -Werror -pedantic -std=c++11
INCFLAGS = -I$(IDIR)
LDFLAGS = -L$(LDIR) -l$(LIBNAME)

$(shell mkdir -p $(DDIR) > /dev/null)
DEPFLAGS = -MT $@ -MMD -MP -MF $(DDIR)/$*.Td
POSTCOMPILE = @mv -f $(DDIR)/$*.Td $(DDIR)/$*.d
SRCS = $(addsuffix .cc, $(basename $(LIBOBJS)))

.PHONY : all lib clean

all : lib

lib : $(addprefix $(ODIR)/, $(LIBOBJS))
	@ar rcs $(LDIR)/lib$(LIBNAME).a $^

$(ODIR)/%.o : $(SDIR)/%.cc $(DDIR)/%.d
	$(CXX) $(DEPFLAGS) $(CXXFLAGS) $(INCFLAGS) -c $< -o $@
	$(POSTCOMPILE)

$(DDIR)/%.d : ;

include $(patsubst %,$(DDIR)/%.d,$(basename $(SRCS)))

clean:
	@rm -rf *~ $(ODIR)/*.o $(LDIR)/*
