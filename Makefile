CXX ?= g++
CXXFLAGS = -Wall -Werror -pedantic -std=c++11

IDIR = ./include
LDIR = ./lib
SDIR = ./src
ODIR = ./build
TDIR = ./test

INCFLAGS = -I$(IDIR)
LDFLAGS = -L$(LDIR) -lbbrcit_kde

# Object files to package into the library. 
LIBOBJS = Kde2d.o ProdKde2d.o ProdAdaKde2d.o fft.o file_io_utils.o

# Unit tests
TESTBINS = test_fft generate_bimodal_gauss

.PHONY : all lib test clean

all : lib test

lib : $(addprefix $(ODIR)/, $(LIBOBJS))
	ar rcs $(LDIR)/libbbrcit_kde.a $^

test : lib $(addprefix $(TDIR)/, $(TESTBINS))

$(TDIR)/% : $(TDIR)/%.cc lib
	$(CXX) $(CXXFLAGS) $(INCFLAGS) $(LDFLAGS) $< -o $@

$(ODIR)/%.o : $(SDIR)/%.cc $(IDIR)/%.h
	$(CXX) $(CXXFLAGS) $(INCFLAGS) -c $< -o $@

clean:
	@rm -f *~ $(ODIR)/*.o $(LDIR)/* $(addprefix $(TDIR)/, $(TESTBINS))
