CC ?= gcc
CXX ?= g++

CXXFLAGS = -Wall -std=c++11

INCFLAGS = -I../
INCFLAGS += -I../data
OBJ = gauss_legendre.o fft.o Kde2d.o ProdKde2d.o data/file_io_utils.o
OBJECTS = $(addprefix ../, $(OBJ))

all : test_gauss test_kde test_cv test_fcv test_fft

test_cv : test_cv.cc
	$(CXX) $(CXXFLAGS) $(INCFLAGS) $^ $(OBJECTS) -o $@

test_fcv : test_fcv.cc
	$(CXX) $(CXXFLAGS) $(INCFLAGS) $^ $(OBJECTS) -o $@

test_fft : test_fft.cc
	$(CXX) $(CXXFLAGS) $(INCFLAGS) $^ $(OBJECTS) -o $@

test_gauss : test_gauss.cc
	$(CXX) $(CXXFLAGS) $(INCFLAGS) $^ $(OBJECTS) -o $@

test_kde : test_kde.cc ../Kde2d.h ../ProdKde2d.h
	$(CXX) $(CXXFLAGS) $(INCFLAGS) $< $(OBJECTS) -o $@

clean:
	@rm -f *~ *.o `find . -perm +100 -type f`
