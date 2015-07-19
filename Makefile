CC ?= gcc
CXX ?= g++

CXXFLAGS = -Wall -std=c++11

UTILDIR = data

INCDIR = $(UTILDIR)
LIBDIR = $(UTILDIR)

all : subsample make_contours prepare_kde_data test_gauss test_kde cv

cv : cv.cc gauss_legendre.o Kde2d.o ProdKde2d.o $(UTILDIR)/file_io_utils.o
	$(CXX) $(CXXFLAGS) $^ -o $@

test_gauss : test_gauss.cc gauss_legendre.o Kde2d.o ProdKde2d.o $(UTILDIR)/file_io_utils.o
	$(CXX) $(CXXFLAGS) $^ -o $@

test_kde : test_kde.o Kde2d.o gauss_legendre.o ProdKde2d.o $(UTILDIR)/file_io_utils.o
	$(CXX) $(CXXFLAGS) $^ -o $@

test_kde.o : test_kde.cc Kde2d.h ProdKde2d.h
	$(CXX) $(CXXFLAGS) -I$(INCDIR) -c $< -o $@

gauss_legendre.o : gauss_legendre.c gauss_legendre.h
	$(CC) -Wall -c $< -o $@

make_contours : make_contours.o Kde2d.o gauss_legendre.o ProdKde2d.o $(UTILDIR)/file_io_utils.o
	$(CXX) $(CXXFLAGS) $^ -o $@

make_contours.o : make_contours.cc Kde2d.h
	$(CXX) $(CXXFLAGS) -I$(INCDIR) -c $< -o $@

Kde2d.o : Kde2d.cc Kde2d.h gauss_legendre.h $(UTILDIR)/file_io_utils.h
	$(CXX) $(CXXFLAGS) -I$(INCDIR) -c $< -o $@

ProdKde2d.o : ProdKde2d.cc ProdKde2d.h Kde2d.h
	$(CXX) $(CXXFLAGS) -I$(INCDIR) -c $< -o $@

prepare_kde_data : prepare_kde_data.o $(UTILDIR)/file_io_utils.o
	$(CXX) $(CXXFLAGS) $^ -o $@

prepare_kde_data.o : prepare_kde_data.cc $(UTILDIR)/file_io_utils.h 
	$(CXX) $(CXXFLAGS) -I$(INCDIR) -c $< -o $@

subsample : subsample.o $(UTILDIR)/file_io_utils.o
	$(CXX) $(CXXFLAGS) $^ -o $@

subsample.o : subsample.cc $(UTILDIR)/file_io_utils.h 
	$(CXX) $(CXXFLAGS) -I$(INCDIR) -c $< -o $@

clean:
	@rm -f *~ *.o `find . -perm +100 -type f`
