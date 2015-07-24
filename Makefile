CC ?= gcc
CXX ?= g++

CXXFLAGS = -Wall -std=c++11

UTILDIR = data

INCDIR = $(UTILDIR)
LIBDIR = $(UTILDIR)

all : cache_kde_score make_contours test_gauss test_kde test_cv cv_scan kde_scan

kde_scan : kde_scan.cc gauss_legendre.o Kde2d.o ProdKde2d.o $(UTILDIR)/file_io_utils.o
	$(CXX) $(CXXFLAGS) -I$(INCDIR) $^ -o $@

cv_scan : cv_scan.cc gauss_legendre.o Kde2d.o ProdKde2d.o $(UTILDIR)/file_io_utils.o
	$(CXX) $(CXXFLAGS) -I$(INCDIR) $^ -o $@

test_cv : test_cv.cc gauss_legendre.o Kde2d.o ProdKde2d.o $(UTILDIR)/file_io_utils.o
	$(CXX) $(CXXFLAGS) -I$(INCDIR) $^ -o $@

test_gauss : test_gauss.cc gauss_legendre.o Kde2d.o ProdKde2d.o $(UTILDIR)/file_io_utils.o
	$(CXX) $(CXXFLAGS) $^ -o $@

test_kde : test_kde.o Kde2d.o gauss_legendre.o ProdKde2d.o $(UTILDIR)/file_io_utils.o
	$(CXX) $(CXXFLAGS) $^ -o $@

test_kde.o : test_kde.cc Kde2d.h ProdKde2d.h
	$(CXX) $(CXXFLAGS) -I$(INCDIR) -c $< -o $@

gauss_legendre.o : gauss_legendre.c gauss_legendre.h
	$(CC) -Wall -c $< -o $@

cache_kde_score : cache_kde_score.o gauss_legendre.o Kde2d.o ProdKde2d.o $(UTILDIR)/file_io_utils.o
	$(CXX) $(CXXFLAGS) $^ -o $@

cache_kde_score.o : cache_kde_score.cc ProdKde2d.h $(UTILDIR)/file_io_utils.o
	$(CXX) $(CXXFLAGS) -I$(INCDIR) -c $< -o $@

make_contours : make_contours.o Kde2d.o gauss_legendre.o ProdKde2d.o $(UTILDIR)/file_io_utils.o
	$(CXX) $(CXXFLAGS) $^ -o $@

make_contours.o : make_contours.cc Kde2d.h
	$(CXX) $(CXXFLAGS) -I$(INCDIR) -c $< -o $@

Kde2d.o : Kde2d.cc Kde2d.h gauss_legendre.h $(UTILDIR)/file_io_utils.h
	$(CXX) $(CXXFLAGS) -I$(INCDIR) -c $< -o $@

ProdKde2d.o : ProdKde2d.cc ProdKde2d.h Kde2d.h
	$(CXX) $(CXXFLAGS) -I$(INCDIR) -c $< -o $@

clean:
	@rm -f *~ *.o `find . -perm +100 -type f`
