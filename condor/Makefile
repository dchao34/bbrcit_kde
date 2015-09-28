CC ?= gcc
CXX ?= g++

CXXFLAGS = -Wall -std=c++11

INCFLAGS = -I../ -I../data
INC = ../gauss_legendre.h ../fft.h ../Kde2d.h ../ProdKde2d.h ../ProdAdaKde2d.h ../data/file_io_utils.h
OBJECTS = $(addprefix ../, gauss_legendre.o Kde2d.o ProdKde2d.o ProdAdaKde2d.o fft.o data/file_io_utils.o)

all : cv_scan cache_kde_score cache_akde_score

cv_scan : cv_scan.o $(OBJECTS)
	$(CXX) $^ -o $@

cv_scan.o : cv_scan.cc $(INC)
	$(CXX) $(CXXFLAGS) $(INCFLAGS) -c $< -o $@

cache_akde_score : cache_akde_score.o $(OBJECTS)
	$(CXX) $^ -o $@

cache_akde_score.o : cache_akde_score.cc $(INC)
	$(CXX) $(CXXFLAGS) $(INCFLAGS) -c $< -o $@

cache_kde_score : cache_kde_score.o $(OBJECTS)
	$(CXX) $^ -o $@

cache_kde_score.o : cache_kde_score.cc $(INC)
	$(CXX) $(CXXFLAGS) $(INCFLAGS) -c $< -o $@

clean:
	@rm -f *~ *.o cache_akde_score cache_kde_score cv_scan
