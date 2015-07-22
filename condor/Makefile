CC ?= gcc
CXX ?= g++

CXXFLAGS = -Wall -std=c++11

INCFLAGS = -I../ -I../data
INC = ../gauss_legendre.h ../Kde2d.h ../ProdKde2d.h ../data/file_io_utils.h
OBJECTS = $(addprefix ../, gauss_legendre.o Kde2d.o ProdKde2d.o data/file_io_utils.o)

all : cache_kde_score

cache_kde_score : cache_kde_score.o $(OBJECTS)
	$(CXX) $^ -o $@

cache_kde_score.o : cache_kde_score.cc $(INC)
	$(CXX) $(CXXFLAGS) $(INCFLAGS) -c $< -o $@

clean:
	@rm -f *~ *.o `find . -perm +100 -type f`
