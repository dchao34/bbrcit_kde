CXX ?= g++

CXXFLAGS = -Wall -std=c++11

all : weighted_down_sample split_generic

weighted_down_sample : weighted_down_sample.o file_io_utils.o
	$(CXX) $(CXXFLAGS) $^ -o $@

weighted_down_sample.o : weighted_down_sample.cc file_io_utils.h 
	$(CXX) $(CXXFLAGS) -c $< -o $@

split_generic : split_generic.o file_io_utils.o
	$(CXX) $(CXXFLAGS) $^ -o $@

split_generic.o : split_generic.cc file_io_utils.h 
	$(CXX) $(CXXFLAGS) -c $< -o $@

file_io_utils.o : file_io_utils.cc file_io_utils.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	@rm -f *~ *.o `find . -perm +100 -type f`
