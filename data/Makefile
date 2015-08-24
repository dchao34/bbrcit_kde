CXX ?= g++

CXXFLAGS = -Wall -std=c++11

BINARIES = dither prepare_kde_data subsample weighted_down_sample split_generic

all : $(BINARIES)

prepare_kde_data : prepare_kde_data.o file_io_utils.o
	$(CXX) $(CXXFLAGS) $^ -o $@

prepare_kde_data.o : prepare_kde_data.cc file_io_utils.h 
	$(CXX) $(CXXFLAGS) -c $< -o $@

dither : dither.o file_io_utils.o
	$(CXX) $(CXXFLAGS) $^ -o $@

dither.o : dither.cc file_io_utils.h 
	$(CXX) $(CXXFLAGS) -c $< -o $@

subsample : subsample.o file_io_utils.o
	$(CXX) $(CXXFLAGS) $^ -o $@

subsample.o : subsample.cc file_io_utils.h 
	$(CXX) $(CXXFLAGS) -c $< -o $@

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
	@rm -f *~ *.o $(BINARIES)
