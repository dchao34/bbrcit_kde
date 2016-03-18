#ifndef CV_UTILS_H__
#define CV_UTILS_H__

#include <vector>
#include <ostream>
#include <string>
#include <utility>

// read weighted 2d points from file `fname`. each row is corresponds to 
// a point. the last column is the weight; otherwise, the value in column `i`
// is the value of coodinate `i`. 
template <typename PointT>
void read_2dpoints(const std::string &fname, std::vector<PointT> &points);


// generate a 2-dimensional grid of points. 
template <typename PointT>
void generate_2dgrid(std::vector<PointT> &grid, 
                     double start_x, double end_x, int steps_x,
                     double start_y, double end_y, int steps_y);
template <>
void generate_2dgrid<std::pair<double,double>>(
    std::vector<std::pair<double,double>> &grid, 
    double start_x, double end_x, int steps_x,
    double start_y, double end_y, int steps_y);

// cross two vectors into a vector of pairs: 
// e.g. [a1, ..., aM] x [b1, bN] => 
//      [(a1,b1), (a1, b2), ... (a1, bN), ... (aM, b1), ..., (aM, bN) ]
template<typename T>
std::vector<std::pair<T,T>>
make_crossed_vector_pairs(const std::vector<T> &v1, const std::vector<T> &v2);

// write 2-dimensional values attached to grid points generated using
// a call to `generate_2dgrid(...)` with arguments `(start|end|steps)_(x|y)`.
void write_2dgrid_values(std::ostream &os, std::vector<double> &values,
                         double start_x, double end_x, int steps_x, 
                         double start_y, double end_y, int steps_y);
void write_2dgrid_values(std::ostream &os, const std::vector<double> &values,
                         double start_x, double end_x, int steps_x, 
                         double start_y, double end_y, int steps_y);
template<typename PointT>
void write_2dgrid_values(std::ostream &os, std::vector<PointT> &point_values,
                         double start_x, double end_x, int steps_x, 
                         double start_y, double end_y, int steps_y);

// write 2-dimensional points in `data` to output stream `os`. 
// each row corresponds to a single point whose value in the 
// `i`th dimension corresponds to the value in column `i`. 
template<typename PointT>
void write_2dscatter_data(std::ostream &os, const std::vector<PointT> &data);


#include "cv_utils_impl.h"


#endif
