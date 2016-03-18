#include <vector>
#include <ostream>

#include "cv_utils.h"

void write_2dgrid_values(std::ostream &os, const std::vector<double> &values,
                         double start_x, double end_x, int steps_x, 
                         double start_y, double end_y, int steps_y) {

  double delta_x = (end_x-start_x)/steps_x;
  double delta_y = (end_y-start_y)/steps_y;
  for (int i = 0; i <= steps_x; ++i) { os << start_x + i * delta_x << " "; } 
  os << std::endl;
  for (int i = 0; i <= steps_y; ++i) { os << start_y + i * delta_y << " "; } 
  os << std::endl;

  for (size_t i = 0; i < values.size(); ++i) {
    if (i % (steps_x+1) == 0 && i) { os << std::endl; }
    os << values[i] << " ";
  }
  os << std::endl;
}


template <>
void generate_2dgrid<std::pair<double,double>>(
    std::vector<std::pair<double,double>> &grid, 
    double start_x, double end_x, int steps_x,
    double start_y, double end_y, int steps_y) {
  
  grid.clear();

  double delta_x = (end_x-start_x)/steps_x;
  double delta_y = (end_y-start_y)/steps_y;
  for (int j = 0; j <= steps_y; ++j) {
    for (int i = 0; i <= steps_x; ++i) {
      grid.push_back({start_x+i*delta_x, start_y+j*delta_y});
    }
  }

  return;
}
