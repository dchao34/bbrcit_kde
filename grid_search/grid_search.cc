#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <string>
#include <chrono>
#include <algorithm>

#include <boost/program_options.hpp>

#include <Kernels/EpanechnikovProductKernel2d.h>
#include <Kernels/GaussianProductKernel2d.h>
#include <KernelDensity.h>

#include "cv_utils.h"
#include "custom_program_option_utils.h"

namespace {
  using KernelType = bbrcit::EpanechnikovProductKernel2d<float>;
  using KernelDensityType = bbrcit::KernelDensity<2,KernelType,double>;
  using DataPointType = KernelDensityType::DataPointType;
}

namespace po = boost::program_options;

void grid_search(const po::variables_map &vm);

// helper function for printing vectors
template<class T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v) {
  std::copy(v.begin(), v.end(), std::ostream_iterator<T>(os, " "));
  return os;
}

int main(int argc, char **argv) {

  try {

    // define program options 
    po::options_description generic("Generic options");
    generic.add_options()
        ("help,h", "produce help message")
    ;

    po::options_description config("Configuration options");
    config.add_options()
        ("rel_tol", po::value<double>(), "relative tolerace for the evaluation error. ")
        ("abs_tol", po::value<double>(), "absolute tolerace for the evaluation error. ")
        ("cuda_device_number", po::value<int>(), "cuda gpu device number used for this session. ")
        ("gpu_block_size", po::value<int>(), "block size for the gpu kernel. ")
        ("refpt_max_leaf_size", po::value<int>(), "maximum leaf size of reference point tree. ")
        ("qgrid_max_leaf_size", po::value<int>(), "maximum leaf size of the query grid tree. ")

        ("input_refpts_fname", po::value<std::string>(), "path to the input reference points. ")
        ("output_scatter_fname", po::value<std::string>(), "path to output matplotlib scatter plot data. ")

        ("skip_cross_validation", po::value<bool>(), "skip cross validation to go straight to evaluation. ")
        ("use_manual_bwgrid", po::value<bool>(), "use manually specified bandwidth grid for grid search. ")
        ("cv_method", po::value<std::string>(), "method to use for cross validation. one of the following: \n"
                                                "lsq_conv: least squares based on convolution kernel. \n"
                                                "lsq_numint: least squares based on numerical integration.")

        ("manual_bwx", po::value<std::string>(), "manual bandwidths in the x-dimension. ")
        ("manual_bwy", po::value<std::string>(), "manual bandwidths in the y-dimension. ")

        ("start_bwx", po::value<double>(), "starting x bandwidth for the automatic bandwidth grid. ")
        ("end_bwx", po::value<double>(), "ending x bandwidth for the automatic bandwidth grid. ")
        ("steps_bwx", po::value<int>(), "bandwidth increments in x for the automatic bandwidth grid. ")
        ("start_bwy", po::value<double>(), "starting y bandwidth for the automatic bandwidth grid. ")
        ("end_bwy", po::value<double>(), "ending y bandwidth for the automatic bandwidth grid. ")
        ("steps_bwy", po::value<int>(), "bandwidth increments in y for the automatic bandwidth grid. ")
        ("output_gridsearch_fname", po::value<std::string>(), "path to output matplotlib grid search data. ")

        ("start_qix", po::value<double>(), "starting x value for the numerical integration grid. ")
        ("end_qix", po::value<double>(), "ending x value for the numerical integration grid. ")
        ("steps_qix", po::value<int>(), "increments in x values for the numerical integration grid. ")
        ("start_qiy", po::value<double>(), "starting y value for the numerical integration grid. ")
        ("end_qiy", po::value<double>(), "ending y value for the numerical integration grid. ")
        ("steps_qiy", po::value<int>(), "increments in y values for the numerical integration grid. ")

        ("use_gridsearch_best", po::value<bool>(), "whether to use the best bandwidth found during grid search. ")
        ("eval_bwx", po::value<double>(), "x bandwidth to evaluate and output for plotting. ")
        ("eval_bwy", po::value<double>(), "y bandwidth to evaluate and output for plotting. ")

        ("start_qx", po::value<double>(), "starting x value for the query grid. ")
        ("end_qx", po::value<double>(), "ending x value for the query grid. ")
        ("steps_qx", po::value<int>(), "increments in x values for the query grid. ")
        ("start_qy", po::value<double>(), "starting y value for the query grid. ")
        ("end_qy", po::value<double>(), "ending y value for the query grid. ")
        ("steps_qy", po::value<int>(), "increments in y values for the query grid. ")

        ("alpha", po::value<double>(), "adaptive kernel sensitivity. ")
        ("output_eval_fname", po::value<std::string>(), "path to output matplotlib data to for plotting kde. ")
        ("output_adaptive_eval_fname", po::value<std::string>(), "path to output matplotlib data to for plotting adaptive kde. ")

    ;

    po::options_description hidden("Hidden options");
    hidden.add_options()
        ("config_file", po::value<std::string>(), "name of a configuration file. ")
    ;

    po::options_description cmdline_options;
    cmdline_options.add(generic).add(config).add(hidden);

    po::options_description config_file_options;
    config_file_options.add(config);

    po::options_description visible;
    visible.add(generic).add(config);

    po::positional_options_description p;
    p.add("config_file", -1);

    // parse program options and configuration file
    po::variables_map vm;
    store(po::command_line_parser(argc, argv).
          options(cmdline_options).positional(p).run(), vm);
    notify(vm);

    if (vm.count("help") || !vm.count("config_file")) {
      std::cout << std::endl;
      std::cout << "Usage: cross_validate [options] config_fname" << std::endl;
      std::cout << visible << "\n";
      return 0;
    }

    std::ifstream fin(vm["config_file"].as<std::string>());
    if (!fin) {
      std::cout << "cannot open config file: ";
      std::cout << vm["config_file"].as<std::string>() << std::endl;
      return 0;
    }

    store(parse_config_file(fin, config_file_options), vm);
    notify(vm);

    // launch grid seach 
    grid_search(vm);


  } catch(std::exception& e) {

    std::cerr << "error: " << e.what() << "\n";
    return 1;

  } catch(...) {

    std::cerr << "Exception of unknown type!\n";
    return 1;
  }


  return 0;
}


void grid_search(const po::variables_map &vm) {

  // 0. setup general utilities
  // --------------------------

  std::ofstream fout;
  std::chrono::high_resolution_clock::time_point start, end;
  std::chrono::duration<double, std::milli> elapsed;
  double rel_tol = vm["rel_tol"].as<double>();
  double abs_tol = vm["abs_tol"].as<double>();

  std::cout << "+ performance parameters: \n" << std::endl;

  std::cout << "  relative tolerance: " << rel_tol << std::endl;
  std::cout << "  absolute tolerance: " << abs_tol << std::endl;
  std::cout << std::endl;

#ifdef __CUDACC__
  int cuda_device_number = vm["cuda_device_number"].as<int>();
  cudaSetDevice(cuda_device_number);
  cudaDeviceProp deviceProp;
  cudaGetDeviceProperties(&deviceProp, cuda_device_number);
  std::cout << "  gpu device number used for this session: ";
  std::cout << cuda_device_number << "\n";
  std::cout << "  device name: " << deviceProp.name << std::endl;

  int gpu_block_size = vm["gpu_block_size"].as<int>();
  std::cout << "  gpu block size: " << gpu_block_size << std::endl;
  std::cout << std::endl;
#endif

  int refpt_max_leaf_size = vm["refpt_max_leaf_size"].as<int>();
  std::cout << "  reference point tree max leaf size: ";
  std::cout << refpt_max_leaf_size << std::endl;

  int qgrid_max_leaf_size = vm["qgrid_max_leaf_size"].as<int>();
  std::cout << "  query grid tree max leaf size: ";
  std::cout << qgrid_max_leaf_size << std::endl;


  std::cout << std::endl;


  // 1. read in reference data
  // -------------------------

  // setup parameters
  std::string input_refpts_fname = vm["input_refpts_fname"].as<std::string>();
  std::string output_scatter_fname = vm["output_scatter_fname"].as<std::string>();

  // read data from file
  std::cout << "+ reading in reference data. \n" << std::endl;

  std::vector<DataPointType> data;
  read_2dpoints(input_refpts_fname, data);

  std::cout << "  => found " << data.size() << " points. \n" << std::endl;

  // write data for making scatter plots in matplotlib
  fout.open(output_scatter_fname);
  write_2dscatter_data(fout, data); 
  fout.close();
  std::cout << "  => wrote results to: " << output_scatter_fname << ". \n" << std::endl;


  // 2. construct the kernel density estimator
  // -----------------------------------------

  // build kernel density estimator
  std::cout << "+ constructing kernel density estimator. \n" << std::endl;

  start = std::chrono::high_resolution_clock::now();
  KernelDensityType kde(data, refpt_max_leaf_size);
  end = std::chrono::high_resolution_clock::now();
  elapsed = end - start;
  std::cout << "  => running time: " << elapsed.count() << " ms. \n" << std::endl;


  // 3. grid search cross validation
  // -------------------------------

  // setup parameters
  bool skip_cross_validation = vm["skip_cross_validation"].as<bool>();
  bool use_manual_bwgrid = vm["use_manual_bwgrid"].as<bool>();
  std::string cv_method = vm["cv_method"].as<std::string>();

  std::string manual_bwx = vm["manual_bwx"].as<std::string>();
  std::string manual_bwy = vm["manual_bwy"].as<std::string>();

  double start_bwx = vm["start_bwx"].as<double>();
  double end_bwx = vm["end_bwx"].as<double>();
  int steps_bwx = vm["steps_bwx"].as<int>();
  double start_bwy = vm["start_bwy"].as<double>();
  double end_bwy = vm["end_bwy"].as<double>();
  int steps_bwy = vm["steps_bwy"].as<int>();

  double start_qix = vm["start_qix"].as<double>();
  double end_qix = vm["end_qix"].as<double>();
  int steps_qix = vm["steps_qix"].as<int>();
  double start_qiy = vm["start_qiy"].as<double>();
  double end_qiy = vm["end_qiy"].as<double>();
  int steps_qiy = vm["steps_qiy"].as<int>();

  std::string output_gridsearch_fname = vm["output_gridsearch_fname"].as<std::string>();

  double best_bwx, best_bwy;
  if (!skip_cross_validation) {
    
    // define recognized cross validation methods 
    std::set<std::string> recognized_cv_methods;
    recognized_cv_methods.insert("lsq_conv");
    recognized_cv_methods.insert("lsq_numint");

    // decide which cv method to use
    if (recognized_cv_methods.find(cv_method) == recognized_cv_methods.end()) {
      throw std::invalid_argument("grid_search: unrecognized cv method. ");
    }

    if (cv_method == "lsq_conv") {
      std::cout << "+ using convolution based least squares cross validation. \n" << std::endl;
    } else {
      std::cout << "+ using numerical integration based least squares cross validation. " << std::endl;
      std::cout << "    grid dimensions: " << steps_qix << " by " << steps_qiy << std::endl;
      std::cout << "    grid bounds: ";
      std::cout << "[" << start_qix << ", " << end_qix << "]" << " x ";
      std::cout << "[" << start_qiy << ", " << end_qiy << "]\n" << std::endl;
    }

    // setup the search grid
    std::vector<std::pair<double, double>> bandwidth_grid;
    if (!use_manual_bwgrid) {
      generate_2dgrid(bandwidth_grid, start_bwx, end_bwx, steps_bwx,
                                      start_bwy, end_bwy, steps_bwy);
    } else {
      std::vector<double> manual_bwx = tokenize<double>(vm["manual_bwx"].as<std::string>());
      std::vector<double> manual_bwy = tokenize<double>(vm["manual_bwy"].as<std::string>());
      bandwidth_grid = make_crossed_vector_pairs(manual_bwx, manual_bwy);
    }

    if (use_manual_bwgrid) {
      std::cout << "  using manual bandwidth grid: \n";
      std::cout << "    grid points: [" << manual_bwx << "] x [" << manual_bwy << "]\n";
    } else {
      std::cout << "  using automatic bandwidth grid: \n";
      std::cout << "    grid dimensions: " << steps_bwx << " by " << steps_bwy << std::endl;
      std::cout << "    grid bounds: ";
      std::cout << "[" << start_bwx << ", " << end_bwx << "]" << " x ";
      std::cout << "[" << start_bwy << ", " << end_bwy << "]\n";
    }
    std::cout << std::endl;

    // grid search
    std::cout << "  grid searching: \n" << std::endl;
    
    start = std::chrono::high_resolution_clock::now();
    
    std::vector<double> cv_scores;
    double curr_cv, best_cv = std::numeric_limits<double>::max();
    for (const auto &p : bandwidth_grid) {
      kde.kernel().set_bandwidths(p.first, p.second);

      if (cv_method == "lsq_conv") {
#ifndef __CUDACC__
        curr_cv = kde.lsq_convolution_cross_validate(rel_tol, abs_tol);
#else
        curr_cv = kde.lsq_convolution_cross_validate(rel_tol, abs_tol, gpu_block_size);
#endif
      } else if (cv_method == "lsq_numint") {
#ifndef __CUDACC__
        curr_cv = lsq_numint_cross_validate(kde, 
                                            start_qix, end_qix, steps_qix,
                                            start_qiy, end_qiy, steps_qiy,
                                            rel_tol, abs_tol, qgrid_max_leaf_size);
#else
        curr_cv = lsq_numint_cross_validate(kde, 
                                            start_qix, end_qix, steps_qix,
                                            start_qiy, end_qiy, steps_qiy,
                                            rel_tol, abs_tol, qgrid_max_leaf_size, 
                                            gpu_block_size);
#endif
      } else {
        assert(false);
      }

      cv_scores.push_back(curr_cv);
      if (curr_cv < best_cv) { 
        best_cv = curr_cv; 
        best_bwx = p.first; best_bwy = p.second; 
      }
      std::cout << "  (" << p.first << ", " << p.second << ") " << curr_cv << std::endl;
    }
    std::cout << std::endl;

    end = std::chrono::high_resolution_clock::now();
    elapsed = end - start;

    std::cout << "  => best bandwidth: ";
    std::cout << "(" << best_bwx << ", " << best_bwy << "), "; 
    std::cout << "cv score = " << best_cv << "\n" << std::endl;

    std::cout << "  => running time per iteration: ";
    std::cout << elapsed.count() / bandwidth_grid.size() << " ms  \n" << std::endl;

    // write cross validation scores for matplotlib
    if (!use_manual_bwgrid) {
      fout.open(output_gridsearch_fname);
      write_2dgrid_values(fout, cv_scores, start_bwx, end_bwx, steps_bwx, start_bwy, end_bwy, steps_bwy); 
      fout.close();
      std::cout << "  => wrote results to: " << output_gridsearch_fname << ". \n" << std::endl;
    }

  }


  // 4. evaluate cross validated density for plotting
  // ------------------------------------------------

  // setup parameters
  bool use_gridsearch_best = vm["use_gridsearch_best"].as<bool>();
  double eval_bwx = vm["eval_bwx"].as<double>();
  double eval_bwy = vm["eval_bwy"].as<double>();

  double start_qx = vm["start_qx"].as<double>();
  double end_qx = vm["end_qx"].as<double>();
  int steps_qx = vm["steps_qx"].as<int>();
  double start_qy = vm["start_qy"].as<double>();
  double end_qy = vm["end_qy"].as<double>();
  int steps_qy = vm["steps_qy"].as<int>();
  std::string output_eval_fname = vm["output_eval_fname"].as<std::string>();

  if (skip_cross_validation) { use_gridsearch_best = false; }

  if (use_gridsearch_best) {
    eval_bwx = best_bwx;
    eval_bwy = best_bwy;
  }


  // generate query grid
  std::vector<DataPointType> grid, queries;
  generate_2dgrid(grid, start_qx, end_qx, steps_qx, start_qy, end_qy, steps_qy);

  std::cout << "+ evaluating cross validated kde over plotting grid. \n" << std::endl;
  std::cout << "  grid dimensions: " << steps_qx << "x" << steps_qy << std::endl;
  std::cout << "  grid bounds: ";
  std::cout << "[" << start_qx << ", " << end_qx << "]" << " x ";
  std::cout << "[" << start_qy << ", " << end_qy << "]" << std::endl;
  std::cout << "  x bandwidth: " << eval_bwx << std::endl;
  std::cout << "  y bandwidth: " << eval_bwy << "\n" << std::endl;

  // evaluate
  kde.kernel().set_bandwidths(eval_bwx, eval_bwy);
  queries = grid;

  start = std::chrono::high_resolution_clock::now();

#ifndef __CUDACC__
  kde.eval(queries, rel_tol, abs_tol, qgrid_max_leaf_size);
#else
  kde.eval(queries, rel_tol, abs_tol, qgrid_max_leaf_size, gpu_block_size);
#endif

  end = std::chrono::high_resolution_clock::now();
  elapsed = end - start;
  std::cout << "  => running time: " << elapsed.count() << " ms. \n" << std::endl;

  // write evaluation results for matplotlib
  fout.open(output_eval_fname);
  write_2dgrid_values(fout, queries, start_qx, end_qx, steps_qx, start_qy, end_qy, steps_qy); 
  fout.close();
  std::cout << "  => wrote results to: " << output_eval_fname << ". \n" << std::endl;


  // 5. adaptive density
  // -------------------

  // setup parameters
  double alpha = vm["alpha"].as<double>();

  // convert 
  std::cout << "+ converting to adaptive density. \n" << std::endl;
  std::cout << "  alpha: " << alpha << std::endl;
  std::cout << std::endl;

  start = std::chrono::high_resolution_clock::now();

#ifndef __CUDACC__
  kde.adapt_density(alpha, rel_tol, abs_tol);
#else
  kde.adapt_density(alpha, rel_tol, abs_tol, gpu_block_size);
#endif

  end = std::chrono::high_resolution_clock::now();
  elapsed = end - start;
  std::cout << "  => running time: " << elapsed.count() << " ms. \n" << std::endl;


  // 6. evaluate adaptive density for plotting
  // -----------------------------------------

  // setup parameters
  std::string output_adaptive_eval_fname = vm["output_adaptive_eval_fname"].as<std::string>();

  // evaluate
  std::cout << "+ evaluating adaptive density over plotting grid. \n" << std::endl;

  queries = grid;

  start = std::chrono::high_resolution_clock::now();

#ifndef __CUDACC__
  kde.eval(queries, rel_tol, abs_tol, qgrid_max_leaf_size);
#else
  kde.eval(queries, rel_tol, abs_tol, qgrid_max_leaf_size, gpu_block_size);
#endif

  end = std::chrono::high_resolution_clock::now();
  elapsed = end - start;
  std::cout << "  => running time: " << elapsed.count() << " ms. \n" << std::endl;

  fout.open(output_adaptive_eval_fname);
  write_2dgrid_values(fout, queries, start_qx, end_qx, steps_qx, start_qy, end_qy, steps_qy); 
  fout.close();
  std::cout << "  => wrote results to: " << output_adaptive_eval_fname << ". \n" << std::endl;

  return;
}
