#include <cstdlib>
#include <iostream>
#include <fstream>
#include <random>

using namespace std;

void print_usage() {
  cout << "usage: ./generate_bimodal_gauss1d sample_size" << endl;
}

int main(int argc, char *argv[]) {

  if (argc != 2) {
    print_usage();
    return EXIT_FAILURE;
  }
  int sample_size = atoi(argv[1]);

  ofstream fout("bimodal_gauss1d.csv");

  random_device rd;
  default_random_engine e(rd());
  uniform_real_distribution<> u(0, 1);
  normal_distribution<> g(0.0, 1.0);

  for (int i = 0; i < sample_size; ++i) {
    double x = g(e);
    if (u(e) < 0.8) {
      x *= 0.3;
      x -= 0.5;
    } else {
      x *= 0.2;
      x += 0.7;
    }
    fout << x << endl;
  }

  return 0;

}
