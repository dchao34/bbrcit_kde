#include <cstdlib>
#include <iostream>
#include <fstream>
#include <random>

using namespace std;

void print_usage() {
  cout << "usage: ./generate_bimodal_gauss2d sample_size" << endl;
}

int main(int argc, char *argv[]) {

  if (argc != 2) {
    print_usage();
    return EXIT_FAILURE;
  }
  int sample_size = atoi(argv[1]);

  ofstream fout("bimodal_gaussian.csv");

  random_device rd;
  default_random_engine e(rd());
  uniform_real_distribution<> u(0, 1);
  normal_distribution<> g(0.0, 1.0);

  for (int i = 0; i < sample_size; ++i) {
    double x = g(e), y = g(e);
    if (u(e) < 0.5) {
      x *= 0.5; y *= 0.1;
      x += 0.1; y += 1.0;
    } else {
      x *= 0.3; y *= 0.25;
      x = 0.8660254 * x + 0.5 * y;
      y = -0.5 * x + 0.8660254 * y;
      x += 0.8; y += 0.4;
    }
    fout << x << " " << y << endl;
  }

  return 0;

}
