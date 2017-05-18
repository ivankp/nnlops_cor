#include <vector>

#include "mat.hh"

using std::cout;
using std::endl;

template <typename T>
sym_mat<double> cov(const T& h, unsigned i) {
  return { h, [i](const auto& bin){ return bin.at(i)-bin[0]; } };
}

int main() {

  std::vector<std::vector<double>> hist {
    {10, 11, 8},
    {100, 102, 98},
    {1000, 1003, 998}
  };

  auto cov1 = cov(hist,1);
  auto cor1 = cor(cov1);

  cout << "cov 1\n" << cov1 << endl;
  cout << "cor 1\n" << cor1.cor << endl;

  auto cov2 = cov(hist,2);
  auto cor2 = cor(cov2);

  cout << "cov 2\n" << cov2 << endl;
  cout << "cor 2\n" << cor2.cor << endl;

  cov1 += cov2;
  auto cor12 = cor(cov1);

  cout << "cov 1 & 2\n" << cov1 << endl;
  cout << "cor 1 & 2\n" << cor12.cor << endl;
}
