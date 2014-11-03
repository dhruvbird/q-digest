#include <iostream>
#include <random>
#include <algorithm>
#include "qdigest.h"

void compare_percentiles(double p,
                         std::vector<int> const &v,
                         qdigest::QDigest const &digest) {
  int ip = (p * 1000 + 9) / 10;
  std::cerr << ip << "th percentile: " << v.at(p * v.size() - 1)
            << " v/s " << digest.percentile(p) << "\n";
}

void test_poisson_distribution(int n, int k, int seed) {
  std::cerr << "<< test_poisson_distribution >>\n";
  std::vector<int> v;
  int number = 1;
  int repeat = 1;
  bool flipped = false;
  for (; v.size() != n; ++number) {
    for (int i = 0; i < repeat && v.size() != n; ++i) {
      v.push_back(number);
    }
    // std::cerr << "number: " << number << ", repeat: " << repeat << "\n";
    // std::cerr << v.size() << "\n";
    if (v.size() <= n/2) {
      if (!flipped) repeat += 3;
    } else {
      if (!flipped) repeat += 3;
      flipped = true;
      repeat -= 3;
    }
    if (repeat < 1) repeat = 2;
  }

  auto w = v;
  std::sort(v.begin(), v.end());
  std::random_shuffle(w.begin(), w.end());

  qdigest::QDigest digest(k);
  for (size_t i = 0; i < w.size(); ++i) {
    digest.insert(w[i], 1);
  }

  std::cerr << "# nodes: " << digest.size() << "\n"
            << "Compression Ratio: " << digest.compression_ratio() << "\n";

  compare_percentiles(0.01, v, digest);
  compare_percentiles(0.02, v, digest);
  compare_percentiles(0.03, v, digest);
  for (double d = 0.05; d < 1.0; d += 0.05) {
    compare_percentiles(d, v, digest);
  }
}

void test_geometric_distribution(int n, int k, int seed) {
  std::cerr << "<< test_geometric_distribution >>\n";
  std::vector<int> v;
  int number = 1;
  int repeat = 1;
  for (; v.size() != n; number += 1, repeat *= 2) {
    for (int i = 0; i < repeat && v.size() != n; ++i) {
      v.push_back(number);
    }
  }
  auto w = v;
  std::random_shuffle(w.begin(), w.end());

  qdigest::QDigest digest(k);
  for (size_t i = 0; i < w.size(); ++i) {
    digest.insert(w[i], 1);
  }

  std::cerr << "# nodes: " << digest.size() << "\n"
            << "Compression Ratio: " << digest.compression_ratio() << "\n";

  compare_percentiles(0.01, v, digest);
  compare_percentiles(0.02, v, digest);
  compare_percentiles(0.03, v, digest);
  for (double d = 0.05; d < 1.0; d += 0.05) {
    compare_percentiles(d, v, digest);
  }
}

void test_random_distribution(int n, int k, int seed) {
  std::cerr << "<< test_random_distribution >>\n";
  srand(seed);
  qdigest::QDigest digest(k);
  std::vector<int> v;

  for (int i = 0; i < n; ++i) {
    auto number = rand() % n;
    v.push_back(number);
    digest.insert(number, 1);
  }
  std::cerr << "# nodes: " << digest.size() << "\n"
            << "Compression Ratio: " << digest.compression_ratio() << "\n";

  std::sort(v.begin(), v.end());

  for (double d = 0.05; d < 1.0; d += 0.05) {
    compare_percentiles(d, v, digest);
  }
}

int main() {
  const int K = 100;
  int seed = 377;
  int N = 65535;

  test_random_distribution(N, K, seed);
  std::cerr << std::endl;
  test_geometric_distribution(N, K, seed);
  std::cerr << std::endl;
  test_poisson_distribution(N, K, seed);
}
