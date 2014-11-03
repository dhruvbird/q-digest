#include <iostream>
#include "qdigest.h"

int main() {
  const int K = 10;
  qdigest::QDigest digest(K);

  digest.insert(10, 1);

  std::cerr << digest << std::endl;

  digest.insert(51, 2);

  std::cerr << digest << std::endl;

  digest.insert(52, 40);
  digest.insert(54, 20);

  std::cerr << digest << std::endl;

  digest.insert(7, 2);

  std::cerr << digest << std::endl;

  digest.insert(9, 2);
  digest.insert(21, 10);

  std::cerr << digest << std::endl;

  std::cerr << "50th percentile is: " << digest.percentile(0.5) << "\n";
  std::cerr << "70th percentile is: " << digest.percentile(0.7) << "\n";
  std::cerr << "2nd percentile is: " << digest.percentile(0.02) << "\n";

  qdigest::QDigest digest2(K);
  digest2.fromString(digest.toString());

  std::cerr << digest << std::endl;
}
