#include <iostream>
#include "qdigest.h"

int main() {
  const int K = 10;
  qdigest::QDigest digest(K);

  digest.insert(10, 1);

  digest.printTree(std::cerr);

  digest.insert(51, 2);

  digest.printTree(std::cerr);

  digest.insert(52, 40);
  digest.insert(54, 20);

  digest.printTree(std::cerr);

  digest.insert(7, 2);

  digest.printTree(std::cerr);

  digest.insert(9, 2);
  digest.insert(21, 10);

  digest.printTree(std::cerr);

  std::cerr << "50th percentile is: " << digest.percentile(0.5) << "\n";
  std::cerr << "70th percentile is: " << digest.percentile(0.7) << "\n";
  std::cerr << "2nd percentile is: " << digest.percentile(0.02) << "\n";

}
