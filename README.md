q-digest
========

An implementation of the Q-Digest Data Structure in C++; documented [here](http://www.cs.virginia.edu/~son/cs851/papers/ucsb.sensys04.pdf).

The Q-Digest lets you efficiently compute approximate percentile of large data sets with integral values. The data structure is space constrained and is parameterized by **(K)**, the *compression factor*. Smaller values of *K* indicate more compression and less accuracy, whereas larger values of *K* indicate lesser compression and more accuracy.

The accuracy of the Q-Digest is reasonable for data that has values which are drawn uniformly from a random distribution of values. It does well on Poisson distributions, and extremely well on geometric distributions. The file [test.cpp](https://github.com/dhruvbird/q-digest/blob/master/test.cpp) shows the differences.

### Public API

```C++
namespace qdigest {
  QDigest::QDigest(size_t K, size_t upper_bound = 1);
  
  QDigest::QDigest(QDigest &&rhs);
  
  void QDigest::swap(QDigest &rhs);
  
  void QDigest::insert(size_t key, unsigned int count);
  
  size_t QDigest::size() const;
  
  bool QDigest::empty() const;
  
  // Returns (# internal nodes) / (sum of counts for all keys)
  double QDigest::compression_ratio() const;
  
  // Returns the 100p'th percentile. p takes values in the range [0..100]
  size_t QDigest::percentile(double p) const;
  
  // Serialize the structure into a string for storage/transmission
  std::string QDigest::toString() const;

  // De-serialize the structure into *this given a serialized representation
  void QDigest::fromString(std::string const &ser);

  // Merge 'rhs' into *this
  void QDigest::merge(QDigest const &rhs);
  QDigest& QDigest::operator+=(QDigest const &rhs);

  // Pretty print the Q-Digest (for debugging)
  std::ostream& operator<<(std::ostream &out, QDigest const &digest);
}
```

You can not directly copy a Q-Digest into another, since copying a Q-Digest is an expensive operation and involves copying/cloning all the internal nodes. However, you can do the following to make a copy if you wish:
```C++
qdigest::QDigest old(100, 10);
old.insert(1, 3);
old.insert(4, 2);

qdigest::QDigest curr(100);
// Since 'curr' started out empty, 'curr' is now a copy of 'old'
curr += old;
```
