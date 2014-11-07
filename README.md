q-digest
========

An implementation of the Q-Digest Data Structure in C++; documented [here](http://www.cs.virginia.edu/~son/cs851/papers/ucsb.sensys04.pdf).

The Q-Digest lets you efficiently compute approximate percentile of large data sets with integral values. The data structure is space constrained and is parameterized by **(K)**, the *compression factor*. Smaller values of *K* indicate more compression and less accuracy, whereas larger values of *K* indicate lesser compression and more accuracy.

The accuracy of the Q-Digest is reasonable for data that has values which are drawn uniformly from a random distribution of values. It does well on Poisson distributions, and extremely well on geometric distributions. The file [test.cpp](https://github.com/dhruvbird/q-digest/blob/master/test.cpp) shows the differences.


