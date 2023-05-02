# pyscan
This is a python wrapper around a large number of anomaly detection algorithms written in c++. These algorithms are from several papers:
* [The Kernel Scan Statistic](https://arxiv.org/abs/1906.09381)
* [Scalable Spatial Scan Statistics for Trajectories](https://arxiv.org/abs/1906.01693)
* [Computing Approximate Statistical Discrepancy](https://arxiv.org/abs/1804.11287)
* [Scalable Spatial Scan Statistics through Sampling](https://dl.acm.org/citation.cfm?id=2996939)
* [Spatial Scan Statistics: Approximations and Performance Studies](http://www.cs.utah.edu/~jeffp/papers/stat-disc-KDD06.pdf)
* [The Hunting of the Bump: On Maximizing Statistical Discrepancy](http://www.cs.utah.edu/~jeffp/papers/stat-disc-SODA06.pdf)

If you are interested in using any of these for comparison studies, finding disease outbreaks, collaberations, etc please reach out and email me. I am happy to be of assistance.

To compile this you will need:

* python python 3.x
* cgal
* gsl
* cmake

## Instructions

```
> git submodule update --init
> mkdir build
> cd build
> cmake ..
> make
``` 
You should then be able to use this library as a standard python module by doing:
import pyscan

### For Mac M1 Chip User
In CMakeLists.txt, you need to comment out 
```
set(CMAKE_CXX_FLAGS_RELEASE "-fPIC -w -O2 -march=native -DNDEBUG")
set(CMAKE_C_FLAGS_RELEASE "-fPIC -w -O2 -march=native -DNDEBUG")
```
and uncomment
```
set(CMAKE_CXX_FLAGS_RELEASE "-fPIC -w -O2 -mcpu=apple-m1 -DNDEBUG")
set(CMAKE_C_FLAGS_RELEASE "-fPIC -w -O2 -mcpu=apple-m1 -DNDEBUG")
```

Make sure use this library on arm64 running Python with the same Python version that built the library. 


## Website
https://michaelmathen.github.io/pyscan/
