# pyscan
This is a python wrapper around a large number of anomaly detection algorithms written in c++. These algorithms are from several papers:
https://arxiv.org/abs/1906.09381
https://arxiv.org/abs/1906.01693
https://arxiv.org/abs/1804.11287
https://dl.acm.org/citation.cfm?id=2996939

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

## Website
https://michaelmathen.github.io/pyscan/
