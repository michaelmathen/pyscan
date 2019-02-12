# pyscan
This is a python wrapper around several anomaly detection algorithms written in c++. These algorithms are from the paper 
Computing Approximate Statistical Discrepancy.
To compile this you will need:

* boost.python
* python 2.7 or python 3.x
* cmake

#Instructions

```
> mkdir build
> cd build
> cmake ..
> make
``` 
You should then be able to use this library as a standard python module by doing:
import pyscan
