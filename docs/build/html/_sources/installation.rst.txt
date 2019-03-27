Installation
==================

You can get pyscan from https://github.com/michaelmathen/pyscan.git. To compile you will need the following dependencies.

*  boost.python
*  python python 3.x
*  cmake

In addition, this library needs a modern version of gcc or clang that supports the c++17 standard. The `Google Test Framework <https://github.com/google/googletest>`_. will be downloaded to build unit tests.::

    git clone https://github.com/michaelmathen/pyscan.git
    cd pyscan
    mkdir build
    cd build
    cmake ..
    make

This will by default do a release build creating a number of files in the build directory. You can add this directory to your PYTHONPATH for instance:

    export PYTHONATH=$PYTHONPATH:/path/to/pyscan/build/ 

You should then be able to use this library as a standard python module by doing 

    import pyscan

You can also just copy and paste the libpyscan.so and pyscan.py file in the build directory into the directory you wish to use pyscan from.

