c++ -O3 -shared -std=c++11 -fPIC    \
-I ../include/pybind11              \
-I NewtonKrylov ../*.cpp            \
`python3-config --cflags --ldflags` \
bindings.cpp -o NewtonKrylov.so
