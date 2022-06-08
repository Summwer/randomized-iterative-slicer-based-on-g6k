# randomized-iterative-slicer-based-on-g6k
I've implemented the randomized iterative slicer algorithm in g6k kernel.

Guidance to implement the randomized iterative slicer algorithm:

```
cd kernel
g++ -O3 -march=native -Wp,-U_FORTIFY_SOURCE -fPIC -Ofast -march=native -ftree-vectorize -funroll-loops -std=c++11 -pthread -Wall -Wextra -fPIC -shared -o libsieving.so sieving.cpp control.cpp fht_lsh.cpp params.cpp cpuperf.cpp randomized_iterative_slicer.cpp
cd ..
export LD_LIBRARY_PATH="/home/cryptothesis/summer/debug/slicer/kernel:/home/cryptothesis/summer/debug/slicer/g6k-fplll/fplll/.lib:/home/cryptothesis/summer/debug/slicer"
g++ -o test test.cpp -L./kernel -lsieving -pthread -lfplll -lgmp -lmpfr
```
