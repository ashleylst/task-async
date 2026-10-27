# task-async
This repository includes the files for testing task-based approach on HPE resources.

For compilation of `jacobi_ompss.cpp` in directory `/opt/ompss/test`:
```
clang++ -O3 -march=native -ffast-math -fompss-2 jacobi_ompss.cpp -lnuma -o jacobi_ompss
```
Consider change `-march=native` to the necessary flags for optimization on the platforms.

To run:
```
./jacobi_ompss --n=64 --NX=512 --NY=512 --nt=100
```
n: block size

NX/NY: matrix dimension

nt: time steps
