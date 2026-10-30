# task-async
This repository includes the files for testing task-based approach on HPE resources.

## 1. Compilation
For compilation of `jacobi_ompss.cpp` in directory `/opt/ompss/test`:
```
clang++ -O3 -march=native -ffast-math -fompss-2 jacobi_ompss.cpp -lnuma -o jacobi_ompss
```
Consider change `-march=native` to the necessary flags for optimization on the platforms.

To run:
```
taskset -c <CORES> ./jacobi_ompss --n=64 --NX=512 --NY=512 --nt=100
```
n: block size

NX/NY: matrix dimension

nt: time steps

## 2. Runtime Settings

If necessary, add the paths pointing to the OmpSs-2 installs, like:
```
export LD_RUN_PATH="$LD_RUN_PATH:$HOME/package/ompss-2/lib"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$HOME/package/ompss-2/lib"
```
Runtime settings can be changed using env var `NANOS6_CONFIG_OVERRIDE`, a common setting I use is
```
export NANOS6_CONFIG_OVERRIDE="loader.report_prefix=#,numa.tracking=on,version.instrument=none"
```
The `loader.report` gives the verbose log for runtime, `numa.tracking` is enabled for NUMA APIs, `version.instrument` controls the instrumentation of the program and default to none. 

For instrumentation, change `none` to `ovni`, which should generates traces in a folder `/ovni` under the working directory.
Use the following to read and join the traces in `.prv` format so we can visualize the execution:
```
<PATH_TO_OVNI_INSTALL>/ovni/bin/ovniemu ./ovni
```
`.prv` traces and `.pcf` config files should be generated in the `./ovni` folder.
