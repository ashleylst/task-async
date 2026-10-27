#include <sstream>
#include <iostream>
#include <vector>
#include <random>
#include <iomanip>
#include <algorithm>
#include "chrono"
#include <getopt.h>
#include <cstdlib>
#include <nanos6.h>
#include <numa.h>
#include <numaif.h>

#define RP 10

void vector_double_define_random(double* phi, int start, int end) {
    std::random_device seed;
    std::mt19937 randomInt(seed());
    std::uniform_real_distribution<> unit(1.0, 10.0);

    for (int i = start; i < end; i++)
        phi[i] = unit(randomInt);

}

//idx does the wrapping here
static inline int idx(int row, int col, int N)
{
    return row * N + col;
}


void print_matrix(int NX, int NY, double *m){
    for (int i = 0; i < NY; i++) {
        for (int j = 0; j < NX; j++) {
            std::cout << m[idx(i, j, NX)] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}


double computeMedian(std::vector<long>& nums) {
    if (nums.empty()) {
        return 0.0;  // Handle empty vector case
    }

    // Make a copy to avoid modifying the original vector
    std::vector<long> sorted = nums;
    std::sort(sorted.begin(), sorted.end());

    size_t n = sorted.size();
    if (n % 2 == 0) {  // Even number of elements
        return (sorted[n/2 - 1] + sorted[n/2]) / 2.0;
    } else {  // Odd number of elements
        return sorted[n/2];
    }
}

void define_uniform(double* phi, int start, int end){
    for (int i = start; i < end; ++i) {
        phi[i] = i;
    }
}

void define_uniform_1(double* phi, double* phi_2, int NX, int NY, int n){
    int num_blocks = std::ceil((NY + n - 1) / n);
    for (int b = 0; b < num_blocks; b++)
    {
#pragma oss task
        {
            int y_end = (b == num_blocks - 1) && (NY % n != 0) ? NY : n * (b + 1);
            for (int i = b * n * NX; i < y_end * NX; ++i) {
                phi[i] = i;
		phi_2[i] = 0;
            }
        }
    }
}

static inline double calculate(double up, double down, double left, double right)
{
    return 0.25 * (up + down + left + right);
}



static inline void jacobi_kernel(int NX, int NY, int n, bool end, double* next,
                                 double *up,
                                 double *mid, double *down)
{
   /* int cpu = sched_getcpu();
    
    if (cpu == -1) {
        perror("sched_getcpu failed");
    }

    // 2. Get the NUMA node for that CPU
    int node = numa_node_of_cpu(cpu);

    std::stringstream ss;
    ss  << "Task id " << id << " is currently running on CPU " << cpu
              << " which belongs to NUMA node " << node << "\n";
    std::cout << ss.str();*/
    //std::osyncstream(std::cout) << "Task id " << id << " is currently running on CPU " << cpu 
    //          << " which belongs to NUMA node " << node << std::endl;

    /// top boundary
    if (up != nullptr)
    {
        for (int x = 1; x < NX - 1; x++) {
            next[idx(0, x, NX)] = calculate(up[idx(n - 1, x, NX)], mid[idx(1, x, NX)],
                                           mid[idx(0, x - 1, NX)], mid[idx(0, x + 1, NX)]);
        }
    }
    /// mid
    int y_end = (end && NY % n != 0) ? NY % n - 1 : n - 1;
    for (int y = 1; y < y_end; y++) {
        for (int x = 1; x < NX - 1; x++) {
            next[idx(y, x, NX)] = calculate(mid[idx(y-1, x, NX)], mid[idx(y+1, x, NX)],
                                           mid[idx(y, x-1, NX)], mid[idx(y, x+1, NX)]);
        }
    }
    /// bottom boundary
    if (down != nullptr)
    {
        for (int x = 1; x < NX - 1; x++) {
            next[idx(n-1, x, NX)] = calculate(mid[idx(n-2, x, NX)], down[idx(0, x, NX)],
                                             mid[idx(n-1, x-1, NX)], mid[idx(n-1, x+1, NX)]);
        }
    }
}

void stencil(int nt, int n, int NX, int NY, double *current, double *next, double **U){
    int num_block = std::ceil((NY + n - 1) / n);
    //int id = 0;

    for (int t = 0; t < nt; t++) {
	if (num_block == 1) {
            // CASE 1: N == n. Only one block, no up/down neighbors.
            // The task covers the entire array: (0; N*N)
#pragma oss task out(next[0]) in(current[0])
            // Call kernel, but adjust parameters since it's the only block
            jacobi_kernel(NX, NY, n, false, &next[0], nullptr, &current[0], nullptr);
	}
	else {
#pragma oss task out(next[0]) in(current[0], current[n*NX])
        jacobi_kernel(NX, NY, n, false, &next[0], nullptr, &current[0], &current[n * NX]);
	//id++;
        for (int j = 1; j < num_block - 1; j++) {
#pragma oss task out(next[j * n * NX]) in(current[(j - 1) * n * NX], current[j * n * NX], current[(j + 1) * n * NX])
            jacobi_kernel(NX, NY, n, false, &next[j * n * NX],
                          &current[(j - 1) * n * NX],
                          &current[j * n * NX],
                          &current[(j + 1) * n * NX]);
	  //  id++;
        }
	if(NY % n != 1){
#pragma oss task out(next[(num_block-1) * n * NX]) in(current[(num_block-2) * n * NX], current[(num_block-1) * n * NX])
        	jacobi_kernel(NX, NY, n, true, &next[(num_block-1) * n * NX], &current[(num_block-2) * n * NX],
                      &current[(num_block-1) * n * NX], nullptr);
	//	id++;
	}
	}
        current = U[(t + 1) % 2];
        next = U[t % 2];	
    }
#pragma oss taskwait
}

void address_map(double* phi, std::size_t size){
    long page_size = numa_pagesize(); 
    for (size_t i = 0; i < size; i += page_size) {
    	int node = -1;
    	void* addr = (char*)phi + i; // Pointer arithmetic in bytes

    	if (get_mempolicy(&node, NULL, 0, addr, MPOL_F_NODE | MPOL_F_ADDR) == 0) {
            // This node ID applies to every double within this page range
            printf("Address %p (offset %zu) is on Node %d\n", addr, i/sizeof(double), node);
    	}
    }
}

int main(int argc, char* argv[]){
    int n = -1, nt = -1, NX = -1, NY = -1;

    // Define long options
    static struct option long_options[] = {
            {"n", required_argument, 0, 'p'},
            {"nt", required_argument, 0, 't'},
            {"NX", required_argument, 0, 'x'},
	    {"NY", required_argument, 0, 'y'},
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
    };

    int opt;
    int option_index = 0;
    while ((opt = getopt_long(argc, argv, "p:t:x:y:h", long_options, &option_index)) != -1) {
        switch (opt) {
            case 'p':
                n = std::atoi(optarg);
                break;
            case 't':
                nt = std::atoi(optarg);
                break;
            case 'x':
                NX = std::atoi(optarg);
                break;
	    case 'y':
                NY = std::atoi(optarg);
                break;
            case 'h':
            default:
                std::cout << "Usage: " << argv[0]
                          << " -p <np> -t <nt> -x <nx>\n"
                          << "   or: " << argv[0]
                          << " --np <np> --nt <nt> --nx <nx>\n";
                return 0;
        }
    }

    // Check required values
    if (n < 0 || nt < 0 || NX < 0 || NY < 0) {
        std::cerr << "Error: Missing required parameters.\n";
        std::cerr << "Use --help for usage.\n";
        return 1;
    }

    nanos6_bitmask_t bitmask;
    nanos6_bitmask_set_wildcard(&bitmask, NUMA_ALL);
    size_t numa_nodes = nanos6_count_setbits(&bitmask);
    size_t size = NX * NY * sizeof(double);
    size_t block_size = size/numa_nodes;
    double* U[2];

    U[0] = (double *)nanos6_numa_alloc_block_interleave(size, &bitmask, block_size);
    U[1] = (double *)nanos6_numa_alloc_block_interleave(size, &bitmask, block_size);
    //U[0] = (double *)malloc(NX * NY * sizeof(double));
    //U[1] = (double *)malloc(NX * NY * sizeof(double));

    double* current = U[0];
    double* next = U[1];

    //address_map(current, size);

   // unsigned int n_threads = std::max(atoi(std::getenv("OMP_NUM_THREADS")), 1);
    std::vector<long> sync(RP, 0.0);
    std::vector<long> exec(RP, 0.0);
    std::vector<long> total(RP, 0.0);

//    for (int w = 0; w < 100; w++)
//    {
//        define_uniform_1(current, next, NX, NY, n);
//    }
//#pragma oss taskwait

    for (int i = 0; i < RP; i++) {
        std::cout << "ROI: " << i << std::endl;

        //std::vector<std::chrono::steady_clock::time_point> task_end(n_threads);
        //std::vector<long> task_duration(n_threads, 0.0);
        //std::vector<long> barrier(n_threads, 0.0);
        std::chrono::steady_clock::time_point taskwait;

        //define_uniform(current, 0, N * N);
        define_uniform_1(current, next, NX, NY, n);
#pragma oss taskwait
	//std::fill(next, next + N * N, 0);
        // Set boundaries to 1.0

	for (std::size_t j = 0; j < NX; j++)
        {
            next[j] = current[j] = 1.0;                                // top boundary
            next[(NY-1) * NX + j] = current[(NY-1) * NX + j] = 1.0;        // bottom boundary
        }
        for (std::size_t j = 0; j < NY; j++)
        {
            next[j * NX] = current[j * NX] = 1.0;                        // left boundary
            next[j * NX + (NX-1)] = current[j * NX + (NX-1)] = 1.0;        // right boundary
        }

	//std::cout << "U[0] mapping: " << std::endl;
	//address_map(current, size);
	//std::cout << "U[1] mapping: " << std::endl;
	//address_map(next, size);

	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  
        stencil(nt, n, NX, NY, current, next, U);                  

        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

        std::cout << "Time = " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count()
                  << std::endl;
        total[i] = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count();

//	for (int i = 0; i < N*N; ++i) {
  //          std::cout << U[nt % 2][i] << " ";
    //    }
        //auto task_finish = *std::max_element(task_end.begin(), task_end.end());
        /*for(int i = 0; i < n_threads; i++)
        {
        std::cout << "Thread " << i << " Task exec time: " << std::setprecision (15) << task_duration[i] << std::endl;
        }*/
        //std::cout << "Task exec time: " << std::setprecision (15) <<
        //std::reduce(task_duration.begin(), task_duration.end()) << std::endl;
        //exec[i] = std::reduce(task_duration.begin(), task_duration.end());

        //std::cout << "Task wait time: " << std::setprecision (15) <<
        //std::chrono::duration_cast<std::chrono::nanoseconds>(taskwait-task_finish).count() << std::endl;

        //std::cout << "Barrier waiting time: " << std::setprecision (15) <<
        //std::reduce(barrier.begin(), barrier.end()) << std::endl;
        //sync[i] = std::chrono::duration_cast<std::chrono::nanoseconds>(taskwait-task_finish).count() +
        //       std::reduce(barrier.begin(), barrier.end());
    }

   // print_matrix(NX, NY, next);
    std::cout << "overall: " << std::endl;
    //std::cout << "median of synchronization " << computeMedian(sync) << std::endl;
    //std::cout << "median of task exec " << computeMedian(exec) << std::endl;
    std::cout << "median of total time " << computeMedian(total) << std::endl;
    
    //free(U[0]);
    //free(U[1]);
    nanos6_numa_free(U[0]);
    nanos6_numa_free(U[1]);
}
