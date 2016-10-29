# CUDA-clustalW
---

  In computational biology, sequence alignment is of priority concern and many methods have been developed to solve sequence alignment-related problems for biological applicatons.
ClustalW is a progressive multiple sequence alignment tool to align a set of sequences by repeatedly aligning pairs of sequences and previously generated alignments. Several algorithms or tools 
have been ported on GPUs with CUDA in computational biology, such as MUMmerGPU, CUDA-MEME, CUDA-BLASTP, and etc. Liu et al. proposed a tool MSA-CUDA to parallelize all three stages 
of ClustalW v2.0.9 processing pipeline by using inter-task parallelization. CUDA ClustalW v1.0 is a GPU version of ClustalW v2.0.11 which is implemented by using intra-task parallelization and 
Synchronous Diagonal Multiple Threads type. Several optimization methods were designed to improve the performance of CUDA ClustalW v1.0. From the experimental results, the CUDA ClustalW 
v1.0 can achieve about 22x speedups for 1000 sequences with length of 1532 in the distance matrix calculation step by comparing to ClustalW v2.0.11 on single-GPU. For the overall execution time, 
the CUDA ClustalW v1.0 can achieve about 33x speedups by comparing to ClustalW v2.0.11 on two-GPUs.


## Version 1.0.0 (March 2013) 

* Linux x86 64-bit
* CPU/GPU coprocess
* Support Multiple GPUs
* NVIDIA CUDA support
* Base on clustalW 2.0

## How to run

1. source code:  

    All in the CPU and GPU directory.

2. How to compile?  

	compile:
  
		$ ./make.sh
    
	clean:  
  
		$ ./make.sh clean

3. Test data:

	Testing data in the "./test_data" directory,
  
	There are eight test data stored as fasta format:

|      number of sequence      |  sequence length|
|----------|:-------------:|------:|
|  100 | 97 |
|    100   |   498 |
| 100 |   1002 |
|  100|			1523    
|1000|			  97
|1000|			 498
|1000|			1002
|1000|			1523
  
4. How to test?

a. CPU:

    $ ./run_cpu.sh [number of sequence] [sequence length]

 Example:

  	$ ./run_cpu.sh 100 498	

b. GPU:

    $ ./run_gpu.sh [number of sequence] [sequence length]

Example:

    $ ./run_gpu.sh 100 498	
		

	


