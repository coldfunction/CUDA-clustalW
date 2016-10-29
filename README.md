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
		number of sequence	sequence length
	(1)		 100			  97
	(2)		 100			 498
	(3)		 100     		1002
	(4)		 100			1523
	(5)	   	1000			  97
	(6)		1000			 498
	(7)		1000			1002
	(8)		1000			1523

4. How to test?
	a. CPU:
		$ ./run_cpu.sh [number of sequence] [sequence length]
		Example:
		$ ./run_cpu.sh 100 498	

	b. GPU:
		$ ./run_gpu.sh [number of sequence] [sequence length]
		Example:
		$ ./run_gpu.sh 100 498	
		

	


