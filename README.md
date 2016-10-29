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
		

	


