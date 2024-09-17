1# SOP_Solver

A parallel Branch and Bound solution to the Sequential Ordering Problem (SOP).

## Background

This code is an implementation of an algorithm developed in the work of Dr. Ghassan Shobaki et. al. in [1] [2] [3] [4] and used as a reference for future developments.

[1] G. Shobaki, J. Jamal. “An Exact Algorithm for the Sequential Ordering Problem and its Application to Switching Energy Minimization in Compilers”. Comput Optim Appl 61, pp. 343–372 (2015). DOI: 10.1007/s10589-015-9725-9.

[2] J. Jamal, G. Shobaki, V. Papapanagiotou, L. M. Gambardella and R. Montemanni. “Solving the Sequential Ordering Problem using Branch and Bound”. 2017 IEEE Symposium Series on Computational Intelligence (SSCI), pp. 1-9. DOI: 10.1109/SSCI.2017.8280805.

[3] T. Gonggiatgul, G. Shobaki, and P. Muyan-Özçelik. “A Parallel Branch and
Bound Algorithm with History-Based Domination and its Application to the
Sequential Ordering Problem”. J. Parallel and Distrib. Comput. 172, pp. 131–143 (2023). DOI: 10.1016/j.jpdc.2022.10.007

[4] G. Shobaki, T. Gonggiatgul, J. Normington and P. Muyan-Ozcelik. “Combining a Parallel Branch-and-Bound Algorithm with a Strong Heuristic to Solve the Sequential Ordering Problem”. In Proc. Third International Workshop on Parallel and Distributed Algorithms for Decision Sciences (PDADS), August 2023. DOI: 10.1145/3605731.3608929.

## References and Dependencies

A reimplementation of https://github.com/tasuka98/PSOP_Solver based upon the referenced work.

The dynamic Hungarian library code is referenced from https://github.com/rod409/SOP, but does not need to be downloaded seperately.

The LKH library code is taken and modified from the [LKH3 Home Page](http://webhotel4.ruc.dk/~keld/research/LKH-3/), but does not need to be downloaded seperately.

The boost library is required before compilation. It can be installed via `sudo apt-get install libboost-all-dev`.

## Usage

### Config File

A config file is supplied with a series of runtime parameters to alter the solver's behavior, including the enabling or disabling of features such as work stealing or thread stopping. There is a seperate config file for each test suite (soplib vs. tsplib) because different parameter values were found to be ideal between the two suites. Following is a summary of each parameter:

- <b>Time_Limit</b>: the maximum amount of time the solver is allowed to use, in seconds. Default 1 hour.
- <b>Global_Pool_Size</b>: the minimum size of the global pool before enumeration can begin. Default 32.
- <b>Restrict_Per</b>: the percentage (0,1) of total available memory that the history table can use before the solver should use the more restrictive algorithm for adding history table entries. Default 0.8.
- <b>History_Depth</b>: when the more restrictive algorithm for adding history table entries is used, always add solutions of this length or less. Default 100.
- <b>Enable (Work Stealing)</b>: whether to enable the work stealing routine that assigns idle threads new work taken from the local pool of another thread. Default 1 (true).
- <b>Enable (Thread Stopping)</b>: whether to enable the thread stopping routine that finds threads working in redundant subspaces based on history table entries and assigns them new, useful work. Default 1 (true).
- <b>Enable (LKH)</b>: whether to allot one thread to run the LKH heuristic instead of branch and bound enumeration. Default 1 (true).
- <b>Enable (Progress Estimation)</b>: whether to enable the progress tracking heuristic, useful for hard instances to approximate how far the solver was able to get before timeout. Default 0 (false).

### Running a Specific Instance

You can run the solver using the following commands, after the boost library is successfully installed:

    make
    ./sop_solver <Instance Location> <Number of Threads> <Config File>

Example: `./sop_solver soplib/R.200.100.1.sop 32 soplib_config.txt`

### Use of the Script to Run Entire Test Suite

Script files are supplied to run many tests in sequence. There is a seperate script for each test suite (soplib vs. tsplib) in order to accomodate their different structure. The script has a number of parameters that can be set to run only a specific subset of instances from that suite, and compiles all of the output into a file `outfile_<type>_raw.log` with a denser summary file `outfile_<type>.log` that includes only the name of each instance, its running time, and the cost of the best solution found.

- <b>startInstance</b>: the first instance to process, skipping any lexicographically earlier; if empty, begins at the first instance. Default "" (start at beginning).
- <b>endInstance</b>: the last instance to process, skipping any lexicographically later; if empty, continues through the entire suite. Default "" (continue until the last).
- <b>dataset</b>: a parameter to determine what portion of the test suite should be processed. Default 0 (entire suite).
  - <i>dataset =&nbsp; 1</i> &nbsp;&nbsp;->&nbsp;&nbsp; exclude hard instances
  - <i>dataset =&nbsp; 0</i> &nbsp;&nbsp;->&nbsp;&nbsp; include all instances
  - <i>dataset = -1</i> &nbsp;&nbsp;->&nbsp;&nbsp; include only hard instances
- <b>num_threads</b>: the total number of threads for the solver to use. Default 32, as any higher resulted in worse results.

You can use the following code to run an entire family of tests:

    make
    bash <Script Name>

Example: `bash run_soplib.sh`

## Multiple History Table

### Overview

To manage the history table efficiently, we use two key variables in the configuration file:

- `Number_of_Buckets`: This determines how many buckets the history table is split into.
- `Bucket_size`: This specifies the depth each bucket should carry. For example, if the size is 25, entries of depth 0 to 25 go into the first bucket, 26 to 50 into the second, and 51+ into the last, assuming `Number_of_Buckets` is set to 3.

### Improving Memory Management with Multiple History Tables

Let's consider a sequence of events happening sequentially, with `Restrict_Per` in the configuration set to 0.9.

- Initialization: Since Number_Of_Buckets is set to 3, we will stop inserting into the last bucket when the memory consumption level reaches 0.7. This is calculated using the formula: `Restrict_Per` - (`Number_of_buckets` - 1) \* 0.1.
- Memory Consumption Thresholds:
  - 0.7: Stop inserting into the last bucket.
  - 0.8: Stop inserting into the second last bucket.
  - 0.9: Free the memory held by the last bucket, then by the second last bucket, and finally stop inserting into the first bucket.

### Example

#### Initial Setup:

- Number_of_Buckets = 3
- Bucket_size = 25
- Restrict_Per = 0.9

#### Memory Management Steps:

- Memory Consumption Reaches 0.7: Stop inserting into the last bucket.
- Memory Consumption Reaches 0.8: Stop inserting into the second last bucket.
- Memory Consumption Reaches 0.9:
  - Free the memory held by the last bucket.
  - Free the memory held by the second last bucket.
  - Stop inserting into the first bucket.

This approach ensures that memory is managed efficiently by gradually restricting insertions and freeing memory as consumption levels increase.

## Version

- `version_1.0.0`: Branch before merging the history changes
- `main`: Latest branch with all the history changes merged into it
