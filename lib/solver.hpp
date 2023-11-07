#ifndef SOLVER_H
    #define SOLVER_H

    #include <iostream>
    #include <fstream>
    #include <sstream>
    #include <iomanip>
    #include <vector>
    //#include <queue>
    #include <deque>
    #include <stack>
    //#include <unordered_map>
    #include <unordered_set>
    #include <math.h>
    #include <cmath>
    #include <limits>
    #include <algorithm>
    #include <cstring>
    #include <numeric>
    #include <functional>
    #include <chrono>
    #include <thread>
    #include <mutex>
    #include <condition_variable>

    #include <boost/container/vector.hpp>
    #include <boost/dynamic_bitset.hpp>

    extern "C" {
        #include "LKH/LKHmain.h"
    }


    #include "history_table.hpp"
    #include "local_pool.hpp"
    #include "graph.hpp"
    #include "hungarian.hpp"
    // #include "active_tree.hpp"
    //#include "precedence.hpp"

    /* A Thread Stop request. */
    struct request_packet {
        int target_last_node;
        int target_depth;
        int target_prefix_cost;
        int target_thread;
        boost::dynamic_bitset<> key;
    };

    /* All the information necessary about the current node in the enumeration tree.
        Since sop_state is a struct, it is always passed by value. */
    struct sop_state {
        std::vector<int> current_path; //the current partial path being considered
        int current_cost; //sum cost of the current_path
        /*  for each node, whether it is included in the current path
            This typing is required because of the special behavior of std::vector<bool> that is implemented as a bitset instead of an iterable C array. */
        boost::container::vector<bool> taken_arr; 
        std::vector<int> depCnt; //for each node, the number of unsatisfied dependencies before it is ready

        //previously called load_info
        int lower_bound = -1; //the lower bound cost of a complete solution beginning with current_path
        int origin_node = -1; //the first node in this path, after the virtual starting node, used to evenly distribute threads between subspaces in solve_parallel
        
        Hungarian hungarian_solver;

        pair<boost::dynamic_bitset<>,int> history_key;
        HistoryNode* cur_parent_hisnode = NULL;
        
        int enumeration_depth = 0; // the depth of the enumeration recursion stack

        // int initial_depth = 0; //depth at which enumeration began once GPQ was initially filled
        // int suffix_cost = 0;
        // unsigned long long current_node_value = -1; //the portion out of ULLONG_MAX of the working tree that is under this node (the partial path represented by this state)
    };

    /* A single B&B solver assigned one to each thread to process possible paths. */
    class solver {
        private:
            int thread_id = -1;         //a number 0 through (thread_total - 1) identifying this thread

            int instance_size = -1;     //number of nodes in the graph, including virtual starting and ending nodes ("real" nodes would be instance_size - 2)
            sop_state problem_state;    //this thread's current state
            // sop_state back_up_state;
            // HistoryNode* current_hisnode;

            // deque<instrct_node> wrksteal_pool;
            // deque<instrct_node> *local_pool = NULL;
            // std::vector<recur_state> recur_stack;

            // Active_Allocator Allocator;
            // Active_Path cur_active_tree;
            // bool abandon_work = false;
            // bool abandon_share = false;
            // bool grabbed = false;
            // int restart_group_id = -1;
            // int mg_id = -1;
            // bool speed_search = false;
            // int lb_curlv = INT_MAX;

            // //Restart
            // int concentrate_lv = 0;

            // //Thread Stopping
            // int stop_depth = -1;
            // int last_node = -1;
            // bool stop_init = false; //INVESTIGATE; might be whether this thread has ever been stopped before

            /* Build graph based on .sop input file specified in filename. */
            void retrieve_input();
            /* Transforms dependency and Hungarian graphs, adding redundant edges from grandparents, great grandparents, etc., and initializes in_degree. */
            void transitive_redundancy();
            /* Computes the transitive closure of the graph. Used for calculating precedence density. */
            size_t transitive_closure(std::vector<std::vector<int>>& isucc_graph);
            /* Performs the Nearest Neighbor Heuristic for the SOP to find an initial solution. */
            vector<int> nearest_neighbor(vector<int>* partial_solution);
            /* Sort the cost graph in descending order of weight. Required for nearest neighbor heuristic. */
            void sort_weight(std::vector<std::vector<edge>>& graph);
            /* Find the highest edge weight in the entire cost graph. Required for Hungarian algorithm. */
            int get_maxedgeweight();
            /* Generate the cost matrix that the Hungarian algorithm uses. */
            vector<vector<int>> get_cost_matrix(int max_edge_weight);

            /* Called from solver::solve, divides work among global pool and each thread, and begins the threads with calls to solver::enumerate. */
            void solve_parallel(); 
            /* Returns true if any sop_state in the container has a depth different than any other, false otherwise. Used for initial splitting in solve_parallel. */
            bool split_level_check(deque<sop_state>* solver_container);

            /* Recursive function that each thread runs to process its assigned spaces of the enumeration tree, checking one node and then its children and their children, etc. */
            void enumerate();
            /* */
            bool enumeration_pre_check(path_node& active_node);
            /*called when pruning a node in enumerate*/
            void prune(int source_node, int taken_node);

            /* Computes a dynamic lower bound based on the previous path with this node added, using the MCPM relaxation. 
                Contains the fix and undue calls internally. 
                src - the number of this node's parent
                dst - the number of the node to be added
                Return - the lower bound computed */
            int dynamic_hungarian(int src, int dst);

            /* Search the history table for previously processed similar paths, and compares the current path to that entry, if found. 
                key - the history key corresponding to the current partial path
                lowerbound - a return variable, which contains the lower bound found in the history table, if a corresponding entry was found
                found - a return variable, true if an entry already existed, false otherwise
                entry - a return variable, a pointer to the history node corresponding to this path
                cost - the cost of the current path
                Return - true if this node still needs to be processed, false if it should be pruned */
            bool history_utilization(Key& key, int cost, int* lowerbound, bool* found, HistoryNode** entry);
            /* Add a new entry to the history table. 
                key - the history key corresponding to the partial path this entry represents
                lower_bound - the lower bound cost of a complete solution beginning with this path
                entry - a return variable, holds a pointer to the entry created, unless NULL is passed
                backtracked - if the subtree under this node has already been fully explored */
            void push_to_history_table(Key& key,int lower_bound,HistoryNode** entry,bool backtracked);

            /* */
            void solver::workload_request();
            /* */
            path_node solver::workstealing();
            /* Build an sop_state based off the information in a path_node. */
            sop_state generate_solver_state(path_node& subproblem);
            /* Build a hungarian solver state based upon the problem_state. Used in generate_solver_state. */
            void regenerate_hungstate();
            //assign_workload???
            //push_to_pool???

            //check_request_buffer
            //thread_stop

            /* Repeatedly run LKH routine. */
            void run_lkh();

        public:
            /* Takes config information and defines all runtime parameters from those strings. */
            void assign_parameter(vector<string> setting);
            /* Primary function that initializes and begins the solver. */
            void solve(string f_name, int thread_num);
    };

    /* These 64 byte structs are necessary for some shared resources in order to reduce cache coherency problems. 
        Since 64 bytes is the size of one cache line, this ensures that each item sits on its own cache line, so
        that the same line won't be accessed and modified by different threads. */
    struct int_64 {
        int val = 0;
        bool explored = false;
        bool padding[59];
    };
    struct bool_64 {
        bool val = false;
        bool padding[63];
    };
    struct unsigned_long_64 {
        unsigned long val = 0;
        bool padding[56];
    };
    struct mutex_64 {
        mutex lck;
        bool padding[24];
    };
    struct lptr_64 {
        deque<path_node> *local_pool = NULL;
        bool padding[56];
    };

#endif
