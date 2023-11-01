#ifndef SOLVER_H
    #define SOLVER_H

    #include "hungarian.hpp"
    // #include "active_tree.hpp"


    // there may be room for OS / hardware specific preformance enchancements 
// but it seems unlikely to provide much benefit
// does not prevent starvation but I doubt it will impact preformance
    class spin_lock { 
            std::atomic_flag locked = ATOMIC_FLAG_INIT ;
        public:
            void lock() {
                while (locked.test_and_set(std::memory_order_acquire)) { ; }
            }
            void unlock() {
                locked.clear(std::memory_order_release);
            }
    };
    
    //used to pick the next thread to be stolen from
    class workstealing_targeter {
        public:
            workstealing_targeter(int size){
                size = size;
            }
            int get(){
                lock.lock();
                int return_value = cur;
                cur = (cur + 1) % size;
                lock.unlock();
                return return_value;
            }
        private:
            int cur = 0;
            int size;
            spin_lock lock;
    };

    /* An edge in the cost graph, from src to dst with cost weight. */
    class edge {
        public:
            int src;
            int dst;
            int weight;
            edge(int source, int destination, int edge_weight): src{source},dst{destination},weight{edge_weight} {}
            bool operator==(const edge &rhs) const {return this->src == rhs.src;}
    };

    /* A node in the graph, containing the information needed to consider this node compared to other children at each level. */
    class node {
        public:
            int n = -1; //node number
            int lb = -1; //lower bound cost of a complete path including the current prefix followed by this node
            int nc = -1; //cost of the edge from the current node to this child
            int partial_cost = -1; //total cost of path to parent
            // unsigned long long current_node_value = -1; //the portion out of ULLONG_MAX of the working tree that is under this node
            // bool node_invalid = false; //investigate
            // bool pushed = false;
            // HistoryNode* his_entry = NULL;
            // Active_Node* act_entry = NULL;

            node(int node_number): n{node_number} {}
            node(int node_number, int lower_bound): n{node_number},lb{lower_bound} {}
    };

    /* Entry in local or global pools containing all the information about a position in the enumeration tree needed to regenerate an sop_state. */
    class path_node {
        public:
            vector<int> sequence;       //the partial path represented by this node in the enumeration tree
            int lower_bound = -1;       //the lower bound cost of a complete path beginning with sequence
            int origin_node = -1;       //the first node in this path, after the virtual starting node, used to evenly distribute threads between subspaces in solve_parallel
            // int parent_lv = -1;
            // bool* invalid_ptr = NULL; //investigate
            // bool deprecated = false; //For Thread Stopping. if this node exists in a redundant subspace and so does not need to be processed
            // unsigned long long current_node_value = -1; //the portion out of ULLONG_MAX of the working tree that is under this node

            //Active_Path partial_active_path;
            //HistoryNode* root_his_node;

            path_node () {}
            path_node (vector<int> partial_path, int lb, int origin) {
                sequence = partial_path;
                lower_bound = lb;
                origin_node = origin;
            }
            // path_node(vector<int> sequence_src,int originate_src, int load_info_src, 
            //             int best_costrecord_src, unsigned long long value, Active_Path temp_active_path, 
            //             HistoryNode* temp_root_hisnode) {
            //     sequence = sequence_src;
            //     lower_bound = load_info_src;
            //     origin_node = originate_src;
            //     parent_lv = best_costrecord_src;
            //     current_node_value = value;
            //     partial_active_path = temp_active_path;
            //     root_his_node = temp_root_hisnode;
            // }
    };

    class local_pool {
        public:
            /*Grabs a node from the shallowest / zero pool*/
            bool pop_from_zero_list(int thread_number, path_node &result_node){
                if (pools[thread_number].size() == 0)
                    return false;
                locks[thread_number].lock();

                if (pools[thread_number].size() == 0 || pools[thread_number].front().empty()){
                    locks[thread_number].unlock();
                    return false;
                }

                result_node = pools[thread_number].front().front();
                pools[thread_number].front().pop_front();

                if (pools[thread_number].front().empty()){
                    pools[thread_number].pop_front();
                }

                locks[thread_number].unlock();
                return true;
            };
            /*Grabs a node from the deepest / active pool*/
            bool pop_from_active_list(int thread_number, path_node &result_node){
                if (pools[thread_number].size() == 0)
                    return false;
                if (pools[thread_number].size() == 1){
                    locks[thread_number].lock();
                }

                if (pools[thread_number].size() == 0 || pools[thread_number].back().empty()){
                    if (pools[thread_number].size() == 1)
                    {
                        locks[thread_number].unlock();
                    }
                    return false;
                }

                result_node = pools[thread_number].back().front();
                pools[thread_number].back().pop_front();

                if (pools[thread_number].size() == 1){
                    locks[thread_number].unlock();
                }
                return true;
            };
            /*Pushes new list to the back of the local pool*/
            void push_list(int thread_number, std::deque<path_node> list){
                locks[thread_number].lock();

                pools[thread_number].push_back(list);

                locks[thread_number].unlock();
            };
            /*removes */
            void pop_active_list(int thread_number){
                locks[thread_number].lock();
                pools[thread_number].pop_back();
                locks[thread_number].unlock();
            };
            bool out_of_work(int thread_number){
                return pools[thread_number].size() == 0;
            };
            local_pool(int thread_count){
                locks = vector<spin_lock>(thread_count);
                pools = vector<std::deque<std::deque<path_node>>>(thread_count);
            }
        private:
            std::vector<spin_lock> locks;
            std::vector<std::deque<std::deque<path_node>>> pools;
    };

    /* All the information necessary about the current node in the enumeration tree.
        Since sop_state is a struct, it is always passed by value. */
    struct sop_state {
        vector<int> current_path; //the current partial path being considered
        int current_cost; //sum cost of the current_path
        /*  for each node, whether it is included in the current path
            This typing is required because of the special behavior of std::vector<bool> that is implemented as a bitset instead of an iterable C array. */
        boost::container::vector<bool> taken_arr; 
        vector<int> depCnt; //for each node, the number of unsatisfied dependencies before it is ready

        //previously called load_info
        int lower_bound = -1; //the lower bound cost of a complete solution beginning with current_path
        int origin_node = -1; //the first node in this path, after the virtual starting node, used to evenly distribute threads between subspaces in solve_parallel
        
        Hungarian hungarian_solver;

        pair<boost::dynamic_bitset<>,int> history_key;
        HistoryNode* cur_parent_hisnode = NULL;
        
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
            // vector<recur_state> recur_stack;

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
            size_t transitive_closure(vector<vector<int>>& isucc_graph);
            /* Performs the Nearest Neighbor Heuristic for the SOP to find an initial solution. */
            vector<int> nearest_neighbor(vector<int>* partial_solution);
            /* Sort the cost graph in descending order of weight. Required for nearest neighbor heuristic. */
            void sort_weight(vector<vector<edge>>& graph);
            /* Find the highest edge weight in the entire cost graph. Required for Hungarian algorithm. */
            int get_maxedgeweight();
            /* Generate the cost matrix that the Hungarian algorithm uses. */
            vector<vector<int>> get_cost_matrix(int max_edge_weight);

            /* Repeatedly run LKH routine. */
            void run_lkh();
            /* Called from solver::solve, divides work among global pool and each thread, and begins the threads with calls to solver::enumerate. */
            void solve_parallel(); 
            /* Returns true if any sop_state in the container has a depth different than any other, false otherwise. Used for initial splitting in solve_parallel. */
            bool split_level_check(deque<sop_state>* solver_container);

            /* Recursive function that each thread runs to process its assigned spaces of the enumeration tree, checking one node and then its children and their children, etc. */
            void enumerate();
            /* */
            bool enumeration_pre_check(path_node& active_node);

            /* */
            int dynamic_hungarian(int src, int dest);

            /* */
            bool history_utilization(pair<boost::dynamic_bitset<>,int>& key,int* lowerbound,bool* found,HistoryNode** entry, int cost);
            /* */
            void push_to_historytable(pair<boost::dynamic_bitset<>,int>& key,int lower_bound,HistoryNode** entry,bool backtracked);

            /* Build an sop_state based off the information in a path_node. */
            sop_state generate_solver_state(path_node& subproblem);
            /* Build a hungarian solver state based upon the problem_state. Used in generate_solver_state. */
            void regenerate_hungstate();

            /* */
            void solver::workload_request();
            /* */
            path_node solver::workstealing();

            /* Initialize local pools. */
            //void local_pool_config(float precedence_density);
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
