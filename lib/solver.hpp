#ifndef SOLVER_H
    #define SOLVER_H

    #include "hungarian.hpp"

    class edge {
        public:
            int src;
            int dest;
            int weight;
            edge(int source, int destination, int edge_weight): src{source},dest{destination},weight{edge_weight} {}
            bool operator==(const edge &rhs) const {return this->src == rhs.src;}
    };

    /* All the information necessary about the current node in the enumeration tree */
    struct sop_state {
        vector<int> current_path;
        int current_cost; //sum cost of the current_path
        vector<int> taken_arr; //a vector of bits for each node in the graph, 1 for taken, 0 for not taken
        vector<int> depCnt; //for each node, the number of unsatisfied dependencies before it is ready

        Hungarian hungarian_solver;

        pair<boost::dynamic_bitset<>,int> hist_table_key;
        HistoryNode* cur_parent_hisnode = NULL;

        // int initial_depth = 0; //depth at which enumeration began once GPQ was initially filled
        // int suffix_cost = 0;
        // ////////////////////
        // int originate = -1;
        // int load_info = -2;
        // ///////////////////
        // unsigned long long current_node_value = -1; //the portion out of ULLONG_MAX of the working tree that is under this node (the partial path represented by this state)
    };

    /* A single B&B solver assigned one to each thread to process possible paths. */
    class solver {
        private:
            int thread_id = -1;         //a number 0 through (thread_total - 1) identifying this thread

            int instance_size = -1;     //number of nodes in the graph, including artificial starting and ending nodes ("real" nodes would be instance_size - 2)
            sop_state problem_state;    //this thread's current state

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
            void solve_parallel(int thread_num, int pool_size); 

            /* Recursive function that each thread runs to process its assigned spaces of the enumeration tree, checking one node and then its children and their children, etc. */
            void enumerate();

            /* Initialize local pools. */
            //void local_pool_config(float precedence_density);
        public:
            /* Takes config information and defines all runtime parameters from those strings. */
            void assign_parameter(vector<string> setting);
            /* Primary function that initializes and begins the solver. */
            void solve(string filename, int thread_num);
    };

#endif
