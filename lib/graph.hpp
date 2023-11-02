/* Basic graph classes, including the more complex path_node. */

#ifndef GRAPH_H
    #define GRAPH_H

    #include <vector>

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
        //TODO: do we even need minimal nodes anymore? It seems that we are only ever using path_nodes
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
            std::vector<int> sequence;       //the partial path represented by this node in the enumeration tree
            int lower_bound = -1;       //the lower bound cost of a complete path beginning with sequence
            int origin_node = -1;       //the first node in this path, after the virtual starting node, used to evenly distribute threads between subspaces in solve_parallel
            // int parent_lv = -1;
            // bool* invalid_ptr = NULL; //investigate
            // bool deprecated = false; //For Thread Stopping. if this node exists in a redundant subspace and so does not need to be processed
            // unsigned long long current_node_value = -1; //the portion out of ULLONG_MAX of the working tree that is under this node

            //Active_Path partial_active_path;
            //HistoryNode* root_his_node;

            path_node () {}
            path_node (std::vector<int> partial_path, int lb, int origin) {
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

#endif