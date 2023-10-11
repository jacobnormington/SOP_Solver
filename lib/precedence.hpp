#ifndef PRECEDENCE_H
#define PRECEDENCE_H

#include <string>
#include <vector>
#include <climits>
#include <queue>
#include <deque>
#include <chrono>
#include <boost/dynamic_bitset.hpp>
#include "solver.hpp"


class transformer {
    private:
        unsigned size = 0;
        vector<boost::dynamic_bitset<>> transitive_closure_map;
        vector<pair<unsigned,unsigned>> independent_pair;
        vector<pair<unsigned,unsigned>> superior_pair;
        vector<pair<boost::dynamic_bitset<>,boost::dynamic_bitset<>>> sucpre_list;
    public:
        transformer(unsigned node_count);
        void generate_superior_pair(vector<vector<edge>>& cost_graph, vector<vector<int>>& ipred_graph, vector<vector<int>>& isucc_graph);
        void transform(vector<vector<edge>>& cost_graph, vector<vector<int>>& ipred_graph, vector<vector<int>>& isucc_graph, vector<vector<edge>>& hung_graph);
        void vaild_indepnodes(vector<vector<edge>>& cost_graph,const int src, const int dest, bool pred, int& max_diff, int& count);
        void transitive_closure(vector<vector<int>>& ipred_graph, vector<vector<int>>& isucc_graph);
        void print_transitive_map();
        bool pred_valid(const int i, const int k);
        bool succ_valid(const int i, const int k);
        bool cycle_detection(vector<vector<int>>& isucc_graph);
};

#endif