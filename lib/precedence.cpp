#include "precedence.hpp"
#include <iostream>
#include <fstream>
#include <stack>
#include <sstream>
#include <cmath>
#include <chrono>
#include <math.h>
#include <limits>
#include <algorithm>
#include <cstring>
#include <numeric>
#include <functional>
#include <vector>
#include <deque>
#include <boost/dynamic_bitset.hpp>

using namespace std;

transformer::transformer(unsigned node_count) {
    size = node_count;
    transitive_closure_map = vector<boost::dynamic_bitset<>>(node_count);
    sucpre_list = vector<pair<boost::dynamic_bitset<>,boost::dynamic_bitset<>>>(node_count);
    for (unsigned i = 0; i < transitive_closure_map.size(); i++) {
        transitive_closure_map[i] = boost::dynamic_bitset<>(node_count);
        sucpre_list[i].first = boost::dynamic_bitset<>(node_count);
        sucpre_list[i].second = boost::dynamic_bitset<>(node_count);
    }
}

void transformer::transform(vector<vector<edge>>& cost_graph, vector<vector<int>>& ipred_graph, vector<vector<int>>& isucc_graph, vector<vector<edge>>& hung_graph) {
    transitive_closure(ipred_graph,isucc_graph);
    generate_superior_pair(cost_graph,ipred_graph,isucc_graph);
    unsigned num_artificaleg = 0;

    while(!superior_pair.empty()) {
        cout << "edge <" << superior_pair.back().first << "->" << superior_pair.back().second << "> added" << endl;
        num_artificaleg++;
        ipred_graph[superior_pair.back().second].push_back(superior_pair.back().first);
        isucc_graph[superior_pair.back().first].push_back(superior_pair.back().second);
        hung_graph[superior_pair.back().second][superior_pair.back().first] = edge(superior_pair.back().second,superior_pair.back().first,-1);
        /*
        if (cycle_detection(isucc_graph)) {
            cout << "error: cycle detected" << endl;
            exit(1);
        }
        */
        independent_pair.clear();
        superior_pair.clear();
        for (unsigned i = 0; i < transitive_closure_map.size(); i++) transitive_closure_map[i].reset();
        transitive_closure(ipred_graph,isucc_graph);
        generate_superior_pair(cost_graph,ipred_graph,isucc_graph);
    }

    cout << "Total number of AE is " << num_artificaleg << endl;
}

bool transformer::cycle_detection(vector<vector<int>>& isucc_graph) {
    deque<int> dfs_queue;
    boost::dynamic_bitset<> visited_node = boost::dynamic_bitset<>(size);
    boost::dynamic_bitset<> finished_node = boost::dynamic_bitset<>(size);
    int selected_node = 0;

    for (unsigned i = 0; i < size; i++) {
        dfs_queue.push_back(i);
        while (!dfs_queue.empty()) {
            selected_node = dfs_queue.back();

            if (visited_node[selected_node]) {
                finished_node[selected_node] = true;
                dfs_queue.pop_back();
            }
            else {
                visited_node[selected_node] = true;
                for (int node : isucc_graph[selected_node]) {
                    if (!visited_node[node] && !finished_node[node]) {
                        dfs_queue.push_back(node);
                    }
                    else if (visited_node[node] && !finished_node[node]) {
                        cout << "detected visited node is " << selected_node << "-->" << node << endl;
                        for (unsigned i = 0; i < dfs_queue.size(); i++) {
                            cout << dfs_queue[i] << ",";
                        }
                        cout << endl;
                        return true;
                    }
                }
            }
        }
        visited_node.reset();
        finished_node.reset();
    }
    return false;
}


void transformer::transitive_closure(vector<vector<int>>& ipred_graph, vector<vector<int>>& isucc_graph) {
    deque<int> dfs_queue;
    boost::dynamic_bitset<> taken_arr = boost::dynamic_bitset<>(size);
    boost::dynamic_bitset<> finished_arr = boost::dynamic_bitset<>(size);
    int selected_node = 0;

    for (unsigned i = 0; i < size; i++) {
        //cout << "current node is " << i << endl;
        transitive_closure_map[i][i] = true;
        dfs_queue.push_back(i);
        while (!dfs_queue.empty()) {
            selected_node = dfs_queue.back();
            dfs_queue.pop_back();
            for (int node : isucc_graph[selected_node]) {
                if (!taken_arr[node]) {
                    taken_arr[node] = true;
                    dfs_queue.push_back(node);
                }
                if (!transitive_closure_map[i][node]) transitive_closure_map[i][node] = true;
                if (!transitive_closure_map[node][i]) transitive_closure_map[node][i] = true;
            }
        }
        taken_arr.reset();
    }

    for (unsigned i = 0; i < size; i++) {
        for (unsigned k = 0; k < size; k++) {
            if (!transitive_closure_map[i][k] && !transitive_closure_map[k][i]) {
                independent_pair.push_back(pair<int,int>(i,k));
            }
        }
    }
    //print_transitive_map();
    return;
}

void transformer::vaild_indepnodes(vector<vector<edge>>& cost_graph,const int src, const int dest, bool pred, int& max_diff, int& count) {
    if (pred) {
        for (unsigned i = 0; i < transitive_closure_map[dest].size(); i++) {
            if (!transitive_closure_map[dest][i] && i != src) {
                count++;
                int diff = cost_graph[i][src].weight - cost_graph[i][dest].weight;
                if (diff > max_diff) max_diff = diff;
            }
        }
    }
    else {
        for (unsigned i = 0; i < transitive_closure_map[dest].size(); i++) {
            if (!transitive_closure_map[dest][i] && i != src) {
                count++;
                int diff = cost_graph[src][i].weight - cost_graph[dest][i].weight;
                if (diff > max_diff) max_diff = diff;
            }
        }
    }
    return;
}

void transformer::generate_superior_pair(vector<vector<edge>>& cost_graph, vector<vector<int>>& ipred_graph, vector<vector<int>>& isucc_graph) {
    for (auto indep_pair : independent_pair) {
        int max_diffA = INT_MIN;
        int max_diffB = INT_MIN;
        int max_diffC = INT_MIN;
        int max_diffD = INT_MIN;
        int special_diff = INT_MIN;

        int A1 = 0;
        int A2 = 0;
        int A3 = 0;
        int A4 = 0;

        special_diff = cost_graph[indep_pair.second][indep_pair.first].weight - cost_graph[indep_pair.first][indep_pair.second].weight;

        for (int node : ipred_graph[indep_pair.first]) {
            A1++;
            int diff = cost_graph[node][indep_pair.second].weight - cost_graph[node][indep_pair.first].weight;
            if (diff > max_diffA) max_diffA = diff;
        }
        vaild_indepnodes(cost_graph,indep_pair.second,indep_pair.first,true,max_diffA,A1);

        for (int node : isucc_graph[indep_pair.first]) {
            A2++;
            int diff = cost_graph[indep_pair.second][node].weight - cost_graph[indep_pair.first][node].weight;
            if (diff > max_diffB) max_diffB = diff;
        }
        vaild_indepnodes(cost_graph,indep_pair.second,indep_pair.first,false,max_diffB,A2);

        for (int node : ipred_graph[indep_pair.second]) {
            A3++;
            int diff = cost_graph[node][indep_pair.first].weight - cost_graph[node][indep_pair.second].weight;
            if (diff > max_diffC) max_diffC = diff;
        }
        vaild_indepnodes(cost_graph,indep_pair.first,indep_pair.second,true,max_diffC,A3);

        for (int node : isucc_graph[indep_pair.second]) {
            A4++;
            int diff = cost_graph[indep_pair.first][node].weight - cost_graph[indep_pair.second][node].weight;
            if (diff > max_diffD) max_diffD = diff;
        }
        vaild_indepnodes(cost_graph,indep_pair.first,indep_pair.second,false,max_diffD,A4);

        if ((max(special_diff,max(max_diffB,max_diffC)) + max_diffA + max_diffD) > 0) {
            //cout << "max diff A,B,C,D = " << max_diffA << "," << max_diffB << "," << max_diffC << "," << max_diffD << endl;
            continue;
        }
        cout << "..................................";
        cout << "max diff A,B,C,D = " << max_diffA << "," << max_diffB << "," << max_diffC << "," << max_diffD << endl;
        cout << "neighbor A,B,C,D = " << A1 << "," << A2 << "," << A3 << "," << A4 << endl;
        cout << "..................................";
        
        superior_pair.push_back(indep_pair);
    }
    independent_pair.clear();
}

void transformer::print_transitive_map() {
    
    for (auto sub_map : transitive_closure_map) {
        for (unsigned i = 0; i < sub_map.size(); i++) {
            cout << sub_map[i] << ",";
        }
        cout << endl;
    }
    
    /*
    cout << "Independent pair array contains " << endl;
    for (auto pair : independent_pair) {
        cout << "[" << pair.first << "," << pair.second << "],";
    }
    cout << "Printing finished" << endl;
    cout << endl;
    */
    return;
}

bool transformer::pred_valid(const int i, const int k) {
    for (unsigned j = 0; j < size; j++) {
        if (sucpre_list[i].first[j] && sucpre_list[i].first[j] != sucpre_list[k].first[j]) return false;
    }
    return true;
}

bool transformer::succ_valid(const int i, const int k) {
    for (unsigned j = 0; j < size; j++) {
        if (sucpre_list[k].second[j] && sucpre_list[k].second[j] != sucpre_list[i].second[j]) return false;
    }
    return true;
}