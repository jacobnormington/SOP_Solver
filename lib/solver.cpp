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
#include <thread>
#include <queue>
#include <deque>
#include <mutex>
#include <bits/stdc++.h>
#include <unordered_map>
#include <condition_variable>

extern "C" {
    #include "LKH/LKHmain.h"
}

#include "solver.hpp"
#include "hash.hpp"
#include "precedence.hpp"

#define TABLE_SIZE 541065431


//////Runtime Parameters (Read Only)/////
    //from command line arguments
    static string filename;                         //name of the sop input file
    static int thread_total = 0;                    //number of threads to use (not counting the LKH thread)

    //from config file
    static int t_limit = 0;                         //time limit, in seconds
    static int global_pool_size = 0;                //Minimum size of the global pool at before enumeration begins
    static int local_depth = 0;                     //Minimum size to maintain in the local pool
    static float inhis_mem_limit = -1;              //0-1, the percentage of memory usage beyond which new entries shouldn't be added to the history table
    static unsigned int inhis_depth = -1;           //after inhis_mem_limit exceeded, will still add an entry if the current depth is less than inhis_depth
    // int exploitation_per;                        //percent of threads that should be devoted to searching already promising subspaces in thread restart, while 1 - exploitation_per percent are devoted to exploring new subspaces
    // group_sample_time                            //period on which to schedule thread restart
    // static int tgroup_ratio = 0;                 //
    static bool enable_workstealing = false;
    static bool enable_threadstop = false;
    static bool enable_lkh = false;
    static bool enable_progress_estimation = false;

    //derived attributes
    static int max_edge_weight = 0;                 //highest weight of any edge in the cost graph
    static float pre_density = 0;                   //number of edges in precedence graph (including derived edges) / the maximum possible
    static int local_pool_size = 0;                 //determined based on presidence density
/////////////////////////////////////////


///////////Shared Resources//////////////
    static vector<vector<edge>> cost_graph;         //n by n matrix with the cost from node i to node j
    static vector<vector<int>> dependency_graph;    //list for each node i of every node j that is dependent on it
    static vector<vector<edge>> in_degree;          //list for each node i of every node j that it is dependent on (in edge format)
    static vector<vector<edge>> hungarian_graph;    //graph in the format that the Hungarian algorithm requires
    static vector<vector<int>> outgoing_graph;

    //static sop_state default_state;                 //state of the solver before enumeration begins
    static vector<Hungarian> initial_hungarian_state;   //for each thread, the Hungarian solver as it began
    static Hash_Map history_table(TABLE_SIZE);
    //GLOBAL POOL
    //LOCAL POOLS

    static vector<int> best_solution;               //the lowest cost solution found so far in any thread
    static int best_cost = INT_MAX;                 //the cost of best_solution

    // std::chrono::time_point<std::chrono::system_clock> start_time_limit;
    // std::chrono::time_point<std::chrono::system_clock> time_point;
/////////////////////////////////////////


///////////Synchronization Variables/////
    //pthread_mutex_t Sol_lock = PTHREAD_MUTEX_INITIALIZER;   //lock for best_solution and its cost
/////////////////////////////////////////


///////////Thread Stopping Variables/////
    //
/////////////////////////////////////////


///////////Work Stealing Variables///////
    //
/////////////////////////////////////////


///////////Thread Restart Variables//////
    //MADE REDUNDANT DO NOT IMPLEMENT
/////////////////////////////////////////


///////////LKH Variables/////////////////
    // thread LKH_thread;
    // int* bestBB_tour = NULL;                        //an array of the best solution not found by LKH, but 1-indexed
    // bool BB_SolFound = false;                       //whether the current best solution was found by B&B, rather than LKH
    // bool BB_Complete = false;
    // bool local_searchinit = true;
    // bool initial_LKHRun = true;
/////////////////////////////////////////


///////////Diagnostic Variables//////////
    //static vector<unsigned_long_64> enumerated_nodes;             //total number of nodes processed by each thread
    //static vector<unsigned long long> estimated_trimmed_percent;  //estimated percentage of entire tree pruned or fully enumerated in each thread, stored as an integer out of ULLONG_MAX
    //static vector<int_64> num_resume;
    //static vector<int_64> num_stop;
    //static vector<double> lp_time;
    //static vector<double> steal_wait;
    //static vector<vector<double>> proc_time;
    //static vector<int> steal_cnt;
    //something to track history entry usage
/////////////////////////////////////////



void solver::assign_parameter(vector<string> setting) {
    t_limit = atoi(setting[0].c_str());
    std::cout << "Time limit = " << t_limit << std::endl;

    global_pool_size = atoi(setting[1].c_str());
    std::cout << "GPQ size = " << global_pool_size << std::endl;
    local_depth = atoi(setting[2].c_str());
    std::cout << "LPQ depth = " << local_depth << std::endl;

    inhis_mem_limit = atof(setting[3].c_str());
    std::cout << "History table mem limit = " << inhis_mem_limit << std::endl;
    inhis_depth = atof(setting[4].c_str());
    std::cout << "History table depth to always add = " << inhis_depth << std::endl;

    // exploitation_per = atof(setting[5].c_str())/float(100);
    // std::cout << "Restart exploitation/exploration ratio is " << exploitation_per << std::endl;
    // group_sample_time = atoi(setting[6].c_str());
    // std::cout << "Group sample time = " << group_sample_time << std::endl;
    // tgroup_ratio = atoi(setting[7].c_str());
    // std::cout << "Number of promising thread per exploitation group = " << tgroup_ratio << std::endl;

    // if (!atoi(setting[8].c_str())) enable_workstealing = false;
    // else enable_workstealing = true;

    // if (!atoi(setting[9].c_str())) enable_threadstop = false;
    // else enable_threadstop = true;

    // if (!atoi(setting[10].c_str())) enable_lkh = false;
    // else enable_lkh = true;

    // if (!atoi(setting[11].c_str())) enable_progress_estimation = false;
    // else enable_progress_estimation = true;

    return;
}

void solver::retrieve_input() {
    ifstream inFile;
    string line;

    // string fd_name;
    // for (unsigned i = 0; i < filename.size(); i++) {
    //     if (i == filename.size() - 4) break;
    //     fd_name += filename[i];
    // }
    // fd_name += "_BB.sop";

    //inFile.open(fd_name);
    inFile.open(filename);
    if (inFile.fail()) {
        cerr << "Error: input file " << filename << " -> " << strerror(errno) << endl;
        exit(EXIT_FAILURE);
    }
    else cout << "Input file is " << filename << endl;

    //Read input one line at a time, storing the matrix
    vector<vector<int>> file_matrix;
    vector<int> file_line;
    bool edge_weight_section = false;
    while (getline(inFile,line)) {
        if (line == "EDGE_WEIGHT_SECTION") { //skip unnecessary headers
            edge_weight_section = true;

            //read the dimension line 
            getline(inFile,line);
            instance_size = stoi(line);
        }
        else if (edge_weight_section) { //read cost/precedence matrix
            stringstream sstream;
            sstream << line;
            string weight;
            int weight_num;
            while (sstream >> weight) {
                stringstream(weight) >> weight_num;
                file_line.push_back(weight_num);
            }
            file_matrix.push_back(file_line);
            file_line.clear();
        }
    }

    if (instance_size != int(file_matrix.size())) {
        cerr << "Error: Instance Size Ambiguous. Expected " << instance_size << ". Was " << file_matrix.size() << endl;
        exit(EXIT_FAILURE);
    }

    cost_graph = vector<vector<edge>>(instance_size);
    hungarian_graph = vector<vector<edge>>(instance_size);
    dependency_graph = vector<vector<int>>(instance_size);

    for (int i = 0; i < (int)file_matrix.size(); i++) {
        int j = 0;
        for (auto edge_weight: file_matrix[i]) {
            if (edge_weight < 0) {
                cost_graph[i].push_back(edge(i,j,file_matrix[j][i]));
                hungarian_graph[i].push_back(edge(i,j,-1));
                dependency_graph[j].push_back(i);
            }
            else { 
                cost_graph[i].push_back(edge(i,j,edge_weight));
                hungarian_graph[i].push_back(edge(i,j,edge_weight));
            }
            j++;
        }
    }

    return;
}

void solver::transitive_redundancy() {
    in_degree = std::vector<vector<edge>>(instance_size);
    for (int i = 0; i < instance_size; i++) {
        for (long unsigned int k = 0; k < dependency_graph[i].size(); k++) {
            int c = dependency_graph[i][k];
            in_degree[c].push_back(edge(i,c,hungarian_graph[i][c].weight));
        }
    }

    for(int i = 0; i < instance_size; ++i) {
        vector<edge> preceding_nodes;
        for (int k = 0; k < (int)dependency_graph[i].size(); k++) preceding_nodes.push_back(edge(i,dependency_graph[i][k],-1));
        unordered_set<int> expanded_nodes;
        for(int j = 0; j < (int)preceding_nodes.size(); ++j){
            vector<edge> st;
            st.push_back(preceding_nodes[j]);
            while(!st.empty()){
                edge dependence_edge = st.back();
                st.pop_back();
                if(dependence_edge.src != i){
                    hungarian_graph[dependence_edge.dest][i].weight = -1;
                    hungarian_graph[i][dependence_edge.dest].weight = -1;
                    expanded_nodes.insert(dependence_edge.dest);
                }

                for(int dest : dependency_graph[dependence_edge.dest]){
                    if(expanded_nodes.find(dest) == expanded_nodes.end()){
                        st.push_back(edge(dependence_edge.dest,dest,-1));
                        expanded_nodes.insert(dest);
                    }
                }
            }
        } 
    }

    for(int i = 0; i < instance_size; ++i) {
        const vector<edge> preceding_nodes = in_degree[i];
        unordered_set<int> expanded_nodes;
        for(int j = 0; j < (int)preceding_nodes.size(); ++j){
            vector<edge> st;
            st.push_back(preceding_nodes[j]);
            while(!st.empty()){
                edge dependence_edge = st.back();
                st.pop_back();
                if(dependence_edge.src != i){
                    hungarian_graph[i][dependence_edge.dest].weight = -1;
                    hungarian_graph[dependence_edge.dest][i].weight = -1;
                    expanded_nodes.insert(dependence_edge.dest);
                }
                for(const edge& e : in_degree[dependence_edge.dest]){
                    if(expanded_nodes.find(e.dest) == expanded_nodes.end()){
                        st.push_back(e);
                        expanded_nodes.insert(e.dest);
                    }
                }
            }
        } 
    }

    return;
}

size_t solver::transitive_closure(vector<vector<int>>& isucc_graph) {
    size_t node_num = 0;
    deque<int> dfs_queue;
    boost::dynamic_bitset<> taken_arr = boost::dynamic_bitset<>(instance_size);
    boost::dynamic_bitset<> finished_arr = boost::dynamic_bitset<>(instance_size);
    vector<boost::dynamic_bitset<>> transitive_closure_map = vector<boost::dynamic_bitset<>>(instance_size);
    for (int i = 0; i < instance_size; i++) transitive_closure_map[i] = boost::dynamic_bitset<>(instance_size);
    int selected_node = 0;

    for (int i = 0; i < instance_size; i++) {
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

    for (int i = 0; i < instance_size; i++) {
        for (int k = 0; k < instance_size; k++) {
            if (transitive_closure_map[i][k]) node_num++;
        }
    }
    //print_transitive_map();
    return node_num;
}

// void solver::local_pool_config(float precedence_density) {
//     //
// }

vector<int> solver::nearest_neighbor(vector<int>* partial_solution) {
    vector<int> solution;
    int current_node;
    int solution_cost = 0;
    bool visit_arr[instance_size];
    int depCnt_arr[instance_size];
    vector<vector<edge>> sorted_costgraph = cost_graph; 
    sort_weight(sorted_costgraph);
    memset(depCnt_arr,0,instance_size*sizeof(int));

    for (int i = 0; i < instance_size; i++) {
        visit_arr[i] = false;
        for (long unsigned int k = 0; k < dependency_graph[i].size(); k++) {
            depCnt_arr[dependency_graph[i][k]]++;
        }
    }

    current_node = 0;
    for (auto node : *partial_solution) {
        visit_arr[node] = true;
        for (long unsigned int i = 0; i < dependency_graph[node].size(); i++) {
            depCnt_arr[dependency_graph[node][i]]--;
        }
        solution.push_back(node);
        solution_cost += cost_graph[current_node][node].weight;
        current_node = node;
    }

    int num = solution.size();
    
    while (num < instance_size) {
        bool taken = false;
        for (auto node: sorted_costgraph[current_node]) {
            if (!visit_arr[node.dest] && !depCnt_arr[node.dest]) {
                current_node = node.dest;
                solution_cost += node.weight;
                solution.push_back(current_node);
                num++;
                visit_arr[node.dest] = true;
                taken = true;
                break;
            }
        }
        if (!taken) {
            std::cout << "Error generating Nearest Neighbor Heuristic" << std::endl;
            std::cout << "current node is " << current_node << std::endl;
            for (auto node: sorted_costgraph[current_node]) {
                cout << node.dest << "," << visit_arr[node.dest] << "," << depCnt_arr[node.dest] << endl;
            }
            exit(EXIT_FAILURE);
        }
        for (long unsigned int i = 0; i < dependency_graph[current_node].size(); i++) {
            depCnt_arr[dependency_graph[current_node][i]]--;
        }
    }

    best_cost = solution_cost;
    return solution;
}

void solver::sort_weight(vector<vector<edge>>& graph) {
    int size = graph.size();
    for (int i = 0; i < size; i++) {
        stable_sort(graph[i].begin(),graph[i].end(),compare_edge);
    }
    return;
}

int solver::get_maxedgeweight() {
    int max = 0;
    for (int i = 0; i < instance_size; i++) {
        for (edge edge : cost_graph[i]) {
            int weight = edge.weight;
            if (weight > max) max = weight;
        }
    }
    return max;
}

vector<vector<int>> solver::get_cost_matrix(int max_edge_weight) {
    vector<vector<int>> matrix(instance_size);
    for(int i = 0; i < instance_size; ++i){
		matrix[i] = vector<int>(instance_size, max_edge_weight*2);
	}

    for (vector<edge> edge_list : hungarian_graph) {
        for (auto edge : edge_list) {
            int i = edge.src;
            int k = edge.dest;
            int weight = edge.weight;
            if (weight != -1 && i != k) {
                matrix[i][k] = weight * 2;
            }
        }
    }

    return matrix;
}

void solver::solve(string f_name, int thread_num) {
    if (thread_num < 1) {
        std::cerr << "Invalid Thread Number Input" << std::endl;
        exit(EXIT_FAILURE);
    }

    if (enable_lkh) thread_total = thread_num - 1;
    else thread_total = thread_num;
    filename = f_name;

    retrieve_input();
    transitive_redundancy();

    //Calculate Precedance Constraint Density
    int total_edges = 0;
    for (unsigned i = 1; i < dependency_graph.size() - 1; i++) {
        total_edges += dependency_graph[i].size() - 1;
    }
    pre_density = float(transitive_closure(dependency_graph)) / float((instance_size)*(instance_size-1));
    std::cout << "Precedance Density = " << pre_density << std::endl;

    //local_pool_config(float((instance_size - 1)*(instance_size-2))/float(total_edges));

    //Find Initial Best Solution
    std::cout << "Nearest Neighbor Heuristic Initialized" << std::endl;
    vector<int> temp_solution;
    temp_solution.push_back(0);
    best_solution = nearest_neighbor(&temp_solution);
    std::cout << "Initial Best Solution Is " << best_cost << std::endl;

    //I DON'T KNOW WHAT THIS IS FOR
    // bestBB_tour = new int[instance_size];
    // for (int i = 0; i < instance_size; i++) bestBB_tour[i] = best_solution[i] + 1;

    //THREAD RESTART
    // promise_Tlimit = (thread_total * exploitation_per) / tgroup_ratio;
    // std::cout << "Maximum exploitation group during thread restart is set to " << promise_Tlimit << std::endl;

    max_edge_weight = get_maxedgeweight();
    problem_state.hungarian_solver = Hungarian(instance_size, max_edge_weight+1, get_cost_matrix(max_edge_weight+1));
    problem_state.hungarian_solver.start();
    problem_state.depCnt = vector<int>(instance_size,0);
    problem_state.taken_arr = vector<int>(instance_size,0);

    for (int i = 0; i < instance_size; i++) {
        for (unsigned k = 0; k < dependency_graph[i].size(); k++) {
            problem_state.depCnt[dependency_graph[i][k]]++;
        }
    }
    // default_state = problem_state; //a copy of problem_state, since structs are passed by value
    std::cout << "Instance size is " << instance_size - 2 << std::endl;

    // thread_load = new load_stats [thread_total];

    history_table.set_up_mem(thread_total,instance_size);

    // abandon_wlk_array = vector<bool_64>(thread_total);
    // num_resume = vector<int_64>(thread_total);
    // num_stop = vector<int_64>(thread_total);
    // glp = vector<lptr_64>(thread_total);
    // lp_lock = vector<mutex_64>(thread_total);
    // busy_arr = vector<atomwrapper<bool>>(thread_total);

    // steal_wait = vector<double>(thread_total,0);
    // lp_time = vector<double>(thread_total,0);
    // proc_time = vector<vector<double>>(thread_total);
    // for (int i = 0; i < thread_total; i++) proc_time[i] = vector<double>(3,0);
    // initial_hungstate = vector<Hungarian>(thread_total);
    // for (int i = 0; i < thread_total; i++) initial_hungstate.push_back(default_state.hungarian_solver);
    // steal_cnt = vector<int>(thread_total,0);
    // enumerated_nodes = vector<unsigned_long_64>(thread_total);
    // estimated_trimmed_percent = vector<unsigned long long>(thread_total,0);

    //if (enable_lkh) LKH_thread = thread(run_lkh);

    auto start_time = chrono::high_resolution_clock::now();
    solve_parallel(thread_total,global_pool_size);
    auto end_time = chrono::high_resolution_clock::now();

    //if (enable_lkh) if (LKH_thread.joinable()) LKH_thread.join();

    auto total_time = chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
    cout << "------------------------" << thread_total << " thread" << "------------------------------" << endl;
    cout << best_cost << "," << setprecision(4) << total_time / (float)(1000000) << endl;
    
    
    ///// BEGIN DIAGNOSTICS /////

    // for (int i = 0; i < thread_total; i++) {
    //     cout << "--------------------------------------------------------" << endl;
    //     cout << "Thread " << i << ": Steal Count = " << steal_cnt[i] << endl;
    //     //cout << "                    Wait Time = " << steal_wait[i] / (float)(1000000) << " s" << endl;
    //     //cout << "                    Local Pool Time = " << lp_time[i] / (float)(1000000) << " s" << endl;
    //     //cout << "                    Load Process Time 0 = " << proc_time[i][0] / (float)(1000000) << " s" << endl;
    //     //cout << "                    Load Process Time 1 = " << proc_time[i][1] / (float)(1000000) << " s" << endl;
    //     //cout << "                    Load Process Time 2 = " << proc_time[i][2] / (float)(1000000) << " s" << endl;
    //     //cout << "--------------------------------------------------------" << endl;
    // }



    // unsigned long total_enodes = 0;
    // for (int i = 0; i < thread_total; i++) total_enodes += enumerated_nodes[i].val;
    // for (int i = 0; i < thread_total; i++) {
    //     cout << "thread " << i << " enumerated nodes = " << enumerated_nodes[i].val << ", " << double(enumerated_nodes[i].val)/double(total_enodes) << "%" << endl;
    // }
    //history_table.print_curmem();
    //cout << "Total enumerated nodes are " << total_enodes << endl;



    //print out progress estimates
    // if (enable_progress_estimation)
    // {
    //     unsigned long long total_progress = 0;
    //     for (int i = 0; i < thread_total; i++)
    //     {
    //         total_progress += estimated_trimmed_percent[i];
    //         //cout << "thread " << i << " trimmed = " << ((double) estimated_trimmed_percent[i])/ULLONG_MAX*100 << "% of the tree" << endl;
    //     }
    //     cout << "Total Progress = " << round(((double) total_progress)/ULLONG_MAX*100) << "%" << endl;
    //     cout << total_progress << " / " << ULLONG_MAX << endl;
    // }
        
    return;
}

// void solver::run_lkh() {
//     while(!BB_Complete) {
//         LKH(&filename[0],initial_LKHRun);
//         if (initial_LKHRun) {
//             initial_LKHRun = false;
//         }
//         BB_SolFound = false;
//     }
//     return;
// }

void solve_parallel(int thread_num, int pool_size) {
    
}





int solver::dynamic_hungarian(int src, int dest) {
    problem_state.hungarian_solver.fix_row(src, dest);
	problem_state.hungarian_solver.fix_column(dest, src);
	problem_state.hungarian_solver.solve_dynamic();
    //number_of_lb++;
    return problem_state.hungarian_solver.get_matching_cost()/2;
}

bool solver::history_utilization(pair<boost::dynamic_bitset<>,int>& key,int* lowerbound,bool* found,HistoryNode** entry, int cost) {
    size_t val = hash<boost::dynamic_bitset<>>{}(key.first);
    int bucket_location = (val + key.second) % TABLE_SIZE;
    HistoryNode* history_node = history_table.retrieve(key,bucket_location);

    if (history_node == NULL) return true;

    *found = true;
    *entry = history_node;
    HistoryContent content = history_node->Entry.load();
    *lowerbound = content.lower_bound;

    if (cost >= content.prefix_cost) return false;

    int target_ID = history_node->active_threadID;
    int imp = content.prefix_cost - cost;
    
    if (!history_node->explored) {
        if (enable_threadstop && active_thread > 0) {
            buffer_lock.lock();
            if (request_buffer.empty() || request_buffer.front().target_thread != target_ID || request_buffer.front().target_depth > (int)problem_state.cur_solution.size()) {
                //num_stop[thread_id].val++;
                request_buffer.push_front({problem_state.cur_solution.back(),(int)problem_state.cur_solution.size(),
                                           content.prefix_cost,target_ID,key.first});
                if (!stop_sig) {
                    stop_cnt = 0;
                    stop_sig = true;
                }
            }
            buffer_lock.unlock();
        }
    }

    if (imp <= content.lower_bound - best_cost) {
        return false;
    }
    history_node->Entry.store({cost,content.lower_bound - imp});
    *lowerbound = content.lower_bound - imp;
    *entry = history_node;
    history_node->explored = false;
    history_node->active_threadID = thread_id;
    
    return true;
}

void solver::push_to_historytable(pair<boost::dynamic_bitset<>,int>& key,int lower_bound,HistoryNode** entry,bool backtracked) {
    if (entry == NULL) history_table.insert(key,problem_state.cur_cost,lower_bound,thread_id,backtracked,problem_state.cur_solution.size());
    else *entry = history_table.insert(key,problem_state.cur_cost,lower_bound,thread_id,backtracked,problem_state.cur_solution.size());
    return;
}









//Compare two edges in the cost graph, a and b by their weight. a > b (therefore "better"), if a.weight < b.weight
bool compare_edge(const edge a, const edge b) { return a.weight < b.weight; }

bool split_sort(const sop_state& src, const sop_state& dest) {
    if (src.cur_solution.size() == dest.cur_solution.size()) return src.load_info < dest.load_info;
    return src.cur_solution.size() < dest.cur_solution.size();
}
bool GPQ_sort(const instrct_node& src, const instrct_node& dest) {
    return src.load_info > dest.load_info;
}
bool shared_pool_sort(const instrct_node& src, const instrct_node& dest) {
    if (src.sequence.size() == dest.sequence.size()) return src.load_info > dest.load_info;
    return src.sequence.size() < dest.sequence.size();
}
bool local_sort(const instrct_node& src, const instrct_node& dest) {
    return src.load_info > dest.load_info;
}
bool bound_sort(const node& src,const node& dest) {
    if (src.lb == dest.lb) return src.nc > dest.nc;
    return src.lb > dest.lb;
}
bool nearest_sort(const node& src,const node& dest) {
    if (src.nc == dest.nc) return src.lb > dest.lb;
    return src.nc > dest.nc;
}