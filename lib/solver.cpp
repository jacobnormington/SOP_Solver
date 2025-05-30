#include "solver.hpp"

extern "C"
{
#include "LKH/LKHmain.h"
}

#define TABLE_SIZE 541065431 // number of buckets in the history table

//////Runtime Parameters (Read Only)/////
// from command line arguments
static string filename;      // name of the sop input file
static int thread_total = 0; // number of threads to use for B&B enumeration (not counting the LKH thread)
static bool limit_insertion = false;

// from config file
static int t_limit = 0;          // time limit, in seconds
static int global_pool_size = 0; // Minimum size of the global pool at before enumeration begins
// static int local_depth = 0;                     //Minimum size to maintain in the local pool
static float inhis_mem_limit = -1; // 0-1, the percentage of memory usage beyond which new entries shouldn't be added to the history table
static float mem_limit = -1;       // 0-1, the percentage of memory usage beyond which we will start blocking new entries into the table
static int number_of_groups = 1;
static int bucket_size = 0;
static bool is_all_table_blocked = false;
// static unsigned int inhis_depth = -1;           //after inhis_mem_limit exceeded, will still add an entry if the current depth is less than inhis_depth
// int exploitation_per;                        //percent of threads that should be devoted to searching already promising subspaces in thread restart, while 1 - exploitation_per percent are devoted to exploring new subspaces
// group_sample_time                            //period on which to schedule thread restart
// static int tgroup_ratio = 0;                 //
static bool enable_workstealing = false;
static bool enable_threadstop = false;
static bool enable_lkh = true;
static bool enable_heuristic = false;
// static bool enable_progress_estimation = false;

// derived attributes
static int max_edge_weight = 0; // highest weight of any edge in the cost graph
static float pre_density = 0;   // number of edges in precedence graph (including derived edges) / the maximum possible
// static int local_pool_size = 0;                 //determined based on presidence density
/////////////////////////////////////////

///////////Shared Resources//////////////
static vector<vector<edge>> cost_graph;      // n by n matrix with the cost from node i to node j
static vector<vector<int>> dependency_graph; // list for each node i of every node j that is dependent on it
static vector<vector<edge>> in_degree;       // list for each node i of every node j that it is dependent on (in edge format) i.e. what nodes must precede it
static vector<vector<edge>> hungarian_graph; // graph in the format that the Hungarian algorithm requires
static vector<vector<int>> outgoing_graph;

static sop_state default_state; // state of a solver before the first node was taken
// static vector<Hungarian> initial_hungarian_state;   //for each thread, the Hungarian solver as it began
static History_Table history_table(TABLE_SIZE);           // the history table
static vector<path_node> global_pool;                     // a global pool of nodes that haven't yet been processed by any thread, only ever take from the back
static local_pool *local_pools;                           // each thread's local pool, with the internal tools to manage them, only ever take from the back
static vector<atomic<unsigned long long>> work_remaining; // used for work stealing, hold an estimate of how much work is left for a thread to do

///////////Synchronization Variables/////
// pthread_mutex_t Sol_lock = PTHREAD_MUTEX_INITIALIZER;   //lock for any updates to best_solution and its cost
static mutex best_solution_lock;
static mutex global_pool_lock; // lock for getting nodes from the global pool
// static mutex Split_lock;
// static mutex asssign_mutex;
// static mutex thread_load_mutex;
// static condition_variable Idel;
// static condition_variable Thread_Stop_Check;
// static condition_variable Resume_State;
// static vector<int> selected_orgin;
// static mutex Select_Mutex;
// static mutex Select_SharedMutex;
// static mutex Resume_Lock;
// static mutex launch_lck;
// pthread_mutex_t Sol_lock = PTHREAD_MUTEX_INITIALIZER;   //lock for any updates to best_solution and its cost

static mutex diagnostics_lock;
static int diagonstics_period = 10;
static float diagonstics_targetTime = 0;
static atomic<bool> time_out(false);  // whether the instance has timed out
static atomic<int> active_threads(0); // the number of threads still working
// static atomic<int> selected_thread (-1);
// static atomic<int> restart_cnt (0);
// static atomic<int> total_restarts (0);
// static atomic<unsigned> idle_counter (0);
// static atomic<size_t> resload_cnt (0);
// static atomic<bool> limit_insert (false);
// static atomic<bool> check_status_safe (true);
// static atomic<bool> resume_success (true);
// static atomic<bool> resume_check (false);
// static atomic<bool> exploit_init (false);

static vector<int> best_solution; // the lowest cost solution found so far in any thread
int best_cost = INT_MAX;          // the cost of best_solution, this is an extern (global) variable shared by LKH

static timer main_timer; // when solve_parallel started (before processing begins, but after all the basic setup, file parsing, etc.)

// std::chrono::time_point<std::chrono::system_clock> time_point;

/////////////////////////////////////////

///////////Thread Stopping Variables/////
static deque<request_packet> request_buffer; //
static mutex buffer_lock;                    // lock for accessing the request_buffer
// static vector<thread_request> thread_requests;
static vector<thread_request> thread_requests(32);

// static mutex pause_lock;                       //
// static mutex ptselct_lock;                     //
static atomic<bool> stop_sig(false); // if any threads are currently being requested to stop
// static vector<long long int> history_table_pruning_success = vector<long long int>(380); // uncomment to track the pruning at each depth and set the size of the instance

static atomic<int> thread_stop_requested(0);
static atomic<int> thread_stop_check(0);
static atomic<int> thread_stopped_successfully(0); // how many threads should stop
/////////////////////////////////////////

///////////Work Stealing Variables///////
//
/////////////////////////////////////////

///////////Thread Restart Variables//////
// MADE REDUNDANT DO NOT IMPLEMENT
/////////////////////////////////////////

///////////LKH Variables/////////////////
thread LKH_thread;
int *bestBB_tour = NULL;  // an array of the best solution not found by LKH, but 1-indexed
bool BB_SolFound = false; // whether the current best solution was found by B&B, rather than LKH
bool BB_Complete = false;
bool local_searchinit = true;
bool initial_LKHRun = true;
pthread_mutex_t Sol_lock = PTHREAD_MUTEX_INITIALIZER;
/////////////////////////////////////////

///////////Diagnostic Variables//////////
static vector<unsigned long long> enumerated_nodes; // total number of nodes processed by each thread
static atomic<int> times_work_stolen;
static atomic<int> steal_misses;
static vector<atomic<int>> steal_attempts = vector<atomic<int>>(32);
static vector<atomic<int>> steal_success = vector<atomic<int>>(32);
static atomic<double> time_workstealing;
static vector<double> steal_times;
static mutex steal_times_lock;
// static vector<unsigned long long> estimated_trimmed_percent;  //estimated percentage of entire tree pruned or fully enumerated in each thread, stored as an integer out of ULLONG_MAX
// TODO: change estimated_trimmed_percent to use unsigned_long_64 (and ULONG_MAX) instead of unsigned long long (and ULLONG_MAX)
// static vector<int_64> num_resume;
// static vector<int_64> num_stop;
// static vector<double> lp_time;
// static vector<double> steal_wait;
// static vector<vector<double>> proc_time;
// static vector<int> steal_cnt;
// something to track history entry usage
/////////////////////////////////////////

/* --------------------- Static Functions -------------------------*/
// Compare two edges in the cost graph, a and b by their weight. a > b (therefore "better"), if a.weight < b.weight
bool compare_edge(const edge a, const edge b) { return a.weight < b.weight; }
// TODO: ensure that this actually means that you are picking the "best" nodes, and it isn't accidentally reverse-sorted
// sort by decreasing lower bound (back is the best)
bool global_pool_sort(const path_node &src, const path_node &dest) { return src.lower_bound > dest.lower_bound; }
// sort by decreasing lower bound (back is the best)
bool local_pool_sort(const path_node &src, const path_node &dest) { return src.lower_bound > dest.lower_bound; }

static double gp_const;     // store the total work in the global pool initially
static double gp_remaining; // variable to store the number of remaining work from global pool

void lkh()
{
    while (!BB_Complete)
    {
        LKH(&filename[0], initial_LKHRun);
        if (initial_LKHRun)
        {
            initial_LKHRun = false;
        }
        BB_SolFound = false;
        // TODO: add entries to the history table corresponding to the LKH solution
    }
    return;
}

void print_diagnostics()
{
    if (main_timer.get_time_seconds() > diagonstics_targetTime)
    {
        diagnostics_lock.lock();
        if (main_timer.get_time_seconds() <= diagonstics_targetTime)
        {
            diagnostics_lock.unlock();
            return;
        }
        cout << main_timer.get_time_seconds() << endl;
        local_pools->print();
        for (int i = 0; i < (int)work_remaining.size(); i++)
        {
            cout << i << ": " << work_remaining[i] << ", ";
        }
        cout << endl;
        cout << active_threads << endl;

        diagonstics_targetTime = main_timer.get_time_seconds() + diagonstics_period;
        diagnostics_lock.unlock();
    }
}
/* ---------------------       END        -------------------------*/
/* --------------------- Static Functions -------------------------*/

void solver::assign_parameter(vector<string> setting)
{
    t_limit = atoi(setting[0].c_str());
    std::cout << "Time limit = " << t_limit << std::endl;

    global_pool_size = atoi(setting[1].c_str());
    std::cout << "GPQ size = " << global_pool_size << std::endl;

    inhis_mem_limit = atof(setting[3].c_str());
    std::cout << "History table mem limit = " << inhis_mem_limit << std::endl;

    // inhis_depth = atof(setting[4].c_str());
    // std::cout << "History table depth to always add = " << inhis_depth << std::endl;

    // exploitation_per = atof(setting[5].c_str())/float(100);
    // std::cout << "Restart exploitation/exploration ratio is " << exploitation_per << std::endl;
    // group_sample_time = atoi(setting[6].c_str());
    // std::cout << "Group sample time = " << group_sample_time << std::endl;
    // tgroup_ratio = atoi(setting[7].c_str());
    // std::cout << "Number of promising thread per exploitation group = " << tgroup_ratio << std::endl;

    if (!atoi(setting[8].c_str()))
        enable_workstealing = false;
    else
        enable_workstealing = true;

    if (!atoi(setting[9].c_str()))
        enable_threadstop = false;
    else
        enable_threadstop = true;

    if (!atoi(setting[10].c_str()))
        enable_lkh = false;
    else
        enable_lkh = true;

    // if (!atoi(setting[11].c_str())) enable_progress_estimation = false;
    // else enable_progress_estimation = true;
    number_of_groups = atoi(setting[12].c_str());
    std::cout << "Number of groups = " << number_of_groups << std::endl;

    bucket_size = atoi(setting[13].c_str());
    if (bucket_size > 0)
        std::cout << "Group size in config = " << bucket_size << std::endl;

    if (number_of_groups - 1 >= 1)
        mem_limit = inhis_mem_limit - (number_of_groups - 1) * 0.1;
    else
    {
        mem_limit = inhis_mem_limit; // if the bucket size is one, the mem limit would be same as config file
        is_all_table_blocked = true; // by setting this to true, we are saying that all the buckets are blocked except first bucket
    }
    std::cout << "Blocking mem limit = " << mem_limit << std::endl;

    return;
}

void print_workdone()
{
    cout << "gp const: " << gp_const << "\n";
    cout << "gp remaining: " << gp_remaining << "\n";

    long double total_work_done = 0; // work done in threads  // thread - work remaining
    long double total_work_rem = 0;

    // Calculate the total remaining work
    for (const auto &work : work_remaining)
    {
        total_work_done += (1 - (static_cast<long double>(work.load()) / ULLONG_MAX));
        total_work_rem += (static_cast<long double>(work.load()) / ULLONG_MAX);
    }

    std::cout << "Total Work done in threads: " << total_work_done << std::endl;
    std::cout << "Total Work remaining in threads: " << total_work_rem << std::endl;

    long double total_global_work_done = 0;
    // Total global work done

    total_global_work_done = gp_const - gp_remaining - 31 + total_work_done;

    std::cout << "Total global work done: " << total_global_work_done << std::endl;
    // Calculate the percentage of work done
    long double percentage_done = (total_global_work_done / gp_const) * 100;
    std::cout << "Percentage of work done: " << percentage_done << "%" << std::endl;
}
void solver::solve(string f_name, int thread_num)
{
    if (thread_num < 1)
    {
        std::cerr << "Invalid Thread Number Input" << std::endl;
        exit(EXIT_FAILURE);
    }

    if (enable_lkh)
        thread_total = thread_num - 1;
    else
        thread_total = thread_num;
    if (global_pool_size < thread_total)
        global_pool_size = thread_total;
    filename = f_name;
    retrieve_input();
    transitive_redundancy();

    // Calculate Precedance Constraint Density
    int total_edges = 0;
    for (unsigned i = 1; i < dependency_graph.size() - 1; i++)
    {
        total_edges += dependency_graph[i].size() - 1;
    }
    pre_density = float(transitive_closure(dependency_graph)) / float((instance_size) * (instance_size - 1));
    std::cout << "Precedance Density = " << pre_density << std::endl;

    // local_pool_config(float((instance_size - 1)*(instance_size-2))/float(total_edges));

    // Find Initial Best Solution
    std::cout << "Nearest Neighbor Heuristic Initialized" << std::endl;
    vector<int> temp_solution;
    temp_solution.push_back(0);
    best_solution = nearest_neighbor(&temp_solution);
    std::cout << "Initial Best Solution Is " << best_cost << std::endl;

    // necessary for LKH
    //  bestBB_tour = new int[instance_size];
    //  for (int i = 0; i < instance_size; i++) bestBB_tour[i] = best_solution[i] + 1;

    // THREAD RESTART
    //  promise_Tlimit = (thread_total * exploitation_per) / tgroup_ratio;
    //  std::cout << "Maximum exploitation group during thread restart is set to " << promise_Tlimit << std::endl;

    max_edge_weight = get_maxedgeweight();
    problem_state.hungarian_solver = Hungarian(instance_size, max_edge_weight + 1, get_cost_matrix(max_edge_weight + 1));
    problem_state.hungarian_solver.start();
    problem_state.depCnt = vector<int>(instance_size, 0);
    problem_state.taken_arr = boost::container::vector<bool>(instance_size, false);
    problem_state.current_cost = 0; // TODO: this shouldn't be necessary, since it is in the struct definition, but it seems to be, investigate

    for (int i = 0; i < instance_size; i++)
    {
        for (unsigned k = 0; k < dependency_graph[i].size(); k++)
        {

            problem_state.depCnt[dependency_graph[i][k]]++;
        }
    }

    // we need to split the instance into 3 equal parts if the group size is 0 in the configuration file
    if (bucket_size == 0)
    {
        bucket_size = instance_size / number_of_groups;
        std::cout << "Group size after splitting instance size in 3 equal part = " << bucket_size << std::endl;
    }
    else if (bucket_size * (number_of_groups - 1) >= instance_size - 2)
    {

        cout << "Exiting the code. Invalid bucket size" << endl;
        exit(EXIT_FAILURE);
    }

    default_state = problem_state; // a copy of problem_state, since structs are passed by value
    std::cout << "Instance size is " << instance_size - 2 << std::endl;

    // thread_load = new load_stats [thread_total];
    local_pools = new local_pool(thread_total);
    history_table.initialize(thread_total, TABLE_SIZE, number_of_groups, bucket_size);
    // thread_requests.resize(thread_total);
    // for (int i = 0; i < thread_total; ++i)
    // {
    //     thread_requests[i] = std::make_unique<thread_request>();
    // }

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
    // initial_hungarian_state = vector<Hungarian>(thread_total);
    // for (int i = 0; i < thread_total; i++) initial_hungarian_state.push_back(default_state.hungarian_solver);
    // steal_cnt = vector<int>(thread_total,0);
    // enumerated_nodes = vector<unsigned_long_64>(thread_total);
    // estimated_trimmed_percent = vector<unsigned long long>(thread_total,0);

    bestBB_tour = new int[instance_size];
    for (int i = 0; i < instance_size; i++)
        bestBB_tour[i] = best_solution[i] + 1;

    if (enable_lkh)
        LKH_thread = thread(lkh);

    auto start_time = chrono::high_resolution_clock::now();
    solve_parallel();
    auto end_time = chrono::high_resolution_clock::now();

    if (enable_lkh)
        if (LKH_thread.joinable())
            LKH_thread.join();

    // DIAGNOSTIC : Enumerated Nodes
    //  unsigned long long enumerated_nodes_sum = 0;
    //  for(int i = 0; i < enumerated_nodes.size(); i++){
    //      enumerated_nodes_sum += enumerated_nodes[i];
    //  }

    auto total_time = chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
    std::cout << "------------------------" << thread_total << " thread"
              << "------------------------------" << std::endl;
    cout << "thread stop requested: " << thread_stop_requested << "\n";
    cout << "thread stop check: " << thread_stop_check << "\n";
    cout << "thread stopped successfully: " << thread_stopped_successfully << "\n";

    for (int i = 0; i < steal_success.size(); i++)
        cout << steal_success[i] << ", ";
    cout << endl;
    for (int i = 0; i < steal_success.size(); i++)
        cout << steal_attempts[i] << ", ";
    cout << endl;
    cout << "total work stolen: " << times_work_stolen << endl;
    cout << "steal misses: " << steal_misses << endl;

    double percent_time_active = (((double)total_time / 1000000) * 32 - time_workstealing) / ((double)total_time / 1000000 * 32);
    cout << "active time: " << percent_time_active << endl;

    std::cout << "best_cost: " << best_cost << "," << setprecision(4) << total_time / (float)(1000000) << std::endl
              << std::endl;

    for (int i = 0; i < steal_times.size(); i++)
    {
        cout << steal_times[i] << endl;
    }

    print_workdone();
    // to count the number of entries at different level in history table and their references
    // history_table.track_entries_and_references();

    /** uncomment to print the pruning happening at each level
    // Print the data at each index
    int total = 0;
    for (size_t i = 0; i < history_table_pruning_success.size(); ++i)
    {
        if (history_table_pruning_success[i] != 0)
        {
            total += history_table_pruning_success[i];
            cout << "Index " << i << ": " << history_table_pruning_success[i] << endl;
        }
    }
    */

    // std::cout << "Best solution path: ";
    // for (int element : best_solution)
    // {
    //     std::cout << element << ", ";
    // }
    // std::cout << std::endl;

    ///// BEGIN DIAGNOSTICS /////

    // int total_cost = 0;
    //  std::cout << "Final solution: " << best_solution[0] << " ";
    //  for (int i = 1; i < (int) best_solution.size(); i++) {
    //      std::cout << best_solution[i] << " ";
    //      int src = best_solution[i-1];
    //      int dst = best_solution[i];
    //      total_cost += cost_graph[src][dst].weight;
    //  }
    // std::cout << std::endl;
    // std::cout << "Cost of final solution is: " << total_cost << ", Reported: " << best_cost << std::endl << std::endl;

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
    //     std::cout << "thread " << i << " enumerated nodes = " << enumerated_nodes[i].val << ", " << double(enumerated_nodes[i].val)/double(total_enodes) << "%" << std::endl;
    // }
    // std::cout << "Total enumerated nodes are " << total_enodes << std::endl << std::endl;

    // history_table.print_curmem();
    // std::cout << std::endl;

    // print out progress estimates
    //  if (enable_progress_estimation)
    //  {
    //      unsigned long long total_progress = 0;
    //      for (int i = 0; i < thread_total; i++)
    //      {
    //          total_progress += estimated_trimmed_percent[i];
    //          //std::cout << "thread " << i << " trimmed = " << ((double) estimated_trimmed_percent[i])/ULLONG_MAX*100 << "% of the tree" << std::endl;
    //      }
    //      std::cout << "Total Progress = " << round(((double) total_progress)/ULLONG_MAX*100) << "%" << std::endl;
    //      std::cout << total_progress << " / " << ULLONG_MAX << std::endl << std::endl;
    //  }

    return;
}

void solver::solve_parallel()
{
    main_timer.restart();

    vector<solver> solvers(thread_total);
    deque<sop_state> *solver_container;
    vector<thread> Thread_manager(thread_total);

    /* Begin splitting nodes of the enumeration tree until GPQ has reached a minimum size. */
    vector<node> ready_list;
    solver_container = new deque<sop_state>();

    // Take the first node (this is the "virtual" node that is always taken first)
    problem_state.current_path.push_back(0);
    problem_state.taken_arr[0] = true;
    for (int vertex : dependency_graph[0])
        problem_state.depCnt[vertex]--;
    for (int i = 0; i < instance_size; i++)
    {
        if (!problem_state.depCnt[i] && !problem_state.taken_arr[i])
        {
            // Push vertices with 0 dependencies into the ready list, unless they are already taken
            ready_list.push_back(node(i));
        }
    }
    // for (int i : problem_state.depCnt)
    // {
    //     std::cout << i << " ";
    // }
    // std::cout << " " << std::endl;

    // Process first generation of nodes (the ready list as defined above)
    sop_state initial_state = problem_state;
    for (auto node : ready_list)
    {

        int taken_node = node.n;
        int cur_node = initial_state.current_path.back();
        initial_state.current_path.push_back(taken_node);
        initial_state.taken_arr[taken_node] = true;
        for (int vertex : dependency_graph[taken_node])
            initial_state.depCnt[vertex]--;
        initial_state.current_cost += cost_graph[cur_node][taken_node].weight; // technically, this is always 0, but required for strict correctness

        initial_state.hungarian_solver.fix_row(cur_node, taken_node);
        initial_state.hungarian_solver.fix_column(taken_node, cur_node);
        initial_state.hungarian_solver.solve_dynamic();
        initial_state.lower_bound = initial_state.hungarian_solver.get_matching_cost() / 2;
        initial_state.origin_node = taken_node;

        /* I DON'T KNOW WHAT THIS IS ACCOUNTING FOR */
        if (pre_density != 0)
            solver_container->push_back(initial_state); // a copy, since sop_state is a struct, passed by value
        else if (instance_size > thread_total)
        {
            global_pool.push_back(path_node(initial_state.current_path, initial_state.lower_bound, initial_state.origin_node));
        }

        initial_state.taken_arr[taken_node] = false;
        for (int vertex : dependency_graph[taken_node])
            initial_state.depCnt[vertex]++;
        initial_state.current_cost -= cost_graph[cur_node][taken_node].weight;
        initial_state.current_path.pop_back();
        initial_state.hungarian_solver.undue_row(cur_node, taken_node);
        initial_state.hungarian_solver.undue_column(taken_node, cur_node);
    }

    // Repeat split operation until enough nodes are produced, but every node must be of the same depth
    while (global_pool.empty() && (solver_container->size() < (size_t)global_pool_size || split_level_check(solver_container)))
    {
        // sort(solver_container->begin(),solver_container->end(),split_sort);
        sop_state target = solver_container->front();
        if (target.current_path.size() == (unsigned)(instance_size - 1))
            break;
        solver_container->pop_front();

        ready_list.clear();
        for (int i = instance_size - 1; i >= 0; i--)
        {
            if (!target.depCnt[i] && !target.taken_arr[i])
                ready_list.push_back(node(i));
        }
        if (ready_list.empty()) // the instance is unsolveable if the ready_list is ever empty
        {
            std::cout << "Precedence Constraints are Unsatisfiable" << std::endl;
            std::cout << "Current path: ";
            for (auto node : target.current_path)
                std::cout << node << " ";
            std::cout << std::endl;
            std::cout << "Ready List is Empty." << std::endl;
            exit(EXIT_SUCCESS);
        }

        for (auto node : ready_list)
        {
            // Push split node back into Global Pool
            int taken_node = node.n;
            int cur_node = target.current_path.back();
            target.current_path.push_back(taken_node);
            target.taken_arr[taken_node] = true;
            for (int vertex : dependency_graph[taken_node])
                target.depCnt[vertex]--;
            target.current_cost += cost_graph[cur_node][taken_node].weight;
            if (target.current_path.size() == (size_t)instance_size && target.current_cost < best_cost)
            { // this shouldn't be possible
                best_solution = target.current_path;
                best_cost = target.current_cost;
            }

            if (target.current_cost < best_cost)
            {
                target.hungarian_solver.fix_row(cur_node, taken_node);
                target.hungarian_solver.fix_column(taken_node, cur_node);
                target.hungarian_solver.solve_dynamic();
                target.lower_bound = target.hungarian_solver.get_matching_cost() / 2;
                solver_container->push_back(target);
                target.hungarian_solver.undue_row(cur_node, taken_node);
                target.hungarian_solver.undue_column(taken_node, cur_node);
            }

            target.current_cost -= cost_graph[cur_node][taken_node].weight;
            for (int vertex : dependency_graph[taken_node])
                target.depCnt[vertex]++;
            target.taken_arr[taken_node] = false;
            target.current_path.pop_back();
        }
    }

    if (global_pool.empty())
    {
        for (auto problem : *solver_container)
        {
            global_pool.push_back(path_node(problem.current_path, problem.lower_bound, problem.origin_node));
        }
    }
    sort(global_pool.begin(), global_pool.end(), global_pool_sort);

    history_table.update_gp_depth(global_pool.at(0).sequence.size());
    gp_const = global_pool.size();
    gp_remaining = global_pool.size();
    // cout << "Setting the parameters to track the work done" << endl;

    delete solver_container;
    /* End Splitting Operation */

    /* Assign subproblems from the GPQ into solvers for each thread. */
    int thread_cnt = 0;
    boost::container::vector<bool> origin_taken_arr = boost::container::vector<bool>(instance_size, false);
    while (thread_cnt < thread_total)
    {                                                                  // continue, even taking duplicates from the same origin, in order to get work for every thread
        fill(origin_taken_arr.begin(), origin_taken_arr.end(), false); // to more evenly distribute work, only one child of each of the first generation should be taken
        for (int i = global_pool.size() - 1; i >= 0; i--)
        { // only ever take from the back of the global_pool, which is the best lower_bound
            if (thread_cnt >= thread_total)
                break;
            unsigned origin = global_pool[i].origin_node;

            if (!origin_taken_arr[origin])
            {
                path_node problem = global_pool[i];
                global_pool.erase(global_pool.begin() + i); // remove from the pool

                solvers[thread_cnt].problem_state = default_state;
                int taken_node = -1;
                int cur_node = problem.sequence.front(); // curr_node is 0
                int size = problem.sequence.size();
                solvers[thread_cnt].problem_state.current_path.push_back(0);
                solvers[thread_cnt].problem_state.taken_arr[cur_node] = true;
                for (int vertex : dependency_graph[cur_node])
                    solvers[thread_cnt].problem_state.depCnt[vertex]--;

                for (int k = 1; k < size; k++)
                { // begin after the first, trivial, node
                    taken_node = problem.sequence[k];
                    solvers[thread_cnt].problem_state.current_path.push_back(taken_node);
                    solvers[thread_cnt].problem_state.taken_arr[taken_node] = true;
                    for (int vertex : dependency_graph[taken_node])
                        solvers[thread_cnt].problem_state.depCnt[vertex]--;
                    solvers[thread_cnt].problem_state.current_cost += cost_graph[cur_node][taken_node].weight;
                    solvers[thread_cnt].problem_state.hungarian_solver.fix_row(cur_node, taken_node);
                    solvers[thread_cnt].problem_state.hungarian_solver.fix_column(taken_node, cur_node);
                    cur_node = taken_node; // TODO: Not sure if this is correct
                }

                solvers[thread_cnt].problem_state.origin_node = problem.origin_node;
                solvers[thread_cnt].problem_state.lower_bound = problem.lower_bound;
                // solvers[thread_cnt].problem_state.initial_depth = solvers[thread_cnt].problem_state.cur_solution.size();
                // solvers[thread_cnt].problem_state.current_node_value = problem.current_node_value; //for progress estimation

                solvers[thread_cnt].thread_id = thread_cnt;
                solvers[thread_cnt].instance_size = instance_size;
                // solvers[thread_cnt].lb_curlv = problem.lower_bound;
                //  solvers[thread_cnt].cur_active_tree = Active_Path(solvers[thread_cnt].problem_state.cur_solution.size());
                //  solvers[thread_cnt].cur_active_tree.set_threadID(thread_cnt, thread_total);

                // solvers[thread_cnt].local_pool = new deque<path_node>();
                // glp[thread_cnt].local_pool = solvers[thread_cnt].local_pool;

                // thread_load[thread_cnt].cur_sol = &solvers[thread_cnt].problem_state.cur_solution;

                boost::dynamic_bitset<> bit_vector(instance_size, false);
                for (auto node : solvers[thread_cnt].problem_state.current_path)
                {
                    bit_vector[node] = true;
                }
                int last_element = solvers[thread_cnt].problem_state.current_path.back();
                solvers[thread_cnt].problem_state.history_key = make_pair(bit_vector, last_element);

                thread_cnt++;
                origin_taken_arr[origin] = true;
            }
        }
    }

    work_remaining = std::vector<std::atomic<unsigned long long>>(thread_cnt); // PROGRESS initializing work remaining vector
    enumerated_nodes = std::vector<unsigned long long>(thread_cnt);
    for (int i = 0; i < thread_cnt; i++)
    {
        work_remaining[i] = ULLONG_MAX;
        enumerated_nodes[i] = 0;
    }

    for (int i = 0; i < thread_total; i++)
    {
        Thread_manager[i] = thread(&solver::enumerate, move(solvers[i]));
        active_threads++;
    }
    for (int i = 0; i < thread_total; i++)
    { // waits until every thread is finished
        if (Thread_manager[i].joinable())
        {
            Thread_manager[i].join();
        }
    }
    active_threads = 0;
    BB_Complete = true;

    if (time_out)
    {
        cout << "instance timed out " << endl;
    }

    return;
}

void solver::enumerate()
{
    while (!time_out)
    {
        // if (main_timer.get_time_seconds() > 3598)
        //     local_pools->print_top_sequence_sizes_end(thread_total);

        // PROGRESS variables
        int ready_node_count = 0;
        int pruned_count = 0;

        deque<path_node> ready_list;
        // bool limit_insertion = false;
        for (int taken_node = 0; taken_node < instance_size; taken_node++)
        {
            if (!problem_state.depCnt[taken_node] && !problem_state.taken_arr[taken_node])
            { // only consider nodes that haven't already been taken, and who have no remaining dependencies
                ready_node_count++;

                // triming
                int source_node = problem_state.current_path.back();
                int lower_bound = -1;
                bool taken = false;

                problem_state.current_path.push_back(taken_node);
                problem_state.history_key.first[taken_node] = true;
                problem_state.history_key.second = taken_node;
                problem_state.current_cost += cost_graph[source_node][taken_node].weight;

                HistoryNode *his_node = NULL;
                // Active_Node* active_node = NULL;

                if (problem_state.current_cost >= best_cost)
                { // backtracking
                    pruned_count++;
                    prune(source_node, taken_node);
                    continue;
                }

                if (problem_state.current_path.size() == (size_t)instance_size)
                { // if you've reached a leaf node (complete solution)
                    if (problem_state.current_cost < best_cost)
                    {
                        best_solution_lock.lock();
                        if (problem_state.current_cost < best_cost)
                        { // make sure it is still true
                            pthread_mutex_lock(&Sol_lock);
                            if (problem_state.current_cost < best_cost)
                            { // second check is for parallel reasons
                                best_cost = problem_state.current_cost;
                                best_solution = problem_state.current_path;
                                // LKH
                                if (enable_lkh)
                                {
                                    for (int i = 0; i < (int)(best_solution.size()); i++)
                                        bestBB_tour[i] = best_solution[i] + 1;
                                    BB_SolFound = true;
                                }

                                // DIAGNOSTIC: best cost
                                std::cout << "Best Cost = " << best_cost << " Found in Thread " << thread_id; // TODO: add toggle
                                std::cout << " at time = " << main_timer.get_time_seconds() << std::endl;
                            }
                            pthread_mutex_unlock(&Sol_lock);
                        }
                        best_solution_lock.unlock();
                    }

                    pruned_count++;
                    prune(source_node, taken_node);
                    continue;
                }

                // true: no entry exist in history table
                // false: someone else is performing better than this node
                // false: the improvement is not worth it (we are adding an entry in the request buffer)
                // true: the improvement is worth it (we are adding an entry in the request buffer)
                bool decision = history_utilization(problem_state.history_key, problem_state.current_cost, &lower_bound, &taken, &his_node);
                if (!taken)
                { // if there is no similar entry in the history table
                    lower_bound = dynamic_hungarian(source_node, taken_node);
                    // TODO_VIKAS: can we check the lower bound with the best cost before inserting into the history table

                    if (!limit_insertion)
                    {
                        if (history_table.get_current_size() < mem_limit * history_table.get_max_size())
                            push_to_history_table(problem_state.history_key, lower_bound, &his_node, false);
                        else if (number_of_groups == 1)
                        {
                            /**
                             * if number of groups is 1, we want to ignore the blocking and the deleting of the history table
                             * by calling free_subtable_memory function, we are blocking the insertion of the history table
                             */

                            history_table.free_subtable_memory(&mem_limit);
                            cout << "Blocking Insertion at time is: " << main_timer.get_time_seconds() << endl;

                            /** to prevent further checking, we set limit_insertion to true */
                            limit_insertion = true;
                        }
                        else
                        {
                            if (!is_all_table_blocked)
                            {
                                if (!history_table.check_and_manage_memory(problem_state.current_path.size(), &mem_limit, &is_all_table_blocked))
                                    push_to_history_table(problem_state.history_key, lower_bound, &his_node, false);
                            }
                            else
                            {
                                if (history_table.get_current_size() >= mem_limit * history_table.get_max_size())
                                {
                                    // local_pools->print_top_sequence_sizes(thread_total);
                                    bool is_space_increased_or_available = history_table.free_subtable_memory(&mem_limit);
                                    if (is_space_increased_or_available)
                                    {
                                        if (problem_state.current_path.size() <= bucket_size)
                                            push_to_history_table(problem_state.history_key, lower_bound, &his_node, false);
                                    }
                                    else
                                    {
                                        cout << "Blocking Insertion at time is: " << main_timer.get_time_seconds() << endl;

                                        /** to prevent further checking, we set limit_insertion to true */
                                        limit_insertion = true;
                                    }
                                }
                                else
                                    push_to_history_table(problem_state.history_key, lower_bound, &his_node, false);
                            }
                        }
                    }
                }
                else if (taken && !decision)
                { // if this path is dominated by another path
                    // tracking the pruning at current depth
                    // history_table_pruning_success[problem_state.current_path.size()]++;
                    pruned_count++;
                    prune(source_node, taken_node);
                    continue;
                }

                if (lower_bound >= best_cost)
                {
                    if (his_node != NULL)
                    {
                        HistoryContent content = his_node->entry.load();
                        // TODO_VIKAS: it will never be greater because prefix_cost is same as the problem_state current_cost
                        if (content.prefix_cost >= problem_state.current_cost)
                            his_node->explored = true;
                    }
                    pruned_count++;
                    prune(source_node, taken_node);
                    continue;
                }
                //"good" node, add it to the ready_list, then reset problem state
                path_node temp(problem_state.current_path, lower_bound, problem_state.origin_node, problem_state.history_key);
                ready_list.push_back(temp);
                problem_state.current_path.pop_back();
                problem_state.current_cost -= cost_graph[source_node][taken_node].weight;
                problem_state.history_key.first[taken_node] = false;
                problem_state.history_key.second = source_node;
            }
        }

        // PROGRESS
        // COMMENT:ready_node_count : number of nodes in the ready list
        unsigned long long next_work_above = problem_state.work_above / ready_node_count;
        unsigned long long remainder = problem_state.work_above % ready_node_count;
        work_remaining[thread_id] -= next_work_above * pruned_count + remainder;

        // COMMENT: Assigning work remaining for each path_node object in ready list
        for (unsigned i = 0; i < ready_list.size(); i++)
        {
            ready_list[i].current_node_value = next_work_above;
        }

        // DIAGNOSTIC: enum_nodes
        // enumerated_nodes[thread_id] += ready_node_count;

        // Sort the ready list and push into local pool
        if (!ready_list.empty())
            std::sort(ready_list.begin(), ready_list.end(), local_pool_sort);
        local_pools->push_list(thread_id, ready_list);

        int lb_liminsert = problem_state.lower_bound; // save lower bound through enumeration for limit insertion in the history table

        /* Begin enumeration. */
        path_node active_node;
        while (local_pools->pop_from_active_list(thread_id, active_node))
        {
            if (enable_threadstop)
            {
                // COMMENT: we are comparing the active_node(thread_id, last_node) with the request_buffer(thread_id, request_buffer)
                // Don't we need other parameters such as depth, prefix_cost or  history_key
                bool prefix_key_matched = false;
                if (check_stop_request(active_node.history_key, active_node.sequence, &prefix_key_matched))
                {
                    work_remaining[thread_id] -= active_node.current_node_value;
                    if (prefix_key_matched)
                        break;
                    else
                        continue;
                }
            }
            if (enumeration_pre_check(active_node))
                continue; // enumeration-time backtracking, and other preprocessing

            /* Take */
            int src = problem_state.current_path.back();  // the number of the predecessor to the node being considered
            int taken_node = active_node.sequence.back(); // this node's number
            problem_state.current_path.push_back(taken_node);
            problem_state.taken_arr[taken_node] = true;
            problem_state.current_cost += cost_graph[src][taken_node].weight;
            for (int vertex : dependency_graph[taken_node])
                problem_state.depCnt[vertex]--;
            problem_state.history_key.first[taken_node] = true;
            problem_state.history_key.second = taken_node;
            problem_state.hungarian_solver.fix_row(src, taken_node);
            problem_state.hungarian_solver.fix_column(taken_node, src);
            // HistoryNode *previous_hisnode = current_hisnode;
            // current_hisnode = history_entry;
            // problem_state.suffix_cost = 0;
            problem_state.enumeration_depth++;
            problem_state.work_above = active_node.current_node_value;

            enumerate();

            /* Untake */
            problem_state.enumeration_depth--;
            // problem_state.suffix_cost????
            // current_hisnode = previous_hisnode;
            problem_state.hungarian_solver.undue_row(src, taken_node);
            problem_state.hungarian_solver.undue_column(taken_node, src);
            problem_state.history_key.first[taken_node] = false;
            problem_state.history_key.second = src;
            for (int vertex : dependency_graph[taken_node])
                problem_state.depCnt[vertex]++;
            problem_state.current_cost -= cost_graph[src][taken_node].weight;
            problem_state.taken_arr[taken_node] = false;
            problem_state.current_path.pop_back();

            if (thread_id == 0)
            { // check if out of time
                if (main_timer.get_time_seconds() > t_limit)
                {
                    time_out = true;
                    active_threads = 0;
                    local_pools->pop_active_list(thread_id);
                    return;
                }
            }
        }
        while (local_pools->pop_from_active_list(thread_id, active_node))
        {
            work_remaining[thread_id] -= active_node.current_node_value;
        }
        local_pools->pop_active_list(thread_id); // TODO: make sure with thread stopping that this is handled properly
        // if (stop_init && (int)problem_state.cur_solution.size() <= stop_depth) {
        //     stop_init = false;
        //     stop_depth = -1;
        //     last_node = -1;
        // }

        // TODO_VIKAS: Shouldn't we pass true in the last parameter, since we the iteration is completed
        if (limit_insertion && history_table.get_current_size() < history_table.get_max_size())
        {
            push_to_history_table(problem_state.history_key, lb_liminsert, NULL, true);
        }

        // TODO: replace above with this, for taking depth into account
        //  if (limit_insertion && history_table.get_current_size() < history_table.get_max_size() && problem_state.current_path.size() >= inhis_depth) { //don't add if it was already added
        //      push_to_history_table(problem_state.history_key,lb_liminsert,NULL,false);
        //  }

        // cur_active_tree.pop_back(stop_init,Allocator);

        if (problem_state.enumeration_depth != 0 || !workload_request()) // only make a workload request from lowest recursion depth, then end if there is no more work
            return;
    }
    return;
}

void solver::retrieve_input()
{
    ifstream inFile;
    string line;

    // string fd_name;
    // for (unsigned i = 0; i < filename.size(); i++) {
    //     if (i == filename.size() - 4) break;
    //     fd_name += filename[i];
    // }
    // fd_name += "_BB.sop";

    // inFile.open(fd_name);
    inFile.open(filename);
    if (inFile.fail())
    {
        cerr << "Error: cannot open input file " << filename << " -> " << strerror(errno) << endl;
        exit(EXIT_FAILURE);
    }
    else
        cout << "Input file is " << filename << endl;

    // Read input one line at a time, storing the matrix
    vector<vector<int>> file_matrix;
    vector<int> file_line;
    bool edge_weight_section = false;
    while (getline(inFile, line))
    {
        if (line == "EDGE_WEIGHT_SECTION")
        { // skip unnecessary headers
            edge_weight_section = true;

            // read the dimension line
            getline(inFile, line);
            instance_size = stoi(line);
        }
        else if (line == "EOF")
        {
            break;
        }
        else if (edge_weight_section)
        { // read cost/precedence matrix
            stringstream sstream;
            sstream << line;
            string weight;
            int weight_num;
            while (sstream >> weight)
            {
                stringstream(weight) >> weight_num;
                file_line.push_back(weight_num);
            }
            file_matrix.push_back(file_line);
            file_line.clear();
        }
    }

    if (instance_size != int(file_matrix.size()))
    {
        std::cerr << "Error: Instance Size Ambiguous. Expected " << instance_size << ". Was " << file_matrix.size() << std::endl;
        for (int i = 0; i < (int)file_matrix.size(); i++)
        {
            for (int j = 0; i < (int)file_matrix[i].size(); i++)
            {
                std::cerr << file_matrix[i][j] << " ";
            }
            std::cerr << std::endl;
        }

        exit(EXIT_FAILURE);
    }

    cost_graph = vector<vector<edge>>(instance_size);
    hungarian_graph = vector<vector<edge>>(instance_size);
    dependency_graph = vector<vector<int>>(instance_size);

    for (int i = 0; i < (int)file_matrix.size(); i++)
    {
        int j = 0;
        for (auto edge_weight : file_matrix[i])
        {
            if (edge_weight < 0)
            {
                cost_graph[i].push_back(edge(i, j, file_matrix[j][i]));
                hungarian_graph[i].push_back(edge(i, j, -1));
                dependency_graph[j].push_back(i);
            }
            else
            {
                cost_graph[i].push_back(edge(i, j, edge_weight));
                hungarian_graph[i].push_back(edge(i, j, edge_weight));
            }
            j++;
        }
    }

    return;
}

void solver::transitive_redundancy()
{
    in_degree = std::vector<vector<edge>>(instance_size);
    for (int i = 0; i < instance_size; i++)
    {
        for (long unsigned int k = 0; k < dependency_graph[i].size(); k++)
        {
            int c = dependency_graph[i][k];
            in_degree[c].push_back(edge(i, c, hungarian_graph[i][c].weight));
        }
    }

    for (int i = 0; i < instance_size; i++)
    {
        vector<edge> preceding_nodes;
        for (int k = 0; k < (int)dependency_graph[i].size(); k++)
            preceding_nodes.push_back(edge(i, dependency_graph[i][k], -1));
        unordered_set<int> expanded_nodes;
        for (int j = 0; j < (int)preceding_nodes.size(); ++j)
        {
            vector<edge> st;
            st.push_back(preceding_nodes[j]);
            while (!st.empty())
            {
                edge dependence_edge = st.back();
                st.pop_back();
                if (dependence_edge.src != i)
                {
                    hungarian_graph[dependence_edge.dst][i].weight = -1;
                    hungarian_graph[i][dependence_edge.dst].weight = -1;
                    expanded_nodes.insert(dependence_edge.dst);
                }

                for (int dest : dependency_graph[dependence_edge.dst])
                {
                    if (expanded_nodes.find(dest) == expanded_nodes.end())
                    {
                        st.push_back(edge(dependence_edge.dst, dest, -1));
                        expanded_nodes.insert(dest);
                    }
                }
            }
        }
    }

    for (int i = 0; i < instance_size; i++)
    {
        const vector<edge> preceding_nodes = in_degree[i];
        unordered_set<int> expanded_nodes;
        for (int j = 0; j < (int)preceding_nodes.size(); ++j)
        {
            vector<edge> st;
            st.push_back(preceding_nodes[j]);
            while (!st.empty())
            {
                edge dependence_edge = st.back();
                st.pop_back();
                if (dependence_edge.src != i)
                {
                    hungarian_graph[i][dependence_edge.dst].weight = -1;
                    hungarian_graph[dependence_edge.dst][i].weight = -1;
                    expanded_nodes.insert(dependence_edge.dst);
                }
                for (const edge &e : in_degree[dependence_edge.dst])
                {
                    if (expanded_nodes.find(e.dst) == expanded_nodes.end())
                    {
                        st.push_back(e);
                        expanded_nodes.insert(e.dst);
                    }
                }
            }
        }
    }

    return;
}

size_t solver::transitive_closure(vector<vector<int>> &isucc_graph)
{
    size_t node_num = 0;
    deque<int> dfs_queue;
    boost::dynamic_bitset<> taken_arr = boost::dynamic_bitset<>(instance_size);
    boost::dynamic_bitset<> finished_arr = boost::dynamic_bitset<>(instance_size);
    vector<boost::dynamic_bitset<>> transitive_closure_map = vector<boost::dynamic_bitset<>>(instance_size);
    for (int i = 0; i < instance_size; i++)
        transitive_closure_map[i] = boost::dynamic_bitset<>(instance_size);
    int selected_node = 0;

    for (int i = 0; i < instance_size; i++)
    {
        // cout << "current node is " << i << endl;
        transitive_closure_map[i][i] = true;
        dfs_queue.push_back(i);
        while (!dfs_queue.empty())
        {
            selected_node = dfs_queue.back();
            dfs_queue.pop_back();
            for (int node : isucc_graph[selected_node])
            {
                if (!taken_arr[node])
                {
                    taken_arr[node] = true;
                    dfs_queue.push_back(node);
                }
                if (!transitive_closure_map[i][node])
                    transitive_closure_map[i][node] = true;
                if (!transitive_closure_map[node][i])
                    transitive_closure_map[node][i] = true;
            }
        }
        taken_arr.reset();
    }

    for (int i = 0; i < instance_size; i++)
    {
        for (int k = 0; k < instance_size; k++)
        {
            if (transitive_closure_map[i][k])
                node_num++;
        }
    }
    // print_transitive_map();
    return node_num;
}

vector<int> solver::nearest_neighbor(vector<int> *partial_solution)
{
    vector<int> solution;
    int current_node;
    int solution_cost = 0;
    bool visit_arr[instance_size];
    int depCnt_arr[instance_size];
    vector<vector<edge>> sorted_costgraph = cost_graph;
    sort_weight(sorted_costgraph);
    memset(depCnt_arr, 0, instance_size * sizeof(int));

    for (int i = 0; i < instance_size; i++)
    {
        visit_arr[i] = false;
        for (long unsigned int k = 0; k < dependency_graph[i].size(); k++)
        {
            depCnt_arr[dependency_graph[i][k]]++;
        }
    }

    current_node = 0;
    for (auto node : *partial_solution)
    {
        visit_arr[node] = true;
        for (long unsigned int i = 0; i < dependency_graph[node].size(); i++)
        {
            depCnt_arr[dependency_graph[node][i]]--;
        }
        solution.push_back(node);
        solution_cost += cost_graph[current_node][node].weight;
        current_node = node;
    }

    int num = solution.size();

    while (num < instance_size)
    {
        bool taken = false;
        for (auto node : sorted_costgraph[current_node])
        { // this is actually taking an edge, since each element of the cost graph is an edge
            if (!visit_arr[node.dst] && !depCnt_arr[node.dst])
            {
                current_node = node.dst;
                solution_cost += node.weight;
                solution.push_back(current_node);
                num++;
                visit_arr[node.dst] = true;
                taken = true;
                break;
            }
        }
        if (!taken)
        {
            std::cout << "Error generating Nearest Neighbor Heuristic" << std::endl;
            std::cout << "current node is " << current_node << std::endl;
            for (auto node : sorted_costgraph[current_node])
            {
                cout << node.dst << "," << visit_arr[node.dst] << "," << depCnt_arr[node.dst] << endl;
            }
            exit(EXIT_FAILURE);
        }
        for (long unsigned int i = 0; i < dependency_graph[current_node].size(); i++)
        {
            depCnt_arr[dependency_graph[current_node][i]]--;
        }
    }

    best_cost = solution_cost;
    return solution;
}

void solver::sort_weight(vector<vector<edge>> &graph)
{
    int size = graph.size();
    for (int i = 0; i < size; i++)
    {
        stable_sort(graph[i].begin(), graph[i].end(), compare_edge);
    }
    return;
}

int solver::get_maxedgeweight()
{
    int max = 0;
    for (int i = 0; i < instance_size; i++)
    {
        for (edge edge : cost_graph[i])
        {
            int weight = edge.weight;
            if (weight > max)
                max = weight;
        }
    }
    return max;
}

vector<vector<int>> solver::get_cost_matrix(int max_edge_weight)
{
    vector<vector<int>> matrix(instance_size);
    for (int i = 0; i < instance_size; ++i)
    {
        matrix[i] = vector<int>(instance_size, max_edge_weight * 2);
    }

    for (vector<edge> edge_list : hungarian_graph)
    {
        for (auto edge : edge_list)
        {
            int i = edge.src;
            int k = edge.dst;
            int weight = edge.weight;
            if (weight != -1 && i != k)
            {
                matrix[i][k] = weight * 2;
            }
        }
    }

    return matrix;
}

bool solver::split_level_check(deque<sop_state> *solver_container)
{
    // unsigned target_level = solver_container->front().current_path.size();
    // for (auto node : *solver_container) {
    //     if (node.current_path.size() != target_level) return true;
    // }
    // return false;
    // this is unnecessary
    // Since nodes are only added to the back, and only taken from the front, you can just check the front vs. back in constant time
    return solver_container->front().current_path.size() != solver_container->back().current_path.size();
}

bool solver::enumeration_pre_check(path_node &active_node)
{
    if (active_node.lower_bound >= best_cost
        // || stop_init
        // || (enable_threadstop && active_node.his_entry != NULL
        //                       && active_node.his_entry->Entry.load().prefix_cost < active_node.partial_cost) //TODO: thread stopping
    )
    {
        // if (enable_progress_estimation) //pruning due to enumeration-time backtracking
        //     estimated_trimmed_percent[thread_id] += active_node.current_node_value; //add the value of this node you are trimming
        // PROGRESS
        work_remaining[thread_id] -= active_node.current_node_value;

        // cur_active_tree.incre_children_cnt(Allocator);
        // if (active_node.his_entry != NULL && active_node.his_entry->active_threadID == thread_id) {
        //     active_node.his_entry->explored = true;
        // }
        // DIAGNOSTIC: view path
        // cout << "at pre check Pruned node " << active_node.sequence.back() << endl;
        return true;
    }
    return false;
}

void solver::prune(int source_node, int taken_node)
{
    // DIAGNOSTIC: view path
    // cout << "Pruning node " << taken_node << endl;
    //  if (enable_progress_estimation)
    //      estimated_trimmed_percent[thread_id] += dest.current_node_value; //add the value of this node you are trimming

    problem_state.current_path.pop_back();
    problem_state.current_cost -= cost_graph[source_node][taken_node].weight;
    problem_state.history_key.first[taken_node] = false;
    problem_state.history_key.second = source_node;
}

int solver::dynamic_hungarian(int src, int dst)
{
    problem_state.hungarian_solver.fix_row(src, dst);
    problem_state.hungarian_solver.fix_column(dst, src);
    problem_state.hungarian_solver.solve_dynamic();
    // number_of_lb++;
    int lb = problem_state.hungarian_solver.get_matching_cost() / 2; // based on the conversion between MCPM and SOP, the actual lower bound is the Hungarian solution / 2
    problem_state.hungarian_solver.undue_row(src, dst);
    problem_state.hungarian_solver.undue_column(dst, src);
    return lb;
}

bool solver::history_utilization(Key &key, int cost, int *lowerbound, bool *found, HistoryNode **entry)
{
    HistoryNode *history_node = history_table.retrieve(key, problem_state.current_path.size());

    if (history_node == NULL)
        return true;

    *found = true;
    *entry = history_node;
    HistoryContent content = history_node->entry.load();
    *lowerbound = content.lower_bound;

    int target_ID = history_node->active_threadID; // find whoever was working in this subspace
    int imp = content.prefix_cost - cost;
    if (cost >= content.prefix_cost)
        return false;

    if (!history_node->explored)
    { // TODO: thread stopping
        if (enable_threadstop && active_threads > 0 && target_ID != thread_id)
        { // then issue thread stop request, since this path is superior
            if (!thread_requests[target_ID].has_request || thread_requests[target_ID].request.target_depth > (int)problem_state.current_path.size())
            {
                thread_requests[target_ID].lock.lock();
                if (!thread_requests[target_ID].has_request || thread_requests[target_ID].request.target_depth > (int)problem_state.current_path.size()) // extra validation
                {
                    thread_stop_requested++;
                    thread_requests[target_ID].request = request_packet(problem_state.current_path.back(), (int)problem_state.current_path.size(),
                                                                        content.prefix_cost, target_ID, key.first);
                    thread_requests[target_ID].has_request = true;
                }
                thread_requests[target_ID].lock.unlock();
            }
        }
    }

    if (imp <= content.lower_bound - best_cost)
    {
        return false;
    }
    history_node->entry.store({cost, content.lower_bound - imp});
    *lowerbound = content.lower_bound - imp;
    *entry = history_node;
    history_node->explored = false;
    history_node->active_threadID = thread_id;

    return true;
}

void solver::push_to_history_table(Key &key, int lower_bound, HistoryNode **entry, bool backtracked)
{
    if (entry == NULL)
        history_table.insert(key, problem_state.current_cost, lower_bound, thread_id, backtracked, problem_state.current_path.size(), instance_size / 3);
    else
        *entry = history_table.insert(key, problem_state.current_cost, lower_bound, thread_id, backtracked, problem_state.current_path.size(), instance_size / 3);
    return;
}

static int stolen_from;
/* BEGIN WORK STEALING*/
// WORKSTEALING
bool solver::workload_request()
{
    if (thread_requests[thread_id].has_request)
    {
        thread_requests[thread_id].lock.lock();
        if (thread_requests[thread_id].has_request)
        {
            thread_requests[thread_id].has_request = false;
        }
        thread_requests[thread_id].lock.unlock();
    }
    if (work_remaining[thread_id] != 0)
    {

        // cout << "ERROR!!! thread " << thread_id << " at workstealing with " << work_remaining[thread_id] << " work remaining" << endl;
    }
    local_pools->set_pool_depth(thread_id, INT32_MAX);
    if (!global_pool.empty())
    {
        global_pool_lock.lock();
        if (!global_pool.empty())
        {
            problem_state = generate_solver_state(global_pool.back());
            local_pools->set_pool_depth(thread_id, 0);
            global_pool.pop_back();
            if (global_pool.empty())
            {
                /**
                 * if global pool is empty, we want to ignore the blocking & deletion of the table
                 * updating the mem_limit to 0.9 (90% of the memory limit)
                 *
                 * updating the number of groups to 1 to treat the history table as a single subspace
                 * and prevent the blocking and the deletion of the history table
                 */
                if (enable_heuristic)
                {
                    cout << "Start treating the history table as a single subspace" << endl;
                    mem_limit = 0.9;
                    number_of_groups = 1;
                }

                cout << "GLOBAL POOL EMPTY at time: " << main_timer.get_time_seconds() << endl;
                gp_remaining = global_pool.size();

                // print_workdone();

            } // updating the variable to track how many of the threads completed the work assigned to them from the primary subspace
            gp_remaining = global_pool.size();
            global_pool_lock.unlock();
            return true;
        }
        global_pool_lock.unlock();
    }

    active_threads--;
    if (enable_workstealing)
    {
        timer t;
        path_node new_node;
        int misses = 0;
        while (true)
        {
            if (active_threads <= 0)
            {
                time_workstealing = time_workstealing + t.get_time_seconds();
                return false;
            }
            int target = local_pools->choose_victim(thread_id, work_remaining, stolen_from);

            if (target == -1)
            {
                misses++;
                stolen_from = 0;
                continue;
            }
            steal_attempts[target]++;
            if (local_pools->pop_from_zero_list(target, new_node, thread_id))
            {
                work_remaining[target] -= new_node.current_node_value;
                problem_state = generate_solver_state(new_node);
                problem_state.work_above = new_node.current_node_value;
                active_threads++;
                times_work_stolen++;
                steal_success[thread_id]++;
                time_workstealing = time_workstealing + t.get_time_seconds();
                steal_misses += misses;
                // steal_times_lock.lock();
                // steal_times.push_back(main_timer.get_time_seconds());
                // steal_times_lock.unlock();
                return true;
            }
            stolen_from = stolen_from | (1 << target);
            misses++;
        }
    }
    return false;
}

sop_state solver::generate_solver_state(path_node &subproblem)
{
    sop_state state = default_state;

    work_remaining[thread_id] = subproblem.current_node_value;
    int cur_node = subproblem.sequence.front();
    int taken_node = -1;
    int size = subproblem.sequence.size();
    state.current_path.push_back(0);
    state.taken_arr[cur_node] = true;
    for (int vertex : dependency_graph[cur_node])
        state.depCnt[vertex]--;

    for (int k = 1; k < size; k++)
    { // begin after the first, trivial, node
        taken_node = subproblem.sequence[k];
        state.current_path.push_back(taken_node);
        state.taken_arr[taken_node] = true;
        for (int vertex : dependency_graph[taken_node])
            state.depCnt[vertex]--;
        state.current_cost += cost_graph[cur_node][taken_node].weight;
        state.hungarian_solver.fix_row(cur_node, taken_node);
        state.hungarian_solver.fix_column(taken_node, cur_node);
        cur_node = taken_node;
    }

    // generate history table key
    boost::dynamic_bitset<> bit_vector(instance_size, false);
    for (auto node : state.current_path)
    {
        bit_vector[node] = true;
    }
    int last_element = state.current_path.back();
    state.history_key = make_pair(bit_vector, last_element);

    // cur_active_tree.generate_path(sequence_node.partial_active_path);
    // cur_active_tree.set_threadID(thread_id, thread_total);
    // if (current_hisnode != NULL && !current_hisnode->explored) {
    //     current_hisnode->active_threadID = thread_id;
    // }

    state.origin_node = subproblem.origin_node;
    state.lower_bound = subproblem.lower_bound;
    // solvers[thread_cnt].problem_state.initial_depth = solvers[thread_cnt].problem_state.cur_solution.size();
    // solvers[thread_cnt].problem_state.current_node_value = problem.current_node_value; //for progress estimation

    // std::cout << "New SOP state" << std::endl;
    // print_state(state);
    // exit(EXIT_FAILURE);

    return state;
}

/* BEGIN THREAD STOPPING */

/* END THREAD STOPPING */

// DIAGNOSTIC
/* BEGIN DIAGNOSTIC FUNCTIONS */

void solver::print_state(sop_state &state)
{
    std::cout << "Path: ";
    for (int i = 0; i < (int)state.current_path.size(); i++)
        std::cout << state.current_path[i] << " ";
    std::cout << std::endl;
    std::cout << "Cost: " << state.current_cost << std::endl;
    std::cout << "Lower Bound: " << state.lower_bound << std::endl;
    std::cout << "Origin: " << state.origin_node << std::endl;
    std::cout << "Key: ";
    for (int i = 0; i < instance_size; i++)
        std::cout << state.history_key.first[i] << " ";
    std::cout << " -- " << state.history_key.second << std::endl;
    std::cout << "Taken: ";
    for (int i = 0; i < (int)state.taken_arr.size(); i++)
        std::cout << state.taken_arr[i] << " ";
    std::cout << std::endl;
    std::cout << "Unfulfilled Dependencies: ";
    for (int i = 0; i < (int)state.depCnt.size(); i++)
        std::cout << state.depCnt[i] << " ";
    std::cout << std::endl;
}

bool solver::check_stop_request(std::pair<boost::dynamic_bitset<>, int> history_key, vector<int> sequence, bool *prefixKeyMatched)
{
    if (thread_requests[thread_id].has_request)
    {
        thread_requests[thread_id].lock.lock();
        if (thread_requests[thread_id].has_request)
        {
            request_packet rp = thread_requests[thread_id].request;
            // if (rp.target_thread == thread_id)
            // {
            if (rp.target_depth <= sequence.size())
            {
                int current_cost = problem_state.current_cost;
                for (int i = 0; i < sequence.size() - rp.target_depth; i++)
                {
                    current_cost -= cost_graph[cost_graph.size() - i - 1][cost_graph.size() - i].weight;
                }
                if (rp.target_last_node == sequence[rp.target_depth - 1] && current_cost >= rp.target_prefix_cost)
                {
                    thread_stop_check++;
                    if (rp.key == history_key.first) // will only occur when the size of the target_depth and sequence size is same
                    {
                        // tracking the pruning at current depth
                        // history_table_pruning_success[sequence.size()]++;
                        thread_stopped_successfully++;
                        thread_requests[thread_id].has_request = false;
                        thread_requests[thread_id].lock.unlock();
                        return true; // Indicate that a stop request was found and handled
                    }
                    else if (rp.key == generate_history_key(sequence, rp.target_depth)) // will only occur when the size of the target_depth and sequence size is different
                    {
                        // tracking the pruning at current depth
                        // history_table_pruning_success[sequence.size()]++;
                        thread_stopped_successfully++;
                        *prefixKeyMatched = true;
                        thread_requests[thread_id].lock.unlock();
                        return true; // Indicate that a stop request was found and handled
                    }
                }
            }
            thread_requests[thread_id].has_request = false;
            thread_requests[thread_id].lock.unlock();
            return false; // Indicate that a stop request was found and handled
            // }
            // else
            // {
            //     cout << "thread id mismatch \n";
            // }
            thread_requests[thread_id].lock.unlock();
        }
    }
    return false;
}

boost::dynamic_bitset<> solver::generate_history_key(const vector<int> &sequence, int depth)
{
    boost::dynamic_bitset<> temp_history_key(instance_size);
    for (int i = 0; i < depth; ++i)
    {
        temp_history_key.set(sequence[i]);
    }
    return temp_history_key;
}

/* END DIAGNOSTIC FUNCTIONS */

// bool shared_pool_sort(const path_node& src, const path_node& dest) {
//     if (src.sequence.size() == dest.sequence.size()) return src.lower_bound > dest.lower_bound;
//     return src.sequence.size() < dest.sequence.size();
// }
// bool bound_sort(const node& src,const node& dest) {
//     if (src.lb == dest.lb) return src.nc > dest.nc;
//     return src.lb > dest.lb;
// }
// bool nearest_sort(const node& src,const node& dest) {
//     if (src.nc == dest.nc) return src.lb > dest.lb;
//     return src.nc > dest.nc;
// }
// DEPRICATED; you don't need to sort the GPQ during splitting since the order doesn't matter yet, you have to split evenly so all nodes are of equal depth
// bool split_sort(const sop_state& src, const sop_state& dest) {
//     if (src.current_path.size() == dest.cur_solution.size()) return src.load_info < dest.load_info;
//     return src.cur_solution.size() < dest.cur_solution.size();
// }
