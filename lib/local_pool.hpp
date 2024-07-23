#ifndef LOCAL_H
    #define LOCAL_H

    #include <deque> //for the local pool structure
    #include <queue>
    #include <iostream>
    #include "synchronization.hpp"
    #include "graph.hpp"

    /* The local pool construct, consisting of pools of nodes for each thread. In 
        each thread, the pool is organized by depth, so that nodes are stolen from 
        only the shallowest part of the pool, and added at the deepest level. */
    class local_pool {
        private:
            int thread_count;
            int current_target = 0;
            std::vector<spin_lock> locks;
            std::vector<std::deque<std::deque<path_node>>> pools;
            std::vector<int> depths;
        public:
            local_pool(int thread_count);
            void add_to_depth_queue(int thread);
            /*Grabs a node from the shallowest / zero pool*/
            bool pop_from_zero_list(int thread_number, path_node &result_node, int stealing_thread);
            /*Grabs a node from the deepest / active pool*/
            bool pop_from_active_list(int thread_number, path_node &result_node);
            bool peek_from_active_list(int thread_number, path_node &result_node);
            /*Pushes new list to the back of the local pool*/
            void push_list(int thread_number, std::deque<path_node> list);
            /*removes */
            void pop_active_list(int thread_number);
            /* Determines if a specific thread's local pool is completely empty. */
            bool out_of_work(int thread_number);
            /* Returns a thread number of the best victim, other than you, for workstealing. 
                thread_number - this thread's number, to ensure you aren't recommended to steal from yourself
                Return - the thread number of the thread to steal from */
            int choose_victim(int thread_number,std::vector<std::atomic<unsigned long long>>& work_remaining, int stolen_from);
            //diagnostic
            int active_pool_size(int thread_number); 
            //diagnostic
            void print();
            //sets the relative depth of the pool
            void set_pool_depth(int thread_id, int depth);
    };

#endif