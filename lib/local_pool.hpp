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
        //TODO: should we instead use a priority queue here, a heap implementation instead of naively just a list of lists?
        private:
            std::vector<spin_lock> locks;
            std::vector<std::deque<std::deque<path_node>>> pools;
            
            //workstealing variables
            int thread_count;
            int current_target = 0;
            spin_lock workstealing_lock;
            spin_lock queue_lock;
            

            struct my_comparator
            {
                // queue elements are vectors so we need to compare those
                bool operator()(std::vector<int> const& a, std::vector<int> const& b) const
                {
                    return a[1] > b[1];
                }   
            };

            std::priority_queue<std::vector<int>, std::vector<std::vector<int>>, my_comparator> depth_queue;


        public:
        std::vector<int> depths;
            local_pool(int thread_count);
            void add_to_depth_queue(int thread);
            /*Grabs a node from the shallowest / zero pool*/
            bool pop_from_zero_list(int thread_number, path_node &result_node, int stealing_thread);
            /*Grabs a node from the deepest / active pool*/
            bool pop_from_active_list(int thread_number, path_node &result_node);
            /*Pushes new list to the back of the local pool*/
            void push_list(int thread_number, std::deque<path_node> list);
            /*removes */
            void pop_active_list(int thread_number);
            /* Determines if a specific thread's local pool is completely empty. */
            bool out_of_work(int thread_number);
            /* Returns a thread number of the best victim, other than you, for workstealing. 
                thread_number - this thread's number, to ensure you aren't recommended to steal from yourself
                Return - the thread number of the thread to steal from */
            int choose_victim(int thread_number,std::vector<std::atomic<unsigned long long>>& work_remaining);

            void set_pool_depth(int thread, int depth);

            int active_pool_size(int thread_number); //diagnostic

            void print();
    };

#endif