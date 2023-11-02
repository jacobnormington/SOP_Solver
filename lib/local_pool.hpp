#ifndef LOCAL_H
    #define LOCAL_H

    #include <deque> //for the local pool structure

    #include "synchronization.hpp"
    #include "graph.hpp"

    //used to pick the next thread to be stolen from
    class workstealing_targeter {
        //TODO: does this actually need to be its own module?
        public:
            /* Setup a new targeter.
                size - the number of threads alloted to B&B enumeration */
            workstealing_targeter(int size){
                size = size;
            }
            /* Determine which thread should be stolen from. */
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

    /* The local pool construct, consisting of pools of nodes for each thread. In 
        each thread, the pool is organized by depth, so that nodes are stolen from 
        only the shallowest part of the pool, and added at the deepest level. */
    class local_pool {
        //TODO: should we instead use a priority queue here, a heap implementation instead of naively just a list of lists?
        private:
            std::vector<spin_lock> locks;
            std::vector<std::deque<std::deque<path_node>>> pools;
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
    };

#endif