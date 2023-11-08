#include "local_pool.hpp"

    bool local_pool::pop_from_zero_list(int thread_number, path_node &result_node){
        if (pools[thread_number].size() <= 1)
            return false;
        locks[thread_number].lock();

        if (pools[thread_number].size() <= 1 || pools[thread_number].front().empty()){
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

    bool local_pool::pop_from_active_list(int thread_number, path_node &result_node){
        if (pools[thread_number].size() == 0)
            return false;
        if (pools[thread_number].size() == 1){
            locks[thread_number].lock();
        }

        if (pools[thread_number].size() == 0 || pools[thread_number].back().empty()){
            if (pools[thread_number].size() == 1){
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

    void local_pool::push_list(int thread_number, std::deque<path_node> list){
        locks[thread_number].lock();

        pools[thread_number].push_back(list);

        locks[thread_number].unlock();
    };

    void local_pool::pop_active_list(int thread_number){
        locks[thread_number].lock();
        pools[thread_number].pop_back();
        locks[thread_number].unlock();
    };

    bool local_pool::out_of_work(int thread_number){
        return pools[thread_number].size() == 0;
    };

    int local_pool::choose_victim(int thread_number){
        workstealing_lock.lock();
        while (current_target != thread_number)
            (current_target + 1) % thread_count;
        int return_value = current_target;
        current_target = (current_target + 1) % thread_count;
        workstealing_lock.unlock();
        return return_value;
    }
