#include "local_pool.hpp"
#include <iostream>

    bool local_pool::pop_from_zero_list(int thread_number, path_node &result_node){
        if (pools[thread_number].size() <= 1)
            return false;
        locks[thread_number].lock();

        while(pools[thread_number].front().empty() && pools[thread_number].size() > 1) pools[thread_number].pop_front();


        if (pools[thread_number].size() <= 1){
            locks[thread_number].unlock();
            return false;
        }
        
        result_node = pools[thread_number].front().back();
        pools[thread_number].front().pop_back();

        if (pools[thread_number].front().empty()){
            pools[thread_number].pop_front();
        }

        locks[thread_number].unlock();
        return true;
    };

    bool local_pool::pop_from_active_list(int thread_number, path_node &result_node){

        if (pools[thread_number].size() == 0)
            return false;
        // if (pools[thread_number].size() == 1){
        //     locks[thread_number].lock();
        // }

        if (pools[thread_number].size() == 0 || pools[thread_number].back().empty()){
            if (pools[thread_number].size() == 1){
                locks[thread_number].unlock();
            }
            return false;
        }

        result_node = pools[thread_number].back().back();
        pools[thread_number].back().pop_back();

        // if (pools[thread_number].size() == 1){
        //     locks[thread_number].unlock();
        // }
        return true;
    };

    void local_pool::push_list(int thread_number, std::deque<path_node> list){
        locks[thread_number].lock();

        pools[thread_number].push_back(list);

        locks[thread_number].unlock();
    };

    void local_pool::pop_active_list(int thread_number){
        locks[thread_number].lock();
        if(pools[thread_number].size() > 0)
            pools[thread_number].pop_back();
        locks[thread_number].unlock();
    };

    bool local_pool::out_of_work(int thread_number){
        return pools[thread_number].size() == 0;
    };

    // int local_pool::choose_victim(int thread_number){
    //     workstealing_lock.lock();
    //     do {
    //         current_target = (current_target + 1) % thread_count;
    //     } while (current_target == thread_number); //skip the requesting thread
    //     int return_value = current_target;
    //     workstealing_lock.unlock();
    //     return return_value;
    // }

    // int local_pool::choose_victim(int thread_number, std::vector<std::atomic<unsigned long long>>& work_remaining){
    //     workstealing_lock.lock();
    //     double max_value = 0;
    //     int max_id = -1;
    //     for(int i = 0; i < thread_count; i++){
    //         if(i == thread_number)
    //             continue;
    //         if(work_remaining[i] > max_value){
    //             max_value = work_remaining[i];
    //             max_id = i;
    //         }
    //     }
    //     workstealing_lock.unlock();
    //     return max_id;
    // }

    int local_pool::choose_victim(int thread_number, std::vector<std::atomic<unsigned long long>>& work_remaining){
        int target = rand() % 30;
        while(target == thread_number) target = rand() % 30;
        return target;
    }

    int local_pool::active_pool_size(int thread_number) { //TODO: this is not strictly necessary
        return pools[thread_number].back().size();
    }