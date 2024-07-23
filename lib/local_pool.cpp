#include "local_pool.hpp"

    local_pool::local_pool(int thread_count){
        locks = std::vector<spin_lock>(thread_count);
        pools = std::vector<std::deque<std::deque<path_node>>>(thread_count);
        depths = std::vector<int>(thread_count);
        this->thread_count = thread_count;
    }

    bool local_pool::pop_from_zero_list(int thread_number, path_node &result_node, int stealing_thread){
        if (pools[thread_number].size() <= 1){
            return false;
        }
        locks[thread_number].lock();

        while(pools[thread_number].front().empty() && pools[thread_number].size() > 1) {
            pools[thread_number].pop_front();
            depths[thread_number]++;    
        }


        if (pools[thread_number].size() <= 1){
            locks[thread_number].unlock();
            return false;
        }
        
        result_node = pools[thread_number].front().back();
        pools[thread_number].front().pop_back();
        depths[stealing_thread] = depths[thread_number] + 1;

        if (pools[thread_number].front().empty()){
            pools[thread_number].pop_front();
            depths[thread_number]++;
        }

        locks[thread_number].unlock();
        return true;
    };

    bool local_pool::pop_from_active_list(int thread_number, path_node &result_node){

        if (pools[thread_number].size() == 0)
            return false;

        if (pools[thread_number].size() == 0 || pools[thread_number].back().empty()){
            if (pools[thread_number].size() == 1){
                locks[thread_number].unlock();
            }
            return false;
        }

        result_node = pools[thread_number].back().back();
        pools[thread_number].back().pop_back();

        return true;
    };

    bool local_pool::peek_from_active_list(int thread_number, path_node &result_node){

        if (pools[thread_number].size() == 0)
            return false;

        if (pools[thread_number].size() == 0 || pools[thread_number].back().empty()){
            if (pools[thread_number].size() == 1){
                locks[thread_number].unlock();
            }
            return false;
        }

        result_node = pools[thread_number].back().back();
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

    int local_pool::choose_victim(int thread_number, std::vector<std::atomic<unsigned long long>>& work_remaining, int stolen_from){
        unsigned long long max_value = 0;
        int max_id = -1;
        bool flag = false;
        for(int i = 0; i < thread_count; i++){
            if((stolen_from & (1 << i)) != 0 || i == thread_number)
                continue;
            locks[i].lock();
            unsigned long long node_value = 0;
            if(pools[i].size() > 1 && pools[i].front().size() != 0) {
                node_value = pools[i].front().back().current_node_value;
            }
            locks[i].unlock();
            if(node_value > max_value){
                max_value = node_value;
                max_id = i;
                continue;
            }
            if(max_value == 0 && (!flag || work_remaining[i] > work_remaining[max_id])){
                max_value = node_value;
                max_id = i;
                flag = true;
            }
        }
        return max_id;
    }

    // int local_pool::choose_victim(int thread_number, std::vector<std::atomic<unsigned long long>>& work_remaining, int stolen_from){
    //     double max_value = -1;
    //     int max_id = -1;
    //     //std::cout << stolen_from << std::endl;
    //     for(int i = 0; i < thread_count; i++){
    //         if((stolen_from & (1 << i)) != 0 )
    //             continue;
    //         if(i == thread_number)
    //             continue;
    //         if(work_remaining[i] / (depths[i] + 1) > max_value){
    //             max_value = work_remaining[i] / (depths[i] + 1);
    //             max_id = i;
    //             continue;
    //         }
    //         // if(work_remaining[i] == max_value && depths[i] < depths[max_id]){
    //         //     max_value = work_remaining[i];
    //         //     max_id = i;
    //         // }
    //     }
    //     return max_id;
    // }

    // int local_pool::choose_victim(int thread_number, std::vector<std::atomic<unsigned long long>>& work_remaining, int a){
    //     int target = rand() % 30;
    //     while(target == thread_number) target = rand() % 30;
    //     return target;
    // }

    int local_pool::active_pool_size(int thread_number) { //TODO: this is not strictly necessary
        return pools[thread_number].back().size();
    }

    void local_pool::print(){
        for(int i = 0; i < pools.size(); i++){
            std::cout << pools[i].size()<<", ";
        }
        std::cout << std::endl;
    }

    void local_pool::set_pool_depth(int thread_id, int depth){
        depths[thread_id] = depth;
    }