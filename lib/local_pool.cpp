#include "local_pool.hpp"

    local_pool::local_pool(int thread_count){
        locks = std::vector<spin_lock>(thread_count);
        pools = std::vector<std::deque<std::deque<path_node>>>(thread_count);
        depths = std::vector<int>(thread_count);
        this->thread_count = thread_count;
        for (int i = 0; i < thread_count; i++){
            std::vector<int> temp;
            temp.push_back(i);
            temp.push_back(0);
            depth_queue.push(temp);
        }
    }

    void local_pool::add_to_depth_queue(int thread){
        queue_lock.lock();
        std::vector<int> temp;
        temp.push_back(thread);
        temp.push_back(depths[thread]);
        depth_queue.push(temp);
        queue_lock.unlock();
    }      

    bool local_pool::pop_from_zero_list(int thread_number, path_node &result_node, int stealing_thread){
        if (pools[thread_number].size() <= 1){
            add_to_depth_queue(thread_number);
            return false;
        }
        locks[thread_number].lock();

        while(pools[thread_number].front().empty() && pools[thread_number].size() > 1) {
            pools[thread_number].pop_front();
            depths[thread_number]++;    
        }


        if (pools[thread_number].size() <= 1){
            locks[thread_number].unlock();
            add_to_depth_queue(thread_number);
            return false;
        }
        
        result_node = pools[thread_number].front().back();
        pools[thread_number].front().pop_back();
        depths[stealing_thread] = depths[thread_number] + 1;

        if (pools[thread_number].front().empty()){
            pools[thread_number].pop_front();
            depths[thread_number]++;
        }
        add_to_depth_queue(thread_number);

    if (pools[thread_number].front().empty())
    {
        pools[thread_number].pop_front();
    }

    locks[thread_number].unlock();
    return true;
};

bool local_pool::pop_from_active_list(int thread_number, path_node &result_node)
{

    if (pools[thread_number].size() == 0)
        return false;
    if (pools[thread_number].size() == 1)
    {
        locks[thread_number].lock();
    }

    if (pools[thread_number].size() == 0 || pools[thread_number].back().empty())
    {
        if (pools[thread_number].size() == 1)
        {
            locks[thread_number].unlock();
        }
        return false;
    }

    result_node = pools[thread_number].back().back();
    pools[thread_number].back().pop_back();

    if (pools[thread_number].size() == 1)
    {
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

    // int local_pool::choose_victim(int thread_number, std::vector<std::atomic<unsigned long long>>& work_remaining){
    //     int best = depths[0];
    //     std::deque<int> set;
    //     for(int i = 0; i < depths.size(); i++ ){
    //         if(pools[i].size() == 0){
    //             //depths[i] = INT32_MAX;
    //             continue;
    //         }
    //         if(depths[i] == best){
    //             set.push_front(i);
    //         }
    //         if(depths[i] < best){
    //             set.clear();
    //             set.push_front(i);
    //         }
    //     }
    //     if(set.size() == 0) return -1;
    //     return set[rand() % set.size()];
    // }

    // int local_pool::choose_victim(int thread_number, std::vector<std::atomic<unsigned long long>>& work_remaining){
    //     int target = rand() % 30;
    //     while(target == thread_number) target = rand() % 30;
    //     return target;
    // }

    int local_pool::choose_victim(int thread_number, std::vector<std::atomic<unsigned long long>>& work_remaining){
        queue_lock.lock();
        if(depth_queue.empty()){ 
            queue_lock.unlock();
            return -1;
        } 
        int temp = depth_queue.top().at(0);
        depth_queue.pop();
        
        queue_lock.unlock();
        return temp;
    }

    void local_pool::set_pool_depth(int thread, int depth){
        depths[thread] = depth;
    }

    int local_pool::active_pool_size(int thread_number) { //TODO: this is not strictly necessary
        return pools[thread_number].back().size();
    }

    void local_pool::print(){
        for(int i = 0; i < pools.size(); i++){
            std::cout << pools[i].size()<<", ";
        }
        std::cout << std::endl;
    }
