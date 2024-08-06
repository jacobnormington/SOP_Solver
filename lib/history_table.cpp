#include "history_table.hpp"

#define MEMORY_RESTRIC 0.95
#define SAMPLE_FREQUENCY 60

Memory_Module::Memory_Module()
{
    bucket_block = new Bucket[BUCKET_BLK_SIZE];
    history_block = (HistoryNode *)malloc(HIS_BLK_SIZE * sizeof(HistoryNode));
    bucket_counter = 0;
    his_node_counter = 0;
}
Memory_Module::~Memory_Module()
{
    // cout << "destructor triggered\n";
    delete[] bucket_block; // Use delete[] to free the array of Buckets
    free(history_block);   // Use free to deallocate memory allocated with malloc
}

Bucket *Memory_Module::get_bucket()
{
    if (bucket_counter >= BUCKET_BLK_SIZE || bucket_block == NULL)
    {
        bucket_block = new Bucket[BUCKET_BLK_SIZE];
        bucket_counter = 0;
    }
    Bucket *bucket = bucket_block + bucket_counter;
    bucket_counter++;
    return bucket;
}

HistoryNode *Memory_Module::retrieve_his_node()
{
    // HistoryNode* node = NULL;

    if (his_node_counter >= HIS_BLK_SIZE || history_block == NULL)
    {
        history_block = (HistoryNode *)malloc(HIS_BLK_SIZE * sizeof(HistoryNode));
        his_node_counter = 0;
    }
    HistoryNode *node = history_block + his_node_counter;
    his_node_counter++;

    return node;
}

History_Table::History_Table(size_t size)
{
    struct sysinfo info;
    if (sysinfo(&info) != 0)
    {
        cout << "can't retrieve system memory info\n";
        exit(EXIT_FAILURE);
    }
    num_buckets = size;
    total_ram = info.totalram;
    max_size = (double)info.freeram * MEMORY_RESTRIC - (size * 4);
    current_size = 0;

    // std::cout << "Max bucket size is " << max_size / 1000000 << " MB" << std::endl;
    // std::cout << "Total Available Memory In OS is " << info.totalram / 1000000 << " MB" << std::endl;

    insert_count = 0;
}

void History_Table::initialize(int thread_num, size_t size, int number_of_groups, int group_size)
{
    num_of_groups = number_of_groups;
    groups_size = group_size;
    block_count.resize(number_of_groups, 0);

    map.resize(number_of_groups);
    memory_allocators.resize(number_of_groups);

    table_lock.resize(number_of_groups);

    blocked_groups.resize(number_of_groups, false);
    is_data_available.resize(number_of_groups, true);

    group_locks = vector<spin_lock>(number_of_groups);

    for (int i = 0; i < number_of_groups; i++)
    {
        map[i].resize(size);
        memory_allocators[i].resize(thread_num);

        table_lock[i] = vector<spin_lock>(size / COVER_AREA + 1);
    }
}

size_t History_Table::get_max_size() { return max_size; }
size_t History_Table::get_current_size() { return current_size; }

unsigned long History_Table::get_free_mem()
{
    struct sysinfo info;
    if (sysinfo(&info) != 0)
    {
        cout << "can't retrieve sys mem info\n";
        exit(1);
    }
    return (double)info.freeram;
}

void History_Table::print_curmem()
{
    std::cout << "Mem exhausted is " << current_size / 1000000 << "MB" << std::endl;
    return;
}

HistoryNode *History_Table::insert(Key &key, int prefix_cost, int lower_bound, unsigned thread_id, bool backtracked, unsigned depth)
{
    int group_index;
    if (depth <= num_of_groups * groups_size)
        group_index = std::ceil(static_cast<double>(depth) / groups_size) - 1;
    else
        group_index = num_of_groups - 1;

    if (!is_data_available[group_index] || blocked_groups[group_index])
    {
        block_count[group_index]++;
        return NULL;
    }

    HistoryNode *node = memory_allocators[group_index][thread_id].retrieve_his_node();

    if (thread_id % 4  == 0)
    {
        insert_count++;
        if (insert_count >= 100000)
        {
            current_size = total_ram - get_free_mem();
            insert_count = 0;
        }
    }

    if (node == NULL)
        return NULL;

    size_t val = hash<boost::dynamic_bitset<>>{}(key.first);
    int bucket = (val + key.second) % num_buckets;

    node->explored = backtracked;
    node->entry.store({prefix_cost, lower_bound});
    node->active_threadID = thread_id;
    // node->depth = depth;
    // node->usage_cnt = 0;

    table_lock[group_index][bucket / COVER_AREA].lock();
    if (!is_data_available[group_index] || blocked_groups[group_index])
    {
        block_count[group_index]++;
        table_lock[group_index][bucket / COVER_AREA].unlock();
        return NULL;
    }
    if (map[group_index][bucket] == NULL)
    {
        map[group_index][bucket] = memory_allocators[group_index][thread_id].get_bucket();
        map[group_index][bucket]->push_back(make_pair(key, node));
        table_lock[group_index][bucket / COVER_AREA].unlock();
        return node;
    }

    map[group_index][bucket]->push_back(make_pair(key, node));
    table_lock[group_index][bucket / COVER_AREA].unlock();
    return node;
}

HistoryNode *History_Table::retrieve(Key &key, int depth)
{
    int group_index;
    if (depth <= num_of_groups * groups_size)
        group_index = std::ceil(static_cast<double>(depth) / groups_size) - 1;
    else
        group_index = num_of_groups - 1;

    if (!is_data_available[group_index])
        return NULL;

    size_t val = hash<boost::dynamic_bitset<>>{}(key.first);
    int bucket = (val + key.second) % num_buckets;

    table_lock[group_index][bucket / COVER_AREA].lock();
    if (!is_data_available[group_index] || map[group_index][bucket] == NULL)
    {
        table_lock[group_index][bucket / COVER_AREA].unlock();
        return NULL;
    }
    else if (map[group_index][bucket]->size() == 1)
    {
        if (key.second == map[group_index][bucket]->begin()->first.second && key.first == map[group_index][bucket]->begin()->first.first)
        {
            table_lock[group_index][bucket / COVER_AREA].unlock();
            return map[group_index][bucket]->begin()->second;
        }
        else
        {
            table_lock[group_index][bucket / COVER_AREA].unlock();
            return NULL;
        }
    }
    if (is_data_available[group_index])
        for (auto iter = map[group_index][bucket]->begin(); iter != map[group_index][bucket]->end(); iter++)
        {
            if (key.second == iter->first.second && key.first == iter->first.first)
            {
                table_lock[group_index][bucket / COVER_AREA].unlock();
                return iter->second;
            }
        }
    table_lock[group_index][bucket / COVER_AREA].unlock();

    return NULL;
}

bool History_Table::check_and_manage_memory(int depth, float *updated_mem_limit, bool *is_all_table_blocked)
{
    int group_index;
    if (depth <= num_of_groups * groups_size)
        group_index = std::ceil(static_cast<double>(depth) / groups_size) - 1;
    else
        group_index = num_of_groups - 1;

    for (int i = memory_allocators.size() - 1; i > 0; --i)
    {
        // cout << "attempted to block the group :" << group_index << "\n";
        if (current_size < *updated_mem_limit * max_size)
            return blocked_groups[group_index];

        if (blocked_groups[i]) // to prevent using unnecessary locks
            continue;

        group_locks[i].lock();
        if (current_size < *updated_mem_limit * max_size)
        {
            group_locks[i].unlock();
            return blocked_groups[group_index];
            // return group_index > i; // if the subtable number is lesser than i, we CAN insert that record by returning FALSE
        }
        if (!blocked_groups[i])
        {
            cout << "Blocking insertion in bucket: " << i + 1 << " when the mem limit is: " << *updated_mem_limit << std::endl;
            blocked_groups[i] = true;
            *is_all_table_blocked = (i == 1);
            *updated_mem_limit += 0.1f;
            cout << "Updated memory limit: " << *updated_mem_limit << std::endl;

            for (int j = memory_allocators.size() - 1; j >= 0; --j)
            {
                cout << blocked_groups[j] << " ";
            }
            cout << "\n";
            group_locks[i].unlock();
            return blocked_groups[i]; // if the subtable number is equal to the blocking group, we CAN"T insert that record by returning TRUE
        }
        group_locks[i].unlock();
    }
    return blocked_groups[group_index];
    ;
}

bool History_Table::free_subtable_memory(float *mem_limit)
{
    for (int i = memory_allocators.size() - 1; i > 0; --i)
    {
        if (current_size < *mem_limit * max_size)
            return true;

        if (blocked_groups[i] && is_data_available[i]) // Ensure the group is blocked and data is there
        {
            // cout << "locking the group\n";
            group_locks[i].lock();
            if (current_size < *mem_limit * max_size)
            {
                group_locks[i].unlock();
                return true;
            }

            // // Free the memory if group entry is blocked and data is available
            if (blocked_groups[i] && is_data_available[i])
            {
                is_data_available[i] = false;
                cout << "starting the deallocation of memory\n";
                current_size = total_ram - get_free_mem();
                std::cout << "current used size before (in bytes): " << current_size << std::endl;
                std::cout << "size of subtable " << i + 1 << " (in bytes): " << sizeof(memory_allocators[i]) << std::endl;

                for (auto &allocator : memory_allocators[i])
                {
                }
                memory_allocators[i].clear(); // Clear the vector, effectively freeing memory

                std::cout << "Freed memory for subtable: " << i + 1 << std::endl;

                current_size = total_ram - get_free_mem();
                std::cout << "current used size after freeing space (in bytes): " << current_size << std::endl;
                for (int k = 0; k < 3; k++)
                {
                    cout << block_count[k] << " ";
                }
                cout << "\n";
                group_locks[i].unlock();
                return true;
            }
            group_locks[i].unlock();
        }
    }
    return false;
}
