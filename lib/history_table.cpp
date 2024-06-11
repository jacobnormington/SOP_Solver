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

Bucket *Memory_Module::get_bucket()
{
    if (bucket_counter == BUCKET_BLK_SIZE || bucket_block == NULL)
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

    table_lock.resize(number_of_groups);
    map.resize(number_of_groups);
    for (int i = 0; i < number_of_groups; i++)
    {
        table_lock[i] = vector<spin_lock>(size / COVER_AREA + 1);
        map[i].resize(size);
    }

    memory_allocators.resize(number_of_groups);
    for (int i = 0; i < number_of_groups; i++)
    {
        memory_allocators[i].resize(thread_num);
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

    HistoryNode *node = memory_allocators[group_index][thread_id].retrieve_his_node();

    if (thread_id == 0)
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

    size_t val = hash<boost::dynamic_bitset<>>{}(key.first);
    int bucket = (val + key.second) % num_buckets;

    table_lock[group_index][bucket / COVER_AREA].lock();
    if (map[group_index][bucket] == NULL)
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
