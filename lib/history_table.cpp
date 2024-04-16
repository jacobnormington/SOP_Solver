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

    table_lock = vector<spin_lock>(size / COVER_AREA + 1);
    map.resize(size);

    insert_count = 0;
}

void History_Table::initialize(int thread_num)
{
    for (int i = 0; i < thread_num; i++)
    {
        memory_allocators.push_back(Memory_Module());
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

// void History_Table::adjust_max_size(int node_count,int GPQ_size) {
//     max_size -= (GPQ_size * (12*node_count + 2*node_count*node_count) * 4);
// }

/* void History_Table::analyze_table() {
    vector<long long unsigned> total_entrynum = vector<long long unsigned>(max_depth,0);
    vector<long long unsigned> total_uses = vector<long long unsigned>(max_depth,0);
    vector<long long unsigned> maximum_uses = vector<long long unsigned>(max_depth,0);
    vector<long long unsigned> minimum_uses = vector<long long unsigned>(max_depth,INT_MAX);

    for (unsigned i = 0; i < size; i++) {
        if (map[i] != NULL) {
            auto iter = map[i]->begin();
            while (iter != map[i]->end()) {
                total_entrynum[iter->second->depth]++;
                total_uses[iter->second->depth] += iter->second->usage_cnt;
                if ((unsigned)iter->second->usage_cnt > maximum_uses[iter->second->depth]) maximum_uses[iter->second->depth] = (unsigned)iter->second->usage_cnt;
                if ((unsigned)iter->second->usage_cnt < minimum_uses[iter->second->depth]) minimum_uses[iter->second->depth] = (unsigned)iter->second->usage_cnt;
                ++iter;
            }
        }
    }

    int lv_delta = ceil(float(max_depth) / 5);
    unsigned starting_lv = 0;
    unsigned ending_lv = lv_delta;

    long long unsigned x_entry = 0;
    long long unsigned x_uses = 0;
    long long unsigned x_maximum = 0;
    long long unsigned x_minimum = 0;
    for (unsigned i = 0; i < 5; i++) {
        for (unsigned j = starting_lv; j < ending_lv; j++) {
            x_entry += total_entrynum[j];
            x_uses += total_uses[j];
            if (maximum_uses[j] > x_maximum) x_maximum = maximum_uses[j];
            if (minimum_uses[j] < x_minimum) x_minimum = minimum_uses[j];
        }
        cout << "---------------------" << i * 20 + 20 << "% [" << starting_lv << "->" << ending_lv <<"]-----------------------" << endl;
        cout << "Total_Entry: [" << x_entry << "]" << endl;
        cout << "Total_Uses: [" << x_uses << "]" << endl;
        cout << "Uses Per Entry: [" << (double)x_uses/(double)x_entry << "]" << endl;
        cout << "Maximum Uses: [" << x_maximum << "]" << endl;
        cout << "Minimum Uses: [" << x_minimum << "]" << endl;
        if (i == 4) cout << "------------------------------------------------" << endl;
        starting_lv = ending_lv;
        ending_lv = ending_lv + lv_delta;
        if (ending_lv >= max_depth) ending_lv = max_depth;
        x_entry = 0;
        x_uses = 0;
        x_maximum = 0;
        x_minimum = 0;
    }
    return;
}
*/

HistoryNode *History_Table::insert(Key &key, int prefix_cost, int lower_bound, unsigned thread_id, bool backtracked, unsigned depth)
{
    HistoryNode *node = memory_allocators[thread_id].retrieve_his_node();

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

    table_lock[bucket / COVER_AREA].lock();
    if (map[bucket] == NULL)
    {
        map[bucket] = memory_allocators[thread_id].get_bucket();
        map[bucket]->push_back(make_pair(key, node));
        table_lock[bucket / COVER_AREA].unlock();
        return node;
    }

    map[bucket]->push_back(make_pair(key, node));
    table_lock[bucket / COVER_AREA].unlock();
    return node;
}

HistoryNode *History_Table::retrieve(Key &key)
{
    size_t val = hash<boost::dynamic_bitset<>>{}(key.first);
    int bucket = (val + key.second) % num_buckets;

    table_lock[bucket / COVER_AREA].lock();
    if (map[bucket] == NULL)
    {
        table_lock[bucket / COVER_AREA].unlock();
        return NULL;
    }
    else if (map[bucket]->size() == 1)
    {
        if (key.second == map[bucket]->begin()->first.second && key.first == map[bucket]->begin()->first.first)
        {
            table_lock[bucket / COVER_AREA].unlock();
            return map[bucket]->begin()->second;
        }
        else
        {
            table_lock[bucket / COVER_AREA].unlock();
            return NULL;
        }
    }

    for (auto iter = map[bucket]->begin(); iter != map[bucket]->end(); iter++)
    {
        if (key.second == iter->first.second && key.first == iter->first.first)
        {
            table_lock[bucket / COVER_AREA].unlock();
            return iter->second;
        }
    }
    table_lock[bucket / COVER_AREA].unlock();

    return NULL;
}
