#include <vector>
#include <unordered_map>
#include <iostream>
#include <list>
#include <mutex>
#include <cmath>
#include <climits>
#include <cerrno>
#include <atomic>
#include <sys/sysinfo.h>
#include <malloc.h>
#include "history.hpp"
#include "hash.hpp"

unsigned long starting_max_size = 0;
#define MEMORY_RESTRIC 0.95
#define SAMPLE_FREQUENCY 60

Mem_Allocator::Mem_Allocator() {
    Bucket_blk = new Bucket[BUCKET_BLK_SIZE];
    history_block = (HistoryNode*)malloc(HIS_BLK_SIZE*sizeof(HistoryNode));
    counter = 0;
    node_counter = 0;
}

Bucket* Mem_Allocator::Get_bucket() {
    if (counter == BUCKET_BLK_SIZE || Bucket_blk == NULL) {
        Bucket_blk = new Bucket[BUCKET_BLK_SIZE];
        counter = 0;
    }
    Bucket* bucket = Bucket_blk + counter;
    counter++;
    return bucket;
}

HistoryNode* Mem_Allocator::retrieve_his_node() {
    HistoryNode* node = NULL;

    if (node_counter >= HIS_BLK_SIZE || history_block == NULL) {
        history_block = (HistoryNode*)malloc(HIS_BLK_SIZE*sizeof(HistoryNode));
        node_counter = 0;
    }
    node = history_block + node_counter;
    node_counter++;

    return node;
}

Hash_Map::Hash_Map(size_t size) {
    struct sysinfo info;
    this->size = size;
	if(sysinfo(&info) != 0){
        cout << "can't retrieve sys mem info\n";
		exit(1);
	}
    starting_max_size = (unsigned long)info.freeram;
    max_size = (double)info.freeram * MEMORY_RESTRIC - (size * 4);
    total_size = info.totalram;
    cur_size = 0;

    //cout << "Max bucket size is " << max_size / 1000000 << " MB" << endl;
    //cout << "Total Available Memory In OS is " << info.totalram / 1000000 << " MB" << endl;

    hash_lock = vector<spinlock>(size/COVER_AREA + 1);
    History_table.resize(size);

    insert_count = 0;
}

void Hash_Map::adjust_max_size(int node_count,int GPQ_size) {
    max_size -= (GPQ_size * (12*node_count + 2*node_count*node_count) * 4);
}

void Hash_Map::set_up_mem(int thread_num,int node_count) {
    for (int i = 0; i < thread_num; i++) {
        Mem_Manager.push_back(Mem_Allocator());
    }
}


size_t Hash_Map::get_max_size() {
    return max_size;
}

size_t Hash_Map::get_cur_size() {
    return cur_size;
}

unsigned long Hash_Map::get_free_mem() {
    struct sysinfo info;
	if(sysinfo(&info) != 0){
        cout << "can't retrieve sys mem info\n";
		exit(1);
	}
    return (double)info.freeram;
}

void Hash_Map::print_curmem() {
    cout << "Mem exhausted is " << cur_size / 1000000 << "MB" << endl;
    return;
}

/*
void Hash_Map::analyze_table() {
    vector<long long unsigned> total_entrynum = vector<long long unsigned>(max_depth,0);
    vector<long long unsigned> total_uses = vector<long long unsigned>(max_depth,0);
    vector<long long unsigned> maximum_uses = vector<long long unsigned>(max_depth,0);
    vector<long long unsigned> minimum_uses = vector<long long unsigned>(max_depth,INT_MAX);

    for (unsigned i = 0; i < size; i++) {
        if (History_table[i] != NULL) {
            auto iter = History_table[i]->begin();
            while (iter != History_table[i]->end()) {
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
    
HistoryNode* Hash_Map::insert(pair<boost::dynamic_bitset<>,int>& item,int prefix_cost,int lower_bound, unsigned thread_id, bool backtracked, unsigned depth) {
    HistoryNode* node = Mem_Manager[thread_id].retrieve_his_node();

    if (thread_id == 0) {
        insert_count++;
        if (insert_count >= 100000) {
            cur_size = total_size - get_free_mem();
            insert_count = 0;
        }
    }

    if (node == NULL) return NULL;

    size_t val = hash<boost::dynamic_bitset<>>{}(item.first);
    int key = (val + item.second) % size;
    
    if (backtracked) node->explored = true;
    else node->explored = false;
    node->Entry.store({prefix_cost,lower_bound});
    node->active_threadID = thread_id;
    //node->depth = depth;
    //node->usage_cnt = 0;

    hash_lock[key/COVER_AREA].lock();
    if (History_table[key] == NULL) {
        History_table[key] = Mem_Manager[thread_id].Get_bucket();
        History_table[key]->push_back(make_pair(item,node));
        hash_lock[key/COVER_AREA].unlock();
        return node;
    }
    
    History_table[key]->push_back(make_pair(item,node));
    hash_lock[key/COVER_AREA].unlock();
    return node;
}

HistoryNode* Hash_Map::retrieve(pair<boost::dynamic_bitset<>,int>& item,int key) {

    hash_lock[key/COVER_AREA].lock();
    if (History_table[key] == NULL) {
        hash_lock[key/COVER_AREA].unlock();
        return NULL;
    }
    else if (History_table[key]->size() == 1) {
        if (item.second == History_table[key]->begin()->first.second && item.first == History_table[key]->begin()->first.first) {
            hash_lock[key/COVER_AREA].unlock();
            return History_table[key]->begin()->second;
        }
        else {
            hash_lock[key/COVER_AREA].unlock();
            return NULL;
        }
    }

    for (auto iter = History_table[key]->begin(); iter != History_table[key]->end(); iter++) {
        if (item.second == iter->first.second && item.first == iter->first.first) {
            hash_lock[key/COVER_AREA].unlock();
            return iter->second;
        }
    }
    hash_lock[key/COVER_AREA].unlock();

    return NULL;
}
