#ifndef HASH_H
#define HASH_H

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
#include <pthread.h>
#include <boost/smart_ptr/detail/spinlock.hpp>
#include <boost/dynamic_bitset.hpp>
#include "history.hpp"

#define BUCKET_BLK_SIZE 81920
#define HIS_BLK_SIZE 81920
#define FREED_SIZE 100000000
#define COVER_AREA 10

typedef pair<pair<boost::dynamic_bitset<>,int>,HistoryNode*> Entry;
typedef list<Entry> Bucket;

struct Recycled_Blk {
    HistoryNode** Node_Blk = NULL;
    int total_counter = -1;
};
 
class spinlock{
    public:
        void lock() {
            while(lck.test_and_set(std::memory_order_acquire)){}
        };
        void unlock() {
            lck.clear(std::memory_order_release);
        };
    private:
        std::atomic_flag lck = ATOMIC_FLAG_INIT;
};

class Mem_Allocator {
    private:
        Bucket* Bucket_blk = NULL;
        unsigned counter;
        HistoryNode* history_block = NULL;
        Recycled_Blk recycled_nodes;
        unsigned node_counter = 0;
        unsigned recycle_nodecnter = 0;
        unsigned max_level = 0;
    public:
        Mem_Allocator();
        Bucket* Get_bucket();
        HistoryNode* retrieve_his_node();
        unsigned long Request_NodeMem();
        bool Check_Reserve();
};

class Hash_Map {
    private:
        atomic<unsigned long> cur_size;
        int insert_count;
        unsigned long total_size = 0;
        unsigned long max_size = 0;
        unsigned size = 0;
        unsigned thread_cnt = 0;
        size_t node_size = 0;
        vector<Bucket*> History_table;
        vector<spinlock> hash_lock;
        vector<Mem_Allocator> Mem_Manager;
    public:
        Hash_Map(size_t size);
        size_t get_max_size();
        size_t get_cur_size();
        unsigned long get_free_mem();
        uint32_t hash_func(pair<vector<bool>,int> item);
        void adjust_max_size(int node_count, int GPQ_size);
        void set_up_mem(int thread_num,int node_count);
        void analyze_table();
        void print_curmem();
        HistoryNode* insert(pair<boost::dynamic_bitset<>,int>& item,int prefix_cost,int lower_bound, unsigned thread_id, bool backtracked, unsigned depth);
        HistoryNode* retrieve(pair<boost::dynamic_bitset<>,int>& item,int key);
};

#endif