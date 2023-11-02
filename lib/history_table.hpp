#ifndef HASH_H
    #define HASH_H

    #include <iostream>

    #include <sys/sysinfo.h>

    // #include <boost/smart_ptr/detail/spinlock.hpp>
    #include <boost/dynamic_bitset.hpp>

    #include "history.hpp"
    #include "synchronization.hpp"

    #define BUCKET_BLK_SIZE 81920
    #define HIS_BLK_SIZE 81920
    #define FREED_SIZE 100000000
    #define COVER_AREA 10

    typedef pair<boost::dynamic_bitset<>,int> Key;  //a key to the history table
    typedef pair<Key,HistoryNode*> Entry;           //a single entry in the history table
    typedef list<Entry> Bucket;                     //a bucket of

    /* Manages memory allocation for the history table in blocks in order to reduce system call overhead. */
    class Memory_Module { 
        private:
            Bucket* bucket_block = NULL;        //current block of buckets
            unsigned bucket_counter;            //counter in the current block of buckets
            HistoryNode* history_block = NULL;  //current block of history nodes
            unsigned his_node_counter = 0;      //counter in the current block of history nodes
            // Recycled_Blk recycled_nodes;
            // unsigned recycled_node_counter = 0;
            // unsigned max_level = 0;
        public:
            Memory_Module();
            Bucket* get_bucket();
            HistoryNode* retrieve_his_node();
            // unsigned long request_nodeMem();
            // bool check_reserve();
    };

    // struct Recycled_Blk {
    //     HistoryNode** node_block = NULL;
    //     int total_counter = -1;
    // };

    /* A thread-safe collection of history entries for previously processed partial paths. 
        For efficiency, the history table is allocated in buckets of many nodes, not individually. */
    class History_Table {
        private: 
            unsigned num_buckets = 0;               //the number of buckets the history table should be stored in
            vector<Bucket*> map;                    //the collection of history nodes
            vector<spin_lock> table_lock;           //a read-write lock for every X adjacent buckets, defined by COVER_AREA
            vector<Memory_Module> memory_allocators;   //a memory allocator for each thread, in order to reduce synchronization overhead

            unsigned long total_ram = 0;            //the total amount of memory in the system, in bytes
            unsigned long max_size = 0;             //the maximum allowed size of the history table, in bytes
            atomic<unsigned long> current_size;     //the current size of the history table, in bytes
            int insert_count;                       //a counter to ensure that, periodically, current_size is updated
        public:
            /* Create a new history table object.
                size - the number of buckets the history table should be stored in */
            History_Table(size_t size);
            /* Initialize the history table memory module (in parallel solver sets up one for each thread). 
                thread_count -  the number of threads alloted to B&B enumeration 
                node_count - the number of nodes in the cost graph */
            void initialize(int thread_count);
            /* Returns the max allowed size of the history table, in bytes. */
            size_t get_max_size();
            /* Returns the current size of the history table, in bytes. */
            size_t get_current_size();
            /* Calculates the amount of free ram, in bytes, available on system. */
            unsigned long get_free_mem();
            /* Print to console the total amount of memory exhausted so far. */
            void print_curmem();
            // void adjust_max_size(int node_count, int GPQ_size);
            //void analyze_table();

            /* Add a new entry to the history table. 
                key - the history key corresponding to the partial path this entry represents
                prefix_cost - the cost of the partial path this entry represents
                lower_bound - the lower bound cost of a complete solution beginning with this path
                thread_id - the thread searching this subspace
                backtracked - whether the subspace under this node has already been explored
                depth - the depth of this node (size of the current partial path)
                Return- a pointer to the node created */
            HistoryNode* insert(Key& key, int prefix_cost,int lower_bound, unsigned thread_id, bool backtracked, unsigned depth);
            /* Find a history table entry based on a specific key. 
                key - the history key corresponding to the partial path this entry represents
                Return- a pointer to the node found, if any */
            HistoryNode* retrieve(Key& key);
    };

#endif