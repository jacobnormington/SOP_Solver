#ifndef ACTIVETREE_H
#define ACTIVETREE_H

#include <atomic>
#include <vector>
#include <queue>
#include <deque>
#include "history.hpp"
#include "active_allocate.hpp"

using namespace std;

class Active_Path {
    private:
        vector<Active_Node*> Path;
        int thread_id;
    public:
        Active_Path();
        Active_Path(int arr_size);
        Active_Node* get_element(int index);
        Active_Node* get_node(int thread_id);
        Active_Node* back();
        HistoryNode* back_history();
        int upward_propagation(Active_Allocator& Allocator);
        size_t get_size();
        bool push_back(int children_num, HistoryNode* his_node, Active_Allocator& Allocator);
        bool incre_children_cnt(Active_Allocator& Allocator);
        void collect_deprecated_nodes(Active_Allocator& Allocator, int depth);
        void generate_path(Active_Path& partial_path);
        bool assign_hislink(HistoryNode* link);
        bool check_deprecation_status(int level);
        bool pop_back(bool out_dated, Active_Allocator& Allocator);
        void set_threadID(int id, int total_id);
};

#endif