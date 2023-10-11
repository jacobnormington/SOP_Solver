#include "active_allocate.hpp"
#include <climits>
#include <vector>

using namespace std;
#define ACTIVE_BLK_SIZE 81920

Active_Allocator::Active_Allocator() {
    Active_Block = new Active_Node[ACTIVE_BLK_SIZE];
    counter = 0;
}

Active_Node* Active_Allocator::assign_node() {
    Active_Node* node = NULL;
    if (Reuse_Node.empty()) {
        if (counter >= ACTIVE_BLK_SIZE) {
            Active_Block = new Active_Node[ACTIVE_BLK_SIZE];
            counter = 0;
        }
        node = Active_Block + counter;
        counter++;
    }
    else {
        node = Reuse_Node.back();
        Reuse_Node.pop_back();
    }
    node->deprecated = false;
    node->total_children_cnt = -1;
    node->cur_children_cnt = -1;
    node->history_link = NULL;
    return node;
}

void Active_Allocator::delete_node(Active_Node* node) {
    if (node != NULL) Reuse_Node.push_back(node);
    return;
}