#ifndef ACTIVEALLOCATE_H
#define ACTIVEALLOCATE_H

#include <vector>
#include "history.hpp"

using namespace std;

class Active_Allocator {
    private:
        int counter;
        Active_Node* Active_Block;
        vector<Active_Node*> Reuse_Node;
    public:
        Active_Allocator();
        Active_Node* assign_node();
        void delete_node(Active_Node* node);
};

#endif