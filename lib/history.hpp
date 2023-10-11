#ifndef HISTORY_H
#define HISTORY_H

#include <iostream>
#include <vector>
#include <atomic>
#include <cstdint>
#include <mutex>

using namespace std;

class HistoryNode;

struct HistoryContent {
	int32_t prefix_cost = -1;
	int32_t lower_bound = -1;
};

struct Active_Node {
    HistoryNode* history_link = NULL;
	int32_t total_children_cnt = 0; //total number of children of this node
    int32_t cur_children_cnt = 0;	//number of children processed so far
    atomic<int16_t> cur_threadID;
	atomic<bool> deprecated;
	mutex nlck;
};

struct HistoryNode {
	atomic<bool> explored; //if the subspace under this node has already been fully explored
	atomic<uint8_t> active_threadID; //the thread that is exploring that subspace
	atomic<HistoryContent> Entry;
};

#endif