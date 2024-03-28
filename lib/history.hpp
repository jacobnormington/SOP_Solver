#ifndef HISTORY_H
#define HISTORY_H

#include <atomic>
#include <cstdint>
#include <mutex>

using namespace std;

struct HistoryContent
{
	int32_t prefix_cost = -1; // the cost of the current path represented by this node
	int32_t lower_bound = -1; // the lower bound cost of a solution starting with this path
};

/* A single node in the history table. */
struct HistoryNode
{
	atomic<bool> explored;			 // if the subspace under this node has already been fully explored
	atomic<uint8_t> active_threadID; // the thread that is exploring this subspace (there can only ever be one, because any others would be stopped)
	atomic<HistoryContent> entry;
};

struct Active_Node
{
	HistoryNode *history_link = NULL;
	int32_t total_children_cnt = 0; // total number of children of this node
	int32_t cur_children_cnt = 0;	// number of children processed so far
	atomic<int16_t> cur_threadID;
	atomic<bool> deprecated; // whether this node is still necessary to consider
	mutex nlck;				 // lock on total and current children counts
};

#endif