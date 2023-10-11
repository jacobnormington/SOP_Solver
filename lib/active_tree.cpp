#include "active_tree.hpp"
#include "solver.hpp"
#include <climits>
#include <vector>
#include <queue>
#include <deque>
#include <mutex>

static const int RECYCLED = -2;
static const int ZOMBIE = -1;

mutex active_map_lock;

Active_Path::Active_Path() {
    Path = vector<Active_Node*>(1,NULL);
}

Active_Path::Active_Path(int arr_size) {
    Path = vector<Active_Node*>(arr_size,NULL);
}

Active_Node* Active_Path::get_element(int index) {
    if (index < int(Path.size()) && index >= 0) return Path[index];
    return NULL;
}

Active_Node* Active_Path::back() {
    if (!Path.empty()) return Path.back();
    return NULL;
}

void Active_Path::set_threadID(int id, int total_id) {
    thread_id = id;
}

bool Active_Path::push_back(int children_num, HistoryNode* his_node, Active_Allocator& Allocator) {

    Active_Node* n = Allocator.assign_node();
    n->total_children_cnt = children_num;
    n->cur_children_cnt = 0;
    n->cur_threadID = thread_id;

    if (his_node != NULL && his_node->active_threadID == thread_id) {
        n->history_link = his_node;
        his_node->explored = false;
    }
    else n->history_link = NULL;
    Path.push_back(n);

    return true;
}

bool Active_Path::pop_back(bool out_dated, Active_Allocator& Allocator) {
    if (Path.back() == NULL) {
        cout << "pop_back error" << endl;
        cout << "size is " << Path.size() << endl;
        exit(1);
    }

    Path.back()->nlck.lock();
    if (Path.back()->total_children_cnt == Path.back()->cur_children_cnt) {
        if (Path.size() > 1 && Path.end()[-2] != NULL) {
            Path.end()[-2]->nlck.lock();
            Path.end()[-2]->cur_children_cnt++;
            if ((Path.end()[-2]->total_children_cnt == Path.end()[-2]->cur_children_cnt) && Path.end()[-2]->cur_threadID == ZOMBIE) {
                //if (Path.end()[-2]->cur_threadID == RECYCLED) cout << "loc 1" << endl; 
                Path.end()[-2]->nlck.unlock();
                upward_propagation(Allocator);
            }
            else Path.end()[-2]->nlck.unlock();
        }
        if (Path.back()->history_link != NULL) {
            Path.back()->history_link->explored = true;
        }
        Path.back()->deprecated = false;
        Path.back()->cur_threadID = RECYCLED;
        Allocator.delete_node(Path.back());
    }
    else {
        Path.back()->cur_threadID = ZOMBIE;
        Path.back()->deprecated = out_dated;
    }
    Path.back()->nlck.unlock();
    
    Path.pop_back();

    return true;
}

bool Active_Path::check_deprecation_status(int level) {
    if (level >= int(Path.size()) || level < 0) return false;
    else if (Path[level] == NULL) return false;

    if (Path[level]->deprecated) return true;

    return false;
}

bool Active_Path::incre_children_cnt(Active_Allocator& Allocator) {
    if (Path.back() == NULL) {
        cout << "children incre error" << endl;
        cout << "size is " << Path.size() << endl;
        exit(1);
    }

    Path.back()->nlck.lock();
    /*
    if (Path.back()->total_children_cnt == Path.back()->cur_children_cnt) {
        cout << "children incre path error" << endl;
        cout << "size is " << Path.size() << endl;
        exit(1);
    }
    */
    Path.back()->cur_children_cnt++;
    if (Path.back()->total_children_cnt == Path.back()->cur_children_cnt && Path.back()->cur_threadID == ZOMBIE) {
        //if (Path.back()->cur_threadID == RECYCLED) cout << "loc 2" << endl; 
        if (Path.size() > 1 && Path.end()[-2] != NULL) {
            Path.end()[-2]->nlck.lock();
            Path.end()[-2]->cur_children_cnt++;
            if (Path.end()[-2]->total_children_cnt == Path.end()[-2]->cur_children_cnt && Path.end()[-2]->cur_threadID == ZOMBIE) {
                //if (Path.end()[-2]->cur_threadID == RECYCLED) cout << "loc 3" << endl; 
                Path.end()[-2]->nlck.unlock();
                upward_propagation(Allocator);
            }
            else Path.end()[-2]->nlck.unlock();
        }
        
        if (Path.back()->history_link != NULL) {
            Path.back()->history_link->explored = true;
        }
        Path.back()->cur_threadID = RECYCLED;
        Allocator.delete_node(Path.back());
        Path.back()->nlck.unlock();
        Path.pop_back();
    }
    else Path.back()->nlck.unlock();

    return false;
}

void Active_Path::generate_path(Active_Path& partial_path) {
    Path = partial_path.Path;
}

bool Active_Path::assign_hislink(HistoryNode* link) {
    if (Path.back() != NULL) Path.back()->history_link = link;
    return false;
}

int Active_Path::upward_propagation(Active_Allocator& Allocator) {
    //cout << "upward_propagation initialized" << endl;
    
    /*
    for (unsigned i = 0; i < Path.size() - 2; i++) {
        if (Path[i] != NULL) {
            Path[i]->nlck.lock();
            if (Path[i]->total_children_cnt == Path[i]->cur_children_cnt) {
                cout << "Collect error!!!!!" << endl;
                cout << "Total children count is " << Path[i]->total_children_cnt << ", and children count is " << Path[i]->cur_children_cnt << endl;
                cout << "Thread num is " << Path[i]->cur_threadID << endl;
                exit(1);
            }        
            Path[i]->nlck.unlock();
        }
    }
    */

    for (int i = Path.size() - 2; i > 0; i--) {
        if (Path[i] == NULL) break;

        Path[i]->nlck.lock();
        if (Path[i]->total_children_cnt == Path[i]->cur_children_cnt && Path[i]->cur_threadID == ZOMBIE) {
            //if (Path[i]->cur_threadID == RECYCLED) cout << "loc 4" << endl; 
            if (Path[i-1] != NULL) {
                Path[i-1]->nlck.lock();
                Path[i-1]->cur_children_cnt++;
                Path[i-1]->nlck.unlock();
            }
            else {
                Path[i]->nlck.unlock();
                return i;
            }
    
            if (Path[i]->history_link != NULL) {
                Path[i]->history_link->explored = true;
                Path[i]->history_link = NULL;
            }
            Path[i]->cur_threadID = RECYCLED;
            Allocator.delete_node(Path[i]);
            Path[i]->nlck.unlock();
            //Path[i] = NULL;
        }
        else {
            Path[i]->nlck.unlock();
            return i;
        }
    }

    return 1;
}