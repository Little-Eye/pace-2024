////////////////////////////////
/// usage : 1.	general trie implementation.
/// 
/// note  : 1.	
////////////////////////////////

#ifndef CN_HUST_GOAL_COMMON_PREFIX_TREE
#define CN_HUST_GOAL_COMMON_PREFIX_TREE


#include "GOAL/Typedef.h"


namespace goal {

// it is recommended to make the `Data` a pointer type.
template<typename Ch, typename Str, typename Data>
struct PrefixTree {
    // return the leaf node for `str`, set the `data` of each new node along the path.
    PrefixTree<Ch, Str, Data>* add(const Str &str, const Data &defaultData) {
        PrefixTree<Ch, Str, Data> *p = this;
        for (auto c = str.begin(); c != str.end(); ++c) {
            auto n = p->tree.find(*c);
            if (n == p->tree.end()) {
                p = &(p->tree[*c]);
                p->data = defaultData;
            } else {
                p = &(n->second);
            }
        }
        return p;
    }

    // return NULL if `str` does not exist.
    PrefixTree<Ch, Str, Data>* find(const Str &str) {
        PrefixTree<Ch, Str, Data> *p = this;
        for (auto c = str.begin(); c != str.end(); ++c) {
            auto n = p->tree.find(*c);
            if (n == p->tree.end()) { return nullptr; }
            p = &(n->second);
        }
        return p;
    }

    // `onNode` can be `std::function<bool(Data&)>` and return true if the sub-tree is cut.
    template<typename OnNode>
    void traverse(OnNode onNode) {
        if (onNode(data)) { return; }
        for (auto n = tree.begin(); n != tree.end(); ++n) { n->second.traverse(onNode); }
    }

    Data data;

protected:
    Map<Ch, PrefixTree<Ch, Str, Data>> tree;
};

}


#endif // CN_HUST_GOAL_COMMON_PREFIX_TREE
