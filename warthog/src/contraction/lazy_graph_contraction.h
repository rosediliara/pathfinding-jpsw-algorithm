#ifndef WARTHOG_LAZY_GRAPH_CONTRACTION_H
#define WARTHOG_LAZY_GRAPH_CONTRACTION_H

// lazy_graph_contraction.h
//
// Applies a contraction operation to each node in a 
// graph. 
//
// Contraction in this case refers to a localised graph
// operation where additional ``shortcut'' edges are added to 
// the graph in such a way that the ``contracted'' node can
// be bypassed entirely during search (i.e. its neighbours are
// connected directly).
//
// In a contraction hierarchy every node is assigned a specific
// ``level'' which is equal to the contraction priority of the
// node. The highest priority nodes (== lowest level) are 
// contracted first and the process continues until every node
// has been contracted.
//
// The node ordering is done in a lazy manner
// using a variety of heuristics.
//
// For more details see:
// [Geisbergerger, Sanders, Schultes and Delling. 
// Contraction Hierarchies: Faster and Simpler Hierarchical 
// Routing in Road Networks. In Proceedings of the 2008
// Workshop on Experimental Algorithms (WEA)]
//
// @author: dharabor
// @created: 2016-01-25
//

#include "ch_data.h"
#include "constants.h"
#include "contraction.h"
#include "flexible_astar.h"
#include "heap.h"
#include "solution.h"
#include "bidirectional_graph_expansion_policy.h"
#include "bidirectional_search.h"

#include <cstdint>
#include <unordered_map>

//using namespace warthog::ch;

namespace warthog
{

namespace ch
{

class ch_pair
{
    public:
        ch_pair() 
            : node_id_(0), cval_(INT32_MAX) { } 

        ch_pair(uint32_t node_id, int32_t cval)
            : node_id_(node_id), cval_(cval) { }

        ch_pair(const ch_pair& other)
        { node_id_ = other.node_id_; cval_ = other.cval_; }

        ~ch_pair() { } 

        uint32_t node_id_; 
        int32_t cval_; // 'value' of contracting this node (lower is better)
};

bool
operator<(const ch_pair& first, const ch_pair& second);

class lazy_graph_contraction 
{
    public:
        lazy_graph_contraction();

        // contract the graph lazily
        //
        // @param c_pct: limit contraction to k% of nodes with 
        // highest priority
        // @param verify_priorities: when true, the priority of each top
        // candidate node is recomputed. if the new value is not better than
        // the second best node on the contraction queue the node is requeued. 
        // process repeats until the top node priority is verified.
        // the hope is that this process yields a better ordering
        bool
        contract(warthog::ch::ch_data* ret, bool verify_priorities=false, uint32_t c_pct=100);

        void
        set_verbose(bool verbose) { this->verbose_ = verbose; }

        bool 
        get_verbose() { return verbose_; } 

        size_t
        mem();

    private:
        warthog::graph::xy_graph* g_;
        std::vector<uint32_t>* order_;
        bool verbose_;
        std::vector<warthog::graph::edge> shortcuts_;

        // node order stuff
        warthog::heap<ch_pair>* heap_;
        warthog::heap_node<ch_pair>* hn_pool_;

        // these objects get recycled across all witness searches
        warthog::solution sol_;

        // track the height of each node in the hierarchy 
        std::vector<uint32_t> height_;

        uint32_t ws_max_expansions_; 
        warthog::zero_heuristic* heuristic_;
        warthog::bidirectional_expander<warthog::ch::bypass_filter>* fexpander_;
        warthog::bidirectional_expander<warthog::ch::bypass_filter>* bexpander_;
        warthog::ch::bypass_filter* c_filter_; // track contractions 
        warthog::apriori_filter* u_filter_; // track neighbours updated
        warthog::bidirectional_search<
            warthog::zero_heuristic,
            warthog::bidirectional_expander<warthog::ch::bypass_filter>>* alg_;

        // metrics
        uint64_t total_expansions_;
        uint64_t total_searches_;
        uint64_t total_lazy_updates_;

        void
        preliminaries(warthog::ch::ch_data*);

        void
        postliminaries();

        uint32_t
        next(bool verify_priorities, uint32_t c_pct);

        double
        witness_search(uint32_t from_id, uint32_t to_id, double via_len, bool resume);

        int32_t
        contract_node(uint32_t node_id, bool metrics_only);

        void
        disconnect_node(uint32_t node_id);
};

}

}

#endif

