#ifndef WARTHOG_BIDIRECTIONAL_GRAPH_EXPANSION_POLICY_H
#define WARTHOG_BIDIRECTIONAL_GRAPH_EXPANSION_POLICY_H

// search/bidirectional_graph_expansion_policy.h
//
// An expansion policy for bidirectional search algorithms.
// When expanding in the forward direction this policy
// generates all outgoing neighbours as successors.
// When expanding in the backward direction this policy
// generates all incoming neighbours as successors.
//
// @author: dharabor
// @created: 2019-08-31
//

#include "apriori_filter.h"
#include "expansion_policy.h"
#include "problem_instance.h"
#include "search_node.h"
#include "xy_graph.h"

#include <vector>

namespace warthog{

class problem_instance;
class search_node;

template< class FILTER >
class bidirectional_expander //: public  expansion_policy
{
    public:
        // @param g: the input contracted graph
        // @param backward: when true successors are generated by following 
        // incoming arcs rather than outgoing arcs (default is outgoing)
        bidirectional_expander(warthog::graph::xy_graph* g, 
                bool backward, FILTER* filter = 0)
        {
            g_ = g;
            backward_ = backward;
            filter_ = filter;

            if(backward_)
            {
                fn_begin_iter_ = 
                    &warthog::bidirectional_expander<FILTER>::get_bwd_begin_iter;
                fn_end_iter_ = 
                    &warthog::bidirectional_expander<FILTER>::get_bwd_end_iter;

                fn_rev_end_iter_ = 
                    &warthog::bidirectional_expander<FILTER>::get_fwd_end_iter;
                fn_rev_begin_iter_ = 
                    &warthog::bidirectional_expander<FILTER>::get_fwd_begin_iter;
            }
            else
            {
                fn_begin_iter_ = 
                    &warthog::bidirectional_expander<FILTER>::get_fwd_begin_iter;
                fn_end_iter_ = 
                    &warthog::bidirectional_expander<FILTER>::get_fwd_end_iter;

                fn_rev_begin_iter_ = 
                    &warthog::bidirectional_expander<FILTER>::get_bwd_begin_iter;
                fn_rev_end_iter_ = 
                    &warthog::bidirectional_expander<FILTER>::get_bwd_end_iter;
            }

            nodepool_.resize(g_->get_num_nodes());
            for(uint32_t i = 0; i < g_->get_num_nodes(); i++)
            {
                nodepool_[i].set_id(i);
            }
        }

        ~bidirectional_expander()
        { }

        inline void 
        expand(warthog::search_node* current, warthog::problem_instance* pi)
        {
            warthog::graph::node* n = g_->get_node(current->get_id());
            begin_ = (this->*fn_begin_iter_)(n);
            end_ = (this->*fn_end_iter_)(n);
        }

        inline void
        first(warthog::search_node*& ret, double& cost)
        {
            for(it_ = begin_; it_ != end_; it_++)
            {
            assert((*it_).node_id_ < g_->get_num_nodes());
            if(!filter_ || !filter_->filter((*it_).node_id_))
            {
                ret = generate((*it_).node_id_); 
                cost = (*it_).wt_;
                it_++;
                return;
            }
            }
            ret = 0;
            cost = warthog::COST_MAX;
        }

        inline void
        next(warthog::search_node*& ret, double& cost)
        {
            for( ; it_ != end_; it_++)
            {
            assert((*it_).node_id_ < g_->get_num_nodes());
            if(!filter_ || !filter_->filter((*it_).node_id_))
            {
                ret = generate((*it_).node_id_); 
                cost = (*it_).wt_;
                it_++;
                return;
            }
            }
            ret = 0;
            cost = warthog::COST_MAX;
        }	

        inline warthog::search_node* 
        generate_start_node(warthog::problem_instance* pi)
        {
            uint32_t s_graph_id = g_->to_graph_id((uint32_t)pi->start_id_);
            if(s_graph_id == warthog::INF32) { return 0; }
            return generate(s_graph_id);
        }

            inline warthog::search_node*
            generate_target_node(warthog::problem_instance* pi)
        {
            uint32_t t_graph_id = g_->to_graph_id((uint32_t)pi->target_id_);
            if(t_graph_id == warthog::INF32) { return 0; }
            return generate(t_graph_id);
        }

        // check if a given search node @param n corresponds to the
        // target. we do this here to decouple the internal 
        // representation of states from the search algorithm which
        // only knows about warthog::search_node objects.
        bool
        is_target(warthog::search_node* n, warthog::problem_instance* pi)
        {
            return n->get_id() == pi->target_id_;
        }

        inline warthog::search_node*
        generate(warthog::sn_id_t nid)
        { return &nodepool_[nid]; }

        inline void
        get_xy(warthog::sn_id_t node_id, int32_t& x, int32_t& y)
        { g_->get_xy((uint32_t)node_id, x, y); }

        inline size_t
        get_num_nodes() 
        { return g_->get_num_nodes(); }

        inline size_t
        mem()
        {
            return 
            sizeof(warthog::search_node)*nodepool_.size() + 
            sizeof(this);
        }

    private:
        bool backward_;
        warthog::graph::xy_graph* g_;

        warthog::graph::edge_iter begin_, end_, it_;
        std::vector<warthog::search_node> nodepool_;
        FILTER* filter_;

        // we use function pointers to iterate over the right set of
        // neighbours (incoming or outgoing) 
        typedef warthog::graph::edge_iter
                (warthog::bidirectional_expander<FILTER>::*chep_get_iter_fn) 
                (warthog::graph::node* n);
        
        // pointers to the neighbours in the direction of the search
        chep_get_iter_fn fn_begin_iter_;
        chep_get_iter_fn fn_end_iter_;

        // pointers to neighbours in the reverse direction to the search
        chep_get_iter_fn fn_rev_begin_iter_;
        chep_get_iter_fn fn_rev_end_iter_;

        inline warthog::graph::edge_iter
        get_fwd_begin_iter(warthog::graph::node* n) 
        { return n->outgoing_begin(); }

        inline warthog::graph::edge_iter
        get_fwd_end_iter(warthog::graph::node* n) 
        { return n->outgoing_end(); }

        inline warthog::graph::edge_iter
        get_bwd_begin_iter(warthog::graph::node* n) 
        { return n->incoming_begin(); }

        inline warthog::graph::edge_iter
        get_bwd_end_iter(warthog::graph::node* n) 
        { return n->incoming_end(); }
};
typedef bidirectional_expander<warthog::apriori_filter> bidirectional_graph_expansion_policy;


}
#endif

