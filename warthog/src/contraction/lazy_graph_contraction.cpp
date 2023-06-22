#include "apriori_filter.h"
#include "euclidean_heuristic.h"
#include "flexible_astar.h"
#include "graph_expansion_policy.h"
#include "heap.h"
#include "helpers.h"
#include "lazy_graph_contraction.h"
#include "zero_heuristic.h"

#include <algorithm>
#include <functional>
#include <utility>

struct shortcut_t
{
    uint32_t from_;
    uint32_t to_;
    uint32_t mid_;
    warthog::graph::edge_cost_t cost_;
};
std::vector<shortcut_t> shortcuts;

bool very_verbose = false;
bool debug=false;

warthog::ch::lazy_graph_contraction::lazy_graph_contraction()
    : verbose_(false) 
{ }

bool
warthog::ch::operator<(const ch_pair& first, const ch_pair& second)
{
    return first.cval_ < second.cval_;
}

bool
warthog::ch::lazy_graph_contraction::contract(
        warthog::ch::ch_data* ret,
        bool verify_priorities, uint32_t c_pct)
{
    warthog::timer mytimer;
    double t_begin = mytimer.get_time_micro();
    std::cerr << "contracting graph " << ret->g_->get_filename() << std::endl;
    std::ifstream ifs(ret->g_->get_filename());
    ret->type_ = warthog::ch::ch_type::UP_DOWN;

    preliminaries(ret);
    if(c_pct < 100)
    { std::cerr << "partial contraction " << "(first "<<c_pct<<"% of nodes only)\n"; }
    uint32_t edges_before = ret->g_->get_num_edges_out();

    // contract nodes in the order produced by ::next
    //uint32_t total_nodes = (uint32_t)((ret->g_->get_num_nodes()*c_pct) / 100);
    double t_last = mytimer.get_time_micro();

    //verbose_ = true;
    //very_verbose = true;
    while(true)
    {
        // contract the highest priority node  +
        // record some search effort metrics
        uint64_t num_expansions = total_expansions_;
        uint64_t num_searches = total_searches_;
        uint64_t num_lazy = total_lazy_updates_;
        mytimer.start();

        uint32_t best_id = next(verify_priorities, c_pct);
        if(best_id == warthog::INF32) { break; }


        uint32_t search_count = total_searches_;

        if(this->verbose_)
        {
            std::cerr 
                << "contracting " << best_id << " pri: " 
                << hn_pool_[best_id].get_element().cval_ 
                //<< " #searches: " << (total_searches_ - search_count)
                << std::endl;
        }

        contract_node(best_id, false);
        disconnect_node(best_id);

        warthog::graph::node* bestnode = g_->get_node(best_id);
        for(uint32_t i = 0; i < bestnode->out_degree() + bestnode->in_degree(); i++)
        {
            uint32_t neighbour_id; 
            if(i < bestnode->in_degree())
            {
                neighbour_id = (bestnode->incoming_begin() + i)->node_id_;
            }
            else
            {
                neighbour_id = (bestnode->outgoing_begin() + (i - bestnode->in_degree()))->node_id_;
            }
            if(!u_filter_->filter(neighbour_id)) { continue; } // already processed
            u_filter_->remove(neighbour_id);

            // re-compute importance metrics
            search_count = total_searches_;
            uint32_t expd = total_expansions_;
            height_[neighbour_id] = std::max(height_[neighbour_id], height_[best_id]+1);
            int32_t cval = contract_node(neighbour_id, true);

            total_lazy_updates_++;

            // update heap 
            if(this->verbose_)
            {
                std::cerr <<  "updating " << neighbour_id << " old_cval " 
                          << hn_pool_[neighbour_id].get_element().cval_ 
                          << " new_cval " << cval 
                          << " expd " << (total_expansions_ - expd)
                          //<< " #searches: " << (total_searches_ - search_count)
                          << std::endl;
            }

            if(cval < hn_pool_[neighbour_id].get_element().cval_)
            {
               hn_pool_[neighbour_id].get_element().cval_ = cval;
               heap_->decrease_key(&hn_pool_[neighbour_id]);
            }
            else
            {
                hn_pool_[neighbour_id].get_element().cval_ = cval;
                heap_->increase_key(&hn_pool_[neighbour_id]);
            }
        }

        mytimer.stop();

        num_expansions = total_expansions_ - num_expansions;
        num_searches = total_searches_ - num_searches;
        num_lazy = total_lazy_updates_ - num_lazy;

        if((mytimer.get_time_micro() - t_last) > 1000000)
        {
            t_last = mytimer.get_time_micro();
            std::cerr
                //<< "\r "
                //<< (int)((order_.size() / (double)g_->get_num_nodes()) * 100)
                << "time " << (int)((t_last - t_begin)/1000000) << "s"
                << " prog " << order_->size()
                << " nid " << best_id 
                << " pri " << hn_pool_[best_id].get_element().cval_
                << " outdeg " << bestnode->out_degree() 
                << " indeg " << bestnode->in_degree() 
                << " #witness " << num_searches
                << " #updates " << num_lazy
                << " #expansions " << num_expansions
                //<< "; micros: " << (int32_t)mytimer.elapsed_time_micro();
                << std::endl;
        }
    }

    // add shortcuts into CH graph
    for(uint32_t i = 0; i < shortcuts.size(); i++)
    {
        shortcut_t tmp = shortcuts.at(i);
        ret->g_->get_node(tmp.from_)->add_outgoing(
                warthog::graph::edge(tmp.to_, tmp.cost_));
        ret->g_->get_node(tmp.to_)->add_incoming(
                warthog::graph::edge(tmp.from_, tmp.cost_));
    }

    // until now we kept track of the _order_ in which nodes were
    // contracted. but for online search we need to know each node's _level_
    // here we convert from one to the other
    warthog::helpers::value_index_swap_array(*ret->level_);

    // sort successors in descending contraction order
    // (the highest-level successor appears first)
    warthog::ch::sort_successors(ret);

    // cleanup
    postliminaries();

    std::cerr
        << "\ngraph, contracted. " << std::endl
        << "time "
        << ((double)(mytimer.get_time_micro() - t_begin)) / 1e6 << " (s)"
        << "; edges before " << edges_before 
        << "; edges after " << ret->g_->get_num_edges_out() << std::endl
        << "total searches " << total_searches_ 
        << "; total expansions (x 1e6) " << ((double)total_expansions_)/1e6 
        << std::endl
        << "lazy updates "<< total_lazy_updates_
        << std::endl;
    return true;
}

void
warthog::ch::lazy_graph_contraction::preliminaries(warthog::ch::ch_data* chd)
{
    //g_ = chd->g_;
    g_ = new warthog::graph::xy_graph(*chd->g_);
    order_ = chd->level_;
    size_t sz_g = g_->get_num_nodes();

    // init algorithm for witness searches
    heuristic_ = new warthog::zero_heuristic();
    c_filter_ = new warthog::ch::bypass_filter();
    u_filter_ = new warthog::apriori_filter(g_->get_num_nodes());

    fexpander_ = new warthog::bidirectional_expander<warthog::ch::bypass_filter>(
                g_, false, c_filter_);
    bexpander_ = new warthog::bidirectional_expander<warthog::ch::bypass_filter>(
                g_, true, c_filter_);
    alg_ = new bidirectional_search< 
                 warthog::zero_heuristic,
                 warthog::bidirectional_expander<warthog::ch::bypass_filter>>
                     (fexpander_, bexpander_, heuristic_);
    ws_max_expansions_ = warthog::INF32;

    // init queue used to track contraction priorities
    heap_ = new warthog::heap<ch_pair>((uint32_t)sz_g, true);
    hn_pool_ = new heap_node<ch_pair>[sz_g];
    

    // initial heights of all the nodes is 0
    height_.resize(sz_g, 0);

    // initialise edge labels with hop numbers (initially, all 1)
    for(uint32_t i = 0; i < g_->get_num_nodes(); i++)
    {
        warthog::graph::node* n = g_->get_node(i);

        std::sort(
                n->outgoing_begin(), n->outgoing_end(), 
                [](warthog::graph::edge& e1, warthog::graph::edge& e2)
                { return e1.node_id_ < e2.node_id_; });

        std::sort(
                n->incoming_begin(), n->incoming_end(), 
                [](warthog::graph::edge& e1, warthog::graph::edge& e2)
                { return e1.node_id_ < e2.node_id_; } );

        for(uint32_t j = 0; j < n->out_degree(); j++)
        {
            warthog::graph::edge& out = *(n->outgoing_begin() + j);
            out.label_ = 1;
        }
        for(uint32_t j = 0; j < n->in_degree(); j++)
        {
            warthog::graph::edge& in = *(n->incoming_begin() + j);
            in.label_ = 1;
        }
    }

    // create an initial ordering by performing on each node
    // a faux contraction operation
    std::cerr << "creating initial contraction order\n";
    total_searches_ = 0;
    total_expansions_ = 0;
    total_lazy_updates_ = 0;

    double pct_counter=0;
    for(uint32_t i = 0;
        i < g_->get_num_nodes(); i++)
    {
        uint32_t exps = total_expansions_;
        int32_t priority = contract_node(i, true);
        hn_pool_[i] = heap_node<ch_pair>(ch_pair(i, priority));
        heap_->push(&hn_pool_[i]);

        if(this->verbose_)
        //if(true)
        {
            std::cerr << "node " << i << " initial_priority " << priority 
                      << " queue_pos " << hn_pool_[i].get_index() 
                      << " exps " << (total_expansions_ - exps) 
                      << std::endl;

        }
        else 
        if( (i / (double)g_->get_num_nodes()) >= pct_counter)
        {
            std::cerr << pct_counter*100 << "%...";
            pct_counter+=0.1;
        }
    }
    std::cerr << "done" << std::endl;
}

void
warthog::ch::lazy_graph_contraction::postliminaries()
{
    delete g_;

    delete [] hn_pool_;
    hn_pool_ = 0;

    delete heap_;
    heap_ = 0;

    delete alg_;
    alg_ = 0;

    delete c_filter_;
    c_filter_ = 0;

    delete u_filter_;
    u_filter_ = 0;
    
    delete heuristic_;
    heuristic_ = 0;

    delete fexpander_;
    fexpander_ = 0;

    delete bexpander_;
    bexpander_ = 0;

}

// identifies the node that should be contracted next NB: nodes are ranked in a
// lazy way; an initial priority is generated (lower is better) but priorities
// are not updated until a node is popped off the heap.  at this stage we
// re-compute the node priority and possibly insert it back into the heap if
// it is no longer the ``best'' candidate
//
// @return the next node id or warthog::INF32 if all nodes have been contracted
uint32_t
warthog::ch::lazy_graph_contraction::next(
        bool verify_priorities, uint32_t c_pct)
{
    // early abort; all nodes contracted
    if(heap_->size() == 0)
    {
        return warthog::INF32;
    }

    // early abort; partial contraction limit reached
    uint32_t pct = (uint32_t)((order_->size() / (warthog::cost_t)g_->get_num_nodes()) * 100);
    if(pct >= c_pct)
    {
        std::cerr << "\npartial contraction finished "
                  << "(processed "<< pct << "% of all nodes)";
        return warthog::INF32;
    }

    // pop the best contraction candidate off the heap
    heap_node<ch_pair>* best = heap_->pop();
    return best->get_element().node_id_;
}

// NB: assumes the via-node is already marked as contracted
// (and will thus not be expanded)
warthog::cost_t
warthog::ch::lazy_graph_contraction::witness_search(
        uint32_t from_id, uint32_t to_id, warthog::cost_t via_len, bool resume)
{
    alg_->set_cost_cutoff(via_len);
    alg_->set_max_expansions_cutoff(ws_max_expansions_);
    alg_->set_time_cutoff(1);
    warthog::graph::xy_graph* g = this->g_;

    // need to specify start + target ids using the identifier
    // that appears in the input file
    warthog::problem_instance pi(
            g->to_external_id(from_id),
            g->to_external_id(to_id));

    sol_.reset();
    sol_.path_.clear();

    if(very_verbose)
    {
        std::cerr << "search: " << from_id << " to " << to_id << std::endl;
    }
    pi.verbose_ = very_verbose;

    // gogogo
    alg_->__get_pathcost(pi, sol_, resume);
    total_expansions_ += sol_.nodes_expanded_;
    total_searches_++;

    return sol_.sum_of_edge_costs_;
}

// contract a node or compute its contraction priority 
// (depending on the value of @param metrics_only)
// @param n: the id of the node to contract
//
// @return the contraction priority (always computed)
int32_t
warthog::ch::lazy_graph_contraction::contract_node( 
        uint32_t node_id, bool metrics_only)
{
    // expansion limit per witness search (as recommended by RoutingKit)
    ws_max_expansions_ = 500; 
    c_filter_->add(node_id);

    // the set of in/out neighbours pairs that we consider as fixed
    // (contracting adds new edges, and we don't iterate over these)
    warthog::graph::node* n = g_->get_node(node_id);
    uint32_t in_deg = n->in_degree();
    uint32_t out_deg = n->out_degree();

    // metrics that track how contracting @param node_id affects the graph
    // (here we follow RoutingKit's implementation and count deletions from 1)
    uint32_t edges_added=0, edges_deleted=1;
    uint32_t hops_added=0, hops_deleted=1;

    // witness searches
    for(uint32_t i = 0; i < in_deg; i++)
    {
        warthog::graph::edge& e_in = *(n->incoming_begin()+i);

        bool resume = false;
        for(uint32_t j = 0; j < out_deg; j++)
        {
            warthog::graph::edge& e_out = *(n->outgoing_begin()+j);
            if(e_in.node_id_ == e_out.node_id_) { continue; }

            // contraction will introduce a shortcut edge only if
            // the via-n path <in, n, out> is the only shortest path
            // NB: we run a resumable bi-dijkstra since the incoming
            // node is common for all outgoing nodes
            warthog::cost_t via_len = e_in.wt_ + e_out.wt_;
            warthog::cost_t cost_cutoff = e_in.wt_ + e_out.wt_;
            uint32_t expd = total_expansions_;
            warthog::cost_t witness_len = 
                witness_search(e_in.node_id_, e_out.node_id_, cost_cutoff, resume);
            resume = true;
            
            if(witness_len > via_len)
            {
                // track various metrics related to the contraction operation
                uint32_t shortcut_hops = e_in.label_ + e_out.label_;
                edges_added += 1;
                hops_added += shortcut_hops;

                if(!metrics_only)
                {
                    if(this->verbose_)
                    {
                        std::cerr 
                            << "shortcut "
                            << e_in.node_id_ << ", "
                            << e_out.node_id_ 
                            //<< " mid " << node_id 
                            << " wt " << via_len 
                            << " hops " << (e_in.label_ + e_out.label_)
                            << " expd " << (total_expansions_ - expd)
                            << std::endl;
                    }

                    // add shortcuts to search space
                    // NB: these shortcuts can appear on witness paths
                    // for other in/out pairs of node @param node_id
                    // (this is very subtle + different from simply estimating 
                    // the contraction priority of a node)
                    g_->get_node(e_in.node_id_)->add_outgoing(
                    warthog::graph::edge(
                        e_out.node_id_,
                        (warthog::graph::edge_cost_t)via_len, 
                        shortcut_hops));

                    g_->get_node(e_out.node_id_)->add_incoming(
                        warthog::graph::edge(
                            e_in.node_id_,
                            (warthog::graph::edge_cost_t)via_len, 
                            shortcut_hops));

                    // also store the shortcuts separately for later
                    shortcuts.push_back({e_in.node_id_, e_out.node_id_, node_id, via_len});
                }
            }
        }
    }

    if(metrics_only) 
    { c_filter_->clear(); }
    else
    { order_->push_back(node_id); }

    // track more metrics related to the contraction operation
    for(uint32_t i = 0; i < in_deg; i++)
    {
        warthog::graph::edge& e_in = *(n->incoming_begin()+i);
        edges_deleted++;
        hops_deleted += e_in.label_;

        // after contraction, need to know which neighbours to update
        if(!metrics_only)
        { u_filter_->add(e_in.node_id_); }
    }
    for(uint32_t j = 0; j < out_deg; j++)
    {
        warthog::graph::edge& e_out = *(n->outgoing_begin()+j);
        edges_deleted++;
        hops_deleted += e_out.label_;

        // after contraction, need to know which neighbours to update
        if(!metrics_only)
        { u_filter_->add(e_out.node_id_); }
    }

    // compute a value to indicate how attractive it is to 
    // contract @param node_id (lower is better)
    // NB: here we use the heuristic from RoutingKit 
    return (int32_t)
        1 + 
        (1000*height_[node_id]) + 
        (1000*edges_added) / edges_deleted +
        (1000*hops_added) / hops_deleted;
}

size_t
warthog::ch::lazy_graph_contraction::mem()
{
    return
        heap_->mem() +
        alg_->mem() +
        sizeof(*hn_pool_)*g_->get_num_nodes() +
        sizeof(this);
}

void
warthog::ch::lazy_graph_contraction::disconnect_node(uint32_t node_id)
{
    warthog::graph::node* n = g_->get_node(node_id);
    for(warthog::graph::edge_iter it = n->outgoing_begin();
            it != n->outgoing_end(); it++)
    {
        uint32_t nei_id = (*it).node_id_;
        warthog::graph::node* nei = g_->get_node(nei_id);
            
        warthog::graph::edge_iter dit = nei->find_edge(node_id, 
                nei->incoming_begin(), nei->incoming_end());
        if(dit != nei->incoming_end())
        { nei->del_incoming(dit); }
    }

    for(warthog::graph::edge_iter it = n->incoming_begin(); 
            it != n->incoming_end(); it++)
    {
        uint32_t nei_id = (*it).node_id_;
        warthog::graph::node* nei = g_->get_node(nei_id);
        warthog::graph::edge_iter dit = nei->find_edge(node_id, 
                nei->outgoing_begin(), nei->outgoing_end());
        if(dit != nei->outgoing_end())
        { nei->del_outgoing(dit); }
    }
}

