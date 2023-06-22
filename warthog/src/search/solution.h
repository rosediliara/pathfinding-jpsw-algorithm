#ifndef WARTHOG_SOLUTION_H
#define WARTHOG_SOLUTION_H

// a wrapper for solutions found by algorithms implemented
// in the warthog library
//
// @author: dharabor
// @created: 2017-05-03
//

#include <vector>
#include <ostream>
#include "constants.h"

namespace warthog
{

class solution
{
    public:
        solution()
        { 
            reset (); 
            path_.reserve(1024); 
        }

        solution(const solution& other) :   
            sum_of_edge_costs_(other.sum_of_edge_costs_), 
            time_elapsed_nano_(other.time_elapsed_nano_), 
            nodes_expanded_(other.nodes_expanded_),
            nodes_touched_(other.nodes_touched_),
            nodes_surplus_(other.nodes_surplus_),
            nodes_reopen_(other.nodes_reopen_),
            heap_ops_(other.heap_ops_),
            path_(other.path_)
        { }

        inline void
        print(std::ostream& out)
        {
            print_metrics(out);
            out << std::endl;
            print_path(out);
            out << std::endl;
        }

        inline void
        print_metrics(std::ostream& out)
        {
            out 
                << "sum_of_edge_costs=" << sum_of_edge_costs_ 
                << "time_elapsed_nano=" << time_elapsed_nano_ 
                << std::endl
                << "nodes expanded=" << nodes_expanded_ 
                << " touched = " << nodes_touched_ 
                << " reopened = " << nodes_reopen_ 
                << " surplus= " << nodes_surplus_
                << " heap-ops= " << heap_ops_;
        }

        inline void
        print_path(std::ostream& out)
        {
            out << "path=";
           for(auto &state : path_) { out << state << " "; }
        }

        void
        reset()
        {
            sum_of_edge_costs_ = warthog::COST_MAX;
            time_elapsed_nano_ = 0;
            nodes_expanded_ = 0; 
            nodes_touched_= 0;
            nodes_surplus_ = 0;
            nodes_reopen_ = 0;
            heap_ops_ = 0;
            path_.clear();
        }

        // metrics
        warthog::cost_t sum_of_edge_costs_;
        double time_elapsed_nano_;
        uint32_t nodes_expanded_;
        uint32_t nodes_touched_;
        uint32_t nodes_surplus_;
        uint32_t nodes_reopen_;
        uint32_t heap_ops_;

        // the sequence of states that comprise 
        // a solution path
        std::vector<warthog::sn_id_t> path_;
};

}

std::ostream& operator<<(std::ostream& str, warthog::solution& sol);

#endif
