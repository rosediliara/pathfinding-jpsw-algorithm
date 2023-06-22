#include "cfg.h"
#include "contraction.h"
#include "dimacs_parser.h"
#include "geography.h"
#include "graph.h"

#include <cmath>
#include <numeric>
#include <iostream>
#include <filesystem>

// Reorder the edges of a node to match their bearing wrt North
//
// This function could be improved by not having to allocate a two vectors at
// each iteration -- e.g., pre-allocate them in 'orient_edges()' and handle
// bounds.
std::vector<size_t>
orient_node_edges(warthog::graph::xy_graph& g, uint32_t n_id,
                  warthog::graph::ECAP_T degree,
                  warthog::graph::edge_iter start)
{
    std::vector<double> angles(degree);
    int32_t nx, ny;
    g.get_xy(n_id, nx, ny);

    for (uint32_t edge_idx = 0; edge_idx < degree; edge_idx++)
    {
        warthog::graph::edge* e = start + edge_idx;
        // Longitude and latitude of the point
        int32_t lat, lng;
        g.get_xy(e->node_id_, lat, lng);
        // Bearing wrt central node
        angles.at(edge_idx) = warthog::geo::get_bearing_xy(nx, ny, lat, lng);
    }

    // Compute the index ordering given all edges' angles
    std::vector<size_t> idx(degree);
    std::iota(idx.begin(), idx.end(), 0);

    // We sort with $\ge$ as we want a *clockwise rotation*
    stable_sort(idx.begin(), idx.end(), [&](size_t i1, size_t i2)
    {
        return angles[i1] > angles[i2];
    });

    return idx;
}

// Orient graph edges where 0 is the first clockwise edge from the North.
void
orient_edges(warthog::graph::xy_graph& g)
{
    std::vector<warthog::graph::edge> save(warthog::graph::ECAP_MAX);

    for (uint32_t i = 0; i < g.get_num_nodes(); i++)
    {
        warthog::graph::node* n = g.get_node(i);
        warthog::graph::edge_iter head = n->outgoing_begin();
        warthog::graph::ECAP_T out_deg = n->out_degree();

        // Although 'graph::node::clear()' does not free the memory, we still
        // need to save the *values* in a temporary location, otherwise we end
        // up erasing them.
        for(warthog::graph::ECAP_T j = 0; j < out_deg; j++)
        {
            save.at(j) = *(head + j);
        }

        n->clear();
        std::vector<size_t> order = orient_node_edges(g, i, out_deg, head);
        for(auto id : order)
        {
            n->add_outgoing(save.at(id));
        }
        assert(n->out_degree() == out_deg);
    }
}

void
help()
{
    std::cerr 
       << "Converts between the graph format used at the 9th DIMACS Implementation\n"
       << "Challenge and the xy graph format used by Warthog.\n"
       << "Usage: ./dimacs2xy --input [dimacs .co file] [dimacs .gr file]\n"
       << "\t--order enforce azimuth ordering of edges in the resulting graph\n";
}

int 
main(int argc, char** argv)
{
    int order = 0;
	// parse arguments
	warthog::util::param valid_args[] = 
	{
		{"core",  required_argument, 0, 2},
		{"input",  required_argument, 0, 2},
		{"order", no_argument, &order, 1},
		{0,  0, 0, 0}
	};

    warthog::util::cfg cfg;
	cfg.parse_args(argc, argv, "-hc:", valid_args);
    
    if(argc < 2)
    {
		help();
        exit(0);
    }

    std::string gr_file = cfg.get_param_value("input");
    std::string co_file = cfg.get_param_value("input");
    if(gr_file == "")
    {
        std::cerr << "err; missing --input [co file] [gr file]\n";
        return EINVAL;
    }

    if(co_file == "")
    {
        std::cerr << "err; missing --input [co file] [gr file]\n";
        return EINVAL;
    }

    // load 
    warthog::dimacs_parser parser(co_file.c_str(), gr_file.c_str());

    // convert
    warthog::graph::xy_graph g_xy;
    warthog::graph::dimacs_to_xy_graph(parser, g_xy);

    // Order edges by azimuth
    if(order)
    {
        orient_edges(g_xy);
    }

    std::filesystem::path gr_path = gr_file;
    std::filesystem::path co_path = co_file;
    // dump
    std::cout 
        << "# ******************************************************************************\n"
        << "# This XY graph created with the tool \"dimacs2xy\" from the Warthog pathfinding\n"
        << "# library (https://bitbucket.org/dharabor/pathfinding).\n"
        << "# The source input files (in DIMACS format) were: \n"
        << "# " << gr_path.filename() << "\n"
        << "# " << co_path.filename() << "\n"
        << "# ******************************************************************************\n"
        << g_xy;
    return 0;
}

