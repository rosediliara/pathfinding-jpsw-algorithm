#define CATCH_CONFIG_RUNNER

#include "catch.hpp"
#include "graph_expansion_policy.h"
#include "xy_graph.h"
#include "cpd_search.h"
#include "dimacs_parser.h"
#include "graph_oracle.h"
#include "cpd_heuristic.h"
#include "cpd_extractions.h"
#include "cast.h"

using namespace std;

int
main(int argv, char* args[])
{
    Catch::Session session;
    int res = session.run(argv, args);
    return res;
}

bool
perturb_edge(
    warthog::graph::xy_graph* g, uint32_t from, uint32_t to, uint32_t coef=2)
{
    warthog::graph::node* n = g->get_node(from);

    for (uint32_t j = 0; j < n->out_degree(); j++)
    {
        warthog::graph::edge* e = (n->outgoing_begin() + j);

        if (e->node_id_ == to)
        {
            e->wt_ *= coef;
            return true;
        }
    }

    return false;
}

SCENARIO("Test CPD A* on a square matrix", "[cpd][square][astar]")
{
    string map_name = "square01.map";
    warthog::graph::xy_graph g;
    warthog::gridmap d(map_name.c_str());

    warthog::graph::gridmap_to_xy_graph(&d, &g, false);
    warthog::simple_graph_expansion_policy expander(&g);
    warthog::cpd::graph_oracle oracle(&g);
    std::string cpd_filename = "square01.cpd";
    std::ifstream ifs(cpd_filename);

    if(ifs.is_open())
    {
        ifs >> oracle;
    }
    else
    {
        std::cerr << "Could not find CPD file '" << cpd_filename << "'\n";
        return;
    }

    warthog::sn_id_t start = 0;
    warthog::sn_id_t goal = 19;
    // For some reason cannot use `warthog::sn_id_t`
    int32_t x, y;
    // Cannot cut corners
    warthog::cost_t cost = warthog::ONE *
        (warthog::DBL_ONE * 6 + warthog::DBL_ROOT_TWO);

    g.get_xy(start, x, y);
    REQUIRE(x == 0);
    REQUIRE(y == 0);
    g.get_xy(goal, x, y);
    REQUIRE(x == 400000);  // (5 - 1) * GRID_TO_GRAPH_SCALE_FACTOR
    REQUIRE(y == 400000);

    GIVEN("A CPD heuristic")
    {
        warthog::cpd_heuristic h(&oracle);
        warthog::pqueue_min open;
        warthog::solution sol;
        warthog::cpd_search<
            warthog::cpd_heuristic,
            warthog::simple_graph_expansion_policy>
                astar(&h, &expander, &open);

        THEN("We should be searching quite fast")
        {
            warthog::problem_instance pi(start, goal, true);

            astar.get_path(pi, sol);

            REQUIRE(sol.sum_of_edge_costs_ == cost);
            REQUIRE(sol.nodes_expanded_ == 0);
        }
    }

    GIVEN("A perturbation")
    {
        warthog::cpd_heuristic h(&oracle);
        warthog::sn_id_t hit = 6;
        warthog::graph::xy_graph p;
        warthog::graph::gridmap_to_xy_graph(&d, &p, false);
        REQUIRE(perturb_edge(&p, hit, 10));
        g.perturb(p);

        THEN("The last unperturbed node is before the perturbation")
        {
            warthog::sn_id_t last;
            warthog::cost_t lb;
            warthog::cost_t ub;
            h.h(start, goal, lb, ub, last);
            REQUIRE(last == hit);
        }

        THEN("The node is cached and other paths have the same")
        {
            warthog::sn_id_t last;
            warthog::cost_t lb;
            warthog::cost_t ub;

            h.h(start, goal, lb, ub, last);
            REQUIRE(last == hit);

            h.h(2, goal, lb, ub, last);
            REQUIRE(last == hit);
        }
    }
}

SCENARIO("Test CPD search on a modified cross.", "[cpd][astar][cross]")
{
    // Load map
    string map_name = "cross01.map";
    warthog::graph::xy_graph g;
    warthog::gridmap d(map_name.c_str());
    warthog::graph::gridmap_to_xy_graph(&d, &g, false);

    // Load CPD
    warthog::cpd::graph_oracle oracle(&g);
    std::string cpd_filename = map_name + ".cpd";
    warthog::cpd_heuristic h(&oracle, 1.0);
    std::ifstream ifs(cpd_filename);

    if(ifs.is_open())
    {
        ifs >> oracle;
    }
    else
    {
        std::cerr << "Could not find CPD file '" << cpd_filename << "'\n";
        return;
    }

    // Search algorithm
    warthog::simple_graph_expansion_policy expander(&g);
    warthog::pqueue_min open;
    warthog::solution sol;
    warthog::cpd_search<
        warthog::cpd_heuristic,
        warthog::simple_graph_expansion_policy>
        astar(&h, &expander, &open);

    // Instance data
    warthog::sn_id_t start = 0;
    warthog::sn_id_t goal = 19;
    // Cannot cut corners
    warthog::cost_t cost = warthog::ONE *
        (warthog::DBL_ONE * 6 + warthog::DBL_ROOT_TWO);
    std::vector<warthog::sn_id_t> optipath = {19, 14, 10, 8, 3, 2, 1, 0};

    GIVEN("The unperturbed graph")
    {
        THEN("The path is correct")
        {
            warthog::problem_instance pi(start, goal, true);

            astar.get_path(pi, sol);

            REQUIRE(sol.path_ == optipath);
            REQUIRE(sol.sum_of_edge_costs_ == cost);
            REQUIRE(sol.nodes_expanded_ <= 1);
        }
    }

    GIVEN("A perturbation not on the optimal path")
    {
        warthog::graph::xy_graph p;
        warthog::graph::gridmap_to_xy_graph(&d, &p, false);
        REQUIRE(perturb_edge(&p, 0, 5));
        g.perturb(p);

        THEN("The search is not affected")
        {
            warthog::problem_instance pi(start, goal, true);

            astar.get_path(pi, sol);

            REQUIRE(sol.path_ == optipath);
            REQUIRE(sol.nodes_expanded_ <= 1);
        }
    }

    GIVEN("A perturbation on the optimal path")
    {
        // Affect (0, 0) -> (0, 1)
        warthog::graph::xy_graph p;
        warthog::graph::gridmap_to_xy_graph(&d, &p, false);
        REQUIRE(perturb_edge(&p, 0, 1));
        g.perturb(p);

        THEN("The search is affected but cost is the same")
        {
            warthog::problem_instance pi(start, goal, true);

            astar.get_path(pi, sol);

            REQUIRE(sol.path_ != optipath);
            REQUIRE(sol.sum_of_edge_costs_ == cost);
            // We actually only expand one node: the start; we early stop when
            // fetching 5 from the queue.
            REQUIRE(sol.nodes_expanded_ == 1);
        }
    }

    GIVEN("A perturbed map")
    {
        // Now we will modify the map enough to warrant search: the diagonal
        // moves at the corner -- which are the only ones that offer
        // alternatives.
        warthog::graph::xy_graph p;
        warthog::graph::gridmap_to_xy_graph(&d, &p, false);
        REQUIRE(perturb_edge(&p, 3, 8));
        REQUIRE(perturb_edge(&p, 11, 16));
        g.perturb(p);

        THEN("The optimal path is now two straight lines")
        {
            warthog::problem_instance pi(start, goal, true);

            astar.get_path(pi, sol);

            REQUIRE(sol.path_ != optipath);
            REQUIRE(sol.path_.size() == (optipath.size() + 1));
            REQUIRE(sol.sum_of_edge_costs_ == warthog::ONE * 8);
        }
       
        WHEN("We have a bit of time to search")
        {
            astar.set_max_ms_cutoff(1);

            THEN("We explore some of the tree")
            {
                warthog::problem_instance pi(start, goal, true);
                astar.get_path(pi, sol);

                REQUIRE(sol.nodes_expanded_ > 0);
            }
        }

        WHEN ("We limit the number of k moves")
        {
            astar.set_max_k_moves(2);

            THEN("We can still explore some")
            {
                warthog::problem_instance pi(start, goal, true);
                astar.get_path(pi, sol);

                // TODO have the actual number of nodes of some sort?
                REQUIRE(sol.nodes_expanded_ > 0);
            }
        }
    }

    GIVEN("A sub-optimal heuristic")
    {
        // 100% quality cutoff
        astar.set_quality_cutoff(1.0);

        WHEN("The graph is perturbed within acceptable range")
        {
            // Modify optimal edge (0, 0) -> (0, 1), but the total cost of the
            // path is within the acceptance factor.
            warthog::graph::xy_graph p;
            warthog::graph::gridmap_to_xy_graph(&d, &p, false);
            REQUIRE(perturb_edge(&p, 0, 1));
            g.perturb(p);

            THEN("The search is not modified")
            {
                warthog::problem_instance pi(start, goal, true);

                astar.get_path(pi, sol);

                REQUIRE(sol.path_ == optipath);
            }
        }
    }
}

SCENARIO("Test reverse CPD A* on a square matrix", "[reverse][square][astar]")
{
    string map_name = "square01.map";
    warthog::graph::xy_graph g;
    warthog::gridmap d(map_name.c_str());
    warthog::graph::gridmap_to_xy_graph(&d, &g, false);
    warthog::simple_graph_expansion_policy expander(&g);
    warthog::cpd::graph_oracle_base<warthog::cpd::REVERSE> oracle(&g);
    // Needs to be created by hand
    std::string cpd_filename = "square01.xy-rev.cpd";
    std::ifstream ifs(cpd_filename);

    REQUIRE(ifs.is_open());
    ifs >> oracle;

    warthog::sn_id_t start = 0;
    warthog::sn_id_t goal = 19;
    warthog::problem_instance pi(start, goal, true);
    // For some reason cannot use `warthog::sn_id_t`
    int32_t x, y;
    // Cannot cut corners
    warthog::cost_t cost = warthog::ONE *
        (warthog::DBL_ONE * 6 + warthog::DBL_ROOT_TWO);

    g.get_xy(start, x, y);
    REQUIRE(x == 0);
    REQUIRE(y == 0);
    g.get_xy(goal, x, y);
    REQUIRE(x == 400000);  // (5 - 1) * GRID_TO_GRAPH_SCALE_FACTOR
    REQUIRE(y == 400000);

    GIVEN("A reverse CPD")
    {
        THEN("We can perform extractions")
        {
            warthog::cpd_extractions_base<warthog::cpd::REVERSE> cpd(
                &g, &oracle);
            warthog::solution sol;

            cpd.get_pathcost(pi, sol);
            REQUIRE(sol.sum_of_edge_costs_ == cost);
        }

        THEN("Search should not expand nodes without congestion")
        {
            warthog::solution sol;
            warthog::pqueue_min open;
            warthog::cpd_heuristic_base<warthog::cpd::REVERSE> h(&oracle);
            warthog::cpd_search<
                warthog::cpd_heuristic_base<warthog::cpd::REVERSE>,
                warthog::simple_graph_expansion_policy>
                    astar(&h, &expander, &open);

            astar.set_max_expansions_cutoff(10);
            astar.get_path(pi, sol);

            REQUIRE(sol.sum_of_edge_costs_ == cost);
            REQUIRE(sol.nodes_expanded_ == 0);
        }
    }
}

SCENARIO("Test table lookup A* on a square matrix", "[table][square][astar]")
{
    string map_name = "square01.map";
    warthog::graph::xy_graph g;
    warthog::gridmap d(map_name.c_str());
    warthog::graph::gridmap_to_xy_graph(&d, &g, false);
    warthog::simple_graph_expansion_policy expander(&g);
    warthog::cpd::graph_oracle_base<warthog::cpd::REV_TABLE> oracle(&g);
    // Needs to be created by hand
    std::string cpd_filename = "square01.xy-table.cpd";
    std::ifstream ifs(cpd_filename);

    REQUIRE(ifs.is_open());
    ifs >> oracle;

    warthog::sn_id_t start = 0;
    warthog::sn_id_t goal = 19;
    warthog::problem_instance pi(start, goal, true);
    // For some reason cannot use `warthog::sn_id_t`
    int32_t x, y;
    // Cannot cut corners
    warthog::cost_t cost = warthog::ONE *
        (warthog::DBL_ONE * 6 + warthog::DBL_ROOT_TWO);

    g.get_xy(start, x, y);
    REQUIRE(x == 0);
    REQUIRE(y == 0);
    g.get_xy(goal, x, y);
    REQUIRE(x == 400000);  // (5 - 1) * GRID_TO_GRAPH_SCALE_FACTOR
    REQUIRE(y == 400000);

    GIVEN("A table CPD")
    {
        THEN("We can perform extractions")
        {
            warthog::cpd_extractions_base<warthog::cpd::REV_TABLE> cpd(
                &g, &oracle);
            warthog::solution sol;

            cpd.get_pathcost(pi, sol);
            REQUIRE(sol.sum_of_edge_costs_ == cost);
        }

        THEN("Search should not expand nodes without congestion")
        {
            warthog::solution sol;
            warthog::pqueue_min open;
            warthog::cpd_heuristic_base<warthog::cpd::REV_TABLE> h(&oracle);
            warthog::cpd_search<
                warthog::cpd_heuristic_base<warthog::cpd::REV_TABLE>,
                warthog::simple_graph_expansion_policy>
                    astar(&h, &expander, &open);

            astar.set_max_expansions_cutoff(10);
            astar.get_path(pi, sol);

            REQUIRE(sol.sum_of_edge_costs_ == cost);
            REQUIRE(sol.nodes_expanded_ == 0);
        }
    }
}
