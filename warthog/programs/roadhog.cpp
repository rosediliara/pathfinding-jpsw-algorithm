// programs/roadhog.cpp
//
// Pulls together a variety of different algorithms for
// routing on road graphs.
//
// @author: dharabor
// @created: 2016-11-24
//

#include "anytime_astar.h"
#include "apex_filter.h"
#include "bb_filter.h"
#include "bb_labelling.h"
#include "bch_search.h"
#include "bch_expansion_policy.h"
#include "bch_bb_expansion_policy.h"
#include "bidirectional_graph_expansion_policy.h"
#include "bidirectional_search.h"
#include "cfg.h"
#include "constants.h"
#include "contraction.h"
#include "cpd_extractions.h"
#include "cpd_heuristic.h"
#include "cpd_search.h"
#include "depth_first_search.h"
#include "dimacs_parser.h"
#include "euclidean_heuristic.h"
#include "fch_bb_expansion_policy.h"
#include "fch_expansion_policy.h"
#include "flexible_astar.h"
#include "graph_expansion_policy.h"
#include "graph_oracle.h"
#include "cpd_graph_expansion_policy.h"
#include "lazy_graph_contraction.h"
#include "xy_graph.h"
#include "solution.h"
#include "timer.h"
#include "workload_manager.h"
#include "zero_heuristic.h"

#include "getopt.h"

#include <fstream>
#include <functional>
#include <iomanip>
#include <memory>
#include <sstream>
#include <unordered_map>

// check computed solutions are optimal
int checkopt = 0;
// print debugging info during search
int verbose = 0;
// display program help on startup
int print_help = 0;

// suppress the header row when printing results? (default: no)
int suppress_header = 0;

long nruns = 1;

void
help()
{
    std::cerr
    << "==> manual <==\n"
    << "This program solves point-to-point pathfinding problems on road networks.\n"
    << "Road networks are specified as xy-graphs and collections of instances are \n"
    << "specified as p2pz. Both formats are similar to that used at the 9th DIMACS \n"
    << "challenge. The main differences are: \n"
    << "(i) a single file format (cf. gr/co files); \n"
    << "(ii) node/arc ids are zero indexed (cf. 1-indexed) and; \n"
    << "(iii) the enforcement of the triangle inequality for all arc costs \n"
    << "(cf. not, as was the case at the competition)\n\n"
    << ""
    << "The following are valid program parameters:\n"
    << "\t--alg [ algorithm name (required) ]\n"
    << "\t--input [ algorithm-specific input files (omit to show options) ] \n"
    << "\t--problem [ ss or p2p problem file (required) ]\n"
    << "\t--verbose (print debug info; omitting this param means no)\n"
    << "\t--nruns [int (repeats per instance; default=" << nruns << ")]\n"
    << "\nRecognised values for --alg:\n"
    << "\tastar, astar-bb, dijkstra, bi-astar, bi-dijkstra\n"
    << "\tbch, bch-astar, bch-bb, fch, fch-bb\n"
    << "\tdfs, cpd, cpd-search\n";
}

void
run_experiments( warthog::search* algo, std::string alg_name,
        warthog::dimacs_parser& parser, std::ostream& out)
{
    std::cerr << "running experiments\n";
    std::cerr << "(averaging over " << nruns << " runs per instance)\n";

    if(!suppress_header)
    {
        std::cout
            << "id\talg\texpanded\ttouched\treopen\tsurplus\theap_ops"
            << "\tnanos\tpcost\tplen\tmap\n";
    }
    uint32_t exp_id = 0;
    for(auto it = parser.experiments_begin();
            it != parser.experiments_end();
            it++)
    {
        warthog::dimacs_parser::experiment exp = (*it);
        warthog::solution sol;
        warthog::sn_id_t start_id = exp.source;
        warthog::sn_id_t target_id = exp.p2p ? exp.target : warthog::SN_ID_MAX;
        warthog::problem_instance pi(start_id, target_id, verbose);
        uint32_t expanded=0, reopen=0, heap_ops=0, touched=0, surplus=0;
        double nano_time = DBL_MAX;
        for(uint32_t i = 0; i < nruns; i++)
        {
            sol.reset();
            algo->get_path(pi, sol);

            expanded += sol.nodes_expanded_;
            reopen += sol.nodes_reopen_;
            heap_ops += sol.heap_ops_;
            touched += sol.nodes_touched_;
            surplus += sol.nodes_surplus_;
            nano_time = nano_time < sol.time_elapsed_nano_
                            ?  nano_time : sol.time_elapsed_nano_;
        }

        out
            << exp_id++ <<"\t"
            << alg_name << "\t"
            << expanded / nruns << "\t"
            << touched / nruns << "\t"
            << reopen / nruns << "\t"
            << surplus / nruns << "\t"
            << heap_ops / nruns << "\t"
            << (long long)sol.sum_of_edge_costs_ << "\t"
            << (int32_t)((sol.path_.size() == 0) ? -1 : (int32_t)(sol.path_.size()-1)) << "\t"
            << parser.get_problemfile()
            << std::endl;
    }
}

void
run_astar(warthog::util::cfg& cfg,
    warthog::dimacs_parser& parser, std::string alg_name)
{
    std::string xy_filename = cfg.get_param_value("input");
    if(xy_filename == "")
    {
        std::cerr << "parameter is missing: --input [xy-graph file]\n";
        return;
    }

    warthog::graph::xy_graph g;
    std::ifstream ifs(xy_filename);
    ifs >> g;

    warthog::simple_graph_expansion_policy expander(&g);
    warthog::euclidean_heuristic h(&g);
    warthog::pqueue_min open;

    warthog::flexible_astar<
        warthog::euclidean_heuristic,
        warthog::simple_graph_expansion_policy,
        warthog::pqueue_min>
            alg(&h, &expander, &open);

    run_experiments(&alg, alg_name, parser, std::cout);
}

void
run_astar_bb(warthog::util::cfg& cfg,
    warthog::dimacs_parser& parser, std::string alg_name)
{
    std::string xy_filename = cfg.get_param_value("input");
    if(xy_filename == "")
    {
        std::cerr << "parameter is missing: --input [xy-graph file]\n";
        return;
    }

    warthog::graph::xy_graph g;
    std::ifstream ifs(xy_filename);
    ifs >> g;
    ifs.close();

    warthog::label::bb_labelling lab(&g);
    std::string label_filename = xy_filename + ".label.bb";

    ifs.open(label_filename.c_str());
    if(ifs.is_open())
    {
        ifs >> lab;
        ifs.close();
    }
    else
    {
        warthog::util::workload_manager workload(g.get_num_nodes());
        workload.set_all_flags(true);
        lab.precompute(&workload);

        warthog::timer t;
        t.start();
        std::cerr << "saving precompute data to "
            << label_filename << "...\n";

        std::ofstream ofs(label_filename,
                std::ios_base::out|std::ios_base::binary);
        ofs << lab;
        if(!ofs.good())
        {
            std::cerr << "\nerror trying to write to file "
                << label_filename << std::endl;
        }
        ofs.close();
        t.stop();
        std::cerr << "done. time " << t.elapsed_time_nano() / 1e9 << " s\n";

    }

    warthog::bb_filter bbf(&lab);
    warthog::graph_expansion_policy<warthog::bb_filter> expander(&g, &bbf);
    warthog::euclidean_heuristic h(&g);
    warthog::pqueue_min open;

    warthog::flexible_astar<
        warthog::euclidean_heuristic,
        warthog::graph_expansion_policy<warthog::bb_filter>,
        warthog::pqueue_min>
            alg(&h, &expander, &open);

    run_experiments(&alg, alg_name, parser, std::cout);
}

void
run_dijkstra(warthog::util::cfg& cfg,
    warthog::dimacs_parser& parser, std::string alg_name )
{
    // load up the graph
    std::string xy_filename = cfg.get_param_value("input");
    if(xy_filename == "")
    {
        std::cerr << "err; require --input [xy-graph file]\n";
        return;
    }
    warthog::graph::xy_graph g;
    std::ifstream ifs(xy_filename);
    ifs >> g;

    warthog::simple_graph_expansion_policy expander(&g);
    warthog::zero_heuristic h;
    warthog::pqueue<warthog::cmp_less_search_node_f_only, warthog::min_q> open;

    warthog::flexible_astar<
        warthog::zero_heuristic,
        warthog::simple_graph_expansion_policy,
        warthog::pqueue<warthog::cmp_less_search_node_f_only, warthog::min_q>>
            alg(&h, &expander, &open);

    run_experiments(&alg, alg_name, parser, std::cout);
}

void
run_bi_astar( warthog::util::cfg& cfg,
        warthog::dimacs_parser& parser, std::string alg_name)
{
    std::string xy_filename = cfg.get_param_value("input");
    if(xy_filename == "")
    {
        std::cerr << "parameter is missing: --input [xy-graph file]\n";
        return;
    }

    warthog::graph::xy_graph g;
    std::ifstream ifs(xy_filename);
    ifs >> g;
    ifs.close();

    warthog::bidirectional_graph_expansion_policy fexp(&g, false);
    warthog::bidirectional_graph_expansion_policy bexp(&g, true);

    warthog::euclidean_heuristic h(&g);
    warthog::bidirectional_search<
        warthog::euclidean_heuristic,
        warthog::bidirectional_graph_expansion_policy>
            alg(&fexp, &bexp, &h);

    run_experiments(&alg, alg_name, parser, std::cout);
}

void
run_bi_dijkstra( warthog::util::cfg& cfg,
     warthog::dimacs_parser& parser, std::string alg_name )
{
    std::string xy_filename = cfg.get_param_value("input");
    if(xy_filename == "")
    {
        std::cerr << "parameter is missing: --input [xy-graph file]\n";
        return;
    }

    warthog::graph::xy_graph g(0, "", true);
    std::ifstream ifs(xy_filename);
    ifs >> g;
    ifs.close();

    warthog::bidirectional_graph_expansion_policy fexp(&g, false);
    warthog::bidirectional_graph_expansion_policy bexp(&g, true);

    warthog::zero_heuristic h;
    warthog::bidirectional_search<
        warthog::zero_heuristic,
        warthog::bidirectional_graph_expansion_policy>
            alg(&fexp, &bexp, &h);

    run_experiments(&alg, alg_name, parser, std::cout);
}


void
run_bch(warthog::util::cfg& cfg,
        warthog::dimacs_parser& parser, std::string alg_name)
{
    std::string chd_file = cfg.get_param_value("input");
    if(chd_file == "")
    {
        std::cerr << "err; missing chd input file\n";
        return;
    }

    warthog::ch::ch_data chd;
    chd.type_ = warthog::ch::UP_ONLY;
    std::ifstream ifs(chd_file.c_str());
    if(!ifs.is_open())
    {
        std::cerr << "err; invalid path to chd input file\n";
        return;
    }

    ifs >> chd;
    ifs.close();

    warthog::bch_expansion_policy fexp(chd.g_);
    warthog::bch_expansion_policy bexp (chd.g_, true);
    warthog::zero_heuristic h;
    warthog::bch_search<
        warthog::zero_heuristic,
        warthog::bch_expansion_policy>
            alg(&fexp, &bexp, &h);

    run_experiments(&alg, alg_name, parser, std::cout);
}

void
run_bch_backwards_only(warthog::util::cfg& cfg, warthog::dimacs_parser& parser,
        std::string alg_name)
{
    std::string chd_file = cfg.get_param_value("input");
    if(chd_file == "")
    {
        std::cerr << "err; missing chd input file\n";
        return;
    }

    warthog::ch::ch_data chd;
    chd.type_ = warthog::ch::UP_ONLY;
    std::ifstream ifs(chd_file.c_str());
    if(!ifs.is_open())
    {
        std::cerr << "err; invalid path to chd input file\n";
        return;
    }

    ifs >> chd;
    ifs.close();

    warthog::bch_expansion_policy bexp (chd.g_, true);
    warthog::zero_heuristic h;
    warthog::pqueue_min open;

    warthog::flexible_astar<
        warthog::zero_heuristic,
        warthog::bch_expansion_policy,
        warthog::pqueue_min>
            alg(&h, &bexp, &open);

    std::cerr << "running experiments\n";
    std::cerr << "(averaging over " << nruns << " runs per instance)\n";

    if(!suppress_header)
    {
        std::cout
            << "id\talg\texpanded\ttouched\treopen\tsurplus\theap_ops"
            << "\tnanos\tpcost\tplen\tmap\n";
    }
    uint32_t exp_id = 0;
    for(auto it = parser.experiments_begin();
            it != parser.experiments_end();
            it++)
    {
        warthog::dimacs_parser::experiment exp = (*it);
        warthog::solution sol;
        warthog::problem_instance pi(exp.source, warthog::INF32, verbose);
        uint32_t expanded=0, reopen=0, heap_ops=0, touched=0, surplus=0;
        double nano_time = DBL_MAX;
        for(uint32_t i = 0; i < nruns; i++)
        {
            sol.reset();
            alg.get_path(pi, sol);

            expanded += sol.nodes_expanded_;
            heap_ops += sol.heap_ops_;
            touched += sol.nodes_touched_;
            reopen += sol.nodes_reopen_;
            surplus += sol.nodes_surplus_;
            nano_time = nano_time < sol.time_elapsed_nano_
                            ?  nano_time : sol.time_elapsed_nano_;
        }

        std::cout
            << exp_id++ <<"\t"
            << alg_name << "\t"
            << expanded / nruns << "\t"
            << touched / nruns << "\t"
            << reopen / nruns << "\t"
            << surplus / nruns << "\t"
            << heap_ops / nruns << "\t"
            << (long long)sol.sum_of_edge_costs_ << "\t"
            << (int32_t)((sol.path_.size() == 0) ? -1 : (int32_t)(sol.path_.size()-1)) << "\t"
            << parser.get_problemfile()
            << std::endl;
    }
}

void
run_bch_astar(warthog::util::cfg& cfg,
              warthog::dimacs_parser& parser, std::string alg_name)
{
    std::string chd_file = cfg.get_param_value("input");
    if(chd_file == "")
    {
        std::cerr << "err; required --input [chd file]\n";
        return;
    }

    warthog::ch::ch_data chd;
    chd.type_ = warthog::ch::UP_ONLY;
    std::ifstream ifs(chd_file.c_str());
    if(!ifs.is_open())
    {
        std::cerr << "err; invalid path to chd input file\n";
        return;
    }

    ifs >> chd;
    ifs.close();

    warthog::euclidean_heuristic h(chd.g_);
    warthog::bch_expansion_policy fexp(chd.g_);
    warthog::bch_expansion_policy bexp (chd.g_, true);
    warthog::bch_search<
        warthog::euclidean_heuristic,
        warthog::bch_expansion_policy>
            alg(&fexp, &bexp, &h);

    run_experiments(&alg, alg_name, parser, std::cout);
}

void
run_bch_bb(warthog::util::cfg& cfg, warthog::dimacs_parser& parser,
        std::string alg_name)
{
    // load up the contraction hierarchy
    std::string chd_file = cfg.get_param_value("input");
    std::string alg_params = cfg.get_param_value("alg");
    if(chd_file == "")
    {
        std::cerr << "err; require --input [chd file]\n";
        return;
    }

    warthog::ch::ch_data chd;
    chd.type_ = warthog::ch::UP_ONLY;
    std::ifstream ifs(chd_file.c_str());
    if(!ifs.is_open())
    {
        std::cerr << "err; invalid path to chd input file\n";
        return;
    }

    ifs >> chd;
    ifs.close();

    // the "cutoff" tells what percentage of nodes from the hierarchy
    // will have exact labels computed by Dijkstra SSSP search.
    // a cutoff of 0.9 for example omits the bottom 10% of nodes
    // in the hierarchy and preprocesses the remaining 90%.
    // With a cutoff of 1, we preprocess all nodes. With a
    // cutoff of 0, we preprocess none and use DFS labels only.
    int pct_dijkstra = 0;
    std::string label_filename = "label.bb-dfs";
    if(alg_params != "")
    {
        pct_dijkstra = std::stoi(alg_params.c_str());
        if(!(pct_dijkstra >= 0 && pct_dijkstra <= 100))
        {
            std::cerr << "dijkstra percentage must be in range 0-100\n";
            return;
        }
    }
    label_filename += "-dijk-";
    label_filename += std::to_string(pct_dijkstra);

    warthog::label::dfs_labelling lab(&chd);

    // load up the edge label data (or else precompute it)
    label_filename =  chd_file + "." + label_filename;
    ifs.open(label_filename);
    if(ifs.is_open())
    {
        ifs >> lab;
        ifs.close();
    }
    else
    {
        std::cerr
            << "err; label file does not exist: "
            << label_filename << std::endl
            << "you could try to generate it with "
            << "--alg fch-bb " << pct_dijkstra << "\n";
        return;

        //warthog::util::workload_manager workload(chd.g_->get_num_nodes());
        //double cutoff = (((double)pct_dijkstra)/100);
        //uint32_t min_level = (uint32_t)(chd.level_->size()*(1-cutoff));
        //for(size_t i = 0; i < chd.g_->get_num_nodes(); i++)
        //{
        //    if(chd.level_->at(i) >= min_level)
        //    { workload.set_flag((uint32_t)i, true); }
        //}

        //lab.precompute(&workload);

        //warthog::timer t;
        //t.start();
        //std::cerr << "saving precompute data to "
        //    << label_filename << "...\n";

        //std::ofstream out(label_filename,
        //        std::ios_base::out|std::ios_base::binary);
        //out << lab;
        //if(!out.good())
        //{
        //    std::cerr << "\nerror trying to write to file "
        //        << label_filename << std::endl;
        //}
        //out.close();
        //t.stop();
        //std::cerr << "done. time " << t.elapsed_time_nano() / 1e9 << " s\n";

    }

    warthog::bch_bb_expansion_policy fexp(&lab, false);
    warthog::bch_bb_expansion_policy bexp (&lab, true);
    warthog::zero_heuristic h;
    warthog::bch_search<
        warthog::zero_heuristic,
        warthog::bch_bb_expansion_policy>
            alg(&fexp, &bexp, &h);

    run_experiments(&alg, alg_name, parser, std::cout);
}

void
run_fch(warthog::util::cfg& cfg, warthog::dimacs_parser& parser,
        std::string alg_name)
{
    std::string chd_file = cfg.get_param_value("input");
    if(chd_file == "")
    {
        std::cerr << "err; missing chd input file\n";
        return;
    }

    // load up the graph
    warthog::ch::ch_data chd;
    chd.type_ = warthog::ch::UP_DOWN;
    std::ifstream ifs(chd_file.c_str());
    if(!ifs.is_open())
    {
        std::cerr << "err; invalid path to chd input file\n";
        return;
    }

    ifs >> chd;
    ifs.close();

    warthog::fch_expansion_policy fexp(&chd);
    warthog::euclidean_heuristic h(chd.g_);
    warthog::pqueue_min open;

    warthog::flexible_astar<
        warthog::euclidean_heuristic,
        warthog::fch_expansion_policy,
        warthog::pqueue_min>
            alg(&h, &fexp, &open);

    // extra metric; how many nodes do we expand above the apex?
    std::function<uint32_t(warthog::search_node*)> fn_get_apex =
    [&chd, &fexp] (warthog::search_node* n) -> uint32_t
    {
        while(true)
        {
            warthog::search_node* p = fexp.generate(n->get_parent());
            if(!p || chd.level_->at(p->get_id()) < chd.level_->at(n->get_id()))
            { break; }
            n = p;
        }
        return chd.level_->at(n->get_id());
    };

    run_experiments(&alg, alg_name, parser, std::cout);
}

void
run_fch_bb(warthog::util::cfg& cfg, warthog::dimacs_parser& parser,
        std::string alg_name)
{
    std::string alg_params = cfg.get_param_value("alg");
    std::string chd_file = cfg.get_param_value("input");
    if(chd_file == "")
    {
        std::cerr << "err; require --input [chd file]\n";
        return;
    }

    // load up the graph
    warthog::ch::ch_data chd;
    chd.type_ = warthog::ch::UP_DOWN;
    std::ifstream ifs(chd_file.c_str());
    if(!ifs.is_open())
    {
        std::cerr << "err; invalid path to chd input file\n";
        return;
    }

    ifs >> chd;
    ifs.close();

    // define the workload

    // the "cutoff" tells what percentage of nodes from the hierarchy
    // will have exact labels computed by Dijkstra SSSP search.
    // a cutoff of 0.9 for example omits the bottom 10% of nodes
    // in the hierarchy and preprocesses the remaining 90%.
    // With a cutoff of 1, we preprocess all nodes. With a
    // cutoff of 0, we preprocess none and use DFS labels only.
    int pct_dijkstra = 0;
    std::string label_filename = "label.bb-dfs";
    if(alg_params != "")
    {
        pct_dijkstra = std::stoi(alg_params.c_str());
        if(!(pct_dijkstra >= 0 && pct_dijkstra <= 100))
        {
            std::cerr << "dijkstra percentage must be in range 0-100\n";
            return;
        }
    }
    label_filename += "-dijk-";
    label_filename += std::to_string(pct_dijkstra);

    warthog::label::dfs_labelling lab(&chd);

    // load up the edge label data (or else precompute it)
    label_filename =  chd_file + "." + label_filename;
    ifs.open(label_filename);
    if(ifs.is_open())
    {
        ifs >> lab;
        ifs.close();
    }
    else
    {
        warthog::util::workload_manager workload(chd.g_->get_num_nodes());
        double cutoff = (((double)pct_dijkstra)/100);
        uint32_t min_level = (uint32_t)(chd.level_->size()*(1-cutoff));
        for(size_t i = 0; i < chd.g_->get_num_nodes(); i++)
        {
            if(chd.level_->at(i) >= min_level)
            { workload.set_flag((uint32_t)i, true); }
        }

        lab.precompute(&workload);

        warthog::timer t;
        t.start();
        std::cerr << "saving precompute data to "
            << label_filename << "...\n";

        std::ofstream out(label_filename,
                std::ios_base::out|std::ios_base::binary);
        out << lab;
        if(!out.good())
        {
            std::cerr << "\nerror trying to write to file "
                << label_filename << std::endl;
        }
        out.close();
        t.stop();
        std::cerr << "done. time " << t.elapsed_time_nano() / 1e9 << " s\n";

    }

    warthog::fch_bb_expansion_policy fexp(&lab);
    warthog::euclidean_heuristic h(chd.g_);
    warthog::pqueue_min open;

    warthog::flexible_astar<
        warthog::euclidean_heuristic,
        warthog::fch_bb_expansion_policy,
        warthog::pqueue_min>
            alg(&h, &fexp, &open);

    // extra metric; how many nodes do we expand above the apex?
    //std::function<uint32_t(warthog::search_node*)> fn_get_apex =
    //[&chd, &fexp] (warthog::search_node* n) -> uint32_t
    //{
    //    while(true)
    //    {
    //        warthog::search_node* p = fexp.generate(n->get_parent());
    //        if(!p || chd.level_->at(p->get_id()) < chd.level_->at(n->get_id()))
    //        { break; }
    //        n = p;
    //    }
    //    return chd.level_->at(n->get_id());
    //};

    run_experiments(&alg, alg_name, parser, std::cout);
}

std::vector<std::pair<unsigned, warthog::graph::edge>>
load_diff(std::string diff_file, warthog::graph::xy_graph& g)
{
  std::vector<std::pair<unsigned, warthog::graph::edge>> edges;
  std::ifstream ifs(diff_file);
  int num;
  ifs >> num;
  edges.resize(num);
  for (int i=0; i<num; i++)
  {
    int u, v, w;
    ifs >> u >> v >> w;
    edges[i].first = u;
    edges[i].second = warthog::graph::edge(v, w);
  }
  return edges;
}

template<warthog::cpd::symbol SYM>
void
run_cpd_search(warthog::util::cfg& cfg,
    warthog::dimacs_parser& parser, std::string alg_name)
{
    warthog::graph::xy_graph g;
    std::ifstream ifs;
    // We first load the xy_graph and its diff as we need them to be *read* in
    // reverse order.
    std::string xy_filename = cfg.get_param_value("input");
    if(xy_filename == "")
    {
        std::cerr << "parameter is missing: --input graph.xy [graph.xy.diff "
                  << "[graph.xy.cpd]]\n";
        return;
    }

    ifs.open(xy_filename);
    if (!ifs.good())
    {
        std::cerr << "Could not open xy-graph: " << xy_filename << std::endl;
        return;
    }

    ifs >> g;
    ifs.close();

    // Check if we have a second parameter in the --input
    std::string diff_filename = cfg.get_param_value("input");
    if (diff_filename == "")
    {
        diff_filename = xy_filename + ".diff";
    }

    ifs.open(diff_filename);
    if (!ifs.good())
    {
        std::cerr <<
            "Could not open diff-graph: " << diff_filename << std::endl;
        return;
    }

    g.perturb(ifs);
    ifs.close();

    // read the cpd
    warthog::cpd::graph_oracle_base<SYM> oracle(&g);
    std::string cpd_filename = cfg.get_param_value("input");
    if(cpd_filename == "")
    {
        cpd_filename = xy_filename + ".cpd";
    }

    ifs.open(cpd_filename);
    if(ifs.is_open())
    {
        ifs >> oracle;
        ifs.close();
    }
    else
    {
        std::cerr << "Could not find CPD file '" << cpd_filename << "'\n";
        return;
    }

    warthog::simple_graph_expansion_policy expander(&g);
    warthog::cpd_heuristic_base<SYM> h(&oracle, 1.0);
    warthog::pqueue_min open;

    warthog::cpd_search<
        warthog::cpd_heuristic_base<SYM>,
        warthog::simple_graph_expansion_policy,
        warthog::pqueue_min>
            alg(&h, &expander, &open);

    // Set options for CPD searh
    std::stringstream ss;
    std::string scale = cfg.get_param_value("fscale");
    std::string tlim = cfg.get_param_value("uslim");
    std::string moves = cfg.get_param_value("kmoves");
    double f_scale;
    uint32_t us_lim;
    uint32_t k_moves;

    if (scale != "")
    {
        ss << scale;
        ss >> f_scale;

        if (f_scale > 0.0)
        {
            alg.set_quality_cutoff(f_scale);
        }
    }

    if (tlim != "")
    {
        ss << tlim;
        ss >> us_lim;

        if (us_lim > 0)
        {
            alg.set_max_us_cutoff(us_lim);
        }
    }

    if (moves != "")
    {
        ss << moves;
        ss >> k_moves;

        if (k_moves > 0)
        {
            alg.set_max_k_moves(k_moves);
        }
    }

    run_experiments(&alg, alg_name, parser, std::cout);
}

template<warthog::cpd::symbol SYM>
void
run_cpd(warthog::util::cfg& cfg,
    warthog::dimacs_parser& parser, std::string alg_name)
{
    std::string xy_filename = cfg.get_param_value("input");
    if(xy_filename == "")
    {
        std::cerr << "parameter is missing: --input [xy-graph file]\n";
        return;
    }


    warthog::graph::xy_graph g;
    std::ifstream ifs(xy_filename);
    ifs >> g;
    ifs.close();

    warthog::cpd::graph_oracle_base<SYM> oracle(&g);
    std::string cpd_filename = cfg.get_param_value("input");
    if(cpd_filename == "")
    {
        cpd_filename = xy_filename + ".cpd";
    }

    ifs.open(cpd_filename);
    if(ifs.is_open())
    {
        ifs >> oracle;
        ifs.close();
    }
    else
    {
        std::cerr << "Could not find CPD file '" << cpd_filename << "'\n";
        return;
    }

    warthog::cpd_extractions_base<SYM> cpd_extract(&g, &oracle);

    run_experiments(&cpd_extract, alg_name, parser, std::cout);
}

void
run_dfs(warthog::util::cfg& cfg,
    warthog::dimacs_parser& parser, std::string alg_name)
{
    std::string xy_filename = cfg.get_param_value("input");
    if(xy_filename == "")
    {
        std::cerr << "parameter is missing: --input [xy-graph file]\n";
        return;
    }

    warthog::graph::xy_graph g;

    std::ifstream ifs(xy_filename);
    ifs >> g;
    ifs.close();

    warthog::simple_graph_expansion_policy expander(&g);
    warthog::zero_heuristic h;
    warthog::pqueue_min open;

    warthog::depth_first_search<
        warthog::zero_heuristic,
        warthog::simple_graph_expansion_policy,
        warthog::pqueue_min>
            alg(&h, &expander, &open);

    run_experiments(&alg, alg_name, parser, std::cout);
}

void
run_cpd_dfs(warthog::util::cfg& cfg,
    warthog::dimacs_parser& parser, std::string alg_name)
{
    std::string xy_filename = cfg.get_param_value("input");
    if(xy_filename == "")
    {
        std::cerr << "parameter is missing: --input [xy-graph file]\n";
        return;
    }

    warthog::graph::xy_graph g;

    std::ifstream ifs(xy_filename);
    ifs >> g;
    ifs.close();

    warthog::cpd::graph_oracle oracle(&g);
    std::string cpd_filename = cfg.get_param_value("input");
    if(cpd_filename == "")
    {
        cpd_filename = xy_filename + ".cpd";
    }

    ifs.open(cpd_filename);
    if(ifs.is_open())
    {
        ifs >> oracle;
    }
    else
    {
        std::cerr << "Could not find CPD file '" << cpd_filename << "'\n";
        return;
    }

    warthog::cpd_graph_expansion_policy expander(&oracle);
    warthog::zero_heuristic h;
    warthog::pqueue_min open;

    warthog::depth_first_search<
        warthog::zero_heuristic,
        warthog::cpd_graph_expansion_policy,
        warthog::pqueue_min>
            alg(&h, &expander, &open);

    run_experiments(&alg, alg_name, parser, std::cout);
}

void
run_dimacs(warthog::util::cfg& cfg)
{
    std::string alg_name = cfg.get_param_value("alg");
    std::string par_nruns = cfg.get_param_value("nruns");
    std::string problemfile = cfg.get_param_value("problem");

    if((alg_name == ""))
    {
        std::cerr << "parameter is missing: --alg\n";
        return;
    }

    if((problemfile == ""))
    {
        std::cerr << "parameter is missing: --problem\n";
        return;
    }

    if(par_nruns != "")
    {
       char* end;
       nruns = strtol(par_nruns.c_str(), &end, 10);
    }

    warthog::dimacs_parser parser;
    parser.load_instance(problemfile.c_str());
    if(parser.num_experiments() == 0)
    {
        std::cerr << "err; specified problem file contains no instances\n";
        return;
    }

    if(alg_name == "dijkstra")
    {
        run_dijkstra(cfg, parser, alg_name);
    }
    else if(alg_name == "astar")
    {
        run_astar(cfg, parser, alg_name);
    }
    else if(alg_name == "astar-bb")
    {
        run_astar_bb(cfg, parser, alg_name);
    }
    else if(alg_name == "bi-dijkstra")
    {
        run_bi_dijkstra(cfg, parser, alg_name);
    }
    else if(alg_name == "bi-astar")
    {
        run_bi_astar(cfg, parser, alg_name);
    }
    else if(alg_name == "bch")
    {
        run_bch(cfg, parser, alg_name);
    }
    else if(alg_name == "bchb")
    {
        run_bch_backwards_only(cfg, parser, alg_name);
    }
    else if(alg_name == "bch-astar")
    {
        run_bch_astar(cfg, parser, alg_name);
    }
    else if(alg_name == "bch-bb")
    {
        run_bch_bb(cfg, parser, alg_name);
    }
    else if(alg_name == "fch")
    {
        run_fch(cfg, parser, alg_name);
    }
    else if(alg_name == "fch-bb")
    {
        run_fch_bb(cfg, parser, alg_name);
    }
    else if(alg_name == "cpd-search")
    {
        run_cpd_search<warthog::cpd::FORWARD>(cfg, parser, alg_name);
    }
    else if(alg_name == "table-search")
    {
        run_cpd_search<warthog::cpd::TABLE>(cfg, parser, alg_name);
    }
    else if(alg_name == "rev-table-search")
    {
        run_cpd_search<warthog::cpd::REV_TABLE>(cfg, parser, alg_name);
    }
    else if(alg_name == "cpd")
    {
        run_cpd<warthog::cpd::FORWARD>(cfg, parser, alg_name);
    }
    else if(alg_name == "rev-cpd")
    {
        run_cpd<warthog::cpd::REVERSE>(cfg, parser, alg_name);
    }
    else if(alg_name == "table")
    {
        run_cpd<warthog::cpd::TABLE>(cfg, parser, alg_name);
    }
    else if(alg_name == "rev-table")
    {
        run_cpd<warthog::cpd::REV_TABLE>(cfg, parser, alg_name);
    }
    else if(alg_name == "dfs")
    {
        run_dfs(cfg, parser, alg_name);
    }
    else if(alg_name == "cpd-dfs")
    {
        run_cpd_dfs(cfg, parser, alg_name);
    }
    else
    {
        std::cerr << "invalid search algorithm\n";
    }
}


int
main(int argc, char** argv)
{
    // parse arguments
    warthog::util::param valid_args[] =
    {
        {"alg",  required_argument, 0, 1},
        {"nruns",  required_argument, 0, 1},
        {"help", no_argument, &print_help, 1},
        {"checkopt",  no_argument, &checkopt, 1},
        {"verbose",  no_argument, &verbose, 1},
        {"noheader",  no_argument, &suppress_header, 1},
        {"input",  required_argument, 0, 1},
        {"problem",  required_argument, 0, 1},
        {"fscale", required_argument, 0, 1},
        {"uslim", required_argument, 0, 1},
        {"kmoves", required_argument, 0, 1},
        {0,  0, 0, 0}
    };

    warthog::util::cfg cfg;
    cfg.parse_args(argc, argv, "-f", valid_args);

    if(argc == 1 || print_help)
    {
        help();
        exit(0);
    }

    run_dimacs(cfg);
}
