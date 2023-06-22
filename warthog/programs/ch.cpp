#include "cfg.h"
#include "bch_expansion_policy.h"
#include "dimacs_parser.h"
#include "graph.h"
#include "lazy_graph_contraction.h"
#include "xy_graph.h"

#include <iostream>
#include <string>

int verbose=false;
int verify=false;
int has_input=0;
warthog::util::cfg cfg;

void
help()
{
    std::cerr
        << "=>  manual <== \n"
        << "This program creates a contraction hierarchy from an xy input graph \n"
        << "(a custom format, similar to that used at the 9th DIMACS challenge)\n";
	std::cerr << "The following are valid program parameters:\n"
    << "\t--input [xy graph file]\n"
	<< "\t--verbose (optional; prints debugging info when compiled with debug symbols)\n"
	<< "\t--verify (optional; verify lazy priorities before contraction.\n"
    << "\t          slow but can produce less shortcut edges)\n";
}

void
contract_graph()
{
    warthog::ch::ch_data chd(true);
    std::string outfile;
    std::string xy_file;

    // load up the input graph/grid
    if(cfg.get_num_values("input") <= 1)
    {
        xy_file = cfg.get_param_value("input");
        std::ifstream ifs(xy_file);
        ifs >> *chd.g_;
        chd.g_->set_filename(xy_file.c_str());
        chd.up_degree_->resize(chd.g_->get_num_nodes(), 0);
    }
    else
    {
        std::cerr << "err; input graph (--input [xy graph file]) not specified\n";
        return;
    }

    // create a new contraction hierarchy with dynamic node ordering
    warthog::ch::lazy_graph_contraction contractor;
    contractor.set_verbose(verbose);
    contractor.contract(&chd, verify);


    // save the result
    outfile = chd.g_->get_filename();
    outfile += ".ch";
    std::cerr << "saving contracted graph to file " << outfile << std::endl;
    std::ofstream ofs(outfile.c_str());
    ofs << chd;
    ofs.close();
    std::cerr << "all done!\n";
}

int main(int argc, char** argv)
{

	// parse arguments
	warthog::util::param valid_args[] =
	{
		{"verbose", no_argument, &verbose, 1},
		{"verify", no_argument, &verify, 1},
		{"input",  required_argument, &has_input, 1},
		{0,  0, 0, 0}
	};
	cfg.parse_args(argc, argv, "abc:d:", valid_args);

    if(!has_input)
    {
		help();
        exit(0);
    }

    contract_graph();
}
