// programs/mapf.cpp
// 
// Pulls together a variefy of algorithms for 
// solving MAPF instances 
//
// @author: dharabor
// @created: 2019-10-29
//

#include "cbs.h"
#include "cbs_ll_heuristic.h"
#include "cfg.h"
#include "constants.h"
#include "flexible_astar.h"
#include "gridmap.h"
#include "gridmap_expansion_policy.h"
#include "cbs_ll_expansion_policy.h"
#include "manhattan_heuristic.h"
#include "mapf/plan.h"
#include "scenario_manager.h"
#include "timer.h"
#include "sipp/sipp_expansion_policy.h"
#include "sipp/jpst_gridmap.h"
#include "sipp/temporal_jps_expansion_policy.h"
#include "getopt.h"

#include <fstream>
#include <functional>
#include <iomanip>
#include <deque>
#include <sstream>
#include <string>
#include <unordered_map>
#include <memory>

#include "time_constraints.h"

// check computed solutions are optimal
int checkopt = 0;
// print debugging info during search
int verbose = 0;
// display program help on startup
int print_help = 0;

void
help()
{
	std::cerr << "valid parameters:\n"
	<< "\t--alg [algorithm name]\n"
	<< "\t--scen [scenario filename]\n"
	<< "\t--plan [plan filename (existing plan describing paths of higher priority agents)]\n"
	<< "\t--verbose (optional)\n"
    << "\nRecognised values for --alg:\n"
    << "\tcbs_ll, jpst, sipp\n";
}


// run SIPP as a prioritised planning algorithm; assume each instance is an
// agent and plan them one after the other (i.e. the agent order is the
// input order that the instances appear in the scenario file. 
// run until completion or failure.
void
run_sipp(warthog::scenario_manager& scenmgr, std::string alg_name, std::string plan_file)
{
    warthog::gridmap gm(scenmgr.get_experiment(0)->map().c_str());
	warthog::manhattan_heuristic heuristic(gm.header_width(), gm.header_height());
    warthog::sipp_gridmap sipp_map(&gm);
	warthog::sipp_expansion_policy expander(&sipp_map);
    warthog::pqueue_min open;

	warthog::flexible_astar<
		warthog::manhattan_heuristic,
	   	warthog::sipp_expansion_policy,
        warthog::pqueue_min> astar(&heuristic, &expander, &open);

    warthog::mapf::plan theplan;

    // this function blocks (== makes obstacles from) the planned paths of 
    // higher priority agents 
    auto add_higher_priority_plan = 
    [&sipp_map, &expander](warthog::solution& sol) -> void
    {
        for(uint32_t i = 0; i < sol.path_.size(); i++)
        {
            int32_t x, y;
            expander.get_xy(sol.path_.at(i), x, y);

            uint32_t xy_id = (uint32_t)(sol.path_.at(i));
            warthog::cost_t timestep = (sol.path_.at(i) >> 32);

            if(xy_id == (uint32_t)sol.path_.back()) 
            { 
                // agents occupy their target location for one 
                // timestep then disappear (DIFFERENT FROM STANDARD MAPF!!)
                sipp_map.add_obstacle( (uint32_t)x, (uint32_t)y, 
                    timestep, timestep+1, warthog::cbs::move::WAIT);

                if(verbose)
                {
                    std::cerr  << " add obstacle: " << xy_id 
                               << " (" << x << ", " << y << ") @ "
                               << "(" << timestep << ", " << timestep+1 << ") "
                               << " dir " << warthog::cbs::WAIT << std::endl;
                }
                break;
            }

            int32_t nx, ny;
            expander.get_xy((uint32_t)sol.path_.at(i+1), nx, ny);
            warthog::cost_t timestep_next = (sol.path_.at(i+1) >> 32);
            assert(nx == x || ny == y);

            // impute the next move from xy locations of nodes on the path
            warthog::cbs::move direction = warthog::cbs::move::WAIT;
            //if(nx == x && ny == y) {  direction = warthog::cbs::move::WAIT; }
            if(nx == x && ny < y) {  direction = warthog::cbs::move::NORTH; }
            if(nx == x && ny > y) {  direction = warthog::cbs::move::SOUTH; }
            if(nx < x && ny == y) {  direction = warthog::cbs::move::WEST; }
            if(nx > x && ny == y) {  direction = warthog::cbs::move::EAST; }

            sipp_map.add_obstacle(
                (uint32_t)x, (uint32_t)y, timestep, timestep_next, direction);

            if(verbose)
            {
                std::cerr  << " add obstacle: " << xy_id << 
                           " (" << x << ", " << y << ") @ "
                           << "(" << timestep << ", " << timestep_next << ") "
                           << " dir " << warthog::cbs::WAIT << std::endl;
            }

            //if(i == (sol.path_.size()-2))
            //{
            //    // after arriving, agents block their target 
            //    // location indefinitely 
            //    start_time = end_time;
            //    end_time = warthog::COST_MAX;
            //    sipp_map.add_obstacle( (uint32_t)nx, (uint32_t)ny, 
            //        end_time, end_time+1, warthog::cbs::move::WAIT);
            //    if(verbose)
            //    {
            //        std::cerr 
            //            << "add obstacle: (" << nx << ", " << ny 
            //            << ") @ (" << start_time << ", " << end_time
            //            << ") dir WAIT" << std::endl;
            //    }
            //    break;
            //}
        }
    };

    warthog::mapf::plan hoplan; 
    if(plan_file != "")
    {
        std::ifstream ifs(plan_file);
        ifs >> hoplan;
        ifs.close();
    }

    // max number of previous agents still moving when we plan the current
    // agent (their paths become dynamic obstacles for the current agent). 
    //uint32_t sz_max_queue = 5;
    //std::deque<warthog::mapf::plan> obstacles;

    std::cout 
        << "id\talg\texpanded\ttouched\treopen\tsurplus\theapops"
        << "\tnanos\tpcost\tplen\tmap\n";
	for(unsigned int i=0; i < scenmgr.num_experiments(); i++)
	{
		warthog::experiment* exp = scenmgr.get_experiment(i);

		uint32_t startid = exp->starty() * exp->mapwidth() + exp->startx();
		uint32_t goalid = exp->goaly() * exp->mapwidth() + exp->goalx();
        warthog::problem_instance pi(startid, goalid, verbose);
        warthog::solution sol;

        //if(i == 621) { pi.verbose_ = true; }
        astar.get_path(pi, sol);
		std::cout
            << i<<"\t" 
            << alg_name << "\t" 
            << sol.nodes_expanded_ << "\t" 
            << sol.nodes_touched_ << "\t"
            << sol.nodes_reopen_ << "\t"
            << sol.nodes_surplus_ << "\t"
            << sol.heap_ops_ << "\t"
            << sol.time_elapsed_nano_ << "\t"
            << sol.sum_of_edge_costs_ << "\t" 
            << (sol.path_.size()-1) << "\t" 
            << scenmgr.last_file_loaded() 
            << std::endl;

         theplan.paths_.push_back(sol);

         // all subsequent agents need to avoid locations on the
         // newly found plan (i.e. we perform prioritised planning)
         // add_higher_priority_plan(sol);

         // load plans of higher priority agents (if any) and block their
         // temporal locations to avoid collisions
         //if(obstacles.size() == sz_max_queue)
         //{
         //    obstacles.pop_front();
         //}
         //obstacles.push_back(theplan);
         //for(auto& myplan : obstacles)
         //{
         //    add_higher_priority_plan(myplan);
         //}

         if(i < hoplan.paths_.size())
         {
             add_higher_priority_plan(hoplan.paths_.at(i));
         }
	}

    // some extra info about sipp's performance
    size_t max_intervals = 0;
    size_t tot_intervals = 0;
    uint32_t map_sz = sipp_map.gm_->header_height() * sipp_map.gm_->header_width();
    for(uint32_t i = 0; i < map_sz; i++)
    {
        max_intervals = sipp_map.get_all_intervals(i).size() > max_intervals ? 
                        sipp_map.get_all_intervals(i).size() : max_intervals;
        tot_intervals += sipp_map.get_all_intervals(i).size();
    }

    std::string tmp_planfile = scenmgr.last_file_loaded() + "." + alg_name + ".plan";
    std::cerr  << "writing plan to " << tmp_planfile << std::endl;
    std::ofstream ofs(tmp_planfile);
    ofs << theplan;
    ofs.close();
	std::cerr << "done. \n";
	std::cerr << "total memory: "<< astar.mem() + scenmgr.mem() << "\n";
	std::cerr << "max sipp intervals per node: "<< max_intervals << "\n";
	std::cerr << "avg sipp intervals per node: "<< tot_intervals / (double)map_sz << "\n";
}

void
run_jpst(warthog::scenario_manager& scenmgr, std::string alg_name, std::string plan_file)
{
    warthog::gridmap gm(scenmgr.get_experiment(0)->map().c_str());
    warthog::manhattan_heuristic heuristic(gm.header_width(), gm.header_height());
    warthog::jpst_gridmap jpst_gm(&gm);
    warthog::temporal_jps_expansion_policy expander(&jpst_gm);
    warthog::pqueue_min open;

    //jpst_gm.t_gm_r_->print(std::cerr);
    //exit(0);

    warthog::flexible_astar<
    	warthog::manhattan_heuristic,
       	warthog::temporal_jps_expansion_policy,
        warthog::pqueue_min> astar(&heuristic, &expander, &open);

    warthog::mapf::plan theplan;

    // this function blocks (== makes obstacles from) the planned paths of 
    // higher priority agents 
    auto add_higher_priority_plan = 
    [&jpst_gm, &expander](warthog::solution& sol) -> void
    {
        for(uint32_t i = 0; i < sol.path_.size(); i++)
        {
            int32_t x, y;
            expander.get_xy(sol.path_.at(i), x, y);

            uint32_t xy_id = (uint32_t)(sol.path_.at(i));
            warthog::cost_t timestep = (sol.path_.at(i) >> 32);

            if(xy_id == (uint32_t)sol.path_.back()) 
            { 
                // agents occupy their target location for one 
                // timestep then disappear (DIFFERENT FROM STANDARD MAPF!!)
                jpst_gm.add_obstacle( (uint32_t)x, (uint32_t)y, 
                    timestep, timestep+1, warthog::cbs::move::WAIT);

                if(verbose)
                {
                    std::cerr  << " add obstacle: " << xy_id 
                               << " (" << x << ", " << y << ") @ "
                               << "(" << timestep << ", " << timestep+1 << ") "
                               << " dir " << warthog::cbs::WAIT << std::endl;
                }

                break;
            }

            int32_t nx, ny;
            expander.get_xy(sol.path_.at(i+1), nx, ny);
            warthog::cost_t timestep_next = (sol.path_.at(i+1) >> 32);
            assert(nx == x || ny == y);

            // impute the next move from xy locations of nodes on the path
            warthog::cbs::move direction = warthog::cbs::move::WAIT;
            //if(nx == x && ny == y) {  direction = warthog::cbs::move::WAIT; }
            if(nx == x && ny < y) {  direction = warthog::cbs::move::NORTH; }
            if(nx == x && ny > y) {  direction = warthog::cbs::move::SOUTH; }
            if(nx < x && ny == y) {  direction = warthog::cbs::move::WEST; }
            if(nx > x && ny == y) {  direction = warthog::cbs::move::EAST; }

            jpst_gm.add_obstacle(
                (uint32_t)x, (uint32_t)y, timestep, timestep_next, direction);

            if(verbose)
            {
                std::cerr  << " add obstacle: (" << x << ", " << y << ") @ (" 
                           << timestep << ", " << timestep_next << ") dir " 
                           << direction << std::endl;
            }

            //if(i == (sol.path_.size()-2))
            //{
            //    // after arriving, agents block their target 
            //    // location indefinitely 
            //    start_time = end_time;
            //    end_time = warthog::COST_MAX;
            //    jpst_gm.add_obstacle( (uint32_t)nx, (uint32_t)ny, 
            //        end_time, end_time+1, warthog::cbs::move::WAIT);
            //    if(verbose)
            //    {
            //        std::cerr 
            //            << "add obstacle: (" << nx << ", " << ny 
            //            << ") @ (" << start_time << ", " << end_time
            //            << ") dir WAIT" << std::endl;
            //    }
            //    break;
            //}
        }
    };

    // load plans of higher priority agents (if any) and block their
    // temporal locations to avoid collisions
    warthog::mapf::plan hoplan; 
    if(plan_file != "")
    {
        std::ifstream ifs(plan_file);
        ifs >> hoplan;
        ifs.close();
    }

    std::cout 
        << "id\talg\texpanded\ttouched\treopen\tsurplus\theapops"
        << "\tnanos\tpcost\tplen\tmap\n";
    for(unsigned int i=0; i < scenmgr.num_experiments(); i++)
    {
    	warthog::experiment* exp = scenmgr.get_experiment(i);

    	uint32_t startid = exp->starty() * exp->mapwidth() + exp->startx();
    	uint32_t goalid = exp->goaly() * exp->mapwidth() + exp->goalx();
        warthog::problem_instance pi(startid, goalid, verbose);
        warthog::solution sol;

        //if(i == 621) { pi.verbose_ = true; }
        astar.get_path(pi, sol);
    	std::cout
            << i<<"\t" 
            << alg_name << "\t" 
            << sol.nodes_expanded_ << "\t" 
            << sol.nodes_touched_ << "\t"
            << sol.nodes_reopen_ << "\t"
            << sol.nodes_surplus_ << "\t"
            << sol.heap_ops_ << "\t"
            << sol.time_elapsed_nano_ << "\t"
            << sol.sum_of_edge_costs_ << "\t" 
            << (sol.path_.size()-1) << "\t" 
            << scenmgr.last_file_loaded() 
            << std::endl;

        theplan.paths_.push_back(sol);

        // we previously loaded up to k reference plans 
        // these plans are added as additional temporal obstacles.
        // We add them one at a time after planning each of 
        // the first i_th agents with i <= k
        if(i < hoplan.paths_.size())
        {
            add_higher_priority_plan(hoplan.paths_.at(i));
        }
    }

    // some extra info about sipp's performance
    size_t max_intervals = 0;
    size_t tot_intervals = 0;
    uint32_t map_sz = jpst_gm.gm_->header_height() * jpst_gm.gm_->header_width();
    for(uint32_t i = 0; i < map_sz; i++)
    {
        max_intervals = jpst_gm.get_all_intervals(i).size() > max_intervals ? 
                        jpst_gm.get_all_intervals(i).size() : max_intervals;
        tot_intervals += jpst_gm.get_all_intervals(i).size();
    }

    std::string tmp_planfile = scenmgr.last_file_loaded() + "." + alg_name + ".plan";
    std::cerr  << "writing plan to " << tmp_planfile << std::endl;
    std::ofstream ofs(tmp_planfile);
    ofs << theplan;
    ofs.close();
    std::cerr << "done. \n";
    std::cerr << "total memory: "<< astar.mem() + scenmgr.mem() << "\n";
    std::cerr << "max sipp intervals per node: "<< max_intervals << "\n";
    std::cerr << "avg sipp intervals per node: "<< tot_intervals / (double)map_sz << "\n";
}

void
run_cbs_ll(warthog::scenario_manager& scenmgr, std::string alg_name)
{
    warthog::gridmap gm(scenmgr.get_experiment(0)->map().c_str());
	warthog::cbs_ll_heuristic heuristic(&gm);
	warthog::cbs_ll_expansion_policy expander(&gm, &heuristic);
    warthog::pqueue_min open;

	warthog::flexible_astar<
		warthog::cbs_ll_heuristic,
	   	warthog::cbs_ll_expansion_policy,
        warthog::pqueue_min>
            astar(&heuristic, &expander, &open);

    warthog::mapf::plan theplan;

    std::cout 
        << "id\talg\texpanded\ttouched\treopen\tsurplus\theapops"
        << "\tnanos\tpcost\tplen\tmap\n";
	for(unsigned int i=0; i < scenmgr.num_experiments(); i++)
	{
		warthog::experiment* exp = scenmgr.get_experiment(i);

		uint32_t startid = exp->starty() * exp->mapwidth() + exp->startx();
		uint32_t goalid = exp->goaly() * exp->mapwidth() + exp->goalx();
        warthog::problem_instance pi(startid, goalid, verbose);
        warthog::solution sol;

        astar.get_path(pi, sol);
		std::cout
            << i<<"\t" 
            << alg_name << "\t" 
            << sol.nodes_expanded_ << "\t" 
            << sol.nodes_touched_ << "\t"
            << sol.nodes_reopen_ << "\t"
            << sol.nodes_surplus_ << "\t"
            << sol.heap_ops_ << "\t"
            << sol.time_elapsed_nano_ << "\t"
            << sol.sum_of_edge_costs_ << "\t" 
            << (sol.path_.size()-1) << "\t" 
            << scenmgr.last_file_loaded()
            << std::endl;

        // the path of the agent now becomes an obstacle for 
        // the next agent. We assume that agents reach their
        // target and then disappear after one unit of time
        for(uint32_t i = 0; i < sol.path_.size(); i++)
        {
            int32_t x, y;
            uint32_t xy_id = (uint32_t)sol.path_.back();
            uint32_t timestep = (uint32_t)(sol.path_.at(i) >> 32);
            expander.get_xy((uint32_t)sol.path_.back(), x, y);
            if((uint32_t)sol.path_.at(i) == (uint32_t)sol.path_.back()) 
            { 
                expander.get_constraint(xy_id)->v_ = true;
                if(verbose)
                {
                    std::cerr << "add obstacle (" << x << ", " << y << ") @ " << i << std::endl;
                }
                break;
            }

            int32_t nx, ny;
            expander.get_xy(sol.path_.at(i+1), nx, ny);

            // block each cell occupied by the agent
            expander.get_constraint(xy_id)->v_ = true;
            if(verbose)
            {
                std::cerr << "add obstacle (" << x << ", " << y << ") @ " 
                          << i << std::endl;
            }
            assert(expander.get_constraint(xy_id)->v_);

            // block any other agent from swapping positions with the agent
            // (i.e. prevent edge collisions) 
            if(nx !=  x || ny != y)
            {
                // compute the opposite direction
                warthog::cbs::move direction = warthog::cbs::WAIT;
                if(nx == x && ny < y) {  direction = warthog::cbs::move::SOUTH; }
                if(nx == x && ny > y) {  direction = warthog::cbs::move::NORTH; }
                if(nx < x && ny == y) {  direction = warthog::cbs::move::EAST; }
                if(nx > x && ny == y) {  direction = warthog::cbs::move::WEST; }

                // place a constraint that forbids moving in the opposite
                // direction to the agent. the constraint is placed on the 
                // tile corresponding to the kth xy location on the path but
                // at one timestep earlier than the arrival of the current agent
                uint32_t n_xy_id = sol.path_.at(i+1);
                warthog::sn_id_t block_id = ((uint64_t)timestep << 32) | n_xy_id;
                expander.get_constraint(block_id)->e_ 
                    |= (uint8_t)(1 << direction);
            }
        }

        theplan.paths_.push_back(sol);

	}
    std::string tmp_planfile = scenmgr.last_file_loaded() + "." + alg_name + ".plan";
    std::cerr  << "writing plan to " << tmp_planfile << std::endl;
    std::ofstream ofs(tmp_planfile);
    ofs << theplan;
    ofs.close();
	std::cerr << "done. total memory: "<< astar.mem() + scenmgr.mem() << "\n";
}



int 
main(int argc, char** argv)
{
	// parse arguments
	warthog::util::param valid_args[] = 
	{
		{"scen",  required_argument, 0, 0},
		{"plan",  required_argument, 0, 0},
		{"alg",  required_argument, 0, 1},
		{"help", no_argument, &print_help, 1},
		{"checkopt",  no_argument, &checkopt, 1},
		{"verbose",  no_argument, &verbose, 1},
		{"format",  required_argument, 0, 1},
		{0,  0, 0, 0}
	};

	warthog::util::cfg cfg;
	cfg.parse_args(argc, argv, "-f", valid_args);

    if(argc == 1 || print_help)
    {
		help();
        exit(0);
    }

    std::string sfile = cfg.get_param_value("scen");
    std::string alg = cfg.get_param_value("alg");
    std::string planfile = cfg.get_param_value("plan");

	if(alg == "" || sfile == "")
	{
        std::cerr << "Err. Must specify a scenario file and search algorithm. "
                  << "Try --help for options.\n";
		exit(0);
	}

	warthog::scenario_manager scenmgr;
	scenmgr.load_scenario(sfile.c_str());

    if(alg == "cbs_ll")
    {
        run_cbs_ll(scenmgr, alg); 
    }
    else if(alg == "sipp")
    {
        run_sipp(scenmgr, alg, planfile);
    }
    else if(alg == "jpst")
    {
        run_jpst(scenmgr, alg, planfile);
    }
    else
    {
        std::cerr << "err; invalid search algorithm: " << alg << "\n";
    }
}
