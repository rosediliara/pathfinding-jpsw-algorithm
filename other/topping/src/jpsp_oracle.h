#ifndef WARTHOG_JPSP_ORACLE_H
#define WARTHOG_JPSP_ORACLE_H

// jpsp_oracle.cpp
// 
// An oracle for identifying jump point successors. This implementation
// is similar to the one used by the JPS+ algorithm.
//
// The oracle is constructed by pre-processing the gridmap and trying to jump
// in every direction from every traversable tile. If there exists a jump
// point in the given direction the id of the jump point is recorded and 
// stored. If there is no jump point a dummy value is stored instead. 
//
// There are eight such ids kept for every traversable tile on the map 
// (one for each possible jump direction). The reason we keep ids instead of 
// boolean values is to allow the search to jump directly to the next jump 
// point.
//
// @author dharabor
// @created: 2017-01-17
//

#include "gridmap.h"
#include "jps.h"
#include "offline_jump_point_locator.h"

#include <cassert>
#include <unordered_map>

namespace warthog
{

class jpsp_oracle
{
    public:
        jpsp_oracle(warthog::gridmap* gm);
        virtual ~jpsp_oracle();

        // return the id of the next jump point from location 
        // (@param x, @param y) when traveling in direction @param d.
        // Given a jump point at location (jp_x, jp_y) its id is
        // calculated as jp_y * [map_width] + jp_x;
        // if there is no next jump point the return value is warthog::INF.
        uint32_t
        next_jump_point(uint32_t x, uint32_t y, warthog::jps::direction d)
        {
            uint32_t loc_id = gm_->to_padded_id(x, y);
            uint32_t jp_id;
            warthog::cost_t num_steps;
            oracle_->jump(d, loc_id, goal_id_, jp_id, num_steps);

            if(jp_id != warthog::INF)
            {
                uint32_t jp_x, jp_y;
                gm_->to_unpadded_xy(jp_id, jp_x, jp_y);
                jp_id = jp_y * gm_->header_width() + jp_x;
            }
            return num_steps;
        }

        // call this at the start of every new pathfinding query;
        // before any calls to ::next_jump_point
        inline void
        set_goal_location(uint32_t x, uint32_t y)
        {
            goal_id_ = gm_->to_padded_id(x, y);
        }

        // sanity check; every jump point is correctly identified
        // by the oracle. similarly, the oracle never says a node 
        // is a jump point when it shouldn't.
        // returns true if the oracle passes the test and false otherwise
        bool
        sanity_check();

    private:
        warthog::offline_jump_point_locator* oracle_;
        warthog::gridmap* gm_;
        uint32_t goal_id_;

        void
        preprocess();
};

}

#endif

