#ifndef WARTHOG_COST_TABLE_H
#define WARTHOG_COST_TABLE_H

// cost_table.h
//
// A utility for mapping terrain types to their weights in weighted maps.
//
// @author: Mark Carlson
// @created: 2022-06-30
//

#include "labelled_gridmap.h"
#include <array>

namespace warthog {

class cost_table
{
public:
    cost_table()
    {
        costs_[0] = 0.0;
        for (int i = 1; i < 256; i++) {
            costs_[i] = 0.0 / 0.0;
        }
    }

    cost_table(const char* filename);

    warthog::cost_t lowest_cost(warthog::vl_gridmap& map);

    warthog::cost_t& operator[](uint8_t index)
    {
        return costs_[index];
    }

private:
    std::array<warthog::cost_t, 256> costs_;
};

}

#endif