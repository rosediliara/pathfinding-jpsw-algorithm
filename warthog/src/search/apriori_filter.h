#ifndef WARTHOG_APRIORI_FILTER_H
#define WARTHOG_APRIORI_FILTER_H

// search/apriori_filter.h
//
// Sometimes during search it is desriable to ignore certain nodes that
// are identified apriori as not useful for any query. This filter helps to 
// achieve the goal by keeping track of a byte-sized label  for each node in a
// discrete graph. During search the filter can be queried about the 
// state of each node and those nodes whose bit is set can be ignored 
// (not generated or expanded).
//
// NB: This class duplicates the functionality of warthog::bitfifield_filter,
// trading space for a little bit of extra speed (no boolean logic ops).
// 
// 
// @author: dharabor
// @created: 2016-07-19

#include "constants.h"
#include "forward.h"
#include <cstdint>
#include <cstring>

namespace warthog
{

class apriori_filter 
{
    public:
        apriori_filter(uint32_t sz)
        {
            filter_sz_ = sz;
            filter_ = new uint8_t[filter_sz_];
        }

        ~apriori_filter()
        { delete [] filter_; }

        inline void
        add(uint32_t node_id)
        {
            filter_[node_id] = true;
        }

        inline void
        remove(uint32_t node_id)
        {
            filter_[node_id] = false;
        }

        inline bool
        filter(uint32_t node_id)
        {
            return filter_[node_id];
        }

        // clear all filter flags
        inline void
        reset_filter()
        {
            for(uint32_t i = 0; i < filter_sz_; i++)
            {
                filter_[i] = false;
            }
        }

        inline size_t
        mem()
        {
            return 
                sizeof(*filter_)*filter_sz_ +
                sizeof(this);
        }

    private:
        uint8_t* filter_;
        uint32_t filter_sz_;
};

}

#endif
