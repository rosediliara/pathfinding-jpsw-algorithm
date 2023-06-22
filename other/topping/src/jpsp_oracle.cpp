#include "jpsp_oracle.h"
#include "online_jump_point_locator.h"

warthog::jpsp_oracle::jpsp_oracle(warthog::gridmap* gm)
    : gm_(gm)
{
    oracle_ = 0;
    preprocess();
}

warthog::jpsp_oracle::~jpsp_oracle()
{
    delete oracle_;
}

void
warthog::jpsp_oracle::preprocess()
{
    oracle_ = new offline_jump_point_locator(gm_);
}

bool
warthog::jpsp_oracle::sanity_check()
{
    bool retval = true;
    warthog::online_jump_point_locator jpl(gm_);
    for(uint32_t y = 0; y < gm_->header_height(); y++)
    {
        for(uint32_t x = 0; x < gm_->header_width(); x++)
        {
            uint32_t loc_id = gm_->to_padded_id(x, y);
            if(!gm_->get_label(loc_id)) { continue; }

            for(uint32_t i = 0; i < 8; i++)
            {
                warthog::jps::direction d = (warthog::jps::direction)(1<<i);
                uint32_t jp_id;
                warthog::cost_t jp_cost;
                jpl.jump(d, loc_id, warthog::INF, jp_id, jp_cost);

                if(jp_id != warthog::INF)
                {
                    uint32_t jp_x, jp_y;
                    gm_->to_unpadded_xy(jp_id, jp_x, jp_y);
                    jp_id = jp_y * gm_->header_width() + jp_x;
                }
                uint32_t oracle_value = this->next_jump_point(x, y, d);
                if(oracle_value != jp_cost)
                {
                    std::cerr 
                        << "sanity fail at "
                        << "("<< x <<", "<< y <<") dir " << d << "; "
                        << " oracle="<<oracle_value << " jpl="<< jp_cost
                        << "\n";
                    return false;
                }
                retval = retval && (oracle_value == jp_cost);
            }
        }
    }
    return retval;
}
