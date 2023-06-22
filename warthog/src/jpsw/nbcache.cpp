#include "nbcache.h"
#include <algorithm>

uint8_t temp[8];

warthog::nbcache::nbcache(warthog::cost_table& costs)
    : local_map_(3, 3), expander_(&local_map_, costs), pqueue_(), cached_()
{
    local_nb_.nw = local_map_.to_padded_id(0, 0);
    local_nb_.n = local_map_.to_padded_id(1, 0);
    local_nb_.ne = local_map_.to_padded_id(2, 0);
    local_nb_.w = local_map_.to_padded_id(0, 1);
    local_nb_.h = local_map_.to_padded_id(1, 1);
    local_nb_.e = local_map_.to_padded_id(2, 1);
    local_nb_.sw = local_map_.to_padded_id(0, 2);
    local_nb_.s = local_map_.to_padded_id(1, 2);
    local_nb_.se = local_map_.to_padded_id(2, 2);
}

void rotate_left(std::array<uint8_t, 9>& cells, int count)
{
    // uint8_t temp[8];
    if (count > 2)
    {
        for (int i = 0; i < 8; i++)
        {
            temp[i] = cells[i];
        }
        for (int i = 0; i < 8; i++)
        {
            cells[(i - count + 8) % 8] = temp[i];
        }
    }
    else
    {
        for (int i = 0; i < count; i++) 
        {
            uint8_t first = cells[0];
            for (int j = 1; j < 8; j++) {
                cells[j - 1] = cells[j];
            }
            cells[7] = first;
        }
    }
    
    
    
    // std::rotate(cells.begin(), cells.begin() + count, cells.end());
}

void rotate(warthog::nb_key& key, warthog::jps::direction going)
{
    switch (going) {
        case warthog::jps::NORTH:
            key.diagonal = false;
            break;
        case warthog::jps::NORTHWEST:
            key.diagonal = true;
            break;
        case warthog::jps::EAST:
            key.diagonal = false;
            rotate_left(key.cells, 2);
            break;
        case warthog::jps::NORTHEAST:
            key.diagonal = true;
            rotate_left(key.cells, 2);
            break;
        case warthog::jps::SOUTH:
            key.diagonal = false;
            rotate_left(key.cells, 4);
            break;
        case warthog::jps::SOUTHEAST:
            key.diagonal = true;
            rotate_left(key.cells, 4);
            break;
        case warthog::jps::WEST:
            key.diagonal = false;
            rotate_left(key.cells, 6);
            break;
        case warthog::jps::SOUTHWEST:
            key.diagonal = true;
            rotate_left(key.cells, 6);
            break;
        default:
            assert(false);
    }
}

uint8_t warthog::nbcache::successors(
        vl_gridmap& map, nbhood_labels& nb, warthog::jps::direction going)
{
    nb_key key = {
        {
            map.get_label(nb.nw),
            map.get_label(nb.n),
            map.get_label(nb.ne),
            map.get_label(nb.e),
            map.get_label(nb.se),
            map.get_label(nb.s),
            map.get_label(nb.sw),
            map.get_label(nb.w),
            map.get_label(nb.h),
        },
        false
    };
    rotate(key, going);

    auto size = cached_.size();
    auto& set = cached_[key];
    if (cached_.size() != size) {
        local_map_.set_label(local_nb_.nw, key.cells[0]);
        local_map_.set_label(local_nb_.n,  key.cells[1]);
        local_map_.set_label(local_nb_.ne, key.cells[2]);
        local_map_.set_label(local_nb_.e,  key.cells[3]);
        local_map_.set_label(local_nb_.se, key.cells[4]);
        local_map_.set_label(local_nb_.s,  key.cells[5]);
        local_map_.set_label(local_nb_.sw, key.cells[6]);
        local_map_.set_label(local_nb_.w,  key.cells[7]);
        local_map_.set_label(local_nb_.h,  key.cells[8]);
        set = calculate_successors(key.diagonal ? 8 : 7);
    }

    switch (going) {
        case warthog::jps::NORTH:
        case warthog::jps::NORTHWEST:
            return set;
        case warthog::jps::EAST:
        case warthog::jps::NORTHEAST:
            return ((set >> 3) & 0b00000001)
                 | ((set >> 2) & 0b00100000)
                 | ((set >> 1) & 0b00010010)
                 | ((set << 1) & 0b10000000)
                 | ((set << 2) & 0b01001100);
        case warthog::jps::SOUTH:
        case warthog::jps::SOUTHEAST:
            return ((set >> 3) & 0b00010000)
                 | ((set >> 1) & 0b00100101)
                 | ((set << 1) & 0b01001010)
                 | ((set << 3) & 0b10000000);
        case warthog::jps::WEST:
        case warthog::jps::SOUTHWEST:
            return ((set >> 2) & 0b00010011)
                 | ((set >> 1) & 0b01000000)
                 | ((set << 1) & 0b00100100)
                 | ((set << 2) & 0b10000000)
                 | ((set << 3) & 0b00001000);
        default:;
    }
    assert(false);
    return 0;
}

uint8_t warthog::nbcache::calculate_successors(uint32_t source)
{
    warthog::problem_instance pi(source);

    warthog::search_node* start = expander_.generate_start_node(&pi);
    pi.start_id_ = start->get_id();
    start->init(pi.instance_id_, warthog::SN_ID_MAX, 0, 0);
    pqueue_.push(start);

    while (pqueue_.size()) {
        warthog::search_node* current = pqueue_.pop();
        current->set_expanded(true);
        expander_.expand(current, &pi);

        warthog::search_node* n;
        warthog::cost_t cost_to_n;

        for (expander_.first(n, cost_to_n); n != 0; expander_.next(n, cost_to_n)) {
            warthog::cost_t gval = current->get_g() + cost_to_n;
            if (n->get_search_number() != current->get_search_number()) {
                n->init(current->get_search_number(), current->get_id(), gval, gval);
                pqueue_.push(n);
            } else if (gval < n->get_g()) {
                n->relax(gval, current->get_id());
                pqueue_.decrease_key(n);
            } else if (gval == n->get_g()) {
                uint32_t h = n->get_id();
                uint32_t p = n->get_parent();
                bool existing_is_ortho =
                    p == h - 1 || p == h - local_map_.width() ||
                    p == h + 1 || p == h + local_map_.width();
                uint32_t c = current->get_id();
                bool new_is_ortho =
                    c == h - 1 || c == h - local_map_.width() ||
                    c == h + 1 || c == h + local_map_.width();
                if (
                    (new_is_ortho && !existing_is_ortho) ||
                    (new_is_ortho == existing_is_ortho && c == local_nb_.h)
                ) {
                    // tiebreak in favor of ortho-last
                    n->set_parent(current->get_id());
                }
            }
        }
    }

    uint8_t successors = warthog::jps::NONE;
    warthog::search_node* n = expander_.generate(local_nb_.nw);
    if (n->get_search_number() == pi.instance_id_ && n->get_parent() == local_nb_.h) {
        successors |= warthog::jps::NORTHWEST;
    }
    n = expander_.generate(local_nb_.n);
    if (n->get_search_number() == pi.instance_id_ && n->get_parent() == local_nb_.h) {
        successors |= warthog::jps::NORTH;
    }
    n = expander_.generate(local_nb_.ne);
    if (n->get_search_number() == pi.instance_id_ && n->get_parent() == local_nb_.h) {
        successors |= warthog::jps::NORTHEAST;
    }
    n = expander_.generate(local_nb_.w);
    if (n->get_search_number() == pi.instance_id_ && n->get_parent() == local_nb_.h) {
        successors |= warthog::jps::WEST;
    }
    n = expander_.generate(local_nb_.e);
    if (n->get_search_number() == pi.instance_id_ && n->get_parent() == local_nb_.h) {
        successors |= warthog::jps::EAST;
    }
    n = expander_.generate(local_nb_.sw);
    if (n->get_search_number() == pi.instance_id_ && n->get_parent() == local_nb_.h) {
        successors |= warthog::jps::SOUTHWEST;
    }
    n = expander_.generate(local_nb_.s);
    if (n->get_search_number() == pi.instance_id_ && n->get_parent() == local_nb_.h) {
        successors |= warthog::jps::SOUTH;
    }
    n = expander_.generate(local_nb_.se);
    if (n->get_search_number() == pi.instance_id_ && n->get_parent() == local_nb_.h) {
        successors |= warthog::jps::SOUTHEAST;
    }
    return successors;
}
