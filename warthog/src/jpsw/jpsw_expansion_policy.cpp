#include "jpsw_expansion_policy.h"
#include "jps.h"

struct jumpcache
{
    uint32_t version_;
    uint32_t target_;
    warthog::cost_t g_;
};

struct warthog::jpsw_extra
{
    warthog::cost_t prospective_g_;
    uint32_t search_number_;
    uint8_t successors_;
    bool reached_orthogonally_;

    uint8_t successor_sets_[8];
    jumpcache jump_cache_[4];
};

warthog::jpsw_expansion_policy::jpsw_expansion_policy(
        nbcache& nbcache, vl_gridmap& map, cost_table& costs)
    : expansion_policy(map.height()*map.width()),
    map_(map), costs_(costs), nbcache_(nbcache)
{
    extra_ = new jpsw_extra[map.width() * map.height()];
    for (uint32_t id = 0; id < map.width()*map.height(); id++) {
        for (int i = 0; i < 8; i++) {
            extra_[id].successor_sets_[i] = warthog::jps::ALL;
        }
        for (int i = 0; i < 4; i++) {
            extra_[id].jump_cache_[i].version_ = 0;
        }
    }

    row_versions_ = new uint32_t[map.height()];
    for (uint32_t i = 0; i < map.height(); i++) {
        row_versions_[i] = 1;
    }
    col_versions_ = new uint32_t[map.width()];
    for (uint32_t i = 0; i < map.width(); i++) {
        col_versions_[i] = 1;
    }
}

warthog::jpsw_expansion_policy::~jpsw_expansion_policy()
{
    delete[] extra_;
    delete[] row_versions_;
    delete[] col_versions_;
}

void warthog::jpsw_expansion_policy::fill_nb_cache()
{
    for (uint32_t id = 0; id < map_.width() * map_.height(); id++) {
        if (!map_.get_label(id)) {
            continue;
        }

        for (int i = 0; i < 8; i++) {
            nbhood_successors(id, (warthog::jps::direction) (1 << i));
        }
    }
}

void warthog::jpsw_expansion_policy::fill_jump_cache()
{
    for (uint32_t id = 0; id < map_.width() * map_.height(); id++) {
        if (!map_.get_label(id)) {
            continue;
        }

        auto row_version = row_versions_[id / map_.width()];
        calculate_jump<0>(id, -1, row_version, [this](uint32_t h) {
            return horizontal_cost(h - 1);
        });
        calculate_jump<1>(id, 1, row_version, [this](uint32_t h) {
            return horizontal_cost(h);
        });

        auto col_version = col_versions_[id % map_.width()];
        calculate_jump<2>(id, -(int) map_.width(), col_version, [this](uint32_t h) {
            return vertical_cost(h - map_.width());
        });
        calculate_jump<3>(id, map_.width(), col_version, [this](uint32_t h) {
            return vertical_cost(h);
        });
    }
}

size_t warthog::jpsw_expansion_policy::mem()
{
    return expansion_policy::mem() + sizeof(*this) + map_.mem() + nbcache_.mem()
            + sizeof(jpsw_extra) * map_.width() * map_.height()
            + sizeof(uint32_t) * (map_.width() + map_.height());
}

warthog::search_node*
warthog::jpsw_expansion_policy::generate_start_node(warthog::problem_instance* pi)
{
    uint32_t max_id = map_.header_width() * map_.header_height();
    if ((uint32_t) pi->start_id_ >= max_id) {
        return 0;
    }
    uint32_t padded_id = map_.to_padded_id((uint32_t) pi->start_id_);
    if (map_.get_label(padded_id) == 0) {
        return 0;
    }

    nbhood_labels nb = nbhood(padded_id);

    // Find initial successor set
    auto& extra = extra_[padded_id];
    extra.search_number_ = pi->instance_id_;
    extra.successors_ = 0;
    extra.prospective_g_ = 0.0;
    extra.reached_orthogonally_ = true;
    if (costs_[map_.get_label(nb.n)]) {
        extra.successors_ |= warthog::jps::NORTH;
        prospect(nb.n, vertical_cost(nb.n), true, pi);
        if (costs_[map_.get_label(nb.e)] && costs_[map_.get_label(nb.ne)]) {
            extra.successors_ |= warthog::jps::NORTHEAST;
            prospect(nb.ne, diagonal_cost(nb.n), false, pi);
        }
        if (costs_[map_.get_label(nb.w)] && costs_[map_.get_label(nb.nw)]) {
            extra.successors_ |= warthog::jps::NORTHWEST;
            prospect(nb.nw, diagonal_cost(nb.nw), false, pi);
        }
    }
    if (costs_[map_.get_label(nb.s)]) {
        extra.successors_ |= warthog::jps::SOUTH;
        prospect(nb.s, vertical_cost(nb.h), true, pi);
        if (costs_[map_.get_label(nb.e)] && costs_[map_.get_label(nb.se)]) {
            extra.successors_ |= warthog::jps::SOUTHEAST;
            prospect(nb.se, diagonal_cost(nb.h), false, pi);
        }
        if (costs_[map_.get_label(nb.w)] && costs_[map_.get_label(nb.sw)]) {
            extra.successors_ |= warthog::jps::SOUTHWEST;
            prospect(nb.sw, diagonal_cost(nb.w), false, pi);
        }
    }
    if (costs_[map_.get_label(nb.e)]) {
        extra.successors_ |= warthog::jps::EAST;
        prospect(nb.e, horizontal_cost(nb.h), true, pi);
    }
    if (costs_[map_.get_label(nb.w)]) {
        extra.successors_ |= warthog::jps::WEST;
        prospect(nb.w, horizontal_cost(nb.w), true, pi);
    }

    return generate(padded_id);
}

warthog::search_node*
warthog::jpsw_expansion_policy::generate_target_node(warthog::problem_instance* pi)
{
    uint32_t max_id = map_.header_width() * map_.header_height();
    if ((uint32_t) pi->target_id_ >= max_id) {
        return 0;
    }
    uint32_t padded_id = map_.to_padded_id((uint32_t) pi->target_id_);
    if (map_.get_label(padded_id) == 0) {
        return 0;
    }
    return generate(padded_id);
}

void
warthog::jpsw_expansion_policy::expand(warthog::search_node* node, warthog::problem_instance* pi)
{
    reset();

    uint32_t id = node->get_id();
    auto nb = nbhood(id);

    uint8_t successors = extra_[id].successors_;

    // do orthogonal directions first so successors are updated for diagonal branch pruning
    if (successors & warthog::jps::NORTH) {
        // check that the successor is not pruned by prospective g
        if (extra_[nb.n].prospective_g_ == node->get_g() + vertical_cost(nb.n)) {
            jump_north(id, nb, node->get_g(), 0.0, pi);
        } else {
            successors &= ~warthog::jps::NORTH;
        }
    }
    if (successors & warthog::jps::WEST) {
        // check that the successor is not pruned by prospective g
        if (extra_[nb.w].prospective_g_ == node->get_g() + horizontal_cost(nb.w)) {
            jump_west(id, nb, node->get_g(), 0.0, pi);
        } else {
            successors &= ~warthog::jps::WEST;
        }
    }
    if (successors & warthog::jps::EAST) {
        // check that the successor is not pruned by prospective g
        if (extra_[nb.e].prospective_g_ == node->get_g() + horizontal_cost(nb.h)) {
            jump_east(id, nb, node->get_g(), 0.0, pi);
        } else {
            successors &= ~warthog::jps::EAST;
        }
    }
    if (successors & warthog::jps::SOUTH) {
        // check that the successor is not pruned by prospective g
        if (extra_[nb.s].prospective_g_ == node->get_g() + vertical_cost(nb.h)) {
            jump_south(id, nb, node->get_g(), 0.0, pi);
        } else {
            successors &= ~warthog::jps::SOUTH;
        }
    }

    // now that orthogonal successors are updated, do the diagonal directions
    if (successors & warthog::jps::NORTHWEST) {
        // check that the successor is not pruned by prospective g
        if (extra_[nb.nw].prospective_g_ == node->get_g() + diagonal_cost(nb.nw)
                && !extra_[nb.nw].reached_orthogonally_) {
            jump_nw(nb, node->get_g(), successors, pi);
        }
    }
    if (successors & warthog::jps::NORTHEAST) {
        // check that the successor is not pruned by prospective g
        if (extra_[nb.ne].prospective_g_ == node->get_g() + diagonal_cost(nb.n)
                && !extra_[nb.ne].reached_orthogonally_) {
            jump_ne(nb, node->get_g(), successors, pi);
        }
    }
    if (successors & warthog::jps::SOUTHWEST) {
        // check that the successor is not pruned by prospective g
        if (extra_[nb.sw].prospective_g_ == node->get_g() + diagonal_cost(nb.w)
                && !extra_[nb.sw].reached_orthogonally_) {
            jump_sw(nb, node->get_g(), successors, pi);
        }
    }
    if (successors & warthog::jps::SOUTHEAST) {
        // check that the successor is not pruned by prospective g
        if (extra_[nb.se].prospective_g_ == node->get_g() + diagonal_cost(nb.h)
                && !extra_[nb.se].reached_orthogonally_) {
            jump_se(nb, node->get_g(), successors, pi);
        }
    }
}

template<int SLOT, typename F>
void warthog::jpsw_expansion_policy::calculate_jump(
    uint32_t start, int delta, uint32_t version, F cost)
{
    auto jump_dist = 0;
    uint32_t last_id = 0;
    auto jump_cost = 0.0;
    auto nb = nbhood(start);
    do {
        jump_dist += 1;
        nb = nbhood(nb.h + delta);
        last_id = nb.h;
        if (!locally_uniform(nb)) {
            break;
        }
        auto& jc = extra_[nb.h].jump_cache_[SLOT];
        if (jc.version_ == version) {
            jump_cost = jc.g_;
            nb = nbhood(jc.target_);
            break;
        }
    } while (true);

    for (int i = 0; i < jump_dist; i++) {
        last_id -= delta;
        jump_cost += cost(last_id);
        extra_[last_id].jump_cache_[SLOT].target_ = nb.h;
        extra_[last_id].jump_cache_[SLOT].g_ = jump_cost;
        extra_[last_id].jump_cache_[SLOT].version_ = version;
    }

    assert(last_id == start);
}

void warthog::jpsw_expansion_policy::jump_west(
        uint32_t from, nbhood_labels nb, double g, double cost, warthog::problem_instance* pi)
{
    auto start = nb.h;
    auto version = row_versions_[start / map_.width()];
    if (extra_[nb.h].jump_cache_[0].version_ != version) {
        calculate_jump<0>(start, -1, version, [this](uint32_t h) {
            return horizontal_cost(h - 1);
        });
    }
    cost += extra_[nb.h].jump_cache_[0].g_;
    nb = nbhood(extra_[nb.h].jump_cache_[0].target_);

    if (start > pi->target_id_ && nb.h < pi->target_id_) {
        // Jump overshoots target, so adjust to target
        nb = nbhood(pi->target_id_);
        cost -= extra_[nb.h].jump_cache_[0].g_;
    }

    if (nb.h == pi->target_id_ || nbhood_successors(nb.h, warthog::jps::WEST) != 0) {
        add_neighbour(generate(nb.h), cost);
        reach(nb.h, warthog::jps::WEST, g + cost, pi);
    }
}

void warthog::jpsw_expansion_policy::jump_east(
        uint32_t from, nbhood_labels nb, double g, double cost, warthog::problem_instance* pi)
{
    auto start = nb.h;
    auto version = row_versions_[start / map_.width()];
    if (extra_[nb.h].jump_cache_[1].version_ != version) {
        calculate_jump<1>(start, 1, version, [this](uint32_t h) {
            return horizontal_cost(h);
        });
    }
    cost += extra_[nb.h].jump_cache_[1].g_;
    nb = nbhood(extra_[nb.h].jump_cache_[1].target_);

    if (start < pi->target_id_ && nb.h > pi->target_id_) {
        // Jump overshoots target, so adjust to target
        nb = nbhood(pi->target_id_);
        cost -= extra_[nb.h].jump_cache_[1].g_;
    }

    if (nb.h == pi->target_id_ || nbhood_successors(nb.h, warthog::jps::EAST) != 0) {
        add_neighbour(generate(nb.h), cost);
        reach(nb.h, warthog::jps::EAST, g + cost, pi);
    }
}

void warthog::jpsw_expansion_policy::jump_north(
        uint32_t from, nbhood_labels nb, double g, double cost, warthog::problem_instance* pi)
{
    auto start = nb.h;
    auto version = col_versions_[start % map_.width()];
    if (extra_[nb.h].jump_cache_[2].version_ != version) {
        calculate_jump<2>(start, -(int) map_.width(), version, [this](uint32_t h) {
            return vertical_cost(h - map_.width());
        });
    }
    cost += extra_[nb.h].jump_cache_[2].g_;
    nb = nbhood(extra_[nb.h].jump_cache_[2].target_);

    if (start % map_.width() == pi->target_id_ % map_.width() && start > pi->target_id_ && nb.h < pi->target_id_) {
        // Jump overshoots target, so adjust to target
        nb = nbhood(pi->target_id_);
        cost -= extra_[nb.h].jump_cache_[2].g_;
    }

    if (nb.h == pi->target_id_ || nbhood_successors(nb.h, warthog::jps::NORTH) != 0) {
        add_neighbour(generate(nb.h), cost);
        reach(nb.h, warthog::jps::NORTH, g + cost, pi);
    }
}

void warthog::jpsw_expansion_policy::jump_south(
        uint32_t from, nbhood_labels nb, double g, double cost, warthog::problem_instance* pi)
{
    auto start = nb.h;
    auto version = col_versions_[start % map_.width()];
    if (extra_[nb.h].jump_cache_[3].version_ != version) {
        calculate_jump<3>(start, map_.width(), version, [this](uint32_t h) {
            return vertical_cost(h);
        });
    }
    cost += extra_[nb.h].jump_cache_[3].g_;
    nb = nbhood(extra_[nb.h].jump_cache_[3].target_);

    if (start % map_.width() == pi->target_id_ % map_.width() && start < pi->target_id_ && nb.h > pi->target_id_) {
        // Jump overshoots target, so adjust to target
        nb = nbhood(pi->target_id_);
        cost -= extra_[nb.h].jump_cache_[3].g_;
    }
    if (nb.h == pi->target_id_ || nbhood_successors(nb.h, warthog::jps::SOUTH) != 0) {
        add_neighbour(generate(nb.h), cost);
        reach(nb.h, warthog::jps::SOUTH, g + cost, pi);
    }
}

void warthog::jpsw_expansion_policy::jump_nw(
        nbhood_labels nb, double g, int successor_set, warthog::problem_instance* pi)
{
    uint32_t from = nb.h;
    double cost = diagonal_cost(nb.nw);
    nb = nbhood(nb.nw);
    while (locally_uniform(nb) && nb.h != pi->target_id_) {
        if (successor_set & warthog::jps::NORTH) {
            jump_north(from, nb, g, cost, pi);
        }
        if (successor_set & warthog::jps::WEST) {
            jump_west(from, nb, g, cost, pi);
        }
        cost += diagonal_cost(nb.nw);
        nb = nbhood(nb.nw);
    }
    add_neighbour(generate(nb.h), cost);
    reach(nb.h, warthog::jps::NORTHWEST, g + cost, pi);
}

void warthog::jpsw_expansion_policy::jump_ne(
        nbhood_labels nb, double g, int successor_set, warthog::problem_instance* pi)
{
    uint32_t from = nb.h;
    double cost = diagonal_cost(nb.n);
    nb = nbhood(nb.ne);
    while (locally_uniform(nb) && nb.h != pi->target_id_) {
        if (successor_set & warthog::jps::NORTH) {
            jump_north(from, nb, g, cost, pi);
        }
        if (successor_set & warthog::jps::EAST) {
            jump_east(from, nb, g, cost, pi);
        }
        cost += diagonal_cost(nb.n);
        nb = nbhood(nb.ne);
    }
    add_neighbour(generate(nb.h), cost);
    reach(nb.h, warthog::jps::NORTHEAST, g + cost, pi);
}

void warthog::jpsw_expansion_policy::jump_sw(
        nbhood_labels nb, double g, int successor_set, warthog::problem_instance* pi)
{
    uint32_t from = nb.h;
    double cost = diagonal_cost(nb.w);
    nb = nbhood(nb.sw);
    while (locally_uniform(nb) && nb.h != pi->target_id_) {
        if (successor_set & warthog::jps::SOUTH) {
            jump_south(from, nb, g, cost, pi);
        }
        if (successor_set & warthog::jps::WEST) {
            jump_west(from, nb, g, cost, pi);
        }
        cost += diagonal_cost(nb.w);
        nb = nbhood(nb.sw);
    }
    add_neighbour(generate(nb.h), cost);
    reach(nb.h, warthog::jps::SOUTHWEST, g + cost, pi);
}

void warthog::jpsw_expansion_policy::jump_se(
        nbhood_labels nb, double g, int successor_set, warthog::problem_instance* pi)
{
    uint32_t from = nb.h;
    double cost = diagonal_cost(nb.h);
    nb = nbhood(nb.se);
    while (locally_uniform(nb) && nb.h != pi->target_id_) {
        if (successor_set & warthog::jps::SOUTH) {
            jump_south(from, nb, g, cost, pi);
        }
        if (successor_set & warthog::jps::EAST) {
            jump_east(from, nb, g, cost, pi);
        }
        cost += diagonal_cost(nb.h);
        nb = nbhood(nb.se);
    }
    add_neighbour(generate(nb.h), cost);
    reach(nb.h, warthog::jps::SOUTHEAST, g + cost, pi);
}

void warthog::jpsw_expansion_policy::reach(
        uint32_t id, warthog::jps::direction direction, double g, warthog::problem_instance* pi)
{
    constexpr auto ORTHO_DIRS =
        warthog::jps::NORTH | warthog::jps::SOUTH | warthog::jps::EAST | warthog::jps::WEST;

    auto& extra = extra_[id];
    auto nb = nbhood(id);
    warthog::search_node* node = generate(id);
    if (pi->instance_id_ != node->get_search_number() || g < node->get_g()) {
        // determine successor set, update prospective g for the reached node
        extra.successors_ = nbhood_successors(id, direction);
        prospect(id, g, (direction & ORTHO_DIRS) != 0, pi);

        // for each direction, update the prospective g value of the node reached
        if (extra.successors_ & warthog::jps::NORTHWEST) {
            prospect(nb.nw, g + diagonal_cost(nb.nw), false, pi);
        }
        if (extra.successors_ & warthog::jps::NORTH) {
            prospect(nb.n, g + vertical_cost(nb.n), true, pi);
        }
        if (extra.successors_ & warthog::jps::NORTHEAST) {
            prospect(nb.ne, g + diagonal_cost(nb.n), false, pi);
        }

        if (extra.successors_ & warthog::jps::WEST) {
            prospect(nb.w, g + horizontal_cost(nb.w), true, pi);
        }
        if (extra.successors_ & warthog::jps::EAST) {
            prospect(nb.e, g + horizontal_cost(nb.h), true, pi);
        }

        if (extra.successors_ & warthog::jps::SOUTHWEST) {
            prospect(nb.sw, g + diagonal_cost(nb.w), false, pi);
        }
        if (extra.successors_ & warthog::jps::SOUTH) {
            prospect(nb.s, g + vertical_cost(nb.h), true, pi);
        }
        if (extra.successors_ & warthog::jps::SOUTHEAST) {
            prospect(nb.se, g + diagonal_cost(nb.h), false, pi);
        }
    } else if (pi->instance_id_ == node->get_search_number() && g == node->get_g()) {
        // when there are multiple canonical optimal paths to a node, its successors are the
        // intersection of the canonical successor sets in each direction.
        extra.successors_ &= nbhood_successors(id, direction);
        extra.reached_orthogonally_ |= (direction & ORTHO_DIRS) != 0;
    }
}

void warthog::jpsw_expansion_policy::prospect(
        uint32_t id, double g, bool ortho, warthog::problem_instance* pi)
{
    auto& extra = extra_[id];
    if (extra.search_number_ != pi->instance_id_ || g < extra.prospective_g_) {
        extra.search_number_ = pi->instance_id_;
        extra.prospective_g_ = g;
        extra.reached_orthogonally_ = ortho;
    } else if (extra.search_number_ == pi->instance_id_ && g == extra.prospective_g_) {
        extra.reached_orthogonally_ |= ortho;
    }
}

int warthog::jpsw_expansion_policy::nbhood_successors(
        uint32_t to, warthog::jps::direction going)
{
    nbhood_labels nb = nbhood(to);
    uint8_t* successor_set;
    switch (going) {
        case warthog::jps::NORTHWEST:
            successor_set = &extra_[nb.h].successor_sets_[0];
            break;
        case warthog::jps::NORTH:
            successor_set = &extra_[nb.h].successor_sets_[1];
            break;
        case warthog::jps::NORTHEAST:
            successor_set = &extra_[nb.h].successor_sets_[2];
            break;
        case warthog::jps::WEST:
            successor_set = &extra_[nb.h].successor_sets_[3];
            break;
        case warthog::jps::EAST:
            successor_set = &extra_[nb.h].successor_sets_[4];
            break;
        case warthog::jps::SOUTHWEST:
            successor_set = &extra_[nb.h].successor_sets_[5];
            break;
        case warthog::jps::SOUTH:
            successor_set = &extra_[nb.h].successor_sets_[6];
            break;
        case warthog::jps::SOUTHEAST:
            successor_set = &extra_[nb.h].successor_sets_[7];
            break;
        default:
            assert(false);
            return 0;
    }
    if (*successor_set == warthog::jps::ALL) {
        *successor_set = nbcache_.successors(map_, nb, going);
    }
    return *successor_set;
}
