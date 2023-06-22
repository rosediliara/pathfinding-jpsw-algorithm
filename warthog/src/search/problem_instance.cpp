#include "problem_instance.h"

uint32_t warthog::problem_instance::instance_counter_ = 1;

std::ostream& operator<<(std::ostream& str, warthog::problem_instance& pi)
{
    pi.print(str);

    return str;
}
