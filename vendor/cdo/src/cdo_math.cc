#include "cdo_math.h"

#include <cmath>
#include <limits>

// "borrowed" from <linux/bitops.h> from linux-2.4 (lib/healpix)
// clang-format off
static unsigned int my_hweight32(unsigned int w) {
    unsigned int res = (w & 0x55555555) + ((w >> 1) & 0x55555555);
    res = (res & 0x33333333) + ((res >> 2) & 0x33333333);
    res = (res & 0x0F0F0F0F) + ((res >> 4) & 0x0F0F0F0F);
    res = (res & 0x00FF00FF) + ((res >> 8) & 0x00FF00FF);
    return (res & 0x0000FFFF) + ((res >> 16) & 0x0000FFFF);
}
// clang-format on

namespace cdo
{

unsigned int
is_power_of_two(unsigned int x)
{
  return (my_hweight32(x) == 1);
}

double
NaN()
{
  return (std::numeric_limits<double>::has_quiet_NaN) ? std::numeric_limits<double>::quiet_NaN() : std::acos(2.0);
}

}  // namespace cdo
