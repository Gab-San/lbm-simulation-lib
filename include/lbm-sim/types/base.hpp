#ifndef __LBM_SIM_TYPES_BASE_HPP
#define __LBM_SIM_TYPES_BASE_HPP

#include <cstdint>
#include <vector>

namespace lbm {
namespace types {

using boundary_t = uint8_t;
using boundary_mask_t = std::vector<boundary_t>;

using dim_t = unsigned short int;

} // namespace types
} // namespace lbm

#endif // __LBM_SIM_TYPES_BASE_HPP
