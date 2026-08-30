/**
 * @file nvcc-compat.hpp
 * @brief Works around an nvcc parse failure in quill's bundled fmt.
 *
 * fmt's 128-bit integer helper declares @c operator~ in a way nvcc does not
 * accept; the definition here supplies it. Active only under @c __CUDACC__,
 * and a no-op for every other compiler.
 */

#ifdef __CUDACC__
#include "quill/bundled/fmt/format.h"

namespace fmtquill {
inline namespace v12 {
namespace detail {
constexpr auto operator~(uint128 n) -> uint128 { return {~n.high(), ~n.low()}; }
} // namespace detail
} // namespace v12
} // namespace fmtquill
#endif
