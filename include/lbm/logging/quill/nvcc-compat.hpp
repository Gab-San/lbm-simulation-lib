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
