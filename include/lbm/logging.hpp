#ifndef __LBM_LOGGING_HPP
#define __LBM_LOGGING_HPP

#include "lbm-sim/core/point.hpp"
#include "lbm-sim/core/vector.hpp"

#include "quill/Logger.h"

#include "quill/DeferredFormatCodec.h"
#include "quill/bundled/fmt/ostream.h"

#include "quill/core/Codec.h"

#include "quill/sinks/Sink.h"

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

template <typename T, unsigned short int dim>
struct fmtquill::formatter<lbm::utils::Point<T, dim>>
    : fmtquill::ostream_formatter {};

template <typename T, unsigned short int dim>
struct quill::Codec<lbm::utils::Point<T, dim>>
    : quill::DeferredFormatCodec<lbm::utils::Point<T, dim>> {};

template <typename T, unsigned short int dim>
struct fmtquill::formatter<lbm::utils::Vector<T, dim>>
    : fmtquill::ostream_formatter {};

template <typename T, unsigned short int dim>
struct quill::Codec<lbm::utils::Vector<T, dim>>
    : quill::DeferredFormatCodec<lbm::utils::Vector<T, dim>> {};

namespace lbm {
namespace logging {
void setup_quill();

std::shared_ptr<quill::Sink> get_console_sink();

quill::Logger *create_or_get_logger(const std::string &logger_name);
} // namespace logging
} // namespace lbm

#endif // __LBM_LOGGING_HPP
