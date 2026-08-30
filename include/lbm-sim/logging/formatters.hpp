/**
 * @file formatters.hpp
 * @brief quill/fmt formatters for the library's own value types.
 *
 * Both Point and Vector already have an @c operator<<, so the formatters
 * derive from @c ostream_formatter and the codecs from
 * @c DeferredFormatCodec: the argument is copied into the log queue and
 * formatted later, on the backend thread.
 *
 * Included only when the quill backend is selected; the ostream backend
 * formats through @c operator<< directly and needs nothing here.
 */

#pragma once

#include "lbm-sim/core/point.hpp"
#include "lbm-sim/core/vector.hpp"

#include "quill/DeferredFormatCodec.h"
#include "quill/bundled/fmt/ostream.h"

#include "quill/core/Codec.h"

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
